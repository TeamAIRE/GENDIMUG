from utils.classes.sorting.cdsprocessor import CdsProcessor
from utils.classes.storage.ogtree import OgTree
from utils.helpers.id_functions import get_prefix


class GffSorter:
    """
    Sorts the gff data and creates multiple storage objects (stored in a ContigStorage object)
    for feature data and ontogenetic links extracted from the corresponding gff3 line,
    with the following transformations:

        - start and stop positions are converted to int
        - the last column of attributes is converted to a dict {"ID": id, "Parent": xxx ...}

        Relational and feature attribute data are stored using these objects:

        -> ContigStorage object carrying data for the current contig in the following objects:
            > Gff: stores each gff line keyed by feature id and with tuple of gff data as values
            > FeatureTypes : dict {str: set of feature ids} for particular set of ids
        -> cds_storage : CdsStorage object that manages CDS fragmentation, renaming of fragments with:
            > CdsProcessor object
        -> contig_description : dict that stores contig data from feature "region" or "chromosome" lines

    Attributes
    ----------
        contig_storage: ContigStorage object, initialized for the provided database and contig ID
        database: str, database name
        contig: str, current contig ID
        cds_storage: CdsStorage object
        contig_description: dict {contig ID: contig data tuple}
        child_parent_dict: dict {child feature: parent feature}
    """
    def __init__(self, contig_storage, cds_storage, contig_description, keep_children=True):
        self.contig_storage = contig_storage
        self.database = contig_storage.database
        self.contig = contig_storage.contig
        self.cds_storage = cds_storage
        self.contig_description = contig_description
        self.keep_children = keep_children
        # temporary storage of the links child to parent (from the Parent flag in the gff attributes)
        self.child_parent_dict = {}

    def sort(self, list_gff_data):
        """ Initializes and builds storage objects.
        Returns ContigStorage, CdsStorage, child_parent_dict and contig_description dict"""

        self.contig_storage.feature_types.initialize("pseudogenes",
                                                     "filtered",
                                                     "additional_features")

        self.cds_storage.contig_cds[self.contig_storage.contig] = set()

        # sort features to extract contig and cds data as well as child-parent links
        self._triage(list_gff_data)

        # Renumber fragmented cds
        if self.keep_children:
            self._add_suffixes_to_cds_parts()

        # Delete obsolete ids
        if self.cds_storage.cds_to_delete and self.keep_children:
            self.contig_storage.gff.remove_set(self.cds_storage.cds_to_delete)

        return self.contig_storage, self.cds_storage, self.child_parent_dict, self.contig_description

    def _triage(self, list_gff_data):
        """ Sorts the various features based on their attributes. Creates a child_parent_dict
        {child feature: parent feature} """
        database = self.contig_storage.database
        for data in list_gff_data:
            self.contig_storage.gff.update(data)
            _, feature, list_positions, _, _, attributes = data.data_tuple
            feature_id = data.id

            if feature in {"region", "chromosome"}:  # marker for contig feature in either refseq/genbank
                # or ensembl GFF, respectively
                for (contig, strand, start, stop) in list_positions:
                    self._update_contig_data(feature_id,
                                             start,
                                             stop,
                                             attributes)

            elif feature.lower() == "cds" and self.keep_children:
                # cds ids will be reformatted in the next step, as there is no 'part' attribute for
                # fragmented cds defined by exonic positions
                self.cds_storage.cds.add(feature_id)

            else:
                # all other features in the gff except contigs and trans-spliced genes
                if self._is_pseudogene(feature_id):
                    self.contig_storage.feature_types.update("pseudogenes", {feature_id})
                self.child_parent_dict = OgTree.update_child_parent_dict(self.child_parent_dict,
                                                                         feature_id,
                                                                         attributes,
                                                                         database=database)

    def _add_suffixes_to_cds_parts(self):
        """ Splits multi-position CDS features into numbered parts matched to their parent exons.
        Single-position CDS features are renamed with a '-1' suffix for consistency.
        Original CDS IDs are marked for deletion """

        cds_from_contig = self.contig_storage.gff.get_features_for_contig_of_type(contig=self.contig,
                                                                                  features={"cds", "cdspart"})
        for cds in cds_from_contig:
            self.cds_storage.contig_cds[self.contig].add(cds)
            cds_data = self.contig_storage.gff.data(cds).data_tuple
            cds_processor = CdsProcessor(cds, cds_data)
            self.child_parent_dict, suffixes = cds_processor.process_cds(self.contig_storage, self.child_parent_dict)
            self.cds_storage.cds_fragments[cds] = sorted(suffixes)
            self.cds_storage.cds_to_delete.add(cds)

    def _update_contig_data(self, contig_id, start, stop, attributes):
        """ Updates contig data in a dictionary {contig ID: (start, stop, name, is_circular)} """

        if self.database != "ensembl":
            if ":" in contig_id:
                contig_id = contig_id.split(":")[0]
        else:
            if ":" in contig_id:
                contig_id = contig_id.split(":")[1]

        name = attributes.get("Name", "chromosome")
        self.contig_description[contig_id] = (start, stop, name, self.is_circular(attributes, database=self.database))

    def _is_pseudogene(self, gene_id: str):
        """ Returns True if the gene is a pseudogene, False otherwise """

        feature: str = self.contig_storage.gff.feature(gene_id)
        attributes: dict[str, str] = self.contig_storage.gff.attributes(gene_id)
        if self.database in {"refseq", "genbank"}:
            if get_prefix(gene_id) == "gene" and (feature == "pseudogene"
                                                  or "pseudogene" in attributes.get("biotype", "")
                                                  or attributes.get("pseudo", "") == "true"):
                return True
        elif self.database == "ensembl":
            if get_prefix(gene_id) == "gene" and (feature == "pseudogene"
                                                  or "pseudogene" in attributes.get("biotype", "")):
                return True
        return False

    @staticmethod
    def is_circular(attributes: dict[str, str], database="refseq") -> bool:
        """ Returns True if the 'Is_circular' key in the attributes dict is 'true', False otherwise"""
        if database in {"refseq", "genbank"}:
            circular: str = attributes.get("Is_circular", "false")
            if circular.lower() == "true":
                return True
        elif database == "ensembl":  # only Vertebrates genomes, only MT is circular
            if attributes["ID"] in {"chromosome:MT", "chromosome:MtDNA"}:
                return True
        return False

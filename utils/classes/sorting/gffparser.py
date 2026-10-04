from utils.classes.storage.gffdata import GffData
from utils.helpers.init_functions import initialize_dict, initialize_parameter


class GffParser:
    """
    Parses the gff data in the GffData object format from the corresponding gff3 line:
        - start and stop positions are converted to int
        - the last column of attributes is converted to a dict {"ID": id, "Parent": xxx ...}
        - fragmented genes (because of trans-splicing) are identified and treated separately

    Attributes
    ----------
        ignored_features: list[str], feature types to filter out
        only_contig: list[str], contig IDs to use
        database: str, refseq (default), genbank or ensembl
    """
    def __init__(self, global_parameters):
        """:param: parameters: GlobalParameters object"""
        if global_parameters.path_ignored_features:
            with open(global_parameters.path_ignored_features, "r") as f:
                self.ignored_features = f.readlines()
        else:
            self.ignored_features = set()
        self.only_contig = initialize_parameter(global_parameters.only_contig, [])
        self.database = global_parameters.database
        self.start_cutoff = global_parameters.start
        self.end_cutoff = global_parameters.end

    def import_annotation(self, path_gff):
        """Wrapper for importing genomic annotation from a gff file. Returns a tuple of dicts:
        contig_data: dict {contig ID: list(GffData objects)},
        contig_genepart: dict of dicts {contig ID: {gene: list geneparts}} """

        print("Loading genome annotation file")
        contig_data, contig_genepart = self._read_gff(path_gff)

        for contig in contig_genepart:
            for gene in contig_genepart[contig]:
                # order lists of geneparts, if any, by suffix
                contig_genepart[contig][gene].sort()

        return contig_data, contig_genepart

    def _read_gff(self, path_gff):
        """ Reads each line of the GFF file to import genomic feature data ,"""

        contig_genepart = {}   # dict {contig: {original gene ID: list(geneparts)}}
        contig_data = {}   # dict {contig: [GffData for each Gff line]}

        with open(path_gff, "r") as gff:
            # parse non-commented lines to populate gff_dict

            for line in gff:
                if line.rstrip() == "##FASTA":
                    break

                elif not line.startswith("#"):
                    # GffData objects are stored in a list for each contig
                    data, contig_genepart = self._parse_genomic_gff_line(line, contig_genepart)
                    if data.id:
                        contig, _, start, end = data.positions[0]
                        if contig in self.only_contig:
                            if data.feature in {"region", "chromosome"}:  # no filtering on positions
                                contig_data = initialize_dict(contig_data, contig, [])
                                contig_data[contig].append(data)
                            else:
                                if self.end_cutoff:
                                    if start >= self.start_cutoff and end <= self.end_cutoff:
                                        contig_data = initialize_dict(contig_data, contig, [])
                                        contig_data[contig].append(data)
                                else:
                                    if start >= self.start_cutoff:
                                        contig_data = initialize_dict(contig_data, contig, [])
                                        contig_data[contig].append(data)

        return contig_data, contig_genepart

    def _parse_genomic_gff_line(self, line, contig_genepart):
        """ Updates the gff_dict for each new line of the gff input with an identifier.
        Returns the feature ID, the corresponding GffData object and an updated contig_genepart dict """

        data = GffData(line, self.database)

        if not (data.feature in self.ignored_features) and data.id:
            if "genepart-" in data.id:
                # build the contig_genepart dict {contig: {original gene ID: list(geneparts)}} for trans-spliced ORFs
                original_id = data.attributes["original name"]
                contig, _, _, _ = data.positions[0]
                if contig in self.only_contig:
                    contig_genepart = initialize_dict(contig_genepart, key=contig, value={})
                    contig_genepart[contig] = initialize_dict(contig_genepart[contig], key=original_id, value=[])
                    contig_genepart[contig][original_id].append(data.id)
            return data, contig_genepart

        return (GffData(("", ("", "", [], "", "", {})), database=self.database, is_tuple=True),
                contig_genepart)  # tuple (empty data.id, empty GffData object, contig_genepart dict)

from utils.classes.storage.gffdata import GffData
from utils.classes.storage.ogtree import OgTree
from utils.classes.sorting.cdscasehandler import CDSCaseHandler
from utils.helpers.id_functions import get_prefix


class CdsProcessor:
    """CDS-specific object that manages the renaming and genomic positions of the CDS to allow fragmentation,
    and reroute CDS ontogenetic edges from rna to exon. Relies on the CdsCaseHandler object to identify particular cases
    of fragmentation/overlapping positions"""

    def __init__(self, cds, cds_data):
        self.cds = cds
        self.cds_data = cds_data
        self.exon_by_coords = {}
        self.parent = ""
        self.suffixes = []
        if self.cds_data:
            _, _, self.list_cds_positions, _, self.list_cds_phases, self.attributes = cds_data
        if "Parent" in self.attributes:
            self.parent = self.attributes.get("Parent", "")
        if self.parent:  # normalizes the transcript prefix
            self.parent.replace("transcript:", "rna-")

    def process_cds(self, contig_storage, child_parent_dict):
        """Main method controlling the flow of corrections"""

        if not self.cds or not self.parent:
            return child_parent_dict, self.suffixes

        # Map each exon (start, end) position, if any, to its ID
        self._build_exon_lookup(contig_storage, child_parent_dict)

        # Modify CDS name and attributes
        modified_cds = CDSCaseHandler(self.cds, self.cds_data, self.exon_by_coords)
        modified_cds.identify_case()
        self.suffixes = modified_cds.suffixes
        fragments = modified_cds.fragments
        data_tuples = modified_cds.data_tuples
        for i, new_id in enumerate(fragments):
            data_tuple = data_tuples[i]
            self._register_cds_part(contig_storage, child_parent_dict, new_id, data_tuple, replace=False)
        return child_parent_dict, self.suffixes

    def _build_exon_lookup(self, contig_storage, child_parent_dict):
        """Returns a dict mapping (start, end) to exon_id for all exons sharing the given parent"""
        sibling_exons = [child for child, parent in child_parent_dict.items()
                         if (parent == self.parent and get_prefix(child) == "exon")]
        for exon in sibling_exons:
            [(_, _, start, end)] = contig_storage.gff.list_positions(exon)
            self.exon_by_coords[(start, end)] = exon

    @staticmethod
    def _register_cds_part(contig_storage, child_parent_dict, new_id, data_tuple, replace=False):
        """Registers a new CDS part in contig_storage and updates child_parent_dict."""
        database = contig_storage.database
        contig_storage.gff.update(GffData((new_id, data_tuple), database=contig_storage.database, is_tuple=True),
                                  replace=replace)
        attributes = data_tuple[-1]
        child_parent_dict = OgTree.update_child_parent_dict(child_parent_dict, new_id, attributes, database=database)
        return child_parent_dict

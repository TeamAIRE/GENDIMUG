from utils.classes.storage.gff import Gff
from utils.classes.storage.intervaldata import IntervalData
from utils.classes.storage.featuretypes import FeatureTypes
from utils.classes.storage.ogtree import OgTree
from utils.classes.storage.contigsequence import ContigSequence


class ContigStorage:
    """
    Contig-specific storage object to centralize the passage of storage dicts and objects to various functions
    Attributes
    ----------
        contig: str, current contig ID
        start_cutoff: int, start position cutoff
        end_cutoff: int, end position cutoff
        database: str, database name for GFF3 format
        gff: Gff object
        positions: IntervalData object
        feature_types: FeatureTypes object
        og_tree: OgTree object
        genome: ContigSequence object
        num_index: dict {feature_id: int}
        edgetype_index: dict {edgetype: int}
        edgeclass_index: dict {edgeclass: int}
    """
    def __init__(self):
        self.contig = ""
        self.start_cutoff = 0
        self.end_cutoff = 0
        self.database = "refseq"
        self.gff = Gff()
        self.positions = IntervalData()
        self.feature_types = FeatureTypes()
        self.og_tree = OgTree()
        self.genome = ContigSequence()
        self.num_index = {}
        self.edgetype_index = {}
        self.edgeclass_index = {}
        self.edge_db = None
        self.dedup_edge_db = None

    def is_contig(self, contig_id):
        """ Assigns the current contig ID to attributes"""
        self.contig = contig_id
        self.positions.is_contig(self.contig)

    def backbone_ids(self):
        """ Returns the set of ids in parent_tree nodes, genic and intergenic sets """
        all_genes = self.feature_types.get_ids("genic")
        all_intergenic = self.feature_types.get_ids("intergenic")
        all_nested_genes = self.feature_types.get_ids("nested_genes")
        return (all_genes - all_nested_genes) | all_intergenic

    def update_positions(self, feature_id):
        """ Updates the self.positions object with the provided feature id, and its list of positions """
        list_positions = self.gff.list_positions(feature_id)
        if list_positions:
            self.positions.update(feature_id, list_positions)

    def update_positions_for_ids(self, iter_ids):
        """ Batch updates the self.backbone_position_dict dict with features from an iterable of ids """
        set_ids = {feature_id for feature_id in iter_ids if feature_id in self.gff.ids()}
        for feature_id in set_ids:
            self.update_positions(feature_id)

    def node_numeric_index(self):
        """ Converts each feature index to a numeric value to decrease the final edge file size """
        index = 0
        for feature_id in self.gff.ids():
            self.num_index[feature_id] = index
            index += 1

from utils.classes.mappers.mapperbase import Mapper
from utils.classes.sorting.nestingsorter import NestingSorter
from utils.helpers.interval_functions import matches_generator
from utils.helpers.update_functions import update_dict_of_sets


class GenicMapper(Mapper):

    def __init__(self, global_parameters, contig_storage, gene_genepart):
        super().__init__(global_parameters=global_parameters, contig_storage=contig_storage)
        self.gene_genepart = gene_genepart

    def map(self):
        """ Infers genic features (genes + pseudogenes) """
        self._genic_backbone()
        return self.contig_storage

    def _genic_backbone(self):
        """ Builds the set of genic features, including fragmented trans-spliced ORFs; identifies nested genes"""
        # build the genic set of ids: parent genes, geneparts, pseudogenes, childless genes
        self._infer_genic_backbone()

        # detect nested genes:
        self._infer_nested_genes()

        # update the key "additional_features" in feature_types, to store all non-genic/children ids
        og_tree_nodes = self.contig_storage.og_tree.all_nodes()
        genic_and_children = self.contig_storage.feature_types.get_ids("genic") | og_tree_nodes
        set_additional = self.contig_storage.gff.ids() - genic_and_children
        self.contig_storage.feature_types.update("additional_features", set_additional)

    def _infer_genic_backbone(self):
        """ Infers the genic backbone of the network by implementing corrections on the GFF data for specific cases:
        - modify_ts_rna_positions: replaces the original (sometimes wrong) genomic positions of trans-spliced RNAs
        by the list of positions from its parent geneparts (allowing changes of contig and strand between parts)
        - include_pseudogenes_as_genes: includes the gene-like orphan features in the network
          (in the "genic" feature set)
        - find_childless_genes: identifies the gene ids that are not linked to any child and returns a set of these ids
          (also added in the "genic" feature set) """
        self.contig_storage.feature_types.initialize("genic")

        # genic categories
        parent_genes = set(gene for gene in self.contig_storage.og_tree.get_genes() if gene not in self.gene_genepart)
        trans_spliced_geneparts = set(genepart for gene in self.gene_genepart
                                      for genepart in self.contig_storage.og_tree.children(gene, "genepart"))
        pseudogenes = self.contig_storage.feature_types.get_ids("pseudogenes")
        childless_genes = self.contig_storage.feature_types.get_ids("childless_genes")

        self.contig_storage.feature_types.update("genic", set.union(*[parent_genes,
                                                                      trans_spliced_geneparts,
                                                                      pseudogenes,
                                                                      childless_genes]))

        self.contig_storage.update_positions_for_ids(self.contig_storage.feature_types.get_ids("genic"))

    def _infer_nested_genes(self):
        """ Detects all genes nested into other genes; the nesting relationships are stored in a temporary dict
        {nested gene: set(nesting genes)} """
        nested_dict = {}
        for strand, interval_dict in self.contig_storage.positions.intervals().items():
            coords, nesting_dict = NestingSorter(interval_dict).remove_nested()
            for nested, _, nesting, _ in matches_generator(nesting_dict, interval_dict, interval_dict):
                nested_dict = update_dict_of_sets(nested_dict, nested, nesting)
        self.contig_storage.feature_types.update("nested_genes", set(nested_dict))
        self.contig_storage.feature_types.update("non_nested",
                                                 self.contig_storage.feature_types.get_ids("genic") - set(nested_dict))

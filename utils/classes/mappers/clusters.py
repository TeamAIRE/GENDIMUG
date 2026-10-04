from utils.classes.mappers.mapperbase import Mapper
from utils.classes.sorting.intervalmatcher import IntervalMatcher
from utils.helpers.init_functions import initialize_dict
from utils.helpers.update_functions import update_dict_of_sets
from utils.helpers.interval_functions import (convert_to_interval_dict,
                                              matches_generator)


class ClusterMapper(Mapper):

    def __init__(self, global_parameters, contig_storage):
        super().__init__(global_parameters=global_parameters, contig_storage=contig_storage)
        self.g_ig_features = {}

    def map(self):
        """ Maps features in clusters relative to backbone features. Returns the g_ig_features dict
        {strand: {genic/intergenic_id: set(feature_ids)}} """
        # dict {strand: set((start, stop, feature_id))}for nested genes and their children features
        nested_genes_dict = self._get_nested_genes_dict()
        # set of backbone feature ids
        backbone = self.contig_storage.backbone_ids()
        for strand in ["+", "-"]:
            nested_gene_interval_dict = convert_to_interval_dict(nested_genes_dict.get(strand, {}))
            backbone_interval_dict = self.contig_storage.positions.get_dict(strand, keep_ids=backbone)

            filtered_ids = set()

            if not self.genic_backbone and self.keep_children:
                filtered_ids = backbone

            elif self.genic_backbone or not self.keep_children:
                parent_genes = self.contig_storage.og_tree.get_genes()
                all_children = set()
                for gene in parent_genes:
                    children = self.contig_storage.og_tree.children(gene)
                    all_children.update(children)
                self.contig_storage.feature_types.update("children", all_children)
                filtered_ids = backbone | all_children

            additional_interval_dict = self.contig_storage.positions.get_dict(strand, filter_ids=filtered_ids)

            if backbone_interval_dict:

                self.g_ig_features = initialize_dict(self.g_ig_features, strand, {})

                # add nested genes features to the additional_interval_dict,
                # so they will be treated as additional features
                if nested_gene_interval_dict:
                    for key_coord, set_nested_gene_ids in nested_gene_interval_dict.items():
                        additional_interval_dict = update_dict_of_sets(additional_interval_dict,
                                                                       key=key_coord,
                                                                       value=set_nested_gene_ids)

                # matching_dict: dict {tuple test position: set of ref tuples with matching positions}
                matching_dict = IntervalMatcher(backbone_interval_dict, additional_interval_dict).match()

                for id_backbone, _, id_additional, _ in matches_generator(matching_dict,
                                                                          backbone_interval_dict,
                                                                          additional_interval_dict):

                    self.g_ig_features[strand] = update_dict_of_sets(self.g_ig_features[strand],
                                                                     key=id_backbone,
                                                                     value=id_additional)
        return self.g_ig_features

    def _get_nested_genes_dict(self):
        """ Builds a strand_dict {strand: {(start, stop): set(feature_ids)} for nested genes and their children """
        set_nested_genes = self.contig_storage.feature_types.get_ids("nested_genes")
        nested_genes_dict = {}
        for gene in set_nested_genes:
            children_set = self.contig_storage.og_tree.children(gene)
            list_positions = self.contig_storage.gff.list_positions(gene)
            nested_genes_dict = self._update_strand_dict_of_sets(nested_genes_dict,
                                                                 gene,
                                                                 list_positions)
            for child in children_set:
                list_positions = self.contig_storage.gff.list_positions(child)
                if list_positions:
                    nested_genes_dict = self._update_strand_dict_of_sets(nested_genes_dict,
                                                                         child,
                                                                         list_positions)
        return nested_genes_dict

    def _update_strand_dict_of_sets(self, strand_dict, feature_id, list_positions):
        """ Updates a strand_dict {strand: set(tuple(start, stop, feature_id))} with a feature id,
        keyed by a position tuple """
        if list_positions:
            list_positions = [(contig, strand, start, stop)
                              for (contig, strand, start, stop) in list_positions if all([contig, strand, start, stop])]
        for (contig, strand, start, stop) in list_positions:
            if contig == self.contig:   # would work even if the feature was fragmented on different contigs
                strand_dict = initialize_dict(strand_dict, strand, {})
                strand_dict = update_dict_of_sets(strand_dict[strand], key=(start, stop), value=feature_id)
        return strand_dict

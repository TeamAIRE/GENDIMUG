from itertools import combinations
from networkx import parse_edgelist, connected_components
from utils.classes.edges.edgebase import EdgeParameters
from utils.classes.storage.gffdata import GffData


class IdenticalEdges(EdgeParameters):

    def __init__(self,
                 global_parameters,
                 contig_storage,
                 set_edges=None,
                 filter_ids=None,
                 keep_ids=None,
                 filter_prefixes=None,
                 keep_prefixes=None):

        super().__init__(set_edges=set_edges,
                         filter_ids=filter_ids,
                         keep_ids=keep_ids,
                         filter_prefixes=filter_prefixes,
                         keep_prefixes=keep_prefixes)

        self.contig_storage = contig_storage

        # additional edge storage
        self.redundant_edges = set()
        self.kept_identical = set()

        # boolean parameters
        self.keep_children = global_parameters.keep_children
        self.keep_redundant = global_parameters.keep_redundant

    def get_edges(self):
        """ Main entry point to detect identical features, and to deduplicate the network. Returns an updated
        ContigStorage object """

        # generate edges between identical features
        print("\n + cis-identical edges...")
        filtered_ids = self.filter_ids
        if not self.keep_children:
            filtered_ids.update(self.contig_storage.feature_types.get_ids("children"))
        self._get_identical_edges(self.contig_storage.gff.get_identical_features(
            filtered_ids=filtered_ids))
        self._get_reverse_complementary_edges(self.contig_storage.gff.get_reverse_complementary_features(
            filtered_ids=filtered_ids))

        # optional deduplication of features with the same positions if they are from the same type
        if not self.keep_redundant:
            # removing redundant nodes of the same type, and nested tandem repeats
            self._deduplicate()
        else:
            # keep all identical edges
            self.contig_storage.edge_db.save_edges(self.edges)

        return self.contig_storage

    def _get_identical_edges(self, identical_features):
        """ Generates all identical edges between identical features """
        if identical_features:
            for set_ids in identical_features.values():
                for edge in list(combinations(set_ids, 2)):  # each combination of two identical features, as a tuple
                    node1, node2 = sorted(edge)  # the tuple is sorted alphanumerically as a list
                    self.edges.add((node1, node2, "cis-identical"))
                    # identical edge is created with nodes in alphanumeric order

    def _get_reverse_complementary_edges(self, rc_features):
        """ Generates all reverse complementary edges """
        if rc_features:
            for set_tuples in rc_features.values():
                for edge in set_tuples:  # each combination of two identical features, as a tuple
                    node1, node2 = sorted(edge)  # the tuple is sorted alphanumerically as a list
                    self.edges.add((node1, node2, "reverse-complementary"))
                    # reverse-complementary edge is created with nodes in alphanumeric order

    def _deduplicate(self):
        """ Deduplicates identical features from nodes and edge pools """
        print("Done.")
        if self.edges:
            print("\nDeduplicating network...")
            self._filter_identical_edges()
            self.contig_storage.edge_db.save_edges(self.kept_identical)

            # identify connected components in the graph of redundant identical features
            graph = self._identical_graph()
            ccs = connected_components(graph)
            conversion = self._compress_identical_component(ccs)

            # filter edges
            all_edges = self.contig_storage.edge_db.load_edges()
            deduplicated_edges = set()
            for node1, node2, edgetype in all_edges:
                deduplicated_edges.add((conversion.get(node1, node1), conversion.get(node2, node2), edgetype))
            print("Done.")

            # save edges
            self.contig_storage.dedup_edge_db.save_edges(deduplicated_edges)
        else:
            # simply assign the edge_db path as dedup_edge_db path
            self.contig_storage.dedup_edge_db.db_path = self.contig_storage.edge_db.db_path

    def _filter_identical_edges(self):
        """ Splits the set of identical edges into two sets: 1) redundant_edges: set of identical edges of the same type
        2) kept_identical: set of other identical edges """
        for n1, n2, edgetype in self.edges:
            if edgetype == "cis-identical":
                data_n1, data_n2 = self.contig_storage.gff.data(n1), self.contig_storage.gff.data(n2)
                if (data_n1.feature == data_n2.feature and data_n1.prefix() == data_n2.prefix()
                        and data_n1.prefix() not in {"gene"}):
                    self.redundant_edges.add((n1, n2, edgetype))
                else:
                    self.kept_identical.add((n1, n2, edgetype))
            else:
                self.kept_identical.add((n1, n2, edgetype))

    def _identical_graph(self):
        """ Builds a networkx graph from the set of identical edges of the same type """
        redundant = {f"{n1},{n2}" for n1, n2, _ in self.redundant_edges}
        return parse_edgelist(redundant, data=False, delimiter=",")

    def _compress_identical_component(self, ccs):
        """ Connected components are represented by set of feature ids. Ids and attributes are concatenated
        (separated by ||) and the new id + attributes are included in the ContigStorage.Gff object. Returns a conversion
        dict {id : new concatenated id} """
        conversion = {}
        for cc in ccs:
            db, feature, list_positions, score, list_phases, fused_attributes = "", "", [], "", [], {}

            feature_ids = sorted(list(cc))
            new_id = "==".join(feature_ids)

            for feature_id in feature_ids:
                db, feature, list_positions, score, list_phases, attributes = self.contig_storage.gff.data_tuple(
                    feature_id)
                for k, v in attributes.items():
                    if k not in fused_attributes:
                        fused_attributes[k] = []
                    fused_attributes[k].append(v)

            for feature_id in feature_ids:
                conversion[feature_id] = new_id

            for k, v in fused_attributes.items():
                if len(set(v)) == 1:
                    new_v = list(set(v))[0]
                else:
                    new_v = "||".join(v)
                fused_attributes[k] = new_v

            data_tuple = (db, feature, list_positions, score, list_phases, fused_attributes)

            # update ContigStorage.Gff object
            self.contig_storage.gff.update(GffData((new_id, data_tuple),
                                                   database=self.contig_storage.database,
                                                   is_tuple=True))
            # single feature IDs and data are removed
            self.contig_storage.gff.remove_set(feature_ids)
        return conversion

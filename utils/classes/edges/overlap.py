from utils.classes.edges.edgefactory import EdgeFactory
from utils.classes.sorting.intervalmatcher import IntervalMatcher
from utils.helpers.interval_functions import is_identical, matches_sets_generator


class Overlap(EdgeFactory):
    """
    Edge factory for overlapping/nesting edges
    Attributes
    ----------
        dict1: interval dict {(start, stop): set(ids)}
        strand1: strand associated to features in dict1, str
        dict2: Optional, interval dict {(start, stop): set(ids)}
        strand2: strand associated to features in dict2, str
        edges: Optional, set of edges
        filter_ids: Optional, set of ids to filter out
        keep_ids: Optional, set of ids to keep despite filters
        filter_prefixes: Optional, set of prefixes to filter out
        keep_prefixes: Optional, set of prefixes to keep despite filters
        trans: bool, whether to compute trans (True) or cis edges (False, default)
    :return: set of edges
    """

    def get_edges(self):
        """
        Computes all nesting/overlapping edges for the provided interval dicts.
        :return: set(edges)
        """
        if not (self.dict1 and self.dict2):
            for edge in self.edges:
                if edge:
                    yield edge

        for node1, node2, edgetype in self._get_edges_by_interval_intersection():
            if self._passes_filter(node1, node2):
                yield node1, node2, edgetype

    def _get_edges_by_interval_intersection(self):
        """
        Computes overlapping/nesting edges from the intersection of dict1 and dict2 intervals.
        :return: set of (node1, node2, edgetype) tuples
        """
        matching_dict = IntervalMatcher(self.dict1, self.dict2).match()

        for feature_id_set1, (start1, stop1), feature_id_set2, (start2, stop2) in matches_sets_generator(matching_dict,
                                                                                                         self.dict1,
                                                                                                         self.dict2):
            if is_identical((start1, stop1), (start2, stop2)):
                continue  # will be linked by cis-identical edges

            indices = {"a": feature_id_set1, "b": feature_id_set2}
            i, j, edgetype = self.get_edge_and_type((start1, stop1, "a"),
                                                    (start2, stop2, "b"),
                                                    self.strand1,
                                                    trans=self.trans)
            if edgetype not in {"cis-overlap", "cis-nested", "trans-overlap", "trans-nested"}:
                continue  # skip any edgetype that is not overlapping

            for node1 in indices[i]:
                for node2 in indices[j]:
                    if node1 != node2:
                        yield node1, node2, edgetype

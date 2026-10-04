from utils.classes.edges.edgebase import EdgeParameters
from utils.helpers.init_functions import initialize_parameter
from utils.helpers.id_functions import get_prefix
from utils.helpers.interval_functions import (is_identical,
                                              is_nested,
                                              is_overlapping)


class EdgeFactory(EdgeParameters):
    """
    Base subclass to generate edges from input lists of intervals list(set of tuples (start, stop, feature_id))
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
    """
    def __init__(self,
                 dict1=None,
                 strand1=None,
                 dict2=None,
                 strand2=None,
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
        self.dict1 = initialize_parameter(dict1, {})
        self.strand1 = initialize_parameter(strand1, "+")
        self.dict2 = initialize_parameter(dict2, dict1)
        self.strand2 = initialize_parameter(strand2, strand1)
        self.trans: bool = self.strand1 and self.strand2 and self.strand1 != self.strand2 and len(
            {self.strand1, self.strand2} & {"+", "-"}) == 2

    @staticmethod
    def _default_constraint(edge):
        """
        Placeholder edge constraint function taking and edge tuple as argument; always returns True
        (except if the edge tuple is not correctly formatted: the default constraint is to be an edge).
        :param edge: tuple(node1: str, node2: str, edgetype: str)
        :return: True
        """
        n1, n2, edgetype = edge
        return n1 and n2 and edgetype

    def _passes_filter(self, id1, id2):
        """
        Returns True if the (id1, id2) pair should be included as an edge.
        :param id1: str, feature ID 1
        :param id2: str, feature ID 2
        :return: bool
        """
        prefix1 = get_prefix(id1)
        prefix2 = get_prefix(id2)

        # 'keep' filters have priority on 'remove' filters; prefix filters have priority on ID filters
        if len({prefix1, prefix2} & self.keep_prefixes) == 2:
            return True
        if len({id1, id2} & self.keep_ids) == 2:
            return True
        if {prefix1, prefix2} & self.filter_prefixes:
            return False
        if {id1, id2} & self.filter_ids:
            return False
        return True

    @staticmethod
    def get_edge_and_type(tup1: tuple[int, int, str],
                          tup2: tuple[int, int, str],
                          strand: str,
                          trans: bool = False) -> tuple[str, str, str] | None:
        """
        Returns an edge and its relationship type as (node1, node2, edgetype), or None if intervals are not valid
        (e.g. if end < start) or no relationship is detected.

        Relationship types (prefixed with "cis-" or "trans-" based on `trans`):
          - cis-identical/reverse-complementary (trans): same start and end
          - cis/trans-nested: one feature fully contained within the other (including size one features)
          - cis/trans-overlap: features partially overlap (excluding size one features)
          - cis-sequential: features are non-overlapping and ordered

        For directional edges (nested, overlap, sequential) the first returned node is the one that comes before/inside
        the other, oriented with respect to `strand`.

        :param tup1:   (start: int, end: int, id: str) of the first feature
        :param tup2:   (start: int, end: int, id: str) of the second feature
        :param strand: "+", "-"
        :param trans:  True for trans-relationships (features on different strands),
                       False (default) for cis (features on the same strand)
        :return:       (node1, node2, edge_type) or None
        """
        (start1, end1, id1), (start2, end2, id2) = tup1, tup2
        prefix = "trans" if trans else "cis"

        # --- no treatment of invalid positions ---
        if end1 < start1 or end2 < start2 or strand not in {"-", "+"}:
            return None

        # --- identical (cis only) ---
        if is_identical((start1, end1), (start2, end2)):
            node1, node2 = sorted((id1, id2))
            if trans:
                return node1, node2, "reverse-complementary"
            else:
                return node1, node2, "cis-identical"

        # --- nested ---
        if is_nested((start1, end1), (start2, end2)):
            if start1 <= start2 and end2 <= end1:  # seq2 inside seq1
                return id2, id1, f"{prefix}-nested"
            if start2 <= start1 and end1 <= end2:  # seq1 inside seq2
                return id1, id2, f"{prefix}-nested"

        # --- overlapping ---
        # for trans edges, the strand is irrelevant

        if is_overlapping((start1, end1), (start2, end2)):
            if strand == "-" and not trans:
                id1, id2 = id2, id1  # flip for minus strand
            if start1 <= start2 and end1 <= end2:
                return id1, id2, f"{prefix}-overlap"
            if start2 <= start1 and end2 <= end1:
                return id2, id1, f"{prefix}-overlap"

        # --- sequential ---
        # for trans edges, the strand is irrelevant

        if start2 >= end1 and start1 < end2:  # seq1 before seq2
            node1, node2 = (id1, id2)
            if not trans:
                node1, node2 = (id2, id1) if strand == "-" else (id1, id2)
            return node1, node2, f"{prefix}-sequential"
        if start1 >= end2 and start2 < end1:  # seq2 before seq1
            node1, node2 = (id2, id1)
            if not trans:
                node1, node2 = (id1, id2) if strand == "-" else (id2, id1)
            return node1, node2, f"{prefix}-sequential"

        return None

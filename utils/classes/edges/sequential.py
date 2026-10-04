from utils.classes.edges.edgefactory import EdgeFactory
from utils.classes.sorting.nestingsorter import NestingSorter
from utils.helpers.init_functions import initialize_parameter


class Sequential(EdgeFactory):
    """
    Edge factory for sequential edges, will only use the first contig_dict.
    Attributes
    ----------
        dict1: interval dict {(start, end): set(ids)}
        strand1: strand associated to features in dict1, str
        dict2: Optional, interval dict {(start, end): set(ids)}
        strand2: strand associated to features in dict2, str
        edges: Optional, set of edges
        filter_ids: Optional, set of ids to filter out
        keep_ids: Optional, set of ids to keep despite filters
        filter_prefixes: Optional, set of prefixes to filter out
        keep_prefixes: Optional, set of prefixes to keep despite filters
        trans: bool, whether to compute trans (True) or cis edges (False, default)
    :return: set of edges
    """

    def get_edges(self, edge_constraint=None):
        """
        Computes all sequential edges for the provided interval dicts.
        :param edge_constraint: Optional callable
        :return: set(edges)
        """
        if not self.trans and self.dict1:
            return self._get_cis_sequential_edges(
                strand=self.strand1,
                interval_dict=self.dict1,
                edge_constraint=edge_constraint
            )

        if self.trans and self.dict1 and self.dict2:
            return self._get_trans_sequential_edges(
                plus_interval_dict=self.dict1,
                minus_interval_dict=self.dict2,
            )

        return set()

    def _get_cis_sequential_edges(self, strand, interval_dict, edge_constraint=None):
        """ Computes and returns a set of cis-sequential edges for the features in the provided interval_dict.
        This is done by removing nested features from the ordered list of feature positions, then by adding edges
        for sequential features """

        edge_constraint = edge_constraint if edge_constraint else self._default_constraint
        coord_sorting = {"+": False, "-": True}

        # remove nested features
        non_nested, _ = NestingSorter(interval_dict).remove_nested()

        if not non_nested or not interval_dict:
            return set()

        non_nested.sort(reverse=coord_sorting[strand])

        for i, (start1, end1) in enumerate(non_nested[:-1]):
            start2, end2 = non_nested[i + 1]
            feat1_set = interval_dict[(start1, end1)]
            feat2_set = interval_dict[(start2, end2)]
            set_nicknames = {"a": feat1_set, "b": feat2_set}

            key1, key2, edgetype = self.get_edge_and_type((start1, end1, "a"), (start2, end2, "b"), strand)
            if edgetype != "cis-sequential":
                continue

            for node1 in set_nicknames[key1]:
                for node2 in set_nicknames[key2]:
                    if edge_constraint((node1, node2, edgetype)) and self._passes_filter(node1, node2):
                        yield node1, node2, edgetype

        return set()

    def _get_trans_sequential_edges(self, plus_interval_dict, minus_interval_dict,  edge_constraint=None):
        """ Computes and returns a set of trans-sequential edges between plus- and minus-strand interval dicts.
        It is done by first identifying disjoint intervals for features from both strands, then only add edges
         for sequential features on different strands """

        edge_constraint = edge_constraint if edge_constraint else self._default_constraint

        plus_interval_dict = initialize_parameter(value=plus_interval_dict, default={})
        minus_interval_dict = initialize_parameter(value=minus_interval_dict, default={})

        filtered_list = self._get_trans_disjoint_intervals(plus_interval_dict, minus_interval_dict)
        strand_index = self._build_strand_position_index(plus_interval_dict, minus_interval_dict)

        for end, start, strand_end, strand_start in filtered_list:
            end_ids = strand_index[strand_end]["end"][end]
            start_ids = strand_index[strand_start]["start"][start]
            for node_end in end_ids:
                for node_start in start_ids:
                    edgetype = "trans-sequential"
                    if (edge_constraint((node_end, node_start, edgetype))
                            and self._passes_filter(node_end, node_start)):
                        yield node_end, node_start, edgetype

    @staticmethod
    def _get_trans_disjoint_intervals(plus_interval_dict, minus_interval_dict):
        """ Returns a list of disjoint interval pairs (end, start, strand_end, strand_start) for convergent (+→-)
        and divergent (-→+) trans-sequential edges candidates """

        plus_starts = {(t[0], "+") for t in plus_interval_dict}
        plus_ends = {(t[1], "+") for t in plus_interval_dict}
        minus_starts = {(t[0], "-") for t in minus_interval_dict}
        minus_ends = {(t[1], "-") for t in minus_interval_dict}

        list_types = {"divergent":  (sorted(plus_starts | minus_ends),   "-", "+"),
                      "convergent": (sorted(plus_ends | minus_starts),  "+", "-")}

        filtered_list = []
        for fused_list, first_strand, second_strand in list_types.values():
            for i in range(1, len(fused_list)):
                position1, strand1 = fused_list[i - 1]
                position2, strand2 = fused_list[i]
                if position1 < position2 and strand1 == first_strand and strand2 == second_strand:
                    filtered_list.append((position1, position2, strand1, strand2))

        return filtered_list

    @staticmethod
    def _build_strand_position_index(plus_interval_dict, minus_interval_dict):
        """ Builds a dict {strand: {start|end: {position: set(ids)}}} used for resolution of trans-sequential
        edge endpoints """

        index = {}
        index_positions = {0: "start", 1: "end"}

        for strand, interval_dict in {"+": plus_interval_dict, "-": minus_interval_dict}.items():
            index.setdefault(strand, {"start": {}, "end": {}})

            for (start, end), set_ids in interval_dict.items():
                for i, position in enumerate((start, end)):
                    key = index_positions[i]
                    index[strand][key].setdefault(position, set()).update(set_ids)

        return index

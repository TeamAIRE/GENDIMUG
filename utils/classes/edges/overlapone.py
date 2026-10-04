from utils.classes.edges.edgebase import EdgeParameters
from utils.helpers.update_functions import update_dict_of_sets
from utils.helpers.positions_functions import get_length_from_positions


class BaseOverlap(EdgeParameters):

    def __init__(self, contig_storage,
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

    def get_edges(self):
        """ Returns edges between features overlapping on one base (not detected by IntervalTree.overlaps
        that uses strict boundaries for intervals) """

        # all features for the current contig for cis-overlaps
        cis_candidates = self.contig_storage.gff.ids() - self.filter_ids

        # size one cis overlaps
        self._get_cis_overlaps(cis_candidates)

        # all features except intergenic for trans-overlaps
        trans_candidates = cis_candidates - (self.contig_storage.feature_types.get_ids("intergenic") | self.filter_ids)

        # size one trans overlaps
        self._get_trans_overlaps(trans_candidates)

        return self.edges

    def _get_cis_overlaps(self, cis_candidates):
        """ Returns edges between features overlapping on a single position in the same orientation (in cis),
        meaning the start position of one is the end position of the other """

        # dict {(strand, start/end): {start/end: set(feature IDs)} }
        cis_dict = self._build_cis_dicts(cis_candidates)

        for strand in ["+", "-"]:

            # identify cases when the start of a feature is the end of another feature
            coincidences = set(cis_dict[(strand, "start")]) & set(cis_dict[(strand, "end")])

            for pos in coincidences:
                end_ids = cis_dict[(strand, "end")][pos]
                start_ids = cis_dict[(strand, "start")][pos]

                # define node sets
                nodes1, nodes2 = end_ids, start_ids
                if strand == "-":
                    nodes1, nodes2 = nodes2, nodes1

                # build edges
                for n1 in nodes1:
                    length1 = get_length_from_positions(self.contig_storage.gff.list_positions(n1))
                    if length1 > 1:
                        for n2 in nodes2:
                            length2 = get_length_from_positions(self.contig_storage.gff.list_positions(n2))
                            if length2 > 1:
                                self.edges.add((n1, n2, "cis-overlap"))

    def _get_trans_overlaps(self, trans_candidates):
        """ Returns edges between features overlapping on a single position in the opposite orientation (in trans),
        meaning the start position of one is the end position of the other."""

        # dict {(strand, start/end): {start/end: set(feature IDs)} }
        cis_dict = self._build_cis_dicts(trans_candidates)

        # dispatch all trans-strand tests
        trans_data = {"convergent": (cis_dict[("+", "end")], cis_dict[("-", "start")]),
                      "divergent": (cis_dict[("+", "start")], cis_dict[("-", "end")])}

        positions = {"convergent": set(cis_dict[("+", "end")]) & set(cis_dict[("-", "start")]),
                     "divergent": set(cis_dict[("+", "start")]) & set(cis_dict[("-", "end")])}

        for trans_type, (nodes1, nodes2) in trans_data.items():
            for pos in positions[trans_type]:
                for n1 in nodes1[pos]:
                    length1 = get_length_from_positions(self.contig_storage.gff.list_positions(n1))
                    if length1 > 1:
                        for n2 in nodes2[pos]:
                            length2 = get_length_from_positions(self.contig_storage.gff.list_positions(n2))
                            if length2 > 1:
                                self.edges.add((n1, n2, "trans-overlap"))

    def _build_cis_dicts(self, candidate_ids):
        """ Builds a # dict {(strand, start/end): {start/end: set(feature IDs)} }
        - keys = tuples (strand, start / end position)
        - values = dict {start / end position: set(ids)} """

        # plus strand
        plus_dict = self.contig_storage.positions.get_dict("+", keep_ids=candidate_ids)
        plus_start_dict, plus_end_dict = self._convert_interval_dict_to_start_and_end_dicts(plus_dict)

        # minus strand
        minus_dict = self.contig_storage.positions.get_dict("-", keep_ids=candidate_ids)
        minus_start_dict, minus_end_dict = self._convert_interval_dict_to_start_and_end_dicts(minus_dict)

        cis_dict = {("+", "start"): plus_start_dict,
                    ("-", "start"): minus_start_dict,
                    ("+", "end"): plus_end_dict,
                    ("-", "end"): minus_end_dict}

        return cis_dict

    @staticmethod
    def _convert_interval_dict_to_start_and_end_dicts(interval_dict):
        """ Converts an interval_dict {(start, end): set(ids)} to two dicts {start: set(ids)} and {end: set(ids)} """
        start_dict, end_dict = {}, {}
        for (start, end), set_ids in interval_dict.items():
            start_dict = update_dict_of_sets(start_dict, key=start, value=set_ids)
            end_dict = update_dict_of_sets(end_dict, key=end, value=set_ids)
        return start_dict, end_dict

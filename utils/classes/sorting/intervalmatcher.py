from intervaltree import IntervalTree
from utils.helpers.update_functions import filter_dict_of_sets, reverse_dict_of_sets


class IntervalMatcher:
    """Computes overlapping interval matches from iter_tuples2 → iter_tuples1, handling size-one intervals
     (start == end) that intervaltree cannot process. Returns a dict mapping each interval in iter_tuples2 to its
     set of overlapping intervals from iter_tuples1 """

    def __init__(self, iter_tuples1=None, iter_tuples2=None):
        self.iter_tuples1 = iter_tuples1
        self.iter_tuples2 = iter_tuples2
        self.tree = None
        self.intervals_1 = None
        self.intervals_2 = None
        self.size_one_1 = None
        self.size_one_2 = None

    def match(self):
        """ Creates the matching dict between position tuples: {(start, end): set of tuples (matching positions)}"""

        self.intervals_1, self.size_one_1 = self._filter_intervals(self.iter_tuples1)
        self.intervals_2, self.size_one_2 = self._filter_intervals(self.iter_tuples2)

        # Standard intervals (start < end) : use interval tree for overlap detection
        matching_dict = self._match_intervals()

        # Size-one intervals in iter_tuples1 (start == end): match against all of iter_tuples2
        # Results are reversed so keys stay in iter_tuples2 space
        if self.size_one_1:
            positions_1 = [start for start, _ in self.size_one_1]
            size_one_matches = self._match_size_one_to_intervals(positions_1, self.intervals_2)
            self._merge_reversed(size_one_matches, matching_dict)

        # Size-one intervals in iter_tuples2: match against all of iter_tuples1
        if self.size_one_2:
            positions_2 = [start for start, _ in self.size_one_2]
            size_one_matches = self._match_size_one_to_intervals(positions_2, self.intervals_1)
            self._merge(size_one_matches, matching_dict)

        return filter_dict_of_sets(matching_dict, cutoff=0)

    @staticmethod
    def _merge(source: dict, target: dict) -> None:
        """ Merges source into target, union of sets for existing keys """
        for key, values in source.items():
            target.setdefault(key, set()).update(values)

    def _merge_reversed(self, source: dict, target: dict) -> None:
        """ Merges the reverse of source into target (swaps keys ↔ values) """
        self._merge(reverse_dict_of_sets(source), target)

    @staticmethod
    def _filter_intervals(iter_intervals):
        """ Removes positions of size 1 (where start == end) in an iterable of interval tuples (start, end),
        as they can not be used within an interval tree. Returns the filtered iterable of intervals and
        the iterable of size-one intervals """
        intervals = set()
        size_one_intervals = set()
        for tup in iter_intervals:
            start, end = tup[0], tup[1]
            (size_one_intervals if start == end else intervals).add((start, end))
        return intervals, size_one_intervals

    def _get_interval_matching_set(self, start, end):
        """ Computes the matches between an interval (start, end) and an interval tree; returns the set of
        matching intervals """
        return {(e.begin, e.end) for e in self.tree.overlap(start, end)}

    def _match_intervals(self):
        """ Returns the dict of sets of position intervals, as tuples (start:int, end:int), matching all intervals
        from a reference list of intervals: matching dict {(start, end): set of tuples (matching positions)} """
        self.tree = IntervalTree.from_tuples(self.intervals_1)
        matching_dict = {}
        for start, end in self.intervals_2:
            matching_dict[(start, end)] = self._get_interval_matching_set(start, end)
        return matching_dict

    @staticmethod
    def _match_size_one_to_intervals(iter_position, iter_tuples, keep_identical: bool = True):
        """ Returns the set of matching dict {(position, position): set of tuples (matching positions)}."""
        list_tuples = list(iter_tuples)
        if not list_tuples:
            return {}

        tree = IntervalTree.from_tuples((start, end + 1, (start, end)) 
                                        if start == end else (start, end, (start, end))
                                        for start, end in list_tuples)

        matching_dict_one = {}
        for pos in iter_position:
            matches = {interval.data for interval in tree.at(pos)}
            if not keep_identical:
                matches -= {(pos, pos)}
            if matches:
                matching_dict_one[(pos, pos)] = matches
        return matching_dict_one

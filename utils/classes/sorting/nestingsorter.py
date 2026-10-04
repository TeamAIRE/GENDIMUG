from intervaltree import IntervalTree
from utils.helpers.init_functions import initialize_parameter


class NestingSorter:
    """
    Sorts nested and nesting features; the remove_nested method returns a list of non-nested position intervals and a
    nesting dict {(start, end): set of nested positions (start, end)}
    Attributes
    ----------
        positions: Iterable of tuples (start: int, end: int) or tuples (start: int, end: int, feature_id:str)
    """
    def __init__(self, positions):
        self.positions = list(positions)

    def remove_nested(self, dict_nesting_nested=None):
        """ Removes all nested positions from an iterable of position tuples (start, end). Returns an updated sorted
        list of non-nested tuples and a filtered dict_nesting_nested. """
        list_intervals = initialize_parameter(self.positions, [])
        dict_nesting_nested = initialize_parameter(dict_nesting_nested, {})

        if not list_intervals:
            return list_intervals, dict_nesting_nested

        unique = list({(t[0], t[1]) for t in list_intervals})
        true_intervals = [(start, end) for start, end in unique if start != end]
        size_one = [(start, end) for start, end in unique if start == end]

        tree = IntervalTree.from_tuples(true_intervals) if true_intervals else IntervalTree()

        for start, end in true_intervals:
            # envelop() returns intervals strictly contained within (start, end), identical intervals are excluded

            children = {(interval.begin, interval.end) for interval in tree.envelop(start, end)
                        if (interval.begin, interval.end) != (start, end)}

            # size-one intervals contained in (start, end)
            children.update((ps, pe) for ps, pe in size_one if start <= ps <= end)

            if children:
                dict_nesting_nested.setdefault((start, end), set()).update(children)

        dict_nesting_nested = self._keep_direct_links(dict_nesting_nested)

        all_nested = set().union(*dict_nesting_nested.values()) if dict_nesting_nested else set()
        non_nested = sorted(t for t in list_intervals if (t[0], t[1]) not in all_nested)

        return non_nested, dict_nesting_nested

    @staticmethod
    def _keep_direct_links(dict_nesting_nested):
        """ Keeps only direct nesting -> nested links in the dict. For each nesting interval,
        removes any nested interval that appears in another nested interval's nested subset.
        Returns the filtered dict """
        filtered = {}
        for nesting, set_nested in dict_nesting_nested.items():
            # Collect everything nested under any nested interval
            set_sub_nested = set().union(*(dict_nesting_nested[nested]
                                           for nested in set_nested
                                           if nested in dict_nesting_nested))
            # Direct children are those not reachable via another child
            filtered[nesting] = set_nested - set_sub_nested
        return filtered

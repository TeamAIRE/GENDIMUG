from utils.helpers.init_functions import initialize_parameter, initialize_dict
from utils.helpers.id_functions import get_prefix


class IntervalData:
    """
    Stores positions in an interval dictionary {contig: {strand: {(start, stop): set(ids)}}}
    Attributes
    ----------
        interval_dict: dict {strand: {(start, stop): set(ids)}}
    """
    def __init__(self):
        self.interval_dict = {}
        self.contig = ""

    def is_contig(self, contig_id):
        """ Assigns the current contig ID to attributes"""
        self.contig = contig_id

    def intervals(self):
        """ Returns the interval_dict """
        return self.interval_dict

    def update(self, feature_id, list_positions):
        """ Updates the interval dict based on the positions of the feature """
        list_positions = self._check_list_positions(list_positions)
        if list_positions:
            for (contig, strand, start, stop) in list_positions:
                if contig == self.contig:
                    start, stop = tuple(sorted([start, stop]))
                    self.interval_dict = initialize_dict(self.interval_dict, strand, {})
                    self.interval_dict[strand] = initialize_dict(self.interval_dict[strand],
                                                                 (start, stop), set())
                    self.interval_dict[strand][(start, stop)].add(feature_id)

    def remove(self, feature_id, list_positions):
        """ Removes a feature from the interval dict based on the positions of the feature """
        # check the list positions is correctly formated
        list_positions = self._check_list_positions(list_positions)
        if list_positions:
            for (contig, strand, start, stop) in list_positions:
                if contig == self.contig:
                    start, stop = tuple(sorted([start, stop]))
                    if self.interval_dict[strand].get((start, stop), set()):
                        self.interval_dict[strand][(start, stop)] = (self.interval_dict[strand][(start, stop)]
                                                                     - {feature_id})
                    if not self.interval_dict[strand].get((start, stop), set()):  # test for empty set
                        self.interval_dict[strand].pop((start, stop), None)

    def remove_set(self, set_feature_ids, list_positions):
        """ Removes a set of feature ids (with the same positions) from the interval dict """
        for feature_id in set_feature_ids:
            self.remove(feature_id, list_positions)

    @staticmethod
    def _check_list_positions(list_positions):
        """ Checks that the list positions is correctly formated"""
        return [(contig, strand, start, stop) for (contig, strand, start, stop) in list_positions
                if all([contig, strand, start, stop])] if list_positions else []

    def get_dict(self, strand, filter_ids=None, keep_ids=None, filter_prefixes=None, keep_prefixes=None):
        """ Returns the interval_dict {(start, stop): set(ids)} corresponding to the provided contig and strand.
        Features can be filtered out by prefixes or by feature ids:
        - provide a set of prefixes to remove with filter_prefixes
        - provide a set of prefixes to keep exclusively with keep_prefixes
        - provide a set of ids to remove with filter_ids
        - provide a set of ids to keep exclusively with keep_ids """
        filter_prefixes = initialize_parameter(value=filter_prefixes, default=set())
        keep_prefixes = initialize_parameter(value=keep_prefixes, default=set())
        filter_ids = initialize_parameter(value=filter_ids, default=set())
        keep_ids = initialize_parameter(value=keep_ids, default=set())
        interval_dict = self.interval_dict.get(strand, {})
        filtered_dict = {}
        for (start, stop), set_ids in interval_dict.items():
            if filter_ids:
                set_ids = set_ids - filter_ids
            if keep_ids:
                set_ids = set_ids & keep_ids
            if filter_prefixes:
                set_ids = {feature_id for feature_id in set_ids if get_prefix(feature_id) not in filter_prefixes}
            if keep_prefixes:
                set_ids = {feature_id for feature_id in set_ids if get_prefix(feature_id) in keep_prefixes}
            if set_ids:
                filtered_dict[(start, stop)] = set_ids
        return filtered_dict

    def get_list(self, strand, filter_ids=None, keep_ids=None, filter_prefixes=None, keep_prefixes=None):
        """ Returns the ordered list of tuples (start, stop, feature_id) corresponding to the provided contig
        and strand.vFeatures can be filtered out by prefixes or by feature ids:
        - provide a set of prefixes to remove with filter_prefixes
        - provide a set of prefixes to keep exclusively with keep_prefixes
        - provide a set of ids to remove with filter_ids
        - provide a set of ids to keep exclusively with keep_ids """
        filter_prefixes = initialize_parameter(value=filter_prefixes, default=set())
        keep_prefixes = initialize_parameter(value=keep_prefixes, default=set())
        filter_ids = initialize_parameter(value=filter_ids, default=set())
        keep_ids = initialize_parameter(value=keep_ids, default=set())
        list_features = []
        interval_dict = self.interval_dict.get(strand, {})
        for (start, stop), set_ids in interval_dict.items():
            if filter_ids:
                set_ids = set_ids - filter_ids
            if keep_ids:
                set_ids = set_ids & keep_ids
            if filter_prefixes:
                set_ids = {feature_id for feature_id in set_ids if get_prefix(feature_id) not in filter_prefixes}
            if keep_prefixes:
                set_ids = {feature_id for feature_id in set_ids if get_prefix(feature_id) in keep_prefixes}
            for feature_id in set_ids:
                list_features.append((start, stop, feature_id))
        return sorted(list_features)

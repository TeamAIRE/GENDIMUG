from numpy import array
from utils.helpers.update_functions import update_dict_of_sets
from utils.helpers.positions_functions import get_min_max_positions


def is_match(coord1: tuple[int, int], coord2: tuple[int, int]) -> bool:
    """Returns True if coord1 and coord2 have overlapping positions, False otherwise"""
    return is_overlapping(coord1, coord2) or is_nested(coord1, coord2) or is_identical(coord1, coord2)


def is_identical(coord1: tuple[int, int], coord2: tuple[int, int]) -> bool:
    """Returns True if coord1 == coord2, False otherwise"""
    return coord1 == coord2


def is_nested(coord1: tuple[int, int], coord2: tuple[int, int]) -> bool:
    """Returns True if one interval is nested into the other (one position can be identical but not both,
     and one interval can be of size one), False otherwise"""
    (start1, end1), (start2, end2) = coord1, coord2
    return (start1 < start2 <= end2 <= end1
            or start1 <= start2 <= end2 < end1
            or start2 < start1 <= end1 <= end2
            or start2 <= start1 <= end1 < end2)


def is_overlapping(coord1: tuple[int, int], coord2: tuple[int, int]) -> bool:
    """Returns True if the intervals are overlapping or nested (one position can be identical but not both,
     size one interval not allowed), False otherwise"""
    (start1, end1), (start2, end2) = coord1, coord2
    return (start1 < start2 < end1 <= end2
            or start1 <= start2 < end1 < end2
            or start1 < start2 <= end1 < end2
            or start2 < start1 < end2 <= end1
            or start2 <= start1 < end2 < end1
            or start2 < start1 <= end2 < end1)


def has_overlapping_intervals(list_positions):
    """
    Returns True if any two intervals in the list overlap.
    :param list_positions: list of (contig, strand, start, end) tuples, where start < end
    :return: True if at least one overlap exists, False otherwise
    """
    sorted_list_positions = sorted(list_positions, key=lambda x: x[2])
    for i in range(len(sorted_list_positions) - 1):
        current_contig, current_strand, _, current_end = sorted_list_positions[i]
        next_contig, next_strand, next_start, _ = sorted_list_positions[i + 1]
        if current_contig == next_contig and current_strand == next_strand:
            if next_start <= current_end:
                return True
    return False


def has_nested_intervals(list_intervals):
    """
    Returns True if there is any nested interval in a list of intervals (excluding identical intervals),
    using numpy arrays of start and end positions.
    :param list_intervals: list of tuples (start: int, end: int, ...)
    :return: bool
    """
    if list_intervals:
        list_intervals.sort()
        list1 = list_intervals[0:-1]
        list2 = list_intervals[1:]
        # numpy arrays of start and end positions:
        s1 = array([e[0] for e in list1])
        s2 = array([e[0] for e in list2])
        e1 = array([e[1] for e in list1])
        e2 = array([e[1] for e in list2])
        return (any(((s1 < s2) & (s2 <= e2) & (e2 <= e1))) or
                any(((s1 <= s2) & (s2 <= e2) & (e2 < e1))) or
                any(((s2 <= s1) & (s1 <= e1) & (e1 < e2))) or
                any(((s2 < s1) & (s1 <= e1) & (e1 <= e2))))
    return False


def get_length_from_overlapping_positions(list_positions):
    """
    Computes the cumulated length of all segments in a list of position tuples
    :param list_positions: list of tuples (contig:str, strand:str, start:int, end:int)
    :return: int
    """
    feature_length = 0
    list_clusters = group_overlapping_intervals(list_positions)
    for cluster in list_clusters:
        if len(cluster) == 1:
            for _, _, start, end in cluster:
                feature_length += abs(end - (start - 1))
        else:  # overlapping intervals due to (ribosomal or other) slippage:
            # the feature length correspond to the length of genomic positions involved
            # (for CDS, may not be multiple of 3)
            min_start, max_end = get_min_max_positions(cluster)
            feature_length += abs(max_end - (min_start - 1))
    return feature_length


def group_overlapping_intervals(list_positions):
    """
    Returns a sorted list of lists, each sublist containing the sorted intervals with overlaps ("clusters").
    :param list_positions: list of tuples (contig: str, strand: str, start: int, end: int)
    :return: list of sublists, each containing overlapping intervals"""
    result = []
    if len(list_positions) == 1:
        return [list_positions]

    elif has_overlapping_intervals(list_positions):

        current_cluster = []
        # sort positions by start:
        list_positions = sorted(list_positions, key=lambda x: x[2])

        for i in range(len(list_positions[:-1])):
            (c1, s1, start1, end1) = list_positions[i]
            current_cluster.append((c1, s1, start1, end1))

            # all positions but the last one
            if i < len(list_positions) - 2:
                (c2, s2, start2, end2) = list_positions[i + 1]
                if not is_overlapping((start1, end1), (start2, end2)):
                    result.append(current_cluster.copy())
                    current_cluster = []

            # last position
            else:
                (c2, s2, start2, end2) = list_positions[i + 1]
                if is_overlapping((start1, end1), (start2, end2)):
                    current_cluster.append((c2, s2, start2, end2))
                    result.append(current_cluster.copy())

                else:
                    result.append(current_cluster.copy())
                    result.append([(c2, s2, start2, end2)])
    return result


def convert_to_interval_dict(interval_set):
    """
    Converts a set of interval tuples (start, stop, id) to an interval_dict {(start, stop): set(ids)},
    while filtering out intergenic regions and feature with same start and stop.
    :param interval_set: set((start, stop, id))
    :return: interval_dict {(start, stop): set(ids)}
    """
    interval_dict = {}
    for start, stop, feature_id in interval_set:
        coord = tuple(sorted([start, stop]))
        interval_dict = update_dict_of_sets(interval_dict, key=coord, value=feature_id)
    return interval_dict


def matches_generator(matching_dict, keys_interval_dict, values_interval_dict):
    """
    Generator looping through matching_dict and interval_dicts to return pairs of matching features
    (ids and start/end positions), which will be used as key-value pair in a storage dict.
    :param matching_dict: dict {(start, end): set of tuples (matching positions)}
    :param keys_interval_dict: dict {(start, end): set(key ids)}
    :param values_interval_dict: dict {(start, end): set(value ids)}
    :return: Generator of tuples of matching (key_feature, (start_k, end_k), value_feature, (start_v, end_v)
    """
    for value_coord, set_value_coord in matching_dict.items():
        for key_coord in set_value_coord:
            for key_id in keys_interval_dict.get(key_coord, set()):
                for value_id in values_interval_dict.get(value_coord, set()):
                    yield key_id, key_coord, value_id, value_coord


def matches_sets_generator(matching_dict, keys_interval_dict, values_interval_dict):
    """
    Generator looping through matching_dict and interval_dicts to return pairs of matching features,
     which will be used as key-value pair in a storage dict.
    :param matching_dict: dict {(start, stop): set of tuples (matching positions)}
    :param keys_interval_dict: dict {(start, stop): set(key ids)}
    :param values_interval_dict: dict {(start, stop): set(value ids)}
    :return: tuple of matching (key_feature1, (start_k, stop_k), value_feature2, (start_v, stop_v)
    """
    for value_coord, set_key_coords in matching_dict.items():
        set_values = values_interval_dict.get(value_coord, set())
        for key_coord in set_key_coords:
            set_keys = keys_interval_dict.get(key_coord, set())
            yield set_keys, key_coord, set_values, value_coord

from utils.helpers.init_functions import initialize_parameter, initialize_dict, initialize_contig_dict


def update_dict_of_sets(dict_of_sets, key, value):
    """
    Returns a dict of sets initialized for the provided key, with the provided value(s) added to the set at that key.
    :param dict_of_sets: dict
    :param key: key to set
    :param value: single hashable value or set of hashable values to add
    :return: updated dict {key: set(value)}
    """
    dict_of_sets = initialize_dict(dict_of_sets, key, set())
    values = value if isinstance(value, set) else {value}
    dict_of_sets[key].update(values)
    return dict_of_sets


def filter_dict_of_sets(dict_of_sets, cutoff=0):
    """
    Filters a dict of sets { a: set(b) } based on the length of set(b)
    :param dict_of_sets: dict { a: set(b) }
    :param cutoff: size limit
    :return: filtered dict_of_sets: dict { a: set(b) } with len(set(b)) > cutoff
    """
    return {k: set_v for k, set_v in dict_of_sets.items() if len(set_v) > cutoff}


def reverse_dict_of_sets(dict_of_sets):
    """
    Reverses a dict of sets { a: set(b) } to a dict {b: set(a)}
    :param dict_of_sets: dict { a: set(b) }
    :return: dict {b: set(a)}
    """
    rev_dict = {}
    for k, set_v in dict_of_sets.items():
        for v in set_v:
            rev_dict = update_dict_of_sets(rev_dict, key=v, value=k)
    return rev_dict


def update_interval_dict_from_positions(feature_id, list_positions, interval_dict=None):
    """
    Updates the interval_dict {contig: {strand: {(start, stop): set of ids} } }
    with the provided id, contig, strand, start and stop positions.
    :param feature_id: str, feature id
    :param list_positions: list of tuples (contig:str, strand:str, start:int, stop:int)
    :param interval_dict: optional, default None, dict {contig: {strand: {(start, stop): set of ids} } }
    :return: updated interval_dict
    """
    interval_dict = initialize_parameter(interval_dict, {})
    for (contig, strand, start, stop) in list_positions:
        interval_dict = initialize_contig_dict(interval_dict, contig, strand, {})
        interval_dict[contig][strand] = initialize_dict(interval_dict[contig][strand], (start, stop), set())
        interval_dict[contig][strand][(start, stop)].add(feature_id)
    return interval_dict


def compatible_contig_dicts(dict1, dict2, contig, strand):
    """
    Tests whether contig and strand are keys in dict1 and dict2.
    :param dict1: dict {contig: {strand: container}}
    :param dict2: dict {contig: {strand: container}}
    :param contig: str, contig id
    :param strand: str, + or -
    :return: bool
    """
    return contig in dict1 and contig in dict2 and strand in dict1[contig] and strand in dict2[contig]

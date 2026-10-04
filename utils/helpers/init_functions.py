def initialize_parameter(value=None, default=None):
    """
    Returns value if not None, otherwise default.
    :param value: any object or None
    :param default: value or variable
    :return: updated parameter
    """
    return value if value is not None else default


def initialize_dict(data_dict, key, value):
    """
    Returns a dict of container (dict, list or set) {key: container} initialized for the provided key,
    or the same dict if the key already exists.
    {key: container}.
    :param data_dict: dict
    :param key: key to set
    :param value: value at the provided key
    :return: dict {key: container}
    """
    data_dict.setdefault(key, value)
    return data_dict


def initialize_contig_dict(contig_dict, contig, strand, container):
    """
    Returns a dict of dicts initialized for the provided contig and strand keys
    {contig: {strand: container}}. The container can be a dict, list or set.
    :param contig_dict: dict
    :param contig: str, contig id
    :param strand: str, '+' or '-'
    :param container: data structure dict(), set(), list() etc...
    :return: dict {contig: {strand: container}}
    """
    contig_dict = initialize_dict(contig_dict, key=contig, value={})
    contig_dict[contig] = initialize_dict(contig_dict[contig], key=strand, value=container)
    return contig_dict

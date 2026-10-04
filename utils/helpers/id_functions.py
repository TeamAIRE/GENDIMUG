def get_prefix(feature_id):
    """
    Returns the prefix of the feature id (in GFF3 format: lowercase prefix + '-')
    :param feature_id: str, feature id
    :return: str, prefix of the feature id or ""
    """
    possible_separators = {":", "-"}
    for c in feature_id:
        if c in possible_separators:
            return feature_id.split(c)[0].lower()
    return None


def get_suffix(feature_id):
    """
    Returns the fragment numeric suffix (as int), if any, of the feature id, or None.
    :param feature_id: str, feature id
    :return: int, fragment suffix of the feature id, or None
    """
    if "-" in feature_id:
        suffix = feature_id.split("-")[-1]
        if suffix.isnumeric():
            return int(suffix)
    return None


def reorder_parts(parts):
    """
    Returns a list of parts ordered according to their last digit
    :param parts: iterable of numbered parts
    :return: ordered list of parts
    """
    filtered_parts = [e for e in parts if get_suffix(e)]
    ordered_parts_dict = {get_suffix(part): part for part in filtered_parts}
    return [ordered_parts_dict[k] for k in sorted(ordered_parts_dict)]


def reformat_ensembl_id(feature_id):
    """
    Returns the reformatted feature ID modified from the original Ensembl feature ID.
    :param feature_id: str, feature ID
    :return: str, reformatted feature ID
    """
    reformated_id = feature_id.replace("transcript:", "rna-").replace(":", "-")
    prefix = reformated_id.split("-")[0].lower()
    return prefix + "-" + "".join(reformated_id.split("-")[1:])

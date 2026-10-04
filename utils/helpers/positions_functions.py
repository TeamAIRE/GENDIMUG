def compatible_positions(list_positions1, list_positions2, same_strand=True):
    """
    Tests whether contigs and strands are the same for both lists
    :param list_positions1: list(tuple(contig, strand, start, end))
    :param list_positions2: list(tuple(contig, strand, start, end))
    :param same_strand: bool, whether the features are on the same strand
    :return: bool
    """
    contig1 = _consensus_position(list_positions1, 0)
    contig2 = _consensus_position(list_positions2, 0)
    if same_strand:
        strand1 = _consensus_position(list_positions1, 1)
        strand2 = _consensus_position(list_positions2, 1)
        if contig1 and contig2 and strand1 and strand2:
            return contig1 == contig2 and strand1 == strand2
    else:
        if contig1 and contig2:
            return contig1 == contig2
    return False


def _consensus_position(list_positions, index):
    """
    Returns the consensus unique value for a list of positions at the provided tuple index.
    (e.g. 0 is contig, 1 is strand, 2 is start and 3 is end), or an empty string if there is no consensus.
    :param list_positions: list(tuple(contig, strand, start, end))
    :param index: int, 0 <= index <= 3
    :return: str (contig, strand) or int (start, end)
    """
    consensus = set()
    for span in list_positions:
        consensus.add(span[index])
    if len(consensus) == 1:
        return list(consensus)[0]
    return ""


def get_length_from_positions(list_positions):
    """
    Computes the cumulated length of all segments in a list of position tuples
    :param list_positions: list of tuples (contig:str, strand:str, start:int, end:int)
    :return: int
    """
    feature_length = 0
    for _, _, start, end in list_positions:
        feature_length += abs(end - (start - 1))
    return feature_length


def get_min_max_positions(list_positions):
    """
    Returns the minimal position and the maximal position from a list of position tuples.
    :param list_positions: list(tuple(contig, strand, start, end))
    :return: min position, max position
    """
    starts = list(start for _, _, start, _ in list_positions)
    ends = list(end for _, _, _, end in list_positions)
    return min(starts), max(ends)


def convert_string_positions_to_list(str_positions):
    """
    Converts a string encoding positions to a list of tuples (contig:str, strand:str, start:int, end:int)
    :param str_positions: str
    :return: list of tuples
    """
    if str_positions:
        positions = [tuple(e.split(";")) for e in str_positions.split("||")]
        return [(contig, strand, int(start), int(end)) for contig, strand, start, end in positions]
    else:
        return ""

def attributes_to_dict(attributes):
    """
    Converts the last column data from ggf3 file to a dict (constructed by splitting the line with ";" and "="
    as separators).
    :param attributes: str, data from the last gff3 column
    :return: attributes dict {key: value} derived from the data in the last column of the gff line
    """
    return {e.rstrip().split("=")[0]: e.rstrip().split("=")[1] for e in attributes.split(";")
            if ("=" in e and len(e.split("=")) >= 2)}


def format_gff_line(line):
    """
    Formats a gff line by converting start and end from str to int and the attributes data to a dict.
    :param line: line read from the gff3 input file
    :return: formatted line
    """
    contig, db, feature, start, end, score, strand, phase, attributes = line.rstrip().split("\t")
    return contig, db, feature, int(start), int(end), score, strand, phase, attributes_to_dict(attributes)

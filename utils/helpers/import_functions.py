from utils.helpers.init_functions import initialize_dict


def import_genome_sequence(path_fasta):
    """
    Wrapper for pickling fasta files.
    :param path_fasta: str, path to the input genomic fasta file
    :return: fasta dict {contig_id: sequence}, path_output_dir, set of contig names
    """
    print("Loading genomic fasta file")
    fasta_dict = _import_fasta(path_fasta)  # whole genome
    contig_names = set(fasta_dict)
    return fasta_dict, contig_names


def _import_fasta(path_fasta):
    """
    Converts a genomic fasta file to a dict {contig_id: sequence}.
    :param path_fasta: path to the input genomic fasta file
    :return: fasta dict {contig ID: sequence}
    """
    tag = ""
    fasta_dict = {}

    with open(path_fasta, "r") as f:

        for line in f:

            if line.startswith(">"):
                tag = line.rstrip().split(" ")[0].replace(">", "")
                fasta_dict = initialize_dict(fasta_dict, tag, [])

            else:
                fasta_dict[tag].append(line.rstrip())

    return {tag: "".join(seq_list) for tag, seq_list in fasta_dict.items()}


def read_file_to_dict(path, sep=",", colindexk=0, colindexv=1, skipfirst=False):
    """
    Read a csv file and creates a dictionary from two columns.
    :param path: path to the csv file
    :param sep: separator, comma by default
    :param colindexk: index of the column for dictionary keys
    :param colindexv: index of the column for dictionary values
    :param skipfirst: boolean, True to skip first line (column names), False by default
    :return: dictionary
    """
    with open(path, "r") as f:
        lines = f.readlines()
        list_lines = [i.rstrip() for i in lines]
    attr_dict = {}
    if skipfirst:
        lines = list_lines[1:]
    else:
        lines = list_lines
    for line in lines:
        attr_dict[line.split(sep)[colindexk]] = line.split(sep)[colindexv]
    return attr_dict

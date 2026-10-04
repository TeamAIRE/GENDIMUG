from os.path import join
from utils.helpers.pickle_functions import load_pickle, save_pickle
from utils.helpers.os_sys_functions import check_output_dir


def save_data(contig_storage, path_output_dir):
    """
    Serializes ContigStorage object to pickle file
    :param contig_storage: ContigStorage object
    :param path_output_dir: str, path to the output directory
    :return: None
    """
    print("(Pickling data file)")
    path_temp_output_dir = join(path_output_dir, "temp_data_storage")
    check_output_dir(path_temp_output_dir)
    filename = f"{contig_storage.contig}_contig_storage.pickle"
    save_pickle(contig_storage, join(path_temp_output_dir, filename))


def path_db_edges(contig, path_output_dir, dir_name, filename):
    """
    Returns the path to the sqlite3 database for edges
    :param contig: str, contig ID
    :param path_output_dir: str, path to the output directory
    :param dir_name: str, name of the directory
    :param filename: str, name of the file
    :return:, str, path to the database file
    """
    path_temp_output_dir = join(path_output_dir, dir_name)
    check_output_dir(path_temp_output_dir)
    filename = f"{contig}_{filename}.db"
    return join(path_temp_output_dir, filename)


def load_data(contig, path_output_dir):
    """
    Deserializes ContigStorage object from pickle file
    :param contig: str, contig ID
    :param path_output_dir: str, path to the output directory
    :return: ContigStorage object
    """
    path_temp_output_dir = join(path_output_dir, "temp_data_storage")
    check_output_dir(path_temp_output_dir)
    filename = f"{contig}_contig_storage.pickle"
    return load_pickle(join(path_temp_output_dir, filename))

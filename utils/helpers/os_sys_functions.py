from sys import exit
from os import makedirs, listdir
from os.path import abspath, exists, join, isfile, isdir


def exit_message():
    exit(f"\nEXIT: no features indicated.")


def check_input_path(path_file: str) -> None:
    """ Check if input file exists
    :param path_file: str, path to the input file
    :return: None
    """
    if path_file:
        if not exists(path_file):
            exit(f"\nEXIT: Unable to locate {path_file}. Please provide a correct path to the input file/folder")


def check_output_dir(path_output_dir: str) -> None:
    """
    Creates output directory if parent directory exists
    :param path_output_dir: str, path to the output directory
    :return: None
    """
    path_output_dir: str = abspath(path_output_dir)
    # Create output directory if parent directory exists
    if path_output_dir and not exists(path_output_dir):
        makedirs(path_output_dir)


def get_list_of_files(dir_path: str) -> list[str]:
    """
    Get all filenames in a given directory
    :param dir_path: path to the directory
    :return: list of file paths
    """
    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]


def get_list_of_dirs(dir_path):
    """
    Get all filenames in a given directory
    :param dir_path: path to the directory
    :return: list of file paths
    """
    return [name for name in listdir(dir_path) if isdir(join(dir_path, name))]

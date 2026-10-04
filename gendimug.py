#!/usr/bin/env python3
from utils.classes.storage.globalparameters import GlobalParameters
from utils.helpers.arg_parser_functions import main_arg_parser
from utils.helpers.time_functions import get_current_time
from wrappers.gendimug_main_wrappers import build_network


def main():
    start_time = get_current_time()
    global_parameters = GlobalParameters()
    global_parameters.import_parameters(main_arg_parser())
    build_network(global_parameters, start_time)


if __name__ == "__main__":
    main()

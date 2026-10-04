from utils.helpers.init_functions import initialize_parameter


class EdgeParameters:
    """ Base class for edge factories"""
    def __init__(self,
                 set_edges=None,
                 filter_ids=None,
                 keep_ids=None,
                 filter_prefixes=None,
                 keep_prefixes=None):
        self.edges = initialize_parameter(set_edges, set())
        self.filter_ids = initialize_parameter(filter_ids, set())
        self.keep_ids = initialize_parameter(keep_ids, set())
        self.filter_prefixes = initialize_parameter(filter_prefixes, set())
        self.keep_prefixes = initialize_parameter(keep_prefixes, set())


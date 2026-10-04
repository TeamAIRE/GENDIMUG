from utils.helpers.init_functions import initialize_dict, initialize_parameter


class FeatureTypes:
    """ Dict of sets storing feature IDs in sets for classification into various feature types """

    def __init__(self):
        self.feature_types = {}

    def initialize(self, *args):
        """ Creates empty sets at the provided keys in feature_sets dict """
        keys = initialize_parameter([*args], [])
        for k in keys:
            self.feature_types = initialize_dict(self.feature_types, k, set())

    def update(self, key, set_ids=None):
        """ Updates a set at the provided key with the ids in the set_ids """
        self.initialize(key)
        self.feature_types[key].update(set_ids)

    def remove(self, key, set_ids=None):
        """ Removes the ids in the set_ids from set at the provided key """
        if key in self.feature_types:
            self.initialize(key)
            self.feature_types[key] = self.feature_types[key] - set_ids

    def get_ids(self, *args, add_ids=None):
        """ Returns a set of feature ids corresponding to the provided keys.
         An optional set of ids to add can be provided """
        keys = initialize_parameter([*args], [])
        set_ids = set()
        for k in keys:
            set_ids.update(set(self.feature_types.get(k, set())))
        return set_ids | initialize_parameter(add_ids, set())

    def get_keys(self):
        """ Returns a set of all keys (feature types) (used for debugging) """
        return set(self.feature_types.keys())

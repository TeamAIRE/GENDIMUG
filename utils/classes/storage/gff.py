from itertools import combinations
from utils.helpers.init_functions import initialize_parameter, initialize_dict
from utils.helpers.update_functions import filter_dict_of_sets, update_dict_of_sets


class Gff:
    """
    Stores, updates and retrieves the gff genomic annotation data (stored in a GffData object)
    in a gff_dict keyed by feature id.

    Attributes
    ----------
        gff_dict: dict {feature_id: GffData object}
    """

    def __init__(self):
        self.gff_dict = {}

    # Internal helpers

    def _get_attr(self, feature_id, attr, default):
        """Returns attr from the GffData object for feature_id, or default if absent."""
        data = self.gff_dict.get(feature_id)
        return getattr(data, attr, default) if data else default

    def _iter_redundancy_candidates(self, filtered_ids):
        """Yields feature IDs that are candidates for redundancy checks. Skips IDs in filtered_ids,
         intergenic ('sequence_feature') and contig features """
        filtered_ids = initialize_parameter(filtered_ids, set())
        valid_features = self._get_all_potentially_redundant_features(filter_types={"region",
                                                                                    "chromosome",
                                                                                    "sequence_feature"})
        for feature_id, data in self.iterate():
            if feature_id not in filtered_ids and data.feature in valid_features:
                yield feature_id

    # Global queries

    def ids(self):
        """ Returns the keys of gff_dict """
        return self.gff_dict.keys()

    def all_data(self):
        """ Returns all GffData objects """
        return self.gff_dict.values()

    def iterate(self):
        """ Returns (feature_id, GffData) pairs for every entry """
        return self.gff_dict.items()

    # Single-feature queries

    def data(self, feature_id):
        """ Returns the GffData object for feature_id, or None """
        return self.gff_dict.get(feature_id)

    def data_tuple(self, feature_id):
        """ Returns the data tuple for feature_id """
        if self.data(feature_id) is not None:
            return self.data(feature_id).data_tuple
        else:
            print(f"{feature_id} not found in the node list!")
            return "", "", [], "", [], {}

    def list_positions(self, feature_id):
        """ Returns the list of (contig, strand, start, stop) tuples for feature_id """
        return self._get_attr(feature_id, "positions", [])

    def list_phases(self, feature_id):
        """ Returns the list of phases for feature_id """
        return self._get_attr(feature_id, "phase", [])

    def feature(self, feature_id):
        """ Returns the feature type string for feature_id """
        return self._get_attr(feature_id, "feature", "")

    def attributes(self, feature_id):
        """ Returns the attributes dict for feature_id."""
        return self._get_attr(feature_id, "attributes", {})

    # Data edition

    def update(self, data, replace=False):
        """ Inserts data into gff_dict, or extends positions/phases if the ID already exists.
        When replace is True the existing entry is overwritten entirely """
        if data.id not in self.gff_dict:
            self.gff_dict[data.id] = data
        elif replace:
            self.gff_dict[data.id] = data
        else:
            existing = self.gff_dict[data.id]
            existing.positions.extend(data.positions)
            existing.phase.extend(data.phase)

    def remove(self, key: str):
        """ Deletes key from gff_dict (.pop() usage avoids a check for presence of the key) """
        self.gff_dict.pop(key, None)

    def remove_set(self, set_keys: set[str]):
        """ Deletes every key in set_keys from gff_dict """
        for key in set_keys:
            self.remove(key)

    # Data subsets

    def ids_of_type(self, feature):
        """ Returns the set of IDs whose feature type is feature"""
        return {feature_id for feature_id, data in self.iterate() if data.feature == feature}

    def get_features_for_contig_of_type(self, contig: str, features: set[str] = None):
        """ Returns IDs on contig whose feature type is in features. If features is empty/None,
        all IDs on contig are returned """
        features = initialize_parameter(features, set())
        return {feature_id for feature_id, data in self.iterate()
                if (not features or data.feature.lower() in features) and self._feature_on_contig(feature_id, contig)}

    def _feature_on_contig(self, feature_id, test_contig):
        """ Returns True when feature_id has at least one position on test_contig """
        return any(contig == test_contig for contig, *_ in self.list_positions(feature_id))

    # Position / phase conversion to str

    def convert_list_positions_to_string(self, feature_id):
        """ Converts position tuples (contig, strand, start, end) to a '||'-joined string """
        parts = [f"{contig};{strand};{start};{end}" for contig, strand, start, end in self.list_positions(feature_id)]
        return "||".join(parts) if parts else ""

    def convert_list_positions_to_string_no_strand(self, feature_id):
        """ Converts position tuples to a '||'-joined string, omitting strand.
        Returns (positions_str, list_of_strands), or ('', []) when strand
        information is missing for any position """
        parts, strands = [], []
        for contig, strand, start, end in self.list_positions(feature_id):
            if strand not in ("+", "-"):
                return "", []  # early-return on missing strand info
            parts.append(f"{contig};{start};{end}")
            strands.append(strand)
        return ("||".join(parts), strands) if parts else ("", [])

    def convert_list_phases_to_string(self, feature_id):
        """ Converts the phase list to a '||'-joined string"""
        phases = self.list_phases(feature_id)
        return "||".join(phases) if phases else ""

    # Redundant feature detection

    def get_identical_features(self, filtered_ids=None):
        """ Returns a dict {positions_str: set_of_ids} for features sharing identical positions.
        Only groups with ≥ 2 members are included """
        identical_dict = {}
        for feature_id in self._iter_redundancy_candidates(filtered_ids):
            positions = self.convert_list_positions_to_string(feature_id)
            if positions:
                identical_dict = update_dict_of_sets(identical_dict, key=positions, value=feature_id)
        return filter_dict_of_sets(identical_dict, cutoff=1)

    def get_reverse_complementary_features(self, filtered_ids=None):
        """ Returns a dict {positions_str: set_of_(id1, id2)_pairs} for RC feature pairs """
        rc_dict = {}
        for feature_id in self._iter_redundancy_candidates(filtered_ids):
            result = self.convert_list_positions_to_string_no_strand(feature_id)
            if not result:
                continue
            positions, list_strand = result
            if positions and list_strand:
                rc_dict = update_dict_of_sets(rc_dict, key=positions, value=(feature_id, tuple(list_strand)))
        return self._filter_rc_dict(filter_dict_of_sets(rc_dict, cutoff=1))

    def _get_all_potentially_redundant_features(self, filter_types=None):
        """ Returns all feature types present in gff_dict, minus filter_types """
        filter_types = initialize_parameter(filter_types, set())
        return {data.feature for data in self.all_data() if data.feature} - filter_types

    @staticmethod
    def _filter_rc_dict(rc_dict):
        """ Keeps only pairs with identical positions but fully opposite strands.
        Input:  {positions: set( (feature_id, tuple_strands) )}
        Output: {positions: set( (node1, node2) )} """
        opposite = {"+": "-", "-": "+"}
        filtered_dict = {}
        for positions, set_tuples in rc_dict.items():
            for (node1, strands1), (node2, strands2) in combinations(set_tuples, 2):
                if all(s1 == opposite[s2] for s1, s2 in zip(strands1, strands2)):
                    filtered_dict = initialize_dict(filtered_dict, positions, set())
                    filtered_dict[positions].add((node1, node2))
        return filtered_dict

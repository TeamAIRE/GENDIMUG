from utils.helpers.init_functions import initialize_dict
from utils.helpers.id_functions import get_suffix
from utils.helpers.update_functions import update_dict_of_sets, reverse_dict_of_sets
from utils.helpers.interval_functions import is_match, group_overlapping_intervals


class CDSCaseHandler:
    """
    CDS-specific object that identifies particular cases of fragmentation/overlapping positions, as described below.
    """
    def __init__(self, cds, cds_data, exon_lookup):
        self.cds = cds
        self.cds_data = cds_data
        self.exon_by_coords = exon_lookup

        # storage of positions, phases and attributes derived from the cds_data tuple
        self.list_cds_positions = []
        self.list_cds_phases = []
        self.position_phase = {}
        self.attributes = {}
        if self.cds_data:
            _, _, self.list_cds_positions, _, self.list_cds_phases, self.attributes = cds_data
        if self.list_cds_positions and self.list_cds_phases:
            self.position_phase = {(pos[2], pos[3]): phase[0]
                                   for pos, phase in zip(self.list_cds_positions, self.list_cds_phases)}

        # storage attributes used depending on the CDS cases
        self.cds_exon_matches = {}
        self.exon_cds_matches = {}
        self.suffixes = []
        self.fragments = []
        self.data_tuples = []
        self.list_position_overlap_clusters = []
        self.list_same_length_clusters = []
        self.list_same_parent_clusters = []
        self.parent = ""

        # Boolean attributes for case identification
        self.exon_data = False
        self.is_fragmented = False
        self.has_overlaps = False
        self.has_slippage = False
        self.has_exon_linked_to_multiple_cdsparts = False

    def identify_case(self):
        """ Main entry point to analyze the current CDS"""
        self._test_cds()
        self._specific_lists()
        self._dispatch()

    def _test_cds(self):
        """Determines whether the cds is fragmented or has overlapping position intervals"""
        self.list_position_overlap_clusters = group_overlapping_intervals(self.list_cds_positions)

        # Test if there are clusters of length > 1, indicating overlapping intervals
        if any([len(cluster) > 1 for cluster in self.list_position_overlap_clusters]):
            self.has_overlaps = True

        # test for fragmented CDS
        if self.has_overlaps:
            if len(self.list_position_overlap_clusters) > 1:
                self.is_fragmented = True
        else:
            if len(self.list_cds_positions) > 1 and not self.has_slippage:
                self.is_fragmented = True

        # test for ribosomal (or other) slippage
        if "slippage" in self.attributes.get("exception", ""):
            self.has_slippage = True

        # test for exonic data
        if self.exon_by_coords:
            self.exon_data = True

        # test for the fragmentation of cds with exon data and slippage
        if self.exon_data:
            # compute matches between cds positions and exon positions, if exon_data
            self._build_matching_dict()
            if self.has_slippage:
                if (not self.has_overlaps
                        and any([len(cds_matches) > 1 for cds_matches in self.exon_cds_matches.values()])):
                    self.is_fragmented = True
                    self.has_exon_linked_to_multiple_cdsparts = True

    def _build_matching_dict(self):
        """ Matches each CDS position to an exon ID"""
        if self.exon_by_coords and self.list_cds_positions:
            for position in self.list_cds_positions:
                self.cds_exon_matches = update_dict_of_sets(self.cds_exon_matches,
                                                            key=position,
                                                            value=self._find_matching_exon(position))
        if all([exon for exon in self.cds_exon_matches.values()]):
            self.exon_cds_matches = reverse_dict_of_sets(self.cds_exon_matches)

    def _specific_lists(self):
        """ Generates lists of clusters by shared parent and by shared length"""
        if self.exon_data and not self.has_overlaps and self.is_fragmented and self.has_slippage:
            # list of clusters by shared parent

            parent_order = []
            clusters = {}
            for position in self.list_cds_positions:
                parent_exon = self._find_matching_exon(position)
                parent_order.append(parent_exon)
                clusters = initialize_dict(clusters, parent_exon, [])
                clusters[parent_exon].append(position)
            parent_order = self._deduplicate_list_in_order(parent_order)
            for parent in parent_order:
                self.list_same_parent_clusters.append(clusters[parent])

        elif not self.exon_data and not self.has_overlaps and self.is_fragmented and self.has_slippage:
            # list of clusters by shared fragment length

            current_length = 0
            cluster = []
            for i, position in enumerate(self.list_cds_positions):
                _, _, start, end = position
                length = abs(end - (start - 1))
                if length == current_length:
                    cluster.append(position)
                else:
                    if cluster:
                        self.list_same_length_clusters.append(cluster.copy())
                    cluster = []
                    current_length = length
                    cluster.append(position)
                # last fragment: append the last current cluster
                if i == len(self.list_cds_positions) - 1:
                    self.list_same_length_clusters.append(cluster.copy())

    def _dispatch(self):
        """ Boolean logic to identify CDS case """
        # key = (self.exon_data, self.has_overlaps, self.is_fragmented, self.has_slippage)
        # Groups of cases that share the same handler
        if self.exon_data:
            if self.has_overlaps:
                self._case_overlapping_exon()           # cases 1 & 2
            elif self.is_fragmented and self.has_slippage:
                self._case_slippage_exon()              # case 3
            elif self.is_fragmented:
                self._case_multi_exon()                 # case 4
            else:
                self._case_single_exon()                # case 5
        else:
            if self.has_overlaps and self.is_fragmented and self.has_slippage:
                self._case_fragmented_prokaryote_slippage()   # case 6
            elif self.is_fragmented and not self.has_slippage:
                self._case_fragmented_prokaryote()            # case 7
            elif self.is_fragmented and self.has_slippage:
                self._case_same_length_clusters()             # case 9
            else:
                self._case_simple_prokaryote()                # case 8

    # processing depending on case
    def _case_overlapping_exon(self):
        """Cases 1 & 2: exon + overlapping intervals, fragmented or not."""
        for i, group in enumerate(self.list_position_overlap_clusters):
            for position in group:
                parent_exon = self._find_matching_exon(position)
                if not parent_exon:
                    continue
                if self.is_fragmented:                   # case 1
                    suffix = get_suffix(parent_exon)
                    extra = {"part": str(suffix), "Parent": parent_exon}
                    self._build_fragment(suffix, group, extra, rename=True)
                else:                                    # case 2
                    suffix = i + 1
                    extra = {"Parent": parent_exon}
                    self._build_fragment(suffix, group, extra)
                break  # first exonic match is enough

    def _case_slippage_exon(self):
        """Case 3: exon + fragmented + slippage, cluster by shared parent."""
        for group in self.list_same_parent_clusters:
            for j, position in enumerate(group):
                parent_exon = self._find_matching_exon(position)
                if not parent_exon:
                    continue
                if not self.has_exon_linked_to_multiple_cdsparts:
                    suffix = get_suffix(parent_exon)
                    extra = {"part": str(suffix), "Parent": parent_exon}
                    self._build_fragment(suffix, group, extra, rename=True)
                    break
                else:
                    suffix = j + 1
                    extra = {"part": str(suffix), "Parent": parent_exon}
                    self._build_fragment(suffix, [position], extra, rename=True)

    def _case_multi_exon(self):
        """Case 4: regular eukaryote multi-exon CDS."""
        for position in self.list_cds_positions:
            parent_exon = self._find_matching_exon(position)
            if parent_exon:
                suffix = get_suffix(parent_exon)
                extra = {"part": str(suffix), "Parent": parent_exon}
                self._build_fragment(suffix, [position], extra, rename=True)

    def _case_single_exon(self):
        """Case 5: regular eukaryote single-exon CDS (with optional slippage)."""
        for position in self.list_cds_positions:
            parent_exon = self._find_matching_exon(position)
            if parent_exon:
                suffix = 1
                self._build_fragment(suffix, self.list_cds_positions, {"Parent": parent_exon})
                break

    def _case_fragmented_prokaryote_slippage(self):
        """Case 6: fragmented prokaryote CDS with slippage."""
        for i, group in enumerate(self.list_position_overlap_clusters):
            suffix = i + 1
            self._build_fragment(suffix, group, {"part": str(suffix)}, rename=True)

    def _case_fragmented_prokaryote(self):
        """Case 7: fragmented prokaryote CDS, no slippage."""
        for i, position in enumerate(self.list_cds_positions):
            suffix = i + 1
            self._build_fragment(suffix, [position], {"Part": str(suffix)}, rename=True)

    def _case_same_length_clusters(self):
        """Case 9: prokaryote, fragmented + slippage, cluster by shared length."""
        for i, group in enumerate(self.list_same_length_clusters):
            suffix = i + 1
            self._build_fragment(suffix, group, {"part": str(suffix)}, rename=True)

    def _case_simple_prokaryote(self):
        """Case 8: regular prokaryote CDS."""
        self._build_fragment(1, self.list_cds_positions, {})

    def _find_matching_exon(self, position: tuple[str, str, int, int]) -> str | None:
        """ Returns the exon ID matching in positions with the current CDS part"""
        _, _, cds_start, cds_end = position
        for (exon_start, exon_end), exon_id in self.exon_by_coords.items():
            if is_match((cds_start, cds_end), (exon_start, exon_end)):
                return exon_id
        return None

    def _get_phases(self, positions):
        """ Returns the matching list of phases for the provided list of positions """
        return [self.position_phase[(s, e)] for _, _, s, e in positions]

    def _build_fragment(self, suffix, positions, extra_attributes, rename=False):
        """Append one fragment to self.fragments / self.suffixes / self.data_tuples."""
        base = self.cds.replace("cds-", "cdspart-") if rename else self.cds
        new_id = f"{base}-{suffix}"
        self.suffixes.append(suffix)
        self.fragments.append(new_id)

        db, feature, _, score, _, _ = self.cds_data
        new_attrs = {**self.attributes, "ID": new_id, **extra_attributes}
        self.data_tuples.append((db, feature, positions, score, self._get_phases(positions), new_attrs))

    @staticmethod
    def _deduplicate_list_in_order(list_with_duplicates):
        """ Preserves the order of elements in the list but removes duplicates"""
        deduplicated_list = []
        seen = set()
        for item in list_with_duplicates:
            if item not in seen:
                deduplicated_list.append(item)
                seen.add(item)
        return deduplicated_list


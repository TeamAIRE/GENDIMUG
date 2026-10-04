from utils.classes.mappers.mapperbase import Mapper
from utils.classes.storage.gffdata import GffData
from utils.classes.sorting.nestingsorter import NestingSorter


class IntergenicMapper(Mapper):

    def map(self):
        """ Infers intergenic and terminal features from non-nested gene positions """
        self._map_intergenic_regions_to_contigs()
        self._export_gff(prefix="intergenic", set_ids=self.contig_storage.feature_types.get_ids("intergenic"))
        return self.contig_storage

    def _map_intergenic_regions_to_contigs(self):
        """ Updates the contig_dict and gff_dict dicts with intergenic features inferred from non-nested ORF coordinates
        and creates a set of intergenic ids """
        self.contig_storage.feature_types.initialize("adjacent_genes", "intergenic")
        for strand in ["+", "-"]:

            # order the features on this strand of this contig
            list_genes = self.contig_storage.positions.get_list(
                strand,
                keep_ids=self.contig_storage.feature_types.get_ids("non_nested"))

            # only keep disjoint sequential positions (i.e. stop1 < start2):
            intergenic_positions = self._sort_intervals(list_genes)

            for i_start, i_stop, id1, id2 in intergenic_positions:
                i_start, i_stop = sorted((i_start, i_stop))
                intergenic_id = f"intergenic-{self.contig}_{strand}_{i_start}_{i_stop}"
                feature = "sequence_feature"
                if strand == "-":
                    id1, id2 = id2, id1
                attributes = {"ID": intergenic_id,
                              "Name": intergenic_id,
                              "Note": "Inter-ORF non-coding region. Positions inferred from ORF coordinates",
                              "locus_tag": f"intergenic|{id1}|{id2}"}
                list_positions = [(self.contig, strand, i_start, i_stop)]
                data_tuple = ("inferred", feature, list_positions, ".", ["."], attributes)

                # update ContigStorage objects
                self.contig_storage.gff.update(
                    GffData((intergenic_id, data_tuple), database=self.database, is_tuple=True))
                self.contig_storage.positions.update(intergenic_id, list_positions)
                self.contig_storage.feature_types.update("intergenic", {intergenic_id})

    def _sort_intervals(self, list_intervals):
        """ Returns the sorted list of non-nested, non-overlapping intervals (start, end, id), identified by
        the surrounding genes. Updates ContigStorage.FeatureTypes with the set of gene pairs that are immediately
        adjacent (key "adjacent_genes") (e.g. no intergenic region between two ORFs)"""

        filtered_list: list[tuple[int, int, str, str]] = []
        list_intervals: list[tuple[int, int, str]]
        list_intervals, _ = NestingSorter(list_intervals).remove_nested()
        for i in range(1, len(list_intervals)):
            # Previous interval end
            previous_start, previous_end, id1 = list_intervals[i - 1]
            # Current interval start
            current_start, current_end, id2 = list_intervals[i]
            if previous_end == current_start and previous_start != current_end:
                # single position overlap, not detected by interval tree
                continue

            elif previous_end < current_start:
                # do not include overlapping genes

                if previous_end + 1 == current_start or previous_end == current_start - 1:
                    # adjacent features, pairs of adjacent features are skipped in the list
                    self.contig_storage.feature_types.update("adjacent_genes", {(id1, id2)})

                else:
                    filtered_list.append((previous_end + 1, current_start - 1, id1, id2))

        return filtered_list

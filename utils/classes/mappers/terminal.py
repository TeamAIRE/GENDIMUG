from utils.classes.mappers.mapperbase import Mapper
from utils.classes.storage.gffdata import GffData


class TerminalFeatureMapper(Mapper):
    """ Creates terminal intergenic region features """

    def __init__(self, global_parameters, contig_storage, contig_description):
        super().__init__(global_parameters, contig_storage)
        contig_start, contig_end, _, circular = contig_description[self.contig]
        self.contig_start = contig_start
        self.contig_end = contig_end
        self.circular = circular
        self.has_start_cutoff = self.start_cutoff > 0
        self.has_end_cutoff = self.end_cutoff > 0

    def map(self):
        """ Computes terminal intergenic regions and registers them in ContigStorage objects.
        Returns a dict mapping strand -> (first_terminal_id, last_terminal_id) for edge building. """
        terminal_anchors = {}

        for strand in ["+", "-"]:

            # gather all the positions of genic features for the current strand
            strand_features = self.contig_storage.positions.get_list(
                strand,
                keep_ids=self.contig_storage.feature_types.get_ids("genic"))

            if strand_features:

                # no first terminal feature if start cutoff
                if not self.has_start_cutoff:
                    first_feature = strand_features[0]
                    end_terminal, _, target = first_feature
                else:
                    end_terminal, target = "", ""

                # no last terminal feature if end cutoff
                if not self.has_end_cutoff:
                    last_feature = strand_features[-1]
                    _, start_terminal, source = last_feature
                else:
                    start_terminal, source = "", ""

                # map terminal features depending on circularity of the contig
                if not self.circular:
                    anchors = self._map_linear(strand, source, target, start_terminal, end_terminal)
                else:
                    anchors = self._map_circular(strand, source, target, start_terminal, end_terminal)

                terminal_anchors[strand] = anchors  # {source, target, terminal_ids}

        return self.contig_storage, terminal_anchors

    def _map_linear(self, strand, source, target, start_terminal, end_terminal) -> dict:
        """ Maps terminal regions for linear contigs """

        name = "terminal_intergenic_region"
        first_id = last_id = None

        # no first terminal feature if start cutoff
        if end_terminal:
            first_start, first_stop = (self.contig_start, end_terminal - 1) if end_terminal != self.contig_start \
                else (0, 0)
        else:
            first_start, first_stop = 0, 0

        # no last terminal feature if end cutoff
        if start_terminal:
            last_start, last_stop = (start_terminal + 1, self.contig_end) if start_terminal != self.contig_end \
                else (0, 0)
        else:
            last_start, last_stop = 0, 0

        if strand == "-":
            first_start, first_stop, last_start, last_stop = last_start, last_stop, first_start, first_stop
            source, target = target, source

        if first_start and first_stop:
            first_id = f"terminal-{self.contig}_{strand}_{first_start}_{first_stop}"
            self._register_terminal_regions(first_id, name, [(self.contig, strand, first_start, first_stop)])

        if last_start and last_stop:
            last_id = f"terminal-{self.contig}_{strand}_{last_start}_{last_stop}"
            self._register_terminal_regions(last_id, name, [(self.contig, strand, last_start, last_stop)])

        return {"source": source, "target": target, "first_id": first_id, "last_id": last_id}

    def _map_circular(self, strand, source, target, start_terminal, end_terminal) -> dict:
        """ Maps terminal regions for circular contigs """

        name = "terminal_intergenic_region"
        only_id = None

        if start_terminal and end_terminal:
            if end_terminal != self.contig_start or start_terminal != self.contig_end:
                if end_terminal == self.contig_start:
                    only_start = start_terminal + 1
                    positions = [(self.contig, strand, only_start, self.contig_end)]
                    only_id = f"terminal-{self.contig}_{strand}_{only_start}_{self.contig_end}"
                elif start_terminal == self.contig_end:
                    only_stop = end_terminal - 1
                    positions = [(self.contig, strand, self.contig_start, only_stop)]
                    only_id = f"terminal-{self.contig}_{strand}_{self.contig_start}_{only_stop}"
                else:
                    only_start, only_stop = start_terminal + 1, end_terminal - 1
                    positions = [(self.contig, strand, only_start, self.contig_end),
                                 (self.contig, strand, self.contig_start, only_stop)]
                    only_id = f"terminal-{self.contig}_{strand}_{only_start}_{only_stop}"

                self._register_terminal_regions(only_id, name, positions)

            if strand == "-":
                source, target = target, source

            return {"source": source, "target": target, "only_id": only_id}

        # no terminal feature if first and last genic features not present
        return {"source": "", "target": "", "only_id": ""}

    def _register_terminal_regions(self, feature_id, feature, list_positions):
        """ Updates ContigStorage objects with data corresponding to inferred extragenic terminal regions """
        correct_list_positions = []
        for (contig, strand, start, stop) in list_positions:
            if stop < start:
                start, stop = stop, start
            correct_list_positions.append((contig, strand, start, stop))

        attributes = {"ID": feature_id, "Name": feature_id,
                      "Note": "Terminal non-coding region. Positions inferred from ORF and contig coordinates"}
        data_tuple = ("inferred", feature, correct_list_positions, ".", ["."] * len(list_positions), attributes)
        self.contig_storage.gff.update(GffData((feature_id, data_tuple), database=self.database, is_tuple=True))
        self.contig_storage.positions.update(feature_id, correct_list_positions)
        self.contig_storage.feature_types.update("intergenic", {feature_id})

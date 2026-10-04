from utils.classes.storage.gffdata import GffData
from utils.classes.sorting.intervalmatcher import IntervalMatcher
from utils.helpers.update_functions import update_interval_dict_from_positions, compatible_contig_dicts
from utils.helpers.interval_functions import matches_generator


class DomainMapper:
    """Global manager of the genomic mapping of proteic domains/motifs."""

    def __init__(self, cds_fragments, cds_motifs_data, cds_with_motifs):
        self.contig_storage = None
        self.cds_fragments = cds_fragments
        self.cds_motifs_data = cds_motifs_data
        self.cds_with_motifs = cds_with_motifs

    def map(self, contig_storage):
        """ Main entry point to map all domains to genomic coordinates and wire parent-child links in the og_tree
        for every CDS on the current contig. Returns and updated ContigStorage object"""
        self.contig_storage = contig_storage

        # process each CDS with identified motifs
        for cds_id in self.cds_with_motifs:
            dict_domains = self.cds_motifs_data.get(cds_id, {})
            suffixes = self.cds_fragments.get(cds_id, [])
            if not suffixes:
                continue

            self._process_cds(cds_id, suffixes, dict_domains)

        # include ontogenetic links CDS fragment -> domain
        self.contig_storage.og_tree.update_child_parent_gene()

        return self.contig_storage

    # Per-CDS processing

    def _process_cds(self, cds_id: str, suffixes: list, dict_domains: dict):
        """Runs the full mapping pipeline for a single CDS."""
        list_cds_parts = self._build_cds_part_names(cds_id, suffixes)
        cds_part_interval_dict = self._build_cds_part_intervals(list_cds_parts)
        domain_interval_dict = self._map_domains_on_genome(list_cds_parts, dict_domains)
        self._link_parts_to_domains(cds_part_interval_dict, domain_interval_dict)

    def _build_cds_part_intervals(self, list_cds_parts: list) -> dict:
        """ Builds interval dict for all CDS parts. Returns {contig: {strand: {(start, stop): set(cds_part ids)}}} """
        interval_dict = {}
        for cds_part in list_cds_parts:
            positions = self.contig_storage.gff.list_positions(cds_part)
            if positions:
                interval_dict = update_interval_dict_from_positions(cds_part, positions, interval_dict=interval_dict)
        return interval_dict

    @staticmethod
    def _build_cds_part_names(cds_id: str, suffixes: list) -> list[str]:
        """Returns the list of CDS part identifiers for the given suffixes."""
        if len(suffixes) > 1:
            cds_id = cds_id.replace("cds-", "cdspart-").replace("CDS-", "cdspart-")
        return [f"{cds_id}-{suffix}" for suffix in suffixes]

    # Domain to genome mapping

    def _map_domains_on_genome(self, list_cds_parts: list, dict_domains: dict) -> dict:
        """ Converts amino-acid domain boundaries to genomic coordinates. Skips the CDS entirely if any part is absent
        from the GFF (incorrect FASTA entry). Returns a domain interval dict
        {contig: {strand: {(start, stop): set(domain ids)}}} """
        domain_interval_dict = {}
        known_ids = self.contig_storage.gff.ids()

        if not all(part in known_ids for part in list_cds_parts):
            return domain_interval_dict  # incorrect CDS -> skip silently

        gdict = self._get_nucleotide_to_genomic_map(list_cds_parts)

        self.contig_storage.feature_types.update("domains", set(dict_domains))

        for domain, (prot, tool, feature, pstart, pstop, pscore, _, pphase, pattributes) in dict_domains.items():
            coding_start = pstart * 3 - 3
            coding_stop = pstop * 3 - 1

            filtered_gdict = self._filter_genomic_dict(gdict, range(coding_start, coding_stop + 1))
            list_genomic_positions = self._split_non_consecutive_positions_as_tuples(sorted(filtered_gdict.values()))

            data_tuple = (tool,
                          feature,
                          list_genomic_positions,
                          pscore,
                          [pphase] * len(list_genomic_positions),
                          pattributes)

            self.contig_storage.gff.update(
                GffData((domain, data_tuple), database=self.contig_storage.database, is_tuple=True))
            domain_interval_dict = update_interval_dict_from_positions(domain,
                                                                       list_genomic_positions,
                                                                       interval_dict=domain_interval_dict)
        return domain_interval_dict

    def _get_nucleotide_to_genomic_map(self, cds_parts: list) -> dict:
        """Returns {nucleotide index in CDS: (genomic position, contig, strand)}."""
        gdict = {}
        c = 0
        for part in cds_parts:
            for contig, strand, start, stop in self.contig_storage.gff.list_positions(part) or []:
                for j in range(start, stop + 1):
                    gdict[c] = (j, contig, strand)
                    c += 1
        return gdict

    @staticmethod
    def _filter_genomic_dict(gdict, range_object):
        """ Filters a genomic dict {coding_position: genomic position} for a provided range of coding positions """
        return {k: v for k in range_object if (v := gdict.get(k))}

    @staticmethod
    def _split_non_consecutive_positions_as_tuples(list_of_positions):
        """ Splits a list of sorted (pos, contig, strand) tuples into a list of (contig, strand, start, stop) spans """
        if not list_of_positions:
            return []
        final = []
        # Seed the first group
        seg_start_pos, seg_contig, seg_strand = list_of_positions[0]
        seg_end_pos = seg_start_pos
        for pos, contig, strand in list_of_positions[1:]:
            contiguous = (abs(pos - seg_end_pos) == 1
                          and contig == seg_contig
                          and strand == seg_strand)
            if contiguous:
                seg_end_pos = pos
            else:
                final.append((seg_contig, seg_strand, seg_start_pos, seg_end_pos))
                seg_start_pos, seg_end_pos, seg_contig, seg_strand = pos, pos, contig, strand
        final.append((seg_contig, seg_strand, seg_start_pos, seg_end_pos))
        return final

    # Ontogenetic tree edges

    def _link_parts_to_domains(self, part_interval_dict: dict, domain_interval_dict: dict):
        """Registers parent (CDS part) -> child (domain) edges in the og_tree."""
        for contig in part_interval_dict:
            for strand in part_interval_dict[contig]:
                if not compatible_contig_dicts(part_interval_dict, domain_interval_dict, contig, strand):
                    continue

                part_dict = part_interval_dict[contig][strand]
                dom_dict = domain_interval_dict[contig][strand]

                matching_dict = IntervalMatcher(part_dict, dom_dict).match()

                for part, _, dom, _ in matches_generator(matching_dict, part_dict, dom_dict):
                    gene_id = self.contig_storage.og_tree.get_key_gene(part)
                    self.contig_storage.og_tree.add_edges(gene_id, [(part, dom)])
                    self._transfer_attributes_to_domain(part, dom)

    def _transfer_attributes_to_domain(self, part, dom):
        cds_attributes = self.contig_storage.gff.attributes(part)
        if cds_attributes.get("exception", ""):
            db, feature, list_positions, score, phases, dom_attributes = self.contig_storage.gff.data_tuple(dom)
            new_dom_attributes = {k: v for k, v in dom_attributes.items()}
            new_dom_attributes["exception"] = cds_attributes["exception"]
            new_tuple = (db, feature, list_positions, score, phases, new_dom_attributes)
            new_data = GffData((dom, new_tuple),
                               database=self.contig_storage.database,
                               is_tuple=True)
            self.contig_storage.gff.update(new_data, replace=True)
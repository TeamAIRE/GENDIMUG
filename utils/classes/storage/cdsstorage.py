from utils.classes.mappers.domain import DomainMapper
from utils.helpers.interval_functions import is_overlapping


class CdsStorage:
    """
    Global CDS-specific storage object to centralize the passage of specific storage objects to various functions
    Attributes
    ----------
        cds: set[str], set of CDS IDs
        cds_to_delete: set[str], set of CDS IDs to remove from the ContigStorage.Gff object
        cds_fragments: dict, {cds: list of cds part suffixes}
        cds_motifs_data: dict of dicts {cds ID: {motif ID: tuple of interproscan GFF3 data}
        contig_cds: dict {contig ID: set[str] of CDS IDS from this contig}
        domain_mapper: DomainMapper object
    """
    def __init__(self):
        self.cds = set()
        self.cds_to_delete = set()
        self.cds_fragments = {}
        self.cds_motifs_data = {}
        self.contig_cds = {}
        self.domain_mapper = None

    def map_domains(self, contig_storage):
        """Computes genomic positions for proteic motifs"""
        contig = contig_storage.contig
        cds_with_motifs = self.get_cds_with_motifs_for_contig(contig)
        self.domain_mapper = DomainMapper(self.cds_fragments, self.cds_motifs_data, cds_with_motifs)
        contig_storage = self.domain_mapper.map(contig_storage)
        return contig_storage

    def get_cds_with_motifs_for_contig(self, contig):
        """Returns the set of CDS associated with proteic motifs for the provided contig"""
        return set(self.cds_motifs_data) & set(self.contig_cds[contig])

    # edges between cds parts
    def get_edges_between_cds_parts(self, contig_storage):
        """ Computes sequential edges between cds parts """
        edges = set()
        for cds in self.contig_cds[contig_storage.contig]:
            suffixes = self.cds_fragments[cds]
            if len(suffixes) > 1:
                for i in range(len(suffixes) - 1):
                    cds_root = cds.replace("cds-", "cdspart-")
                    cds_up = f"{cds_root}-{suffixes[i]}"
                    positions_up = contig_storage.gff.list_positions(cds_up)
                    cds_down = f"{cds_root}-{suffixes[i + 1]}"
                    positions_down = contig_storage.gff.list_positions(cds_down)
                    for (_, _, start_up, end_up) in positions_up:
                        for (_, _, start_down, end_down) in positions_down:
                            if not is_overlapping((start_up, end_up), (start_down, end_down)):
                                # overlaps are computed later
                                edges.add((cds_up, cds_down, "cis-sequential"))
        return edges

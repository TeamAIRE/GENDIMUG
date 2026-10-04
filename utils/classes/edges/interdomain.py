from utils.classes.edges.intracluster import EdgeIntra
from utils.helpers.interval_functions import convert_to_interval_dict


class InterdomainEdges:

    def __init__(self, contig_storage, cds_storage):
        self.contig_storage = contig_storage
        self.cds_storage = cds_storage
        self.contig = self.contig_storage.contig

    def get_edges(self):
        # compute edges between domains
        edges_between_domains = self._get_inter_domain_edges()
        self.contig_storage.edge_db.save_edges(edges_between_domains)

    def _get_inter_domain_edges(self):
        """ Computes the edges between domains for all cds """
        edges = set()
        cds_with_motifs = self.cds_storage.get_cds_with_motifs_for_contig(self.contig)
        for cds in cds_with_motifs:
            dom_dict = self.cds_storage.cds_motifs_data[cds]
            if dom_dict:
                domain_set = set()
                for domain, (prot, tool, feature, pstart, pend, pscore, pstrand, pphase,
                             pattributes) in dom_dict.items():
                    domain_set.add((pstart, pend, domain))
                domain_interval_dict = convert_to_interval_dict(domain_set)
                new_edges = EdgeIntra(dict1=domain_interval_dict).get_edges()
                for edge in new_edges:
                    edges.add(edge)
        return edges

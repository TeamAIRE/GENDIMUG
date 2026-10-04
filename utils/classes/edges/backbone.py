from utils.classes.edges.sequential import Sequential
from utils.helpers.id_functions import get_prefix
from utils.helpers.init_functions import initialize_parameter


class EdgeBackbone:
    """ Computes backbone edges, genic-intergenic if gene == False (default), gene-only otherwise """
    def __init__(self, contig_storage, contig_description=None, gene=False):
        self.intergenic = not gene
        self.contig_storage = contig_storage
        self.contig_description = initialize_parameter(contig_description, {})
        self.edges = set()
        self.min_max = {}

    def get_edges(self):
        """ Main entry point for backbone edge building """

        # update ContigStorage.positions with all other mapped backbone features
        # (excluding features in the set "additional features"):
        all_mapped_ids = set(self.contig_storage.gff.ids()) - self.contig_storage.feature_types.get_ids(
            "additional_features")
        self.contig_storage.update_positions_for_ids(all_mapped_ids)

        if self.intergenic:
            # build intergenic edges with only non-nested gene-like features
            self._get_genic_intergenic_edges()

        else:
            # create gene-only backbone edges
            self._get_genic_sequential_edges()

        # add adjacent gene edges (no intergenic region)
        self._get_adjacent_genes_edges()
        return self.edges

    def _get_genic_intergenic_edges(self):
        """ Saves non-nested gene - intergenic relationships as edges ( + edgetype) """
        print(" + (gene <> intergenic/terminal feature) edges...")
        for intergenic_id in self.contig_storage.feature_types.get_ids("intergenic"):
            if get_prefix(intergenic_id) != "terminal":
                attributes = self.contig_storage.gff.attributes(intergenic_id)
                gene1, gene2 = attributes["locus_tag"].split("|")[1:]
                for n1, n2 in [(gene1, intergenic_id), (intergenic_id, gene2)]:
                    self.edges.add((n1, n2, "cis-sequential"))

    def _get_genic_sequential_edges(self):
        """ Saves gene (without the nested genes) - gene relationships as ( + edge type) edges """
        print(" + (gene <> gene) sequential edges...")
        non_nested_genes = self.contig_storage.feature_types.get_ids("non_nested")
        for strand in ["+", "-"]:
            interval_dict = self.contig_storage.positions.get_dict(strand, keep_ids=non_nested_genes)
            self._get_min_max(strand, interval_dict)
            sequential_edges = Sequential(dict1=interval_dict,
                                          strand1=strand).get_edges()
            for edge in sequential_edges:
                self.edges.add(edge)
            # terminal edges are added if the contig is circular
            self._compute_closing_genic_edges(strand)

    def _get_adjacent_genes_edges(self):
        """ Saves adjacent genes relationships as ( + sequential) edges """
        adjacent = self.contig_storage.feature_types.get_ids("adjacent_genes")
        for id1, id2 in adjacent:
            [(_, strand, _, _)] = self.contig_storage.gff.list_positions(id1)
            if strand == "+":
                self.edges.add((id1, id2, "cis-sequential"))
            elif strand == "-":
                self.edges.add((id2, id1, "cis-sequential"))

    def _get_min_max(self, strand, interval_dict):
        """ Saves min and max gene positions in the interval_dict for the current strand """
        start_positions = sorted(interval_dict.keys())
        if start_positions:
            first_start, first_end = start_positions[0]
            first_gene = list(interval_dict[(first_start, first_end)])[0]
        else:
            first_start, first_end, first_gene = 0, 0, ""
        end_positions = sorted(interval_dict.keys(), key=lambda x: x[1])
        if end_positions:
            last_start, last_end = end_positions[-1]
            last_gene = list(interval_dict[(last_start, last_end)])[0]
        else:
            last_start, last_end, last_gene = 0, 0, ""
        self.min_max[strand] = ((first_start, first_end, first_gene), (last_start, last_end, last_gene))

    def _compute_closing_genic_edges(self, strand):
        """ Determines the edges between extremity features and terminal regions in circular contigs,
        depending on strand """
        if self.contig_description:
            contig_start, contig_stop, _, circular = self.contig_description[self.contig_storage.contig]
            if circular:  # create terminal edges to have a circular network
                (first_start, first_end, first_gene), (last_start, last_end, last_gene) = self.min_max[strand]
                if first_gene and last_gene:
                    node1, node2 = (last_gene, first_gene) if strand == "+" else (first_gene, last_gene)
                    self.edges.add((node1, node2, "cis-sequential"))

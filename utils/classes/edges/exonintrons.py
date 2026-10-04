from utils.helpers.id_functions import reorder_parts, get_suffix


class EdgeIntrons:
    """
    Edge factory for intron inference and computation of parent-exon ontogenetic edges and exon-intron edges.
    Attributes:
    ----------
        contig_storage: ContigStorage object
        gene_genepart: dict {gene ID : list of geneparts}
        edges: set of edges as tuples (node1, node2, edgetype)
        trans_spliced_exons: set[str] of trans-spliced exon IDs
    """

    def __init__(self, contig_storage, gene_genepart):
        self.contig_storage = contig_storage
        self.gene_genepart = gene_genepart
        self.edges = set()
        self.trans_spliced_exons = set()

    def get_edges(self):
        """ Returns exon-intron edges (or exon-exon edges when trans-splicing occurs) """
        # Build exon-intron edges
        self._get_exon_intron_edges()
        # add trans-splicing edges
        self._get_exon_trans_splicing_edges()
        return self.edges

    def _get_exon_intron_edges(self):
        """ Builds exon - intron edges """
        # identify exons for trans-spliced ORFs
        for gene in self.gene_genepart:
            self.trans_spliced_exons.update(self.contig_storage.og_tree.children(gene, {"exon"}))

        for intron_id in self.contig_storage.feature_types.get_ids("intronic"):

            # edge previous exon - intron
            exon_before = intron_id.replace("intron", "exon")
            nb = get_suffix(exon_before)
            if nb:
                self.edges.add((exon_before, intron_id, "cis-sequential"))

                # edge intron - following exon
                exon_after = f"{'-'.join(exon_before.split('-')[:-1])}-{nb + 1}"  # increment the exon count
                self.edges.add((intron_id, exon_after, "cis-sequential"))

    def _get_exon_trans_splicing_edges(self):
        """ Builds edges between trans-spliced exons (edgetype == trans-splicing) """

        for gene in self.gene_genepart:

            set_rna = self.contig_storage.og_tree.children(gene, {"rna"})
            if set_rna:
                for rna in set_rna:
                    exons = self.contig_storage.og_tree.filtered_children_from(gene, rna, {"exon"})
                    self._add_edges_trans_splicing(exons)

            else:
                exons = self.contig_storage.og_tree.children(gene, {"exon"})
                self._add_edges_trans_splicing(exons)

    def _add_edges_trans_splicing(self, list_exons):
        """ Adds trans-splicing exon-exon edges if there are more than one exon """
        if len(list_exons) > 1:
            exons = reorder_parts(list_exons)
            for i in range(len(exons) - 1):
                self.edges.add((exons[i], exons[i + 1], "trans-splicing"))

from utils.classes.mappers.mapperbase import Mapper
from utils.classes.storage.gffdata import GffData
from utils.helpers.update_functions import filter_dict_of_sets
from utils.helpers.id_functions import reorder_parts


class IntronMapper(Mapper):

    def __init__(self, global_parameters, contig_storage, gene_genepart):
        super().__init__(global_parameters=global_parameters, contig_storage=contig_storage)
        self.gene_genepart = gene_genepart
        self.parent_to_exons = {}
        self.filtered_parent_to_exons = {}
        self.gene_rna = {}

    def map(self):
        """ Main entry to infer introns: positions, ontogenetic relationships"""

        # get all exon parents in the og_tree dict. Trans-spliced ORFs are excluded.
        self._map_parent_exons_links()

        # Inference of intronic regions from exons coordinates and storage/update of ContigStorage objects
        if self.filtered_parent_to_exons:
            # map intronic positions
            self._map_introns_to_contig()
            # optional export of a GFF3 intron file
            self._export_gff(prefix="intronic", set_ids=self.contig_storage.feature_types.get_ids("intronic"))

        return self.contig_storage

    def _map_parent_exons_links(self):
        """ Collects all RNAs in the ContigStorage.OgTree object as a dict
        {rna: {sets of children by type}} (in eukaryotes usually) or {gene: {sets of children by type}} (prokaryotes).
        The trans-spliced ORFs identified in the gene_genepart dict are excluded """
        for gene_id in self.contig_storage.og_tree.get_genes():
            if gene_id not in self.gene_genepart:  # excludes the trans-spliced RNAs
                self.gene_rna[gene_id] = self.contig_storage.og_tree.children(gene_id,
                                                                              {"rna"})

        for gene_id, set_rna in self.gene_rna.items():

            if set_rna:  # eukaryotic gff annotation, exons are linked to RNAs
                for rna_id in set_rna:
                    self.parent_to_exons[rna_id] = self.contig_storage.og_tree.filtered_children_from(gene_id,
                                                                                                      rna_id,
                                                                                                      {"exon"})

            else:  # prokaryotic gff annotation, exons are linked to genes
                self.parent_to_exons[gene_id] = self.contig_storage.og_tree.children(gene_id,
                                                                                     {"exon"})

            # to infer introns, keep the RNAs with more than one exon
            self.filtered_parent_to_exons = filter_dict_of_sets(self.parent_to_exons, cutoff=1)

    def _map_introns_to_contig(self):
        """ Updates the ContigStorage objects with intronic features inferred from exon coordinates.
        Trans-spliced ORFs are excluded from intron mapping """

        self.contig_storage.feature_types.initialize("intronic")

        for parent_id, set_exons in self.filtered_parent_to_exons.items():
            gene_id = self.contig_storage.og_tree.get_key_gene(parent_id)

            if gene_id:
                new_links = []  # storage of new parent_id > intron_id links to update the og_tree
                exons = reorder_parts(set_exons)  # converts the set of exon ids to an ordered list of exon ids

                for i in range(len(exons) - 1):

                    # naming and positions of intronic feature
                    exon_before, exon_after = exons[i], exons[i + 1]
                    list_positions_before, list_positions_after = (self.contig_storage.gff.list_positions(exon_before),
                                                                   self.contig_storage.gff.list_positions(exon_after))

                    if list_positions_before and list_positions_after:
                        [(contig, strand, start_before, end_before)] = list_positions_before
                        [(_, _, start_after, end_after)] = list_positions_after
                        if start_before and end_before and start_after and end_after:
                            i_start, i_stop = self._get_positions_from_borders((start_before, end_before),
                                                                               (start_after, end_after),
                                                                               strand)
                            intron_id = f"{exon_before.replace('exon', 'intron')}"
                            new_links.append((parent_id, intron_id))

                            # update feature_sets
                            self.contig_storage.feature_types.update("intronic", {intron_id})

                            # update the Gff object
                            locus_tag = gene_id
                            attributes = {"ID": intron_id,
                                          "Parent": parent_id.replace("rna-", "").replace("gene-", ""),
                                          "Name": intron_id,
                                          "Note": "Positions inferred from exons coordinates",
                                          "locus_tag": locus_tag}
                            data_tuple = ("inferred", "intron", [(contig, strand, i_start, i_stop)], ".", ["."],
                                          attributes)
                            self.contig_storage.gff.update(GffData((intron_id, data_tuple),
                                                                   database=self.contig_storage.database,
                                                                   is_tuple=True))

                            # update IntervalData object
                            self.contig_storage.update_positions(intron_id)

                # add ontogenetic edges for all introns in the OgTree object
                self.contig_storage.og_tree.add_edges(gene_id, new_links)

        # update the child to parent gene dict with introns from the created ontogenetic edges in the OgTree object
        self.contig_storage.og_tree.update_child_parent_gene()

    @staticmethod
    def _get_positions_from_borders(coord1: tuple[int, int], coord2: tuple[int, int], strand: str):
        """ Computes the coordinates of a feature from the borders of adjacent features and the strand information"""
        (start1, end1), (start2, end2) = coord1, coord2
        if strand == "-":
            end1, start2 = end2, start1
        return end1 + 1, start2 - 1

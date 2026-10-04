from utils.helpers.init_functions import initialize_parameter
from utils.helpers.id_functions import get_prefix
from utils.helpers.positions_functions import (get_min_max_positions,
                                               compatible_positions)


class EdgeAttributes:
    """
    Manages the addition of edge attributes and node / edge encoding for a set of edges (node1, node2, edgetype)
    Attributes
    ----------
        edges: set of tuples (node1, node2, edgetype)
    """

    def __init__(self,
                 edges: set[tuple[str, str, str]] | None = None
                 ):
        self.edges = initialize_parameter(edges, set())
        self.contig_description = {}
        self.encode = False

    def add_attributes(self, contig_storage, contig_description, encode: bool = False):
        """
        Computes the edgeclass, overlap_length, and distance for each edge in self.edges.
        Returns a set of str formatted edges.
        :param contig_storage: ContigStorage object
        :param contig_description: dict {contig ID: data tuple (start, end, name, bool is_circular)}
        :param encode: bool, whether to encode numerically the nodes, edge types and edge classes
        :return: set(str)
        """
        # populate attributes
        self.contig_description = contig_description
        self.encode = encode

        # initialize final edges set and optional encoding indices
        final_edges = set()
        et_i, ec_i = 0, 0

        # decorate each edge
        for n1, n2, edgetype in self.edges:

            same_strand = False if edgetype.startswith("trans-") or edgetype.startswith("reverse-") else True

            edgeclass = self._get_edgeclass(n1, n2, edgetype)

            # optional encoding of edgetype and edgeclass
            contig_storage, et_i, ec_i = self._increment_encoding_index(contig_storage, edgetype, et_i, edgeclass, ec_i)

            # calculation of overlap_length and distance edge attributes
            list_positions1 = contig_storage.gff.list_positions(n1)
            list_positions2 = contig_storage.gff.list_positions(n2)
            overlap_length, distance = self._get_overlap_and_distance(edgetype,
                                                                      list_positions1,
                                                                      list_positions2,
                                                                      same_strand=same_strand)

            # final edge formatting
            tuple_edge = (n1, n2, edgetype, edgeclass, overlap_length, distance)
            edge = self._conditional_encoding(contig_storage, tuple_edge)
            final_edges.add(edge)

            # CORRECTION for rare trans-spliced ORFS:
            # add a nesting spatial edge from gene part to trans-spliced mRNA
            # (only if the list of positions for the gene part
            # is shorter than the list of positions for the trans-spliced mRNA)
            # this is triggered by a 'genepart_parent_rna' edgeclass

            if edgeclass == "genepart_parent_rna" and len(list_positions1) < len(list_positions2):
                ts_edgetype = "cis-nested" if same_strand else "trans-nested"
                ts_edgeclass = self._get_edgeclass(n1, n2, ts_edgetype)
                # calculation of overlap_length and distance edge attributes
                ts_overlap_length, ts_distance = self._get_overlap_and_distance(edgetype,
                                                                                list_positions1,
                                                                                list_positions2,
                                                                                same_strand=same_strand)
                # final additional edge formatting
                tuple_edge = (n1, n2, ts_edgetype, ts_edgeclass, ts_overlap_length, ts_distance)
                # optional encoding of edgetype and edgeclass
                contig_storage, et_i, ec_i = self._increment_encoding_index(contig_storage,
                                                                            ts_edgetype,
                                                                            et_i,
                                                                            ts_edgeclass,
                                                                            ec_i)

                edge = self._conditional_encoding(contig_storage, tuple_edge)
                final_edges.add(edge)

        return final_edges

    def _get_overlap_and_distance(self, edgetype, list_positions1, list_positions2, same_strand=False):
        """ Returns overlap length and distance for the current edge """
        overlap_length = 0
        distance = 0
        if edgetype in {"cis-sequential", "trans-sequential", "trans-splicing"}:
            distance = self._compute_distance(list_positions1, list_positions2,
                                              same_strand=same_strand)
        else:
            overlap_length = self._get_overlap_length_from_positions(list_positions1, list_positions2,
                                                                     same_strand=same_strand)
        return overlap_length, distance

    def _conditional_encoding(self, contig_storage, tuple_edge):
        """ Returns the edge as a str, with optional encoding (if self.encode==True)"""
        n1, n2, edgetype, edgeclass, overlap_length, distance = tuple_edge
        if self.encode:
            alias1, alias2 = contig_storage.num_index[n1], contig_storage.num_index[n2]
            et_alias = contig_storage.edgetype_index[edgetype]
            ec_alias = contig_storage.edgeclass_index[edgeclass]
            return f"{alias1},{alias2},{et_alias},{ec_alias},{overlap_length},{distance}\n"
        else:
            return f"{n1},{n2},{edgetype},{edgeclass},{overlap_length},{distance}\n"

    def _increment_encoding_index(self, contig_storage, edgetype, et_i, edgeclass, ec_i):
        """"""
        if self.encode:
            if edgetype not in contig_storage.edgetype_index:
                contig_storage.edgetype_index[edgetype] = et_i
                et_i += 1
            if edgeclass not in contig_storage.edgeclass_index:
                contig_storage.edgeclass_index[edgeclass] = ec_i
                ec_i += 1
        return contig_storage, et_i, ec_i

    @staticmethod
    def _get_edgeclass(n1, n2, edgetype):
        """
        Returns the class of the edge (n1, n2, edgetype)
        :param n1: str, node feature1
        :param n2: str, node feature2
        :param edgetype: str, type of edge
        :return: str, class of edge
        """
        return f"{get_prefix(n1)}_{edgetype}_{get_prefix(n2)}"

    # computation of distance between features
    def _compute_distance(self, list_positions1, list_positions2, same_strand=True):
        """
        Returns the distance in positions between two nodes.
        :param list_positions1: list of tuples (contig, strand, start, end) for node 1
        :param list_positions2: list of tuples (contig, strand, start, end) for node 2
        :param same_strand: bool
        :return: int
        """
        contigs = set(c for c, _, _, _ in list_positions1 + list_positions2)
        if len(contigs) != 1:  # features on different contig (very rare)
            distance = "none"
        else:
            contig = list(contigs)[0]
            if contig in self.contig_description:
                _, max_position, _, circular = self.contig_description[contig]
                if circular and max_position:  # circular contig
                    distance = self._get_distance_circular(list_positions1,
                                                           list_positions2,
                                                           max_position,
                                                           same_strand=same_strand)
                else:  # linear contig
                    distance = self._get_distance(list_positions1,
                                                  list_positions2,
                                                  same_strand=same_strand)
            else:  # contig has no data (should never happen)
                distance = "none"
        return distance

    @staticmethod
    def _get_distance(list_positions1, list_positions2, same_strand=True):
        """
        Returns the distance between two features depending on their lists of positions
        :param list_positions1: list(tuple(contig, strand, start, end))
        :param list_positions2: list(tuple(contig, strand, start, end))
        :param same_strand: bool, whether the features are on the same strand
        :return: int
        """
        if compatible_positions(list_positions1, list_positions2, same_strand):
            min1, max1 = get_min_max_positions(list_positions1)
            min2, max2 = get_min_max_positions(list_positions2)
            for end_max, start_min in [(max1, min2), (max2, min1)]:
                if end_max < start_min:
                    return start_min - end_max - 1
        return 0

    def _get_distance_circular(self, list_positions1, list_positions2, max_position, same_strand=True):
        """
        Returns the distance between two features depending on their lists of positions
        :param list_positions1: list(tuple(contig, strand, start, end))
        :param list_positions2: list(tuple(contig, strand, start, end))
        :param same_strand: bool, whether the features are on the same strand
        :return: int
        """
        min_distance = [self._get_distance(list_positions1, list_positions2, same_strand=same_strand)]
        if compatible_positions(list_positions1, list_positions2, same_strand):
            min1, max1 = get_min_max_positions(list_positions1)
            min2, max2 = get_min_max_positions(list_positions2)
            for end_max, start_min in [(max1, min2), (max2, min1)]:
                if end_max > start_min:
                    min_distance.append((max_position - end_max) + start_min)
            return min(min_distance)
        return 0

    # computation of overlap length
    @staticmethod
    def _get_overlap_length(coord1: tuple[int, int], coord2: tuple[int, int]) -> int:
        """
        Calculates the length of the overlap between two positions intervals (inclusive)
        :param coord1: tuple of int (start1, end1) with start1 <= end1
        :param coord2: tuple of int (start2, end2) with start2 <= end2
        :return: int, length of the overlap
        """
        (start1, end1), (start2, end2) = coord1, coord2
        if not (start2 > end1 or start1 > end2):
            start_o, end_o = max(start1, start2), min(end1, end2)
            return end_o - (start_o - 1)
        return 0

    def _get_overlap_length_from_positions(self, list_positions1, list_positions2, same_strand=True):
        """
        Calculates the length of the overlap between two list positions
        :param list_positions1: list(tuple(contig, strand, start, end))
        :param list_positions2: list(tuple(contig, strand, start, end))
        :return: int, length of the overlap
        """
        overlap_length = 0
        for contig1, strand1, start1, end1 in list_positions1:
            for contig2, strand2, start2, end2 in list_positions2:
                if compatible_positions([(contig1, strand1, start1, end1)],
                                        [(contig2, strand2, start2, end2)],
                                        same_strand=same_strand):
                    overlap_length += self._get_overlap_length((start1, end1), (start2, end2))
        return overlap_length

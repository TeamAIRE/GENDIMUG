from utils.classes.edges.edgefactory import EdgeFactory
from utils.classes.sorting.nestingsorter import NestingSorter
from utils.helpers.init_functions import initialize_parameter


class EdgeIntra(EdgeFactory):
    """ Edge factory for intra-cluster edges. Returns the spatial edges in a cluster of features,
    organizing features by levels of nesting """

    def get_edges(self):
        """ Computes all nesting/overlapping edges for the provided interval dicts """
        if not self.dict1:
            for edge in self.edges:
                if edge:
                    yield edge
                else:
                    return None

        # sort the positions of all features associated to the current cluster
        positions = sorted([position for position in self.dict1])
        positions, dict_nesting_nested = NestingSorter(positions).remove_nested()

        # sequential or overlapping only:
        self.edges.update(self._get_inter_feature_edges(strand="+",
                                                        dict_interval=self.dict1,
                                                        list_positions=positions,
                                                        allowed_edgetypes={"cis-sequential", "cis-overlap"})
                          )

        # if there are nested features, compute their sequential/overlapping spatial relationships
        if dict_nesting_nested:
            self.edges.update(self._get_nested_inter_feature_edges(strand="+",
                                                                   dict_interval=self.dict1,
                                                                   dict_nesting_nested=dict_nesting_nested)
                              )

        if self.edges:
            for node1, node2, edgetype in self.edges:
                if self._passes_filter(node1, node2):
                    yield node1, node2, edgetype

        return None

    def _edges_from_interval_pair(self, strand=None,
                                  dict_interval=None,
                                  interval1=None,
                                  interval2=None,
                                  allowed_edgetypes=None):
        """Yields (feat1, feat2, edgetype) for all feature pairs at two position intervals."""

        # default values
        strand = initialize_parameter(strand, "+")
        dict_interval = initialize_parameter(dict_interval, {})
        interval1 = initialize_parameter(interval1, ())
        interval2 = initialize_parameter(interval2, ())
        allowed_edgetypes = initialize_parameter(allowed_edgetypes, set())
        feat1_set = dict_interval.get(interval1, set())
        feat2_set = dict_interval.get(interval2, set())

        if interval1 and interval2:
            start1, end1 = interval1
            start2, end2 = interval2

            # get edgetype for allowed_edgetypes edges and yield edges
            for feat1 in feat1_set:
                for feat2 in feat2_set:
                    node1, node2, edgetype = self.get_edge_and_type((start1, end1, feat1),
                                                                    (start2, end2, feat2),
                                                                    strand)
                    if edgetype in allowed_edgetypes:
                        yield node1, node2, edgetype

        return None

    def _get_nested_inter_feature_edges(self, strand=None, dict_interval=None, dict_nesting_nested=None):
        """ Computes edges between nested and overlapping features """

        # default values
        strand = initialize_parameter(strand, "+")
        dict_interval = initialize_parameter(dict_interval, {})
        dict_nesting_nested = initialize_parameter(dict_nesting_nested, {})
        inter_feature_edges = []

        # dispatch intervals to produce edges depending on nesting level

        for interval1, set_position2 in dict_nesting_nested.items():
            if interval1 in dict_interval:  # nesting interval
                new_edges_nest = set()

                for interval2 in set_position2:
                    if interval2 in dict_interval:  # nested interval

                        # produce the corresponding cis-nested edges between features
                        for edge in self._edges_from_interval_pair(strand=strand,
                                                                   dict_interval=dict_interval,
                                                                   interval1=interval1,
                                                                   interval2=interval2,
                                                                   allowed_edgetypes={"cis-nested"}):
                            new_edges_nest.add(edge)

                # additional edges (cis-sequential or cis-overlap) between nested features
                list_nested = sorted(set_position2)
                new_edges_seq_over = self._get_inter_feature_edges(strand=strand,
                                                                   dict_interval=dict_interval,
                                                                   list_positions=list_nested,
                                                                   allowed_edgetypes={"cis-sequential", "cis-overlap"})

                # add produced edges to the list of edges
                inter_feature_edges.append(new_edges_nest)
                inter_feature_edges.append(new_edges_seq_over)

        if inter_feature_edges:

            return set.union(*inter_feature_edges)

        return set()

    def _get_inter_feature_edges(self, strand=None, dict_interval=None, list_positions=None, allowed_edgetypes=None):
        """ Computes edges between consecutive features; allowed_edgetypes: set of str, allowed edge types
        (e.g. "cis-sequential", "cis-overlap"...) """

        # default values
        strand = initialize_parameter(strand, "+")
        dict_interval = initialize_parameter(dict_interval, {})
        list_positions = initialize_parameter(list_positions, [])
        allowed_edgetypes = initialize_parameter(allowed_edgetypes, set())
        new_edges = set()

        order = {"+": False, "-": True}

        # compute inter_feature edges from the sorted list of positions
        list_positions.sort(reverse=order[strand])
        for i, interval1 in enumerate(list_positions[:-1]):
            interval2 = list_positions[i + 1]
            if interval1 in dict_interval and interval2 in dict_interval:
                for edge in self._edges_from_interval_pair(strand=strand,
                                                           dict_interval=dict_interval,
                                                           interval1=interval1,
                                                           interval2=interval2,
                                                           allowed_edgetypes=allowed_edgetypes):
                    new_edges.add(edge)

        return new_edges

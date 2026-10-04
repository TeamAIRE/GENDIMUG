class TerminalEdges:
    """Builds terminal edges from pre-registered terminal features."""

    def map(self, terminal_anchors: dict, circular: bool) -> set:
        """ Builds edges from the anchor dict produced by TerminalFeatureMapper."""
        edges = set()
        for strand, anchors in terminal_anchors.items():
            if anchors:
                if circular:
                    edges.update(self._edges_circular(anchors))
                else:
                    edges.update(self._edges_linear(anchors))
        return edges

    @staticmethod
    def _edges_linear(anchors: dict) -> set:
        """ Returns a set of terminal edges for linear contigs """
        edges = set()
        if anchors["first_id"]:
            edges.add((anchors["first_id"], anchors["target"], "cis-sequential"))
        if anchors["last_id"]:
            edges.add((anchors["source"], anchors["last_id"], "cis-sequential"))
        return edges

    @staticmethod
    def _edges_circular(anchors: dict) -> set:
        """ Returns a set of terminal edges for circular contigs """
        edges = set()
        if anchors["only_id"]:
            edges.add((anchors["source"], anchors["only_id"], "cis-sequential"))
            edges.add((anchors["only_id"], anchors["target"], "cis-sequential"))
        return edges

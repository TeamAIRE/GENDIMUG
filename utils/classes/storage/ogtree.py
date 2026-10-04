from networkx import DiGraph, connected_components, bfs_tree
from utils.helpers.id_functions import get_prefix, reformat_ensembl_id
from utils.helpers.init_functions import initialize_dict


class OgTree:
    """
    Ontogenetic Parent-children tree, stored as networkx Directed Graphs, better than simple trees
    since they allow cycles.
    Attributes
    ----------
        tree: dict {id: nx DiGraph}
        child_to_parent_dict: dict {id: parent}, allowing to connect any feature to its parent gene
        genepart_to_gene: dict {genepart: original gene id}, allowing to connect any genepart feature to its parent
        trans-spliced gene
    """

    def __init__(self):
        self.tree = {}
        self.child_to_parent_dict = {}
        self.genepart_to_gene = {}

    def build(self, child_parent_dict, gene_genepart):
        """ Builds a tree from the child_parent_dict """
        # transfer the child_parent_dict to parent_tree attributes
        self.child_to_parent_dict = child_parent_dict
        # build the parent tree graph
        gff_edges = set((parent, child) for child, parent in self.child_to_parent_dict.items())
        genetic_graph = self.add_edges_to_network(gff_edges)
        self.populate(genetic_graph, gene_genepart)

    @staticmethod
    def add_edges_to_network(iterable_oriented_edges, network=None):
        """ Grows a networkx DiGraph (if provided) or initializes a new DiGraph with the iterable of tuples
         (node1, node2) """
        if not network:
            network = DiGraph()
        network.add_edges_from(iterable_oriented_edges)
        return network

    def add_edges(self, key, iterable_oriented_edges):
        """ Grows a networkx DiGraph (if provided) or initializes a new DiGraph with the iterable of tuples
         (node1, node2) """
        self.tree = initialize_dict(self.tree, key, DiGraph())
        self.tree[key].add_edges_from(iterable_oriented_edges)

    def populate(self, parent_tree_network, gene_genepart):
        """ Populate the genepart_gene and the parent_tree_dict using edges from the parent_tree_network
        networkx DiGraph and data from gene_genepart, to account for fragmented genes (trans-splicing exception) """
        for gene, set_genepart in gene_genepart.items():
            for genepart in set_genepart:
                self.genepart_to_gene[genepart] = gene
        for cc in connected_components(parent_tree_network.to_undirected()):
            genes = [e for e in cc if get_prefix(e) == "gene"]
            geneparts = [e for e in cc if get_prefix(e) == "genepart"]
            # for trans-spliced genes, correct with geneparts later
            if genes and not geneparts:
                for gene in genes:
                    self.tree[gene] = parent_tree_network.subgraph(cc).copy()
            if geneparts and not genes:
                gene = self.genepart_to_gene[geneparts[0]]
                self.tree[gene] = parent_tree_network.subgraph(cc).copy()

    def remove_edges(self, key, iterable_oriented_edges):
        """ Removes the oriented edges, as an iterable of tuples (node1, node2), from a provided networkx DiGraph """
        if self.tree.get(key, None):
            self.tree[key].remove_edges_from(iterable_oriented_edges)

    def remove_nodes(self, key, iterable_nodes):
        """ Removes the nodes, and associated edges, as an iterable of node ids, from a provided networkx DiGraph """
        if self.tree.get(key, None):
            self.tree[key].remove_nodes_from(iterable_nodes)
            for node in iterable_nodes:
                if node in self.child_to_parent_dict:
                    del self.child_to_parent_dict[node]

    def all_edges(self):
        """ Returns all directed edges from all gene trees, without attributes """
        all_edges = set()
        all_trees = set(self.tree.values())
        for tree in all_trees:
            all_edges.update(set(tree.edges()))
        return all_edges

    def get_ontogenetic_edges(self):
        """ Returns all ontogenetic edges """
        print("\n + ontogenetic edges...")
        ontogenetic_edges = set()
        tree_edges = self.all_edges()
        for (parent, child) in tree_edges:
            ontogenetic_edges.add((parent, child, "parent"))
        return ontogenetic_edges

    def all_nodes(self):
        """ Returns all nodes from all gene trees, without attributes """
        all_nodes = set()
        all_trees = set(self.tree.values())
        for tree in all_trees:
            all_nodes.update(set(tree.nodes()))
        return all_nodes

    def gene_tree_from(self, gene):
        """ Returns the networkx DiGraph gene tree associated to the provided gene """
        return self.tree.get(gene, DiGraph())

    def get_genes(self):
        """ Returns the set of keys of all gene trees """
        return self.tree.keys()

    def parent_genes(self, gene):
        """ Returns the set of gene or genepart ids at the root of the gene tree associated to the provided gene """
        return set(e for e in self.tree.get(gene, DiGraph()).nodes() if get_prefix(str(e)) in {"gene", "genepart"})

    def children(self, gene, set_prefixes=None):
        """ Returns the set of children ids with the provided prefix, if any, in the gene tree associated
        to the provided gene """
        if gene:
            if gene in self.tree:
                if not set_prefixes:
                    return set(e for e in self.tree.get(gene, DiGraph()).nodes() if e != gene)
                else:
                    return set(e for e in self.tree.get(gene, DiGraph()).nodes() if e != gene
                               and get_prefix(str(e)) in set_prefixes)
        return set()

    def children_from(self, gene, feature_id):
        """ Returns the set of all children ids for the provided feature id in the gene tree associated
        to the provided gene """
        return set(n for n in bfs_tree(self.tree.get(gene, DiGraph()), feature_id) if n != feature_id)

    def filtered_children_from(self, gene, feature_id, set_prefixes):
        """ Returns the set of all children ids for the provided feature id with the provided prefixes,
        in the gene tree associated to the provided gene """
        return set(e for e in self.children_from(gene, feature_id) if get_prefix(e) in set_prefixes)

    def update_child_parent_gene(self):
        """ Updates the child_to_parent_dict by iterating through the tree;
        Warning: if a feature has several parents in the graph (i.e. gene parts parenting an mRNA
        after trans-splicing...) those links will be missing has the child will only be connected to the parent
        as described in the GFF file """
        for gene, gene_tree in self.tree.items():
            for node in gene_tree.nodes():
                if node != gene and node not in self.child_to_parent_dict:
                    self.child_to_parent_dict[node] = gene

    def get_key_gene(self, feature_id):
        """ Returns the key gene_id for the provided feature_id """
        return self.child_to_parent_dict.get(feature_id, "")

    def iterate(self):
        """ Generates tuples (gene_id, gene_tree) for each gene_id in the tree """
        return self.tree.items()

    @staticmethod
    def update_child_parent_dict(child_parent_dict, feature_id, attributes, database="refseq"):
        """ Updates parent-child links in a dictionary {id_child: id_parent} """
        if attributes.get("Parent", ""):
            if database in {"refseq", "genbank"}:
                child_parent_dict[feature_id] = attributes.get("Parent", "")
            elif database == "ensembl":
                child_parent_dict[feature_id] = reformat_ensembl_id(attributes.get("Parent", ""))
        return child_parent_dict

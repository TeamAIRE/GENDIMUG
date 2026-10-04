from utils.classes.edges.overlap import Overlap
from utils.classes.edges.sequential import Sequential
from utils.helpers.init_functions import initialize_parameter
from utils.helpers.id_functions import get_prefix


def get_backbone_trans_edges(contig_storage, filter_ids=None, keep_ids=None):
    """
    Computes trans-overlapping/nesting and sequential edges for the current contig.
    Intergenic regions are excluded from the pool of features.
    :param contig_storage: ContigStorage object
    :param filter_ids: set(str)
    :param keep_ids: set(str)
    :return: set of edges
    """
    filter_ids = initialize_parameter(filter_ids, set())
    keep_ids = initialize_parameter(keep_ids, set())
    intergenic = contig_storage.feature_types.get_ids("intergenic")
    all_edges = set()
    all_edges = _backbone_trans_overlap_edges(contig_storage,
                                              all_edges,
                                              filter_ids=intergenic | filter_ids)

    all_edges = _backbone_trans_sequential_edges(contig_storage,
                                                 all_edges,
                                                 filter_ids=intergenic | filter_ids)
    contig_storage.edge_db.save_edges(all_edges)


def _backbone_trans_overlap_edges(contig_storage, all_edges, filter_ids=None, keep_ids=None):
    """
    Computes the trans (between strands) feature overlapping relationships for the backbone features.
    :param contig_storage: ContigStorage object
    :param all_edges: set of edges as tuples (node1, node2, edgetype)
    :param filter_ids: set(str)
    :param keep_ids: set(str)
    :return: updated set of edges
    """
    filter_ids = initialize_parameter(filter_ids, set())
    keep_ids = initialize_parameter(keep_ids, set())
    filter_ids = contig_storage.feature_types.get_ids("filtered", add_ids=filter_ids)

    # strand-specific interval dicts without the filtered features:
    plus_dict = contig_storage.positions.get_dict("+", filter_ids=filter_ids, keep_ids=keep_ids)
    minus_dict = contig_storage.positions.get_dict("-", filter_ids=filter_ids, keep_ids=keep_ids)
    if plus_dict and minus_dict:
        edges_trans_overlap = Overlap(dict1=plus_dict,
                                      strand1="+",
                                      dict2=minus_dict,
                                      strand2="-").get_edges()
        for edge in edges_trans_overlap:
            all_edges.add(edge)
    return all_edges


def _backbone_trans_sequential_edges(contig_storage, all_edges, filter_ids=None, keep_ids=None):
    """
    Computes the trans (between strands) feature sequential relationships for the backbone features.
    :param contig_storage: ContigStorage object
    :param all_edges: set of edges as tuples (node1, node2, edgetype)
    :param filter_ids: set(str)
    :param keep_ids: set(str)
    :return: updated set of edges
    """
    filter_ids = initialize_parameter(filter_ids, set())
    keep_ids = initialize_parameter(keep_ids, set())

    # strand-specific interval dicts with all features:
    plus_dict = contig_storage.positions.get_dict("+", filter_ids=filter_ids, keep_ids=keep_ids)
    minus_dict = contig_storage.positions.get_dict("-", filter_ids=filter_ids, keep_ids=keep_ids)
    if plus_dict and minus_dict:
        # Detecting trans-sequential links, excluding intergenic_region features
        edges_trans_sequential = Sequential(dict1=plus_dict,
                                            strand1="+",
                                            dict2=minus_dict,
                                            strand2="-").get_edges()
        for edge in edges_trans_sequential:
            all_edges.add(edge)
    return all_edges


def get_all_cis_overlapping_edges(contig_storage, filter_ids=None, keep_ids=None):
    """
    Computes cis-overlapping/nesting edges for all contigs and strands.
    :param contig_storage: ContigStorage object
    :param filter_ids: set(str), feature IDs to ignore
    :param keep_ids: set(str), feature IDs to keep
    return None
    """
    all_edges = set()
    filter_ids = contig_storage.feature_types.get_ids("filtered", add_ids=filter_ids)
    keep_ids = initialize_parameter(keep_ids, set())
    for strand in ["+", "-"]:
        interval_dict_all_features = contig_storage.positions.get_dict(strand, filter_ids=filter_ids)
        interval_dict_backbone_without_domains = contig_storage.positions.get_dict(strand,
                                                                                   filter_ids=filter_ids,
                                                                                   filter_prefixes={"dom",
                                                                                                    "intergenic",
                                                                                                    "terminal"})


        # compute overlaps
        backbone_overlap_edges = Overlap(dict1=interval_dict_all_features,
                                         strand1=strand,
                                         dict2=interval_dict_backbone_without_domains,
                                         strand2=strand).get_edges()
        for edge in backbone_overlap_edges:
            all_edges.add(edge)
        contig_storage.edge_db.save_edges(all_edges)


def get_additional_sequential_edges(global_parameters,
                                    contig_storage,
                                    genic_intergenic_mapping,
                                    filter_ids=None,
                                    keep_ids=None):
    """
    Computes all sequential edges for additional features (including those with nested genes).
    Those are generally computed whithin genic /intergenic clusters, except when features overlap multiple clusters.
    :param global_parameters: GlobalParameters object
    :param contig_storage: ContigStorage object
    :param genic_intergenic_mapping: dict {strand: {genic/intergenic_id: set(feature_ids)} }
    :param filter_ids: set(str)
    :param keep_ids: set(str)
    :return: None
    """
    all_edges = set()
    filtered_ids = contig_storage.feature_types.get_ids("filtered", add_ids=filter_ids)
    keep_ids = initialize_parameter(keep_ids, set())
    backbone = contig_storage.backbone_ids()

    for strand in ["+", "-"]:

        if strand not in genic_intergenic_mapping:
            continue

        # 1) one pass over all features for backbone feature to additional feature links

        edge_constraint = None
        interval_dict = contig_storage.positions.get_dict(strand, filter_ids=filtered_ids)#, keep_ids=keep_ids)

        if not global_parameters.genic_backbone:
            edge_constraint = _force_genic_nongenic_links  # This avoids excessive linkage between genic features,
            # bypassing the gene - intergenic links
            # not needed in gene-only networks

        sequential_edges = Sequential(dict1=interval_dict,
                                      strand1=strand).get_edges(edge_constraint=edge_constraint)
        for edge in sequential_edges:
            all_edges.add(edge)

        # 2) build a reverse lookup {feature_id -> cluster_key} for this strand, then compute additional-only edges
        #  and keep only intra-cluster ones

        feature_to_cluster = {}
        for key, set_ids in genic_intergenic_mapping[strand].items():
            all_add_ids = set_ids - backbone
            if len(all_add_ids) > 1:
                for add_id in all_add_ids:
                    feature_to_cluster.setdefault(add_id, set()).add(key)
        if not feature_to_cluster:
            continue

        # Single get_dict / Sequential call covering all additional features at once
        set_add_ids = set(feature_to_cluster)
        interval_dict = contig_storage.positions.get_dict(strand, keep_ids=set_add_ids)
        candidate_edges = Sequential(dict1=interval_dict, strand1=strand).get_edges()

        # Retain only edges whose two endpoints fall in the same cluster
        # (their cluster sets have a non-empty intersection)
        for edge in candidate_edges:
            if feature_to_cluster.get(edge[0], set()) & feature_to_cluster.get(edge[1], set()):
                all_edges.add(edge)

    contig_storage.edge_db.save_edges(all_edges)


def _force_genic_nongenic_links(edge):
    """
    Custom constraint imposed on provided edge, returns True if the constraint is fulfilled.
    Here the constraint is for the edges to be between genic and non-genic features.
    :param edge: tuple(node1: str, node2: str, edgetype: str)
    :return: bool
    """
    genic_prefixes = {"gene", "genepart", "rna", "exon", "intron", "cds", "cdspart", "dom"}
    node1, node2, edgetype = edge
    return ((get_prefix(node1) in genic_prefixes and get_prefix(node2) not in genic_prefixes) or
            (get_prefix(node2) in genic_prefixes and get_prefix(node1) not in genic_prefixes))

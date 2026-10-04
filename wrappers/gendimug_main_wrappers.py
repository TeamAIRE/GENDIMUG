# storage

from utils.classes.storage.cdsstorage import CdsStorage
from utils.classes.storage.contigstorage import ContigStorage

# sorting
from utils.classes.sorting.gffparser import GffParser
from utils.classes.sorting.gffsorter import GffSorter

# mappers
from utils.classes.mappers.genic import GenicMapper
from utils.classes.mappers.terminal import TerminalFeatureMapper
from utils.classes.mappers.intergenic import IntergenicMapper
from utils.classes.mappers.external import ExternalMapper
from utils.classes.mappers.introns import IntronMapper
from utils.classes.mappers.clusters import ClusterMapper

# edges
from utils.classes.edges.backbone import EdgeBackbone
from utils.classes.edges.exonintrons import EdgeIntrons
from utils.classes.edges.overlapone import BaseOverlap
from utils.classes.edges.backboneterminal import TerminalEdges
from utils.classes.edges.identical import IdenticalEdges
from utils.classes.edges.interdomain import InterdomainEdges

# export
from utils.classes.export.edgedatabase import EdgeDatabase

# helpers
from utils.helpers.time_functions import total_time
from utils.helpers.import_functions import import_genome_sequence
from utils.helpers.temporary_save import save_data, load_data, path_db_edges
from utils.helpers.export_functions import write_data, summary

# wrappers
from wrappers.ogtree_building import parental_links
from wrappers.external_import import (import_interproscan_data,
                                      import_repeatmodeler_data,
                                      import_tandemrepeatfinder_data,
                                      import_custom_gff)
from wrappers.complete_features import include_other_genomic_features
from wrappers.spatial_edges import (get_backbone_trans_edges,
                                    get_all_cis_overlapping_edges,
                                    get_additional_sequential_edges)


def build_network(global_parameters, start_time):
    """
    Network building main wrapper.
    :param global_parameters: GlobalParameters object
    :param start_time: time.time() object
    :return: None
    """
    # 1) loading genomic data : sequence and annotation

    # Initialization of whole genome storage structures
    counters = {"nodes": {}, "edges": {}}
    contig_description = {}
    cds_storage = CdsStorage()

    # User info
    global_parameters.print_paths()

    # Import sequence file
    fasta_dict, contig_names = import_genome_sequence(path_fasta=global_parameters.path_fna)
    if not global_parameters.only_contig:
        global_parameters.only_contig = list(contig_names)

    # Import GFF3 annotation file
    parser = GffParser(global_parameters)
    contig_data, contig_gene_fragments = parser.import_annotation(path_gff=global_parameters.path_gff)

    # 2) build networks based on genomic GFF data only, contig by contig, then store the created objects
    contig_store = {}  # only populated in the non-serialize path

    for contig in global_parameters.only_contig:

        fasta_seq, list_gff_data = _get_contig_data(contig,
                                                    fasta_dict,
                                                    contig_data)

        contig_storage, cds_storage, contig_description = _network_from_gff(global_parameters,
                                                                            cds_storage,
                                                                            contig,
                                                                            contig_description,
                                                                            list_gff_data,
                                                                            contig_gene_fragments)
        contig_storage.genome.sequence = fasta_seq

        if global_parameters.serialize:
            save_data(contig_storage, global_parameters.path_output_dir)
        else:
            contig_store[contig] = contig_storage

    # 3) optional global (all contigs) external imports
    if global_parameters.has_domains and global_parameters.keep_children:
        cds_storage = import_interproscan_data(path_ipr=global_parameters.path_ipr,
                                               cds_storage=cds_storage,
                                               database=global_parameters.database)

    int_repeats_storage, tandem_repeats_storage, list_additional_data = {}, {}, []
    if global_parameters.external:
        print(f"\n{'-' * 50}\nEXTERNAL IMPORTS")

        int_repeats_storage = import_repeatmodeler_data(path_rm=global_parameters.path_rm)
        tandem_repeats_storage = import_tandemrepeatfinder_data(path_trf=global_parameters.path_trf,
                                                                contig_names=contig_names)
        list_additional_data = import_custom_gff(path_custom=global_parameters.path_custom)

    # 4) enrich networks with external data and write output for each contig
    for contig in global_parameters.only_contig:
        print(f"\n{'-' * 50}\nCurrent contig: {contig}")

        # recover stored contig data
        if global_parameters.serialize:
            contig_storage = load_data(contig, global_parameters.path_output_dir)
        else:
            contig_storage = contig_store.pop(contig)

        # Include external data
        contig_storage = ExternalMapper(global_parameters,
                                        contig_storage,
                                        cds_storage,
                                        int_repeats_storage,
                                        tandem_repeats_storage,
                                        list_additional_data).map()

        # Include additional features if the network simplification mode allows it
        if global_parameters.external:
            contig_storage = include_other_genomic_features(global_parameters, contig_storage)

        print("\n REMAINING EDGES")
        # Build ontogenetic edges
        if global_parameters.has_ontogenetic_edges:
            ontogenetic_edges = contig_storage.og_tree.get_ontogenetic_edges()
            contig_storage.edge_db.save_edges(ontogenetic_edges)

        # Build remaining edges between additional features/backbone features
        _other_edges(global_parameters, contig_storage, cds_storage)
        
        # Export data and populate network summary counters
        counters["nodes"][contig], counters["edges"][contig] = write_data(global_parameters,
                                                                          contig_storage,
                                                                          contig_description)

    # Final clean up and export of counter data
    global_parameters.remove_temp_directories()
    summary(global_parameters, counters)
    total_time(start_time)


def _network_from_gff(global_parameters, cds_storage, contig, contig_description, list_gff_data, contig_gene_fragments):
    """
    Collects data and builds edges from genomic annotation file (GFF3 format)
    :param global_parameters: GlobalParameters object
    :param cds_storage: CdsStorage object, CDS-specific storage object
    :param contig: str, current contig ID
    :param contig_description: dict[str, list[GffData]], dict {contig ID: list of GffData objects}
    :param list_gff_data: list of GffData objects
    :param contig_gene_fragments: dict[str, dict[str, list[str]]], dict of dicts {contig ID: {gene: list geneparts}}
    :return: ContigStorage object, updated CdsStorage object, list of sets of edges as tuples (node1, node2, edgetype),
    updated contig_description
    """
    has_genic_backbone = global_parameters.genic_backbone
    has_spatial_edges = global_parameters.has_spatial_edges
    complete_network = global_parameters.complete_network
    keep_redundant = global_parameters.keep_redundant
    keep_children = global_parameters.keep_children

    print(f"\n{'-'*50}\nCurrent contig: {contig}")
    print("\nEDGES FROM GENOMIC ANNOTATION")

    # initialize a ContigStorage object
    contig_storage = ContigStorage()
    contig_storage.is_contig(contig)
    contig_storage.database = global_parameters.database

    contig_storage, cds_storage, child_parent_dict, contig_description = GffSorter(contig_storage,
                                                                                   cds_storage,
                                                                                   contig_description,
                                                                                   keep_children).sort(list_gff_data)
    db_path = path_db_edges(contig=contig,
                            path_output_dir=global_parameters.path_output_dir,
                            dir_name="temp_edges_db",
                            filename="edges")
    contig_storage.edge_db = EdgeDatabase(db_path)

    contig_storage.start_cutoff = global_parameters.start
    contig_storage.end_cutoff = global_parameters.end

    if not keep_redundant:

        dedup_db_path = path_db_edges(contig=contig,
                                      path_output_dir=global_parameters.path_output_dir,
                                      dir_name="temp_edges_db",
                                      filename="dedup_edges")
        contig_storage.dedup_edge_db = EdgeDatabase(dedup_db_path)

    gene_genepart = contig_gene_fragments.get(contig, {})

    contig_storage = parental_links(contig_storage, child_parent_dict, gene_genepart)

    contig_storage = GenicMapper(global_parameters, contig_storage, gene_genepart).map()

    edges = set()

    if has_spatial_edges:
        if not has_genic_backbone:  # backbone is made of genic and intergenic features

            print("\nInfering intergenic regions...")

            # mapping of strand-specific intergenic and terminal features

            contig_storage, terminal_anchors = TerminalFeatureMapper(global_parameters=global_parameters,
                                                                     contig_storage=contig_storage,
                                                                     contig_description=contig_description).map()
            contig_storage = IntergenicMapper(global_parameters=global_parameters, contig_storage=contig_storage).map()

            print("\n + (gene <> intergenic) edges...")
            print(" + (gene <> gene) edges...")

            # building backbone intergeinc/genic/terminal edges
            _, _, _, circular = contig_description[contig_storage.contig]
            if not global_parameters.start and not global_parameters.end:
                terminal_edges = TerminalEdges().map(terminal_anchors, circular)
            else:
                terminal_edges = set()
            edges = EdgeBackbone(contig_storage).get_edges() | terminal_edges

        else:

            # backbone is made of genic features only,
            # contig_description is needed to build edges closing the network for circular contigs
            print("\n + (gene <> gene) edges...")
            edges = EdgeBackbone(contig_storage, contig_description, gene=True).get_edges()

    if complete_network or not has_spatial_edges:

        print("\nInfering intronic regions...")
        contig_storage = IntronMapper(global_parameters, contig_storage, gene_genepart).map()

    if complete_network:

        intron_edges = EdgeIntrons(contig_storage, gene_genepart).get_edges()
        edges.update(intron_edges)
        print("\n + (exon <> intron) edges...")

        print("\nCDS fragmentation...")
        cds_part_edges = cds_storage.get_edges_between_cds_parts(contig_storage)
        edges.update(cds_part_edges)
        print("\n + (cds part <> cds part) edges...")

    print("\nDone.")

    contig_storage.edge_db.save_edges(edges)

    return contig_storage, cds_storage, contig_description


def _other_edges(global_parameters, contig_storage, cds_storage):
    """
    Wrapper managing the addition of additional spatial edges between genomic features and additional features
    :param global_parameters: GlobalParameters object
    :param contig_storage: ContigStorage object
    :param cds_storage: CdsStorage object, global storage of cds-specific data
    :return: None
    """
    has_spatial_edges = global_parameters.has_spatial_edges
    complete_network = global_parameters.complete_network
    has_domains = global_parameters.has_domains
    keep_children = global_parameters.keep_children
    genic_backbone = global_parameters.genic_backbone

    if has_domains and complete_network:

        print("\n + (protein motif <> protein motif) edges...")
        InterdomainEdges(contig_storage, cds_storage).get_edges()

    if has_spatial_edges:

        # Mapping of features by backbone clusters:
        # dict {strand: {genic or intergenic ID defining a cluster: set of overlapping features}}
        cluster_mapping = ClusterMapper(global_parameters, contig_storage).map()
        # add trans_overlap and trans_sequential edges between features in different strands
        # (using a filter on intergenic ids to reduce the number of edges)
        print("\n + trans-strand backbone edges...")

        # to create simplified networks, children features can be filtered out (keep_children == False)
        # and genic ids can be kept (genic_backbone == True)

        children = set()
        keep_ids = set()
        if not keep_children:
            children = contig_storage.feature_types.get_ids("children")
            # remove trans-spliced geneparts from the set of children to remove
            geneparts = {feature_id for feature_id in children if feature_id.startswith("genepart")}
            children = children - geneparts
        if genic_backbone:
            keep_ids = contig_storage.feature_types.get_ids("genic")
        get_backbone_trans_edges(contig_storage,
                                 filter_ids=children,
                                 keep_ids=keep_ids)

        # get all cis-nesting/overlapping edges between backbone + additional features
        # (except between domains as their edges are based on proteic positions)
        print("\n + cis-strand overlap edges...")
        get_all_cis_overlapping_edges(contig_storage,
                                      filter_ids=children,
                                      keep_ids=keep_ids)
        # (using a filter on intergenic ids to reduce the number of edges)
        print("\n + additional sequential edges...")

        intergenic = contig_storage.feature_types.get_ids("intergenic")
        get_additional_sequential_edges(global_parameters,
                                        contig_storage,
                                        cluster_mapping,
                                        filter_ids=intergenic | children,
                                        keep_ids=keep_ids)

        # detect remaining one base overlaps

        edges_one_base_overlap = BaseOverlap(contig_storage,
                                             filter_ids=children,
                                             keep_ids=keep_ids).get_edges()

        contig_storage.edge_db.save_edges(edges_one_base_overlap)
        # manage edges between identical features
        contig_storage = IdenticalEdges(global_parameters, contig_storage).get_edges()

        if not keep_children:
            contig_storage.gff.remove_set(children)

    else:
        # in ontogenetic-only network, no collapsing of redundant features, the edge db path is transferred
        contig_storage.dedup_edge_db.db_path = contig_storage.edge_db.db_path


def _get_contig_data(contig, fasta_dict, contig_data):
    """Uniform interface to retrieve (fasta_seq, list_gff_data) regardless of serialize mode."""
    return fasta_dict[contig], contig_data[contig]

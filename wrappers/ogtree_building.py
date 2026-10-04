from utils.classes.storage.gffdata import GffData
from utils.helpers.id_functions import get_prefix, get_suffix, reorder_parts
from utils.helpers.interval_functions import is_match


def parental_links(contig_storage, child_parent_dict, gene_genepart):
    """
    Wrapper that computes the parent_child links and populate their associated dictionaries in the ContigStorage object.
    :param contig_storage: ContigStorage object
    :param child_parent_dict: dict {id_child: id_parent}
    :param gene_genepart: dict {gene ID: list of geneparts}
    :return: updated contig_storage
    """
    # populate the parent-children tree dict
    contig_storage = _build_og_tree(contig_storage, child_parent_dict, gene_genepart)
    # infer the exon to cds parts parent-child links from genomic positions
    contig_storage = _infer_parent_links_exon_cds(contig_storage)
    # identify genes with no child
    contig_storage = _find_childless_genes(contig_storage)
    # build gene_centric {parent gene <-> children} dicts, as genes form the backbone used to:
    # - define intergenic features
    # - cluster intragenic features
    # backbone genes will be all roots (prefixes 'gene' and 'genepart') from the parent tree
    return contig_storage


def _build_og_tree(contig_storage, child_parent_dict, gene_genepart):
    """
    Builds a tree from the parent_links dict
    :param contig_storage: ContigStorage object
    :param child_parent_dict: dict {id_child: id_parent}
    :param gene_genepart: dict {gene ID: list of geneparts}
    :return: updated contig_storage with parent_tree: dict {parent_type: {parent_id: {child_type: set(child_id)}}}
     and updated feature_sets with categories "additional_features"
    """
    # build the parent tree graph
    contig_storage.og_tree.build(child_parent_dict, gene_genepart)
    # modify graphs for trans-spliced orfs:
    contig_storage = _correct_trans_spliced_parent_links(contig_storage, gene_genepart)
    # modify trans-spliced rna positions:
    contig_storage = _get_trans_spliced_rna_positions(contig_storage, gene_genepart)
    return contig_storage


def _correct_trans_spliced_parent_links(contig_storage, gene_genepart):
    """
    Modifies the parental links between trans-spliced genes and their direct children (rna or exon),
    by replacing the gene by gene parts.
    :param contig_storage: ContigStorage object
    :param gene_genepart: dict {gene: list(geneparts)}
    :return: updated og_tree
    """
    for gene, geneparts in gene_genepart.items():
        new_edges = set()
        removed = set()
        edges = contig_storage.og_tree.gene_tree_from(gene)
        for parent, child in list(edges.edges()):
            if parent == gene:
                removed.add((parent, child))
                if get_prefix(child) == "exon":
                    suffix = get_suffix(child)
                    parent = parent.replace("gene", "genepart") + f"-{suffix}"
                    new_edges.add((parent, child))
                elif get_prefix(child) in {"rna"}:
                    for genepart in geneparts:
                        parent = genepart
                        new_edges.add((parent, child))
        contig_storage.og_tree.add_edges(gene, new_edges)
        contig_storage.og_tree.remove_edges(gene, removed)
    return contig_storage


def _get_trans_spliced_rna_positions(contig_storage, gene_genepart):
    """
    Corrects the trans-spliced RNA positions to a list of genepart positions.
    :param contig_storage: ContigStorage object
    :param gene_genepart: dict {gene: list(geneparts)}
    :return: updated contig_storage
    """
    for gene, list_geneparts in gene_genepart.items():
        rnas = contig_storage.og_tree.children(gene, set_prefixes={"rna"})
        for rna in rnas:
            nested = False
            db, feature, old_positions, score, list_phases, attributes = contig_storage.gff.data(rna).data_tuple
            list_geneparts = reorder_parts(list_geneparts)
            genepart_positions = [contig_storage.gff.list_positions(genepart)[0] for genepart in list_geneparts]
            genepart_phases = [contig_storage.gff.list_phases(genepart)[0] for genepart in list_geneparts]
            contig, strand, start, stop = old_positions[0]
            for genepart_position in genepart_positions:
                new_contig, new_strand, new_start, new_stop = genepart_position
                if contig == new_contig and (strand == new_strand or strand == "?"):
                    if new_start <= start <= stop <= new_stop:  # nesting
                        nested = True
                        break
            if not nested:
                # RNA takes the positions  of geneparts
                list_positions = genepart_positions.copy()
                list_phases = genepart_phases
            else:
                # RNA is nested into a genepart (e.g. ncRNA) so they keep their original position
                list_positions = old_positions.copy()
            data_tuple = db, feature, list_positions, score, list_phases, attributes
            contig_storage.gff.update(GffData((rna, data_tuple), database=contig_storage.database, is_tuple=True),
                                      replace=True)
    return contig_storage


def _infer_parent_links_exon_cds(contig_storage):
    """
    Computes parent-child links between exons and cds/cdsparts
    :param contig_storage: ContigStorage object
    :return: updated contig_storage
    """
    # GFF3 (eukaryotic) files from the ncbi have parental links from rnas to both exons and cds
    for gene in contig_storage.og_tree.get_genes():
        gene_tree = contig_storage.og_tree.gene_tree_from(gene)
        rnas = set(n for n in gene_tree.nodes() if get_prefix(n) in {"rna"})
        if rnas:
            for rna in rnas:  # prokaryotic genomes will not have these links
                rna_exons = set((n1, n2) for (n1, n2) in gene_tree.edges() if n1 == rna and get_prefix(n2) == "exon")
                rna_cds = set((n1, n2) for (n1, n2) in gene_tree.edges()
                              if n1 == rna and get_prefix(n2) in {"cds", "cdspart"})
                contig_storage.og_tree.remove_edges(gene, rna_cds)
                list_cds = reorder_parts([tup[1] for tup in rna_cds])
                list_exons = reorder_parts([tup[1] for tup in rna_exons])
                if list_cds:
                    new_links = []
                    list_intervals_exons = [(start, stop, exon)
                                            for exon in list_exons
                                            for _, _, start, stop in contig_storage.gff.list_positions(exon)
                                            ]
                    list_intervals_cds = [(start, stop, cds)
                                          for cds in list_cds
                                          for _, _, start, stop in contig_storage.gff.list_positions(cds)
                                          ]
                    if len(list_intervals_exons) == 1 and len(list_intervals_cds) == 1:
                        new_links.append((list_intervals_exons[0][2], list_intervals_cds[0][2]))
                    else:
                        for (start1, stop1, exon) in list_intervals_exons:
                            for (start2, stop2, cds) in list_intervals_cds:
                                if is_match((start1, stop1), (start2, stop2)):
                                    new_links.append((exon, cds))
                    contig_storage.og_tree.add_edges(gene, new_links)
    return contig_storage


def _find_childless_genes(contig_storage):
    """
    Identifies the gene ids that are not linked to any child and stores the set of these ids
    in contig_storage.feature_types['childless_genes'].
    :param contig_storage: ContigStorage object
    :return: updated contig_storage
    """
    all_parent_genes = {gene_id for gene_id in contig_storage.og_tree.get_genes() if gene_id.startswith("gene-")}
    all_genes = {gene_id for gene_id in contig_storage.gff.ids() if (contig_storage.gff.feature(gene_id) == "gene"
                                                                     and get_prefix(gene_id) != "genepart")}
    childless_genes = all_genes - all_parent_genes
    contig_storage.feature_types.update("childless_genes", childless_genes)
    return contig_storage

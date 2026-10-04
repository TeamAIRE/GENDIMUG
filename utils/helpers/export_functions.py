from os.path import join
from utils.classes.export.edgeattributes import EdgeAttributes
from utils.helpers.id_functions import get_prefix
from utils.helpers.interval_functions import get_length_from_overlapping_positions
from utils.helpers.os_sys_functions import check_output_dir


def summary(global_parameters, counters):
    """
    Saves a summary of the produced network statistics
    :param global_parameters: GlobalParameters Object
    :param counters: dict {"edges": {contig: int}, "nodes": {contig: int}}
    :return: None
    """
    path_output_dir = global_parameters.path_output_dir
    path_summary = join(path_output_dir, f"summary_{global_parameters.root}.tsv")
    counter_edges = counters["edges"]
    counter_nodes = counters["nodes"]
    total_nodes = sum(list(counter_nodes.values()))
    total_edges = sum(list(counter_edges.values()))
    with open(path_summary, "w") as out:
        out.write("contig\tnodes\tedges\n")
        for contig in sorted(counter_edges.keys()):
            nb_nodes = counter_nodes[contig]
            nb_edges = counter_edges[contig]
            out.write(f"{contig}\t{nb_nodes}\t{nb_edges}\n")
        out.write(f"\ntotal\t{total_nodes}\t{total_edges}\n")


def capture_unmapped_domains(contig_storage, path_output_dir, unmapped_domains):
    """
    Optional export of unmapped protein motifs present in the interproscan file, if any.
    :param contig_storage: ContigStorage object
    :param path_output_dir: path to output directory
    :param unmapped_domains: set of unmapped domain IDs
    :return: None
    """
    contig = contig_storage.contig
    out_dir = join(path_output_dir, "unmapped_domains")
    check_output_dir(out_dir)
    filename = f"unmapped_domains_for_{contig}.txt"
    path = join(out_dir, filename)
    with open(path, "w") as out:
        for dom in sorted(unmapped_domains):
            out.write(f"{dom}\n")


def _capture_domains_non_multiple_of_three(contig_storage, path_output_dir, inconsistent_domains):
    """
    Optional export of unmapped protein motifs present in the interproscan file.
    :param contig_storage: ContigStorage object
    :param path_output_dir: path to output directory
    :param inconsistent_domains: set of domain IDs whose positions don't add up to a length multiple of 3
    :return: None
    """
    contig = contig_storage.contig
    out_dir = join(path_output_dir, "domains_of_length_not_multiple_of_3")
    check_output_dir(out_dir)
    filename = f"inconsistent_domains_mapped_for_{contig}.tsv"
    path = join(out_dir, filename)
    with open(path, "w") as out:
        out.write("motif\tlist of positions\ttotal length (bp)\tknown exception\n")
        for data in sorted(inconsistent_domains):
            dom, positions, length, exception = data
            out.write(f"{dom}\t{positions}\t{length}\t{exception}\n")


def write_data(global_parameters, contig_storage, contig_description):
    """
    Wrapper to export edge and node data to output files.
    :param global_parameters: GlobalParameters Object
    :param contig_storage: ContigStorage object
    :param contig_description: dict {contig ID: contig data tuple}
    :return: int, int (number of nodes, number of edges)
    """
    if global_parameters.encode:
        contig_storage.node_numeric_index()
    nb_edges, nodes = _write_network(global_parameters, contig_storage, contig_description)
    nb_nodes = _write_annotation(global_parameters, contig_storage, nodes)
    return nb_nodes, nb_edges


def _write_network(global_parameters, contig_storage, contig_description):
    """
    Saves edge-associated data in a cytoscape-readable csv format. Adds the length of each feature (feature_length)
    :param global_parameters: GlobalParameters Object
    :param contig_storage: ContigStorage object
    :param contig_description: dict {contig ID: contig data tuple}
    :return: int, the total number of edges and set(nodes)
    """
    total = 0
    col_names = {True: "node1,node2,et_index,ec_index,overlap_length,distance\n",
                 False: "node1,node2,edgetype,edgeclass,overlap_length,distance\n"}
    contig = contig_storage.contig
    with open(global_parameters.get_path_network(contig), "w") as f:
        f.write(col_names[global_parameters.encode])
        # add edge attributes
        set_edges = set(contig_storage.edge_db.load_edges() if global_parameters.keep_redundant
                        else contig_storage.dedup_edge_db.load_edges())
        nodes = _get_all_nodes(set_edges)
        set_edges = EdgeAttributes(set_edges).add_attributes(contig_storage,
                                                             contig_description,
                                                             encode=global_parameters.encode)
        total = len(set_edges)
        print(f"\nSaving {total} edges...")
        for edge in set_edges:
            f.write(edge)
    print("Done.")
    return total, nodes


def _get_all_nodes(set_edges):
    """ Returns a set of nodes present in edges """
    set_nodes = set()
    for n1, n2, _ in set_edges:
        set_nodes.update({n1, n2})
    return set_nodes


def _write_annotation(global_parameters, contig_storage, set_nodes):
    """
    Saves node-associated data in a cytoscape-readable csv format. Adds the length of each feature (feature_length)
    and the % GC if a fasta file was provided.
    :param global_parameters: GlobalParameters Object
    :param contig_storage: ContigStorage object
    :return: int, the number of annotations written to file
    """
    chunksize = 10000
    contig = contig_storage.contig
    inconsistent_domains = set()
    # lines to be written to files are written in chunks of size 'chunksize'
    # to reduce the I/O load and speed up the process

    path_output_annot = global_parameters.get_path_node_annotation(contig)
    categories = _compile_categories(contig_storage)
    with open(path_output_annot, "w") as f:
        col_names = {True: f"node,feature_id,feature,positions,phases,{','.join(categories)}\n",
                     False: f"node,feature,positions,phases,{','.join(categories)}\n"}
        f.write(col_names[global_parameters.encode])

        c = 0
        total = 0
        text = set()
        for feature_id in set_nodes:
            db, feature, list_positions, score, list_phases, attributes = contig_storage.gff.data_tuple(feature_id)
            if feature != "region":  # contigs are not nodes in the network
                c += 1
                total += 1

                # create a str that lists all genomic segments in order, separated by pipes
                positions = contig_storage.gff.convert_list_positions_to_string(feature_id)
                phase = contig_storage.gff.convert_list_phases_to_string(feature_id)
                attributes["feature_length"] = get_length_from_overlapping_positions(list_positions)
                attributes["percent_GC"] = contig_storage.genome.get_gc(list_positions)
                line_start = f"{feature_id},{feature},{positions},{phase}"
                if global_parameters.encode:
                    alias = contig_storage.num_index[feature_id]
                    line_start = f"{alias},{line_start}"

                # add values in categories order
                vals = []
                for cat in categories:
                    vals.append(attributes.get(cat, "none"))
                vals = [str(e).replace(",", "|") for e in vals]
                line = f"{line_start},{','.join(vals)}\n"
                text.add(line)

                # test for the size of the list of lines relative to chunksize, if it is of the chunksize, write to file
                # and reset the list
                if c % chunksize == 0:
                    f.write("".join(text))
                    c = 0
                    text = set()

            # this is a test for protein encoding regions, whose length is expected to be multiple of 3 (codons)
            # this is sometimes not the case, due to e.g. ribosomal slippage or low quality sequence.
            # such inconsistent domains are captured below:
            if get_prefix(feature_id) == "dom":
                if attributes["feature_length"] % 3:
                    exception = attributes.get("exception", "")
                    # if different from 0, meaning it's not multiple of 3
                    inconsistent_domains.add((feature_id, positions, attributes["feature_length"], exception))

        print(f"\nSaving {total} node annotations...")
        # final chunk is written to file
        f.write("".join(text))

        # encoding of edge attributes in separate index files
        if global_parameters.encode:
            _export_edge_encoding(global_parameters, contig_storage, contig)

        # signal probable error regarding cutoff parameters
        if global_parameters.start and global_parameters.end:
            if global_parameters.start >= global_parameters.end:
                print("\nWARNING: the chosen start and end cutoff parameters forbid network construction.\n")

        print("Done.")
    # Export detected inconsistent domains
    if inconsistent_domains:
        _capture_domains_non_multiple_of_three(contig_storage, global_parameters.path_output_dir, inconsistent_domains)

    return total


def _compile_categories(contig_storage):
    """
    Compiles all the keys of all attributes dicts in gff_dict, to be used as annotation column names
    :param contig_storage: ContigStorage object
    :return: sorted deduplicated list of categories
    """
    categories = {"feature_length", "percent_GC"}
    for data in contig_storage.gff.all_data():
        categories.update(set(data.attributes))
    return sorted(categories)


def _export_edge_encoding(global_parameters, contig_storage, contig):
    """
    Writes the edge attributes index code in output files
    :param global_parameters: GlobalParameters object
    :param contig_storage: ContigStorage object
    :param contig: contig ID, str
    :return: None
    """
    path_output_edgetype_index = global_parameters.get_path_edgetype_index(contig)
    with open(path_output_edgetype_index, "w") as f:
        f.write("et_index,edgetype\n")
        for k, v in contig_storage.edgetype_index.items():
            f.write(f"{v},{k}\n")
    path_output_edgeclass_index = global_parameters.get_path_edgeclass_index(contig)
    with open(path_output_edgeclass_index, "w") as f:
        f.write("ec_index,edgeclass\n")
        for k, v in contig_storage.edgeclass_index.items():
            f.write(f"{v},{k}\n")

from utils.helpers.id_functions import get_prefix
from utils.helpers.export_functions import capture_unmapped_domains


def include_other_genomic_features(global_parameters, contig_storage):
    """
    Wrapper controlling the inclusion of additional gff features not used in the backbone to the network
    :param global_parameters: GlobalParameters object
    :param contig_storage: ContigStorage object
    :return: updated contig_storage
    """
    all_features = {feature_id for feature_id in contig_storage.gff.ids()
                    if contig_storage.gff.feature(feature_id) != "region"}

    mapped_features = contig_storage.og_tree.all_nodes() | contig_storage.feature_types.get_ids("genic",
                                                                                                "intergenic",
                                                                                                "repeats",
                                                                                                "filtered")
    additional_features = contig_storage.feature_types.get_ids("additional_features")

    # identify and export domains that were not mapped (usually indicating position discrepancies
    # between the protein and genomic annotation files)

    unmapped_domains = _any_unmapped_domain(global_parameters, contig_storage, all_features, mapped_features)
    contig_storage.feature_types.update("filtered", unmapped_domains)

    # identify remaining features and add them to the ContigStorage.IntervalData object
    # exclude region/chromosome = contigs (they would be linked to every feature) and unmapped domains

    filtered_other_features = {feature_id for feature_id in additional_features
                               if contig_storage.gff.feature(feature_id) not in {"region", "chromosome"}}

    remaining_features = filtered_other_features - unmapped_domains

    # add any remaining feature to the ContigStorage.IntervalData object
    if remaining_features:
        for feature_id in remaining_features:
            contig_storage.positions.update(feature_id, contig_storage.gff.list_positions(feature_id))

    return contig_storage


def _any_unmapped_domain(global_parameters, contig_storage, all_features, mapped_features):
    """
    Exports the list of protein domains which were not mapped, usually because of discrepancies between protein
     (amino-acid) and CDS (base pairs) sequences.
    :param global_parameters: GlobalParameters object
    :param contig_storage: ContigStorage object
    :param all_features: set of all features
    :param mapped_features: set of mapped features
    :return: None
    """
    unmapped_domains = set(feature_id for feature_id in all_features
                           if (feature_id not in mapped_features and get_prefix(feature_id) == "dom"))

    if unmapped_domains:
        capture_unmapped_domains(contig_storage=contig_storage,
                                 path_output_dir=global_parameters.path_output_dir,
                                 unmapped_domains=unmapped_domains)
    return unmapped_domains

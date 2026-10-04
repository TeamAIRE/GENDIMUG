from os.path import join, exists
from utils.helpers.init_functions import initialize_contig_dict


class Mapper:
    """ Base mapper class"""
    def __init__(self, global_parameters, contig_storage):

        # transfer of GlobalParameters attributes
        self.generate_gff = global_parameters.generate_gff
        self.path_gff_dir = global_parameters.path_gff_dir
        self.root = global_parameters.root
        self.keep_redundant = global_parameters.keep_redundant
        self.has_domains = global_parameters.has_domains
        self.external = global_parameters.external
        self.complete_network = global_parameters.complete_network
        self.keep_children = global_parameters.keep_children
        self.genic_backbone = global_parameters.genic_backbone
        self.has_spatial_edges = global_parameters.has_spatial_edges

        # transfer of ContigStorage attributes
        self.contig_storage = contig_storage
        self.contig = contig_storage.contig
        self.database = contig_storage.database
        self.start_cutoff = contig_storage.start_cutoff
        self.end_cutoff = contig_storage.end_cutoff

    def _positions_fit(self, list_positions):
        """ Tests whether a list of positions fits into the start and end position cutoffs"""
        res = []
        if not self.end_cutoff:
            for _, _, start, end in list_positions:
                res.append(start >= self.start_cutoff)
            return all(res)
        else:
            for _, _, start, end in list_positions:
                res.append(start >= self.start_cutoff and end <= self.end_cutoff)
            return all(res)

    def _export_gff(self, prefix, set_ids):
        """ Optional export of intronic features in the GFF3 format"""
        if self.generate_gff:
            gff_path = join(self.path_gff_dir, f"{prefix}__{self.root}.gff")
            self._write_gff(path_output_gff=gff_path, set_ids=set_ids)

    def _write_gff(self, path_output_gff: str, set_ids: set[str]) -> None:
        """
        Saves data to an output file in the gff3 format, ordered by contig, start and end.
        :param path_output_gff: path to the output gff3 file
        :param set_ids: set of ids (str)
        :return: None
        """
        reordering = {}
        # lines to be written to files are written in chunks of size 'chunksize' to reduce the I/O load
        # and speed up the process
        chunksize = 10000
        text = ""

        for feature_id in set_ids:
            if feature_id in self.contig_storage.gff.ids():
                db, feature, list_positions, score, phases, attributes = self.contig_storage.gff.data_tuple(feature_id)
                attributes["ID"] = feature_id
                for index, (contig, strand, start, end) in enumerate(list_positions):
                    if len(list_positions) > 1:
                        attributes["part"] = index + 1
                    current_attributes = ";".join([f"{k}={v}" for k, v in attributes.items()])
                    reordering = initialize_contig_dict(reordering, contig, start, {})
                    if end not in reordering[contig][start]:
                        reordering[contig][start][end] = set()
                    phase = phases[index]
                    reordering[contig][start][end].add(f"{contig}\t{db}\t{feature}\t{start}\t{end}"
                                                       f"\t{score}\t{strand}\t{phase}\t{current_attributes}\n")

        already_started = False
        if exists(path_output_gff):
            already_started = True

        with open(path_output_gff, "a") as f:
            if not already_started:
                f.write("##gff-version 3\n")
            c = 0
            for contig in sorted(reordering):
                for start in sorted(reordering[contig]):
                    for end in sorted(reordering[contig][start]):
                        data_set = reordering[contig][start][end]
                        for data in data_set:
                            c += 1
                            text += data
                        if c % chunksize == 0:
                            f.write(text)
                            c = 0
                            text = ""
            f.write(text)

from utils.classes.mappers.mapperbase import Mapper
from utils.classes.storage.gffdata import GffData
from utils.classes.sorting.nestingsorter import NestingSorter
from utils.helpers.update_functions import reverse_dict_of_sets


class ExternalMapper(Mapper):

    def __init__(self, global_parameters, contig_storage, cds_storage, int_repeats, tandem_repeats, list_data):
        super().__init__(global_parameters=global_parameters, contig_storage=contig_storage)
        self.cds_storage = cds_storage
        self.int_repeats = int_repeats
        self.tandem_repeats = tandem_repeats
        self.list_data = list_data

    def map(self):
        if self.has_domains and (self.complete_network or not self.has_spatial_edges):
            # only networks using the domain data
            self.map_interproscan()
        if self.external:  # only networks using the additional features
            self.map_interspersed_repeats()
            self.map_tandem_repeats()
            self.map_custom_features()
        return self.contig_storage

    def map_interproscan(self):
        """ Includes interproscan data into the network """

        # compute the genomic positions of domains based on those of the cds
        # update the Gff object with these positions and OgTree object with matching cds parts and domains

        self.contig_storage = self.cds_storage.map_domains(self.contig_storage)

        # update the current positions dict with all proteic domains
        all_domains = self.contig_storage.feature_types.get_ids("domains")
        for domain in all_domains:
            self.contig_storage.positions.update(domain, self.contig_storage.gff.list_positions(domain))

        # generate an optional GFF3 file for domains
        self._export_gff(prefix="proteic", set_ids=all_domains)

        return self.contig_storage

    def map_interspersed_repeats(self):
        """ Includes repeatmodeler data into the network """

        # update ContigStorage objects for interspersed repeats
        if self.contig in self.int_repeats and self.int_repeats[self.contig]:
            for transp_id, data in self.int_repeats[self.contig].items():
                _, _, list_positions, _, _, _ = data
                if self._positions_fit(list_positions):
                    self.contig_storage.gff.update(GffData((transp_id, data),
                                                           database=self.database,
                                                           is_tuple=True))
                    self.contig_storage.positions.update(transp_id, list_positions)
                    self.contig_storage.feature_types.update("repeats", {transp_id})

            # generate an optional GFF3 file for interspersed repeats
            self._export_gff(prefix="interspersed_repeats_RepeatModeler", set_ids=set(self.int_repeats[self.contig]))

        return self.contig_storage

    def map_tandem_repeats(self):
        """ Includes Tandem Repeat Finder (TRF) data in the network """

        # update ContigStorage objects for tandem repeats
        if self.contig in self.tandem_repeats and self.tandem_repeats[self.contig]:
            for repeat_id, data in self.tandem_repeats[self.contig].items():
                if repeat_id not in self.contig_storage.gff.ids():
                    _, _, list_positions, _, _, _ = data
                    if self._positions_fit(list_positions):
                        self.contig_storage.gff.update(
                            GffData((repeat_id, data), database=self.database, is_tuple=True))
                        self.contig_storage.positions.update(repeat_id, list_positions)
                        self.contig_storage.feature_types.update("repeats", {repeat_id})

            # generate an optional GFF3 file for tandem repeats
            self._export_gff(prefix="tandem_repeats_TRF_TRAP", set_ids=set(self.tandem_repeats[self.contig]))

        # optional removal of nested repeats        
        if not self.keep_redundant:
            self._remove_nested_satellites()

        return self.contig_storage

    def _remove_nested_satellites(self):
        """ Removes short motif repeats nested into larger motif repeats """
        interval_dict_satellites = self.contig_storage.positions.get_dict("+", keep_prefixes={"sat"})
        list_satellites = self.contig_storage.positions.get_list("+", keep_prefixes={"sat"})
        _, nesting_dict = NestingSorter(list_satellites).remove_nested()
        nested_nesting_dict = reverse_dict_of_sets(nesting_dict)
        for start, end in nested_nesting_dict:
            list_positions = [(self.contig, "+", start, end)]
            nested_satellites = interval_dict_satellites[(start, end)]
            self.contig_storage.feature_types.update("filtered", nested_satellites)
            self.contig_storage.gff.remove_set(self.contig_storage.feature_types.get_ids("filtered"))
            self.contig_storage.positions.remove_set(nested_satellites, list_positions)

    def map_custom_features(self):
        """ Includes custom gff data as additional features """

        # update ContigStorage objects for custom features
        for data_dict in self.list_data:
            data = data_dict.get(self.contig, [])
            if data:
                for data_tuple in data:
                    _, _, list_positions, _, _, attributes = data_tuple
                    if self._positions_fit(list_positions):
                        feature_id = attributes.get("ID", "")
                        if feature_id:
                            self.contig_storage.gff.update(GffData((feature_id, data_tuple),
                                                                   database=self.database,
                                                                   is_tuple=True))
                            self.contig_storage.positions.update(feature_id, list_positions)
                            self.contig_storage.feature_types.update("additional_features", {feature_id})

        return self.contig_storage

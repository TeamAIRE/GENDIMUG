from utils.helpers.gff_functions import format_gff_line
from utils.helpers.id_functions import reformat_ensembl_id


class GffData:
    """ Stores the gff data from the parsing of a GFF line or a tuple.
    :param args: arguments; str if tuple == False; (str, Tuple) if tuple == True
    :param is_tuple: bool, if True, args is a tuple of GFF3 data, arguments should be a tuple (feature_id, data_tuple);
    if False (default), args is the line read from the GFF3 input file corresponding to the description
    of a feature (str) """

    def __init__(self, args, database, is_tuple=False):
        if is_tuple:
            self.id, self.data_tuple = args
        else:
            self.id, self.data_tuple = self.parse_gff_line(args, database)
        self.db, self.feature, self.positions, self.score, self.phase, self.attributes = self.data_tuple

    def prefix(self):
        """ Returns the prefix of the feature id, or '' """
        for sep in ["-", ":", "_"]:
            if sep in self.id:
                return self.id.split(sep)[0].lower()
        return ""

    def suffix(self):
        """ Returns the suffix of the feature id, or '' """
        if "-" in self.id:
            suffix = self.id.split("-")[-1]
            if suffix.isnumeric():
                return int(suffix)
        return ""

    def parse_gff_line(self, line, database):
        """ Returns a tuple (feature ID, data tuple) from a GFF3 line """
        contig, db, feature, start, end, score, strand, phase, attributes = format_gff_line(line)
        # only import features with an ID tag
        feature_id = self.reformat_feature_id(feature, attributes, database)
        if not feature_id:
            return "", ("", "", [], "", [], {})
        # refseq/genbank trans-splicing annotations
        if feature_id and feature == "gene" and "exception" in attributes and "part" in attributes:
            # trans-splicing flag in refseq/genbank gff files:
            # detects trans-spliced genes with the attribute 'exception: trans-splicing',
            # and renames the gene parts with the prefix 'genepart' and the part number as suffix '-n'.
            if attributes["exception"] == "trans-splicing" and not feature_id.startswith("genepart-"):
                attributes["original name"] = feature_id
                feature_id = f"{feature_id.replace('gene-', 'genepart-')}-{attributes['part']}"
                attributes["ID"] = feature_id
            # convert contig, strand, start, end information into a list of positions
            # [(contig, strand, start, end)] and the phase information into a list of phases
            # positions will thus have the same format for single-segment and multi-segment features
            # this format allows to have segments on different contigs and strands (e.g. trans-splicing)
        tuple_data = (db, feature, [(contig, strand, start, end)], score, [phase], attributes)
        return feature_id, tuple_data

    def reformat_feature_id(self, feature, attributes, database="refseq"):
        """ Returns the reformatted feature ID from a dict of attributes based on the type of feature,
        uniformizing database specific formatting for network construction """
        if feature not in {"region", "chromosome"}:
            if feature == "exon":
                feature_id = self.reformat_exon(attributes, database)
            else:
                feature_id = self.reformat_general_id(attributes, database)
        else:
            feature_id: str = self.get_contig_id(attributes, database)
        return feature_id

    @staticmethod
    def get_contig_id(attributes, database="refseq"):
        """ Returns the contig ID from a dict of attributes based on database specific formatting """
        index: int = 0
        if database == "ensembl":
            index = 1
        return attributes["ID"].split(":")[index]

    @staticmethod
    def reformat_exon(attributes, database="refseq"):
        """ Returns the exon ID from a dict of attributes based on database specific formatting """
        if database in {"refseq", "genbank"}:
            return attributes["ID"]
        elif database == "ensembl":
            parent = attributes["Parent"]
            suffix = attributes["rank"]
            return parent.replace("transcript:", "exon-") + "-" + suffix
        return None

    @staticmethod
    def reformat_general_id(attributes, database="refseq"):
        """ Returns the reformatted feature ID from a dict of attributes based on database specific formatting """
        if database == "ensembl":
            if "ID" in attributes:
                return reformat_ensembl_id(attributes["ID"])
        return attributes.get("ID", None)

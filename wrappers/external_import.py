from os.path import join
from utils.helpers.init_functions import initialize_dict
from utils.helpers.gff_functions import attributes_to_dict, format_gff_line
from utils.helpers.os_sys_functions import get_list_of_files


def import_interproscan_data(path_ipr, cds_storage, database):
    """
    Imports InterproScan data from a gff3 input file and returns a dict {cds: {interproscan_id: gff line}}.
    :param path_ipr: path to the input Interproscan GFF3 file
    :param cds_storage: CdsStorage object, CDS-specific storage object
    :param database: database name
    :return: cds_motifs_dict {cds ID: {domain ID: data tuple from interproscan GFF file}}
    """
    if path_ipr:
        print("\nImporting proteic features...")
        total = 0
        with open(path_ipr, "r") as f:
            for line in f:
                if line == "##FASTA\n":
                    break
                if not line.startswith("#"):
                    cds, tool, feature, pstart, pstop, score, strand, phase, attributes = line.rstrip().split("\t")
                    cds_id = "cds-" + cds
                    if database == "ensembl":
                        cds_id = "cds-" + cds.split(".")[0]  # removes the version suffix
                    cds_storage.cds_motifs_data = initialize_dict(cds_storage.cds_motifs_data, cds_id, {})
                    if feature == "protein_match":
                        pstart, pstop = int(pstart), int(pstop)
                        attributes = attributes_to_dict(attributes)
                        if "Dbxref" in attributes:
                            domain_id = attributes["Dbxref"].replace("InterPro:", "") + f"_{pstart}_{pstop}"
                            domain_id = f"dom-{cds}_" + domain_id.replace('"', '')
                        else:
                            domain_id = f"dom-{cds}_{attributes['Name']}_{pstart}_{pstop}"
                        if domain_id not in cds_storage.cds_motifs_data[cds_id]:
                            cds_storage.cds_motifs_data[cds_id][domain_id] = (
                                cds_id, tool, feature, pstart, pstop, score, strand, phase, attributes)
                            total += 1
        print(f"{total} proteic features imported from InterproScan.")
    return cds_storage


def import_repeatmodeler_data(path_rm):
    """
    Extracts interspersed repeats from the stockholm file produced by RepeatModeler; gives a unique ID to each repeat:
    id from repeatmodeler + contig + start + stop;  Updates gff_dict with interspersed repeat data from repeatModeler
    analysis
    :param path_rm: path to the input RepeatModeler Stockholm file
    :return: updated contig_storage, set of repeat ids
    """
    int_repeats = {}
    if path_rm:
        print("\nImporting interspersed repeats...")
        total = 0  # counter for features
        with open(path_rm, "r") as f:
            for line in f:
                if "# STOCKHOLM 1.0" in line:
                    feature_id = next(f)
                    feature_id = feature_id.rstrip().replace("#=GF ID    ", "")
                    next(f)
                    descr = next(f)
                    descr = descr.rstrip().replace("#=GF TP    ", "")
                    for i in range(4):  # skip 4 lines
                        next(f)
                elif "//" not in line:
                    if not line.startswith("#"):
                        coords = line.split("    ")[0]
                        if "//" not in coords:
                            contig = coords.split(":")[0]
                            int_repeats = initialize_dict(int_repeats, contig, {})
                            start, stop = coords.split(":")[1].split("-")
                            strand = "+"
                            if int(start) > int(stop):
                                strand = "-"
                                temp = start
                                start = stop
                                stop = temp
                            transp_id = f"intrep-{feature_id}_{contig}_{strand}_{start}_{stop}"
                            list_positions = [(contig, strand, int(start), int(stop))]
                            data = (
                                "RepeatModeler", "dispersed_repeat", list_positions,
                                ".", ["."], {"ID": transp_id, "contig": contig, "Note": descr})
                            int_repeats[contig][transp_id] = data
                            # contig_storage.gff.update(GffData((transp_id, data), is_tuple=True))
                            # contig_storage.positions.update(transp_id, list_positions)
                            total += 1
        print(f"{total} Interspersed Repeats imported from RepeatModeler.")
    return int_repeats


def import_tandemrepeatfinder_data(path_trf, contig_names):
    """
    Extracts tandem repeats (on the + strand) from the gff output files produced by TRF+TRAP,
     and compute their reverse-complementary (on the - strand)
    :param path_trf: path to the input TRAP gff folder
    :param contig_names: set of contig IDs
    :return: dict {contig ID: {repeat ID: repeat data}}
    """
    tandem_repeats = {}
    if path_trf:
        print("\nImporting tandem repeats...")
        total = 0  # counter for features
        list_files = [file for file in get_list_of_files(path_trf) if file.endswith(".gff")]
        for file in list_files:
            contig = get_contig_from_filename(file, contig_names)
            path = join(path_trf, file)
            with (open(path, "r") as f):
                for line in f:
                    if "\t" in line:
                        contig, db, feature, start, stop, score, strand, phase, attributes = parse_trftrap_gff_line(
                            line,
                            contig)
                        tandem_repeats = initialize_dict(tandem_repeats, contig, {})
                        repeat_id = attributes["ID"]
                        list_positions = [(contig, strand, start, stop)]
                        data = (db, feature, list_positions, score, [phase], attributes)
                        tandem_repeats[contig][repeat_id] = data
                        total += 1
        print(f"{total} Tandem Repeats imported from Tandem Repeat Finder.")
    return tandem_repeats


def parse_trftrap_gff_line(line, contig):
    """
    Converts each line from the trf/trap gff format file to gff3-like format:
    - start and stop will be converted to int.
    - attributes will be converted to dict
    - a unique ID will be generated for each repeat using: copy number x repeat unit size + contig + start + stop
    :param line: line read from a trf/trap output gff file
    :param contig: current contig/contig
    :return: tuple of all gff3 format data
    """
    line = line.split("\t")
    _, db, feature, start, stop, score, strand, phase, attributes = line
    attributes = [e.rstrip() for e in attributes.split('" "') if ("=" in e or ";" in e)]
    attr_dict = {}
    for e in attributes:
        if "=" in e:
            k, v = e.split(" = ")
            attr_dict[k.replace('"', '')] = v.replace('"', '')
        elif ";" in e:
            data = e.split(" ; ")
            source = data[0].replace('"', '')
            attr_dict["Note"] = source
            for d in data[1:]:
                k, v = d.split(" ")
                attr_dict[k.replace('"', '')] = v.replace('"', '')
    repeat_id = f"sat-{attr_dict['copy number']}x{attr_dict['repeat unit size']}_{contig}_{strand}_{start}_{stop}"
    attr_dict["ID"] = repeat_id
    return contig, db, feature, int(start), int(stop), score, strand, phase, attr_dict


def get_contig_from_filename(filename, contig_names):
    """
    Parses a filename for contig id and returns the contig id.
    :param filename: str
    :param contig_names: set {id_contig}
    :return: contig id or None if not found
    """
    for contig in contig_names:
        if contig in filename:
            file_data = set(filename.split("_"))
            contig_data = set(contig.split("_"))
            ref_len = len(contig_data)
            if len(contig_data & file_data) == ref_len:
                return contig
    return None


def import_custom_gff(path_custom):
    """
    Imports custom gff data as additional features.
    :param path_custom: path to the input file containing a list of paths to custom GFF3 annotation files.
    :return: list of dicts {contig: list(data tuple)}
    """
    list_data = []
    if path_custom:
        print("\nImporting custom gff3 files...")
        total = 0  # counter for features
        with open(path_custom, "r") as f:
            gff_paths = [path.rstrip() for path in f.readlines()]
            for path in gff_paths:
                data = {}
                with open(path, "r") as g:
                    lines = g.readlines()
                    for line in lines:
                        if not line.startswith("#"):
                            contig, db, feature, start, stop, score, strand, phase, attributes = format_gff_line(line)
                            if "ID" in attributes:
                                data = initialize_dict(data, contig, [])
                                data_tuple = (db, feature, [(contig, strand, start, stop)], score, [phase], attributes)
                                data[contig].append(data_tuple)
                                total += 1
                list_data.append(data)
        n = len(list_data)
        final = "file"
        if n > 1:
            final = "files"
        print(f"{total} Additional features imported from {n} GFF3 input {final}.")
    return list_data

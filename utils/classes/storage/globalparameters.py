from os.path import join, basename, exists
from shutil import rmtree
from utils.helpers.os_sys_functions import check_input_path, check_output_dir, get_list_of_dirs


class GlobalParameters:
    """
    Object storing and passing paths information and command line parameters to various functions:
    Attributes
    ----------
        simplify: int, either 0 (create the complete network), 1 (create a network without children features)
                    or 2 (create a gene-only backbone with additional features), 3 (create a gene-only network)
        only_contig: list[str], list of contig IDs to convert to networks

        keep_redundant: bool, flag to indicate whether to deduplicate (False, default) or keep redundant features (True)
        encode: bool, flag to indicate whether to encode feature IDs, edge types and edge classes as numbers
        generate_gff: bool, flag to indicate whether to generate GFF3 files for inferred features

        path_gff: str, path to the input genomic gff file
        path_output_dir: str, path to the output directory
        path_gfn_dir: str, path to the output network directory
        path_gff_dir: str, path to the output GFF3 files directory
        path_fna: str, path to the input genomic fasta file
        path_custom: str, path to the multiline list of paths to additional gff files
        path_rm: str, path to the input stockholm file from RepeatModeler
        path_trf: str, path to the input gff directory from TRF/TRAP
        path_ipr: str,  path to the input gff3 file from Interproscan
        path_ignored_features: str, path to a text file listing feature types to exclude from the network
        root: str, root file name
        database: str, either refseq, genbank or ensembl
    """
    def __init__(self):

        # network parameters
        self.simplify: int = 0  # default value, complete network
        self.only_contig: list[str] = []  # default value, contig IDs appended from command line parameter -smp
        # or from the complete contig IDs list (default)
        self.start: int = 0
        self.end: int = 0

        # boolean parameters
        self.keep_redundant: bool = False
        self.encode: bool = False
        self.generate_gff: bool = False
        self.serialize: bool = False
        self.has_domains: bool = False
        self.external: bool = False

        # boolean parameters for network simplification (default values correspond to a complete network)
        self.has_spatial_edges: bool = True
        self.has_ontogenetic_edges: bool = True
        self.complete_network: bool = True
        self.keep_children: bool = True
        self.genic_backbone: bool = False

        # path data
        self.path_gff: str = ""
        self.path_output_dir: str = ""
        self.path_gfn_dir: str = ""
        self.path_gff_dir: str = ""
        self.path_fna: str = ""
        self.path_custom: str = ""
        self.path_rm: str = ""
        self.path_trf: str = ""
        self.path_ipr: str = ""
        self.path_ignored_features: str = ""
        self.root: str = ""

        # database GFF3 format
        self.database: str = "refseq"  # default GFF3 format

    def import_parameters(self, args):
        """ Parses the command line arguments stored in the argparse object args """

        self.path_gff = args.path_gff
        self.path_output_dir = args.path_output_dir
        self.path_custom = args.path_custom_gffs
        self.path_fna = args.path_fna
        self.path_rm = args.path_rm  # optional
        self.path_trf = args.path_trf_trap  # optional
        self.path_ipr = args.path_ipr  # optional
        self.path_ignored_features = args.filter  # optional
        self.keep_redundant = args.keep_redundant
        self.encode = args.encode
        self.generate_gff = args.generate_gff
        self.serialize = True

        # define network export subdirectory
        self.path_gfn_dir = join(self.path_output_dir, "genomic_networks")

        # define optional GGF3 output files export subdirectory
        if self.generate_gff:
            self.path_gff_dir = join(self.path_output_dir, "generated_GFF3_annotation")

        # define GGF3 format of genome annotations
        if args.database in {"refseq", "genbank", "ensembl"}:
            self.database = args.database

        # define the type of network
        self.simplify = args.simplify if args.simplify in range(0, 6) else 0  # default value, complete network

        # define the list of contigs
        if args.only_contig:
            self.only_contig = args.only_contig.split(",")

        self.start = args.start if args.start else 0
        self.end = args.end if args.end else 0

        # compute root file name
        if args.root_name:
            self.root = args.root_name
        else:
            filename = basename(self.path_gff)
            self.root = ".".join(filename.split(".")[:-1])
        self._modify_root_name()

        # check and create the required directories
        self.check_paths()

        # test whether there are external imports
        self.has_domains = True if self.path_ipr else False
        self.external = any([self.path_custom, self.path_rm, self.path_trf])

        # adjust simplification boolean parameters

        if self.simplify == 5:
            self.has_spatial_edges = False
        if self.simplify in {1, 2, 3, 4}:
            self.has_ontogenetic_edges = False
        if self.simplify in {2, 3, 4, 5}:
            self.complete_network = False
        if self.simplify in {2, 3, 4}:
            self.keep_children = False
        if self.simplify in {4, 5}:
            self.external = False  # even if paths are provided, they will not be used for simplify 4 and 5
        if self.simplify in {3, 4}:
            self.genic_backbone = True

    def check_paths(self):
        """ Performs a check on the existence of files and directories for each path """

        # mandatory and optional input paths
        for path in [self.path_fna,
                     self.path_gff,
                     self.path_custom,
                     self.path_ipr,
                     self.path_rm,
                     self.path_trf,
                     self.path_ignored_features]:
            check_input_path(path)

        # output directories
        paths = [self.path_output_dir]
        if self.generate_gff:
            paths.append(self.path_gff_dir)
        for path in paths:
            check_output_dir(path)

    def get_path_network(self, contig):
        """Returns the path to save the output network file for the current contig"""
        return self._output_path(f"gfn__{self.root}_{contig}.csv", contig)

    def get_path_node_annotation(self, contig):
        """Returns the path to save the output node annotation file for the current contig"""
        return self._output_path(f"node_annotation__{self.root}_{contig}.csv", contig,
                                 subdirectory=f"annotation_{contig}")

    def get_path_edgetype_index(self, contig):
        """Returns the path to save the output node annotation file for the current contig"""
        return self._output_path(f"edgetype_index__{self.root}_{contig}.csv", contig,
                                 subdirectory=f"annotation_{contig}")

    def get_path_edgeclass_index(self, contig):
        """Returns the path to save the output node annotation file for the current contig"""
        return self._output_path(f"edgeclass_index__{self.root}_{contig}.csv", contig,
                                 subdirectory=f"annotation_{contig}")

    def _output_path(self, output_file, contig, subdirectory=None):
        """ Returns the path to the output directory for a given contig"""
        path_output_dir = self.path_output_dir
        path_output_subdirectory = join(path_output_dir, contig, subdirectory) if subdirectory \
            else join(path_output_dir, contig)
        check_output_dir(path_output_subdirectory)
        return join(path_output_subdirectory, output_file)

    def _modify_root_name(self):
        """ Builds the final root name depending on cutoff parameters"""
        root_extension = ""
        root_start = "start"
        root_end = "end"
        if self.start > 0:
            root_start = str(self.start)
        if self.end > 0:
            root_end = str(self.end)
        if root_start != "start" or root_end != "end":
            root_extension = f"_(from_{root_start}_to_{root_end})"
        if root_extension:
            self.root = self.root + root_extension

    def print_paths(self):
        """ Prints user information relative to paths """
        print(f"\nInput FASTA file: {self.path_fna}\n"
              f"\nInput GFF file: {self.path_gff}\n"
              f"\nOutput directory: {self.path_output_dir}\n"
              f"\nRoot name for output files: {self.root}\n")

    def remove_temp_directories(self) -> None:
        """ Cleanup method removing all temporary subdirectories """
        if exists(self.path_output_dir):
            sub_dirs = get_list_of_dirs(self.path_output_dir)
            temp_sub_dirs = [d for d in sub_dirs if d.startswith("temp_")]
            for temp_dir in temp_sub_dirs:
                rmtree(join(self.path_output_dir, temp_dir))

import argparse


class SmartFormatter(argparse.HelpFormatter):
    """ Allows new line character use in the help text"""
    def _fill_text(self, text, width, indent):
        return text.replace('\n', f'\n{indent}').strip()

def main_arg_parser():
    """.
    Converts command line arguments to parameters in an argparse object
    :return: object containing the parameters
    """
    parser = argparse.ArgumentParser(
        description="This script computes a genomic features network (GFN, multiDiGraph) and a node annotation file\n"
                    "for each contig in a genome, as well as optional gff3 annotation files (parameter: --generate_gff)\n"
                    "for inferred or integrated additional features.\n"
                    "It uses a genomic features file (gff3 format), and a genomic sequence file (fasta format) as its\n"
                    "only mandatory inputs. These should be formatted either in the Refseq, Genbank or Ensembl databases\n"
                    "style (the script uses Refseq style by default, use the --database parameter to change it).\n\n"
                    "The import of external annotation sources, such as output files from InterproScan, RepeatModeler\n"
                    "or Tandem repeat Finder analyses, is natively implemented. Custom annotation in the gff3 format\n"
                    "can also be imported.\n\n"
                    "Numeric encoding of feature IDs, edge types and classes can be used to decrease network size \n"
                    "on disk (parameter: --encode).\n\n"
                    "Filtering is implemented as well to: \n"
                    " (1) select specific contigs (parameter: --only_contig) \n"
                    " (2) specify a start and an end to only create network edges over a specific position interval\n"
                    "     (will apply to all selected contigs; parameters: --start, --end) \n"
                    " (3) adjust network simplification level (parameter: --simplify) \n"
                    " (4) adjust removal of redundant features (parameter: --keep_redundant) \n"
                    " (5) filter out specific types of features (parameter: --filter).",
                    formatter_class=SmartFormatter)
    # mandatory parameters
    parser.add_argument("-gff", "--path_gff", type=str, help="path to the input gff annotation file")
    parser.add_argument("-fna", "--path_fna", type=str, default="",
                        help="Path to the input genomic fasta file "
                             "(used to compute percent GC for each genomic feature)")
    # options
    parser.add_argument("-cus", "--path_custom_gffs", type=str, default="",
                        help="Path to an input file containing a list of paths to custom gff3 files, one path per line")
    parser.add_argument("-rm", "--path_rm", type=str, default="",
                        help="Optional: path to the input RepeatModeler "
                             "Stockholm file (interspersed repeats/transposon mapping)")
    parser.add_argument("-trf", "--path_trf_trap", type=str, default="",
                        help="Optional: path to the folder containing the input TRF/TRAP "
                             "gff files (tandem repeats mapping)")
    parser.add_argument("-ipr", "--path_ipr", type=str, default="",
                        help="Optional: path to the input InterproScan gff3 file")
    parser.add_argument("-o", "--path_output_dir", type=str,
                        help="path to the output directory (default: a new folder in the script directory"
                             ": ./genomic_features_network)",
                        default="./genomic_features_network")
    parser.add_argument("-n", "--root_name", type=str,
                        help="Optional: custom root name for output files "
                             "(default: root name of the provided gff file)", default="")
    parser.add_argument("-db", "--database", type=str,
                        help="Genomic GFF Database name "
                             "(default: refseq; supported: refseq, genbank, ensembl)", default="")
    parser.add_argument("-only", "--only_contig", type=str,
                        help="Optional: which contigs to use to build networks, use the same contig IDs as in "
                             "the Fasta genome sequence file, separated by a comma ',' (default: all contigs) ",
                        default="")
    parser.add_argument("-start", "--start", type=int, help="start position on the chromosome "
                                                            "(default: position 0); for multi-contig filtering,"
                                                            " the same range will be applied to all selected contigs",
                        default=0)
    parser.add_argument("-end", "--end", type=int, help="end position on the chromosome "
                                                        "(default: max position); for multi-contig filtering,"
                                                        " the same range will be applied to all selected contigs",
                        default=0)
    parser.add_argument("-filt", "--filter", type=str,
                        help="Optional: path to a text file listing feature types to exclude from the network."
                             " (These feature types cannot be filtered:"
                             " gene; rna; exon; intron; cds; sequence_feature) (default: No filtering)", default="")
    parser.add_argument("-smp", "--simplify", type=int,
                        help="Optional: outputs a simplified network, takes a value of: "
                             "0 (default): complete spatial and ontogenetic network;"
                             " 1: complete spatial network (no ontogenetic edges);"
                             " 2: spatial backbone without children features (with additional features);"
                             " 3: spatial gene-only backbone network (with additional features); "
                             " 4: spatial gene-only network (only made of genes and pseudogenes); "
                             " 5: ontogenetic network only (redundant features i.e. same type and same positions,"
                             " are kept)",
                        default=0, choices=list(range(0, 6)))
    parser.add_argument("-kr", "--keep_redundant", action='store_true',
                        help="Optional flag: if used, the script will not collapse feature nodes of the same type "
                             "with identical positions into a single node (default: deduplication of features).")
    parser.add_argument("-enc", "--encode", action='store_true',
                        help="Optional flag: if used, the script will encode node names, edge types and classes"
                             " as numbers to reduce file size"
                             " (recommended for large datasets, default: no encoding)")
    parser.add_argument("-gen", "--generate_gff", action='store_true',
                        help="Optional flag: if used, the script will generate additional GFF3 files"
                             " for inferred features (e.g. introns, intergenic regions, repeats..."
                             " default: no generation of GFF files)")
    return parser.parse_args()

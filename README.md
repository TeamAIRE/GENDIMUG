GenDiMug

A generic inclusive multigraph framework for genomic feature organization

GenDiMug is a Python-based framework that transforms standard Genomic annotations (GFF3) and sequence (Fasta) into complex Directed Multigraphs. Unlike linear genome browsers, GenDiMug models the genome as a network of ontogenetic (parent-child) and spatial (cis/trans) relationships, providing a multi-scale view of genomic architecture as well as a network based interface to combine with other -omics derived graphs.



Capabilities

•	Automated feature inference: Automatically derives features often missing from annotations, including intergenic regions, introns, and terminal regions (with specific handling for circular contigs/plasmids). Optional export of this inferences as GFF3 annotation files.
•	Support for fragmented features: Native support for features mapping to multiple intervals, such as CDS, trans-spliced mRNAs and multi-exon protein motifs.
•	Dual-relationship modeling:
•	Ontogenetic edges: Directed "Parent → Child" relationships between features (Gene -> (RNA) -> (Exon) -> CDS -> Protein Domain). Note that RNA/Exon informations are usually provided in eukaryotic genome annotations, but only occasionally in prokaryotic genome annotations.
•	Spatial edges: Connectivity based on genomic positions of features, including cis-(on the same strand) and trans- (trans-strand) sequential, overlap, nested, identical/reverse-complementary relationships.
•	High-performance overlap detection: Utilizes Interval Trees for efficient identification of overlapping and nested features across large-scale genomic datasets.

Multigraph Architecture

GenDiMug constructs a network composed of 12 node types and 9 edge types, allowing for highly granular filtering and topological analysis. The edge orientation is dictated by relative genomic positions, in ascending order for features annotated on the ‘+’ strand and in decreasing order on the ‘-‘ strand. The ontogenetic relationships are derived from the GFF3 ‘Parent’ feature attributes, and rewired to obtain the hierarchical order Gene -> (RNA) -> (Exon) -> CDS -> Protein Domain.

Node & Edge Taxonomy
Category	Types
Genomic Nodes	Gene, Pseudogene, RNA, Exon, Intron, CDS, Domain, Intergenic, Terminal, Tandem Repeat, Interspersed Repeat, Custom/other Features.
Cis-Spatial Edges	Sequential, Overlap, Nested, Identical.
Trans-Spatial Edges	Sequential, Overlap, Nested, Reverse-Complementary.
Hierarchical Edges	Ontogenetic (Parent -> Child).

GenDiMug workflow (complete network)

1.	Feature collection & inference
Collects ORF/Gene positions to infer intergenic regions and terminal ends. For circular contigs, terminal nodes are unified to bridge the graph ends. Parent -> child relationships are retrieved from feature attributes.
2.	Backbone construction
Creates a linear spatial backbone per strand, linking ORFs (genes and pseudogenes) and intergenic regions via cis-sequential/overlap edges. Additional nested ORFs are then layered onto this backbone.
3.	Intra-ORF edges
Exons and Introns are linked via cis-sequential edges. CDS fragments are connected to their corresponding coding exons. Optional additional elements are integrated (protein motifs, DNA repeats…).
4.	Cis/Trans-relationship mapping
Computes spatial relationships between features on the same and opposite strands (e.g., trans-overlap, reverse-complementary) to complete the multigraph.
5.	Optional deduplication of redundant features
Features of the same type with the same genomic positions are represented by a single node. Nested DNA repeats are removed.


Quick Start

Requirements

•	Python > 3.10
Additional python packages:
•	Intervaltree (https://github.com/chaimleib/intervaltree,   https://pypi.org/project/intervaltree/)
•	BioPython (https://biopython.org/,    https://pypi.org/project/biopython/)
•	Networkx (https://networkx.org/,    https://pypi.org/project/networkx/)


Basic Usage

	Python command line:
	python3 gendimug.py 

	Main parameters
	-gff [PATH TO GFF3 ANNOTATION FILE] 
	-fna [PATH TO FASTA GENOME SEQUENCE FILE] 

	Optional parameters:

	Optional features input
	-cus , --path_custom_gffs [PATH TO INPUT CUSTOM GFF3 FORMAT FEATURE FILE]
	-rm, --path_rm [PATH TO INPUT STOCKHOLM FORMAT REPEATMODELER RESULT FILE] 
	-trf, --path_trf_trap [PATH TO INPUT TANDEM REPEAT FINDER/TRAP GFF RESULT DIRECTORY]
	-ipr, --path_ipr [PATH TO INPUT INTERPROSCAN GFF3 RESULT FILE]
	
	Optional output parameters
	-o, --path_output_dir [PATH TO OUTPUT DIRECTORY]
	-n, --root_name [CUSTOM ROOT NAME FOR OUTPUT FILES]
	
	Optional modifiers
	-db, --database [GENOME DATABASE: refseq/genbank/ensembl]
	-smp, --simplify [0: 5]
	
	Optional flags
	-kr, --keep_redundant : when used, the script does not deduplicate features
	-enc, --encode:        when used, feature IDs, edge types and classes are numerically encoded
	-gen, --generate_gff: when used, generates GFF3 annotation files for inferred features
	
	Optional filters
	-only, --only_contig [CONTIG IDS LIST, comma-separated] : Build network only for the indicated contigs
	-start, --start [CUTOFF START POSITION] : Cutoff start position for all selected contigs
	-end, --end [CUTOFF END POSITION] : Cutoff end position for all selected contigs
	-filt, --filter [PATH TO A TEXT FILE LISTING FEATURE TYPES TO EXCLUDE FROM THE NETWORK]


Output files

GenDiMug outputs a network edge file (prefixed by ‘gfn__’ as genome features network) in csv format for each contig in a separate directory named by contig ID, as in the genome annotation file. The node annotation file (and edgetype/edgeclass numeric indices files, if using the ‘—encode’ parameter) are saved in the subdirectory ‘annotations’ in each contig specific directory. A summary file (number of nodes and edges by contig, total number of nodes and edges) is saved in the general output directory.

Data Attributes
Each element in the generated graph carries rich metadata:
•	Edges: Include physical distance (bp), overlap length, and specific edge classes (e.g., exon_cis-sequential_intron).
•	Nodes: Include original GFF3 data, computed local GC content, and a list of genomic position intervals (handling fragmentation).


Simplified networks

5 options are available using the –simplify parameter :

o	(0) Complete spatial and ontogenetic genome features network
o	(1) Complete spatial genome features network
o	(2) Spatial genome features network without children features (includes additional features)
o	(3) Spatial genome features network without children features, with a gene-only backbone (includes additional features)
o	(4) Gene-only spatial network without children features (no additional features)
o	(5) Ontogenetic network only (redundant features are kept)

Curation data for protein motifs/domains

When importing protein motifs or domains from InterproScan, the script attempts to map domains to their encoding genomic positions. Discrepancies between protein and CDS sequences can lead to failure of mapping some domains, which will be captured in a supplementary ‘unmapped domains’ output directory.

Additionally, the script checks that mapped domain genomic positions add up to a length that is a multiple of 3, and captures any domain that fails that check in the ‘domains_of_length_not_multiple_of_3’ supplementary output directory, with exception data( e.g. “ribosomal slippage”, “low quality sequence”…), if available in the CDS attributes in the genomic annotation file.

Network visualization

The csv edge files produced by GenDiMug can be directly imported into graph visualization softwares. The maximum network size that can be managed by the software depends on available memory; however it is always possible to generate a smaller network for visualization purposes using the ‘--start’ and ‘--end’ parameters. 
Here is the procedure to visualize genome features networks in Cytoscape (https://cytoscape.org) :

Network file import (no encoding)

file : gfn__[...].csv
File > Import > Network from File

Node annotation file import

file : node_annotation__[...].csv
File > Import > Table from File
Import Data as: Node Table Columns
Key Column for Network: 'node'

With numerical encoding activated (‘--encode’ parameter)

Numerically encoded network file import

file : gfn__[...].csv
File > Import > Network from File

Numerically encoded node annotation file import

file : node_annotation__[...].csv
File > Import > Table from File
Import Data as: Node Table Columns
Key Column for Network: 'node'

Numerically encoded edgetype annotation file import

file : edgetype_index__[...].csv
File > Import> Table from File
Import Data as: Edge Table Columns
Key Column for Network: 'et_index'

Numerically encoded edgeclass annotation file import

file : edgeclass_index__[...].csv
File > Import> Table from File
Import Data as: Edge Table Columns
Key Column for Network: 'ec_index'

Use the node and edge style options to customize the rendering of nodes and edges.

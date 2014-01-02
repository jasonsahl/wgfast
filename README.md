####Header####
The whole genome focused array SNP typing (WG-FAST) pipeline
written by: Jason Sahl
Email: jasonsahl at gmail dot com

####overview####
The goal of WG-FAST is to phylogenetically genotype an unknown
sample in the context of a well studied pathogen.  This sample
can either be a metagenomics dataset, a metatranscriptomics dataset,
or a single isolate sequencing dataset

####Limitations####
Currently, WG-FAST can only accept paired-end sequence data, with
a heavy emphasis on the Illumina platform.

####External dependencies####
WG-FAST relies on:

1.  RaxML - tested version is 7.7.8 (PTHREADS), must be in your path as 'raxmlHPC-PTHREADS'
citation: 'Stamatakis A. RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics. 2006;22(21):2688-90'
2.  Samtools - tested version is 0.1.19, must be in your path as 'samtools'
citation: 'Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, Genome Project Data Processing S. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-9'
3.  BWA - tested version is 0.7.5a, must be in your path as 'bwa'
citation: 'Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXivorg. 2013(arXiv:1303.3997 [q-bio.GN])'
4.  GATK - tested version is 2.7.2, included with WG-FAST, but also requires Python 1.7+
citation: 'McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research. 2010;20(9):1297-303'
5.  Picard tools - tested version is 1.79, included with WG-FAST
6.  DendroPy - tested version is 3.12.0, must be in PYTHONPATH
7.  BioPython - must be in PYTHONPATH

####Required input####
1.  Directory of sequence reads.  The reads must be named according to Illumina Hi-Seq or Mi-Seq conventions.
The reads must be in the Illumina 1.9+ FastQ format.  If you have old Illumina FastQ encodings, they must be
converted before running WG-FAST.
2.  SNP matrix.  The easiest way to generate this is by using Sniper (link to REPO).
3.  Phylogeny.  A script is included with WG-FAST that can generate an appropriate phylogeny
from the SNP matrix.
4.  Reference genome in FASTA format.  This should be the same FASTA that was used to generate your Sniper
SNP matrix.

####Installation####
1.  To install, enter the wgfast directory and type the command below.  If you don't have suo
privileges, install to base directory using --prefix:
sudo python setup.py install
2.  Open the script (wgfast.py) with a text editor and change the path to your WG-FAST installation directory.
For example:
WGFAST_PATH="/Users/jsahl/wgfast"
3.  To verify your installation, enter the wgfast directory and type the command below.  If everything
is working correctly, all tests should pass:
python tests/test_all_functions.py

####Test data####
1.  A sample SNP matrix is included with WG-FAST:
test_data/nasp_sample_matrix.tsv
2.  A sample phylogeny that is compatible with WG-FAST:
test_data/nasp_raxml.tree
3.  To reduce the size of the repository, raw sequence reads are not included.  Raw reads
from E. coli, downloadable from Genbank, can be used to test the method

####Output printed to screen####
1.  number of SNPs in genome X : these are actual polymorphic positions
2.  number of discarded SNPs in genome X : these are positions that were thrown away
because they didn't pass the coverage or proportion filters (see below)
3.  number of callable positions in genome X : these are all callable positions in your
query genome, both mono and polymorphic
4.  Insertion likelihood values: values from RAxML.  The last column is the number of possible
insertion nodes.  The lower the number, the more likely the placement
5.  maximum subsample distance between Reference and X: Only seen with re-sampling routine.
This shows the maximum distance observed


####Output files####
1.  nasp_matrix_with_unknowns.txt - The original matrix with unknowns pasted on.  However,
non-SNP information has been stripped from this matrix
2.  coverage_out.txt - Average depth of coverage of raw reads against the reference
3.  transformed.tree - Tree compatible with FigTree.  Your unknown samples will be shown
in red
4.  all_patristic_distances.txt - Pairwise patristic distances between all samples in your
analysis
5.  *.subsample.distances.txt - Sub-sampled distances when the re-sampling sub-routine is
invoked
6.  classification_results.txt - RAxML insertion likelihoods for each unknown into the tree


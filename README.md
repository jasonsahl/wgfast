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

###External dependencies###
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

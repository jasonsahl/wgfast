## The whole genome focused array SNP typing (WG-FAST) pipeline  

### Updated 7/16/2019  

#### Citation:  
Jason Sahl, James Schupp, David Rasko, Rebecca Colman, Jeffrey Foster, Paul Keim (2015).
Phylogenetically typing bacterial strains from partial SNP genotypes observed from direct
sequencing of clinical specimen metagenomic data  

#### Contact:  
Please address queries, conserns, improvements to jasonsahl at gmail dot com  

#### What does *WG-FAST* do?  

WG-FAST was designed as a tool to phylogenetically genotype unknown samples, even those
with extremely low read coverage, in the context of a well-studied dataset.

#### What does WG-FAST not do?  

WG-FAST is not intended to identify new SNPs in a dataset (i.e replace a de novo analysis).
If too many samples are processed with WG-FAST, a phylogenetic discovery bias can exist.

#### Installation  
See Readme

#### Unit tests  
-To test that all functions are working correctly in WG-FAST, type:  
```python tests/test_all_functions.py```

#### Required input files  
1. **Directory of sequence reads ("-d")**. The reads must be named according to Illumina HiSeq
or MiSeq conventions. Reads must be in the Illumina 1.9+ FastQ format. If you have old
Illumina FASTQ encodings, they must be converted before running WG-FAST. **Important: names
must not have periods ".", brackets "[]", or other weird characters "=:" in the header.**  
2. **Directory of reference files ("-r")**. This directory should only contain the following
files:  
a. **SNP matrix (must end in ".tsv")**. The easiest way to generate this is by using NASP
(https://github.com/TGenNorth/NASP). If other SNP matrix formats are used, they must conform
to having hte first column including (contig::coordinate) and the column following the SNP calls
must be (#SNPcall). For the sub-sampling routine to complete, a genome must be present in your matrix
that is called 'Reference'. **Important: sample names must not have periods in the header.**  
b. **Phylogeny (must end in .tree)**. A script is included with WG-FAST that can generate an appropriate
phylogeny for a NASP matrix (see below). This script also generates a 'Parameters' file, which can be used
with WG-FAST and cuts down on the computational time required for each subsequent run. The best way to guarantee downstream
compatibility is to generate the phylogeny with RAxML, which is used by the wgfast_prep.py script described below.  
c. **Reference genome in FASTA format (must end in .fasta)**. This should be the same FASTA that was used to call SNPs with NASP.  
d. **RAxML paramters file (must end in .PARAMS) (optional)**. If provided, RAxML will run faster because the likelihood has already
been calculated. Make sure that you use the conda version of RAxML for full compatibility.  

#### Complete list of arguments:  

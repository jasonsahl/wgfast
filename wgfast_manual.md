## The whole genome focused array SNP typing (*WG-FAST*) pipeline  

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

#### What does *WG-FAST* not do?  

*WG-FAST* is not intended to identify new SNPs in a dataset (i.e replace a de novo analysis).
If too many samples are processed with WG-FAST, a phylogenetic discovery bias can exist.

#### Installation  
See Readme

#### Unit tests  
-To test that all functions are working correctly in *WG-FAST*, type:  
```python tests/test_all_functions.py```

#### Required input files  
1. **Directory of sequence reads ("-d")**. The reads must be named according to Illumina HiSeq
or MiSeq conventions. Reads must be in the Illumina 1.9+ FastQ format. If you have old
Illumina FASTQ encodings, they must be converted before running *WG-FAST*. **Important: names
must not have periods ".", brackets "[]", or other weird characters "=:" in the header.**  
2. **Directory of reference files ("-r")**. This directory should only contain the following
files:  
a. **SNP matrix (must end in ".tsv")**. The easiest way to generate this is by using NASP
(https://github.com/TGenNorth/NASP). If other SNP matrix formats are used, they must conform
to having hte first column including (contig::coordinate) and the column following the SNP calls
must be (#SNPcall). For the sub-sampling routine to complete, a genome must be present in your matrix
that is called 'Reference'. **Important: sample names must not have periods in the header.**  
b. **Phylogeny (must end in .tree)**. A script is included with *WG-FAST* that can generate an appropriate
phylogeny for a NASP matrix (see below). This script also generates a 'Parameters' file, which can be used
with *WG-FAST* and cuts down on the computational time required for each subsequent run. The best way to guarantee downstream
compatibility is to generate the phylogeny with RAxML, which is used by the wgfast_prep.py script described below.  
c. **Reference genome in FASTA format (must end in .fasta)**. This should be the same FASTA that was used to call SNPs with NASP.  
d. **RAxML paramters file (must end in .PARAMS) (optional)**. If provided, RAxML will run faster because the likelihood has already
been calculated. Make sure that you use the conda version of RAxML for full compatibility.  

#### Complete list of arguments:  

-h, --help  show this help message and exist  
-r, REFERENCE_DIR, path to reference directory [REQUIRED]  
-d, READ_DIRECTORY, path to directory of fastq.gz files [REQUIRED]  
-p, PROCESSORS, # of processors to use - defaults to 2  
-c, COVERAGE, minimum SNP coverage required to be called, defaults to 3x  
-o, PROPORTION, proportion of alleles to be corrected, defaults to 0.9  
-k, KEEP, keep temporary files? Defaults to F (T or F)  
-s, SUBSAMPLE, run subsample routine? Defaults to T (T or F)  
-n, SUBNUMS, number of subsamples to process, defaults to 100  
-g, DOC, run depth of coverage on all files? Defaults to T (T or F)  
-e, TMP_DIR, temporary directory path for GATK, default = '/tmp'  
-z, insertion method (MP or ML), defaults to ML, MP not well tested  
-f, How close does a subsample have to be from true placement, Defaults to 0.1 (Float)  
-y, Only run sub-sample routine and exit? Defaults to F (T or F)  
-j, MODEL, which model to run with raxml (GTRGAMMA, ASC_GTRGAMMA)  

#### Test Data:  

1. To give *WG-FAST* a try to make sure everything is installed correctly, check out the test_data
directory.  
2. The following command was run from a “run” directory created within the WG-FAST
installation:  
```python ../wgfast.py -r ../test_data -d ../test_data/reads/ -p 4 -c 1 -s F```  

#### What is *WG-FAST* doing?  

*WG-FAST* is a pipeline that wraps other tools together to place the samples into a phylogenetic
context, based on a well-characterized dataset. The order of functions conducted by *WG-FAST*
include:  
1. Adapters are trimmed using BBduk. A list of adapters is included in the “bin” directory
(“illumina_adapters_all.txt”) within the *WG-FAST* distribution. This file includes all of the
standard Illumina adapters that I could find, plus additional ones that we use in our
laboratory. If you have different adapters than these, add the sequences into this file when
you run *WG-FAST*. The minimum length of sequences to keep is hard coded as 50, but I
could change this value to reflect different lengths if needed.  
2. A dictionary is created from the Reference fasta with Picard Tools and the reference fasta is
indexed with SAMtools.  
3. Reads are paired into a single sample based on their names, although single end reads are
also supported. *WG-FAST* assumes that names be something like: "S1_R1_001.fastq.gz” or
“R1_001.fastq.gz”. If the name pairings aren’t recognized, they will be run as single ended.  
4. Reads are mapped to the reference with minimap2, using default settings. SNPs are then
called with the Haplotype caller method in GATK v4.  
5. If selected, coverage across the reference is calculated with the “DepthOfCoverage” method
in SAMtools.  
6. SNPs are compared against the SNPs in your initial SNP matrix. If the SNP is missing from
your VCF file, a “-“ is inserted. Furthermore, if a position fails either the user-defined
thresholds for depth or proportion, the position is replaced with a “-“. The number of
discarded SNPs, or those that didn’t pass the filters, is reported.  
7. The new matrix is converted into a multi-fasta and processed with the evolutionary placement
algorithm in RAxML v8, placing the unknown into a user-provided phylogeny. The resulting
tree is converted into Nexus, such that the unknowns can be easily visualized (Red) in
FigTree (http://tree.bio.ed.ac.uk/software/figtree/).  
8. For the optional sub-sampling routine, the following functions are performed:  
**a.** For each query, the two closest genomes are identified, based on lowest pairwise
patristic distances calculated by DendroPy  
**b.** Two matrices are created for each query, one for each of the nearest neighbors. Each
of these contains the name of the query and the name of the neighbor, ending in
“tmp.matrix”. If you generate these files (and keep them), then run *WG-FAST* in “-y T”
mode, these files will be skipped and will not need to be created again.  
**c.** Each matrix is then converted into a FASTA file and a parameters file is created with
RAxML. The thought process here is that once the parameters files are created, they
can be used for future comparisons. If you run *WG-FAST* with “-y T”, you can re-use
these parameters files and won’t need to generate again. This method uses the
PTHREADS method to take advantage of multiple processors. Currently, only four
separate parameters files can be generated concurrently; this value is hardcoded into
the script. This is due to the large computational demand required by this method, but
could be made into a tunable parameter, if needed.  
**d.** The initial tree is pruned of the genome to be re-inserted. Each sub-sampled neighbor
is then inserted back into the pruned phylogeny, using the GTRGAMMA method, to
stay consistent with the parameters file. The patristic distance is then calculated
between the neighbor and the “Reference”. This is done because if a NASP-formatted
matrix is generated, the reference FASTA file will always be called “Reference”. The
“true” distance is also identified based on the original SNP matrix using all of the
original source data. The sub-sampled distance is then divided by the true distance in
order to get a ratio of placement replication. If this value falls within a user-defined
threshold, then the placement is considered to be correct. For each unknown, the
number of incorrect placements is divided by the total number of placements in order
to get an idea of how stable the placement is in that region of the tree at the level of
SNP density provided by the unknown. The resulting p-value can then be used to
assess the robustness of a placement, based on significance values chosen by the
user.  

Try running the script with the subsample routine with;  
```python ../wgfast.py -r ../test_data -d ../test_data/reads/ -p 4 –c 1```

#### Output printed to screen:  

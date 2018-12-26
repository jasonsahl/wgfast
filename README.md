# WG-FAST
The whole genome focused array SNP typing (WG-FAST) pipeline
written by: Jason Sahl
Email: jasonsahl at gmail dot com

for a more comprehensive overview, look at the Manual

## Overview
The goal of WG-FAST is to phylogenetically genotype an unknown
sample in the context of a well studied pathogen.  This sample
can either be a metagenomics dataset, a metatranscriptomics dataset,
or a single isolate sequencing dataset

## Installation
Tested with a Python 3.6 environment.  
To create environment
```
conda create -n [NAME] python=3.6
```
Activate environment
```
source activate [NAME]
```
Install WGFAST
```
conda install wgfast -c bioconda -c Fofanov -c conda-forge
```

Run test data from github repository
`wgfast -r path/to/test_data -d path/to/test_data/reads/ -p 4 -c 1 -s F`

## WGFAST setup
Running WGFAST requires a path to a directory of reference files `-r` that **must only contain 4 files with the file extensions in each example.**
1. **SNP Matrix** (example: bestsnps.tsv) generated using [NASP](https://github.com/TGenNorth/NASP). 
2. **Phylogeny** (example: nasp_raxml.tree).  A RAxML maximum likelihood tree can be generated from a NASP matrix using `wgfast_prep`. 
3. **Reference genome** in Fasta format (reference.fasta). This should be the same Fasta that was used to call SNPs with NASP.
4. **RAxML parameters file** (nasp.PARAMS). This is also provided when running `wgfast_prep`. If provided, RAxML will run faster because the likelihood has already been calculated. Make sure that the version of RAxML that was used to make this file is the same version that is used in current run or you will receive an error. (this should not be a problem if using the same install of WGFAST for both setup and running.

#### Usage `wgfast_prep`
```
$ wgfast_prep

Must provide matrix.

Usage: wgfast_prep [options]

Options:
  -h, --help            show this help message and exit
  -m MATRIX, --snp_matrix=MATRIX
                        path to NASP snp_matrix [REQUIRED]
  -o MODEL, --model=MODEL
                        model for RAxML, defaults to ASC_GTRGAMMA, can also be
                        GTRGAMMA
  -p PROCESSORS, --processors=PROCESSORS
                        number of processors to use with GTRGAMMA, defaults to
                        4
```

## WGFAST automatic setup
`wgfastdb` will automate the download and curation of bacterial genomes using [NCBITK](https://github.com/andrewsanchez/NCBITK) and [GenBankQC](https://github.com/andrewsanchez/GenBankQC), respectively. It will then generate all of the WGFAST required files using the curated genomes.  

`wgfastdb` can be run using a species or list of species from the command line or by passing a config file (Note: the species name must match exactly a species directory at ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/). The config file should contain comma separated values and must contain each species name and the ID of the reference genome that should be used. Other parameters can be customized for each individual species or the column can be left out if all default values should be used. If passing values at the command line, each of the parameters must be lists that are the same length as the list of species passed. If a single value is passed for a parameter, it will be applied to all species. 

Example Config:
`wgfastdb.cfg`
```
species,reference,unknowns,contigs,distance,assembly_size
Clostridium_botulinum,GCA_000017045,200,3.0,3.0,3.0
Rickettsia_prowazekii,GCA_001602215,200,3.0,3.0,3.0
Francisella_tularensis,GCA_000833335,200,3.0,3.0,3.0
Burkholderia_pseudomallei,GCA_000011545,200,3.0,3.0,3.0
Mycobacterium_tuberculosis,GCA_000016145,200,3.0,3.0,3.0
Burkholderia_mallei,GCA_000011705,200,3.0,3.0,3.0
Chlamydia_psittaci,GCA_000298435,200,3.0,3.0,3.0
Escherichia_coli,GCA_000005845,200,3.0,3.0,3.0
Coxiella_burnetii,GCA_000007765,200,3.0,3.0,3.0
Salmonella_enterica,GCA_001558355,200,3.0,3.0,3.0
Campylobacter_jejuni,GCA_000148705,200,3.0,3.0,3.0
Yersinia_enterocolitica,GCA_000834195,200,3.0,3.0,3.0
Listeria_monocytogenes,GCA_002101275,200,3.0,3.0,3.0
Yersinia_pestis,GCA_000009065,200,3.0,3.0,3.0
Bacillus_anthracis,GCA_000008445,200,3.0,3.0,3.0
```

#### Step 1: Download sequences
If running `wgfastdb` for the first time in a new directory (`PATH`), the sequences for the species listed will be downloaded, formated, and unzipped into `PATH/genomes/SPECIES_NAME/*.fasta`. If you are updating an existing database, make sure to call `wgfastdb` from the root `PATH`, not the `genomes` directory. By default, if the collection already exists, `wgfastdb` will check for updates and make any changes; to skip this step and use the collection as is, pass `--no_update` and `--no_assembly_update`. 

#### Step 2: Sequence curation
Next a quality control step is performed on the collection of genomes. This step has 4 parameters that control which genomes pass the filters: `--unknowns`, `--contigs`, `--assembly_size`, and `--distance`. See usage below or [GenBankQC](https://github.com/andrewsanchez/GenBankQC) for more details.  The fasta files that passed the filters are located at `PATH/genomes/SPECIES_NAME/qc/200-3.0-3.0-3.0/passed/*.fasta`, where `200-3.0-3.0-3.0` are the parameter values set for `--unknowns`, `--contigs`, `--assembly_size`, and `--distance`. 

#### Step 3: Matrix and Tree Generation
Using only the sequences in `PATH/genomes/SPECIES_NAME/qc/200-3.0-3.0-3.0/passed/*.fasta`, a NASP matrix is created using the genome that matches the `--reference` id (`'GCA_XXXXXXXXX*.fasta`) as the reference. This reference must exist in `PATH/genomes/SPECIES_NAME/*.fasta`. When complete, the matrix will be used to generate the RAxML .tree and .PARAMS file. All three files will be placed at `PATH/wgfast/SPECIES_NAME/` using the names listed above, with all intermediate files removed. 

#### Running
When `wgfastdb` is called, the download and curation processes are called as subprocesses of the main process, and the `threads` parameter will be passed and used by each serial curate (step 2) subprocess call. The processes in step 3 are run as a snakemake file that will spawn new processes for each step. Snakemake parameters can be passed at the command line and a snakemake config file can define resources to be used by each rule of the snake file. The `--jobs` or `--cores` snakemake parameter is important to tell snakemake the maximum number of cores/jobs to use.  

Each step is broken up into a number of different Snakemake rules. Each rule can be submitted as a job on a cluster using commands for the appropriate scheduler. This is accomplished by creating a cluster configuration json file that outlines the resources that should be applied to each rule. There can also be a default setting which will be the fall-back position for rules that are not provided a specific configuration. Any of the variables defined in the json can be accessed when defining the cluster script. The values in the cluster configuration file can be accessed using the cluster keyword when passed to the `--cluster` argument.

An example using slurm:
`srun --time 8:00:00 --mem 20000 wgfastdb ./ --config wgfastdb.cfg --threads 8 --cluster-config cluster.cfg --cluster "sbatch --mem {cluster.mem} --time {cluster.time} --job-name (cluster.job-name} --output ./Logs/{cluster.job-name}_%A.log --cpus {cluster.cpus}" --jobs 24` 

For each job, cluster.mem, cluster.time, and cluster.job-name values will be pulled from the cluster.cfg file (see example below) under the name of the rule that is being executed. Certain jobs are allocated a maximum number of threads but it will not use more than specificed by `--jobs/--cores`. If not specified, Snakemake assumes that only one core is available and will therefore scale back the number of threads to one.

In this case, requesting 24 jobs means that snakemake will run up to 24 jobs at a time for rules that have a max number of one thread. For rules like nasp_matrix that can run a max of 8 threads, only 3 jobs can be run at the same time.

The `srun --time 8:00:00 --mem 2000` at the beginning of the line defines the resources for running the main `wgfastdb` process which includes the download and curate steps.

`cluster.cfg`
```
{
    "__default__" :
    {
        "job-name": "wgfastdb",
        "account" : "username",
        "time" : "20:00",
        "mem" : 2000
    },
    "remove_duplicates" :
    {
        "job-name": ""remove_duplicates",
        "time" : "10:00",
    },
    "frankenfastas" :
    {
        "job-name": "frankenfastas",
        "time": "5:00",
    },
    "nasp_matrix" :
    {
        "job-name": "nasp_matrix",
        "time" : "1:00:00",
        "mem" : 60000
    },
    "validate_reference" :
    {
        "job-name": "validate_reference",
    },
    "matrix_to_fasta" :
    {
        "job-name" : "matrix_to_fasta",
        "time" : "30:00",
    },
    "raxml_tree" :
    {
        "job-name": "raxml_tree",
        "time": "00:30:00",
    },
    "raxml_params" :
    {
        "job-name": "raxml_params",
        "time" : "1:00:00",
        "mem" : 30000
    }
}
```





#### Usage `wgfastdb`

```
$ wgfastdb
usage: wgfastdb [-h] [--no_update] [--no_assembly_update]
                (--species SPECIES [SPECIES ...] | --config CONFIG)
                [--unknowns UNKNOWNS [UNKNOWNS ...]]
                [--contigs CONTIGS [CONTIGS ...]]
                [--assembly_size ASSEMBLY_SIZE [ASSEMBLY_SIZE ...]]
                [--distance DISTANCE [DISTANCE ...]]
                [--reference REFERENCE [REFERENCE ...]] [--log LOG]
                [--threads THREADS]
                PATH

Database setup for WGFAST

optional arguments:
  -h, --help            show this help message and exit
  --log LOG, -l LOG     Set log file path (default: ./wgfastdb.log)
  --threads THREADS, -t THREADS
                        Number of worker threads to spawn for curation.
                        (default: 4)

Download:
  Download and update sequences from NCBI

  PATH                  Path to existing fasta database files or location
                        where the files should be downloaded.
  --no_update           Do not sync your collection with the latest assembly
                        versions (default: False)
  --no_assembly_update  Do not download the latest assembly summary and
                        taxonomy dump and use your local copies. (default:
                        False)
  --species SPECIES [SPECIES ...], -s SPECIES [SPECIES ...]
                        List of species to build database. The species name
                        must match exactly a species directory at
                        ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/
                        (default: None)
  --config CONFIG, -cfg CONFIG
                        Path to config table (.csv). Table MUST include a
                        "species" column and may include a column "reference"
                        for reference genomes. Curation parameters can also be
                        set in the table using column headers: "unknowns",
                        "contigs", "assembly_size", and "distance". If a
                        parameter cell is left blank it will be replaced with
                        the default value or a value passed from the command
                        line (default: None)

Curate:
  Curate genomes

  --unknowns UNKNOWNS [UNKNOWNS ...], -n UNKNOWNS [UNKNOWNS ...]
                        Maximum number of unknown bases (not A, T, C, G) for
                        curation. If more than one value is passed, the list
                        must be the same length as the number of species.
                        Otherwise the same value is applied to all species.
                        (default: [200])
  --contigs CONTIGS [CONTIGS ...], -c CONTIGS [CONTIGS ...]
                        Acceptable deviations from median number of contigs
                        for curation. If more than one value is passed, the
                        list must be the same length as the number of species.
                        Otherwise the same value is applied to all species
                        (default: [3.0])
  --assembly_size ASSEMBLY_SIZE [ASSEMBLY_SIZE ...], -a ASSEMBLY_SIZE [ASSEMBLY_SIZE ...]
                        Acceptable devations from median assembly size for
                        curation. If more than one value is passed, the list
                        must be the same length as the number of species.
                        Otherwise the same value is applied to all species
                        (default: [3.0])
  --distance DISTANCE [DISTANCE ...], -d DISTANCE [DISTANCE ...]
                        Acceptable deviations from median MASH distances for
                        curation. If more than one value is passed, the list
                        must be the same length as the number of species.
                        Otherwise the same value is applied to all species
                        (default: [3.0])

Tree:
  SNP matrix and tree generation

  --reference REFERENCE [REFERENCE ...], -r REFERENCE [REFERENCE ...]
                        Define which genome to use as reference by providing
                        accession number (GCA_XXXXXXXXX). This list should be
                        the same length as the number of species. This is
                        REQUIRED from the command line or in the config file.
                        (default: None)
```

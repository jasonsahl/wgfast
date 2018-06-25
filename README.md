#### WG-FAST ####
The whole genome focused array SNP typing (WG-FAST) pipeline
written by: Jason Sahl
Email: jasonsahl at gmail dot com

for a more comprehensive overview, look at the Manual

#### overview ####
The goal of WG-FAST is to phylogenetically genotype an unknown
sample in the context of a well studied pathogen.  This sample
can either be a metagenomics dataset, a metatranscriptomics dataset,
or a single isolate sequencing dataset

#### Installation ####
Tested with a Python 3.6 environment
`conda install wgfast -c bioconda -c tara_furstenau`
Run test data from github repository
`wgfast -r path/to/test_data -d path/to/test_data/reads/ -p 4 -c 1 -s F`


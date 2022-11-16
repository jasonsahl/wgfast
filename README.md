#### *WG-FAST*
The whole genome focused array SNP typing (*WG-FAST*) pipeline
written by: Jason Sahl
Email: jasonsahl at gmail dot com

for a more comprehensive overview, look at the Manual

#### overview
The goal of *WG-FAST* is to phylogenetically genotype an unknown
sample in the context of a well studied pathogen.  This sample
can either be a metagenomics dataset, a metatranscriptomics dataset,
or a single isolate sequencing dataset

#### Installation
1.  The easiest way to install is through conda:  

```conda create -n wgfast python=3.6```  
```source activate wgfast```  
```conda install -c bioconda gatk4 picard raxml samtools bbmap dendropy minimap2 biopython```  
#You might need to install Bioconda with: pip install Biopython  

2. Download the wgfast github repository, install:  
```git clone https://github.com/jasonsahl/wgfast.git```  
```python setup.py build```  
```python setup.py install```

3. Open the script (wgfast.py) with a text editor and change the path to your *WG-FAST* installation directory.  
For example:  
WGFAST_PATH="/Users/jsahl/wgfast"  

4. To verify your installation, enter the wgfast directory and type the command below.  If everything
is working correctly, all tests should pass:  
```python tests/test_all_functions.py```  

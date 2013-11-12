#!/usr/bin/env python

"""Categorizes unknown SNPs.
Dependencies include:

1.  GATK - tested version is
2.  samtools - tested version is
3.  bwa - tested version is [optinal]
4.  novoalign - tested version is [optional]
4.  Picard tools - tested version is
5.  mothur -tested version is
6.  raxmlHPC - tested version is

Input is a SNP matrix: currently, this matrix
can be generated with NASP

"""

from optparse import OptionParser
import subprocess
import os
import sys
from wg_fast.util import *
import errno
from igs.utils import logging as log_isg

WGFAST_PATH="/Users/jsahl/wg_fast"
sys.exec.append("%s" % WGFAST_PATH)
GATK_PATH=WGFAST_PATH+"/bin/GenomeAnalysisTK.jar"
PICARD_PATH=WGFAST_PATH+"/bin/CreateSequenceDictionary.jar"

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastqs cannot be found"
        sys.exit()

def test_filter(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported.  Only select from T and F"
        sys.exit()

def main(aligner,matrix,tree,reference,directory,processors,coverage,proportion,keep,subsample,subnums):
    start_dir = os.getcwd()
    log_isg.logPrint('testing the paths of all dependencies')
    ap=os.path.abspath("%s" % start_dir)
    """need to check the paths for all dependencies"""
    aa = subprocess.call(['which', 'raxmlHPC-PTHREADS'])
    if aa == 0:
        pass
    else:
        print "RAxML must be in your path as raxmlHPC-PTHREADS"
        sys.exit()
    ab = subprocess.call(['which', 'samtools'])
    if ab == 0:
        pass
    else:
        print "samtools must be in your path"
        sys.exit()
    ac = subprocess.call(['which', 'bwa'])
    if ac == 0:
        pass
    else:
        print "bwa must be in your path"
        sys.exit()
    ad = subprocess.call(['which', 'mothur'])
    if ad == 0:
        pass
    else:
        print "mothur must be in your path"
        sys.exit()
    """done checking for dependencies"""
    ref_path=os.path.abspath("%s" % reference)
    dir_path=os.path.abspath("%s" % directory)
    try:
        os.makedirs('%s/scratch' % ap)
    except OSError, e:
        if e.errno != errno.EEXIST:raise
    """copy stuff into the scratch directory"""
    subprocess.check_call("cp %s %s/scratch/reference.fasta" % (ref_path, ap), shell=True)
    """index reference file.  GATK appears to do this incorrectly"""
    subprocess.check_call("samtools faidx %s/scratch/reference.fasta" % (ap), shell=True)
    subprocess.check_call("bwa index %s/scratch/reference.fasta > /dev/null 2>&1" % (ap), shell=True)
    ref_name=get_seq_name(reference)
    reduced=ref_name.replace(".fasta","")
    """creates dict file with picard tools.  In testing, GATK does this incorrectly"""
    try:
        os.system("java -jar %s R=%s/scratch/reference.fasta O=%s/scratch/reference.dict > /dev/null 2>&1" % (PICARD_PATH, ap, ap))
    except:
        print "dict wasn't created"
    fileSets=read_file_sets(dir_path)
    ref_coords = grab_matrix_coords(matrix)
    run_loop(fileSets, dir_path,"%s/scratch/reference.fasta" % ap , processors, GATK_PATH, ref_coords, coverage, proportion)
    """will subsample based on the number of SNPs reported by the following function"""
    used_snps=find_used_snps()
    outnames=merge_vcfs(matrix)
    for name in outnames:
        for k,v in used_snps.iteritems():
            if name==k:
                log_isg.logPrint("number of usable SNPs in genome %s = %s" % (k,v))
    subprocess.check_call("paste ref.list *.tmp.matrix > merged.vcf", shell=True)
    subprocess.check_call("rm -rf ref.list *.tmp.matrix", shell=True)
    merge_matrix(matrix, "merged.vcf")
    matrix_to_fasta("combined.matrix")
    os.system("mv combined.matrix %s/nasp_matrix.with_unknowns.txt" % ap)
    true_dists=dist_seqs("all.fasta", outnames)
    """add in code for the use of multiple aligners"""
    #if "bwa" == aligner:
    #    subprocess.check_call("bwa index %s > /dev/null 2>&1" % ref_path, shell=True)
    #    run_bwa(reference, read1, read2, processors, name)
    #elif "novoalign" == aligner:
    #    """stuff for novoalign"""
    #else:
    #    print "invalid aligner option; choose from bwa or novoalign"   
    """end of code to add in"""
    os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
    try:
        run_raxml("out.fasta", tree, processors)
    except:
        print "raxml encountered an error and unknown couldn't be added"
    if subsample=="T":
        dist_sets=find_two()
        log_isg.logPrint("running subsample routine")
        subsample_snps("nasp_matrix.with_unknowns.txt", dist_sets, used_snps,subnums)
        process_temp_matrices()
        compare_subsample_results(true_dists)
    else:
        log_isg.logPrint("all done")
    if keep == "T":
        pass
    else:
        subprocess.check_call("rm all.dist all.fasta mothur.* raxml.log raxml.out merged.vcf out.fasta* *tmp.matrix renamed.dist", shell=True)
        for outname in outnames:
            subprocess.check_call("rm %s.bam* %s.vcf* %s.filtered.vcf %s.sam.log" % (outname,outname,outname,outname), shell=True)
            os.chdir("%s" % ap)
            subprocess.check_call("rm -rf scratch", shell=True)
    log_isg.logPrint("all done")
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-a", "--aligner", dest="aligner",
                      help="aligner to use, either bwa or novoalign, defaults to bwa",
                      action="store", type="string", default="bwa")
    parser.add_option("-m", "--snp_matrix", dest="matrix",
                      help="path to NASP snp_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-t", "--tree", dest="tree",
                      help="path to input tree [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-r", "--reference_fasta", dest="reference",
                      help="path to reference fasta [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-d", "--directory", dest="directory",
                      help="path to directory of fastq files [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-p", "--processors", dest="processors",
                      help="# of processors to use - defaults to 2",
                      default="2", type="int")
    parser.add_option("-c", "--coverage", dest="coverage",
		      help="minimum SNP coverage required to be called a SNP - defaults to 3",
                      default="3", type="int")
    parser.add_option("-o", "--proportion", dest="proportion",
		      help="proportion of alleles to be called a SNP, defaults to 0.9",
                      default="0.9", type="float")
    parser.add_option("-k", "--keep", dest="keep",
                      help="keep temp files?  Defaults to F",
                      action="callback", callback=test_filter, type="string", default="F")
    parser.add_option("-s", "--subsample", dest="subsample",
                      help="Run subsample routine? Defaults to T",
                      action="callback", callback=test_filter, type="string", default="T")
    parser.add_option("-n", "--subnums", dest="subnums",
                      help="number of subsamples to process, defaults to 100",
                      action="store", type="int", default="100")

    options, args = parser.parse_args()
    
    mandatories = ["matrix", "tree", "reference", "directory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.aligner,options.matrix,options.tree,options.reference,options.directory,
         options.processors,options.coverage,options.proportion,options.keep,options.subsample,
         options.subnums)

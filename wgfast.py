#!/usr/bin/env python

"""WG-FAST
written by Jason Sahl
correspondence: jasonsahl@gmail.com

Dependencies include:

1.  GATK - tested version is 2.7.2.  This version
    requires Java 1.7.
2.  samtools - tested version is 0.1.19
3.  bwa - tested version is 0.7.5a
4.  Picard tools - tested version is 1.79
5.  raxmlHPC - tested version is 7.7.8
6.  DendroPy - only required for sub-sampling
    routine.  Tested version is 3.12.0

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

#WGFAST_PATH="/home/jsahl/tools/wgfast"
WGFAST_PATH="/Users/jsahl/wgfast"
sys.path.append("%s" % WGFAST_PATH)
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

def main(matrix,tree,reference,directory,processors,coverage,proportion,keep,subsample,subnums):
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
    print "citation: 'Stamatakis A. RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics. 2006;22(21):2688-90'"
    print "citation: 'Berger SA, Krompass D, Stamatakis A. Performance, accuracy, and Web server for evolutionary placement of short sequence reads under maximum likelihood. Syst Biol. 2011;60(3):291-302'"
    ab = subprocess.call(['which', 'samtools'])
    if ab == 0:
        pass
    else:
        print "samtools must be in your path"
        sys.exit()
    print "citation: 'Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, Genome Project Data Processing S. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-9'"
    ac = subprocess.call(['which', 'bwa'])
    if ac == 0:
        pass
    else:
        print "bwa must be in your path"
        sys.exit()
    print "citation: 'Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXivorg. 2013(arXiv:1303.3997 [q-bio.GN])'"
    print "Also uses GATK"
    print "citation: 'McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research. 2010;20(9):1297-303'"
    #done checking for dependencies
    print ""
    log_isg.logPrint('WG-FAST pipeline starting')
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
    """write reduced matrix with only the SNP data"""
    write_reduced_matrix(matrix)
    ref_name=get_seq_name(reference)
    reduced=ref_name.replace(".fasta","")
    """creates dict file with picard tools.  In testing, GATK does this incorrectly"""
    try:
        os.system("java -jar %s R=%s/scratch/reference.fasta O=%s/scratch/reference.dict > /dev/null 2>&1" % (PICARD_PATH, ap, ap))
    except:
        print "dict wasn't created"
    fileSets=read_file_sets(dir_path)
    ref_coords = grab_matrix_coords(matrix)
    run_loop(fileSets, dir_path,"%s/scratch/reference.fasta" % ap , processors, GATK_PATH, ref_coords, coverage, proportion, matrix, ap)
    """will subsample based on the number of SNPs reported by the following function"""
    used_snps=find_used_snps()
    outnames=grab_names()
    for name in outnames:
        for k,v in used_snps.iteritems():
            if name==k:
                log_isg.logPrint("number of callable positions in genome %s = %s" % (k,v))
    subprocess.check_call("paste *.tmp.matrix > merged.vcf", shell=True)
    """deletes temporary files that could be confused later on"""
    subprocess.check_call("rm -rf *.tmp.matrix", shell=True)
    subprocess.check_call("paste temp.matrix merged.vcf > combined.matrix", shell=True)
    matrix_to_fasta("combined.matrix")
    os.system("mv combined.matrix %s/nasp_matrix.with_unknowns.txt" % ap)
    """this fixes the SNP output to conform with RAxML"""
    os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
    try:
        run_raxml("out.fasta", tree, processors, "classification_results.txt")
        transform_tree("tree_including_unknowns_noedges.tree")
        print ""
        log_isg.logPrint("Insertion likelihood values:")
        parse_likelihoods("classification_results.txt")
        print ""
    except:
        print "raxml encountered an error and unknown couldn't be added"
    calculate_pairwise_tree_dists("tree_including_unknowns_noedges.tree","all_patristic_distances.txt")
    if subsample=="T":
        """find true distance of each genome to the reference"""
        os.system("sort -g -k 6 all_patristic_distances.txt > tmp_patristic_distances.txt")
        final_sets, distances=find_two("tmp_patristic_distances.txt", outnames)
        get_closest_dists(final_sets, distances, outnames)
        log_isg.logPrint("running subsample routine")
        subsample_snps("nasp_matrix.with_unknowns.txt", final_sets, used_snps, subnums)
        os.system("sed 's/QUERY___//g' tree_including_unknowns_noedges.tree > tmp.tree")
        process_temp_matrices(final_sets, "tmp.tree", processors, "all_patristic_distances.txt")
        compare_subsample_results(outnames)
    else:
        log_isg.logPrint("all done")
    if keep == "T":
        pass
    else:
        try:
            subprocess.check_call("rm all.dist all.fasta mothur.* raxml.log raxml.out merged.vcf out.fasta* *tmp.matrix renamed.dist temp.matrix tmp.tree tmp_patristic_distances.txt", shell=True, stderr=open(os.devnull, 'w'))
        except:
            pass
        for outname in outnames:
            try:
                subprocess.check_call("rm %s.bam* %s.vcf* %s.filtered.vcf %s.sam.log %s.closest.two.txt %s_coverage*" % (outname,outname,outname,outname,outname,outname), shell=True, stderr=open(os.devnull, 'w'))
            except:
                pass
            os.chdir("%s" % ap)
            subprocess.check_call("rm -rf scratch", shell=True)
    log_isg.logPrint("all done")
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
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

    main(options.matrix,options.tree,options.reference,options.directory,
         options.processors,options.coverage,options.proportion,options.keep,options.subsample,
         options.subnums)

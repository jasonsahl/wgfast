#!/usr/bin/env python

"""

WG-FAST
written by Jason Sahl
correspondence: jasonsahl@gmail.com

"""

from optparse import OptionParser
import subprocess
import os
import sys
import errno
import glob
try:
    from wg_fast.util import *
    from igs.utils import logging as log_isg
except:
    print "wgfast path needs to be modified in the wgfast.py file"
    sys.exit()

"""modify line below to reflect your installation directory"""
WGFAST_PATH="/Users/jasonsahl/tools/wgfast"
if os.path.exists(WGFAST_PATH):
    sys.path.append("%s" % WGFAST_PATH)
else:
    print "your WGFAST path is not correct.  Edit the path in wgfast.py and try again"
    sys.exit()
GATK_PATH=WGFAST_PATH+"/bin/GenomeAnalysisTK.jar"
PICARD_PATH=WGFAST_PATH+"/bin/CreateSequenceDictionary.jar"
ADD_GROUPS=WGFAST_PATH+"/bin/AddOrReplaceReadGroups.jar"
TRIM_PATH=WGFAST_PATH+"/bin/trimmomatic-0.30.jar"
    
def main(matrix,tree,reference,directory,parameters,processors,coverage,proportion,keep,subsample,subnums,doc,tmp_dir,insertion_method,fudge,only_subs,model):
    #start_dir = os.getcwd()
    ref_path=os.path.abspath("%s" % reference)
    dir_path=os.path.abspath("%s" % directory)
    #check for binary dependencies
    log_isg.logPrint('testing the paths of all dependencies')
    ap=os.path.abspath("%s" % os.getcwd())
    aa = subprocess.call(['which', 'raxmlHPC-SSE3'])
    if aa == 0:
        pass
    else:
        print "RAxML must be in your path as raxmlHPC-PTHREADS"
        sys.exit()
    print "*citation: 'Stamatakis A. RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics. 2006;22(21):2688-90'"
    print "*citation: 'Berger SA, Krompass D, Stamatakis A. Performance, accuracy, and Web server for evolutionary placement of short sequence reads under maximum likelihood. Syst Biol. 2011;60(3):291-302'"
    ab = subprocess.call(['which', 'samtools'])
    if ab == 0:
        pass
    else:
        print "samtools must be in your path"
        sys.exit()
    print "*citation: 'Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, Genome Project Data Processing S. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-9'"
    ac = subprocess.call(['which', 'bwa'])
    if ac == 0:
        pass
    else:
        print "bwa must be in your path"
        sys.exit()
    print "*citation: 'Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXivorg. 2013(arXiv:1303.3997 [q-bio.GN])'"
    print "Patristic distances calculated with DendroPy"
    print "*citation: 'Sukumaran J, Holder MT. DendroPy: a Python library for phylogenetic computing. Bioinformatics. 2010;26(12):1569-71. Epub 2010/04/28. doi: 10.1093/bioinformatics/btq228. PubMed PMID: 20421198'"
    print "Also uses GATK for variant calling"
    print "*citation: 'McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research. 2010;20(9):1297-303'"
    print "Uses trimmomatic for read trimming"
    print "*citation: Bolger A.M., Lohse M., Usadel B. Trimmomatic: A flexible trimmer for Illumina Sequence Data.  Bioinformatics. 2014.  Doi:10.1093/bioinformatics/btu170"
    print ""
    #done checking for dependencies"""
    log_isg.logPrint('WG-FAST pipeline starting')
    try:
        os.makedirs('%s/scratch' % ap)
    except OSError, e:
        if e.errno != errno.EEXIST:raise
    #copy reference into the scratch directory, where all of the work will take place
    subprocess.check_call("cp %s %s/scratch/reference.fasta" % (ref_path, ap), shell=True)
    #index reference file.  GATK appears to do this incorrectly"""
    subprocess.check_call("samtools faidx %s/scratch/reference.fasta" % (ap), shell=True)
    subprocess.check_call("bwa index %s/scratch/reference.fasta > /dev/null 2>&1" % (ap), shell=True)
    #write reduced matrix with only the SNP data"""
    write_reduced_matrix(matrix)
    ref_name=get_seq_name(reference)
    reduced=ref_name.replace(".fasta","")
    """creates dict file with picard tools.  In testing, GATK does this incorrectly"""
    try:
        os.system("java -jar %s R=%s/scratch/reference.fasta O=%s/scratch/reference.dict > /dev/null 2>&1" % (PICARD_PATH, ap, ap))
    except:
        print "dict wasn't created"
        sys.exit()
    if only_subs == "T":
        pass
    else:
        fileSets=read_file_sets(dir_path)
        ref_coords = get_all_snps(matrix)
        run_loop(fileSets, dir_path,"%s/scratch/reference.fasta" % ap , processors, GATK_PATH, ref_coords, coverage, proportion, matrix, ap,doc,tmp_dir,ADD_GROUPS,TRIM_PATH,WGFAST_PATH)
    """will subsample based on the number of SNPs reported by the following function"""
    used_snps=find_used_snps()
    outnames=grab_names()
    for name in outnames:
        for k,v in used_snps.iteritems():
            if name==k:
                log_isg.logPrint("number of callable positions in genome %s = %s" % (k,v))
    if only_subs == "T":
        try:
            os.system("rm RAxML*")
        except:
            pass
        pass
    else:
        """need to change this to a python function"""
        subprocess.check_call("paste *.tmp.matrix > merged.vcf", shell=True)
        """deletes temporary files that could be confused later on"""
        subprocess.check_call("rm -rf *.tmp.matrix", shell=True)
        subprocess.check_call("paste temp.matrix merged.vcf > combined.matrix", shell=True)
        matrix_to_fasta("combined.matrix", "all.fasta")
        os.system("mv combined.matrix %s/nasp_matrix.with_unknowns.txt" % ap)
        """this fixes the SNP output to conform with RAxML"""
        os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
        if insertion_method == "ML":
            suffix = run_raxml("out.fasta", tree,"out.classification_results.txt", "V", parameters, model, "out")
        elif insertion_method == "MP":
            try:
                suffix = run_raxml("out.fasta", tree,"out.classification_results.txt", "y", parameters, model, "out")
            except:
                print "problem with MP, moving to ML"
                suffix = run_raxml("out.fasta", tree, "out.classification_results.txt", "V", parameters, model, "out")
        else:
            pass
        transform_tree("%s.tree_including_unknowns_noedges.tree" % suffix)
        print ""
        if insertion_method == "ML":
            log_isg.logPrint("Insertion likelihood values:")
            parse_likelihoods("out.classification_results.txt")
            print ""
        else:
            pass
        if insertion_method == "ML":
            calculate_pairwise_tree_dists("%s.tree_including_unknowns_noedges.tree" % suffix,"all_patristic_distances.txt")
        else:
            print "patrisitic distances can't always be calculated with parsimony"
            pass
    if subsample=="T":
        aa = subprocess.call(['which', 'raxmlHPC-PTHREADS-SSE3'])
        if aa == 0:
            pass
        else:
            print "for sub-sample routine, RAxML must be in your path as raxmlHPC-PTHREADS-SSE3"
            sys.exit()
        print "*citation: 'Stamatakis A. RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics. 2006;22(21):2688-90'"
        if insertion_method == "MP":
            print "subsample method is not compatible with MP"
            pass
        else:
            try:
                os.system("sort -g -k 6 all_patristic_distances.txt | sed 's/://g' > tmp_patristic_distances.txt")
            except:
                print "all_patrisitic_distances.txt must be in your analysis directory!"
                sys.exit()
            final_sets, distances=find_two_new("tmp_patristic_distances.txt", outnames)
            results = get_closest_dists_new(final_sets,outnames)
            log_isg.logPrint("running subsample routine, forcing GTRGAMMA model")
            thread_list = []
            files_and_temp_names = [(list(f)) for f in final_sets]
            allsnps = get_all_snps(matrix)
            def _perform_workflow(data):
                f = data
                subsample_snps_dev("temp.matrix", f, used_snps, subnums, allsnps)
            results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))
            temp_matrices = glob.glob(os.path.join(ap, "*tmp.matrix"))
            final_matrices = []
            for matrix in temp_matrices:
                final_matrices.append(matrix.replace("%s/" % ap,""))
            sample_sets = {}
            for matrix in final_matrices:
                entries = matrix.split(".")
                if entries[0] in sample_sets:
                    sample_sets[entries[0]].append(entries[2])
                else:
                    sample_sets[entries[0]]=[entries[2]]
            log_isg.logPrint('creating PARAMS file')
            for k,v in sample_sets.iteritems():
                uniques= []
                [uniques.append(item) for item in v if item not in uniques]
                def _perform_workflow(data):
                    create_params_files(k, uniques, tree, "temp.matrix", final_sets, processors)
                results = set(p_func.pmap(_perform_workflow,
                                      sample_sets,
                                      num_workers=processors))
            log_isg.logPrint('adding unknowns to tree')
            def _perform_workflow_2(data):
                z = data
                process_temp_matrices_dev(final_sets, z, tree, processors, "all_patristic_distances.txt", "V", parameters, model)
            results = set(p_func.pmap(_perform_workflow_2,
                                      final_matrices,
                                      num_workers=processors))
            compare_subsample_results(outnames,distances,fudge)
    else:
        pass
    if keep == "T":
        pass
    else:
        try:
            subprocess.check_call("rm all.dist all.fasta raxml.log raxml.out merged.vcf out.fasta* *tmp.matrix renamed.dist temp.matrix tmp.tree tmp_patristic_distances.txt", shell=True, stderr=open(os.devnull, 'w'))
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
    parser.add_option("-x", "--parameters_file", dest="parameters",
                      help="parameters for RAxML insertion, defaults to NULL",
                      action="store", type="string", default="NULL")
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
    parser.add_option("-g", "--doc", dest="doc",
                      help="run depth of coverage on all files?  Defaults to T",
                      action="callback", callback=test_filter, type="string", default="T")
    parser.add_option("-e", "--temp", dest="tmp_dir",
                      help="temporary directory for GATK analysis, defaults to /tmp",
                      action="store", type="string", default="/tmp")
    parser.add_option("-z", "--method", dest="insertion_method",
                      help="method to insert unknown genomes (MP or ML), defaults to ML",
                      action="callback", type="string", callback=test_methods, default="ML")
    parser.add_option("-f", "--fudge_factor", dest="fudge",
                      help="How close does a subsample have to be from true placement?  Defaults to 0.1",
                      action="store", type="float", default="0.1")
    parser.add_option("-y", "--only_subs", dest="only_subs",
                      help="Only run sub-sample routine and exit? Defaults to F",
                      action="callback", callback=test_filter, type="string", default="F")
    parser.add_option("-j", "--model", dest="model",
                      help="which model to run with raxml, GTRGAMMA, ASC_GTRGAMMA, GTRCAT, ASC_GTRCAT",
                      action="callback", callback=test_models, type="string", default="ASC_GTRGAMMA")

    options, args = parser.parse_args()
    
    mandatories = ["matrix", "tree", "reference", "directory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix,options.tree,options.reference,options.directory,options.parameters,
         options.processors,options.coverage,options.proportion,options.keep,options.subsample,
         options.subnums,options.doc,options.tmp_dir,options.insertion_method,options.fudge,
         options.only_subs,options.model)

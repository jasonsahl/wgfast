#!/usr/bin/env python

"""
WG-FAST
written by Jason Sahl
correspondence: jasonsahl@gmail.com
"""
from __future__ import print_function
from optparse import OptionParser
import subprocess
import os
import sys
import errno
import glob

"""modify line below to reflect your installation directory"""
WGFAST_PATH="/Users/jasonsahl/tools/wgfast"

if os.path.exists(WGFAST_PATH):
    sys.path.append("%s" % WGFAST_PATH)
else:
    print("your WG-FAST path is not correct.  Edit the path in wgfast.py and try again")
    sys.exit()
try:
    from wg_fast.util import *
except:
    print("wgfast path needs to be modified in the wgfast.py file")
    sys.exit()

#Need to update with new GATK4
"""I should no longer need this argument"""
#GATK_PATH=WGFAST_PATH+"/bin/GenomeAnalysisTK.jar"
#TODO: See if I can remove this
#PICARD_PATH=WGFAST_PATH+"/bin/CreateSequenceDictionary.jar"
ADD_GROUPS=WGFAST_PATH+"/bin/AddOrReplaceReadGroups.jar"

def main(reference_dir,read_directory,processors,coverage,proportion,keep,subsample,
    subnums,doc,tmp_dir,fudge,only_subs,model,gatk_method):
    ref_path=os.path.abspath("%s" % reference_dir)
    dir_path=os.path.abspath("%s" % read_directory)
    """Test to make sure all required files are present"""
    tree = "".join(glob.glob(os.path.join(ref_path, "*.tree")))
    tree_list = glob.glob(os.path.join(ref_path, "*.tree"))
    if len(tree_list) == 0:
        print("You need to provide a tree in your reference directory ending in '.tree'")
        sys.exit()
    elif len(tree_list) == 1:
        pass
    elif len(tree_list) > 1:
        print("More than one tree ending in '.tree' found in your reference directory...exiting")
        sys.exit()
    matrix = "".join(glob.glob(os.path.join(ref_path, "*.tsv")))
    matrix_list = glob.glob(os.path.join(ref_path, "*.tsv"))
    if len(matrix_list) == 0:
        print("You need to provide a NASP formatted matrix ending in '.tsv'")
        sys.exit()
    elif len(matrix_list) > 1:
        print("More than one file in your reference directory ending in '.tsv'...exiting")
        sys.exit()
    reference = "".join(glob.glob(os.path.join(ref_path, "*.fasta")))
    reference_list = glob.glob(os.path.join(ref_path, "*.fasta"))
    if len(reference_list) == 0:
        print("You must provide a REFERENCE FASTA file")
        sys.exit()
    elif len(reference_list) > 1:
        print("More than one file in your reference directory ending in '.fasta'...exiting")
        sys.exit()
    """Get the reference information. Outfile is named $name.tmp.txt"""
    get_seq_length(reference, "ref")
    subprocess.check_call("tr ' ' '\t' < ref.tmp.txt > ref.genome_size.txt", shell=True)
    if len(reference_list) == 0:
        print("You need to provide a reference FASTA file in your reference directory, ending in 'fasta'")
        sys.exit()
    elif len(reference_list) > 1:
        print("Multiple files found in your reference directory ending in 'fasta'...exiting")
        sys.exit()
    try:
        parameters = "".join(glob.glob(os.path.join(ref_path, "*.PARAMS")))
        parameters_list = glob.glob(os.path.join(ref_path, "*.PARAMS"))
        if len(parameters_list)>1:
            print("More than one RAxML parameters file found in your reference directory...exiting")
            sys.exit()
        elif len(parameters_list) == 0:
            parameters = "NULL"
    except:
        parameters = "NULL"
    #check for binary dependencies
    log_isg.logPrint('testing the paths of all dependencies')
    ap=os.path.abspath("%s" % os.getcwd())
    aa = subprocess.call(['which', 'raxmlHPC-SSE3'])
    if aa == 0:
        pass
    else:
        print("RAxML must be in your path as raxmlHPC-SSE3")
        sys.exit()
    print("*citation: 'Stamatakis, A. RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics (2014).'")
    print("*citation: 'Berger SA, Krompass D, Stamatakis A. Performance, accuracy, and Web server for evolutionary placement of short sequence reads under maximum likelihood. Syst Biol. 2011;60(3):291-302'")
    ab = subprocess.call(['which', 'samtools'])
    if ab == 0:
        pass
    else:
        print("samtools must be in your path")
        sys.exit()
    print("*citation: 'Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, Genome Project Data Processing S. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-9'")
    ac = subprocess.call(['which', 'bwa'])
    if ac == 0:
        pass
    else:
        print("bwa must be in your path")
        sys.exit()
    """This is the new test for bbduk.sh"""
    ac = subprocess.call(['which', 'bbduk.sh'])
    if ac == 0:
        pass
    else:
        print("bbduk needs to be in your path as bbduk.sh")
        sys.exit()
    ac = subprocess.call(['which', 'gatk'])
    if ac == 0:
        pass
    else:
        print("GATK needs to be in your path as gatk")
        sys.exit()
    print("*citation: 'Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXivorg. 2013(arXiv:1303.3997 [q-bio.GN])'")
    print("Patristic distances calculated with DendroPy")
    print("*citation: 'Sukumaran J, Holder MT. DendroPy: a Python library for phylogenetic computing. Bioinformatics. 2010;26(12):1569-71. Epub 2010/04/28. doi: 10.1093/bioinformatics/btq228. PubMed PMID: 20421198'")
    print("Also uses GATK for variant calling")
    print("*citation: 'McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research. 2010;20(9):1297-303'")
    print("Uses BioPython for FASTA parsing")
    print("*citation :Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3")
    print("")
    #done checking for dependencies"""
    log_isg.logPrint('WG-FAST pipeline starting')
    log_isg.logPrint("WG-FAST was invoked with the following parameters:")
    print("-m %s \\" % "".join(matrix))
    print("-t %s \\" % "".join(tree))
    print("-r %s \\" % "".join(reference))
    print("-d %s \\" % read_directory)
    print("-x %s \\" % "".join(parameters))
    print("-p %s \\" % processors)
    print("-c %s \\" % coverage)
    print("-o %s \\" % proportion)
    print("-k %s \\" % keep)
    print("-s %s \\" % subsample)
    print("-n %s \\" % subnums)
    print("-g %s \\" % doc)
    print("-e %s \\" % tmp_dir)
    print("-f %s \\" % fudge)
    print("-y %s \\" % only_subs)
    print("-j %s \\" % model)
    print("-q %s" % gatk_method)
    print("-------------------------")
    try:
        os.makedirs('%s/scratch' % ap)
    except:
        os.system("rm -rf %s/scratch" % ap)
        os.makedirs('%s/scratch' % ap)
    check_input_files(matrix, reference)
    ########Real work starts here############
    #copy reference into the scratch directory, where all of the work will take place
    if only_subs == "T":
        pass
    else:
        subprocess.check_call("cp %s %s/scratch/reference.fasta" % (reference, ap), shell=True)
        #index reference file.  GATK appears to do this incorrectly"""
        subprocess.check_call("samtools faidx %s/scratch/reference.fasta" % ap, shell=True)
        #TODO: replace with minimap2
        subprocess.check_call("bwa index %s/scratch/reference.fasta > /dev/null 2>&1" % ap, shell=True)

        """creates dict file with picard tools.  In testing, GATK does this incorrectly"""
        try:
            #TODO: Check to see if I can remove this step
            #os.system("java -jar %s R=%s/scratch/reference.fasta O=%s/scratch/reference.dict > /dev/null 2>&1" % (PICARD_PATH, ap, ap))
            os.system("gatk CreateSequenceDictionary -R=%s/scratch/reference.fasta -O=%s/scratch/reference.dict > /dev/null 2>&1" % (ap, ap))
        except:
            print("dict wasn't created")
            sys.exit()
    if os.path.isfile("temp.matrix"):
        pass
    else:
        #write reduced matrix with only the SNP data"""
        write_reduced_matrix(matrix)
    ref_name=get_seq_name(reference)
    if only_subs == "T":
        pass
    else:
        fileSets=read_file_sets(dir_path)
        if len(fileSets) == 0:
            print("No usable file sets found...exiting")
            try:
                os.system("rm temp.matrix")
            except:
                pass
            sys.exit()
        else:
            ref_coords = get_all_snps(matrix)
            log_isg.logPrint("Loop starting")
            #run_loop_dev(fileSets,dir_path,"%s/scratch/reference.fasta" % ap,processors,GATK_PATH,
            #ref_coords,coverage,proportion,matrix,ap,doc,tmp_dir,ADD_GROUPS,WGFAST_PATH,gatk_method)
            run_loop_dev(fileSets,dir_path,"%s/scratch/reference.fasta" % ap,processors,ref_coords,
            coverage,proportion,matrix,ap,doc,tmp_dir,WGFAST_PATH,gatk_method)
    """will subsample based on the number of SNPs reported by the following function"""
    if "T" in doc:
        os.system("cat *breadth.txt > breadth_over_%sx_out.txt" % coverage)
        os.system("cat *sum_cov.txt> coverage_out.txt")
    else:
        pass
    used_snps=find_used_snps()
    #Outnames is required for the sub-sampling routine, even with -y T
    outnames=grab_names()
    for name in outnames:
        for k,v in used_snps.items():
            if name == k:
                log_isg.logPrint("number of callable positions in genome %s = %s" % (k,v))
    if only_subs == "T":
        try:
            #Starts with a clean slate, to replace with new EPA algorithm
            os.system("rm RAxML*")
        except:
            pass
        pass
    else:
        create_merged_vcf()
        subprocess.check_call("paste temp.matrix merged.vcf > combined.matrix", shell=True)
        matrix_to_fasta("combined.matrix", "all.fasta")
        os.system("mv combined.matrix %s/nasp_matrix.with_unknowns.txt" % ap)
        """this fixes the SNP output to conform with RAxML"""
        os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
        suffix = run_raxml("out.fasta", tree,"out.classification_results.txt", "V", parameters, model, "out")
        transform_tree("%s.tree_including_unknowns_noedges.tree" % suffix)
        print("")
        log_isg.logPrint("Insertion likelihood values:")
        parse_likelihoods("out.classification_results.txt")
        print("")
        calculate_pairwise_tree_dists("%s.tree_including_unknowns_noedges.tree" % suffix,"all_patristic_distances.txt")
    if subsample=="T":
        aa = subprocess.call(['which', 'raxmlHPC-PTHREADS-SSE3'])
        if aa == 0:
            pass
        else:
            print("for sub-sample routine, RAxML must be in your path as raxmlHPC-PTHREADS-SSE3")
            sys.exit()
        print("*citation: 'Stamatakis A. RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics. 2006;22(21):2688-90'")
        try:
            if os.path.isfile("tmp_patristic_distances.txt"):
                pass
            else:
                os.system("sort -g -k 6 all_patristic_distances.txt | sed 's/://g' > tmp_patristic_distances.txt")
        except:
            print("all_patrisitic_distances.txt must be in your analysis directory!")
            sys.exit()
        final_sets, distances=find_two_new("tmp_patristic_distances.txt", outnames)
        results = get_closest_dists_new(final_sets,outnames)
        log_isg.logPrint("running subsample routine, forcing GTRGAMMA model")
        """mpshell on this function"""
        allsnps = get_all_snps(matrix)
        subsample_snps_2(final_sets,used_snps,subnums,allsnps,processors,"temp.matrix")
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
        new_sample_dicts = {}
        for k,v in sample_sets.iteritems():
            uniques = []
            [uniques.append(item) for item in v if item not in uniques]
            new_sample_dicts.update({k:uniques})
        log_isg.logPrint('creating PARAMS file')
        if os.path.isfile("*PARAMS"):
            pass
        else:
            create_params_files_dev(new_sample_dicts,tree,"temp.matrix",final_sets,processors)
        try:
            """Must make sure that remove previous RAxML files"""
            subprocess.check_call("rm RAxML*", shell=True, stderr=open(os.devnull, 'w'))
        except:
            pass
        """final_matrices does indeed have all of the temp matrices loaded"""
        log_isg.logPrint('adding unknowns to tree')
        #TODO: enter progress bar instead of printing out lots of text
        process_temp_matrices_2(final_sets,final_matrices,tree,processors,"all_patristic_distances.txt", "V", parameters, model)
        print("-------------------------")
        compare_subsample_results(outnames,distances,fudge)
    else:
        pass
    #Clean up temporary files
    if keep == "T":
        pass
    else:
        try:
            subprocess.check_call("rm ref.tmp.txt ref.genome_size.txt all.dist all.fasta raxml.log raxml.out merged.vcf out.fasta* *tmp.matrix renamed.dist tmp.tree temp.matrix tmp_patristic_distances.txt out* RAxML_portableTree*jplace *.unpaired.fastq.gz", shell=True, stderr=open(os.devnull, 'w'))
        except:
            pass
        for outname in outnames:
            try:
                subprocess.check_call("rm %s* RAxML_log* RAxML_info*" % outname, shell=True, stderr=open(os.devnull, 'w'))
            except:
                pass
            os.chdir("%s" % ap)
            subprocess.check_call("rm -rf scratch", shell=True)
    log_isg.logPrint("all done")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", "--reference_directory", dest="reference_dir",
                      help="path to reference file directory [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-d", "--read_directory", dest="read_directory",
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
    parser.add_option("-g", "--doc", dest="doc",
                      help="run depth of coverage on all files?  Defaults to T",
                      action="callback", callback=test_filter, type="string", default="T")
    #TODO: think about removing this argument
    parser.add_option("-e", "--temp", dest="tmp_dir",
                      help="temporary directory for GATK analysis, defaults to /tmp",
                      action="store", type="string", default="/tmp")
    parser.add_option("-f", "--fudge_factor", dest="fudge",
                      help="How close does a subsample have to be from true placement?  Defaults to 0.1",
                      action="store", type="float", default="0.1")
    parser.add_option("-y", "--only_subs", dest="only_subs",
                      help="Only run sub-sample routine and exit? Defaults to F",
                      action="callback", callback=test_filter, type="string", default="F")
    parser.add_option("-j", "--model", dest="model",
                      help="which model to run with raxml, GTRGAMMA, ASC_GTRGAMMA, GTRCAT, ASC_GTRCAT",
                      action="callback", callback=test_models, type="string", default="ASC_GTRGAMMA")
    parser.add_option("-q", "--gatk_method", dest="gatk_method",
                      help="How to call GATK? Defaults to EMIT_ALL_CONFIDENT_SITES, can be EMIT_ALL_SITES",
                      action="callback", callback=test_gatk, type="string", default="EMIT_ALL_CONFIDENT_SITES")

    options, args = parser.parse_args()

    mandatories = ["reference_dir","read_directory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.reference_dir,options.read_directory,
         options.processors,options.coverage,options.proportion,options.keep,options.subsample,
         options.subnums,options.doc,options.tmp_dir,options.fudge,
         options.only_subs,options.model,options.gatk_method)

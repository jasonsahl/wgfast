import os
import re
import logging
import glob
import subprocess
from subprocess import Popen
from Bio import SeqIO
from Bio import Phylo
import glob
from igs.threading import functional as p_func
from igs.utils import logging as log_isg
from operator import itemgetter
import random

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def get_readFile_components(full_file_path):
    """function adapted from:
    https://github.com/katholt/srst2"""
    (file_path,file_name) = os.path.split(full_file_path)
    m1 = re.match("(.*).gz",file_name)
    ext = ""
    if m1 != None:
        ext = ".gz"
        file_name = m1.groups()[0]
    (file_name_before_ext,ext2) = os.path.splitext(file_name)
    full_ext = ext2+ext
    return(file_path,file_name_before_ext,full_ext)

def read_file_sets(dir_path):        
    """match up pairs of sequence data, adapted from
    https://github.com/katholt/srst2"""
    fileSets = {} 
    forward_reads = {}
    reverse_reads = {} 
    num_paired_readsets = 0
    num_single_readsets = 0
    for infile in glob.glob(os.path.join(dir_path, "*.fastq.gz")):
        (file_path,file_name_before_ext,full_ext) = get_readFile_components(infile)
        m=re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
        if m==None:
            m=re.match("(.*)("+"_R1"+")(_.*)$",file_name_before_ext)
            if m!=None:
                (baseName,read) = m.groups()[0], m.groups()[1]
                forward_reads[baseName] = infile
            else:
                m=re.match("(.*)("+"_R2"+")(_.*)$",file_name_before_ext)
                if m!=None:
                    (baseName,read) = m.groups()[0], m.groups()[1]
                    reverse_reads[baseName] = infile
                else:
                    print "Could not determine forward/reverse read status for input file "
        else:
            baseName, read  = m.groups()[0], m.groups()[3]
            if read == "_R1":
                forward_reads[baseName] = infile
            elif read == "_R2":
                reverse_reads[baseName] = infile
            else:
                print "Could not determine forward/reverse read status for input file "
                fileSets[file_name_before_ext] = infile
                num_single_readsets += 1
    for sample in forward_reads:
        if sample in reverse_reads:
            fileSets[sample] = [forward_reads[sample],reverse_reads[sample]] # store pair
            num_paired_readsets += 1
        else:
            fileSets[sample] = [forward_reads[sample]] # no reverse found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + forward_reads[sample])
    for sample in reverse_reads:
        if sample not in fileSets:
            fileSets[sample] = reverse_reads[sample] # no forward found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + reverse_reads[sample])
                                
    if num_paired_readsets > 0:
        logging.info('Total paired readsets found:' + str(num_paired_readsets))        
    if num_single_readsets > 0:
        logging.info('Total single reads found:' + str(num_single_readsets))

    return fileSets

def run_loop(fileSets, dir_path, reference, processors, gatk, ref_coords, coverage, proportion, matrix):
    files_and_temp_names = [(str(idx), list(f))
                            for idx, f in fileSets.iteritems()]
    def _perform_workflow(data):
        """idx is the sample name, f is the file dictionary"""
        idx, f = data
        run_bwa(reference, f[0], f[1], processors, idx)
        process_sam("%s.sam" % idx, idx)
        run_gatk(reference, processors, idx, gatk)
        #filter_vcf("%s.vcf.out" % idx, ref_coords, idx)
        #parse_vcf("%s.vcf.filtered" % idx, coverage, proportion, idx)
        process_vcf("%s.vcf.out" % idx, ref_coords, coverage, proportion, idx)
        make_temp_matrix("%s.filtered.vcf" % idx, matrix, idx)
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))

def bwa(reference,read1,read2,sam_file, processors, log_file='',**my_opts):
    """controller for bwa, currently only works with bwa mem"""
    mem_arguments = ['bwa', 'mem', '-v', '2', '-M', '-t', '%s' % processors]
    for opt in my_opts.items():
        mem_arguments.extend(opt)
    if "null" in read2:
        mem_arguments.extend([reference,read1])
    else:
        mem_arguments.extend([reference,read1,read2]) 
    if log_file:
       try:
           log_fh = open(log_file, 'w')
       except:
           print log_file, 'could not open'
    else:
        log_fh = PIPE
    try:
        sam_fh = open(sam_file, 'w')
    except:
        print sam_file, 'could not open'
 
    bwa = Popen(mem_arguments, stderr=log_fh, stdout=sam_fh)
    bwa.wait()

def run_bwa(reference, read1, read2, processors, name):
    """launces bwa. Adds in read_group for compatability with GATK"""
    read_group = '@RG\tID:%s\tSM:vac6wt\tPL:ILLUMINA\tPU:vac6wt' % name
    bwa(reference,read1, read2,"%s.sam" % name, processors, log_file='%s.sam.log' % name,**{'-R':read_group}) 

def process_sam(in_sam, name):
    subprocess.check_call("samtools view -h -b -S %s > %s.1.bam 2> /dev/null" % (in_sam, name), shell=True)
    subprocess.check_call("samtools view -u -h -F4 -o %s.2.bam %s.1.bam > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools view -h -b -q1 -F4 -o %s.3.bam %s.2.bam > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools sort %s.3.bam %s > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools index %s.bam > /dev/null 2>&1" % name, shell=True)
    subprocess.check_call("rm %s.1.bam %s.2.bam %s.3.bam %s" % (name, name, name, in_sam), shell=True)

def run_gatk(reference, processors, name, gatk):
    """gatk controller"""
    args = ['java', '-jar', '%s' % gatk, '-T', 'UnifiedGenotyper',
            '-R', '%s' % reference, '-nt', '%s' % processors, '-S', 'silent',
            '-mbq', '17', '-ploidy', '1', '-out_mode', 'EMIT_ALL_CONFIDENT_SITES',
            '-stand_call_conf', '30', '-stand_emit_conf', '30', '-I', '%s.bam' % name,
            '-rf', 'BadCigar']
    try:
        vcf_fh = open('%s.vcf.out' % name, 'w')
    except:
        log_isg.logPrint('could not open vcf file')
    try:
        log_fh = open('%s.vcf.log' % name, 'w')
    except:
        log_isg.logPrint('could not open log file')
    try:
        gatk_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
        gatk_run.wait()
    except:
        log_isg.logPrint("GATK encountered problems and did not run")

def filter_vcf(vcf, ref_coords, name):
    """filters a VCF for coordinates only present in the provided
    NASP SNP matrix"""
    vcf_in = open(vcf, "U")
    vcf_out = open("%s.vcf.filtered" % name, "w")
    outdata = [ ]
    for line in vcf_in:
        """skip VCF comment lines"""
        if line.startswith("#"):
            print >> vcf_out, line,
        else:
            fields=line.split()
            """skip INFO lines put in by GATK"""
            if "INFO" not in fields[0]:
                """must adapt for NASP format"""
                merged_fields=fields[0]+"::"+fields[1]
                try:
                    """checks to see what's in the NASP coord list"""
                    if merged_fields in ref_coords:
                        print >> vcf_out, line,
                        #outdata.append(merged_fields)
                    else:
                        pass
                except:
                    """if no coords are in your ref list, there is likely a problem"""
                    print "Are you sure you have the correct reference?"
                    sys.exit()
            else:
                continue
    vcf_in.close()
    vcf_out.close()
    return outdata
    
def parse_vcf(vcf, coverage, proportion, name):
    """finds SNPs that pass user-defined thresholds
    for coverage and proportion"""
    vcf_in = open(vcf, "U")
    vcf_out = open("%s.filtered.vcf" % name, "w")
    outdata = []
    for line in vcf_in:
        if line.startswith('#'):
           pass
        else:
            fields=line.split()
            """for GATK, a period signifies a reference call.
            First we want to look at the situation where this is
            not the case"""
            if "." != fields[4]:
                snp_fields=fields[9].split(':')
                if int(len(snp_fields))>2:
                    prop_fields=snp_fields[1].split(',')
                    if int(snp_fields[2])>=coverage:
                        if int(prop_fields[1])/int(snp_fields[2])>=float(proportion):
                            print >> vcf_out, fields[0]+"::"+fields[1],fields[4]+"\n",
                            outdata.append(fields[0]+"::"+fields[1]+"::"+fields[4])
                        else:
                            print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                            outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
                        """if problems are encountered, throw in a gap.  Could be too conservative"""
                    else:
                        print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                        outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
                else:
                    pass
            elif "." == fields[4]:
                nosnp_fields=fields[7].split(';')
                cov_fields=nosnp_fields[1].replace("DP=","")
                try:
                    if int(cov_fields)>=coverage:
                        print >> vcf_out, fields[0]+"::"+fields[1],fields[3]+"\n",
                        outdata.append(fields[0]+"::"+fields[1]+"::"+fields[3])
                    else:
                        print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                        outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
                except:
                     print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                     outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
            else:
                print "error in vcf file found!"
                sys.exit()
    vcf_in.close()
    vcf_out.close()
    return outdata

def process_vcf(vcf, ref_coords, coverage, proportion, name):
    """finds SNPs that pass user-defined thresholds
    for coverage and proportion"""
    vcf_in = open(vcf, "U")
    vcf_out = open("%s.filtered.vcf" % name, "w")
    outdata = []
    for line in vcf_in:
        if line.startswith('#'):
           pass
        elif line.startswith('INFO'):
            pass
        else:
            fields=line.split()
            """for GATK, a period signifies a reference call.
            First we want to look at the situation where this is
            not the case"""
            merged_fields=fields[0]+"::"+fields[1]
            if merged_fields in ref_coords:
                if "." != fields[4]:
                    snp_fields=fields[9].split(':')
                    if int(len(snp_fields))>2:
                        prop_fields=snp_fields[1].split(',')
                        if int(snp_fields[2])>=coverage:
                            if int(prop_fields[1])/int(snp_fields[2])>=float(proportion):
                                print >> vcf_out, fields[0]+"::"+fields[1],fields[4]+"\n",
                                outdata.append(fields[0]+"::"+fields[1]+"::"+fields[4])
                            else:
                                print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                                outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
                            """if problems are encountered, throw in a gap.  Could be too conservative"""
                        else:
                            print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                            outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
                    else:
                        pass
                elif "." == fields[4]:
                    nosnp_fields=fields[7].split(';')
                    cov_fields=nosnp_fields[1].replace("DP=","")
                    try:
                        if int(cov_fields)>=coverage:
                            print >> vcf_out, fields[0]+"::"+fields[1],fields[3]+"\n",
                            outdata.append(fields[0]+"::"+fields[1]+"::"+fields[3])
                        else:
                            print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                            outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
                    except:
                        print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                        outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
                else:
                    print "error in vcf file found!"
                    sys.exit()
            else:
                pass
    vcf_in.close()
    vcf_out.close()
    return outdata

def sort_information(x):
    fields = x.split("::")
    return int(fields[1])

def merge_vcfs(matrix):
    curr_dir= os.getcwd()
    all_ids= [ ]
    names=[ ]
    in_matrix = open(matrix, "U")
    """these are all of the screened SNPs"""
    matrix_ids=[ ]
    firstLine = in_matrix.readline()
    for line in in_matrix:
        mfields=line.split()
        matrix_ids.append(mfields[0])
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered.vcf")):
        for line in open(infile, "U"):
            fields=line.split()
            all_ids.append(fields[0])
            names.append(infile)
    combined=matrix_ids+all_ids
    nr=[x for i, x in enumerate(combined) if x not in combined[i+1:]]
    nr_sorted=sorted(nr, key=sort_information)
    #open("ref.list", "a").write("\n")
    #for x in nr_sorted:
    #    open("ref.list", "a").write("%s\n" % x)
    outnames=[]
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered.vcf")):
        name=get_seq_name(infile)
        reduced=name.replace(".filtered.vcf","")
        outnames.append(reduced)
        open("%s.tmp.matrix" % name, 'a').write('%s\n' % name)
        value_dict={}
        new_dicts={}
        for line in open(infile, "U"):
            fields=line.split()
            value_dict.update({fields[0]:fields[1]})
        for k,v in value_dict.iteritems():
            if k in nr_sorted: new_dicts.update({k:v})
            else: new_dicts.update({k:"-"})
        for x in nr_sorted:
            if x not in value_dict.keys():new_dicts.update({x:"-"})
        for key in sorted(new_dicts.iterkeys(), key=sort_information):
            open("%s.tmp.matrix" % name, 'a').write("%s\n" % new_dicts[key])
    return outnames
    in_matrix.close()

def merge_matrix(matrix, merged_vcf):
    curr_dir= os.getcwd()
    in_matrix = open(matrix, "U")
    out_matrix = open("combined.matrix", "a")
    in_vcf = open(merged_vcf, "U")
    firstLine = in_matrix.readline()
    first_fields=firstLine.split()
    last=first_fields.index("#SNPcall")
    query_names=[ ]
    vcf_first=in_vcf.readline()
    vcf_first_fields=vcf_first.split()
    print >> out_matrix, "\t".join(first_fields[:last]),"\t","\t".join(vcf_first_fields),"\t","\t".join(first_fields[last:])
    for inline in open(matrix,"U"):
        matrix_fields=inline.split()
        for myline in open(merged_vcf, "U"):
            vcf_fields=myline.split()
            if matrix_fields[0] == vcf_fields[0]:
                print >> out_matrix,"\t".join(matrix_fields[:last]),"\t","\t".join(vcf_fields[1:]),"\t","\t".join(matrix_fields[last:])
            else:
                pass
    in_matrix.close()
    out_matrix.close()
                
def matrix_to_fasta(matrix_in):
    reduced = [ ]
    out_fasta = open("all.fasta", "w")
    for line in open(matrix_in,"U"):
        fields = line.split()
        reduced.append(fields[1:])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])
    out_fasta.close()

def dist_seqs(fasta_in, outnames):
    os.system('mothur "#dist.seqs(fasta=%s, calc=nogaps)" > /dev/null 2>&1' % fasta_in)
    reduced = [ ]
    true_dists = ()
    for name in outnames:
        mydict={}
        outfile = open("%s.distances.txt" % name, "w")
        for line in open("all.dist", "U"):
            if line.startswith('%s' % name):
		fields = line.split()
		str1 = "".join(fields[1])
		str2 = "".join(fields[2])
		mydict.update({str1:str2})
        temp=sorted(mydict.items(), key=itemgetter(1))
        for x in temp:
            print >> outfile,"\t".join(x)
        reduced = temp[:5]
        print "closest genome from %s" % name,"\t","distance"
        for x in reduced:
            print "\t".join(x)
        print ""
        for y in reduced[:2]:
            #print "true distance of %s to %s = %s" % (name,y[0],y[1])
            y_fixed=y[0].replace('__','::')
            y_fixed_2=y_fixed.replace('.filtered.vcf','')
            true_dists=((name,y_fixed_2,y[1]),)+true_dists
        outfile.close()
    return true_dists

def find_two():
    curr_dir= os.getcwd()
    dist_sets = ()
    for infile in glob.glob(os.path.join(curr_dir, "*.distances.txt")):
        with open(infile) as f:
            newlist=[]
            temp = get_seq_name(infile)
            reduced = temp.replace(".distances.txt","")
            for line in f.readlines()[:2]:
                fields=line.split()
                new_fields=[]
                for x in fields:
                    new_fields.append(x.replace('__','::'))
                newlist.append(new_fields[0])
            dist_sets=((reduced,newlist[0],newlist[1]),)+dist_sets
    return dist_sets
                

def run_raxml(fasta_in, tree, processors):
    args = ['raxmlHPC-PTHREADS', '-T', '%s' % processors, '-f', 'v',
	     '-s', '%s' % fasta_in, '-m', 'GTRGAMMA', '-n', 'out', '-t',
	     '%s' % tree, '>', '/dev/null 2>&1']
    try:
        vcf_fh = open('raxml.out', 'w')
    except:
        log_isg.logPrint('could not open raxml file')
    try:
        log_fh = open('raxml.log', 'w')
    except:
        log_isg.logPrint('could not open log file')
        
    log_isg.logPrint("inserting sequence into tree")
    try:
        raxml_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
        raxml_run.wait()
        log_isg.logPrint("sequence(s) inserted into tree")
    except:
        log_isg.logPrint("sequence(s) were not inserted into tree!!!!!")
    os.system("sed 's/\[[^]]*\]//g' RAxML_labelledTree.out > tree_including_unknowns_noedges.tree")
    subprocess.check_call("mv RAxML_labelledTree.out tree_including_unknowns_edges.tree" , shell=True)
    subprocess.check_call("rm RAxML_*", shell=True)

def grab_matrix_coords(matrix):
    coords = [ ]
    my_matrix = open(matrix, "U")
    firstLine = my_matrix.readline()
    for line in my_matrix:
	fields = line.split()
	coords.append(fields[0])
    return coords
    my_matrix.close()

def subsample_snps(matrix, dist_sets, used_snps,subnums):
    allSNPs = [ ]
    for line in open(matrix, "U"):
        if line.startswith("LocusID"):
            pass
        else:
            fields=line.split()
            allSNPs.append(fields[0])
    for k,v in used_snps.iteritems():
        kept_snps=random.sample(set(allSNPs), int(v))
        for z in dist_sets:
            if z[0]==k:
                for x in range(1,subnums):
                    outfile = open("%s.%s.%s.tmp.matrix" % (k,x,z[1]), "w")
                    in_matrix=open(matrix,"U")
                    firstLine = in_matrix.readline()
                    print >> outfile, firstLine,
                    first_fields = firstLine.split()
                    gindex=first_fields.index(z[1])
                    for line in in_matrix:
                        matrix_fields=line.split()
                        if matrix_fields[0] in kept_snps:
                            print >> outfile, line,
                        else:
                            print >> outfile, "\t".join(matrix_fields[:gindex-1])+"\t"+"-"+"\t"+"\t".join(matrix_fields[gindex:])+"\n",
                for x in range(1,subnums):
                    outfile_2 = open("%s.%s.%s.tmp.matrix" % (k,x,z[2]), "w")
                    in_matrix=open(matrix,"U")
                    firstLine = in_matrix.readline()
                    print >> outfile_2, firstLine,
                    first_fields = firstLine.split()
                    gindex=first_fields.index(z[2])
                    for line in in_matrix:
                        matrix_fields=line.split()
                        if matrix_fields[0] in kept_snps:
                            print >> outfile_2, line,
                        else:
                            print >> outfile_2, "\t".join(matrix_fields[:gindex-1])+"\t"+"-"+"\t"+"\t".join(matrix_fields[gindex:])+"\n",
def find_used_snps():
    curr_dir= os.getcwd()
    used_SNPs = {}
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered.vcf")):
        name=get_seq_name(infile)
        reduced=name.replace('.filtered.vcf', '')
        good_snps=[]
        for line in open(infile, "U"):
            fields=line.split()
            if fields[1] != "-":
                good_snps.append("1")
        used_SNPs[str(reduced)] = int(len(good_snps))
    return used_SNPs

def process_temp_matrices():
   curr_dir= os.getcwd()
   for infile in glob.glob(os.path.join(curr_dir, "*tmp.matrix")):
       name=get_seq_name(infile)
       split_fields=name.split(".")
       outfile=open("%s.%s.subsample.distances.txt" % (split_fields[0],split_fields[2]), "a")
       try:
           matrix_to_fasta(infile)
       except:
           print "problem encountered with file ",infile
       os.system('mothur "#dist.seqs(fasta=all.fasta, calc=nogaps)" > /dev/null 2>&1')
       os.system('sed "s/.filtered.vcf//g" all.dist > renamed.dist')
       for line in open("renamed.dist", "U"):
           if line.startswith("%s" % split_fields[0]):
               fields=line.split()
               new_fields=[]
               for x in fields:
                   new_fields.append(x.replace('__','::'))
               if new_fields[0]==split_fields[0] and new_fields[1]==split_fields[2]:
                   print >> outfile, new_fields[2]+"\n",
               else:
                   pass

def compare_subsample_results(true_dists):
    curr_dir= os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, "*.subsample.distances.txt")):
        all_dists=[ ]
        dists_greater_than_true=[ ]
        dists_equal_to_true=[ ]
        dists_less_than_true=[ ]
        for line in open(infile, "U"):
            all_dists.append(line)
        name=get_seq_name(infile)
        split_fields=name.split(".")
        max_dist=max(all_dists)
        print "maximum subsample distance between %s and %s = %s" % (split_fields[0],split_fields[1],max_dist),
        for x in true_dists:
            if x[0] == split_fields[0] and x[1] == split_fields[1]:
                for dists in all_dists:
                    if float(x[2])>float(dists):
                        dists_greater_than_true.append("1")
                    elif float(x[2])==float(dists):
                        dists_equal_to_true.append("1")
                    else:
                        dists_less_than_true.append("1")
            else:
                pass
        greaters = int(len(dists_greater_than_true))
        equals = int(len(dists_equal_to_true))
        lessers = int(len(dists_less_than_true))
        for x in true_dists:
            if x[0] == split_fields[0] and x[1] == split_fields[1]:
                print "True distance between %s and %s = %s" % (split_fields[0],split_fields[1],x[2])
        print "Subsample distances between %s and %s greater than true value = %s" % (split_fields[0],split_fields[1],greaters)
        print "Subsample distances between %s and %s equal to true value = %s" % (split_fields[0],split_fields[1],equals)
        print "Subsample distances between %s and %s less than true value = %s" % (split_fields[0],split_fields[1],lessers)    
        print ""

def transform_tree(tree):
    infile = open(tree, "U")
    tree_string = []
    for line in infile:
        tree_string.append(line)
    infile.close()
    mytree = Phylo.read(tree, 'newick')
    tree_names = [ ]
    for clade in mytree.find_clades():
        if clade.name:
            tree_names.append(clade.name)
    outfile = open("transformed.tree", "w")
    print >> outfile, "#NEXUS"+"\n",
    print >> outfile, "begin taxa;"+"\n",
    print >> outfile, "\t"+"dimensions ntax="+str(len(tree_names))+";"+"\n",
    print >> outfile, "\t"+"taxlabels"+"\n",
    for clade in mytree.find_clades():
        if clade.name:
            tree_names.append(clade.name)
    nr=[x for i, x in enumerate(tree_names) if x not in tree_names[i+1:]]
    for tree_name in nr:
        if "QUERY__" in tree_name:
            print >> outfile, "\t"+"'%s'" % tree_name+"[&!color=#-3407821]"+"\n",
        else:
            print >> outfile, "\t"+tree_name+"\n",
    print >> outfile, ";"+"\n",
    print >> outfile, "end;"+"\n",
    print >> outfile, ""+"\n",
    print >> outfile, "begin trees;"+"\n",
    for x in tree_string:
        print >> outfile, "\t"+"tree tree_1 = [&R] ",str(x)+"\n",
    print >> outfile, "end;"+"\n",

def write_reduced_matrix(matrix):
    in_matrix = open(matrix, "U")
    outfile = open("temp.matrix", "w")
    firstLine = in_matrix.readline()
    first_fields=firstLine.split()
    last=first_fields.index("#SNPcall")
    print >> outfile, "\t".join(first_fields[:last])+"\n",
    for line in in_matrix:
        fields = line.split()
        print >> outfile, "\t".join(fields[:last])
    outfile.close()
    in_matrix.close()
        
def make_temp_matrix(vcf, matrix, name):
    in_matrix = open(matrix, "U")
    """these are all of the screened SNPs"""
    matrix_ids=[ ]
    firstLine = in_matrix.readline()
    for line in in_matrix:
        mfields=line.split()
        matrix_ids.append(mfields[0])
    open("%s.tmp.matrix" % name, 'a').write('%s\n' % name)
    value_dict={}
    new_dicts={}
    nr_sorted=sorted(matrix_ids, key=sort_information)
    for line in open(vcf, "U"):
        fields=line.split()
        value_dict.update({fields[0]:fields[1]})
    for k,v in value_dict.iteritems():
        if k in nr_sorted: new_dicts.update({k:v})
        else: new_dicts.update({k:"-"})
    for x in nr_sorted:
        if x not in value_dict.keys():new_dicts.update({x:"-"})
    for key in sorted(new_dicts.iterkeys(), key=sort_information):
        open("%s.tmp.matrix" % name, 'a').write("%s\n" % new_dicts[key])

def grab_names():
    curr_dir= os.getcwd()
    outnames = [ ]
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered.vcf")):
        name=get_seq_name(infile)
        reduced=name.replace(".filtered.vcf","")
        outnames.append(reduced)
    return outnames


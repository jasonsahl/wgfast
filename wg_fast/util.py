from __future__ import division
import os
import re
import logging
import glob
import subprocess
import sys
from subprocess import Popen
try:
    from Bio import SeqIO
    from Bio import Phylo
except:
    print "BioPython is not in your PYTHONPATH, but needs to be"
    sys.exit()
try:
    import dendropy
    from dendropy import treecalc
except:
    print "dendropy is not installed, but needs to be"
    sys.exit()
import glob
from igs.threading import functional as p_func
from igs.utils import logging as log_isg
from operator import itemgetter
import random
import threading

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

def process_coverage(name):
    """function required in loop"""
    curr_dir= os.getcwd()
    outfile = open("coverage_out.txt", "a")
    coverage_dict = {}
    try:
        infile = open("%s_coverage.sample_summary" % name, "U")
    except:
        print "infile cannot be used"
        sys.exit()
    for line in infile:
        fields = line.split()
        if fields[0] == name:
            coverage_dict.update({fields[0]:fields[2]})
    if len(coverage_dict)>=1:
        for k,v in coverage_dict.iteritems():
            print >> outfile,k,v+"\n",
    else:
        raise TypeError("dictionary is empty")
    return coverage_dict

def run_loop(fileSets, dir_path, reference, processors, gatk, ref_coords, coverage, proportion, matrix,ap,doc,tmp_dir):
    files_and_temp_names = [(str(idx), list(f))
                            for idx, f in fileSets.iteritems()]
    lock = threading.Lock()
    def _perform_workflow(data):
        """idx is the sample name, f is the file dictionary"""
        idx, f = data
        if len(f)>1:
            run_bwa(reference, f[0], f[1], processors, idx)
        else:
            run_bwa(reference, f[0], "NULL", processors, idx)
        process_sam("%s.sam" % idx, idx)
        run_gatk(reference, processors, idx, gatk, tmp_dir)
        os.system("echo %s.bam > %s.bam.list" % (idx,idx))
        if "T" == doc:
            lock.acquire()
            os.system("java -Djava.io.tmpdir=%s -jar %s -R %s/scratch/reference.fasta -T DepthOfCoverage -o %s_coverage -I %s.bam.list -rf BadCigar > /dev/null 2>&1" % (tmp_dir,gatk,ap,idx,idx))
            lock.release()
            process_coverage(idx)
        else:
            pass
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
    if "NULL" in read2:
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
    """launches bwa. Adds in read_group for compatability with GATK"""
    read_group = '@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tPU:%s' % (name,name,name)
    bwa(reference,read1, read2,"%s.sam" % name, processors, log_file='%s.sam.log' % name,**{'-R':read_group}) 

def process_sam(in_sam, name):
    """samtools runs to remove multiple mapped reads and unmapped reads"""
    subprocess.check_call("samtools view -h -b -S %s > %s.1.bam 2> /dev/null" % (in_sam, name), shell=True)
    subprocess.check_call("samtools view -u -h -F4 -o %s.2.bam %s.1.bam > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools view -h -b -q1 -F4 -o %s.3.bam %s.2.bam > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools sort %s.3.bam %s > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools index %s.bam > /dev/null 2>&1" % name, shell=True)
    subprocess.check_call("rm %s.1.bam %s.2.bam %s.3.bam %s" % (name, name, name, in_sam), shell=True)

def run_gatk(reference, processors, name, gatk, tmp_dir):
    """gatk controller"""
    args = ['java', '-Djava.io.tmpdir=%s' % tmp_dir, '-jar', '%s' % gatk, '-T', 'UnifiedGenotyper',
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

def process_vcf(vcf, ref_coords, coverage, proportion, name):
    """finds SNPs that pass user-defined thresholds
    for coverage and proportion"""
    vcf_in = open(vcf, "U")
    vcf_out = open("%s.filtered.vcf" % name, "w")
    outdata = []
    good_snps = [ ]
    mixed_snps = [ ]
    ref_set = set(ref_coords)
    for line in vcf_in:
        if line.startswith('#'):
           pass
        elif line.startswith('Java'):
            pass
        elif line.startswith('/tmp'):
            pass
        elif line.startswith('Try'):
            pass
        elif line.startswith('INFO'):
            pass
        elif line in ['\n', '\r\n']:
            pass
        else:
            fields=line.split()
            """for GATK, a period signifies a reference call.
            First we want to look at the situation where this is
            not the case"""
            try:
                merged_fields=fields[0]+"::"+fields[1]
                if merged_fields in ref_set:
                    if "." != fields[4]:
                        snp_fields=fields[9].split(':')
                        if int(len(snp_fields))>2:
                            prop_fields=snp_fields[1].split(',')
                            if int(snp_fields[2])>=coverage:
                                if int(prop_fields[1])/int(snp_fields[2])>=float(proportion):
                                    print >> vcf_out, fields[0]+"::"+fields[1],fields[4]+"\n",
                                    outdata.append(fields[0]+"::"+fields[1]+"::"+fields[4])
                                    good_snps.append("1")
                                else:
                                    print >> vcf_out, fields[0]+"::"+fields[1],"-"+"\n",
                                    outdata.append(fields[0]+"::"+fields[1]+"::"+"-")
                                    mixed_snps.append("1")
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
            except:
                pass
    print "number of SNPs in genome %s = " % name, len(good_snps)
    print "number of discarded SNPs in genome %s = " % name, len(mixed_snps)
    vcf_in.close()
    vcf_out.close()
    return outdata
    
def sort_information(x):
    try:
        fields = x.split("::")
        return int(fields[1])
    except:
        raise TypeError("problem encountered parsing fields")

def matrix_to_fasta(matrix_in):
    """function to convert a SNP matrix to a multi-fasta file"""
    reduced = [ ]
    out_fasta = open("all.fasta", "w")
    redux = [ ]
    for line in open(matrix_in,"U"):
        fields = line.split()
        reduced.append(fields[1:])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])
        redux.append(">"+str(x[0])+"".join(x[1:3]))
    out_fasta.close()
    return redux

def run_raxml(fasta_in, tree, processors, out_class_file, insertion_method, parameters):
    if "NULL" == parameters:
        args = ['raxmlHPC-PTHREADS', '-T', '%s' % processors, '-f', '%s' % insertion_method,
	     '-s', '%s' % fasta_in, '-m', 'GTRGAMMA', '-n', 'out', '-t',
	     '%s' % tree, '>', '/dev/null 2>&1']
    else:
        args = ['raxmlHPC-PTHREADS', '-T', '%s' % processors, '-f', '%s' % insertion_method,
	     '-s', '%s' % fasta_in, '-m', 'GTRGAMMA', '-n', 'out', '-R', parameters, '-t',
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
    try:
        subprocess.check_call("cat RAxML_classificationLikelihoodWeights.out >> %s" % out_class_file, shell=True)
    except:
        pass
    subprocess.check_call("rm RAxML_*", shell=True)

def grab_matrix_coords(matrix):
    """retrieve all of the coordinates from a
    NASP formatted SNP matrix"""
    coords = [ ]
    my_matrix = open(matrix, "U")
    firstLine = my_matrix.readline()
    for line in my_matrix:
	fields = line.split()
	coords.append(fields[0])
    return coords
    my_matrix.close()

def subsample_snps(matrix, dist_sets, used_snps, subnums):
    """get a list of all possible positions, depending
    on those positions in the original matrix"""
    allSNPs = [ ]
    for line in open(matrix, "U"):
        if line.startswith("LocusID"):
            pass
        else:
            fields=line.split()
            allSNPs.append(fields[0])
    for k,v in used_snps.iteritems():
        for z in dist_sets:
            if len(z) == 0:
                pass
            else:
                if z[0]==k:
                    for x in range(1,int(subnums)+1):
                        kept_snps=random.sample(set(allSNPs), int(v))
                        outfile = open("%s.%s.%s.tmp.matrix" % (k,x,z[1]), "w")
                        in_matrix=open(matrix,"U")
                        firstLine = in_matrix.readline()
                        print >> outfile, firstLine,
                        first_fields = firstLine.split()
                        fixed_fields = []
                        for x in first_fields:
                            fixed_fields.append(re.sub('[:,]', '', x))
                        gindex=fixed_fields.index(z[1])
                        for line in in_matrix:
                            matrix_fields=line.split()
                            if matrix_fields[0] in kept_snps:
                                print >> outfile, line,
                            else:
                                print >> outfile, "\t".join(matrix_fields[:gindex])+"\t"+"-"+"\t"+"\t".join(matrix_fields[gindex+1:])+"\n",
                        outfile.close()
    return allSNPs

def find_used_snps():
    """report how many SNPs were used in a given sample.  This is
    then used for the sub-sampling routine"""
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

def branch_lengths_2_decimals(str_newick_tree):
    """replaces branch lengths in scientific notation with decimals"""
    colon_s = 0
    comma_back_paren_s = 0
    num = ''
    new_tree = ''
    for count, char in enumerate(str_newick_tree):
        if char == ':':
            colon_s = count
            continue
        if char in (')', ','):
            comma_back_paren_s = 1
            num = '%f' % float(num)
            new_tree += ":" + num
            colon_s = 0
            num = ''
        if colon_s != 0:
            num = num + char
        if colon_s == 0:
            new_tree += char
    new_tree = new_tree.strip('\'').strip('\"').strip('\'') + ";"
    return new_tree

def process_temp_matrices(dist_sets, tree, processors, patristics, insertion_method, parameters):
    curr_dir= os.getcwd()
    os.system("rm tree_including_unknowns_noedges.tree")
    for infile in glob.glob(os.path.join(curr_dir, "*tmp.matrix")):
        """the genome names are parsed out of the tmp.matrices"""
        name=get_seq_name(infile)
        split_fields=name.split(".")
        outfile=open("%s.%s.subsample.distances.txt" % (split_fields[0],split_fields[2]), "a")
        name_fixed = []
        name_fixed.append(re.sub('[:,]', '', split_fields[2]))
        tmptree = open("tmpx.tree", "w")
        """The genomes in the dist_sets will be pruned in
        the subsample routine"""
        to_prune = []
        for x in dist_sets:
            if x[0] == split_fields[0]: 
                to_prune.append(x[1])
        """names will be fixed if they contain characters
        that are not accepted by downstream applications"""
        to_prune_fixed=[]
        for x in to_prune:
            to_prune_fixed.append(re.sub('[:,]', '', x))
        """dendropy is used here to import the tree and prune the taxa"""
        tree_full = dendropy.Tree.get_from_path(tree,schema="newick",preserve_underscores=True)
        tree_full.prune_taxa_with_labels(to_prune_fixed)
        """dendropy uses scientific notation, which needs to be converted
        into decimals"""
        final_tree = branch_lengths_2_decimals(tree_full.as_string("newick"))
        print >> tmptree, final_tree
        tmptree.close()
        """A new tree file is created, changing the dendropy
        format into a more classical Newick format"""
        tmptree2 = open("tmpxz.tree", "w")
        for line in open("tmpx.tree", "U"):
            if line.startswith("[&U]"):
                fields = line.split()
                fixed_fields = [ ]
                for x in fields:
                    fixed_fields.append(x.replace("'",""))
                print >> tmptree2, fixed_fields[1]
            else:
                pass
        tmptree2.close()
        try:
            matrix_to_fasta(infile)
        except:
            print "problem converting matrix to fasta"
        """if problems in the tree names are found, they are removed by the system command"""
        os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
        """raxml is now used to insert the pruned genomes back into the tree"""
        run_raxml("out.fasta", "tmpxz.tree", processors, "subsampling_classifications.txt", insertion_method, parameters)
        """dendropy is used to calculate pairwise patristic distances"""
        calculate_pairwise_tree_dists("tree_including_unknowns_noedges.tree", "resampling_distances.txt")
        """parse the results from raxml and save the results to the subsamples file"""
        for line in open("resampling_distances.txt","U"):
            resample_fields = line.split()
            myid = re.sub("[:']", "",resample_fields[4])
            fixedid = myid.replace("QUERY___","")
            newid = re.sub("[:']","",resample_fields[2])
            fixedid2 = newid.replace("QUERY___","")
            if resample_fields[2] == "'Reference'" and fixedid in name_fixed:
                print >> outfile, "resampled distance between Reference and %s = %s" % (fixedid, resample_fields[5])
            elif resample_fields[4] == "'Reference':" and fixedid2 in name_fixed:
                print >> outfile, "resampled distance between Reference and %s = %s" % (fixedid2, resample_fields[5])
            else:
                pass
            
        os.system("rm all.fasta tmpxz.tree out.fasta tmpx.tree tree_including_unknowns_noedges.tree resampling_distances.txt")
        
def compare_subsample_results(outnames,distances,fudge):
    curr_dir= os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, "*.subsample.distances.txt")):
        name=get_seq_name(infile)
        split_fields=name.split(".")
        genomes_used = [ ]
        all_dists=[ ]
        dists_greater_than_true=[ ]
        dists_equal_to_true=[ ]
        dists_less_than_true=[ ]
        """Here, I get the list of the genomes that are being analyzed"""
        for line in open(infile, "U"):
            fields = line.split()
            all_dists.append(float(fields[7]))
        try:
            max_dist=max(all_dists)
            print ""
            print "maximum subsample distance between %s and %s = %.2f" % (fields[3],fields[5],float(max_dist)),"\n",
        except:
            print "problem found in input file: ", infile
        true_dists = [ ]
        for distance in distances:
            if distance[1] == split_fields[1]:
                true_dists.append(distance[2])
        for all_dist in all_dists:
            if "%.2f" % float(all_dist)> "%.2f" % (float(true_dists[0])+float(fudge)):
                dists_greater_than_true.append("1")
            elif "%.2f" % float(all_dist)<"%.2f" % (float(true_dists[0])-float(fudge)):
                dists_less_than_true.append("1")
            else:
                dists_equal_to_true.append("1")
        greaters = int(len(dists_greater_than_true))
        equals = int(len(dists_equal_to_true))
        lessers = int(len(dists_less_than_true))
        print "True distance between Reference and %s = %.2f" % (split_fields[1],float(true_dists[0]))
        print "Sample: %s" % split_fields[0]
        print "Subsample distances between Reference and %s greater than true value = %s" % (split_fields[2],greaters)
        print "Subsample distances between Reference and %s equal to true value = %s" % (split_fields[2],equals)
        print "Subsample distances between Reference and %s less than true value = %s" % (split_fields[2],lessers)    
        p = (greaters+lessers)/(greaters+lessers+equals)
        print "Placement p value = %.3f" % float(p)
        
def transform_tree(tree):
    """converts a Newick tree into a Nexus-formatted
    tree that can be visualized with FigTree"""
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
    infile.close()
    outfile.close()

def write_reduced_matrix(matrix):
    """This function takes a NASP formatted
    SNP matrix and writes a temporary matrix
    that can be easily combined with temporary files"""
    in_matrix = open(matrix, "U")
    outfile = open("temp.matrix", "w")
    outdata = [ ]
    firstLine = in_matrix.readline()
    first_fields=firstLine.split()
    last=first_fields.index("#SNPcall")
    print >> outfile, "\t".join(first_fields[:last])+"\n",
    for line in in_matrix:
        fields = line.split()
        print >> outfile, "\t".join(fields[:last])
        outdata.append(len(fields[:last]))
    outfile.close()
    in_matrix.close()
    return outdata
        
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
    for line in open(vcf, "U"):
        fields=line.split()
        value_dict.update({fields[0]:fields[1]})
    for k,v in value_dict.iteritems():
        new_dicts.update({k:v})
    for x in matrix_ids:
        if x not in value_dict.keys():new_dicts.update({x:"-"})
    variety = [ ]
    for x in matrix_ids:
        if x in new_dicts:
            variety.append(new_dicts.get('%s' % x))
    if "A" or "T" or "G" or "C" in variety:
        for x in matrix_ids:
            if x in new_dicts:
                open("%s.tmp.matrix" % name, 'a').write("%s\n" % new_dicts.get('%s' % x))
    else:
        print "sample %s had no usable positions!!!" % name
    in_matrix.close()
    return new_dicts

def grab_names():
    """These names will be used for future iterations"""
    curr_dir= os.getcwd()
    outnames = [ ]
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered.vcf")):
        name=get_seq_name(infile)
        reduced=name.replace(".filtered.vcf","")
        outnames.append(reduced)
    return outnames

def parse_likelihoods(infile):
    """This function parses the likelihoods output from RAxML"""
    my_in = open(infile, "U")
    like_dict = {}
    for line in my_in:
        fields = line.split()
        try:
            like_dict[fields[0]].append(fields[2])
        except KeyError:
            like_dict[fields[0]] = [fields[2]]
    print "sample_name","\t","insertion_likelihood","\t","number of potential insertion nodes"
    for k,v in like_dict.iteritems():
        print k+"\t"+v[0]+"\t"+str(len(v))
    my_in.close()

def calculate_pairwise_tree_dists(intree, output):
    """uses dendropy function to calculate all pairwise distances between tree"""
    tree = dendropy.Tree.get_from_path(intree, "newick", preserve_underscores=True)
    outfile = open("%s" % output, "w")
    distances = treecalc.PatristicDistanceMatrix(tree)
    distances_sets = [ ]
    for i, t1 in enumerate(tree.taxon_set):
        for t2 in tree.taxon_set[i+1:]:
            distances_sets.append(distances(t1, t2))
    try:
        for i, t1 in enumerate(tree.taxon_set):
            for t2 in tree.taxon_set[i+1:]:
                print >> outfile, "Distance between '%s' and '%s': %s" % (t1.label, t2.label, distances(t1, t2))
    except:
        print "problem iterating through tree.  Tree is empty or not Newick format"
    outfile.close()
    return distances_sets

def get_closest_dists_new(final_sets, outnames):
    results = []
    for final_set in final_sets:
        if len(final_set) == 0:
            pass
        results.append(final_set[1]+final_set[2])
    return results
            
def find_two_new(infile,outnames):
    """find two closest genomes to each query genome,
    return the names and distances (sorted), for the
    two genomes"""
    distances = ()
    output_tuples = ()
    for outname in outnames:
        outname_tuple = ()
        for line in open(infile, "U"):
            fields=line.split()
            new_fields=[ ]
            for x in fields:
                new_fields.append(re.sub('[:,]', '', x))
            final_fields=[ ]
            for y in new_fields:
                final_fields.append(y.replace("QUERY___",""))
            """We need the distance of all samples, compared to the reference"""
            if "Reference" in final_fields[2].replace("'",""):
                distances=((final_fields[2].replace("'",""),final_fields[4].replace("'",""),final_fields[5].replace("'","")),)+distances
            elif "Reference" in final_fields[4].replace("'",""):
                distances=((final_fields[4].replace("'",""),final_fields[2].replace("'",""),final_fields[5].replace("'","")),)+distances
            if final_fields[4].replace("'","") == outname:
                if "Reference" not in final_fields[2].replace("'","") and final_fields[2].replace("'","") not in outnames:
                    outname_tuple=((final_fields[4].replace("'",""),final_fields[2].replace("'",""),final_fields[5].replace("'","")),)+outname_tuple
            elif final_fields[2].replace("'","") == outname:
                if "Reference" not in final_fields[4].replace("'","") and final_fields[4].replace("'","") not in outnames:
                    outname_tuple=((final_fields[2].replace("'",""),final_fields[4].replace("'",""),final_fields[5].replace("'","")),)+outname_tuple
        for result in sorted(outname_tuple,key=lambda x: float(x[2]))[:2]:
            output_tuples=((result),)+output_tuples
    return output_tuples,distances

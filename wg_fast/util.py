from __future__ import division
import os,os.path
import re
import logging
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
try:
    from igs.threading import functional as p_func
    from igs.utils import logging as log_isg
except:
    print "your PYTHONPATH needs to include the wg-fast repository"
    sys.exit()
from operator import itemgetter
import threading
import collections
import random

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

def test_methods(option, opt_str, value, parser):
    if "ML" in value:
        setattr(parser.values, option.dest, value)
    elif "MP" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported.  Only select from MP or ML"
        sys.exit()

def test_models(option, opt_str, value, parser):
    if "GTRGAMMA" in value:
        setattr(parser.values, option.dest, value)
    elif "ASC_GTRGAMMA" in value:
        setattr(parser.values, option.dest, value)
    elif "GTRCAT" in value:
        setattr(parser.values, option.dest, value)
    elif "ASC_GTRCAT" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "substitution model is not supported"
        sys.exit()

def get_seq_name(in_fasta):
    """used for renaming the sequences - tested"""
    return os.path.basename(in_fasta)

def get_readFile_components(full_file_path):
    """function adapted from:
    https://github.com/katholt/srst2 - tested"""
    (file_path,file_name) = os.path.split(full_file_path)
    m1 = re.match("(.*).gz",file_name)
    ext = ""
    if m1 is not None:
        ext = ".gz"
        file_name = m1.groups()[0]
    (file_name_before_ext,ext2) = os.path.splitext(file_name)
    full_ext = ext2+ext
    return file_path,file_name_before_ext,full_ext

def read_file_sets(dir_path):        
    """match up pairs of sequence data, adapted from
    https://github.com/katholt/srst2 will be tough to test
    with variable names and read paths"""
    fileSets = {} 
    forward_reads = {}
    reverse_reads = {} 
    num_paired_readsets = 0
    num_single_readsets = 0
    for infile in glob.glob(os.path.join(dir_path, "*.fastq.gz")):
        (file_path,file_name_before_ext,full_ext) = get_readFile_components(infile)
        m=re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
        if m is None:
            m=re.match("(.*)("+"_R1"+")(_.*)$",file_name_before_ext)
            if m is not None:
                (baseName,read) = m.groups()[0], m.groups()[1]
                forward_reads[baseName] = infile
            else:
                m=re.match("(.*)("+"_R2"+")(_.*)$",file_name_before_ext)
                if m is not None:
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
    """function required in loop - tested"""
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
    infile.close()
    outfile.close()

def run_loop(fileSets, dir_path, reference, processors, gatk, ref_coords, coverage, proportion, matrix,ap,doc,tmp_dir,picard,trim_path,wgfast_path):
    files_and_temp_names = [(str(idx), list(f)) for idx, f in fileSets.iteritems()]
    lock = threading.Lock()
    def _perform_workflow(data):
        """idx is the sample name, f is the file dictionary"""
        idx, f = data
        if len(f)>1:
            """paired end sequences - won't work for old, short sequences"""
            args=['java','-jar','%s' % trim_path,'PE', '-threads', '%s' % processors,
                  '%s' % f[0], '%s' % f[1], '%s.F.paired.fastq.gz' % idx, 'F.unpaired.fastq.gz',
	          '%s.R.paired.fastq.gz' % idx, 'R.unpaired.fastq.gz', 'ILLUMINACLIP:%s/bin/illumina_adapters_all.fasta:2:30:10' % wgfast_path,
	          'MINLEN:50']
            try:
                vcf_fh = open('%s.trimmomatic.out' % idx, 'w')
            except:
                log_isg.logPrint('could not open trimmomatic file')
            try:
                log_fh = open('%s.trimmomatic.log' % idx, 'w')
            except:
                log_isg.logPrint('could not open log file')
            if os.path.isfile("%s.F.paired.fastq.gz" % idx):
                pass
            else:
                try:
                    trim = Popen(args, stderr=vcf_fh, stdout=log_fh)
                    trim.wait()
                except:
                    log_isg.logPrint('problem enountered trying to run trimmomatic')
            if os.path.isfile("%s_renamed_header.bam" % idx):
                pass
            else:
                run_bwa(reference, '%s.F.paired.fastq.gz' % idx, '%s.R.paired.fastq.gz' % idx, processors, idx)
        else:
            """single end support"""
            args=['java','-jar','%s' % trim_path,'SE', '-threads', '%s' % processors,
                  '%s' % f[0], '%s.single.fastq.gz' % idx, 'ILLUMINACLIP:%s/bin/illumina_adapters_all.fasta:2:30:10' % wgfast_path,
	          'MINLEN:50']
            try:
                vcf_fh = open('%s.trimmomatic.out' % idx, 'w')
            except:
                log_isg.logPrint('could not open trimmomatic file')
            try:
                log_fh = open('%s.trimmomatic.log' % idx, 'w')
            except:
                log_isg.logPrint('could not open log file')
            if os.path.isfile("%s.single.fastq.gz" % idx):
                pass
            else:
                try:
                    trim = Popen(args, stderr=vcf_fh, stdout=log_fh)
                    trim.wait()
                except:
                    log_isg.logPrint("problem encountered with trimmomatic")
            if os.path.isfile("%s_renamed_header.bam" % idx):
                pass
            else:
                run_bwa(reference, '%s.single.fastq.gz' % idx, "NULL", processors, idx)
        if os.path.isfile("%s_renamed_header.bam" % idx):
            pass
        else:
            process_sam("%s.sam" % idx, idx)
            """inserts read group information, required by new versions of GATK"""
            os.system("java -jar %s INPUT=%s.bam OUTPUT=%s_renamed_header.bam SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=illumina RGSM=%s RGPU=name CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT > /dev/null 2>&1" % (picard,idx,idx,idx,idx,idx))
            os.system("samtools index %s_renamed_header.bam" % idx)
        run_gatk(reference, processors, idx, gatk, tmp_dir)
        if "T" == doc:
            lock.acquire()
            os.system("echo %s_renamed_header.bam > %s.bam.list" % (idx,idx))
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
    """controller for bwa, currently only works with bwa mem - untested because they are mostly system calls"""
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
    try:
        sam_fh = open(sam_file, 'w')
    except:
        print sam_file, 'could not open'
 
    bwa = Popen(mem_arguments, stderr=log_fh, stdout=sam_fh)
    bwa.wait()

def run_bwa(reference, read1, read2, processors, name):
    """launches bwa. Adds in read_group for compatability with GATK - untested"""
    read_group = '@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tPU:%s' % (name,name,name)
    bwa(reference,read1, read2,"%s.sam" % name, processors, log_file='%s.sam.log' % name,**{'-R':read_group}) 

def process_sam(in_sam, name):
    """samtools runs to remove multiple mapped reads and unmapped reads - untested, system calls"""
    subprocess.check_call("samtools view -h -b -S %s > %s.1.bam 2> /dev/null" % (in_sam, name), shell=True)
    subprocess.check_call("samtools view -u -h -F4 -o %s.2.bam %s.1.bam > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools view -h -b -q1 -F4 -o %s.3.bam %s.2.bam > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools sort %s.3.bam %s > /dev/null 2>&1" % (name, name), shell=True)
    subprocess.check_call("samtools index %s.bam > /dev/null 2>&1" % name, shell=True)
    subprocess.check_call("rm %s.1.bam %s.2.bam %s.3.bam %s" % (name, name, name, in_sam), shell=True)

def run_gatk(reference, processors, name, gatk, tmp_dir):
    """gatk controller, mbq used to be set to 17, but was recently changed - untested, system call"""
    args = ['java', '-Djava.io.tmpdir=%s' % tmp_dir, '-jar', '%s' % gatk, '-T', 'UnifiedGenotyper',
            '-R', '%s' % reference, '-nt', '%s' % processors, '-S', 'silent',
            '-ploidy', '1', '-out_mode', 'EMIT_ALL_CONFIDENT_SITES',
            '-stand_call_conf', '30', '-stand_emit_conf', '30', '-I', '%s_renamed_header.bam' % name,
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
    for coverage and proportion - tested"""
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
    """simple sort - tested"""
    try:
        fields = x.split("::")
        return int(fields[1])
    except:
        raise TypeError("problem encountered parsing fields")

def matrix_to_fasta(matrix_in, outfile):
    """function to convert a SNP matrix to a multi-fasta file - tested"""
    reduced = [ ]
    out_fasta = open(outfile, "w")
    redux = [ ]
    for line in open(matrix_in,"U"):
        newline = line.strip()
        fields = newline.split()
        reduced.append(fields[1:])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])
        redux.append(">"+str(x[0])+"".join(x[1:3]))
    out_fasta.close()
    return redux

def run_raxml(fasta_in, tree, out_class_file, insertion_method, parameters, model, suffix):
    """untested function, system calls"""
    if "NULL" == parameters:
        args = ['raxmlHPC-SSE3', '-f', '%s' % insertion_method,
	     '-s', '%s' % fasta_in, '-m', '%s' % model, '-n', '%s' % suffix, '-t',
	     '%s' % tree, '--no-bfgs', '>', '/dev/null 2>&1']
    else:
        args = ['raxmlHPC-SSE3', '-f', '%s' % insertion_method,
	     '-s', '%s' % fasta_in, '-m', '%s' % model, '-n', '%s' % suffix, '-R', parameters, '-t',
	     '%s' % tree, '--no-bfgs', '>', '/dev/null 2>&1']
    try:
        vcf_fh = open('%s.raxml.out' % suffix, 'w')
    except:
        log_isg.logPrint('could not open raxml file')
    try:
        log_fh = open('%s.raxml.log' % suffix, 'w')
    except:
        log_isg.logPrint('could not open log file')
        
    log_isg.logPrint("inserting sequence into tree")
    try:
        raxml_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
        raxml_run.wait()
        log_isg.logPrint("sequence(s) inserted into tree")
    except:
        log_isg.logPrint("sequence(s) were not inserted into tree!!!!!")
    os.system("sed 's/\[[^]]*\]//g' RAxML_labelledTree.%s > %s.tree_including_unknowns_noedges.tree" % (suffix, suffix))
    subprocess.check_call("mv RAxML_labelledTree.%s %s_tree_including_unknowns_edges.tree" % (suffix, suffix) , shell=True)
    try:
        subprocess.check_call("cat RAxML_classificationLikelihoodWeights.%s >> %s" % (suffix, out_class_file), shell=True)
    except:
        pass
    os.system("rm RAxML_*.%s" % suffix)
    return suffix

def subsample_snps(matrix, dist_sets, used_snps, subnums):
    """get a list of all possible positions, depending
    on those positions in the original matrix - tested"""
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
                        for y in first_fields:
                            fixed_fields.append(re.sub('[:,]', '', y))
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
    then used for the sub-sampling routine - tested"""
    curr_dir= os.getcwd()
    used_SNPs = {}
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered.vcf")):
        name=get_seq_name(infile)
        reduced=name.replace('.filtered.vcf', '')
        good_snps=[]
        for line in open(infile, "U"):
            fields=line.split()
            try:
                if fields[1] != "-":
                    good_snps.append("1")
                else:
                    pass
            except:
                raise TypeError("abnormal number of fields observed")
        used_SNPs[str(reduced)] = int(len(good_snps))
    return used_SNPs

def branch_lengths_2_decimals(str_newick_tree):
    """replaces branch lengths in scientific notation with decimals = tested"""
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

def fasta_to_tab(fasta):
    """tested"""
    infile = open(fasta, "rU")
    outfile = open("out.tab", "w")
    for record in SeqIO.parse(infile, "fasta"):
        """this list is just for testing,
        and is ok if it's overwritten for each
        fasta"""
        for_test = []
        print >> outfile, record.id, record.seq
        for_test.append(record.id)
        for_test.append(str(record.seq))
    infile.close()
    outfile.close()
    return for_test

def tab_to_fasta(new_tab):
    """tested"""
    infile = open(new_tab, "rU")
    outfile = open("out.fasta", "w")
    for line in infile:
        to_test = []
        fields = line.split()
        print >> outfile, ">"+fields[0]
        print >> outfile, fields[1].upper()
        to_test.append(fields[0])
        to_test.append(fields[1].upper())
    infile.close()
    outfile.close()
    return to_test
    
def tab_to_matrix(tab):
    """tested"""
    reduced = [ ]
    out_matrix = open("tab_matrix", "w")
    for line in open(tab):
        tmp_list = []
        fields = line.split()
        tmp_list.append(fields[0])
        for nucs in fields[1]:
            tmp_list.append(nucs.upper())
        reduced.append(tmp_list)
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_matrix, "\t".join(x)
    out_matrix.close()
    return test
        
def filter_alignment(tab):
    """currently untested, but needs to be"""
    outfile = open("tab.filtered", "w")
    infile = open(tab, "U")
    firstLine = infile.readline()
    print >> outfile, firstLine,
    for line in infile:
        valid_fields = []
        fields = line.split()
        for field in fields:
            """skip fields that might be present when missing data are included"""
            if field != "-" and field != "N" and field != "X":
                valid_fields.append(field)
            else:
                pass
        counter=collections.Counter(valid_fields)
        values=counter.values()
        values.sort(key=int)
        if len(values)>int(1):
            print >> outfile, line,
        else:
            pass
    outfile.close()
    infile.close()

def raxml_calculate_base_tree(in_fasta, model, name):
    args = ['raxmlHPC-SSE3', '-f', 'd', '-p', '12345',
	     '-s', '%s' % in_fasta, '-m', '%s' % model, '-n', '%s' % name, "--no-bfgs",
	     '>', '/dev/null 2>&1']
    try:
        vcf_fh = open('raxml.out', 'w')
    except:
        log_isg.logPrint('could not open raxml file')
    try:
        log_fh = open('raxml.log', 'w')
    except:
        log_isg.logPrint('could not open log file')
    try:
        raxml_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
        raxml_run.wait()
    except:
        print "could not infer base pruned tree"
        sys.exit()
    
def file_to_fasta(matrix, out_fasta):
    """almost identical to matrix_to_fasta. Not tested"""
    reduced = [ ]
    out_matrix = open(out_fasta, "w")
    for line in open(matrix, "U"):
        fields = line.strip().split()
        reduced.append(fields)
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_matrix, ">"+str(x[0])
        print >> out_matrix, "".join(x[1:])
    out_matrix.close()

def prune_fasta(to_prune, infile, outfile):
    my_in = open(infile, "U")
    my_out = open(outfile, "w")
    seqrecords = [ ]
    ids = [ ]
    for record in SeqIO.parse(my_in, "fasta"):
        if record.id not in to_prune:
            seqrecords.append(record)
            ids.append(record.id)
    SeqIO.write(seqrecords, my_out, "fasta")
    my_in.close()
    my_out.close()
    return ids
    
    
    
def remove_invariant_sites(in_fasta, out_fasta):
    """only keep invarint sites, doesn't need testing"""
    fasta_to_tab(in_fasta)
    tab_to_matrix("out.tab")
    filter_alignment("tab_matrix")
    file_to_fasta("tab.filtered", out_fasta)

def process_temp_matrices(dist_sets, tree, processors, patristics, insertion_method, parameters, model):
    """needs testing"""
    curr_dir= os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, "*tmp.matrix")):
        """the genome names are parsed out of the tmp.matrices"""
        name=get_seq_name(infile)
        split_fields=name.split(".")
        outfile=open("%s.%s.subsample.distances.txt" % (split_fields[0],split_fields[2]), "a")
        name_fixed = []
        name_fixed.append(re.sub('[:,]', '', split_fields[2]))
        if os.path.isfile("%s.tree" % (split_fields[0]+split_fields[2])):
            pass
        else:
            tmptree = open("%s.tmp.tree" % ((split_fields[0]+split_fields[2])), "w")
            """The genomes in the dist_sets will be pruned in
            the subsample routine"""
            to_prune = []
            for x in dist_sets:
                if x[0] == split_fields[0]: 
                    if x[1] == split_fields[2] or x[2] == split_fields[2]:
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
            tmptree2 = open("%s.tree" % (split_fields[0]+split_fields[2]), "w")
            for line in open("%s.tmp.tree" % (split_fields[0]+split_fields[2]), "U"):
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
            matrix_to_fasta(infile, "%s.fasta" % (split_fields[0]+split_fields[2]))
        except:
            print "problem converting matrix to fasta"
            sys.exit()
        """if problems in the tree names are found, they are removed by the system command"""
        os.system("sed 's/://g' %s.fasta | sed 's/,//g' > %s_in.fasta" % ((split_fields[0]+split_fields[2]),(split_fields[0]+split_fields[2])))
        prune_fasta(to_prune, "%s_in.fasta" % (split_fields[0]+split_fields[2]), "%s_pruned.fasta" % (split_fields[0]+split_fields[2]))
        """Unknown genomes are added to the tree.  If the parameters file has already been created, don't create it again"""
        if os.path.isfile("%s-PARAMS" % (split_fields[0]+split_fields[2])):
            try:
                #model is now hardcoded as GTRGAMMA; this makes it possible to use multiple processors
                run_raxml("%s_in.fasta" % (split_fields[0]+split_fields[2]), "%s.tree" % (split_fields[0]+split_fields[2]), "%s.subsampling_classifications.txt" % (split_fields[0]+split_fields[2]), insertion_method, "%s-PARAMS" % (split_fields[0]+split_fields[2]), "GTRGAMMA", "%s" % (split_fields[0]+split_fields[2]))
            except:
                print "problem adding unknown sequences to the following tree: %s" % (split_fields[0]+split_fields[2])
                os.system("rm RAx*")
                pass
        else:
            try:
                subprocess.check_call("raxmlHPC-SSE3 -f e -m GTRGAMMA -s %s_pruned.fasta -t %s.tree -n %s-PARAMS --no-bfgs > /dev/null 2>&1" % ((split_fields[0]+split_fields[2]),(split_fields[0]+split_fields[2]),(split_fields[0]+split_fields[2])), shell=True)
                os.system("mv RAxML_binaryModelParameters.%s-PARAMS %s-PARAMS" % ((split_fields[0]+split_fields[2]),(split_fields[0]+split_fields[2])))
                run_raxml("%s_in.fasta" % (split_fields[0]+split_fields[2]), "%s.tree" % (split_fields[0]+split_fields[2]), "%s.subsampling_classifications.txt" % (split_fields[0]+split_fields[2]), insertion_method, "%s-PARAMS" % (split_fields[0]+split_fields[2]), "GTRGAMMA", "%s" % (split_fields[0]+split_fields[2]))
            except:
                print "problem adding unknown sequences to the following tree: %s" % (split_fields[0]+split_fields[2])
                os.system("rm RAx*")
                pass
        try:
            calculate_pairwise_tree_dists("%s.tree_including_unknowns_noedges.tree" % (split_fields[0]+split_fields[2]), "%s.resampling_distances.txt" % (split_fields[0]+split_fields[2]))
        except:
            pass
        for line in open("%s.resampling_distances.txt" % (split_fields[0]+split_fields[2]),"U"):
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
        outfile.close()
        
def compare_subsample_results(outnames,distances,fudge):
    """needs testing"""
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
    tree that can be visualized with FigTree - needs testing"""
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
    that can be easily combined with temporary files - tested"""
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
    """these are all of the screened SNPs - tested"""
    matrix_ids=[ ]
    firstLine = in_matrix.readline()
    for line in in_matrix:
        mfields=line.split()
        matrix_ids.append(mfields[0])
    value_dict={}
    new_dicts={}
    with open(vcf, "U") as my_file:
        first_char = my_file.read(1)
        if first_char:
            my_file.seek(0)
            for line in my_file:
                fields=line.split()
                value_dict.update({fields[0]:fields[1]})
    if len(value_dict)>=1:
        for k,v in value_dict.iteritems():
            new_dicts.update({k:v})
        for x in matrix_ids:
            if x not in value_dict.keys():new_dicts.update({x:"-"})
    variety = [ ]
    for x in matrix_ids:
        if x in new_dicts:
            variety.append(new_dicts.get('%s' % x))
    if len(variety)>=1:
        if "A" or "T" or "G" or "C" in variety:
            open("%s.tmp.matrix" % name, 'a').write('%s\n' % name)
            for x in matrix_ids:
                if x in new_dicts:
                    open("%s.tmp.matrix" % name, 'a').write("%s\n" % new_dicts.get('%s' % x))
        else:
            print "sample %s had no usable positions!!!" % name
    else:
        pass
    in_matrix.close()
    return new_dicts

def grab_names():
    """These names will be used for future iterations - tested"""
    curr_dir= os.getcwd()
    outnames = [ ]
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered.vcf")):
        name=get_seq_name(infile)
        reduced=name.replace(".filtered.vcf","")
        outnames.append(reduced)
    return outnames

def parse_likelihoods(infile):
    """This function parses the likelihoods output from RAxML - tested"""
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
    return like_dict

def calculate_pairwise_tree_dists(intree, output):
    """uses dendropy function to calculate all pairwise distances between tree - tested"""
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
    """tested"""
    results = []
    for final_set in final_sets:
        if len(final_set) == 0:
            pass
        results.append(final_set[1]+final_set[2])
    return results
            
def find_two_new(infile,outnames):
    """find two closest genomes to each query genome,
    return the names and distances (sorted), for the
    two genomes - tested"""
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

def subsample_snps_dev(matrix, final_set, used_snps, subnums, allsnps):
    """needs testing"""
    for k,v in used_snps.iteritems():
        if final_set[0]==k:
            for x in range(1,int(subnums)+1):
                kept_snps=random.sample(set(allsnps), int(v))
                solids = set(kept_snps)
                outfile = open("%s.%s.%s.tmp.matrix" % (k,x,final_set[1]), "w")
                in_matrix=open(matrix,"U")
                firstLine = in_matrix.readline()
                print >> outfile, firstLine,
                first_fields = firstLine.split()
                fixed_fields = []
                for y in first_fields:
                    fixed_fields.append(re.sub('[:,]', '', y))
                gindex=fixed_fields.index(final_set[1])
                for line in in_matrix:
                    matrix_fields=line.split()
                    if matrix_fields[0] in solids:
                        print >> outfile, line,
                    else:
                        print >> outfile, "\t".join(matrix_fields[:gindex])+"\t"+"-"+"\t"+"\t".join(matrix_fields[gindex+1:])+"\n",
                in_matrix.close()
                outfile.close()
        else:
            pass

def get_all_snps(matrix):
    """tested"""
    allSNPs = [ ]
    for line in open(matrix, "U"):
        if line.startswith("LocusID"):
            pass
        else:
            fields=line.split()
            allSNPs.append(fields[0])
    return allSNPs

def create_params_files(id, to_prune_set, full_tree, full_matrix, dist_sets, processors):
    """not currently tested, but needs to be"""
    if int(processors)<=2:
        my_processors = 2
    else:
        my_processors = int(int(processors)/2)
    for item in to_prune_set:
        new_name = str(id)+str(item)
        tmptree = open("%s.tmp.tree" % new_name, "w")
        to_prune = []
        for x in dist_sets:
            if x[0] == id: 
                if x[1] == item or x[2] == item:
                    to_prune.append(x[1])
        to_prune_fixed=[]
        for x in to_prune:
            to_prune_fixed.append(re.sub('[:,]', '', x))
        tree_full = dendropy.Tree.get_from_path(full_tree,schema="newick",preserve_underscores=True)
        tree_full.prune_taxa_with_labels(to_prune_fixed)
        final_tree = branch_lengths_2_decimals(tree_full.as_string("newick"))
        print >> tmptree, final_tree
        tmptree.close()
        tmptree2 = open("%s.tree" % new_name, "w")
        for line in open("%s.tmp.tree" % new_name, "U"):
            if line.startswith("[&U]"):
                fields = line.split()
                fixed_fields = [ ]
                for x in fields:
                    fixed_fields.append(x.replace("'",""))
                print >> tmptree2, fixed_fields[1]
            else:
                pass
        tmptree2.close()
        """result is a pruned tree that is ready for RAxML"""
        matrix_to_fasta(full_matrix, "%s.fasta" % new_name)
        os.system("sed 's/://g' %s.fasta | sed 's/,//g' > %s_in.fasta" % (new_name, new_name))
        prune_fasta(to_prune, "%s_in.fasta" % new_name, "%s_pruned.fasta" % new_name)
        try:
            subprocess.check_call("raxmlHPC-PTHREADS-SSE3 -T %s -f e -m GTRGAMMA -s %s_pruned.fasta -t %s.tree -n %s-PARAMS --no-bfgs > /dev/null 2>&1" % (my_processors, new_name, new_name, new_name), shell=True)
            os.system("mv RAxML_binaryModelParameters.%s-PARAMS %s-PARAMS" % (new_name, new_name))
        except:
            pass
    
def process_temp_matrices_dev(dist_sets, sample, tree, processors, patristics, insertion_method, parameters, model):
    """not currently tested, but needs to be"""
    name=get_seq_name(sample)
    split_fields=name.split(".")
    outfile=open("%s.%s.subsample.distances.txt" % (split_fields[0],split_fields[2]), "a")
    name_fixed = []
    name_fixed.append(re.sub('[:,]', '', split_fields[2]))
    to_prune = []
    for x in dist_sets:
        if x[0] == split_fields[0]: 
            if x[1] == split_fields[2] or x[2] == split_fields[2]:
                to_prune.append(x[1])
    to_prune_fixed=[]
    for x in to_prune:
        to_prune_fixed.append(re.sub('[:,]', '', x))
    full_context = split_fields[0]+split_fields[1]+split_fields[2]
    tree_full = dendropy.Tree.get_from_path(tree,schema="newick",preserve_underscores=True)
    tree_full.prune_taxa_with_labels(to_prune_fixed)
    tmptree = open("%s.tmp.tree" % full_context, "w")
    final_tree = branch_lengths_2_decimals(tree_full.as_string("newick"))
    print >> tmptree, final_tree
    tmptree.close()
    tmptree2 = open("%s.tree" % full_context, "w")
    for line in open("%s.tmp.tree" % full_context, "U"):
        if line.startswith("[&U]"):
            fields = line.split()
            fixed_fields = [ ]
            for x in fields:
                fixed_fields.append(x.replace("'",""))
            print >> tmptree2, fixed_fields[1]
        else:
            pass
    tmptree2.close()
    matrix_to_fasta(sample, "%s.fasta" % full_context)
    os.system("sed 's/://g' %s.fasta | sed 's/,//g' > %s_in.fasta" % (full_context, full_context))
    run_raxml("%s_in.fasta" % full_context, "%s.tree" % full_context, "%s.subsampling_classifications.txt" % full_context, insertion_method, "%s-PARAMS" % (split_fields[0]+split_fields[2]), "GTRGAMMA", "%s" % full_context)
    calculate_pairwise_tree_dists("%s.tree_including_unknowns_noedges.tree" % full_context, "%s.resampling_distances.txt" % full_context)
    for line in open("%s.resampling_distances.txt" % full_context,"U"):
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
    outfile.close()
        

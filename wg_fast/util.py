from __future__ import division
import os,os.path
import re
import logging
import subprocess
import sys
from subprocess import Popen
import time
try:
    from Bio import SeqIO
    from Bio import Phylo
except:
    print("BioPython is not in your PYTHONPATH, but needs to be")
    sys.exit()
try:
    import dendropy
    from dendropy import treecalc
except:
    print("dendropy is not installed, but needs to be")
    sys.exit()
import glob
from operator import itemgetter
import threading
import collections
import random

#logPrint stuff
OUTSTREAM = sys.stdout
ERRSTREAM = sys.stderr
DEBUG = False

def logPrint(msg, stream=None):
    if stream is None:
        stream = OUTSTREAM
    stream.write('LOG: %s - %s\n' % (timestamp(), removeRecursiveMsg(msg)))
    stream.flush()

def errorPrint(msg, stream=None):
    if stream is None:
        stream = ERRSTREAM

    stream.write('ERROR: %s - %s\n' % (timestamp(), removeRecursiveMsg(msg)))
    stream.flush()

def debugPrint(fmsg, stream=None):
    """In this case msg is a function, so the work is only done if debugging is one"""
    if DEBUG:
        if stream is None:
            stream = ERRSTREAM

        stream.write('DEBUG: %s - %s\n' % (timestamp(), removeRecursiveMsg(fmsg())))
        stream.flush()

def timestamp():
    return time.strftime('%Y/%m/%d %H:%M:%S')

def removeRecursiveMsg(msg):
    """
    This takes a message and if it starts with something that looks like
    a message generated with these tools it chops it off.  Useful if using
    one of these logging functions to print output from a program using
    the same logging functions
    """
    if msg.startswith('ERROR: ') or msg.startswith('DEBUG: ') or msg.startswith('LOG: '):
        return msg.split(' - ', 1)[1]
    else:
        return msg

def mp_shell(func, params, numProc):
    from multiprocessing import Pool
    p = Pool(numProc)
    out = p.map(func, params)
    p.terminate()
    return out

def report_stats(results,name,output):
    outfile = open(output, "w")
    total_size = []
    mapped_size = []
    with open(results) as infile:
        for line in infile:
            newline = line.strip()
            fields = newline.split()
            """this just makes sure that the file looks correct"""
            if len(fields)==3:
                total_size.append(float(fields[1]))
                mapped_size.append(float(fields[2]))
            else:
                print("coverage file for %s is malformed" % name)
                print("-------------------------")
    try:
        total_summed = sum(total_size)
        total_mapped = sum(mapped_size)
        mapped_value = (total_mapped/total_summed)*100
        outfile.write(str(name)+"\t"+str.format('{0:.4f}',mapped_value)+"\n")
        outfile.close()
    except:
        pass

def merge_files_by_column(column,file_1,file_2,out_file):
    """Takes 2 file and merge their columns based on the column. It is assumed
    that line ordering in the files do not match, so we read both files into memory
    and join them"""
    join_map = {}
    with open(file_1) as my_file_1:
        for line in my_file_1:
            line.strip()
            row = line.split()
            column_value = row.pop(column)
            join_map[column_value] = row
    with open(file_2) as my_file_2:
        for line in my_file_2:
            line.strip()
            row = line.split()
            column_value = row.pop(column)
            if column_value in join_map:
                join_map[column_value].extend(row)
    fout = open(out_file, 'w')
    for k, v in join_map.items():
        fout.write('\t'.join([k] + v) + '\n')
    fout.close()

def sum_coverage(coverage,cov,name):
    outfile = open("%s.amount_covered.txt" % name, "w")
    another_outfile = open("%s.sum_covered.txt" % name, "w")
    all = []
    my_dict = {}
    cov_dict = {}
    with open(coverage) as my_file:
        for line in my_file:
            fields=line.split()
            fields = map(lambda s: s.strip(), fields)
            all.append(fields)
    for x, y in all:
        """Here we're only counting it if it is above the given coverage threshold"""
        try:
            cov_dict[x].append(int(y))
        except KeyError:
            cov_dict[x] = [int(y)]
        if int(y)>=int(cov):
           try:
               my_dict[x].append(y)
           except KeyError:
               my_dict[x] = [y]
        else:
               pass
    for k,v in my_dict.items():
        outfile.write(str(k)+"\t"+str(len(v))+"\n")
    for k,v in cov_dict.items():
        another_outfile.write(str(k)+"\t"+str(sum(v))+"\n")
    outfile.close()
    another_outfile.close()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print("directory of fastqs cannot be found")
        sys.exit()

def test_filter(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("option not supported.  Only select from T and F")
        sys.exit()

def test_methods(option, opt_str, value, parser):
    if "ML" in value:
        setattr(parser.values, option.dest, value)
    elif "MP" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("option not supported.  Only select from MP or ML")
        sys.exit()

def test_models(option, opt_str, value, parser):
    if "GTRGAMMA" in value:
        setattr(parser.values, option.dest, value)
    elif "ASC_GTRGAMMA" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("substitution model is not supported")
        sys.exit()

def get_seq_length(ref, name):
    """uses BioPython in order to calculated the length of
    each fasta entry in the reference fasta"""
    outfile = open("%s.tmp.txt" % name, "w")
    for record in SeqIO.parse(open(ref), "fasta"):
        outfile.write(str(record.id)+"\t"+str(len(record.seq))+"\n")
    outfile.close()

def remove_column(temp_file, name):
    outfile = open("%s.coverage.out" % name, "w")
    my_fields = []
    with open(temp_file) as my_file:
        for line in my_file:
            fields=line.split()
            del fields[1]
            my_fields.append(fields)
    for x in my_fields:
        outfile.write("\t".join(x))
        outfile.write("\n")
    outfile.close()

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
            #m=re.match("(.*)("+"_R1"+")(_.*)$",file_name_before_ext)
            m=re.match("(.*)("+"_R1"+")(.*)$",file_name_before_ext)
            if m is not None:
                (baseName,read) = m.groups()[0], m.groups()[1]
                forward_reads[baseName] = infile
            else:
                #m=re.match("(.*)("+"_R2"+")(_.*)$",file_name_before_ext)
                m=re.match("(.*)("+"_R2"+")(.*)$",file_name_before_ext)
                if m is not None:
                    (baseName,read) = m.groups()[0], m.groups()[1]
                    reverse_reads[baseName] = infile
                else:
                    print("Could not determine forward/reverse read status for input file %s" % infile)
        else:
            baseName, read  = m.groups()[0], m.groups()[3]
            if read == "_R1":
                forward_reads[baseName] = infile
            elif read == "_R2":
                reverse_reads[baseName] = infile
            else:
                print("Could not determine forward/reverse read status for input file")
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

def get_sequence_length(fastq_in):
    from itertools import islice
    from gzip import GzipFile
    with GzipFile("%s" % fastq_in) as file:
        head = list(islice(file, 2))
    return len(head[1])

def _perform_workflow_run_loop_dev(data):
    """sample ID"""
    idx = data[0]
    """path to data"""
    f = data[1]
    dir_path = data[2]
    reference = data[3]
    ref_coords = data[4]
    coverage = data[5]
    proportion = data[6]
    matrix = data[7]
    scratch_dir = data[8]
    doc = data[9]
    tmp_dir = data[10]
    """This is now the PATH to bbduk.sh"""
    wgfast_path = data[11]
    processors = data[12]
    ploidy = data[13]
    if os.path.isfile("%s.tmp.xyx.matrix" % idx):
        logPrint("%s.tmp.xyx.matrix exists, skipping" % idx)
    else:
        """This means that the data is paired end"""
        logPrint("processing %s" % idx)
        if len(f)>1:
            if os.path.isfile("%s.F.paired.fastq.gz" % idx):
                pass
            else:
                length = int(get_sequence_length(f[0])/2)
                try:
                    subprocess.check_call("bbduk.sh in=%s in2=%s ref=%s/bin/illumina_adapters_all.fasta out=%s.F.paired.fastq.gz out2=%s.R.paired.fastq.gz minlen=%s overwrite=true > /dev/null 2>&1" % (f[0],f[1],wgfast_path,idx,idx,length), shell=True)
                except:
                    print("Read trimmer did not finish correctly")
                    sys.exit()
            #Skip this step, which can be long, if the bam file already exists
            if os.path.isfile("%s_renamed.bam" % idx):
                pass
            else:
                subprocess.check_call("minimap2 -ax sr %s/reference.fasta %s.F.paired.fastq.gz %s.R.paired.fastq.gz | samtools sort -l 0 -@ %s - | samtools view -F 4 -Su -o %s_renamed.bam -" % (scratch_dir,idx,idx,processors,idx),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        else:
            """Single end read support"""
            length = int(get_sequence_length(f[0])/2)
            try:
                subprocess.check_call("bbduk.sh -Xmx2g in=%s ref=%s/bin/illumina_adapters_all.fasta out=%s.single.fastq.gz minlen=%s overwrite=true ignorebadquality" % (f[0],wgfast_path,idx,length), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
            except:
                print("Read trimmer did not finish correctly")
                sys.exit()
            subprocess.check_call("minimap2 -ax sr %s/reference.fasta %s.single.fastq.gz | samtools sort -l 0 -@ %s - | samtools view -F 4 -Su -o %s_renamed.bam -" % (scratch_dir,idx,processors,idx),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        """inserts read group information, required by new versions of GATK"""
        try:
            subprocess.check_call("picard AddOrReplaceReadGroups I=%s_renamed.bam O=%s_renamed_header.bam SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=illumina RGSM=%s RGPU=name CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT" % (idx,idx,idx,idx,idx),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        except:
            time.sleep(5)
            subprocess.check_call("picard AddOrReplaceReadGroups I=%s_renamed.bam O=%s_renamed_header.bam SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=illumina RGSM=%s RGPU=name CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT" % (idx,idx,idx,idx,idx),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        try:
            subprocess.check_call("samtools index %s_renamed_header.bam > /dev/null 2>&1" % idx,stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        except:
            print("problem indexing file with samtools..exiting")
            sys.exit()
        #TODO: be able to change the ploidy, which will affect the number of calls
        subprocess.check_call("gatk HaplotypeCaller -R %s/reference.fasta -I %s_renamed_header.bam -O %s.test.vcf -ERC BP_RESOLUTION -ploidy %s" % (scratch_dir,idx,idx,ploidy), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        subprocess.check_call("gatk GenotypeGVCFs -O %s.vcf.out -R %s/reference.fasta -V %s.test.vcf --include-non-variant-sites" % (idx,scratch_dir,idx), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        if "T" == doc:
            subprocess.check_call("samtools depth -aa %s_renamed_header.bam > %s.coverage" % (idx,idx), shell=True)
            remove_column("%s.coverage" % idx, idx)
            sum_coverage("%s.coverage.out" % idx, coverage, idx)
            #If the amount_covered files are empty, WG-FAST will hang here. I put in this check and need to test
            if os.path.getsize("%s.amount_covered.txt" % idx) == 0:
                pass
            else:
                merge_files_by_column(0,"ref.genome_size.txt", "%s.amount_covered.txt" % idx, "%s.results.txt" % idx)
                merge_files_by_column(0,"ref.genome_size.txt", "%s.sum_covered.txt" % idx, "%s.cov.results.txt" % idx)
                report_stats("%s.results.txt" % idx, idx, "%s_breadth.txt" % idx)
                report_stats("%s.cov.results.txt" % idx, idx, "%s_sum_cov.txt" % idx)
        else:
            pass
        #filtered.vcf will be created in this function
        good_calls = process_vcf("%s.vcf.out" % idx, ref_coords, coverage, proportion, idx)
        if int(good_calls) > 0:
            make_temp_matrix("%s.filtered.vcf" % idx, matrix, idx)
        else:
            print("sample %s had no SNP calls and will not be inserted into the tree" % idx)
            print("-------------------------")

def run_loop_dev(fileSets,dir_path,reference,processors,ref_coords,coverage,proportion,
    matrix,scratch_dir,doc,tmp_dir,wgfast_path,ploidy):
    files_and_temp_names = []
    for idx, f in fileSets.items():
        files_and_temp_names.append([idx,f,dir_path,reference,ref_coords,coverage,proportion,
                                        matrix,scratch_dir,doc,tmp_dir,wgfast_path,processors,ploidy])
    mp_shell(_perform_workflow_run_loop_dev,files_and_temp_names,processors)

def process_vcf(vcf, ref_coords, coverage, proportion, name):
    """finds SNPs that pass user-defined thresholds
    for coverage and proportion - needs to look at tests"""
    vcf_out = open("%s.filtered.vcf" % name, "w")
    good_snps = []
    mixed_snps = []
    mixed_refs = []
    ref_set = set(ref_coords)
    with open(vcf) as vcf_in:
        for line in vcf_in:
            newline=line.strip("\n")
            if newline.startswith('#'):
               pass
            elif newline.startswith('Java'):
                pass
            elif newline.startswith('/tmp'):
                pass
            elif newline.startswith('Try'):
                pass
            elif newline.startswith('INFO'):
                pass
            elif newline in ['\n', '\r\n']:
                pass
            else:
                fields=newline.split()
                """for GATK, a period signifies a reference call.
                First we want to look at the situation where this is
                not the case"""
                merged_fields=fields[0]+"::"+fields[1]
                if merged_fields in ref_set:
                    #This indicates that the position is a SNP
                    if "." != fields[4] and len(fields[4]) == 1 and len(fields[3])==1 and fields[8]!="GT:AD:PGT:PID:PS":
                        if fields[6] == "LowQual":
                            pass
                        else:
                            snp_fields=fields[9].split(':')
                            if int(len(snp_fields))>2:
                                prop_fields=snp_fields[1].split(',')
                                fixed_coverage = int(float(snp_fields[2]))
                                if fixed_coverage>=coverage:
                                    if int(prop_fields[1])/int(snp_fields[2])>=float(proportion):
                                        vcf_out.write(fields[0]+"::"+fields[1]+"\t"+fields[4]+"\n")
                                        good_snps.append("1")
                                    else:
                                        #Changed out a gap character with an N
                                        vcf_out.write(fields[0]+"::"+fields[1]+"\t"+"N"+"\n")
                                        mixed_snps.append("1")
                                    """if problems are encountered, throw in a gap.  Could be too conservative"""
                                else:
                                    vcf_out.write(fields[0]+"::"+fields[1]+"\t"+"N"+"\n")
                            else:
                                pass
                    #This is if the position is reference
                    elif "." in fields[4] and len(fields[4])==1 and len(fields[3])==1:
                        if fields[6] == "LowQual":
                            pass
                        else:
                            if len(fields) == 10:
                                if "DP" in fields[7]:
                                    nosnp_fields=fields[7].split(';')
                                    #This will provide the coverage
                                    if "DP" in nosnp_fields[0]:
                                        cov_fields=nosnp_fields[0].replace("DP=","")
                                    elif "DP" in nosnp_fields[1]:
                                        cov_fields=nosnp_fields[1].replace("DP=","")
                                    fixed_coverage = int(float(cov_fields))
                                    if fixed_coverage>=coverage:
                                        vcf_out.write(fields[0]+"::"+fields[1]+"\t"+fields[3]+"\n")
                                    else:
                                        vcf_out.write(fields[0]+"::"+fields[1]+"\t"+"N"+"\n")
                                        mixed_refs.append("1")
                                else:
                                    vcf_out.write(fields[0]+"::"+fields[1]+"\t"+"N"+"\n")
                                    mixed_refs.append("1")
                    #If it can't determine the status of the position, add an N
                    else:
                        vcf_out.write(fields[0]+"::"+fields[1]+"\t"+"N"+"\n")
    vcf_out.close()
    print("number of SNPs in genome %s = %s" % (name, str(len(good_snps))))
    print("number of discarded SNPs in genome %s = %s" % (name, str(len(mixed_snps))))
    print("number of discarded Reference positions in genome due to no coverage %s = %s" % (name, str(len(mixed_refs))))
    print("-------------------------")
    return len(good_snps)

def sort_information(x):
    """simple sort - tested"""
    try:
        fields = x.split("::")
        return int(fields[1])
    except:
        raise TypeError("problem encountered parsing fields")

def matrix_to_fasta(matrix_in, outfile):
    """function to convert a SNP matrix to a multi-fasta file - tested"""
    reduced = []
    out_fasta = open(outfile, "w")
    redux = []
    with open(matrix_in) as my_file:
        for line in my_file:
            newline = line.strip()
            fields = newline.split()
            reduced.append(fields[1:])
    test=map(list, zip(*reduced))
    for x in test:
        out_fasta.write(">"+str(x[0])+"\n")
        out_fasta.write("".join(x[1:]))
        out_fasta.write("\n")
        redux.append(">"+str(x[0])+"".join(x[1:3]))
    out_fasta.close()
    return redux

def run_raxml(fasta_in, tree, out_class_file, insertion_method, parameters, model, suffix):
    """untested function, system calls"""
    if "NULL" == parameters:
        if "ASC_GTRGAMMA" == model:
            args = ['raxmlHPC-SSE3', '-f', '%s' % insertion_method,
	         '-s', '%s' % fasta_in, '-m', '%s' % model, '-n', '%s' % suffix, '-t',
	         '%s' % tree, '--asc-corr=lewis', '--no-bfgs', '>', '/dev/null 2>&1']
        else:
            args = ['raxmlHPC-SSE3', '-f', '%s' % insertion_method,
	         '-s', '%s' % fasta_in, '-m', '%s' % model, '-n', '%s' % suffix, '-t',
	         '%s' % tree, '--no-bfgs', '>', '/dev/null 2>&1']
    else:
        if "ASC_GTRGAMMA" == model:
            args = ['raxmlHPC-SSE3', '-f', '%s' % insertion_method,
	         '-s', '%s' % fasta_in, '-m', '%s' % model, '-n', '%s' % suffix, '-R', parameters, '-t',
	         '%s' % tree, '--asc-corr=lewis', '--no-bfgs', '>', '/dev/null 2>&1']
        else:
            args = ['raxmlHPC-SSE3', '-f', '%s' % insertion_method,
	         '-s', '%s' % fasta_in, '-m', '%s' % model, '-n', '%s' % suffix, '-R', parameters, '-t',
	         '%s' % tree, '--no-bfgs', '>', '/dev/null 2>&1']
    try:
        vcf_fh = open('%s.raxml.out' % suffix, 'w')
    except:
        logPrint('could not open raxml file')
    try:
        log_fh = open('%s.raxml.log' % suffix, 'w')
    except:
        logPrint('could not open log file')
    #TODO: See why this doesn't catch many RAxML errors
    try:
        raxml_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
        raxml_run.wait()
    except:
        logPrint("sequence(s) were not inserted into tree!!!!!")
    try:
        os.system("sed 's/\[[^]]*\]//g' RAxML_labelledTree.%s > %s.tree_including_unknowns_noedges.tree" % (suffix, suffix))
        subprocess.check_call("mv RAxML_labelledTree.%s %s_tree_including_unknowns_edges.tree" % (suffix, suffix) , shell=True)
    except:
        logPrint("error encountered running RAxML...check log file")
        sys.exit()
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
    with open(matrix) as my_file:
        for line in my_file:
            if line.startswith("LocusID"):
                pass
            else:
                fields=line.split()
                allSNPs.append(fields[0])
    for k,v in used_snps.items():
        for z in dist_sets:
            if len(z) == 0:
                pass
            else:
                if z[0]==k:
                    for x in range(1,int(subnums)+1):
                        kept_snps=random.sample(set(allSNPs), int(v))
                        outfile = open("%s.%s.%s.tmp.matrix" % (k,x,z[1]), "w")
                        in_matrix=open(matrix)
                        firstLine = in_matrix.readline()
                        outfile.write(firstLine)
                        outfile.write("\n")
                        first_fields = firstLine.split()
                        fixed_fields = []
                        for y in first_fields:
                            fixed_fields.append(re.sub('[:,]', '', y))
                        gindex=fixed_fields.index(z[1])
                        for line in in_matrix:
                            matrix_fields=line.split()
                            if matrix_fields[0] in kept_snps:
                                outfile.write(line)
                                outfile.write("\n")
                            else:
                                outfile.write("\t".join(matrix_fields[:gindex])+"\t"+"-"+"\t"+"\t".join(matrix_fields[gindex+1:])+"\n")
                        in_matrix.close()
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
        with open(infile) as my_file:
            for line in my_file:
                fields=line.split()
                try:
                    if fields[1] != "N":
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
    outfile = open("out.tab", "w")
    with open(fasta) as my_fasta:
        for record in SeqIO.parse(my_fasta,"fasta"):
            """this list is just for testing,
            and is ok if it's overwritten for each
            fasta"""
            for_test = []
            outfile.write(str(record.id)+"\t"+str(record.seq)+"\n")
            for_test.append(record.id)
            for_test.append(str(record.seq))
    outfile.close()
    return for_test

def tab_to_fasta(new_tab):
    """tested"""
    outfile = open("out.fasta", "w")
    with open(new_tab) as infile:
        for line in infile:
            to_test = []
            fields = line.split()
            outfile.write(">"+fields[0]+"\n")
            outfile.write(fields[1].upper())
            to_test.append(fields[0])
            to_test.append(fields[1].upper())
    outfile.close()
    return to_test

def tab_to_matrix(tab):
    """tested"""
    reduced = []
    out_matrix = open("tab_matrix", "w")
    with open(tab) as my_file:
        for line in my_file:
            tmp_list = []
            fields = line.split()
            tmp_list.append(fields[0])
            for nucs in fields[1]:
                tmp_list.append(nucs.upper())
            reduced.append(tmp_list)
    test=list(map(list,zip(*reduced)))
    to_return = []
    for x in test:
        out_matrix.write("\t".join(x))
        out_matrix.write("\n")
        to_return.append(list(x))
    out_matrix.close()
    return to_return

def filter_alignment(tab):
    """tested"""
    outfile = open("tab.filtered", "w")
    with open(tab) as infile:
        firstLine = infile.readline()
        test_fields = []
        outfile.write(firstLine)
        outfile.write("\n")
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
            #sorted_values = values.sort(key=int)
            sorted_values = sorted(values)
            if len(sorted_values)>int(1):
                outfile.write(line)
                outfile.write("\n")
            else:
                pass
            test_fields.append(valid_fields)
    outfile.close()
    return test_fields

def raxml_calculate_base_tree(in_fasta, model, name):
    """not tested, all system calls"""
    args = ['raxmlHPC-SSE3', '-f', 'd', '-p', '12345',
	     '-s', '%s' % in_fasta, '-m', '%s' % model, '-n', '%s' % name, "--no-bfgs",
	     '>', '/dev/null 2>&1']
    try:
        vcf_fh = open('raxml.out', 'w')
    except:
        logPrint('could not open raxml file')
    try:
        log_fh = open('raxml.log', 'w')
    except:
        logPrint('could not open log file')
    try:
        raxml_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
        raxml_run.wait()
    except:
        print("could not infer base pruned tree")
        sys.exit()

def file_to_fasta(matrix, out_fasta):
    """almost identical to matrix_to_fasta. Not tested"""
    reduced = [ ]
    out_matrix = open(out_fasta, "w")
    with open(matrix) as my_file:
        for line in my_file:
            fields = line.strip().split()
            reduced.append(fields)
    test=map(list, zip(*reduced))
    for x in test:
        out_matrix.write(">"+str(x[0])+"\n")
        out_matrix.write("".join(x[1:]))
        out_matrix.write("\n")
    out_matrix.close()

def prune_fasta(to_prune, infile, outfile):
    """tested"""
    my_out = open(outfile, "w")
    seqrecords = [ ]
    ids = [ ]
    with open(infile) as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
            if record.id not in to_prune:
                seqrecords.append(record)
                ids.append(record.id)
    SeqIO.write(seqrecords, my_out, "fasta")
    my_out.close()
    return ids

def remove_invariant_sites(in_fasta, out_fasta):
    """only keep invarint sites, all functions are tested"""
    fasta_to_tab(in_fasta)
    tab_to_matrix("out.tab")
    filter_alignment("tab_matrix")
    file_to_fasta("tab.filtered", out_fasta)

def compare_subsample_results(outnames,distances,fudge):
    """needs testing"""
    curr_dir= os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, "*.subsample.distances.txt")):
        name=get_seq_name(infile)
        split_fields=name.split(".")
        genomes_used = []
        all_dists=[]
        dists_greater_than_true=[]
        dists_equal_to_true=[]
        dists_less_than_true=[]
        """Here, I get the list of the genomes that are being analyzed"""
        try:
            with open(infile) as my_file:
                for line in my_file:
                    fields = line.split()
                    all_dists.append(float(fields[7]))
        except:
           print("problem parsing input file: %s" % infile)
        if len(all_dists)>=1:
            max_dist=max(all_dists)
            print("")
            print("maximum subsample distance between %s and %s = %.2f" % (fields[3],fields[5],float(max_dist)))
        else:
            print("problem grabbing distances- make sure that subsample.distances.txt files aren't empty")
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
        try:
            greaters = int(len(dists_greater_than_true))
            equals = int(len(dists_equal_to_true))
            lessers = int(len(dists_less_than_true))
            print("True distance between Reference and %s = %.2f" % (split_fields[1],float(true_dists[0])))
            print("Sample: %s" % split_fields[0])
            print("Subsample distances between Reference and %s greater than true value = %s" % (split_fields[2],greaters))
            print("Subsample distances between Reference and %s equal to true value = %s" % (split_fields[2],equals))
            print("Subsample distances between Reference and %s less than true value = %s" % (split_fields[2],lessers))
            p = (greaters+lessers)/(greaters+lessers+equals)
            print("Placement p value = %.3f" % float(p))
            print("-------------------------")
        except:
            pass

def transform_tree(tree):
    """converts a Newick tree into a Nexus-formatted
    tree that can be visualized with FigTree - needs testing"""
    #infile = open(tree, "U")
    tree_string = []
    with open(tree) as infile:
        for line in infile:
            tree_string.append(line)
    #infile.close()
    mytree = Phylo.read(tree, 'newick')
    tree_names = [ ]
    for clade in mytree.find_clades():
        if clade.name:
            tree_names.append(clade.name)
    outfile = open("transformed.tree", "w")
    outfile.write("#NEXUS"+"\n")
    outfile.write("begin taxa;"+"\n")
    outfile.write("\t"+"dimensions ntax="+str(len(tree_names))+";"+"\n")
    outfile.write("\t"+"taxlabels"+"\n")
    for clade in mytree.find_clades():
        if clade.name:
            tree_names.append(clade.name)
    nr=[x for i, x in enumerate(tree_names) if x not in tree_names[i+1:]]
    for tree_name in nr:
        if "QUERY__" in tree_name:
            outfile.write("\t"+"'%s'" % tree_name+"[&!color=#-3407821]"+"\n")
        else:
            outfile.write("\t"+tree_name+"\n")
    outfile.write(";"+"\n")
    outfile.write("end;"+"\n")
    outfile.write(""+"\n")
    outfile.write("begin trees;"+"\n")
    for x in tree_string:
        outfile.write("\t"+"tree tree_1 = [&R] "+str(x)+"\n")
    outfile.write("end;"+"\n")
    infile.close()
    outfile.close()

def write_reduced_matrix(matrix):
    """This function takes a NASP formatted
    SNP matrix and writes a temporary matrix
    that can be easily combined with temporary files - tested"""
    in_matrix = open(matrix)
    outfile = open("temp.matrix", "w")
    outdata = [ ]
    firstLine = in_matrix.readline()
    first_fields=firstLine.split()
    last=first_fields.index("#SNPcall")
    outfile.write("\t".join(first_fields[:last]))
    outfile.write("\n")
    for line in in_matrix:
        fields = line.split()
        outfile.write("\t".join(fields[:last]))
        outfile.write("\n")
        outdata.append(len(fields[:last]))
    outfile.close()
    in_matrix.close()
    return outdata

def make_temp_matrix(vcf, matrix, name):
    """these are all of the screened SNPs - tested"""
    matrix_ids=[]
    with open(matrix) as in_matrix:
        firstLine = in_matrix.readline()
        for line in in_matrix:
            newline = line.strip("\n")
            mfields=newline.split()
            matrix_ids.append(mfields[0])
    value_dict={}
    new_dicts={}
    with open(vcf) as my_file:
        for line in my_file:
            newline = line.strip("\n")
            if newline.startswith("#"):
                pass
            else:
                fields=newline.split()
                value_dict.update({fields[0]:fields[1]})
                #value_dict could only contain Ns
    if len(value_dict)>=1:
        """here is where value_dict_set is populated"""
        keys = []
        for k,v in value_dict.items():
            keys.append(k)
        value_dict_set = set(keys)
        new_dicts = value_dict
        for x in matrix_ids:
            if x not in value_dict_set:new_dicts.update({x:"N"})
    else:
        print("no usable information in vcf file, did you use the correct reference?")
    #new_dicts contains all calls, good and bad
    value_dict = {}
    """variety will contain a complete set of SNPs"""
    variety = []
    #This ensures that the values are in order
    for x in matrix_ids:
        for k,v in new_dicts.items():
            if x == k:
                variety.append(v)
    variety_set = set(variety)
    #Changed this from a to w on 12/5/19
    out_matrix = open("%s.tmp.xyx.matrix" % name,"w")
    if len(variety)>=1:
        if "A" or "T" or "G" or "C" in value_dict_set:
            out_matrix.write('%s\n' % name)
            for x in matrix_ids:
                if x in new_dicts:
                    out_matrix.write("%s\n" % new_dicts.get('%s' % x))
        else:
            print("sample %s had no usable positions!!!" % name)
    else:
        print("sample %s has problems" % name)
        print("-------------------------")
    out_matrix.close()
    value_dict_set = []
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
    like_dict = {}
    with open(infile) as my_file:
        for line in my_file:
            fields = line.split()
            try:
                like_dict[fields[0]].append(fields[2])
            except KeyError:
                like_dict[fields[0]] = [fields[2]]
    print("sample_name"+"\t"+"insertion_likelihood"+"\t"+"number of potential insertion nodes")
    for k,v in like_dict.items():
        print(k+"\t"+v[0]+"\t"+str(len(v)))
    return like_dict

def calculate_pairwise_tree_dists(intree, output):
    """uses dendropy function to calculate all pairwise distances between tree - tested"""
    tree = dendropy.Tree.get_from_path(intree, "newick", preserve_underscores=True)
    outfile = open("%s" % output, "w")
    distances = tree.phylogenetic_distance_matrix()
    distances_sets = [ ]
    for i,t1 in enumerate(tree.taxon_namespace):
        for t2 in tree.taxon_namespace[i+1:]:
            distances_sets.append(distances(t1, t2))
    try:
        for i, t1 in enumerate(tree.taxon_namespace):
            for t2 in tree.taxon_namespace[i+1:]:
                outfile.write("Distance between '%s' and '%s': %s\n" % (t1.label, t2.label, distances(t1, t2)))
    except:
        print("problem iterating through tree. Tree is empty or not Newick format")
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
        with open(infile) as my_file:
            for line in my_file:
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

def _perform_workflow_subsample_snps_dev(data):
    final_set = data[0]
    """final_set is a tuple of sample, NN, value"""
    used_snps = data[1]
    """used_snps is a dictionary"""
    subnums = data[2]
    """subnums is an integer"""
    allsnps = data[3]
    """allsnps is a list of all SNPs used"""
    matrix = data[4]
    for k,v in used_snps.items():
        if final_set[0]==k:
            for x in range(1,int(subnums)+1):
                kept_snps=random.sample(set(allsnps), int(v))
                solids = set(kept_snps)
                if os.path.isfile("%s.%s.%s.tmp.matrix" % (k,x,final_set[1])):
                    pass
                else:
                    outfile = open("%s.%s.%s.tmp.matrix" % (k,x,final_set[1]), "w")
                    with open(matrix) as my_matrix:
                        firstLine = my_matrix.readline()
                        outfile.write(firstLine)
                        first_fields = firstLine.split()
                        fixed_fields = []
                        for y in first_fields:
                            fixed_fields.append(re.sub('[:,]', '', y))
                        gindex=fixed_fields.index(final_set[1])
                        for line in my_matrix:
                            matrix_fields=line.split()
                            if matrix_fields[0] in solids:
                                outfile.write(line)
                            else:
                                outfile.write("\t".join(matrix_fields[:gindex])+"\t"+"-"+"\t"+"\t".join(matrix_fields[gindex+1:])+"\n")
                    outfile.close()
        else:
            pass

def subsample_snps_2(final_sets,used_snps,subnums,allsnps,processors,matrix):
    files_and_temp_names = []
    for f in final_sets:
        files_and_temp_names.append([f,used_snps,subnums,allsnps,matrix])
    mp_shell(_perform_workflow_subsample_snps_dev, files_and_temp_names, processors)

def subsample_snps_dev(matrix, final_set, used_snps, subnums, allsnps):
    """needs testing"""
    for k,v in used_snps.items():
        if final_set[0]==k:
            for x in range(1,int(subnums)+1):
                kept_snps=random.sample(set(allsnps), int(v))
                solids = set(kept_snps)
                if os.path.isfile("%s.%s.%s.tmp.matrix" % (k,x,final_set[1])):
                    pass
                else:
                    outfile = open("%s.%s.%s.tmp.matrix" % (k,x,final_set[1]), "w")
                    in_matrix=open(matrix)
                    firstLine = in_matrix.readline()
                    outfile.write(firstLine)
                    first_fields = firstLine.split()
                    fixed_fields = []
                    for y in first_fields:
                        fixed_fields.append(re.sub('[:,]', '', y))
                    gindex=fixed_fields.index(final_set[1])
                    for line in in_matrix:
                        matrix_fields=line.split()
                        if matrix_fields[0] in solids:
                            outfile.write(line)
                        else:
                            outfile.write("\t".join(matrix_fields[:gindex])+"\t"+"-"+"\t"+"\t".join(matrix_fields[gindex+1:])+"\n")
                    in_matrix.close()
                    outfile.close()
        else:
            pass

def get_all_snps(matrix):
    """tested"""
    allSNPs = [ ]
    with open(matrix) as my_matrix:
        for line in my_matrix:
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
        if os.path.isfile("%s-PARAMS" % new_name):
            continue
        else:
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
            tmptree.write(final_tree)
            tmptree.close()
            tmptree2 = open("%s.tree" % new_name, "w")
            with open("%s.tmp.tree" % new_name) as my_file:
                for line in my_file:
                    if line.startswith("[&U]"):
                        fields = line.split()
                        fixed_fields = [ ]
                        for x in fields:
                            fixed_fields.append(x.replace("'",""))
                        tmptree2.write(fixed_fields[1])
                    else:
                        pass
            tmptree2.close()
            """result is a pruned tree that is ready for RAxML"""
            matrix_to_fasta(full_matrix, "%s.fasta" % new_name)
            os.system("sed 's/://g' %s.fasta | sed 's/,//g' > %s_in.fasta" % (new_name, new_name))
            prune_fasta(to_prune, "%s_in.fasta" % new_name, "%s_pruned.fasta" % new_name)
            try:
                subprocess.check_call("rm RAxML*%s-PARAMS" % new_name, shell=True, stderr=open(os.devnull, 'w'))
            except:
                pass
            #forcing GTRGAMMA
            try:
                subprocess.check_call("raxmlHPC-PTHREADS-SSE3 -T %s -f e -m GTRGAMMA -s %s_pruned.fasta -t %s.tree -n %s-PARAMS > /dev/null 2>&1" % (my_processors, new_name, new_name, new_name), shell=True)
                os.system("mv RAxML_binaryModelParameters.%s-PARAMS %s-PARAMS" % (new_name, new_name))
            except:
                continue

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
    #full_context includes sample, near neighbor, and replicate
    full_context = split_fields[0]+split_fields[1]+split_fields[2]
    new_name = split_fields[0]+split_fields[2]
    if os.path.isfile("%s.tree_including_unknowns_noedges.tree" % full_context):
        print("tree already present, skipping")
        pass
    else:
        tree_full = dendropy.Tree.get_from_path(tree,schema="newick",preserve_underscores=True)
        tree_full.prune_taxa_with_labels(to_prune_fixed)
        tmptree = open("%s.tmp.tree" % full_context, "w")
        final_tree = branch_lengths_2_decimals(tree_full.as_string("newick"))
        tmptree.write(final_tree)
        tmptree.close()
        tmptree2 = open("%s.tree" % full_context, "w")
        #for line in open("%s.tmp.tree" % full_context, "U"):
        with open("%s.tmp.tree" % full_context) as my_file:
            for line in my_file:
                if line.startswith("[&U]"):
                    fields = line.split()
                    fixed_fields = [ ]
                    for x in fields:
                        fixed_fields.append(x.replace("'",""))
                    tmptree2.write(fixed_fields[1])
                else:
                    pass
        tmptree2.close()
        matrix_to_fasta(sample, "%s.fasta" % full_context)
        os.system("sed 's/://g' %s.fasta | sed 's/,//g' > %s_in.fasta" % (full_context, full_context))
        if os.path.isfile("%s-PARAMS" % new_name):
            try:
                run_raxml("%s_in.fasta" % full_context, "%s.tree" % full_context, "%s.subsampling_classifications.txt" % full_context, insertion_method, "%s-PARAMS" % new_name, "GTRGAMMA", "%s" % full_context)
            except:
                pass
        else:
            try:
                subprocess.check_call("raxmlHPC-PTHREADS-SSE3 -T %s -f e -m GTRGAMMA -s %s_pruned.fasta -t %s.tree -n %s-PARAMS --no-bfgs > /dev/null 2>&1" % (processors, new_name, new_name, new_name), shell=True)
                os.system("mv RAxML_binaryModelParameters.%s-PARAMS %s-PARAMS" % (new_name, new_name))
                run_raxml("%s_in.fasta" % full_context, "%s.tree" % full_context, "%s.subsampling_classifications.txt" % full_context, insertion_method, "%s-PARAMS" % (split_fields[0]+split_fields[2]), "GTRGAMMA", "%s" % full_context)
            except:
                pass
    if os.path.isfile("%s.resampling_distances.txt" % full_context):
        pass
    else:
        try:
            calculate_pairwise_tree_dists("%s.tree_including_unknowns_noedges.tree" % full_context, "%s.resampling_distances.txt" % full_context)
            for line in open("%s.resampling_distances.txt" % full_context,"U"):
                resample_fields = line.split()
                myid = re.sub("[:']", "",resample_fields[4])
                fixedid = myid.replace("QUERY___","")
                newid = re.sub("[:']","",resample_fields[2])
                fixedid2 = newid.replace("QUERY___","")
                if resample_fields[2] == "'Reference'" and fixedid in name_fixed:
                    outfile.write("resampled distance between Reference and %s = %s\n" % (fixedid, resample_fields[5]))
                elif resample_fields[4] == "'Reference':" and fixedid2 in name_fixed:
                    outfile.write("resampled distance between Reference and %s = %s\n" % (fixedid2, resample_fields[5]))
                else:
                    pass
        except:
            pass

def check_input_files(matrix,reference):
    ref_names = []
    with open(reference) as my_ref:
        for record in SeqIO.parse(my_ref,"fasta"):
            ref_names.append(record.id)
    with open(matrix) as f:
        for line in f.readlines()[:2]:
            if line.startswith("LocusID"):
                pass
            else:
                fields = line.split()
                name_fields=fields[0].split("::")
                if name_fields[0] in ref_names:
                    pass
                else:
                    print("The IDs in your Reference don't match the names in your SNP matrix! Please fix and re-start...exiting...")
                    sys.exit()

def create_merged_vcf():
    out_file = open("merged.vcf", "w")
    start_dir = os.getcwd()
    lists = []
    for infile in glob.glob(os.path.join(start_dir, "*.tmp.xyx.matrix")):
        data = open(infile).read().splitlines()
        lists.append(data)
    test=map(list, zip(*lists))
    for x in test:
        out_file.write("\t".join(x))
        out_file.write("\n")
    out_file.close()

def _perform_workflow_create_params(data):
    id = data[0]
    to_prune_set = data[1]
    full_tree = data[2]
    full_matrix = data[3]
    dist_sets = data[4]
    processors = data[5]
    if int(processors)<=2:
        my_processors = 2
    else:
        my_processors = int(int(processors)/2)
    for item in to_prune_set:
        new_name = str(id)+str(item)
        if os.path.isfile("%s-PARAMS" % new_name):
            pass
        else:
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
            tmptree.write(final_tree)
            tmptree.close()
            tmptree2 = open("%s.tree" % new_name, "w")
            with open("%s.tmp.tree" % new_name) as my_file:
                for line in my_file:
                    line = line.replace("'","")
                    tmptree2.write(line)
            tmptree2.close()
            """result is a pruned tree that is ready for RAxML"""
            matrix_to_fasta(full_matrix, "%s.fasta" % new_name)
            os.system("sed 's/://g' %s.fasta | sed 's/,//g' > %s_in.fasta" % (new_name, new_name))
            prune_fasta(to_prune,"%s_in.fasta" % new_name,"%s_pruned.fasta" % new_name)
            try:
                subprocess.check_call("rm RAxML*%s-PARAMS" % new_name, shell=True, stderr=open(os.devnull, 'w'))
            except:
                pass
            #forcing GTRGAMMA
            try:
                subprocess.check_call("raxmlHPC-PTHREADS-SSE3 -T %s -f e -m GTRGAMMA -s %s_pruned.fasta -t %s.tree -n %s-PARAMS > /dev/null 2>&1" % (my_processors, new_name, new_name, new_name), shell=True)
                os.system("mv RAxML_binaryModelParameters.%s-PARAMS %s-PARAMS" % (new_name, new_name))
            except:
                print("problem creating PARAMS file for %s.tree, exiting" % new_name)
                sys.exit()

def create_params_files_dev(new_sample_dicts,tree,matrix,final_sets,processors):
    to_run = []
    for k,v in new_sample_dicts.items():
        to_run.append([k,v,tree,matrix,final_sets,processors])
    mp_shell(_perform_workflow_create_params,to_run,processors)

def process_temp_matrices_2(final_sets,final_matrices,tree,processors,patristic_distances, V, parameters, model):
    to_run = []
    for matrix in final_matrices:
        to_run.append([final_sets,matrix,tree,processors,patristic_distances,V,parameters,model])
    mp_shell(_perform_workflow_temp_matrices, to_run, processors)

def _perform_workflow_temp_matrices(data):
    dist_sets = data[0]
    sample = data[1]
    tree = data[2]
    processors = data[3]
    patristics = data[4]
    insertion_method = data[5]
    parameters = data[6]
    model = data[7]
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
    #full_context includes sample, near neighbor, and replicate
    full_context = split_fields[0]+split_fields[1]+split_fields[2]
    new_name = split_fields[0]+split_fields[2]
    if os.path.isfile("%s.tree_including_unknowns_noedges.tree" % full_context):
        print("tree already present, skipping")
        pass
    else:
        tree_full = dendropy.Tree.get_from_path(tree,schema="newick",preserve_underscores=True)
        tree_full.prune_taxa_with_labels(to_prune_fixed)
        tmptree = open("%s.tmp.tree" % full_context, "w")
        final_tree = branch_lengths_2_decimals(tree_full.as_string("newick"))
        tmptree.write(final_tree)
        tmptree.close()
        tmptree2 = open("%s.tree" % full_context, "w")
        with open("%s.tmp.tree" % full_context) as my_file:
            for line in my_file:
                line = line.replace("'","")
                tmptree2.write(line)
        tmptree2.close()
        matrix_to_fasta(sample, "%s.fasta" % full_context)
        os.system("sed 's/://g' %s.fasta | sed 's/,//g' > %s_in.fasta" % (full_context, full_context))
        if os.path.isfile("%s-PARAMS" % new_name):
            try:
                run_raxml("%s_in.fasta" % full_context, "%s.tree" % full_context, "%s.subsampling_classifications.txt" % full_context, insertion_method, "%s-PARAMS" % new_name, "GTRGAMMA", "%s" % full_context)
            except:
                pass
        else:
            try:
                subprocess.check_call("raxmlHPC-PTHREADS-SSE3 -T %s -f e -m GTRGAMMA -s %s_pruned.fasta -t %s.tree -n %s-PARAMS --no-bfgs > /dev/null 2>&1" % (processors, new_name, new_name, new_name), shell=True)
                os.system("mv RAxML_binaryModelParameters.%s-PARAMS %s-PARAMS" % (new_name, new_name))
                run_raxml("%s_in.fasta" % full_context, "%s.tree" % full_context, "%s.subsampling_classifications.txt" % full_context, insertion_method, "%s-PARAMS" % (split_fields[0]+split_fields[2]), "GTRGAMMA", "%s" % full_context)
            except:
                pass
    if os.path.isfile("%s.resampling_distances.txt" % full_context):
        pass
    else:
        try:
            calculate_pairwise_tree_dists("%s.tree_including_unknowns_noedges.tree" % full_context, "%s.resampling_distances.txt" % full_context)
            for line in open("%s.resampling_distances.txt" % full_context,"U"):
                resample_fields = line.split()
                myid = re.sub("[:']", "",resample_fields[4])
                fixedid = myid.replace("QUERY___","")
                newid = re.sub("[:']","",resample_fields[2])
                fixedid2 = newid.replace("QUERY___","")
                if resample_fields[2] == "'Reference'" and fixedid in name_fixed:
                    outfile.write("resampled distance between Reference and %s = %s\n" % (fixedid, resample_fields[5]))
                elif resample_fields[4] == "'Reference':" and fixedid2 in name_fixed:
                    outfile.write("resampled distance between Reference and %s = %s\n" % (fixedid2, resample_fields[5]))
                else:
                    pass
        except:
            pass
    outfile.close()

def qc_files(fasta,tree):
    #first grab the IDs from the FASTA file:
    fasta_ids = []
    with open(fasta) as my_fasta:
        for record in SeqIO.parse(fasta,"fasta"):
            fasta_ids.append(record.id)
    #Now parse the names from the tree:
    tree_ids = []
    mytree = Phylo.read(tree,'newick')
    for clade in mytree.find_clades():
        if clade.name:
            tree_ids.append(clade.name)
    tree_set = set(tree_ids)
    fasta_set = set(fasta_ids)
    if len(fasta_set) == len(tree_set):
        print("new samples didn't get added correctly...exiting")
        sys.exit()
    else:
        pass

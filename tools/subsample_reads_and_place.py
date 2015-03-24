#!/usr/bin/env python

"""From a SNP matrix, sub-sample
SNPs from a given genome, then find
the SNP level where the unknown genome
is placed correctly >90% of the time"""

from optparse import OptionParser
import sys
import random
import re
import os
import subprocess
import dendropy
import collections

try:
    from dendropy import treecalc
except:
    print "Dendropy needs to be installed for this script to run"
    sys.exit()
from subprocess import Popen
try:
    from Bio import SeqIO
except:
    print "BioPython needs to be installed for this script to run"
    sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def subsample_snps(matrix, name, start):
    """get a list of all possible positions, depending
    on those positions in the original matrix.  Similar
    to method in the main script"""
    allSNPs = [ ]
    for line in open(matrix, "U"):
        if line.startswith("LocusID"):
            pass
        else:
            fields=line.split()
            allSNPs.append(fields[0])
    kept_snps=random.sample(set(allSNPs), int(start))
    outfile = open("%s.%s.tmp.matrix" % (name,start), "w")
    in_matrix=open(matrix,"U")
    firstLine = in_matrix.readline()
    print >> outfile, firstLine,
    first_fields = firstLine.split()
    last=first_fields.index("#SNPcall")
    mygenome=first_fields.index(name)
    fixed_fields = []
    for x in first_fields[:last]:
        fixed_fields.append(re.sub('[:,]', '', x))
    gindex = [ ]
    for x in first_fields[:last]:
        gindex.append(first_fields.index(x))
    for line in in_matrix:
        matrix_fields=line.split()
        if matrix_fields[0] in kept_snps:
            print >> outfile, line,
        else:
            print >> outfile, "\t".join(matrix_fields[:mygenome])+"\t"+"-"+"\t"+"\t".join(matrix_fields[mygenome+1:])+"\n",
    in_matrix.close()
    outfile.close()
    return last

def matrix_to_fasta(matrix_in, name, last):
    """converts a SNP matrix to fasta format.
    Again, slightly different output compared to tested
    function"""
    reduced = [ ]
    out_fasta = open("%s.fasta" % name, "w")
    for line in open(matrix_in):
        fields = line.split("\t")
        reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])
    out_fasta.close()

def branch_lengths_2_decimals(str_newick_tree):
    """replaces branch lengths in scientific notation with decimals
    Identical to tested function in main script"""
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

def insert_sequence(in_fasta, tree, fixed_name, parameters, processors):
    args = ['raxmlHPC-PTHREADS-SSE3', '-f', 'V',
            '-s', '%s' % in_fasta, '-m', 'GTRGAMMA', '-n', '%s' % fixed_name, '-t',
            '%s' % tree, '-T', '%s' % processors, '-R', '%s' % parameters, '--no-bfgs', '>', '/dev/null 2>&1']
    try:
        vcf_fh = open('raxml.out', 'w')
    except:
        print 'could not open raxml file'
    try:
        log_fh = open('raxml.log', 'w')
    except:
        print 'could not open log file'
    try:
        raxml_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
        raxml_run.wait()
        os.system("sed 's/\[[^]]*\]//g' RAxML_labelledTree.%s > %s.tree_including_unknowns_noedges.tree" % (fixed_name,fixed_name))
        subprocess.check_call("rm RAxML*.%s" % fixed_name, shell=True)
    except:
        print "sequence(s) were not inserted into tree!!!!!"
    
def prune_tree(fixed_name,tree):
    tree_full = dendropy.Tree.get_from_path(tree,schema="newick",preserve_underscores=True)
    tree_full.prune_taxa_with_labels(["%s" % fixed_name])
    final_tree = branch_lengths_2_decimals(tree_full.as_string("newick"))
    tmptree = open("%s.tmpx.tree" % fixed_name, "w")
    print >> tmptree, final_tree
    tmptree.close()
    tmptree2 = open("%s.tmpxz.tree" % fixed_name, "w")
    for line in open("%s.tmpx.tree" % fixed_name, "U"):
        if line.startswith("[&U]"):
            fields = line.split()
            fixed_fields = [ ]
            for x in fields:
                fixed_fields.append(x.replace("'",""))
            print >> tmptree2, fixed_fields[1]
        else:
            pass
    tmptree2.close()
    
def calculate_pairwise_tree_dists(intree, output):
    tree = dendropy.Tree.get_from_path(intree, "newick", preserve_underscores=True)
    outfile = open("%s" % output, "w")
    distances = treecalc.PatristicDistanceMatrix(tree)
    for i, t1 in enumerate(tree.taxon_set):
        for t2 in tree.taxon_set[i+1:]:
            print >> outfile, "Distance between '%s' and '%s': %s" % (t1.label, t2.label, distances(t1, t2))
    outfile.close()

def parse_distances(distance_file,fixed_name):
    true_value = []
    infile = open(distance_file, "U")
    for line in infile:
        if len(fixed_name)>1:
            for name in fixed_name:
                fields = line.split()
                if fields[2] == "'%s'" % name and fields[4] == "'Reference':":
                    true_value.append(fields[5])
                elif fields[2] =="'Reference'" and fields[4] =="'%s':" % name:
                    true_value.append(fields[5])
                else:
                    pass
        else:
            fields = line.split()
            if fields[2] == "'%s'" % ''.join(fixed_name) and fields[4] == "'Reference':":
                true_value.append(fields[5])
            elif fields[2] =="'Reference'" and fields[4] =="'%s':" % ''.join(fixed_name):
                true_value.append(fields[5])
            else:
                pass
            
    infile.close()
    return true_value

def get_name_by_ID(in_fasta, ID, out_fasta):
    infile = open(in_fasta, "U")
    output_handle = open(out_fasta, "w")
    seqrecords=[ ]
    for record in SeqIO.parse(infile, "fasta"):
        if record.id == ID:
            seqrecords.append(record)
    SeqIO.write(seqrecords, output_handle, "fasta") 
    infile.close()
    output_handle.close()

def rename_fasta(in_fasta, name, out_fasta):
    infile = open(in_fasta, "U")
    output_handle = open(out_fasta, "w")
    seqrecords=[ ]
    for record in SeqIO.parse(infile, "fasta"):
        seqrecords.append(record.seq)
    print >> output_handle,">"+name+"\n",
    for seqrecord in seqrecords:
        print >> output_handle,seqrecord
    infile.close()
    output_handle.close()

def get_field_index(matrix_in):
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    return last

def remove_sequence(in_fasta,name,out_fasta):
    seqrecords = []
    infile = open(in_fasta, "U")
    output_handle = open(out_fasta, "w")
    for record in SeqIO.parse(infile, "fasta"):
        if record.id != name:
            seqrecords.append(record)
    SeqIO.write(seqrecords, output_handle, "fasta") 
    infile.close()
    output_handle.close()

def main(matrix,tree,name,start,step,end,processors,iterations,deviation):
    aa = subprocess.call(['which', 'raxmlHPC-PTHREADS-SSE3'])
    if aa == 0:
        pass
    else:
        print "RAxML must be in your path as raxmlHPC-PTHREADS-SSE3"
        sys.exit()
    """get starting information"""
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    matrix_path = os.path.abspath("%s" % matrix)
    tree_path = os.path.abspath("%s" % tree)
    """done with getting starting information"""
    os.system("mkdir %s/%s.tmp" % (start_path,name))
    fixed_name = []
    #Changing directory into temporary one
    os.chdir("%s/%s.tmp" % (start_path,name))
    os.system("sed 's/://g' %s | sed 's/,//g' > REF.matrix" % matrix_path)
    fixed_name.append(re.sub('[:,]', '', name))
    calculate_pairwise_tree_dists(tree_path, "%s.all_snps_patristic_distances.txt" % "".join(fixed_name))
    last=get_field_index(matrix_path)
    matrix_to_fasta("REF.matrix","REF",last)
    remove_sequence("REF.fasta", "".join(fixed_name), "REF_pruned.fasta")
    true_value = parse_distances("%s.all_snps_patristic_distances.txt" % "".join(fixed_name),fixed_name)
    outfile = open("%s.results.out" % ''.join(fixed_name), "w")
    prune_tree(''.join(fixed_name),tree_path)
    print "creating parameters file"
    subprocess.check_call("raxmlHPC-PTHREADS-SSE3 -f e -m GTRGAMMA -s REF_pruned.fasta -t %s.tmpxz.tree -n PARAMS --no-bfgs -T %s > /dev/null 2>&1" % ("".join(fixed_name),processors) , shell=True)
    subprocess.check_call("mv RAxML_binaryModelParameters.PARAMS %s.PARAMS" % "".join(fixed_name), shell=True)
    print "starting loop"
    for i in range(start, end+1, step):
        hits = []
        for j in range(1,iterations+1):
            last=subsample_snps("REF.matrix", "".join(fixed_name), i)
            os.system("sed 's/://g' %s.%s.tmp.matrix | sed 's/,//g' > %s.%s.tmp.fixed.matrix" % ("".join(fixed_name),i,"".join(fixed_name),i))
            matrix_to_fasta("%s.%s.tmp.fixed.matrix" % ("".join(fixed_name),i), "%s.%s" % ("".join(fixed_name),i), last)
            get_name_by_ID("%s.%s.fasta" % ("".join(fixed_name),i), ''.join(fixed_name), "%s.%s.%s.tmp.fasta" % ("".join(fixed_name),i,j))
            tmp_name = ''.join(fixed_name)+str(j)
            rename_fasta("%s.%s.%s.tmp.fasta" % ("".join(fixed_name),i,j), tmp_name,"%s.%s.%s.zzyzz.fasta" % ("".join(fixed_name),i,j))
        os.system("cat %s.*.zzyzz.fasta REF_pruned.fasta > %s.joined.fasta" % ("".join(fixed_name),"".join(fixed_name)))
        os.system("rm %s.*.tmp.fasta %s.*.zzyzz.fasta" % ("".join(fixed_name),"".join(fixed_name)))
        insert_sequence("%s.joined.fasta" % "".join(fixed_name), "%s.tmpxz.tree" % "".join(fixed_name), ''.join(fixed_name), "%s.PARAMS" % "".join(fixed_name), processors)
        calculate_pairwise_tree_dists("%s.tree_including_unknowns_noedges.tree" % "".join(fixed_name),"%s.all_patristic_distances.txt" % "".join(fixed_name))
        os.system("cp %s.tree_including_unknowns_noedges.tree %s.%s.%s.tree" % ("".join(fixed_name),"".join(fixed_name),i,j))
        query_names = []
        for j in range(1,iterations+1):
            query_names.append("QUERY___"+"".join(fixed_name)+str(j))
        subsampled_values = parse_distances("%s.all_patristic_distances.txt" % "".join(fixed_name),query_names)
        for value in subsampled_values:
            if (float(value)/float(''.join(true_value)))<(1+float(deviation)):
                hits.append("1")
            else:
                pass
        print >> outfile, i, len(hits)
        if int(len(hits))>=int(0.95*iterations):
            print "optimal value is for %s is %s" % ("".join(fixed_name),i)
            print >> outfile, "optimial value is for %s is %s" % ("".join(fixed_name),i)
            break
    os.system("mv %s.results.out %s" % (''.join(fixed_name), start_path))
    os.chdir("%s" % start_path)
    os.system("rm -rf %s/%s.tmp" % (start_path,name))
        
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-m", "--snp_matrix", dest="matrix",
                      help="path to NASP snp_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-t", "--tree", dest="tree",
                      help="path to input tree [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-n", "--name", dest="name",
                      help="name of genome to test [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-s", "--start", dest="start",
                      help="starting number of SNPs to sample, defaults to 50",
                      action="store", type="int", default="50")
    parser.add_option("-p", "--step", dest="step",
                      help="step through SNPs at this level, defaults to 100",
                      action="store", type="int", default="100")
    parser.add_option("-e", "--end", dest="end",
                      help="ending number of SNPs to sample, defaults to 100000",
                      action="store", type="int", default="100000")
    parser.add_option("-o", "--processors", dest="processors",
                      help="number of processors to use with RAxML, defaults to 4",
                      action="store", type="int", default="4")
    parser.add_option("-i", "--iterations", dest="iterations",
                      help="number of iterations at each level, defaults to 10",
                      action="store", type="int", default="10")
    parser.add_option("-d", "--deviation", dest="deviation",
                      help="deviation from 1, to determine correct placement, defaults to 0.05",
                      action="store", type="float", default="0.05")
    
    options, args = parser.parse_args()
    
    mandatories = ["matrix", "tree", "name"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix,options.tree,options.name,options.start,options.step,options.end,
         options.processors,options.iterations,options.deviation)

#!/usr/bin/env python

"""From a SNP matrix, sub-sample
SNPs from a given genome, then find
the SNP level where the unknown genome
is placed correctly 95% of the time"""

from optparse import OptionParser
import sys
import random
import re
import os
import subprocess
import dendropy
from dendropy import treecalc
from subprocess import Popen
from Bio import SeqIO

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def subsample_snps(matrix, name, start):
    """get a list of all possible positions, depending
    on those positions in the original matrix"""
    allSNPs = [ ]
    for line in open(matrix, "U"):
        if line.startswith("LocusID"):
            pass
        else:
            fields=line.split()
            allSNPs.append(fields[0])
    kept_snps=random.sample(set(allSNPs), int(start))
    outfile = open("%s.tmp.matrix" % start, "w")
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
    """converts an ISG matrix to fasta format"""
    reduced = [ ]
    out_fasta = open("%s.fasta" % name, "w")
    for line in open(matrix_in):
        fields = line.split("\t")
        reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])

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

def insert_sequence(in_fasta, tree, processors):
    args = ['raxmlHPC-PTHREADS', '-T', '%s' % processors, '-f', 'v',
	     '-s', '%s' % in_fasta, '-m', 'GTRGAMMA', '-n', 'out', '-t',
	     '%s' % tree, '>', '/dev/null 2>&1']
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
    except:
        print "sequence(s) were not inserted into tree!!!!!"
    os.system("sed 's/\[[^]]*\]//g' RAxML_labelledTree.out > tree_including_unknowns_noedges.tree")
    subprocess.check_call("rm RAxML_*", shell=True)
    
def prune_tree(fixed_name,tree):
    tree_full = dendropy.Tree.get_from_path(tree,schema="newick",preserve_underscores=True)
    tree_full.prune_taxa_with_labels(fixed_name)
    final_tree = branch_lengths_2_decimals(tree_full.as_string("newick"))
    tmptree = open("tmpx.tree", "w")
    print >> tmptree, final_tree
    tmptree.close()
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
        for name in fixed_name:
            fields = line.split()
            if fields[2] == "'%s'" % name and fields[4] == "'Reference':":
                true_value.append(fields[5])
            elif fields[2] =="'Reference'" and fields[4] =="'%s':" % name:
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
    matrix=open(matrix_in, "rU")
    firstLine = open(matrix_in).readline()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    return last

def main(matrix,tree,name,start,step,end,processors,iterations):
    fixed_name = []
    os.system("sed 's/://g' %s | sed 's/,//g' > REF.matrix" % matrix)
    fixed_name.append(re.sub('[:,]', '', name))
    calculate_pairwise_tree_dists(tree, "all_snps_patristic_distances.txt")
    last=get_field_index(matrix)
    matrix_to_fasta("REF.matrix","REF",last)
    true_value = parse_distances("all_snps_patristic_distances.txt",fixed_name)
    for i in range(start, end+1, step):
        hits = []
        for j in range(1,iterations):
            last=subsample_snps(matrix, name, i)
            os.system("sed 's/://g' %s.tmp.matrix | sed 's/,//g' > %s.tmp.fixed.matrix" % (i,i))
            matrix_to_fasta("%s.tmp.fixed.matrix" % i, i, last)
            get_name_by_ID("%s.fasta" % i, ''.join(fixed_name), "%s.%s.tmp.fasta" % (i,j))
            tmp_name = ''.join(fixed_name)+str(j)
            rename_fasta("%s.%s.tmp.fasta" % (i,j), tmp_name,"%s.%s.zzyzz.fasta" % (i,j))
        prune_tree(fixed_name,tree)
        os.system("cat *.zzyzz.fasta REF.fasta > input.fasta")
        os.system("rm *.tmp.fasta *.zzyzz.fasta")
        insert_sequence("input.fasta", "tmpxz.tree", processors)
        calculate_pairwise_tree_dists("tree_including_unknowns_noedges.tree","all_patristic_distances.txt")
        query_names = []
        for j in range(1,iterations):
            query_names.append("QUERY___"+"".join(fixed_name)+str(j))
        subsampled_values = parse_distances("all_patristic_distances.txt",query_names)
        for value in subsampled_values:
            if float(value)/float(''.join(true_value))<1.02 and float(value)/float(''.join(true_value))>0.98:
                hits.append("1")
        print i, len(hits)
        if int(len(hits))>=95:
            print "optimial value is %s" % i
            break
        
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
                      help="number of iterations at each level",
                      action="store", type="int", default="10")
    
    options, args = parser.parse_args()
    
    mandatories = ["matrix", "tree", "name"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix,options.tree,options.name,options.start,options.step,options.end,options.processors,options.iterations)

    

#!/usr/bin/env python

"""from a NASP matrix, creates the tree
necessary for wg-fast to run"""

from optparse import OptionParser
import subprocess
import os
import sys

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def get_field_index(matrix_in):
    matrix=open(matrix_in, "rU")
    firstLine = open(matrix_in).readline()
    fixed_line = firstLine.strip()
    first_fields = fixed_line.split("\t")
    last=first_fields.index("#SNPcall")
    return last

def matrix_to_fasta(matrix_in, last):
    """converts an ISG matrix to fasta format"""
    reduced = [ ]
    out_fasta = open("all.fasta", "w")
    for line in open(matrix_in):
        fields = line.split("\t")
        reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])

def main(matrix):
    """determines whether or not raxml is in your path"""
    ab = subprocess.call(['which', 'raxmlHPC-SSE3'])
    if ab == 0:
        pass
    else:
        print "RAxML must be in your path as raxmlHPC-PTHREADS"
        sys.exit()
    last=get_field_index(matrix)
    matrix_to_fasta(matrix, last)
    os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
    subprocess.check_call("raxmlHPC-SSE3 -f d -p 12345 -m ASC_GTRGAMMA -s out.fasta -n nasp > /dev/null 2>&1", shell=True)
    subprocess.check_call("raxmlHPC-SSE3 -f e -m ASC_GTRGAMMA -s out.fasta -t RAxML_bestTree.nasp -n PARAMS > /dev/null 2>&1", shell=True)
    subprocess.check_call("mv RAxML_bestTree.nasp nasp_raxml.tree", shell=True)
    subprocess.check_call("mv RAxML_binaryModelParameters.PARAMS nasp.PARAMS", shell=True)
    subprocess.check_call("rm RAxML_* out.fasta all.fasta", shell=True)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--snp_matrix", dest="matrix",
                      help="path to NASP snp_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    options, args = parser.parse_args()
    
    mandatories = ["matrix"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix)

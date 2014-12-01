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

def test_models(option, opt_str, value, parser):
    if "GTRGAMMA" in value:
        setattr(parser.values, option.dest, value)
    elif "ASC_GTRGAMMA" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "substitution model is not supported"
        sys.exit()

def get_field_index(matrix_in):
    """untested function"""
    matrix=open(matrix_in, "rU")
    firstLine = open(matrix_in).readline()
    fixed_line = firstLine.strip()
    first_fields = fixed_line.split("\t")
    last=first_fields.index("#SNPcall")
    matrix.close()
    return last

def matrix_to_fasta(matrix_in, last):
    """converts a NASP matrix to fasta format.
    Similar to tested function in main script,
    but slightly different output"""
    reduced = [ ]
    out_fasta = open("all.fasta", "w")
    for line in open(matrix_in, "U"):
        fields = line.split("\t")
        reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])
    out_fasta.close()

def main(matrix, model):
    """determines whether or not raxml is in your path"""
    ab = subprocess.call(['which', 'raxmlHPC-SSE3'])
    if ab == 0:
        pass
    else:
        print "RAxML must be in your path as raxmlHPC-SSE3"
        sys.exit()
    last=get_field_index(matrix)
    matrix_to_fasta(matrix, last)
    os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
    if model == "ASC_GTRGAMMA":
        subprocess.check_call("raxmlHPC-SSE3 -f d -p 12345 -m %s -s out.fasta -n nasp --asc-corr=lewis --no-bfgs > /dev/null 2>&1" % model, shell=True)
        subprocess.check_call("raxmlHPC-SSE3 -f e -m %s -s out.fasta -t RAxML_bestTree.nasp -n PARAMS --asc-corr=lewis --no-bfgs > /dev/null 2>&1" % model, shell=True)
    else:
        subprocess.check_call("raxmlHPC-SSE3 -f d -p 12345 -m %s -s out.fasta -n nasp --no-bfgs > /dev/null 2>&1" % model, shell=True)
        subprocess.check_call("raxmlHPC-SSE3 -f e -m %s -s out.fasta -t RAxML_bestTree.nasp -n PARAMS --no-bfgs > /dev/null 2>&1" % model, shell=True)
    subprocess.check_call("mv RAxML_bestTree.nasp nasp_raxml.tree", shell=True)
    subprocess.check_call("mv RAxML_binaryModelParameters.PARAMS nasp.PARAMS", shell=True)
    subprocess.check_call("rm RAxML_* out.fasta all.fasta", shell=True)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--snp_matrix", dest="matrix",
                      help="path to NASP snp_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-o", "--model", dest="model",
                      help="model for RAxML, defaults to ASC_GTRGAMMA, can also be GTRGAMMA",
                      action="callback", callback=test_models, type="string", default="ASC_GTRGAMMA")
    options, args = parser.parse_args()
    
    mandatories = ["matrix"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix,options.model)

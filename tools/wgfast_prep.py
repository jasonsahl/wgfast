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
        print('%s file cannot be opened' % option)
        sys.exit()

def test_models(option, opt_str, value, parser):
    if "GTRGAMMA" in value:
        setattr(parser.values, option.dest, value)
    elif "ASC_GTRGAMMA" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("substitution model is not supported")
        sys.exit()

def get_field_index(matrix_in):
    """untested function"""
    with open(matrix_in) as my_matrix:
        firstLine = open(matrix_in).readline().rstrip()
        first_fields = firstLine.split("\t")
        last=first_fields.index("#SNPcall")
    return last

def matrix_to_fasta(matrix_in, last):
    """converts a NASP matrix to fasta format.
    Similar to tested function in main script,
    but slightly different output"""
    reduced = [ ]
    out_fasta = open("all.fasta", "w")
    with open(matrix_in) as my_matrix:
        for line in my_matrix:
            fields = line.split("\t")
            reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        out_fasta.write(">"+str(x[0])+"\n")
        out_fasta.write("".join(x[1:])+"\n")
    out_fasta.close()

def main(matrix,model,processors,algorithm):
    """determines whether or not raxml is in your path"""
    if algorithm == "raxml-ng":
        ab = subprocess.call(['which', 'raxml-ng'])
        if ab == 0:
            pass
        else:
            print("RAxML must be in your path as raxml-ng")
            sys.exit()
    elif algorithm == "raxml-HPC":
        ab = subprocess.call(['which', 'raxmlHPC-PTHREADS-SSE3'])
        if ab == 0:
            pass
        else:
            print("RAxML must be in your path as raxmlHPC-PTHREADS-SSE3")
            sys.exit()
    last=get_field_index(matrix)
    matrix_to_fasta(matrix, last)
    #Prep the creation of the FASTA file, removing odd characters
    os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
    if model == "ASC_GTRGAMMA":
        subprocess.check_call("raxmlHPC-SSE3 -f d -p 12345 -m %s -s out.fasta -n nasp --asc-corr=lewis --no-bfgs > /dev/null 2>&1" % model, stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        subprocess.check_call("raxmlHPC-SSE3 -f e -m %s -s out.fasta -t RAxML_bestTree.nasp -n PARAMS --asc-corr=lewis --no-bfgs > /dev/null 2>&1" % model, stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    else:
        if algorithm == "raxml-HPC":
            subprocess.check_call("raxmlHPC-PTHREADS-SSE3 -T %s -f d -p 12345 -m %s -s out.fasta -n nasp --no-bfgs > /dev/null 2>&1" % (processors,model), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
            subprocess.check_call("raxmlHPC-PTHREADS-SSE3 -T %s -f e -m %s -s out.fasta -t RAxML_bestTree.nasp -n PARAMS --no-bfgs > /dev/null 2>&1" % (processors,model), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        elif algorithm == "raxml-ng":
            subprocess.check_call("raxml-ng --msa out.fasta --model GTR+G --threads %s --prefix nasp" % processors,stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    if algorithm == "raxml-HPC":
        subprocess.check_call("mv RAxML_bestTree.nasp nasp_raxml.tree", shell=True)
        subprocess.check_call("mv RAxML_binaryModelParameters.PARAMS nasp.PARAMS", shell=True)
        subprocess.check_call("rm RAxML_* out.fasta all.fasta", shell=True)
    else:
        subprocess.check_call("mv nasp.raxml.bestTree nasp_raxml.tree", stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        subprocess.check_call("rm nasp.raxml.startTree out.fasta all.fasta", stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    print("Model used: %s" % model)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--snp_matrix", dest="matrix",
                      help="path to NASP snp_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-o", "--model", dest="model",
                      help="model for RAxML, can be ASC_GTRGAMMA or GTRGAMMA [DEFAULT]",
                      action="callback", callback=test_models, type="string", default="GTRGAMMA")
    parser.add_option("-p", "--processors", dest="processors",
                      help="number of processors to use with GTRGAMMA, defaults to 2",
                      action="store", type="int", default="2")
    parser.add_option("-a", "--algorithm", dest="algorithm",
                      help="algorithm to use, either raxml-ng or raxml-HPC",
                      action="store", type="string", default="raxml-ng")
    options, args = parser.parse_args()

    mandatories = ["matrix"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.matrix,options.model,options.processors,options.algorithm)

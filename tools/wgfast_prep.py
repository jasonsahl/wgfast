#!/usr/bin/env python

"""from a NASP matrix, creates the tree
necessary for wg-fast to run"""

from optparse import OptionParser
import subprocess
from wg_fast.util import matrix_to_fasta
import os
import sys

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(matrix):
    """determines whether or not raxml is in your path"""
    ab = subprocess.call(['which', 'raxmlHPC-PTHREADS'])
    if ab == 0:
        pass
    else:
        print "RAxML must be in your path as raxmlHPC-PTHREADS"
        sys.exit()
    matrix_to_fasta(matrix)
    os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
    subprocess.check_call("raxmlHPC-PTHREADS -f d -p 12345 -m GTRGAMMA -s out.fasta -n nasp -T 4 > /dev/null 2>&1", shell=True)
    subprocess.check_call("mv RAxML_bestTree.nasp nasp_raxml.tree", shell=True)
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

#!/usr/bin/env python

"""From a NASP matrix, subsample
SNPs at a given level, calculates
distance matrices, and performs a Mantel
test with the Pearson correlation"""

import optparse
import sys
import subprocess
import random
import os
import glob
from wg_fast.util import test_file

def subsample_snps(matrix, snps, iterations):
    """get a list of all possible positions, depending
    on those positions in the original matrix - tested?"""
    allSNPs = [ ]
    for line in open(matrix, "U"):
        if line.startswith("LocusID"):
            pass
        else:
            fields=line.split()
            allSNPs.append(fields[0])
    for x in range(1,int(iterations)+1):
        kept_snps=random.sample(set(allSNPs), snps)
        outfile = open("%s.%s.tmp.matrix" % (snps, x), "w")
        in_matrix=open(matrix,"U")
        firstLine = in_matrix.readline()
        print >> outfile, firstLine,
        first_fields = firstLine.split()
        for line in in_matrix:
             matrix_fields=line.split()
             if matrix_fields[0] in kept_snps:
                  print >> outfile, line,
    outfile.close()
                    
def matrix_to_fasta(matrix, prefix):
    """converts a NASP matrix to fasta format.  Includes
    SNP call field, which is different than the function
    included for the main WG-FAST pipeline"""
    reduced = [ ]
    out_fasta = open("%s" % prefix, "w")
    firstLine = open(matrix).readline()
    first_fields = firstLine.split()
    last=first_fields.index("#SNPcall")
    for line in open(matrix, "U"):
        fields = line.split()
        reduced.append(fields[1:last])
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_fasta, ">"+str(x[0])
        print >> out_fasta, "".join(x[1:])
    out_fasta.close()

def compare_matrices(start_dir, ref_matrix, processors):
    for infile in glob.glob(os.path.join(start_dir, "*.tmp.matrix.xyzzy")):
        subprocess.check_call(['mothur',
                               '#dist.seqs(fasta=%s, output=lt, processors=%s)' % (infile,processors),'>/dev/null 2>&1'])
    for distfile in glob.glob(os.path.join(start_dir, "*.tmp.matrix.phylip.dist")):
        subprocess.check_call(['mothur',
        '#mantel(phylip1=%s, phylip2=%s, method=pearson)' % (distfile, ref_matrix), '>/dev/null 2>&1'])

def process_results(start_dir):
    outfile = open("results.txt", "w")
    for infile in glob.glob(os.path.join(start_dir, '*.tmp.matrix.phylip.mantel')):
        for line in open(infile, "U"):
            if line.startswith("Mantel"):
                pass
            else:
                fields = line.split()
                print >> outfile, fields[0]
    outfile.close()
            
def main(matrix,snps,iterations,processors):
    start_dir = os.getcwd()
    ac = subprocess.call(['which', 'mothur'])
    if ac == 0:
        pass
    else:
        print "mothur must be in your path"
        sys.exit()
    print "citation: Schloss PD, Westcott SL, Ryabin T, Hall JR, Hartmann M, Hollister EB, Lesniewski RA, Oakley BB, Parks DH, Robinson CJ, Sahl JW, Stres B, Thallinger GG, Van Horn DJ, Weber CF. Introducing mothur: Open-Source, Platform-Independent, Community-Supported Software for Describing and Comparing Microbial Communities. Appl Environ Microbiol. 2009;75(23):7537-41"
    matrix_to_fasta(matrix, "reference_matrix.fasta")
    subprocess.check_call(['mothur',
                           '#dist.seqs(fasta=reference_matrix.fasta, output=lt, processors=%s)' % processors,'>','/dev/null 2>&1'])
    ref_matrix = "reference_matrix.phylip.dist"
    subsample_snps(matrix, snps, iterations)
    for infile in glob.glob(os.path.join(start_dir, '*.tmp.matrix')):
        matrix_to_fasta(infile, "%s.xyzzy" % infile)
    compare_matrices(start_dir, ref_matrix, processors)
    process_results(start_dir)
    os.system("rm *.tmp.matrix* mothur.* reference_matrix*")
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-m", "--matrix", dest="matrix",
                      help="/path/to/NASP matrix [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-s", "--snps", dest="snps",
                      help="number of SNPs to subsample [REQUIRED]",
                      type="int", action="store")
    parser.add_option("-i", "--iterations", dest="iterations",
                      help="number of iterations to perform, defaults to 100",
                      type="int", default="100",action="store")
    parser.add_option("-p", "--processors", dest="processors",
                      help="number of processors to use, defaults to 4",
                      type="int", default="4",action="store")
    options, args = parser.parse_args()
    
    mandatories = ["matrix", "snps"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix,options.snps,options.iterations,options.processors)
    

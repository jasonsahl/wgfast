#!/usr/bin/env python
"""Tool for WGFAST database curation"""

import sys
import os
import optparse
import subprocess
import glob

def mp_shell(func, params, numProc):
    from multiprocessing import Pool
    p = Pool(numProc)
    out = p.map(func, params)
    p.terminate()
    return out

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print("%s directory cannot be found" % option)
        sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def get_field_index(matrix_in):
    firstLine = open(matrix_in).readline().rstrip()
    first_fields = firstLine.split("\t")
    last=first_fields.index("#SNPcall")
    return last

def process_fasta_dir(fasta_dir):
    file_list = []
    for infile in glob.glob(os.path.join(fasta_dir, '*.fasta')):
        file_list.append(infile)
    return file_list

def test_out_dir_contents(file_list,out_dir):
    genome_test = []
    if os.path.exists("%s/bestsnp.tsv" % out_dir):
        nasp_headers = []
        genome_names = []
        with open("%s/bestsnp.tsv" % out_dir) as f:
            for line in f.readlines()[:1]:
                fields = line.split()
                last=fields.index("#SNPcall")
                for field in fields[1:last]:
                    nasp_headers.append(field)
        for file in file_list:
            name_path = file.split("/")
            genome_names.append(name_path[-1].replace(".fasta",""))
        for name in genome_names:
            if name in nasp_headers:
                pass
            else:
                genome_test.append("1")
    else:
        print("Can't find bestsnp.tsv in your wgfast directory")
        sys.exit()
    return len(genome_test)

def _nasp_externals_wf(data):
    tn = data[0]
    f = data[1]
    ref = data[2]
    f_base = os.path.basename(f)
    f_redux = f_base.replace(".fasta","")
    if os.path.exists("%s.frankenfasta" % f_redux):
        pass
    else:
        subprocess.check_call("convert_external_genome --reference %s --external %s" % (ref,f), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)

def make_frankenfastas(fasta_path,ref_path,processors):
    files = []
    files_and_temp_names = []
    for infile in glob.glob(os.path.join(fasta_path, '*.fasta')):
        files.append(infile)
    for idx,f in enumerate(files):
        files_and_temp_names.append([str(idx),os.path.join(fasta_path, f),ref_path])
    mp_shell(_nasp_externals_wf, files_and_temp_names, processors)

def matrix_to_fasta(matrix_in,last):
    """converts a NASP matrix to fasta format.
    Similar to tested function in main script,
    but slightly different output"""
    reduced = []
    out_fasta = open("all.fasta", "w")
    with open(matrix_in) as my_file:
        for line in my_file:
            fields = line.split("\t")
            reduced.append(fields[1:last])
    test=map(list,zip(*reduced))
    for x in test:
        out_fasta.write(">"+str(x[0])+"\n"+"".join(x[1:])+"\n")
    out_fasta.close()

def create_tree_files(matrix):
    last=get_field_index(matrix)
    matrix_to_fasta(matrix, last)
    os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
    subprocess.check_call("raxmlHPC-SSE3 -f d -p 12345 -m ASC_GTRGAMMA -s out.fasta -n nasp --asc-corr=lewis --no-bfgs", stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    subprocess.check_call("raxmlHPC-SSE3 -f e -m ASC_GTRGAMMA -s out.fasta -t RAxML_bestTree.nasp -n PARAMS --asc-corr=lewis --no-bfgs", stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    subprocess.check_call("mv RAxML_bestTree.nasp nasp_raxml.tree", shell=True)
    subprocess.check_call("mv RAxML_binaryModelParameters.PARAMS nasp.PARAMS", shell=True)
    subprocess.check_call("rm RAxML_* out.fasta all.fasta", shell=True)

def main(fasta_dir, reference, wgfast_dir, processors):
    dependencies = ['nasp','find_duplicates','raxmlHPC-SSE3','convert_external_genome']
    for dependency in dependencies:
        ab = subprocess.call(['which', '%s' % dependency])
        if ab == 0:
            pass
        else:
            print("%s isn't in your path, but needs to be" % dependency)
            sys.exit()
    start_dir = os.getcwd()
    #Create a temporary directory where the work will be done
    if os.path.exists("%s/tmp_workdir" % start_dir):
        os.system("rm -rf %s/tmp_workdir" % start_dir)
        os.system("mkdir %s/tmp_workdir" % start_dir)
    else:
        os.system("mkdir %s/tmp_workdir" % start_dir)
    ap_start = os.path.abspath(start_dir)
    fasta_path=os.path.abspath(fasta_dir)
    ref_path=os.path.abspath(reference)
    out_path=os.path.abspath(wgfast_dir)
    ref_name=os.path.basename(ref_path)
    ref_base = ref_name.replace(".fasta","")
    #Check if output directory exists
    #First, get a list of all FASTA files that need to be analyzed
    file_list = process_fasta_dir(fasta_path)
    nasp_value = test_out_dir_contents(file_list,out_path)
    if nasp_value > 0:
        print("%s genome(s) not in your NASP matrix: running NASP..." % nasp_value)
        os.system("mv %s/bestsnp.tsv %s/bestsnp.tsv.bk" % (out_path,out_path))
        #Step 1 is to identify duplicates...requires script in path, outfile is reference.delta
        #change directory to the work directory
        os.chdir("%s/tmp_workdir" % start_dir)
        subprocess.check_call("find_duplicates --reference %s" % ref_path, stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        make_frankenfastas(fasta_path,ref_path,processors)
        subprocess.check_call("nasptool matrix --num-threads %s --reference-fasta %s --reference-dups %s/tmp_workdir/duplicates.txt %s/tmp_workdir/*.frankenfasta" % (processors,ref_path,ap_start,ap_start),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        print("Moving new nasp matrix into place")
        os.system("mv %s/tmp_workdir/bestsnp.tsv %s/bestsnp.tsv" % (start_dir,out_path))
        print("Creating tree files")
        create_tree_files("%s/bestsnp.tsv" % out_path)
        print("moving tree files into place")
        os.system("mv %s/tmp_workdir/nasp_raxml.tree %s/tmp_workdir/nasp.PARAMS %s" % (start_dir,start_dir,out_path))
        #Clean up
        nasp_value = test_out_dir_contents(file_list,out_path)
        if nasp_value > 0:
            print("files were not updated..check input and try again")
        else:
            print("Your WG-FAST directory has been updated")
        os.chdir(start_dir)
        os.system("rm -rf %s/tmp_workdir" % start_dir)
    else:
        print("Your WG-FAST directory is up to date")

if __name__ == "__main__":
    try:
        main(snakemake.params[0],
        snakemake.params[1],
        snakemake.params[2],
        snakemake.threads)
    except NameError:
        usage="usage: %prog [options]"
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-d", "--fasta_dir", dest="fasta_dir",
                        help="/path/to/FASTA dir [REQUIRED]",
                        type="string", action="callback", callback=test_dir)
        parser.add_option("-r", "--reference", dest="reference",
                        help="reference FASTA to use [REQUIRED]",
                        type="string", action="callback", callback=test_file)
        parser.add_option("-w", "--wgfast_dir", dest="wgfast_dir",
                        help="wgfast reference directory [REQUIRED]",
                        type="string", action="callback", callback=test_dir)
        parser.add_option("-p", "--processors", dest="processors",
                        help="number of processors, defaults to 2",
                        type="int", action="store", default=2)
        options, args = parser.parse_args()

        mandatories = ["fasta_dir","reference","wgfast_dir"]

        for m in mandatories:
            if not getattr(options, m, None):
                print("\nMust provide %s.\n" %m)
                parser.print_help()
                exit(-1)

        main(
            options.fasta_dir,
            options.reference,
            options.wgfast_dir,
            options.processors)
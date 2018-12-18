#!/usr/bin/env python
"""Tool for WGFAST database curation"""

import sys
import os
import argparse
import subprocess
import glob
import shutil
from multiprocessing import Pool

from tempfile import TemporaryDirectory
from wgfastdb.parse_util import types, config_logging, fasta


def mp_shell(func, params, numProc):
    p = Pool(numProc)
    out = p.map(func, params)
    p.terminate()
    return out

def get_field_index(matrix_in):
    firstLine = open(matrix_in).readline().rstrip()
    first_fields = firstLine.split("\t")
    last = first_fields.index("#SNPcall")
    return last

def process_fasta_dir(fasta_dir):
    file_list = []
    for infile in glob.glob(os.path.join(fasta_dir, '*.fasta')):
        file_list.append(infile)
    return file_list

def test_out_dir_contents(file_list, out_dir):
    run_nasp = False
    matrix_file = "{}/bestsnp.tsv".format(out_dir)
    if os.path.exists(matrix_file):
        nasp_headers = []
        genome_names = []
        with open(matrix_file) as f:
            for line in f.readlines()[:1]:
                fields = line.split()
                last = fields.index("#SNPcall")
                for field in fields[1:last]:
                    nasp_headers.append(field)
        for file in file_list:
            name_path = file.split("/")
            genome_names.append(name_path[-1].replace(".fasta", ""))
        for name in genome_names:
            if name not in nasp_headers:
                print(
                    "{} genome not in your NASP matrix: running NASP...".format(
                        name
                    ))
                run_nasp = True
        if run_nasp:
            logger.info("Moving existing bestsnp.tsv to bestsnp.tsv.bk")
            os.rename(
                "{}/bestsnp.tsv".format(out_dir),
                "{}/bestsnp.tsv.bk".format(out_dir))
    else:
        run_nasp = True
        logger.info("No bestsnp.tsv in your wgfast directory: running NASP...")
    return run_nasp

def _nasp_externals_wf(data):
    tn = data[0]
    f = data[1]
    ref = data[2]
    f_base = os.path.basename(f)
    f_redux = f_base.replace(".fasta","")
    if os.path.exists("%s.frankenfasta" % f_redux):
        pass
    else:
        subprocess.check_call(
            "convert_external_genome --reference %s --external %s" % (ref, f),
            shell=True)

def make_frankenfastas(fasta_path, ref_path, processors):
    files = []
    files_and_temp_names = []
    for infile in glob.glob(os.path.join(fasta_path, '*.fasta')):
        files.append(infile)
    for idx, f in enumerate(files):
        files_and_temp_names.append(
            [str(idx), os.path.join(fasta_path, f), ref_path])
    mp_shell(_nasp_externals_wf, files_and_temp_names, processors)

def matrix_to_fasta(matrix_in, last):
    """converts a NASP matrix to fasta format.
    Similar to tested function in main script,
    but slightly different output"""
    reduced = []
    with open(matrix_in) as my_file:
        for line in my_file:
            fields = line.split("\t")
            reduced.append(fields[1:last])
    test = map(list, zip(*reduced))
    with open("all.fasta", "w") as out_fasta:
        for x in test:
            out_fasta.write(
                fasta(str(x[0]), "".join(x[1:])))

def remove_files(*args):
    for arg in args:
        os.remove(arg)

def create_tree_files(matrix, wrkdir):
    logger.info("Generating tree files")
    last = get_field_index(matrix)
    matrix_to_fasta(matrix, last)
    os.system("sed 's/://g' all.fasta | sed 's/,//g' > out.fasta")
    subprocess.check_call(
        "raxmlHPC-SSE3 -f d -p 12345 -m ASC_GTRGAMMA -s out.fasta "
        "-n nasp --asc-corr=lewis --no-bfgs -w {}".format(wrkdir),
        shell=True)
    subprocess.check_call(
        "raxmlHPC-SSE3 -f e -m ASC_GTRGAMMA -s out.fasta "
        "-t RAxML_bestTree.nasp -n PARAMS --asc-corr=lewis --no-bfgs "
        "-w {}".format(wrkdir),
        shell=True)
    os.rename("RAxML_bestTree.nasp", "nasp_raxml.tree")
    os.rename(
        "RAxML_binaryModelParameters.PARAMS",
        "nasp.PARAMS"
    )
    remove_files("out.fasta", "all.fasta", *glob.glob("RAxML_*"))
    logger.info("Finished generating tree files")

def nasp(ref_path, fasta_path, tmp_dir, processors):
    # check for duplicates
    logger.info("Finding duplicates in reference: {}".format(ref_path))
    subprocess.check_call(
            "find_duplicates --reference {}".format(ref_path), shell=True)
    logger.info("Generating frankenfastas")
    make_frankenfastas(fasta_path, ref_path, processors)
    logger.info("Building nasp matrices")
    subprocess.check_call(
            "nasp matrix --num-threads {0} --reference-fasta {1} "
            "--reference-dups {2}/duplicates.txt "
            "{2}/*.frankenfasta".format(
                processors, ref_path, tmp_dir),
                shell=True)
    logger.info("Finished running NASP")

def check_dependencies():
    dependencies = [
        'nasp','find_duplicates','raxmlHPC-SSE3','convert_external_genome']
    for dependency in dependencies:
        if shutil.which(dependency) is None:
            problem = "{} is not in your path, but needs to be".format(
                dependency)
            logger.error(problem)
            raise RuntimeError(problem)



def main(fasta_path, ref_path, out_path, processors, log):
    logger = config_logging(log, level="INFO")
    check_dependencies()
    # First, get a list of all FASTA files that need to be analyzed
    file_list = process_fasta_dir(fasta_path)
    run_nasp = test_out_dir_contents(file_list, out_path)
    with TemporaryDirectory() as temp_dir:   
        os.chdir(temp_dir)
        print("Tempdir", temp_dir)

        if run_nasp:
        # Step 1 is to identify duplicates...requires script in path,
        # outfile is reference.delta

            nasp(ref_path, fasta_path, temp_dir, processors)
            try:
                shutil.move(
                    "{}/bestsnp.tsv".format(temp_dir),
                    "{}/bestsnp.tsv".format(out_path))

                logger.info("Moving new nasp matrix into place")
            except FileNotFoundError:
                logger.error(
                    "bestsnp.tsv was not generated .. check input and try again")
                raise

        create_tree_files("{}/bestsnp.tsv".format(out_path), temp_dir)

        try:
            shutil.move(
                "{}/nasp_raxml.tree".format(temp_dir),
                out_path)
            shutil.move(
                "{}/nasp.PARAMS".format(temp_dir),
                out_path
            )
            logger.info("moving tree files into place")
        except FileNotFoundError:
            logger.error(
                "Tree files were not generated .. check input and try again")
        
        logger.info("Moving reference file")
        shutil.move(ref_path, os.path.join(out_path, "reference.fasta"))
        logger.info("Your WG-FAST directory has been updated")

def parse(argv=None):
    if argv is None:
        argv = sys.argv
    print(argv)

    parser = argparse.ArgumentParser(
        prog="wgfast_tree",
        description="Generate tree and matrix files for WG-FAST",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-d", "--fasta_path", type=types['path_type'], required=True,
        help="/path/to/FASTA directory")
    parser.add_argument(
        "-r", "--ref_path", type=types['file_type'], required=True,
        help="reference FASTA to use.")
    parser.add_argument(
        "-w", "--out_path", type=types['outpath_type'] , required=True,
        help="wgfast reference directory")
    parser.add_argument(
        "-p", "--processors", type=types['positive_int'], default=2,
        help="number of processors")
    parser.add_argument(
        "-l", "--log", type=types['outfile_type'], default=None,
        help="Log file path")
    
    args = parser.parse_args(argv)
    return vars(args)

if __name__ == "__main__":
    print(__file__)
    try:
        argv = list(map(str, [
            "--fasta_path", snakemake.input[0],
            "--ref_path", snakemake.input[1],
            "--out_path", snakemake.params[0],
            "--processors", snakemake.threads,
            "--log", snakemake.log[0]]))

    except NameError:
        argv = sys.argv
    
    finally:
        ARGS = parse(argv)
        main(**ARGS)
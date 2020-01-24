import argparse
import sys
import datetime
import os
import json
import logging
import click
import subprocess as sp
import pandas as pd
import numpy as np
from collections import namedtuple

from wgfastdb.parse_util import types, error, warn, config_logging


DEFAULT_LOG_FNAME = "wgfastdb.log"
LOG = logging.getLogger(__name__)

def parse():
    parser = argparse.ArgumentParser(
        prog="wgfastdb",
        description="Database setup for WGFAST",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    download = parser.add_argument_group(
        'Download', 'Download and update sequences from NCBI'
    )
    
    download.add_argument(
        "path", type=types['outpath_type'], default="./",
        metavar="PATH",
        help="""
        Path to existing fasta database files or 
        location where the files should be downloaded.
        """
    )

    download.add_argument(
        "--no_update", action='store_true', default=False,
        help="""
        Do not sync your collection with the latest assembly versions
        """
    )
    download.add_argument(
        "--no_assembly_update", action="store_true", default=False,
        help="""
        Do not download the latest assembly summary and taxonomy
        dump and use your local copies.
        """
    )

    download.add_argument(
        "--download_only", action="store_true", default=False,
        help="""
        Run only the download step.
        """
    )

    group = download.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--species", '-s', nargs='+', type=str,
        help="""
        List of species to build database. The species name
        must match exactly a species directory at 
        ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/
        """
    )
    group.add_argument(
        "--config", '-cfg', type=types['file_type'],
        help="""
        Path to config table (.csv). 
        Table MUST include a "species" column and may 
        include a column "reference" for reference genomes.
        Curation parameters can also be set in the table using column headers: 
        "unknowns", "contigs", "assembly_size", and "distance". 
        If a parameter cell is left blank it will be replaced with the 
        default value or a value passed from the command line
        """
    )

    curate = parser.add_argument_group(
        'Curate', description="Curate genomes"
    )

    curate.add_argument(
        "--curate_only", action="store_true", default=False,
        help="Run only the curate step"
    )
    
    curate.add_argument(
        "--unknowns", '-n', type=types['nonneg_int'], nargs='+', default=[200],
        help="""
        Maximum number of unknown bases (not A, T, C, G) for curation. 
        If more than one 
        value is passed, the list must be the same length as the number of
        species. Otherwise the same value is applied to all species."""
    )

    curate.add_argument(
        "--contigs", '-c', type=types['nonneg_float'], nargs='+', default=[3.0],
        help="""
        Acceptable deviations from median number of contigs for curation. 
        If more than one 
        value is passed, the list must be the same length as the number of 
        species. Otherwise the same value is applied to all species"""
    )

    curate.add_argument(
        "--assembly_size", '-a', type=types['nonneg_float'], nargs='+',
        default=[3.0],
        help="""
        Acceptable devations from median assembly size for curation.
        If more than one 
        value is passed, the list must be the same length as the number of 
        species. Otherwise the same value is applied to all species"""
    )

    curate.add_argument(
        "--distance", '-d', type=types['nonneg_float'], nargs='+',
        default=[3.0],
        help="""
        Acceptable deviations from median MASH distances for curation.
        If more than one 
        value is passed, the list must be the same length as the number of 
        species. Otherwise the same value is applied to all species"""
    )

    tree = parser.add_argument_group(
        'Tree', description="SNP matrix and tree generation"
    )
    tree.add_argument(
        "--reference", '-r', type=str, nargs='+',
        help="""
        Define which genome to use as reference by providing accession number
        (GCA_XXXXXXXXX). This list should be the same length as the number of
        species. This is REQUIRED from the command line or in the config file.
        """
    )

    parser.add_argument(
        "--log", '-l', type=types['outfile_type'],
        default=os.path.join("./", DEFAULT_LOG_FNAME),
        help="Set log file path"
    )

    return parser


# params from args to work on
PARAMS = ["unknowns", "contigs", "assembly_size", "distance", "reference"]
PARAM_TYPES = [
    types['nonneg_int'], types['nonneg_float'],
    types['nonneg_float'], types['nonneg_float'], str]

def get_user_provided_params(argv):
    user_params = []
    for param in PARAMS:
        if "--{}".format(param) in argv or "-{}".format(param[0]) in argv:
            user_params.append(param)
    return user_params
            
def validate_types_in_config(column_data, types):
    validated = []
    for value in column_data:
        validated.append(types(value))
    return validated

def check_and_adjust_length_of_arg_list(args_dict, n_species):
    param_data = {}
    for param in PARAMS:
        if param in args_dict:
            if args_dict[param] is None:
                continue
            elif (len(args_dict[param]) > 1):
                if len(args_dict[param]) != n_species:
                    problem = """
                        The number of values for parameter {0} must be the same
                        as the number of species ({1})
                        """.format(param, n_species)
                    LOG.exception(problem)
                    raise ValueError(problem)
                else:
                    param_data[param] = args_dict[param]
                    # param_data.append(args_dict[param])
            else:
                param_data[param] = args_dict[param] * n_species
    return pd.DataFrame(param_data, columns=PARAMS)

def get_arg_df_from_cml(args, n_species):
    args_dict = vars(args)
    df = check_and_adjust_length_of_arg_list(args_dict, n_species)
    df['species'] = args_dict['species']
    df = df.set_index('species')
    return df
     
def get_arg_df_from_config(args):
    try:
        df = pd.read_csv(args.config, index_col="species", comment='#')
    except ValueError:
        LOG.exception(
            """Required 'species' column is not 
            present in config file: {}""".format(args.config))
        raise
    except Exception:
        LOG.exception(
            "Problem reading config file: {}".format(args.config)
        )
        raise
    return df


def combine_cml_and_cfg_params(cfg, cml, user_params):
    combined_data = []
    for param, param_type in zip(PARAMS, PARAM_TYPES):
        # check if value was overidden from the command line
        if param in user_params:
            # return value from cml df
            combined_data.append(list(cml[param]))
            continue
        try:
            if param in cfg:
                # add column from config if it exists and fill any
                # missing values with default cml args
                combined_data.append(
                    validate_types_in_config(
                        cfg[param].fillna(cml[param].iloc[0]),
                    param_type))
            else:
                # add cml args if column is not present in config
                combined_data.append(
                    list(cml[param])
                )
        except KeyError:
            # this should only occur if reference params are not provided
            # because there are no default values
            LOG.exception("Reference information must be provided")
            raise
    return pd.DataFrame(
        np.array(combined_data).T, columns=PARAMS, index=cfg.index)


def get_arg_df(args, user_params):
    if args.species is not None:
        cml_arg_data = get_arg_df_from_cml(args, len(args.species))
        try:
            cml_arg_data['reference']
        except KeyError:
            LOG.exception("Reference information is required")
            raise
        return cml_arg_data
    else:
        cfg_arg_data = get_arg_df_from_config(args)
        cml_arg_data = get_arg_df_from_cml(args, len(cfg_arg_data.index))
        return combine_cml_and_cfg_params(
            cfg_arg_data, cml_arg_data, user_params)


def unlock_snake(snakefile, work_path, params, snakemake_args):
    cmd = (
            'snakemake --snakefile {snakefile} -d {outpath} '
            '--config params={params}'.format(
                snakefile=snakefile,
                outpath=work_path,
                params=params,
                ).split())
    cmd += snakemake_args
    sp.run(cmd, check=True)

def run_snake(snakefile, work_path, params, snakemake_args):
    cmd = (
            'snakemake --snakefile {snakefile} -d {outpath} '
            '--config params={params}'.format(
                snakefile=snakefile,
                outpath=work_path,
                params=params,
                ).split())
    cmd += snakemake_args
    sp.run(cmd, check=True)
    

def download_sequences(species, paths, update, update_assembly):
    message = "Downloading sequences for species: {0}\nWriting to {1}".format(
        ", ".join(species), paths.genomes)
    LOG.info(message)
    click.secho(message, err=False, fg='green')
    cmd = ["ncbitk", update, update_assembly,
           paths.genomes] + species
    LOG.info("Running command: {}".format(" ".join(cmd)))
    sp.check_output(cmd)
    message = "Finished downloading sequences to: {}".format(
        paths.genomes)
    LOG.info(message)
    click.secho(message, err=False, fg="green")


def curate(species, paths, params, threads):
    params = json.loads(params)
    for spec in species:
        message = (
            "Running genbankqc on {0} using the "
            "following parameters {1} and using {2} threads".format(
                spec, " ".join(
                    ["{0}={1}".format(k, v) for k, v in params[spec].items()]),
                    threads))
        LOG.info(message)
        click.secho(message, err=False, fg='green')

        sp.call(
            list(map(str, ["genbankqc", "species",
            os.path.join(paths.genomes, spec),
            "-n", params[spec]['unknowns'],
            "-c", params[spec]['contigs'],
            "-s", params[spec]['assembly_size'],
            "-d", params[spec]['distance'],
            "--processes", threads])))




def run(params_json, args, snakemake_args):
    paths = make_dirs(args.path)
    snek_path = os.path.dirname(__file__)
    nasp_snek = os.path.join(snek_path, "nasp_snek")
    curate_snek = os.path.join(snek_path, "curate_snek")
    if "--unlock" in snakemake_args:
        unlock_snake(nasp_snek, paths.wrk_dir, params_json, snakemake_args)
    else:
        update = "--no-update" if args.no_update else "--update"
        update_assembly = (
            "--local-assembly" if args.no_assembly_update
            else "--update-assembly")
        species = list(json.loads(params_json).keys())
        if not args.curate_only:
            download_sequences(species, paths, update, update_assembly)
        # curate(species, paths, params_json, args.threads)
        if not args.download_only:
            LOG.info("Running Curation")
            run_snake(
                curate_snek,
                paths.wrk_dir, params_json, snakemake_args
            )
            LOG.info("Generating matrix and tree")
            run_snake(
                nasp_snek,
                paths.wrk_dir, params_json, snakemake_args + ["--use-conda"]
            )

PATHS = namedtuple('paths', ['genomes', 'logs', 'wgfast', 'wrk_dir'])

def make_dirs(outpath):
    genomes = os.path.join(outpath, "genomes")
    logs = os.path.join(outpath, "logs")
    wgfast = os.path.join(outpath, "wgfast")
    wrk_dir = outpath
    os.makedirs(genomes, exist_ok=True)
    os.makedirs(logs, exist_ok=True)
    os.makedirs(wgfast, exist_ok=True)
    return PATHS(genomes=genomes, logs=logs, wgfast=wgfast, wrk_dir=wrk_dir)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = parse()
    # Return help if no command is passed
    if len(argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(0)
    user_params = get_user_provided_params(argv)
    args, snakemake_args = parser.parse_known_args(argv[1:])
    config_logging(args.log)
    arg_df = get_arg_df(args, user_params)
    LOG.info("Parameters set to:\n{}".format(arg_df))
    run(arg_df.to_json(orient="index"), args, snakemake_args)

if __name__ == "__main__":
    main()


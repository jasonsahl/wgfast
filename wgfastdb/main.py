import argparse
import sys
import datetime
import os
import logging
import subprocess as sp
import pandas as pd

from wgfastdb.parse_util import types, error, warn


DEFAULT_LOG_FNAME = wgfastdb.log

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

    dowload.add_argument(
        "--update", action='store_true', default=False,
        help="""
        Sync your collection with the latest assembly versions
        """
    )
    download.add_argument(
        "--update-assembly", action="store_true", default=False,
        help="""
        Download the latest assembly summary and taxonomy
        dump or use your local copies.
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

    parser.add_argument(
        "--threads", '-t', type=types['positive_int'], default=4,
        help="Number of worker threads to spawn."
    )

    parser.add_argument(
        'snakemake_args', nargs=argparse.REMAINDER,
        help="Additional snakemake arguments"
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

def config_logging(file_name, level="INFO"):
    """Configure logging
    Arguments:
        level (str): Logging level
    """
    logging.basicConfig(
        filename=file_name,
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=getattr(logging, level),
        filemode='w',
        format='%(asctime)s %(levelname)s: [%(name)s] %(message)s')
    return logging.getLogger(__name__)




def check_and_adjust_length_of_arg_list(args_dict, n_species):
    param_data = []
    for param in PARAMS:
        if param in args_dict:
            if (len(args_dict[param]) > 1):
                if len(args_dict[param]) != n_species:
                    problem = """
                        The number of values for parameter {0} must be the same
                        as the number of species ({1})
                        """.format(param, n_species)
                    logger.exception(problem)
                    raise ValueError(problem)
                else:
                    param_data.append(args_dict[param])
            else:
                param_data.append(args_dict[param] * n_species)
    return pd.DataFrame(param_data, columns=PARAMS)

def get_arg_df_from_cml(args):
    args_dict = vars(args)
    df = check_and_adjust_length_of_arg_list(args_dict, len(args.species))
    df['species'] = args_dict['species']
    df = df.set_index('species')
    return df
     
def get_arg_df_from_config(args):
    try:
        df = pd.read_csv(args.config, index_col="species")
    except ValueError:
        logger.exception(
            """Required 'species' column is not 
            present in config file: {}""".format(args.config))
        raise
    except Exception:
        logger.exception(
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
            logger.exception("Reference information must be provided")
            raise
    return pd.DataFrame(combined_data, columns=PARAMS, index=cfg['species'])
        


    

def get_arg_df(args, user_params):
    if 'species' in args:
        cml_arg_data = get_arg_data_from_cml(args)
        try:
            cml_arg_data['reference']
        except KeyError:
            logger.exception("Reference information is required")
            raise
        return cml_arg_data
    else:
        cfg_arg_data = get_arg_data_from_config(args)
        cml_arg_data = get_arg_data_from_cml(
            args, len(cfg_arg_data['species']))
        return combine_cml_and_cfg_params(
            cfg_arg_data, cml_arg_data, user_params)
        
def run(params_json, outpath, threads, snakemake_args):
    script = os.path.join(os.path.dirname(__file__), "wgfastdb_tree.py")
    update = "--update" if args.update else "--no-update"
    update_assembly = (
        "--update-assembly" if args.update_assembly else "--local-assembly")
    cmd = 'snakemake --config path={outpath} threads={threads} '
            'params={params} script={script} update={update}, '
            'update_assembly={update_assembly}'.format(
                outpath=outpath,
                threads=threads,
                params=params_json,
                script=script,
                update=update,
                update_assembly=update_assembly).split()
    cmd += snakemake_args
    try:
        sp.call(cmd, check=True)
    except ( KeyboardInterrupt, sp.CalledProcessError) as e:
                warn("Unlocking directory after failed snakemake")
                sp.run(cmd + ["--unlock"], check=True )
                error(e)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = parse()
    # Return help if no command is passed
    if len(argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(0)
    user_params = get_user_provided_params(argv)
    args, snakemake_args = parser.parse_known_args(argv)
    logger = config_logging(args.log)
    arg_df = get_arg_df(args, user_params)
    logger.info("Parameters set to:\n{}".format(arg_df))
    run(arg_df.to_json, args.path, args.threads, snakemake_args)

if __name__ == "__main__":
    main()

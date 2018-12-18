import argparse
import os
import click
import textwrap
import logging
import sys
from glob import glob

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

def warn(msg):
    '''Prints a warning message to stderr.'''
    text = quote(msg, quote="!> ", nl=False)
    click.secho(text, err=True, fg='red')


def quote(text, width=72, quote="", nl=True):
    if not text:
        return ""
    out = ""
    for line in text.split("\n"):
        sublines = textwrap.wrap(line, width=width, replace_whitespace=False)
        sublines = [quote + l for l in sublines]  # if l.strip()]
        out += "\n".join(sublines) + "\n"
    return out

class WGFASTError(Exception):
    pass

def error(msg, exception=True, wrap=True):
    '''Prints an error message to logger and stderr and exits.'''
    if exception:
        raise WGFASTError(msg)
    else:
        msg = quote(msg, width=75) if wrap else msg
        click.secho(msg, fg='red', err=True)
        sys.exit(1)

def path_type(input_path):
    """Path_type is for paths that should already exist.
    If path exists, returns absolute path, otherwise raise
    ArgumentTypeError"""
    if not os.path.isdir(input_path):
        raise argparse.ArgumentTypeError("Not a valid path")
    return os.path.abspath(input_path)


def project_dir_type(input_path):
    """Creates a project directory if one does not exist and
    Throws PermissionError if there are no permissions to create
    directory. Returns absolute path to directory."""
    input_path = os.path.abspath(input_path)
    if not os.path.isdir(input_path):
        try:
            os.makedirs(input_path)
        except PermissionError:
            error(
                "No permission to make directory: {}".format(input_path))
    return input_path


def project_dir_type_cd(input_path):
    """Creates a project directory if one does not exist and
    Throws PermissionError if there are no permissions to create
    directory. Returns absolute path to directory. 
    Changes cwd to input_path"""
    project_dir_type(input_path)
    if os.getcwd() != input_path:
        os.chdir(input_path)
    return input_path


def path_list_type(input_paths):
    """Returns list of absolute paths using glob"""
    glob_list = []
    for path in input_paths:
        glob_list += glob(path)
    if not glob_list:
        raise argparse.ArgumentTypeError("Not a valid path")
    return list({os.path.abspath(path) for path in glob_list})


def outpath_type(input_path):
    """Outpath_type creates a directory if one does not exist.
    Throws PermissionError if there are no permissions to create
    directory. If path already exists and it is not empty, a warning
    is issued. Returns absolute path to directory."""
    input_path = os.path.abspath(input_path)
    try:
        os.makedirs(input_path)
    except PermissionError:
        error(
            "No permission to make directory: {}".format(input_path))
    except OSError:
        warn("Directory already exists: {}. ".format(input_path))
        if os.listdir(input_path):
            warn("Files in {} may be overwritten!".format(input_path))
    return input_path


def file_type(input_file):
    """File type is for input files that
    should already exist.
    Check if input_file is a file.
    If yes, return absolute path to
    file else throw ArgumentTypeError
    exception"""
    input_file = os.path.abspath(input_file)
    if not os.path.isfile(input_file):
        raise argparse.ArgumentTypeError(
            "Not a valid file path: {}".format(input_file))
    return input_file


def outfile_type(input_file):
    """Checks that path to file exists,
    if it doesn't, the path is created,
    returns abs path to file. PermissionError
    exception when there is no permission to 
    create directory"""
    input_file = os.path.abspath(input_file)
    path = os.path.dirname(input_file)
    if not os.path.isdir(path):
        try:
            os.makedirs(path)
        except PermissionError:
            error(
                "No permission to create file: {}".format(input_file)
            )
    return input_file

def write_handle_type(input_file):
    """Returns a handle to an outfile"""
    return open(outfile_type(input_file), 'w')


def read_handle_type(input_file):
    """Returns a handle to an infile"""
    return open(file_type(input_file), 'r') 


def nonneg_int(input_val):
    """Make a non negative int type for argparse"""
    try:
        input_val = int(input_val)
        if input_val < 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid non-negative integer")
    return input_val


def proportion(input_val):
    """make proportion type for argparse"""
    try:
        input_val = float(input_val)
        if input_val > 1 or input_val < 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a proportion")
    return input_val


def positive_int(input_val):
    """Make a positive int type for argparse."""
    try:
        input_val = int(input_val)
        if input_val <= 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a positive integer")
    return input_val

def positive_float(input_val):
    """ Make a positive float type for argparse."""
    try:
        input_val = float(input_val)
        if input_val <= 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a positive float")
    return input_val

def nonneg_float(input_val):
    """ Make a non-negative float type for argparse."""
    try:
        input_val = float(input_val)
        if input_val < 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid non-negative float")
    return input_val

def fasta(header, seq):
    return ">{0}\n{1}\n".format(header, seq)

types = {
    'int': int,
    'str': str,
    'float': float,
    'bool': bool,
    'file_type': file_type,
    'outfile_type': outfile_type,
    'path_type': path_type,
    'outpath_type': outpath_type,
    'positive_int': positive_int,
    'nonneg_int': nonneg_int,
    'proportion': proportion,
    'project_dir_type': project_dir_type,
    'write_handle_type': write_handle_type,
    'read_handle_type': read_handle_type,
    'path_list_type': path_list_type,
    'positive_float': positive_float,
    'nonneg_float': nonneg_float
    }
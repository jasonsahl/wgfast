import os
import re
import click
import logbook

from collections import namedtuple

from genbankqc import Genbank
from genbankqc import Genome
from genbankqc import Species


class CLIGroup(click.Group):
    def parse_args(self, ctx, args):
        try:
            if args[0] in self.commands:
                if len(args) == 1 or args[1] not in self.commands:
                    args.insert(0, '')
        except IndexError:
            pass
        super(CLIGroup, self).parse_args(ctx, args)


@click.group(invoke_without_command=True, no_args_is_help=True, cls=CLIGroup)
@click.pass_context
@click.argument('path', type=click.Path(), required=False)
def cli(ctx, path):
    # TODO: Option for basic info about PATH
    """
    Assess the integrity of your genomes through automated analysis of
    species-based statistics and metadata.
    """
    genbank = Genbank(path)
    _ctx = namedtuple('ctx', ['genbank', 'assembly_summary'])
    ctx.obj = _ctx(genbank=genbank, assembly_summary=genbank.assembly_summary)
    logbook.set_datetime_format("local")
    log_dir = os.path.join(path, ".logs")
    log_file = os.path.join(log_dir, "qc.log")
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)
    handler = logbook.TimedRotatingFileHandler(log_file, backup_count=10)
    handler.push_application()
    if ctx.invoked_subcommand is None:
        genbank.qc()


@cli.command()
@click.pass_obj
@click.argument('path', type=click.Path(exists=True, file_okay=False))
@click.option('--unknowns', '-n',
              type=int, default=200,
              help='Maximum number of unknown bases (not A, T, C, G)')
@click.option('--contigs', '-c',
              type=float, default=3.0,
              help='Acceptable deviations from median number of contigs')
@click.option('--assembly_size', '-s',
              type=float, default=3.0,
              help='Acceptable deviations from median assembly size')
@click.option('--distance', '-d',
              type=float, default=3.0,
              help='Acceptable deviations from median MASH distances')
@click.option('--all',
              type=float,
              help='Acceptable deviations for all metrics')
@click.option('--metadata', is_flag=True,
              help='Get metadata for genome at PATH',)
@click.option('--processes', type=int, default=1,
              help="Number of processes")
def species(ctx, path, unknowns, contigs, assembly_size, distance, all,
            metadata, processes):
    """Run commands on a single species"""
    kwargs = {"max_unknowns": unknowns,
              "contigs": contigs,
              "assembly_size": assembly_size,
              "mash": distance,
              "assembly_summary": ctx.assembly_summary,
              "processes": processes}
    logbook.set_datetime_format("local")
    log_dir = os.path.join(path, ".logs")
    log_file = os.path.join(log_dir, "qc.log")
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)
    handler = logbook.TimedRotatingFileHandler(log_file, backup_count=10)
    handler.push_application()
    species = Species(path, **kwargs)
    species.qc()
    if metadata:
        species.metadata()


@cli.command()
@click.pass_obj
@click.argument('path', type=click.Path(exists=True, dir_okay=False))
@click.option('--metadata', is_flag=True,
              help='Get metadata for genome at PATH')
def genome(ctx, path, metadata):
    """
    Get information about a single genome.
    """

    genome = Genome(path, ctx.assembly_summary)
    if metadata:
        click.echo(genome.metadata)


@cli.command()
@click.pass_obj
@click.argument('path', type=click.Path(exists=True, dir_okay=False))
def log_stats(ctx, path):
    """
    Summarize basic stats for a given log file.
    """
    log_file = os.path.join(ctx.genbank.path, '.logs', path)
    not_enough_genomes = [0, "Not enough genomes"]
    completed_metadata_command = [0, "Completed metadata command"]
    already_complete = [0, "Already complete"]
    tree_already_complete = [0, "Tree already complete"]
    generated_stats = [0, "Generated stats"]
    qc_completed = [0, "qc command completed"]
    stats = [log_file, not_enough_genomes,
             completed_metadata_command,
             already_complete,
             tree_already_complete,
             generated_stats,
             qc_completed]
    with open(log_file) as f:
        for line in f:
            if re.search(not_enough_genomes[1], line):
                not_enough_genomes[0] += 1
            elif re.search(completed_metadata_command[1], line):
                completed_metadata_command[0] += 1
            elif re.search(already_complete[1], line):
                already_complete[0] += 1
            elif re.search(tree_already_complete[1], line):
                tree_already_complete[0] += 1
            elif re.search(generated_stats[1], line):
                generated_stats[0] += 1
            elif re.search(qc_completed[1], line):
                qc_completed[0] += 1
    for i in stats:
        click.echo(i[1])
        click.echo(i[0])

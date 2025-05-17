import argparse
import sys
import os
import numpy as np

from rich_argparse import ArgumentDefaultsRichHelpFormatter

from . import __version__
from . import index, fetch, profile, rebuild
from .log import log

## customize formatter
ArgumentDefaultsRichHelpFormatter.styles['argparse.prog'] = 'default'
ArgumentDefaultsRichHelpFormatter.styles['argparse.default'] = 'grey50'
ArgumentDefaultsRichHelpFormatter.styles['argparse.metavar'] = 'grey50'
ArgumentDefaultsRichHelpFormatter.styles['argparse.groups'] = '#DF8080'
ArgumentDefaultsRichHelpFormatter.styles['argparse.args'] = '#85C2DE'


def parser_index(parser):
    parser_index = parser.add_parser(
        name='index',
        help='Update an existing reference database or construct a new one using MAGs.',
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False)

    parser_index.add_argument(
        dest='files',
        nargs='+',
        metavar='file',
        help='Input fasta <*.fa|*.fna|*.fasta> file(s), gzip optional <*.gz>.')

    required = parser_index.add_argument_group('required arguments')
    required.add_argument(
        '-o',
        '--outdir',
        required=True,
        metavar='DIR',
        help='Output directory.')

    optional = parser_index.add_argument_group('optional arguments')
    optional.add_argument(
        '--compress',
        action='store_true',
        help='Compress the database by creating a tar archive.')

    optional.add_argument(
        '-d',
        '--database',
        metavar='DIR',
        help='Database directory. If given then extend it, otherwise build a new one.')

    optional.add_argument(
        '-a',
        '--taxonomy',
        metavar='FILE',
        help="Path to GTDB-Tk's <gtdbtk.bac120.summary.tsv> or GTDB's <bac120_taxonomy.tsv>. If not given then ignore taxonomy mapping.")

    optional.add_argument(
        '-t',
        '--threads',
        metavar='INT',
        type=int,
        default=os.cpu_count(),
        help=f'Number of threads.')

    additional = parser_index.add_argument_group('additional arguments - sketching')
    additional.add_argument(
        '-k',
        '--kmer',
        dest='k',
        metavar='INT',
        type=int,
        default=31,
        help='K-mer size.')

    additional.add_argument(
        '-s',
        '--scale',
        dest='s',
        metavar='INT',
        type=int,
        default=250,
        help='Scale for downsampling.')

    additional.add_argument(
        '-w',
        '--window',
        dest='w',
        metavar='INT',
        type=int,
        default=25000,
        help='Window size for k-mer grouping.')

    parser_index.add_argument('-v', '--version', action='version', version=__version__, help=argparse.SUPPRESS)
    parser_index.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)
    parser_index.set_defaults(func=index)

def parser_fetch(parser):
    parser_fetch = parser.add_parser(
        name='fetch',
        help='Download the pre-built GTDB database from Zenodo.',
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False)

    required = parser_fetch.add_argument_group('required arguments')
    required.add_argument(
        '-o',
        '--outdir',
        required=True,
        metavar='DIR',
        help='Output directory.')

    parser_fetch.add_argument('-v', '--version', action='version', version=__version__, help=argparse.SUPPRESS)
    parser_fetch.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)
    parser_fetch.set_defaults(func=fetch)

def parser_profile(parser):
    parser_profile = parser.add_parser(
        name='profile',
        help='Profile bacteria growth dynamics from fasta|fastq samples.',
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False)

    parser_profile.add_argument(
        dest='files',
        nargs='+',
        metavar='file',
        help='Input fasta <*.fa|*.fasta> or fastq <*.fq|*.fastq> file(s), gzip optional <*.gz>.')

    required = parser_profile.add_argument_group('required arguments')
    required.add_argument(
        '-o',
        '--outdir',
        required=True,
        metavar='DIR',
        help='Output directory.')

    required.add_argument(
        '-d',
        '--database',
        required=True,
        metavar='DIR',
        help='Database directory.')

    optional = parser_profile.add_argument_group('optional arguments')
    optional.add_argument(
        '--force',
        action='store_true',
        help='Force counting.')

    optional.add_argument(
        '--single',
        action='store_true',
        help='Files are single-end. If not given then merge forward|reward files with <_(1|2)>, <_(R1|R2)> or <_(fwd|rev)>.')

    optional.add_argument(
        '-t',
        '--threads',
        metavar='INT',
        type=int,
        default=os.cpu_count(),
        help=f'Number of threads.')

    additional = parser_profile.add_argument_group('additional arguments - parsing')
    additional.add_argument(
        '-x',
        '--min-dept',
        metavar='FLOAT',
        type=float,
        default=5,
        help='Min. median (window) depth of sketches.')

    additional.add_argument(
        '-y',
        '--max-disp',
        metavar='FLOAT',
        type=float,
        default=np.inf,
        help="Max. median (window) dispersion of sketches' counts.")

    additional.add_argument(
        '-z',
        '--min-frac',
        metavar='FLOAT',
        type=float,
        default=0.75,
        help="Min. fraction of reference genomes' windows covered by sketches.")

    additional.add_argument(
        '-c',
        '--min-cont',
        metavar='FLOAT',
        type=float,
        default=0.25,
        help="Min. containment of reference genomes' sketches after reassignment of shared k-mers.")

    additional = parser_profile.add_argument_group('additional arguments - fitting')
    additional.add_argument(
        '-n',
        '--components',
        metavar='INT',
        type=int,
        default=5,
        help='Number of mixture components.')

    additional.add_argument(
        '-i',
        '--max-iter',
        metavar='INT',
        type=int,
        default=np.inf,
        help='Terminal condition for EM - max. number of iterations.')

    additional.add_argument(
        '-e',
        '--tol',
        metavar='FLOAT',
        type=float,
        default=1e-5,
        help='Terminal condition for EM - tolerance.')

    parser_profile.add_argument('-v', '--version', action='version', version=__version__, help=argparse.SUPPRESS)
    parser_profile.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)
    parser_profile.set_defaults(func=profile)

def parser_rebuild(parser):
    parser_rebuild = parser.add_parser(
        name='rebuild', 
        help='Rebuild a reference database using GTDB representative bacterial genomes.',
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False)

    required = parser_rebuild.add_argument_group('required arguments')
    required.add_argument(
        '-o',
        '--outdir',
        required=True,
        metavar='DIR',
        help='Output directory.')

    optional = parser_rebuild.add_argument_group('optional arguments')
    optional.add_argument(
        '-t',
        '--threads',
        metavar='INT',
        type=int,
        default=os.cpu_count(),
        help=f'Number of threads.')

    additional = parser_rebuild.add_argument_group('additional arguments - sketching')
    additional.add_argument(
        '-k',
        '--kmer',
        dest='k',
        metavar='INT',
        type=int,
        default=31,
        help='K-mer size.')

    additional.add_argument(
        '-s',
        '--scale',
        dest='s',
        metavar='INT',
        type=int,
        default=250,
        help='Scale for downsampling.')

    additional.add_argument(
        '-w',
        '--window',
        dest='w',
        metavar='INT',
        type=int,
        default=25000,
        help='Window size for k-mer grouping.')

    parser_rebuild.add_argument('-v', '--version', action='version', version=__version__, help=argparse.SUPPRESS)
    parser_rebuild.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)
    parser_rebuild.set_defaults(func=rebuild)

def cli():
    parser = argparse.ArgumentParser(
        prog='pilea',
        description=f'Pilea v{__version__}: profiling bacterial growth dynamics from metagenomes with sketching',
        formatter_class=ArgumentDefaultsRichHelpFormatter
    )
    parser.add_argument('-v', '--version', action='version', version=__version__)
    subparsers = parser.add_subparsers(title='subcommands')

    parser_index(subparsers)
    parser_fetch(subparsers)
    parser_profile(subparsers)
    parser_rebuild(subparsers)

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    args.func(**{key: val for key, val in vars(args).items() if key != 'func'})

if __name__ == '__main__':
    cli()

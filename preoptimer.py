import argparse
import sys
from pathlib import Path

from src import PreoptimerError
from src.parameterization import Parameterization
from src.preoptimization import Preoptimization


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
    return number


def load_arguments():
    parser = argparse.ArgumentParser(
        description='PROPTIMUS PREOPTIMER: a parameterizable tool for setting bond lengths '
                    'and bond angles in protein structures.')

    required = parser.add_argument_group('required arguments')
    required.add_argument('--mode',
                          type=str,
                          required=True,
                          choices=['preoptimize', 'parameterize'],
                          help='Mode of the program. Mode "preoptimize" sets bond lengths and bond angles '
                               'according to the parameters. Mode "parameterize" creates these parameters '
                               'from a set of protein structures.')
    required.add_argument('--parameters',
                          type=str,
                          required=True,
                          metavar='FILE',
                          help='JSON file with parameters. In mode "preoptimize", the parameters are read from it '
                               '(available parameters are in the "parameters" directory), in mode "parameterize", '
                               'the calculated parameters are saved to it (the file must not exist).')

    preoptimize = parser.add_argument_group('arguments for --mode preoptimize')
    preoptimize.add_argument('--input_PDB_file',
                             type=str,
                             metavar='FILE',
                             help='PDB file with structure, which should be preoptimized.')
    preoptimize.add_argument('--output_PDB_file',
                             type=str,
                             metavar='FILE',
                             help='PDB file to store preoptimized structure (the file must not exist).')
    preoptimize.add_argument('--only_hydrogens',
                             action='store_true',
                             help='Preoptimize only hydrogen atoms, positions of all other atoms are kept.')

    parameterize = parser.add_argument_group('arguments for --mode parameterize')
    parameterize.add_argument('--proteins_directory',
                              type=str,
                              metavar='DIR',
                              help='Directory with protein structures (*.pdb files) used to calculate the parameters.')
    parameterize.add_argument('--workers',
                              type=positive_int,
                              metavar='N',
                              help='Number of worker processes (default: number of CPUs).')
    args = parser.parse_args()
    if args.mode == 'preoptimize':
        missing = [a for a in ('input_PDB_file', 'output_PDB_file') if getattr(args, a) is None]
        if missing:
            parser.error('mode preoptimize requires: ' + ', '.join('--' + a for a in missing))
        for name in ('input_PDB_file', 'parameters'):
            if not Path(getattr(args, name)).is_file():
                parser.error(f'--{name} "{getattr(args, name)}" does not exist.')
        if Path(args.output_PDB_file).exists():
            parser.error(f'--output_PDB_file "{args.output_PDB_file}" already exists.')
    elif args.mode == 'parameterize':
        if args.proteins_directory is None:
            parser.error('mode parameterize requires: --proteins_directory')
        if not Path(args.proteins_directory).is_dir():
            parser.error(f'--proteins_directory "{args.proteins_directory}" does not exist or is not a directory.')
        if Path(args.parameters).exists():
            parser.error(f'--parameters "{args.parameters}" already exists.')
    return args


def main():
    args = load_arguments()
    try:
        if args.mode == 'preoptimize':
            Preoptimization(args.input_PDB_file, args.parameters, args.output_PDB_file,
                            only_hydrogens=args.only_hydrogens).run()

        elif args.mode == 'parameterize':
            Parameterization(args.proteins_directory, args.parameters, workers=args.workers).run()
    except PreoptimerError as error:
        sys.exit(f"error: {error}")


if __name__ == '__main__':
    main()

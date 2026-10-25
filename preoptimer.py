import argparse
from pathlib import Path

from src.protein import Protein
from src.make_html_graphs import make_html_graphs
from src.parameterization import Parameterization


def load_arguments():
    print("\nParsing arguments... ", end="", flush=True)
    parser = argparse.ArgumentParser(
        description='PROPTIMUS PREOPTIMER: a parameterizable tool for setting bond lengths, '
                    'bond angles and dihedral angles in protein structures.')

    required = parser.add_argument_group('required arguments')
    required.add_argument('--mode',
                          type=str,
                          required=True,
                          choices=['preoptimize', 'parameterize'],
                          help='Mode of the program. Mode "preoptimize" sets bond lengths, bond angles '
                               'and dihedral angles according to the parameters. Mode "parameterize" '
                               'creates these parameters from a set of molecules.')

    preoptimize = parser.add_argument_group('arguments for --mode preoptimize')
    preoptimize.add_argument('--input_PDB_file',
                             type=str,
                             metavar='FILE',
                             help='PDB file with structure, which should be preoptimised.')
    preoptimize.add_argument('--output_PDB_file',
                             type=str,
                             metavar='FILE',
                             help='PDB file to store preoptimized structure.')
    preoptimize.add_argument('--parameters',
                             type=str,
                             metavar='FILE',
                             help='JSON file with parameters. Available parameters are in the "parameters" directory.')
    preoptimize.add_argument('--only_hydrogens',
                             action='store_true',
                             help='Preoptimize only hydrogen atoms, positions of all other atoms are kept.')
    preoptimize.add_argument('--noopt',
                             action='store_true',
                             help='Skip the short AMBER optimization (OpenMM) run before the preoptimization. '
                                  'Without --only_hydrogens it optimizes all atoms except CA, '
                                  'with --only_hydrogens only hydrogens.')

    parameterize = parser.add_argument_group('arguments for --mode parameterize')
    parameterize.add_argument('--proteins_directory',
                              type=str,
                              metavar='DIR',
                              help='Directory with proteins used to calculate the parameters.')
    parameterize.add_argument('--data_dir',
                              type=str,
                              metavar='DIR',
                              help='Directory to store the results.'
                                   'The calculated parameters are saved to parameters.json and HTML graphs of all '
                                   'atom types of all residues to the "graphs" subdirectory.')
    parameterize.add_argument('--workers',
                              type=int,
                              metavar='N',
                              help='Number of worker processes (default: number of CPUs).')
    parameterize.add_argument('--open_browser',
                              action='store_true',
                              help='Open the HTML graphs in the default web browser when they are ready.')
    args = parser.parse_args()
    if args.mode == 'preoptimize':
        missing = [a for a in ('input_PDB_file', 'output_PDB_file', 'parameters') if getattr(args, a) is None]
        if missing:
            parser.error('mode preoptimize requires: ' + ', '.join('--' + a for a in missing))
        for name in ('input_PDB_file', 'parameters'):
            if not Path(getattr(args, name)).is_file():
                parser.error(f'--{name} "{getattr(args, name)}" does not exist.')
        if Path(args.output_PDB_file).exists():
            parser.error(f'--output_PDB_file "{args.output_PDB_file}" already exists.')
    elif args.mode == 'parameterize':
        missing = [a for a in ('proteins_directory', 'data_dir') if getattr(args, a) is None]
        if missing:
            parser.error('mode parameterize requires: ' + ', '.join('--' + a for a in missing))
        if not Path(args.proteins_directory).is_dir():
            parser.error(f'--proteins_directory "{args.proteins_directory}" does not exist or is not a directory.')
        if Path(args.data_dir).exists():
            parser.error(f'--data_dir "{args.data_dir}" already exists.')
    print("ok")
    return args


def main():
    args = load_arguments()
    if args.mode == 'preoptimize':
        protein = Protein(args.input_PDB_file)
        protein.preoptimize(args.parameters, only_hydrogens=args.only_hydrogens,
                             amber_optimization=not args.noopt)
        protein.save(args.output_PDB_file)

    elif args.mode == 'parameterize':
        data_dir = Path(args.data_dir)
        data_dir.mkdir(parents=True)
        parameterization = Parameterization(args.proteins_directory, str(data_dir / 'parameters.json'),
                                            workers=args.workers)
        parameterization.run()
        make_html_graphs(parameterization.data, parameterization.parameters, str(data_dir / 'graphs'),
                         open_browser=args.open_browser)
        parameterization.find_outliers()


if __name__ == '__main__':
    main()


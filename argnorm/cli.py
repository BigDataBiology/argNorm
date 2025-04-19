import argparse
import sys
import textwrap

from .argnorm_version import __version__, is_release
from .atomicwrite import atomic_write
from . import lib


def main():
    """
    Major function to run when running `argnorm` in shell.
    """
    parser = argparse.ArgumentParser(
        description=('argNorm normalizes ARG annotation results from '
                     'different tools and databases to the same ontology, '
                     'namely ARO (Antibiotic Resistance Ontology).'),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\

        If argNorm is helpful in a scientific publication, please cite:

            argNorm: normalization of antibiotic resistance gene annotations to the Antibiotic Resistance Ontology (ARO)
            by Svetlana Ugarcina Perovic, Vedanth Ramji, Hui Chong, Yiqian Duan, Finlay Maguire, Luis Pedro Coelho
            in Bioinformatics (2025) https://doi.org/10.1093/bioinformatics/btaf173
         '''))

    parser.add_argument('--version', '-v',
                        action='version',
                        version=(__version__ if is_release else f'{__version__}-dev')
                        )
    parser.add_argument('tool', type=str.lower,
                        choices=['argsoap', 'abricate', 'deeparg', 'resfinder', 'amrfinderplus', 'groot', 'hamronization'],
                        help='tool (required): The bioinformatics tool used for ARG annotation.')
    parser.add_argument('--db', type=str.lower,
                        choices=lib.DATABASES,
                        help='--db (mostly optional): The database used alongside the ARG annotation tool. This is\
                            only required if abricate or groot is used as a tool. Please refer here for more information\
                                on --db: https://github.com/BigDataBiology/argNorm?tab=readme-ov-file#--db-optional')
    parser.add_argument('-i', '--input', type=str,
                        help='-i (required): The path to the ARG annotation result which needs to be normalized.')
    parser.add_argument('--hamronization_skip_unsupported_tool', action='store_true',
                        help="--hamronization_skip_unsupported_tool (optional): skip rows with unsupported tools\
                            for hamronization outputs. argNorm be default will raise an exception if unsupported\
                             tool is found in hamronization. Use this if you only want argNorm to raise a warning.")
    parser.add_argument('-o', '--output', type=str,
                        help="-o (required): The path to the output file where you would like to store argNorm's results")
    args = parser.parse_args()

    if args.output is None:
        sys.stderr.write('Please specify an output file using `-o` or `--output`\n')
        sys.exit(2)

    if args.tool in ['groot', 'abricate'] and args.db == None:
        sys.stderr.write('Please specify a database using `--db` when using groot or abricate\n')
        sys.exit(2)

    # We only import the normalize function when the user actually wants to run the program
    # This makes running `argnorm -h` much faster because it avoids importing slow modules (e.g. pandas)
    from .normalize import normalize
    result = normalize(args.input,
            tool=args.tool,
            database=args.db,
            skip_on_unsupported_tool=args.hamronization_skip_unsupported_tool
        )

    prop_unmapped = result.ARO.isna().mean()
    print(f'{args.output}:', f'{1 - prop_unmapped:.1%} ARGs mapped.')
    with atomic_write(args.output, mode='w', overwrite=True) as out:
        out.write(f'# argNorm version: {__version__}\n')
        result.to_csv(out, sep='\t', index=False)

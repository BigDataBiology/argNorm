import argparse
import sys
from .argnorm_version import __version__
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
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tool', type=str.lower,
                        choices=['argsoap', 'abricate', 'deeparg', 'resfinder', 'amrfinderplus', 'groot', 'hamronization'],
                        help='The tool you used to do ARG annotation.')
    parser.add_argument('--db', type=str.lower,
                        choices=lib.DATABASES,
                        help='The database you used to do ARG annotation.')
    parser.add_argument('-i', '--input', type=str,
                        help='The annotation result you have')
    parser.add_argument('--hamronized', action='store_true',
                        help='Use hamronization as a tool instead')
    parser.add_argument('--hamronization_skip_unsupported_tool', action='store_true', help="Skip rows with unsupported tools for hamronization outputs")
    parser.add_argument('-o', '--output', type=str,
                        help='The file to save normalization results')
    args = parser.parse_args()

    if args.output is None:
        sys.stderr.write('Please specify an output file using `-o` or `--output`\n')
        sys.exit(2)

    if args.hamronized:
        sys.stderr.write('Upgrade to use hamronization as a tool instead of a flag\n')
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

    prop_unmapped = ((result.ARO == 'ARO:nan').sum() + result.ARO.isna().sum()) / result.shape[0]
    print(f'{args.output}:', f'{round(1 - prop_unmapped, 3):.2%} ARGs mapped.')
    with atomic_write(args.output, mode='w', overwrite=True) as out:
        out.write(f'# argNorm version: {__version__}\n')
        result.to_csv(out, sep='\t', index=False)

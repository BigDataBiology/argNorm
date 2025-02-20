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

    parser.add_argument('-o', '--output', type=str,
                        help='The file to save normalization results')
    args = parser.parse_args()

    if args.hamronized:
        sys.stderr.write('Upgrade to use hamronization as a tool instead of a flag\n')
        sys.exit(2)

    # We only import the normalize function when the user actually wants to run the program
    # This makes running `argnorm -h` much faster because it avoids importing slow modules (e.g. pandas)
    from .normalize import normalize
    result = normalize(args.input,
            tool=args.tool,
            database=args.db,
        )

    prop_unmapped = ((result.ARO == 'ARO:nan').sum() + result.ARO.isna().sum()) / result.shape[0]
    print(f'{args.output}:', f'{round(1 - prop_unmapped, 3):.2%} ARGs mapped.')
    with atomic_write(args.output, mode='w', overwrite=True) as out:
        out.write(f'# argNorm version: {__version__}\n')
        result.to_csv(out, sep='\t', index=False)

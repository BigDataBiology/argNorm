import argparse

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
                        choices=['argsoap', 'abricate', 'deeparg', 'resfinder', 'amrfinderplus'],
                        help='The tool you used to do ARG annotation.')
    parser.add_argument('--db', type=str.lower,
                        choices=['sarg', 'ncbi', 'resfinder', 'deeparg', 'megares', 'argannot', 'resfinderfg'],
                        help='The database you used to do ARG annotation.')
    parser.add_argument('--hamronized', action='store_true',
                        help='Use this if the input is hamronized (processed using the hAMRonization tool)')
    parser.add_argument('-i', '--input', type=str,
                        help='The annotation result you have')
    parser.add_argument('-o', '--output', type=str,
                        help='The file to save normalization results')
    args = parser.parse_args()

    # We only import the normalize function when the user actually wants to run the program
    # This makes running `argnorm -h` much faster because it avoids importing slow modules (e.g. pandas)
    from .normalize import normalize
    result = normalize(args.input,
            tool=args.tool,
            database=args.db,
            is_hamronized=args.hamronized
        )

    prop_unmapped = ((result.ARO == 'ARO:nan').sum() + result.ARO.isna().sum()) / result.shape[0]
    print(f'{args.output}:', f'{round(1 - prop_unmapped, 3):.2%} ARGs mapped.')
    result.to_csv(args.output, sep='\t', index=False)

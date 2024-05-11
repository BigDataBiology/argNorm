import os
import pandas as pd
try:
    pd.options.mode.copy_on_write = True
except pd.errors.OptionError:
    pass
from .drug_categorization import confers_resistance_to, drugs_to_drug_classes
from .lib import get_aro_mapping_table
from .lib import MAPPING_TABLE_ARO_COL, TARGET_ARO_COL

# Column headings for drug categorization output
CONFERS_RESISTANCE_TO_COL = 'confers_resistance_to'
RESISTANCE_TO_DRUG_CLASSES_COL = 'resistance_to_drug_classes'

class BaseNormalizer:
    """
    Inherit this class and customize subclass methods to implement the normalization of new databases/formats.
    """

    def __init__(self, database=None, is_hamronized=False) -> None:
        self.database = database
        self.is_hamronized = is_hamronized

    def run(self, input_file : str):
        """
        Main normalization pipeline.
        """
        original_annot = self.load_input(input_file)
        input_genes = self.get_input_ids(original_annot)

        aro_table = get_aro_mapping_table(self.database)
        aro_table.set_index(self.preprocess_ref_genes(
            aro_table.index
        ), inplace=True)
        mapping = aro_table[MAPPING_TABLE_ARO_COL].to_dict()
        original_annot[TARGET_ARO_COL] = input_genes.map(mapping)

        # Drug categorization
        original_annot[CONFERS_RESISTANCE_TO_COL] = original_annot[TARGET_ARO_COL].map(
                lambda a: ','.join(confers_resistance_to(a)),
                na_action='ignore')
        original_annot[RESISTANCE_TO_DRUG_CLASSES_COL] = original_annot[TARGET_ARO_COL].map(
                lambda a: ','.join(drugs_to_drug_classes(confers_resistance_to(a))),
                na_action='ignore')

        return original_annot

    def preprocess_ref_genes(self, ref_genes):
        """
        Customize this when ref gene and input gene can not exactly match.
        """
        return ref_genes


    def load_input(self, input_file):
        """
        Customize this when it fails to parse the input data.
        """
        return pd.read_csv(input_file, sep='\t')


class ARGSOAPNormalizer(BaseNormalizer):
    def __init__(self, database=None, is_hamronized=False) -> None:
        database = 'sarg'
        super().__init__(database, is_hamronized)


    def get_input_ids(self, itable):
        return itable['reference_accession' if self.is_hamronized else 1]

    def load_input(self, input_file):
        if self.is_hamronized:
            return pd.read_csv(input_file, sep='\t')
        else:
            return pd.read_csv(input_file, sep='\t', header=None)


class DeepARGNormalizer(BaseNormalizer):
    def __init__(self, database=None, is_hamronized=False) -> None:
        database = 'deeparg'
        super().__init__(database, is_hamronized)

    def get_input_ids(self, itable):
        return itable['gene_name' if self.is_hamronized else 'best-hit']

class ResFinderNormalizer(BaseNormalizer):
    def __init__(self, database=None, is_hamronized=False) -> None:
        database = 'resfinder'
        super().__init__(database, is_hamronized)

    def get_input_ids(self, itable):
        return itable['gene_symbol' if self.is_hamronized else 'Accession no.']

    def preprocess_ref_genes(self, ref_genes):
        return ref_genes.str.split('_').str[0 if self.is_hamronized else -1]


class AMRFinderPlusNormalizer(BaseNormalizer):
    def __init__(self, database=None, is_hamronized=False) -> None:
        database = 'ncbi'
        super().__init__(database, is_hamronized)

    def get_input_ids(self, itable):
        return itable['gene_symbol' if self.is_hamronized else 'Accession of closest sequence']

    def preprocess_ref_genes(self, ref_genes):
        return ref_genes.str.split('|').str[5 if self.is_hamronized else 1]


class AbricateNormalizer(BaseNormalizer):
    def __init__(self, database=None, is_hamronized=False) -> None:
        if database not in ['ncbi', 'deeparg', 'resfinder', 'sarg', 'megares', 'argannot', 'resfinderfg']:
            raise Exception(f'{database} is not a supported database.')
        if not is_hamronized and database in ['sarg', 'resfinderfg']:
            raise Exception(f'{database} is not a supported database for raw files.')

        super().__init__(database, is_hamronized)

    def get_input_ids(self, itable):
        if self.is_hamronized:
            col = dict(
                ncbi='gene_symbol',
                deeparg='gene_name',
                resfinder='gene_symbol',
                sarg='gene_symbol',
                megares='reference_accession',
                argannot='reference_accession',
                resfinderfg='gene_name'
            )[self.database]
        else:
            col = dict(
                ncbi='GENE',
                deeparg='best-hit',
                resfinder='GENE',
                megares='ACCESSION',
                argannot='ACCESSION'
            )[self.database]
        if self.database == 'resfinderfg':
            return itable[col].str.split('|').str[1]
        return itable[col]


    def preprocess_argannot_ref_genes(self, ref_gene):
        split_str = ref_gene.split(':')
        if not str(split_str[2][0]).isnumeric() and not '-' in split_str[2]:
            return ':'.join([split_str[1], split_str[3]])
        return ':'.join(ref_gene.split(':')[1:3])

    def preprocess_ref_genes(self, ref_genes):
        process_funcs_by_db = dict(
            ncbi=lambda x: x.split('|')[5],
            deeparg=lambda x: x,
            resfinder=lambda x: '_'.join(x.split('_')[:-1]),
            sarg=lambda x: x,
            megares=lambda x: x.split('|')[0],
            argannot=self.preprocess_argannot_ref_genes,
            resfinderfg=lambda x: x.split('|')[1]
        )
        return ref_genes.map(process_funcs_by_db[self.database])


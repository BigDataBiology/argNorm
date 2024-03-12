import os
import pandas as pd
import warnings

from .drug_categorization import confers_resistance_to, drugs_to_drug_classes

ORIGINAL_ID_COL = 'Original ID'
MAPPING_TABLE_ARO_COL = 'ARO'
TARGET_ARO_COL = 'ARO'

# Column headings for drug categorization output
CONFERS_RESISTANCE_TO_COL = 'confers_resistance_to'
RESISTANCE_TO_DRUG_CLASSES_COL = 'resistance_to_drug_classes'

_ROOT = os.path.abspath(os.path.dirname(__file__))

def is_number(num):
    """
    Required for checking aro mappings to discern between numbers and other
    string identifiers.
    """
    try:
        float(num)
    except ValueError:
        return False

    return True

def get_data_path(path, getting_manual_curation):
    """
    Gets mapping tables and manual curation tables.
    Giving 'True' as argument after 'path' will get manual curation table.
    Else mapping table will be returned.
    """
    if getting_manual_curation:
        return os.path.join(_ROOT, 'data/manual_curation', path)

    return os.path.join(_ROOT, 'data', path)

class BaseNormalizer:
    """
    Inherit this class and customize subclass methods to implement the normalization of tools.
    """

    def __init__(self, database=None, is_hamronized=False, mode=None, uses_manual_curation=True) -> None:
        self.tool = ''
        self.database = database
        self.mode = mode
        self.is_hamronized = is_hamronized
        self.uses_manual_curation = uses_manual_curation
        self._set_input_gene_col()

    def run(self, input_file : str):
        """
        Main normalization pipeline.
        """
        original_annot = self.load_input(input_file)
        input_genes = self.preprocess_input_genes(
            original_annot[self._input_gene_col].str.lower()
        )
        aro_table = self.get_aro_mapping_table()
        aro_table.set_index(self.preprocess_ref_genes(
            aro_table[ORIGINAL_ID_COL].str.lower()
        ), inplace=True)
        mapping = aro_table[MAPPING_TABLE_ARO_COL].to_dict()
        original_annot[TARGET_ARO_COL] = input_genes.map(mapping)

        # Drug categorization
        original_annot[CONFERS_RESISTANCE_TO_COL] = self.initial_drug_categorization(original_annot[TARGET_ARO_COL])
        original_annot[RESISTANCE_TO_DRUG_CLASSES_COL] = self.final_drug_categorization(original_annot[TARGET_ARO_COL])

        return original_annot

    def initial_drug_categorization(self, aro_list):
        result = []
        for aro in aro_list:
            result.append(", ".join(confers_resistance_to(aro)))

        return result

    def final_drug_categorization(self, aro_list):
        drugs = []
        for aro in aro_list:
            drugs.append(confers_resistance_to(aro))

        drug_classes = []
        for drug in drugs:
            drug_classes.append(", ".join(drugs_to_drug_classes(drug)))

        return drug_classes

    def preprocess_ref_genes(self, ref_genes):
        """
        Customize this when ref gene and input gene can not exactly match.
        """
        return ref_genes

    def preprocess_input_genes(self, input_genes):
        """
        Customize this when ref gene and input gene can not exactly match.
        """
        return input_genes


    def _set_input_gene_col(self):
        """
        Always adapt this method to the input data format.
        """
        self._input_gene_col = ''

    def get_aro_mapping_table(self):
        """
        Don't customize this unless you're using your own (not package built-in) reference data.
        """
        df = pd.read_csv(get_data_path(f'{self.tool}_{self.database}_{self.mode}_ARO_mapping.tsv', False), sep='\t', index_col=0)

        if self.uses_manual_curation:
            gene_identifier = 'Original ID'

            if self.database == 'ncbi':
                manual_curation_fname = 'ncbi_manual_curation.tsv'
            elif self.database == 'resfinder':
                manual_curation_fname = 'resfinder_manual_curation.tsv'
            else:
                manual_curation_fname = f'{self.tool}_{self.database}_{self.mode}_manual_curation.tsv'
            manual_curation = pd.read_csv(get_data_path(manual_curation_fname, True), sep='\t')

            aro_nan_indices = [(list(df[gene_identifier]).index(manual_curation.loc[i, gene_identifier])) for i in range(manual_curation.shape[0])]

            for i in range(len(aro_nan_indices)):
                df.loc[aro_nan_indices[i], 'ARO'] = manual_curation.loc[i, 'ARO Replacement']

                if self.tool != 'argsoap' and self.mode != 'orfs':
                    df.loc[aro_nan_indices[i], 'Gene Name in CARD'] = manual_curation.loc[i, 'Gene Name in CARD']
            if self.tool != 'argsoap' or self.mode != 'orfs':
                df[TARGET_ARO_COL] = df[TARGET_ARO_COL].map(lambda a: f'ARO:{int(float(a)) if is_number(a) == True else a}')
        else:
            if self.tool != 'argsoap' or self.mode != 'orfs':
                df[TARGET_ARO_COL] = df[TARGET_ARO_COL].map(lambda a: f'ARO:{int(a) if a == a else "nan"}') # a == a checks that a is not nan

        return df

    def load_input(self, input_file):
        """
        Customize this when it fails to parse the input data.
        """
        return pd.read_csv(input_file, sep='\t')


class ARGSOAPNormalizer(BaseNormalizer):
    def __init__(self, database=None, is_hamronized=False, mode=None, uses_manual_curation=True) -> None:
        if not database:
            warnings.warn('No `database` specified. Will try using SARG.')
            database = 'sarg'
        elif database != 'sarg':
            warnings.warn('The `database` is not supported. Will try using SARG instead.')
            database = 'sarg'
        super().__init__(database, is_hamronized, mode, uses_manual_curation)
        self.tool = 'argsoap'


    def _set_input_gene_col(self):
        if self.is_hamronized and self.mode == 'reads':
            self._input_gene_col = 'reference_accession'
        elif self.is_hamronized and self.mode == 'orfs':
            self._input_gene_col = 'gene_name'
        elif not self.is_hamronized and self.mode == 'reads':
            self._input_gene_col = 1
        elif not self.is_hamronized and self.mode == 'orfs':
            self._input_gene_col = 0
        else:
            self._raise_incorrect_mode_error()

    def load_input(self, input_file):
        if self.is_hamronized:
            return pd.read_csv(input_file, sep='\t')
        else:
            return pd.read_csv(input_file, sep='\t', header=None)

    def preprocess_input_genes(self, input_genes):
        if self.is_hamronized and self.mode == 'reads':
            return input_genes
        elif self.is_hamronized and self.mode == 'orfs':
            return input_genes
        elif not self.is_hamronized and self.mode == 'reads':
            return input_genes
        elif not self.is_hamronized and self.mode == 'orfs':
            return input_genes.apply(lambda x: x.split('_train_msa')[0])
        else:
            self._raise_incorrect_mode_error()

    def preprocess_ref_genes(self, ref_genes):
        if self.is_hamronized and self.mode == 'reads':
            return ref_genes
        elif self.is_hamronized and self.mode == 'orfs':
            return ref_genes.apply(lambda x: x.replace("'", '_').replace('-', '_'))
        elif not self.is_hamronized and self.mode == 'reads':
            return ref_genes
        elif not self.is_hamronized and self.mode == 'orfs':
            return ref_genes.apply(lambda x: x.replace("'", '_').replace('-', '_'))
        else:
            self._raise_incorrect_mode_error()

    def _raise_incorrect_mode_error(self):
        raise ValueError('Please specify correct mode for your input.')


class DeepARGNormalizer(BaseNormalizer):

    def __init__(self, database=None, is_hamronized=False, mode=None, uses_manual_curation=True) -> None:
        if mode:
            warnings.warn('`mode` is not relavant for DeepARG and will be ignored.')
            mode = 'both'
        else:
            warnings.warn('`mode` is not specified. Will use default setting "both".')
            mode = 'both'
        if not database:
            warnings.warn('No `database` specified. Will try using DeepARG.')
            database = 'deeparg'
        elif database != 'deeparg':
            warnings.warn('The `database` is not supported. Will try using DeepARG instead.')
            database = 'deeparg'
        super().__init__(database, is_hamronized, mode, uses_manual_curation)
        self.tool = 'deeparg'

    def _set_input_gene_col(self):
        """
        Always adapt this method to the input data format.
        """
        if self.is_hamronized:
            self._input_gene_col = 'gene_name'
        else:
            self._input_gene_col = 'best-hit'


class ResFinderNormalizer(BaseNormalizer):

    def __init__(self, database=None, is_hamronized=False, mode=None, uses_manual_curation=True) -> None:
        if mode:
            warnings.warn('`mode` is not relavant for ResFinder and will be ignored.')
            mode = 'both'
        else:
            warnings.warn('`mode` is not specified. Will use default setting "both".')
            mode = 'both'
        if not database:
            warnings.warn('No `database` specified. Will try using ResFinder.')
            database = 'resfinder'
        elif database != 'resfinder':
            warnings.warn('The `database` is not supported. Will try using ResFinder instead.')
            database = 'resfinder'
        super().__init__(database, is_hamronized, mode, uses_manual_curation)
        self.tool = 'resfinder'

    def _set_input_gene_col(self):
        """
        Always adapt this method to the input data format.
        """
        if self.is_hamronized:
            self._input_gene_col = 'gene_symbol'
        else:
            self._input_gene_col = 'Accession no.'

    def preprocess_ref_genes(self, ref_genes):
        if self.is_hamronized:
            return ref_genes.apply(lambda x: x.split('_')[0])
        else:
            return ref_genes.apply(lambda x: x.split('_')[-1])


class AMRFinderPlusNormalizer(BaseNormalizer):

    def __init__(self, database=None, is_hamronized=False, mode=None, uses_manual_curation=True) -> None:
        if mode:
            warnings.warn('`mode` is not relavant for AMRFinderPlus and will be ignored.')
            mode = 'both'
        else:
            warnings.warn('`mode` is not specified. Will use default setting "both".')
            mode = 'both'
        if not database:
            warnings.warn('No `database` specified. Will try using NCBI.')
            database = 'ncbi'
        elif database != 'ncbi':
            warnings.warn('The `database` is not supported. Will try using NCBI instead.')
            database = 'ncbi'
        super().__init__(database, is_hamronized, mode, uses_manual_curation)
        self.tool = 'amrfinderplus'

    def _set_input_gene_col(self):
        """
        Always adapt this method to the input data format.
        """
        if self.is_hamronized:
            self._input_gene_col = 'gene_symbol'
        else:
            self._input_gene_col = 'Accession of closest sequence'

    def preprocess_ref_genes(self, ref_genes):
        if self.is_hamronized:
            return ref_genes.apply(lambda x: x.split('|')[5])
        else:
            return ref_genes.apply(lambda x: x.split('|')[1])


class AbricateNormalizer(BaseNormalizer):

    def __init__(self, database=None, is_hamronized=False, mode=None, uses_manual_curation=True) -> None:
        if mode:
            warnings.warn('`mode` is not relavant for Abricate and will be ignored.')
            mode = 'both'
        else:
            warnings.warn('`mode` is not specified. Will use default setting "both".')
            mode = 'both'
        super().__init__(database, is_hamronized, mode, uses_manual_curation)
        self.tool = 'abricate'

    def _set_input_gene_col(self):
        if self.is_hamronized:
            gene_col_by_db = dict(
            ncbi='gene_symbol',
            deeparg='gene_name',
            resfinder='gene_symbol',
            sarg='gene_symbol',
            megares='reference_accession',
            argannot='reference_accession'
        )
        else:
            gene_col_by_db = dict(
            ncbi='GENE',
            deeparg='best-hit',
            resfinder='GENE',
            megares='ACCESSION',
            argannot='ACCESSION'
        )
        self._input_gene_col = gene_col_by_db[self.database]

    def preprocess_input_genes(self, input_genes):
        process_funcs_by_db = dict(
            ncbi=lambda x: x,
            deeparg=lambda x: x,
            resfinder=lambda x: x,
            sarg=lambda x: x,
            megares=lambda x: x,
            argannot=lambda x: x
        )
        return input_genes.apply(process_funcs_by_db[self.database])

    def preprocess_ref_genes(self, ref_genes):
        process_funcs_by_db = dict(
            ncbi=lambda x: x.split('|')[5],
            deeparg=lambda x: x,
            resfinder=lambda x: '_'.join(x.split('_')[:-1]),
            sarg=lambda x: x,
            megares=lambda x: x.split('|')[0],
            argannot=lambda x: x.split('~~~')[-1]
        )
        return ref_genes.apply(process_funcs_by_db[self.database])


"""Source code for the Normalizer class"""

import os
import pandas as pd
import warnings

MAPPING_TABLE_ARO_COL = 'ARO'
TARGET_ARO_COL = 'ARO'

_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'data', path)


class BaseNormalizer:
    """
    Inherit this class and customize subclass methods to implement the normalization of tools.
    """

    def __init__(self, database=None, is_hamronized=False, mode=None) -> None:
        self.tool = ''
        self.database = database
        self.mode = mode
        self.is_hamronized = is_hamronized
        self._set_input_gene_col()
        self._set_ref_gene_and_aro_cols()

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
            aro_table[self.ref_gene_col].str.lower()
        ), inplace=True)
        mapping = aro_table[MAPPING_TABLE_ARO_COL].to_dict()
        original_annot[TARGET_ARO_COL] = input_genes.map(mapping)
        return original_annot


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


    def _set_ref_gene_and_aro_cols(self):
        """
        Customize this when the reference data format is different from the default (e.g. for sarg orfs mode).
        """
        self.ref_gene_col = 'Original ID'

    def _set_input_gene_col(self):
        """
        Always adapt this method to the input data format.
        """
        self._input_gene_col = ''

    def get_aro_mapping_table(self):
        """
        Don't customize this unless you're using your own (not package built-in) reference data.
        """
        df = pd.read_csv(get_data(f'{self.tool}_{self.database}_{self.mode}_ARO_mapping.tsv'), sep='\t', index_col=0)
        if self.tool != 'argsoap' or self.mode != 'orfs':
            df[TARGET_ARO_COL] = 'ARO:' + df[TARGET_ARO_COL].astype(str).apply(lambda x: x.split('.')[0])
        return df

    def load_input(self, input_file):
        """
        Customize this when it fails to parse the input data.
        """
        return pd.read_csv(input_file, sep='\t')


    
class ARGSOAPNormalizer(BaseNormalizer):
    def __init__(self, database=None, is_hamronized=False, mode=None) -> None:
        if not database:
            warnings.warn('No `database` specified. Will try using SARG.')
            database = 'sarg'
        elif database != 'sarg':
            warnings.warn('The `database` is not supported. Will try using SARG instead.')
            database = 'sarg'
        super().__init__(database, is_hamronized, mode)
        self.tool = 'argsoap'
    
    def _set_ref_gene_and_aro_cols(self):
        if self.mode == 'reads':
            self.ref_gene_col = 'Original ID'
        elif self.mode == 'orfs':
            self.ref_gene_col = 'Categories_in_database'
        else:
            self._raise_incorrect_mode_error()

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

    def __init__(self, database=None, is_hamronized=False, mode=None) -> None:
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
        super().__init__(database, is_hamronized, mode)
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

    def __init__(self, database=None, is_hamronized=False, mode=None) -> None:
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
        super().__init__(database, is_hamronized, mode)
        self.tool = 'resfinder'
    
    def _set_input_gene_col(self):
        """
        Always adapt this method to the input data format.
        """
        if self.is_hamronized:
            self._input_gene_col = ''  # TODO add here.
        else:
            self._input_gene_col = 'Accession no.'

    def preprocess_ref_genes(self, ref_genes):
        return ref_genes.apply(lambda x: x.split('_')[-1])


class AMRFinderPlusNormalizer(BaseNormalizer):

    def __init__(self, database=None, is_hamronized=False, mode=None) -> None:
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
        super().__init__(database, is_hamronized, mode)
        self.tool = 'amrfinderplus'
    
    def _set_input_gene_col(self):
        """
        Always adapt this method to the input data format.
        """
        if self.is_hamronized:
            self._input_gene_col = ''  # TODO add this.
        else:
            self._input_gene_col = 'Accession of closest sequence'

    def preprocess_ref_genes(self, ref_genes):
        return ref_genes.apply(lambda x: x.split('|')[1])

    


class AbricateNormalizer(BaseNormalizer):

    def __init__(self, database=None, is_hamronized=False, mode=None) -> None:
        if mode:
            warnings.warn('`mode` is not relavant for Abricate and will be ignored.')
            mode = 'both'
        else:
            warnings.warn('`mode` is not specified. Will use default setting "both".')
            mode = 'both'
        super().__init__(database, is_hamronized, mode)
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

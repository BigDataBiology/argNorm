"""Source code for the Normalizer class"""

import os
import pandas as pd


_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'data', path)

# supported_dbs = ['ncbi', 'resfinder', 'sarg', 'deeparg', 'megares', 'argannot']


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
        self._set_tmp_gene_col()
        self._set_output_aro_col()
        self._set_ref_gene_and_aro_cols()

    def run(self, input_file):
        """
        Main normalization pipeline.
        """
        original_annot = self.load_input(input_file)
        original_annot[self._tmp_gene_col] = self.preprocess_input_genes(
            original_annot[self._input_gene_col].str.lower()
        )
        aro_table = self.get_aro_mapping_table()
        ref_genes = self.preprocess_ref_genes(
            aro_table[self.ref_gene_col].str.lower()
        )
        mapping = self._make_hashmap(ref_genes, aro_table[self.ref_aro_col])
        original_annot[self._aro_col] = original_annot[self._tmp_gene_col].map(mapping)
        return original_annot.drop(columns=self._tmp_gene_col)
        

    def preprocess_ref_genes(self, ref_genes):
        """
        Customize this
        """ 
        return ref_genes

    def preprocess_input_genes(self, input_genes):
        """
        Customize this
        """ 
        return input_genes

    def _set_tmp_gene_col(self):
        """
        Customize this when the default col name causes conflicted with the input
        """   
        self._tmp_gene_col = 'tmp_gene'

    def _set_output_aro_col(self):   
        """
        Customize this when the default col name causes conflicted with the input
        """   
        self._aro_col = 'ARO'

    def _set_ref_gene_and_aro_cols(self):
        """
        Customize this when the reference data format is different from the default (e.g. for sarg orfs mode).
        """
        self.ref_gene_col = 'Original ID'
        self.ref_aro_col = 'ARO'

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
            df[self._aro_col] = 'ARO:' + df[self._aro_col].astype(str).apply(lambda x: x.split('.')[0])
        return df

    def load_input(self, input_file):
        """
        Customize this when it fails to parse the input data.
        """
        return pd.read_csv(input_file, sep='\t')

    def _make_hashmap(self, keys, vals):
        return dict(zip(keys, vals))

    
class ARGSOAPNormalizer(BaseNormalizer):
    def __init__(self, database=None, is_hamronized=False, mode=None) -> None:
        super().__init__(database, is_hamronized, mode)
        self.tool = 'argsoap'
        self.database = 'sarg'
    
    def _set_ref_gene_and_aro_cols(self):
        if self.mode == 'reads':
            self.ref_gene_col = 'Original ID'
        elif self.mode == 'orfs':
            self.ref_gene_col = 'Categories_in_database'
        else:
            self._raise_incorrect_mode_error()
        self.ref_aro_col = 'ARO'

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
        super().__init__(database, is_hamronized, mode)
        self.tool = 'deeparg'


class AbricateNormalizer(BaseNormalizer):

    def __init__(self, database=None, is_hamronized=False, mode=None) -> None:
        super().__init__(database, is_hamronized, mode)
        self.tool = 'abricate'


class HamronizedNormalizer:
    """Major class to perform the normalization"""

    def __init__(self, db):
        self.database = db
        self._set_mapping_table()
        self._aro_col = 'ARO'
        self._set_gene_name_cols()
        self._set_input_gene_name_process_funcs()
        self._set_AROmap_gene_name_process_funcs()
        self._set_input_loading_funcs()
        self._temp_gene_name_col = 'temp_gene_name'

    def _set_mapping_table(self):  ## TODO update this.
        suffix = '_ARO_mapping.tsv'
        self.mapping_tables = dict(
            ncbi=get_data('ncbi' + suffix),
            deeparg=get_data('deeparg' + suffix),
            resfinder=get_data('resfinder' + suffix),
            sarg=get_data('sarg' + suffix),
            megares=get_data('megares' + suffix),
            argannot=get_data('argannot' + suffix)
        )

    def _set_gene_name_cols(self):
        self._gene_name_cols = dict(
            ncbi='gene_symbol',
            deeparg='gene_name',
            resfinder='gene_symbol',
            sarg='gene_symbol',
            megares='reference_accession',
            argannot='reference_accession'
        )
        
    def _set_input_gene_name_process_funcs(self):
        self._preprocess_input_gene_name_funcs = dict(
            ncbi=lambda x: x,
            deeparg=lambda x: x,
            resfinder=lambda x: x,
            sarg=lambda x: x,
            megares=lambda x: x,
            argannot=lambda x: x
        )

    def _set_AROmap_gene_name_process_funcs(self):
        self._preprocess_AROmap_gene_name_funcs = dict(
            ncbi=lambda x: x.split('|')[5],
            deeparg=lambda x: x,
            resfinder=lambda x: '_'.join(x.split('_')[:-1]),
            sarg=lambda x: x,
            megares=lambda x: x.split('|')[0],
            argannot=lambda x: x.split('~~~')[-1]
        )

    def _set_input_loading_funcs(self):
        fun = lambda file: pd.read_csv(file, sep='\t', index_col=0)
        self._input_loading_funcs = dict(
            ncbi=fun,
            deeparg=fun,
            resfinder=fun,
            sarg=fun,
            megares=fun,
            argannot=fun
        )

    def run(self, input_file):
        """Performs the normalization"""
        original_annot = self._load_input(input_file)
        original_annot = self._preprocess_input(original_annot)
        ##
        mapping_table = self.load_mapping_table()
        names = self._preprocess_gene_names(mapping_table['Original ID']).str.lower()
        aros = 'ARO:' + mapping_table[self._aro_col].astype(str).apply(lambda x: x.split('.')[0])
        hashmap = self._make_hashmap(names, aros)
        ##
        original_annot[self._aro_col] = original_annot[self._temp_gene_name_col].\
            str.lower().map(hashmap)
        return original_annot.drop(columns=[self._temp_gene_name_col])

    def _preprocess_input(self, input_annot):
        gene_name_col = self._gene_name_cols[self.database]
        single_gene_name_process = self._preprocess_input_gene_name_funcs[self.database]
        input_annot[self._temp_gene_name_col] = input_annot[gene_name_col].\
            apply(single_gene_name_process)
        return input_annot

    def _preprocess_gene_names(self, gene_names):
        return gene_names.apply(self._preprocess_AROmap_gene_name_funcs[self.database])

    def load_mapping_table(self):
        file = self.mapping_tables[self.database]
        return pd.read_csv(file, sep='\t', index_col=0)

    def _load_input(self, file):
        return self._input_loading_funcs[self.database](file)

    def _make_hashmap(self, keys, vals):
        return dict(zip(keys, vals))


class RawNormalizer(HamronizedNormalizer):

    def __init__(self, db):
        HamronizedNormalizer.__init__(self, db=db)

    def _set_gene_name_cols(self):
        self._gene_name_cols = dict(
            ncbi='GENE',
            deeparg='best-hit',
            resfinder='GENE',
            megares='ACCESSION',
            argannot='ACCESSION'
        )

    def _set_input_loading_funcs(self):
        HamronizedNormalizer._set_input_loading_funcs(self)
        self._input_loading_funcs['sarg'] = lambda x: x # TODO to be updated

    def run(self, input_file):
        return HamronizedNormalizer.run(self, input_file=input_file)

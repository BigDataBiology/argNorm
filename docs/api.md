# argNorm API

## argnorm.lib: core normalization functionality

### argnorm.lib.DATABASES

A list of supported databases.

### argnorm.lib.map_to_aro(): maps a given gene in a specific database to an ARO term.

#### Parameters
* gene (str): The original ID of the gene as mentioned in source database.
* database (str): name of database. Can be: argannot, deeparg, megares, ncbi, resfinderfg, sarg, groot-db, groot-core-db, groot-argannot, groot-resfinder, and groot-card

#### Returns
* pronto.term.Term: A pronto term with the ARO number of input gene. ARO number can be accessed using 'id' attribute and gene name can be accessed using 'name' attribute. 

> Note: if ARO mapping is doesn't exist, None is returned.

#### Example

```
from argnorm.lib import map_to_aro

# Mapping the `ARR-2_1_HQ141279` gene from the `resfinder` database to the ARO
print(map_to_aro('ARR-2_1_HQ141279', 'resfinder'))

# Mapping the `argannot~~~(Bla)cfxA4~~~AY769933:1-966` gene in `groot` using the `groot-argannot` database
print(map_to_aro('argannot~~~(Bla)cfxA4~~~AY769933:1-966', 'groot-argannot'))
```

### argnorm.lib.get_aro_mapping_table(): gets ARO mapping table for a specific database

#### Parameters 
* database (str): name of database. Can be: argannot, deeparg, megares, ncbi, resfinderfg, sarg, groot-db, groot-core-db, groot-argannot, groot-resfinder, and groot-card

#### Returns
* pandas.DataFrame: A pandas dataframe with ARGs mapped to AROs.

#### Example
```
# Getting the ARO mapping table for the `resfinder` database
from argnorm.lib import get_aro_mapping_table
print(get_aro_mapping_table('resfinder'))
```

## argnorm.drug_categorization: drug categorization functions

### argnorm.drug_categorization.confers_resistance_to(): returns a list of the drugs/antibiotics to which a gene confers resistance to.

#### Parameters
* aro_num (str): ARO number. Needs to be in the form 'ARO:number'.

#### Returns
* list[str]: A list with ARO number of the drugs/antibiotics to which the input gene confers resistance to.

#### Example

```
# Getting a list of all the drugs to which tet(X)(ARO:3000205) confers resistance to
from argnorm.drug_categorization import confers_resistance_to
print(confers_resistance_to('ARO:3000205'))
```

### argnorm.drug_categorization.drugs_to_drug_classes(): returns a list of categories of drug classes
> Example: cephem and penam are categorized as beta_lactam antibiotics.

#### Parameters
* list[str]: A list with ARO numbers of the drugs/antibiotics to which a gene confers resistance to.

#### Returns
* (list[str]): A list containing the ARO numbers of the drug classes of each drug given as input to the function. Order of the ARO numbers of the drug classes corresponds to the order in which the drugs were given to the function in the drugs_list.

#### Example

```
# Getting a list of the drug categories of the drugs to which tet(X) (ARO:3000205) confers resistance to
from argnorm.drug_categorization import drugs_to_drug_classes, confers_resistance_to
print(drugs_to_drug_classes(confers_resistance_to('ARO:3000205')))
print(drugs_to_drug_classes(['ARO:0000030', 'ARO:0000051', 'ARO:0000069', 'ARO:3000152', 'ARO:3000528', 'ARO:3000667', 'ARO:3000668']))
```

## argnorm.normalizers: normalizers for specific tools

Normalizers classes for specific tools which normalize ARG annotation outputs. Same functionality as CLI.

All normalizers have 1 optional parameter:
* database (str): name of database. Can be: argannot, deeparg, megares, ncbi, resfinderfg, sarg, groot-db, groot-core-db, groot-card, groot-argannot, and groot-resfinder.

> Note: the database parameter only needs to be specified for AbricateNormalizer and GrootNormalizer. ncbi, deeparg, resfinder, sarg, megares, argannot, resfinderfg are the supported databases for AbricateNormalizer and groot-db, groot-core-db, groot-argannot, groot-resfinder, and groot-card are the supported databases for GrootNormalizer.

Available normalizers:
* argnorm.normalizers.ARGSOAPNormalizer
* argnorm.normalizers.DeepARGNormalizer
* argnorm.normalizers.ResFinderNormalizer
* argnorm.normalizers.AMRFinderPlusNormalizer
* argnorm.normalizers.AbricateNormalizer
* argnorm.normalizers.GrootNormalizer
* argnorm.normalizers.HamronizationNormalizer

### Methods

#### run(): return a pandas.DataFrame with the ARO and drug categorization columns added to the initial ARG annotation output dataframe.

Parameters: 
* input_file (str): path to ARG annotation output on which normalization & drug categorization is to be performed.

Returns: 
* pandas.DataFrame: A pandas dataframe with with ARO and drug categorization columns added to dataframe from `input_file`.

### Example 1: using ResFinderNormalizer
Download the sample data from [here](https://raw.githubusercontent.com/BigDataBiology/argNorm/main/examples/raw/resfinder.resfinder.orfs.tsv) in a folder called `argnorm_normalizers_tutorial`.

```
mkdir argnorm_normalizers_tutorial
cd argnorm_normalizers_tutorial
wget https://raw.githubusercontent.com/BigDataBiology/argNorm/main/examples/raw/resfinder.resfinder.orfs.tsv
```

Save the following piece of Python code in the `argnorm_normalizers_tutorial` folder, and run the script.

```
from argnorm.normalizers import ResFinderNormalizer

resfinder_normalizer = ResFinderNormalizer(is_hamronized=False)
resfinder_normalizer.run('./resfinder.resfinder.orfs.tsv').to_csv('./resfinder.resfinder.orfs.normed.tsv', sep='\t')
```

This will create a file called `resfinder.resfinder.orfs.normed.tsv` with ARO mappings and drug categorization.

### Example 2: using HamronizationNormalizer

Download the sample data [here](https://raw.githubusercontent.com/BigDataBiology/argNorm/7ee9d74c9fa51956ecb7706fa979cc0696ae305d/examples/hamronized/abricate.resfinderfg.tsv), and store it in a folder called `argnorm_normalizers_tutorial`.

```
mkdir argnorm_normalizers_tutorial
cd argnorm_normalizers_tutorial
wget https://raw.githubusercontent.com/BigDataBiology/argNorm/7ee9d74c9fa51956ecb7706fa979cc0696ae305d/examples/hamronized/abricate.resfinderfg.tsv
```

Save the following piece of Python code in the `argnorm_normalizers_tutorial` folder, and run the script.

```
from argnorm.normalizers import HamronizationNormalizer

normalizer = HamronizationNormalizer()
normalizer.run('./abricate.resfinderfg.tsv').to_csv('./abricate.resfinderfg.normed.tsv', sep='\t')
```

This will create a file called `abricate.resfinderfg.normed.tsv` with ARO mappings and drug categorization.
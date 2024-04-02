# argNorm: Normalize ARG Annotations to the ARO

[![Python package](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml/badge.svg)](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml)
[![Downloads](https://pepy.tech/badge/argNorm)](https://pepy.tech/project/argNorm)
![](https://img.shields.io/badge/status-alpha-red?style=flat)

## What is argNorm?

argNorm is a tool to normalize antiobiotic resistance genes (ARGs) by mapping them to the [antibiotic resistance ontology (ARO)](https://obofoundry.org/ontology/aro.html) created by the [CARD database](https://card.mcmaster.ca/).

argNorm also enhances antibiotic resistance gene annotations by providing drug categorization of the drugs that antibiotic resistance genes confer resistance to.

> #### *Note:*
> *This is a beta-quality implementation (subject to changes and some bugs may remain), but you're welcome to try it and provide feedback in the [Issues Page](https://github.com/BigDataBiology/argNorm/issues).*

## Why argNorm?

Right now, many tools exist for annotating ARGs in genomes and metagenomes. However, each tool will have distinct output formats.

The [hAMRonization package](https://github.com/pha4ge/hAMRonization) can normalize file formats, but each tool will use different names/identifiers (_e.g._, `TetA` or `TETA` or `tet(A)` or `tet-A` are all different ways to spell the same gene name).

For a small number of isolate genomes, a human user can manually evaluate the outputs.
However, in metagenomics, especially for large-scale projects, this becomes infeasible.
Thus, `argNorm` normalizes the _output vocabulary_ of ARG annotation tools by mapping them to the same ontology (ARO).

## Drug Categorization

Besides performing normalization, argNorm also provides drug categorization of drugs that antibiotic resistance genes confer resistance to.

For example, the `PBP2b` (`ARO:3003042`) gene confers resistance to the drug class `amoxicillin`. `amoxicillin` is then categorized into a broader category of `beta lactam antibiotic`.

argNorm provides support for this, and adds the `confers_resistance_to` and `resistance_to_drug_classes` columns to ARG annotations.

The `confers_resistance_to` column will contain ARO numbers of all the drug classes that a gene provides resistance to (`ARO:0000064` for `amoxicillin` in the previous example).

The `resistance_to_drug_classes` column will contain ARO numbers of the broader categories of the drug classes in the `confers_resistance_to` column (`ARO:3000007` for `beta lactam antibiotic` in the previous example).

![argNorm Workflow](./images/argnorm_workflow.svg)

## Supported Tools

- [DeepARG](https://bench.cs.vt.edu/deeparg) (v1.0.2)
- [ARGs-OAP](https://galaxyproject.org/use/args-oap/) (v3)
- [ABRicate](https://github.com/tseemann/abricate) (v1.0.1) with NCBI (v3.6), ResFinder (v4.1.11), MEGARes (v2.0), ARG-ANNOT (v5), ResFinderFG (v2)
- [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) (v4.0)
- [AMRFinderPlus](https://github.com/ncbi/amr) (v3.10.30)

## Installation

argNorm can be installed using pip using any one of the commands shown:

```
pip install argnorm

pip install https://github.com/BigDataBiology/argNorm/archive/refs/heads/main.zip
```

## Tutorial Video

[![argNorm Tutorial](https://markdown-videos-api.jorgenkh.no/url?url=https%3A%2F%2Fyoutu.be%2Fvx8MCQ7gDLs)](https://youtu.be/vx8MCQ7gDLs)

## Quick Example 1: argNorm as a Command Line Tool

Here is a quick demo of running argNorm on the command line.

### Step 0: Install argNorm & Create Working Directory

```
pip install argnorm

argnorm -h
```

`argnorm -h` or `argnorm --help` will display all the available options to run argNorm with.

```
> argnorm -h

usage: argnorm [-h] [--db {sarg,ncbi,resfinder,deeparg,megares,argannot}] [--hamronized] [-i INPUT] [-o OUTPUT] {argsoap,abricate,deeparg,resfinder,amrfinderplus}

argNorm normalizes ARG annotation results from different tools and databases to the same ontology, namely ARO (Antibiotic Resistance Ontology).

positional arguments:
  {argsoap,abricate,deeparg,resfinder,amrfinderplus}
                        The tool you used to do ARG annotation.

options:
  -h, --help            show this help message and exit
  --db {sarg,ncbi,resfinder,deeparg,megares,argannot}
                        The database you used to do ARG annotation.
  --hamronized          Use this if the input is hamronized (processed using the hAMRonization tool)
  -i INPUT, --input INPUT
                        The annotation result you have
  -o OUTPUT, --output OUTPUT
                        The file to save normalization results
```

Create a folder called `argNorm_tutorial` and store the downloaded data file in it. Navigate into the `argNorm_tutorial` folder.

### Step 1: Download Input Data

argNorm adds columns to the outputs of ARG annotation tools.

For this example, we'll run argNorm on the ARG annotation output from the [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) tool (with ResFinder as the reference database for the tool).

Click [here](https://raw.githubusercontent.com/BigDataBiology/argNorm/main/examples/raw/resfinder.resfinder.orfs.tsv) to download the input data.

If you are on Linux:
```
wget https://raw.githubusercontent.com/BigDataBiology/argNorm/main/examples/raw/resfinder.resfinder.orfs.tsv
```

### Step 2: Run argNorm

Here is a basic outline of most argNorm commands:

```bash
argnorm [tool] -i [original_annotation.tsv] -o [argnorm_result.tsv] [--hamronized]
```

Here, `tool` refers to the ARG annotation tool used (ResFinder in this case). `original_annotation.tsv` is the path to the input data and `argnorm_result.tsv` is the path to output file where the resulting table from argNorm will be stored. `--hamronized` is an option to indicate if the input data is a result of using the [hAMRonization package](https://github.com/pha4ge/hAMRonization). In our example, the input data is not a result of using the hAMRonization package.

To run argNorm on our input data, use this command in your terminal:

```
argnorm resfinder -i ./resfinder.resfinder.orfs.tsv -o ./resfinder.resfinder.orfs.normed.tsv
```

The argNorm result will be stored in the file `resfinder.resfinder.orf.normed.tsv`.

For more examples of using argNorm as a command line tool, click [here](cli_examples.md).

## Quick Example 2: argNorm as a Library

Here is a quick demo of using argNorm as a Python library.

### Code

Save this piece of Python code to a file called `argnorm_tutorial.py`

```
# Import from argNorm
from argnorm.general import map_to_aro
from argnorm.drug_categorization import confers_resistance_to, drugs_to_drug_classes

# Creating a list of input genes to be mapped to the ARO
list_of_input_genes = ['sul1_2_U12338', 'sul1_3_EU855787', 'sul2_1_AF542061']

# The database from which the `list_of_input_genes` was created is the SARG database
database = 'sarg'

# Looping through `list_of_input` genes and mapping each gene to the ARO
# Storing each ARO mapping in the `list_of_aros` list
list_of_aros = []
for gene in list_of_input_genes:
    list_of_aros.append(map_to_aro(gene, 'sarg'))
print(list_of_aros)

# Looping through `list_of_aros` and finding the drugs to which the each ARO confers resistance to
# Storing each drug in the `list_of_drugs` list
list_of_drugs = []
for aro in list_of_aros:
    list_of_drugs.append(confers_resistance_to(aro))
print(list_of_drugs)

# Looping through `list_of_drugs` and finding the superclass/drug class of each drug
# Storing each superclass/drug class in the `list_of_drug_classes` list
list_of_drug_classes = []
for drug in list_of_drugs:
    list_of_drug_classes.append(drugs_to_drug_classes(drug))
print(list_of_drug_classes)
```

### Explanation

To use argNorm as a library, we must first import it in our Python file:

```
from argnorm.general import map_to_aro
from argnorm.drug_categorization import confers_resistance_to, drugs_to_drug_classes
```

The `general` module contains the function `map_to_aro` which will return the ARO number of a particular antibiotic resistance gene. The `drug_categorization` module contains the functions `confers_resistance_to` and `drugs_to_drug_classes`. The `confers_resistance_to` function returns the drugs to which a gene confers resistance to. The `drugs_to_drug_classes` function returns the drug class to which a specific drug belongs.

The `map_to_aro` function takes two arguments: `gene` and `database`. `gene` is the name of an antibiotic resistance gene. `database` is the database from which `gene` is taken from.

In this example, a list of genes from the SARG database is used, and `map_to_aro` maps each gene in the list to an ARO term. These ARO terms are also stored in the `list_of_aros` list:

```
database = 'sarg'

list_of_aros = []
for gene in list_of_input_genes:
    list_of_aros.append(map_to_aro(gene, 'sarg'))
print(list_of_aros)
```

Once a list of AROs is created for each gene, the `confers_resistance_to` function can be used on each ARO to create a list of drugs to which each gene/ARO confers resistance to:

```
list_of_drugs = []
for aro in list_of_aros:
    list_of_drugs.append(confers_resistance_to(aro))
print(list_of_drugs)
```

Now, each drug in the `list_of_drugs` can be categorized into a broader drug category using the `drugs_to_drug_classes` function:

```
list_of_drug_classes = []
for drug in list_of_drugs:
    list_of_drug_classes.append(drugs_to_drug_classes(drug))
print(list_of_drug_classes)
```

For more examples of argNorm as a library, click [here](lib_examples.md).

## API Reference

Click [here](api_reference) to view the API reference
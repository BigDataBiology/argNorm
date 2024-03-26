# argNorm

[![Python package](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml/badge.svg)](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml)
[![Downloads](https://pepy.tech/badge/argNorm)](https://pepy.tech/project/argNorm)
![](https://img.shields.io/badge/status-alpha-red?style=flat)


## What is argNorm?
argNorm is a tool to normalize antibiotic resistance genes (ARGs) by mapping them to the
[antibiotic resistance ontology (ARO)](https://obofoundry.org/ontology/aro.html) by CARD. It also provides drug categorization of drugs that antibiotic resistance genes confer resistance to.

## Why argNorm?

Right now, many tools exist for annotating ARGs in genomes and metagenomes. However, each tool will have distinct output formats.

The [hAMRonization](https://github.com/pha4ge/hAMRonization) package can normalize file formats, but each tool will use different names/identifiers (_e.g._, `TetA` or `TETA` or `tet(A)` or `tet-A` are all different ways to spell the same gene name).

For a small number of isolate genomes, a human user can still quickly evaluate the outputs.
However, for metagenomics, especially for large-scale projects, this becomes infeasible.
Thus, `argNorm` normalizes the _output vocabulary_ of these tools by mapping all tools to the same ontology (ARO).

### Note:
*This is a beta-quality implementation (subject to changes and some bugs may remain), but you're welcomed to try it and provide feedback.*

*We welcome your feedback on the [Issues Page](https://github.com/BigDataBiology/argNorm/issues).*

## Tutorial Video

[![argNorm Tutorial](https://markdown-videos-api.jorgenkh.no/url?url=https%3A%2F%2Fyoutu.be%2Fvx8MCQ7gDLs)](https://youtu.be/vx8MCQ7gDLs)

## Installation
argNorm can be installed using pip as shown.
```bash
pip install argnorm
```

## Supported tools

Note that CARD RGI already uses ARO, thus there is no need to use argNorm.

- [DeepARG](https://bench.cs.vt.edu/deeparg) (v1.0.2)
- [ARGs-OAP](https://galaxyproject.org/use/args-oap/) (v3)
- [ABRicate](https://github.com/tseemann/abricate) (v1.0.1) with NCBI (v3.6), ResFinder (v4.1.11), MEGARes (v2.0), ARG-ANNOT (v5), ResFinderFG (v2)
- [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) (v4.0)
- [AMRFinderPlus](https://github.com/ncbi/amr) (v3.10.30)

## Basic usage

The only positional argument required is `tool` which can be:
- `deeparg`
- `argsoap`
- `abricate`
- `resfinder`
- `amrfinderplus`

The available options are:
- `-h` or `--help`: shows available options and exits.
- `--db`: database used to perform ARG annotation. Supported databases are:
    - SARG (`sarg`)
    - NCBI (`ncbi`)
    - ResFinder (`resfinder`)
    - DeepARG (`deeparg`)
    - MEGARes (`megares`)
    - ARG-ANNOT (`argannot`)
- `--hamronized`: use this if the input is hamronized by [hAMRonization](https://github.com/pha4ge/hAMRonization)
- `-i` or `--input`: path to the annotation result
- `-o` or `--output`: the file to save normalization results

Use `argnorm -h` or `argnorm --help` to see available options.

```bash
>argnorm -h
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

Here is a basic outline of a command using argNorm.

```bash
argnorm [tool] -i [original_annotation.tsv] -o [annotation_result_with_aro.tsv]
```

## Examples

### ARGs-OAP

```bash
argnorm argsoap -i examples/raw/args-oap.sarg.reads.tsv -o outputs/raw/args-oap.sarg.reads.tsv

argnorm argsoap -i examples/hamronized/args-oap.sarg.reads.tsv -o outputs/hamronized/args-oap.sarg.reads.tsv --hamronized
```

### DeepARG

```bash
argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o outputs/raw/deeparg.deeparg.orfs.tsv

argnorm deeparg -i examples/hamronized/deeparg.deeparg.orfs.tsv -o outputs/hamronized/deeparg.deeparg.orfs.tsv --hamronized
```

### ABRicate

#### Hamronized
```bash
argnorm abricate --db ncbi -i examples/hamronized/abricate.ncbi.tsv -o outputs/hamronized/abricate.ncbi.tsv --hamronized
argnorm abricate --db megares -i examples/hamronized/abricate.megares.tsv -o outputs/hamronized/abricate.megares.tsv --hamronized
argnorm abricate --db argannot -i examples/hamronized/abricate.argannot.tsv -o outputs/hamronized/abricate.argannot.tsv --hamronized
argnorm abricate --db resfinder -i examples/hamronized/abricate.resfinder.tsv -o outputs/hamronized/abricate.resfinder.tsv --hamronized
```

#### Raw
```bash
argnorm abricate --db ncbi -i examples/raw/abricate.ncbi.tsv -o outputs/raw/abricate.ncbi.tsv
argnorm abricate --db megares -i examples/raw/abricate.megares.tsv -o outputs/raw/abricate.megarest.tsv
argnorm abricate --db argannot -i examples/raw/abricate.argannot.tsv -o outputs/raw/abricate.argannot.tsv
```

### ResFinder

#### Hamronized
```bash
argnorm resfinder -i examples/hamronized/resfinder.resfinder.orfs.tsv -o outputs/hamronized/resfinder.resfinder.orfs.tsv --hamronized
argnorm resfinder -i examples/hamronized/resfinder.resfinder.reads.tsv -o outputs/hamronized/resfinder.resfinder.reads.tsv --hamronized
```

#### Raw
```bash
argnorm resfinder -i examples/raw/resfinder.resfinder.orfs.tsv -o outputs/raw/resfinder.resfinder.orfs.tsv
argnorm resfinder -i examples/raw/resfinder.resfinder.reads.tsv -o outputs/raw/resfinder.resfinder.reads.tsv
```

### AMRFinderPlus
```bash
argnorm amrfinderplus -i examples/raw/amrfinderplus.ncbi.orfs.tsv -o outputs/raw/amrfinderplus.ncbi.orfs.tsv

argnorm amrfinderplus -i examples/hamronized/amrfinderplus.ncbi.orfs.tsv -o outputs/hamronized/amrfinderplus.ncbi.orfs.tsv
```

## Drug Categorization

Besides performing normalization, argNorm also provides drug categorization of drugs that antibiotic resistance genes confer resistance to.

For example, the `PBP2b` (`ARO:3003042`) gene confers resistance to the drug class `amoxicillin`. `amoxicillin` is then categorized into a broader category of `beta lactam antibiotic`.

argNorm provides support for this, and adds the `confers_resistance_to` and `resistance_to_drug_classes` columns to ARG annotations.

The `confers_resistance_to` column will contain ARO numbers of all the drug classes that a gene provides resistance to (`ARO:0000064` for `amoxicillin` in the previous example).

The `resistance_to_drug_classes` column will contain ARO numbers of the broader categories of the drug classes in the `confers_resistance_to` column (`ARO:3000007` for `beta lactam antibiotic` in the previous example).

## Authors

- Vedanth Ramji [vedanth.ramji@gmail.com](mailto:vedanth.ramji@gmail.com)*
- Hui Chong [huichong.me@gmail.com](mailto:huichong.me@gmail.com)
- Svetlana Ugarcina Perovic [svetlana.ugarcina@gmail.com](mailto:svetlana.ugarcina@gmail.com)
- Finlay Maguire [finlaymaguire@gmail.com](mailto:finlaymaguire@gmail.com)
- [Luis Pedro Coelho](https://luispedro.org) [luis@luispedro.org](mailto:luis@luispedro.org)

*: current maintainer
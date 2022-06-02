# argNorm

[![Python package](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml/badge.svg)](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml)
[![Pylint](https://github.com/BigDataBiology/argNorm/actions/workflows/pylint.yml/badge.svg)](https://github.com/BigDataBiology/argNorm/actions/workflows/pylint.yml)
[![Downloads](https://pepy.tech/badge/argNorm)](https://pepy.tech/project/argNorm)
![](https://img.shields.io/badge/status-alpha-red?style=flat) 
<!-- ![](https://img.shields.io/github/license/BigDataBiology/argNorm?style=flat) -->

Fast antibiotic resistance gene (ARG) normalization by mapping to the [antibiotic resistance ontology (ARO)](https://obofoundry.org/ontology/aro.html) by CARD.

This is a very-first implementation (**not ready for production**), but you're welcomed to try it and provide feedback to make it better. 

We welcome your feedback on the [Issues Page](https://github.com/BigDataBiology/argNorm/issues). 

## Supported databases

- [x] [deeparg](https://bitbucket.org/gusphdproj/deeparg-largerepo/src/master/database/v2/features.fasta)
- [ ] [sarg](https://smile.hku.hk/SARGs/static/images/Ublastx_stageone2.3.tar.gz)
- [x] [ncbi](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt)
- [x] [argannot](https://github.com/tseemann/abricate/tree/master/db/argannot)
- [x] [megares](https://github.com/tseemann/abricate/tree/master/db/megares)
- [x] [resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db)


## Installation

```bash
pip install argnorm
```

## Basic usage

Use `argnorm -h` too see available options.

```bash
argnorm [database] -i [original_annotation.tsv] -o [annotation_result_with_aro.tsv]
```

## Examples

```bash
argnorm deeparg -i examples/deeparg.deeparg.tsv -o tmp
argnorm megares -i examples/abricate.megares.tsv -o tmp
argnorm argannot -i examples/abricate.argannot.tsv -o tmp
argnorm resfinder -i examples/abricate.resfinder.tsv -o tmp
argnorm ncbi -i examples/abricate.ncbi.tsv -o tmp
```

## Maintainer

|   Name    | Email                 | Organization                                                 |
| :-------: | --------------------- | ------------------------------------------------------------ |
| Hui Chong | huichong.me@gmail.com | Research Assistant, Big Data Biology Lab, Fudan University |
| Svetlana Ugarcina Perovic | svetlana.ugarcina@gmail.com | Postdoc Researcher, Big Data Biology Lab, Fudan University |
| Luis Pedro Coelho | luis@luispedro.org | Principle Investigator, Big Data Biology Lab, Fudan University |

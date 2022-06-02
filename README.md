# argNorm

[![Python package](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml/badge.svg)](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml)
[![Pylint](https://github.com/BigDataBiology/argNorm/actions/workflows/pylint.yml/badge.svg)](https://github.com/BigDataBiology/argNorm/actions/workflows/pylint.yml)
[![Downloads](https://pepy.tech/badge/argNorm)](https://pepy.tech/project/argNorm)
![](https://img.shields.io/badge/status-alpha-red?style=flat) 
<!-- ![](https://img.shields.io/github/license/BigDataBiology/argNorm?style=flat) -->

Fast ARG normalization by mapping to the ARO ontology.

This is a very-first implementation (**not ready for production**), but you're welcomed to try it and provide feedback to make it better. 

We recieve feedback on [GitHub Issue](https://github.com/BigDataBiology/argNorm/issues). 

## Supported databases

- [x] deeparg
- [ ] sarg
- [x] ncbi
- [x] argannot
- [x] megares
- [x] resfinder


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
| Svetlana Ugarcina | svetlana.ugarcina@gmail.com | Postdoc Researcher, Big Data Biology Lab, Fudan University |
| Luis Pedro Coelho | luis@luispedro.org | Principle Investigator, Big Data Biology Lab, Fudan University |
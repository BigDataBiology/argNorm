# argNorm

[![Python package](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml/badge.svg)](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml)
[![Downloads](https://pepy.tech/badge/argNorm)](https://pepy.tech/project/argNorm)
![](https://img.shields.io/badge/status-alpha-red?style=flat) 
<!-- ![](https://img.shields.io/github/license/BigDataBiology/argNorm?style=flat) -->

Fast antibiotic resistance gene (ARG) normalization by mapping to the [antibiotic resistance ontology (ARO)](https://obofoundry.org/ontology/aro.html) by CARD.

This is a very-first implementation (**not ready for production**), but you're welcomed to try it and provide feedback to make it better. 

We welcome your feedback on the [Issues Page](https://github.com/BigDataBiology/argNorm/issues). 

## Supported tools

- [x] [DeepARG](https://bench.cs.vt.edu/deeparg)
- [x] [ARGs-OAP](https://galaxyproject.org/use/args-oap/)
- [ ] [Abricate](https://github.com/tseemann/abricate) (work in progress)

<!-- Hamronized output

- [x] [deeparg](https://bitbucket.org/gusphdproj/deeparg-largerepo/src/master/database/v2/features.fasta)
- [x] [sarg](https://smile.hku.hk/SARGs/static/images/Ublastx_stageone2.3.tar.gz)
- [x] [ncbi](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt)
- [x] [argannot](https://github.com/tseemann/abricate/tree/master/db/argannot)
- [x] [megares](https://github.com/tseemann/abricate/tree/master/db/megares)
- [x] [resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db) -->

<!-- Raw output

- [x] [deeparg](https://bitbucket.org/gusphdproj/deeparg-largerepo/src/master/database/v2/features.fasta)
- [ ] [sarg](https://smile.hku.hk/SARGs/static/images/Ublastx_stageone2.3.tar.gz)
- [x] [ncbi](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt)
- [x] [argannot](https://github.com/tseemann/abricate/tree/master/db/argannot)
- [x] [megares](https://github.com/tseemann/abricate/tree/master/db/megares)
- [ ] [resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db) -->

## Installation

```bash
pip install argnorm
```

## Basic usage

Use `argnorm -h` to see available options.

```bash
argnorm [database] -i [original_annotation.tsv] -o [annotation_result_with_aro.tsv]
```

## Examples

```bash
argnorm argsoap --mode reads -i examples/raw/args-oap.sarg.reads.tsv -o tmp
argnorm argsoap --mode reads -i examples/hamronized/args-oap.sarg.reads.tsv -o tmp --hamronized
argnorm argsoap --mode orfs -i examples/raw/args-oap.sarg.orfs.tsv -o tmp
argnorm argsoap --mode orfs -i examples/hamronized/args-oap.sarg.orfs.tsv -o tmp --hamronized

argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o tmp
argnorm deeparg -i examples/hamronized/deeparg.deeparg.orfs.tsv -o tmp --hamronized
```

<!-- Handling hamronized inputs.

```bash
argnorm deeparg -i examples/hamronized/deeparg.deeparg.orfs.tsv -o tmp
argnorm megares -i examples/hamronized/abricate.megares.tsv -o tmp
argnorm argannot -i examples/hamronized/abricate.argannot.tsv -o tmp
argnorm resfinder -i examples/hamronized/abricate.resfinder.tsv -o tmp
argnorm ncbi -i examples/hamronized/abricate.ncbi.tsv -o tmp
```

Handling raw inputs (not hamronized by hAMRonization)

```bash
argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o tmp --raw
argnorm megares -i examples/raw/abricate.megares.tsv -o tmp --raw
argnorm argannot -i examples/raw/abricate.argannot.tsv -o tmp --raw
# argnorm resfinder -i examples/raw/abricate.resfinder.tsv -o tmp --raw
argnorm ncbi -i examples/raw/abricate.ncbi.tsv -o tmp --raw
``` -->

<!-- ## TODO

- [ ] Add aro categories in outputs, see reference file at: ../quick_amr_db_harmonisation/dbs/aro_categories.tsv.
- [ ] Add sarg support by including the reference file (../quick_amr_db_harmonisation/Ublastx_stageone2.3/DB/structure_20181107.LIST), but how to deal with multiple ARO matches of each query result.

Note: Sarg reads mode and orfs mode output formats are different.  -->

## Maintainers

|   Name    | Email                 | Organization                                                 |
| :-------: | --------------------- | ------------------------------------------------------------ |
| Hui Chong | huichong.me@gmail.com | Research Assistant, Big Data Biology Lab, Fudan University |
| Svetlana Ugarcina Perovic | svetlana.ugarcina@gmail.com | Postdoc Researcher, Big Data Biology Lab, Fudan University |
| Luis Pedro Coelho | luis@luispedro.org | Principle Investigator, Big Data Biology Lab, Fudan University |
| Finlay Maguire | finlaymaguire@gmail.com | Assistant Professor, Department of Community Health and Epidemiology, Dalhousie University |
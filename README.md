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
- [x] [ABRicate](https://github.com/tseemann/abricate) with NCBI, ResFinder, MEGARes, ARG-ANNOT
- [x] [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/)
- [x] [AMRFinderPlus](https://github.com/ncbi/amr)
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
argnorm [tool] -i [original_annotation.tsv] -o [annotation_result_with_aro.tsv]
```

## Examples

```bash
argnorm argsoap --mode reads -i examples/raw/args-oap.sarg.reads.tsv -o tmp
argnorm argsoap --mode reads -i examples/hamronized/args-oap.sarg.reads.tsv -o tmp --hamronized
argnorm argsoap --mode orfs -i examples/raw/args-oap.sarg.orfs.tsv -o tmp
argnorm argsoap --mode orfs -i examples/hamronized/args-oap.sarg.orfs.tsv -o tmp --hamronized

argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o tmp
argnorm deeparg -i examples/hamronized/deeparg.deeparg.orfs.tsv -o tmp --hamronized

argnorm abricate --db ncbi -i examples/raw/abricate.ncbi.tsv -o tmp
argnorm abricate --db megares -i examples/raw/abricate.megares.tsv -o tmp
argnorm abricate --db argannot -i examples/raw/abricate.argannot.tsv -o tmp

argnorm abricate --db ncbi -i examples/hamronized/abricate.ncbi.tsv -o tmp --hamronized
argnorm abricate --db megares -i examples/hamronized/abricate.megares.tsv -o tmp --hamronized
argnorm abricate --db argannot -i examples/hamronized/abricate.argannot.tsv -o tmp --hamronized
argnorm abricate --db resfinder -i examples/hamronized/abricate.resfinder.tsv -o tmp --hamronized

argnorm resfinder -i examples/raw/resfinder.resfinder.orfs.tsv -o outputs/raw/resfinder.resfinder.orfs.tsv
argnorm resfinder -i examples/raw/resfinder.resfinder.reads.tsv -o outputs/raw/resfinder.resfinder.reads.tsv

argnorm amrfinderplus -i examples/raw/amrfinderplus.ncbi.orfs.tsv -o outputs/raw/amrfinderplus.ncbi.orfs.tsv

```

## Maintainers

|   Name    | Email                 | Organization                                                 |
| :-------: | --------------------- | ------------------------------------------------------------ |
| Hui Chong | huichong.me@gmail.com | Research Assistant, Big Data Biology Lab, Fudan University |
| Svetlana Ugarcina Perovic | svetlana.ugarcina@gmail.com | Postdoc Researcher, Big Data Biology Lab, Fudan University |
| Finlay Maguire | finlaymaguire@gmail.com | Assistant Professor, Department of Community Health and Epidemiology, Dalhousie University |
| [Luis Pedro Coelho](https://luispedro.org) | luis@luispedro.org | Principal Investigator, Big Data Biology Lab, Fudan University |

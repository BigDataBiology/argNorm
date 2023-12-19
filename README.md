# argNorm

[![Python package](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml/badge.svg)](https://github.com/BigDataBiology/argNorm/actions/workflows/python-package.yml)
[![Downloads](https://pepy.tech/badge/argNorm)](https://pepy.tech/project/argNorm)
![](https://img.shields.io/badge/status-alpha-red?style=flat) 
<!-- ![](https://img.shields.io/github/license/BigDataBiology/argNorm?style=flat) -->

Antibiotic resistance gene (ARG) normalization by mapping to the [antibiotic resistance ontology (ARO)](https://obofoundry.org/ontology/aro.html) by CARD.

Right now, many tools exist for annotating ARGs in genomes and metagenomes.
However, each tool will output in its own format.
The [hAMRonization](https://github.com/pha4ge/hAMRonization) package can normalize file formats, but each tool will use different names/identifiers (_e.g._, `TetA` or `TETA` or `tet(A)` or `tet-A` are all different ways to spell the same gene name).
For a small number of isolate genomes, a human user can still quickly evaluate the outputs.
However, for metagenomics, especially for large-scale projects, this becomes infeasible.
Thus, `argNorm` normalizes the _output vocabulary_ of these tools by mapping all tools to the same ontology (ARO).

This is a beta-quality implementation (subject to changes and some bugs may remain), but you're welcomed to try it and provide feedback.

We welcome your feedback on the [Issues Page](https://github.com/BigDataBiology/argNorm/issues). 

## Supported tools

Note that CARD RGI already uses ARO, thus there is no need to use argNorm.

- [x] [DeepARG](https://bench.cs.vt.edu/deeparg) (v1.0.2)
- [x] [ARGs-OAP](https://galaxyproject.org/use/args-oap/) (v2.3)
- [x] [ABRicate](https://github.com/tseemann/abricate) (v1.0.1) with NCBI (v3.6), ResFinder (v4.1.11), MEGARes (v2.0), ARG-ANNOT (v5)
- [x] [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) (v4.0)
- [x] [AMRFinderPlus](https://github.com/ncbi/amr) (v3.10.30)
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

## Authors

- Vedanth Ramji [vedanth.ramji@gmail.com](mailto:vedanth.ramji@gmail.com)*
- Hui Chong [huichong.me@gmail.com](mailto:huichong.me@gmail.com)
- Svetlana Ugarcina Perovic [svetlana.ugarcina@gmail.com](mailto:svetlana.ugarcina@gmail.com)
- Finlay Maguire [finlaymaguire@gmail.com](mailto:finlaymaguire@gmail.com)
- [Luis Pedro Coelho](https://luispedro.org) [luis@luispedro.org](mailto:luis@luispedro.org)

*: current maintainer

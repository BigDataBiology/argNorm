#!/bin/bash

echo "Testing ARGs-OAP normalization..."
argnorm argsoap --mode reads -i examples/raw/args-oap.sarg.reads.tsv -o outputs/raw/args-oap.sarg.reads.tsv
argnorm argsoap --mode reads -i examples/hamronized/args-oap.sarg.reads.tsv -o outputs/hamronized/args-oap.sarg.reads.tsv --hamronized
argnorm argsoap --mode orfs -i examples/raw/args-oap.sarg.orfs.tsv -o outputs/raw/args-oap.sarg.orfs.tsv
argnorm argsoap --mode orfs -i examples/hamronized/args-oap.sarg.orfs.tsv -o outputs/hamronized/args-oap.sarg.orfs.tsv --hamronized

echo "Testing DeepARG normalization..."
argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o outputs/raw/deeparg.deeparg.orfs.tsv
argnorm deeparg -i examples/hamronized/deeparg.deeparg.orfs.tsv -o outputs/hamronized/deeparg.deeparg.orfs.tsv --hamronized

echo "Testing Abricate normalization..."
argnorm abricate --db ncbi -i examples/raw/abricate.ncbi.tsv -o outputs/raw/abricate.ncbi.tsv
argnorm abricate --db megares -i examples/raw/abricate.megares.tsv -o outputs/raw/abricate.megares.tsv
argnorm abricate --db argannot -i examples/raw/abricate.argannot.tsv -o outputs/raw/abricate.argannot.tsv
# argnorm abricate --db resfinder -i examples/raw/abricate.resfinder.tsv -o outputs/raw/abricate.resfinder.tsv

argnorm abricate --db ncbi -i examples/hamronized/abricate.ncbi.tsv -o outputs/hamronized/abricate.ncbi.tsv --hamronized
argnorm abricate --db megares -i examples/hamronized/abricate.megares.tsv -o outputs/hamronized/abricate.megares.tsv --hamronized
argnorm abricate --db argannot -i examples/hamronized/abricate.argannot.tsv -o outputs/hamronized/abricate.argannot.tsv --hamronized
argnorm abricate --db resfinder -i examples/hamronized/abricate.resfinder.tsv -o outputs/hamronized/abricate.resfinder.tsv --hamronized

argnorm resfinder -i examples/raw/resfinder.resfinder.orfs.tsv -o outputs/raw/resfinder.resfinder.orfs.tsv
argnorm resfinder -i examples/raw/resfinder.resfinder.reads.tsv -o outputs/raw/resfinder.resfinder.reads.tsv

argnorm amrfinderplus -i examples/raw/amrfinderplus.ncbi.orfs.tsv -o outputs/raw/amrfinderplus.ncbi.orfs.tsv

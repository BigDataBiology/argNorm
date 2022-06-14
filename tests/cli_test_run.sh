#!/bin/bash

echo "Testing ARGs-OAP normalization..."
argnorm argsoap --mode reads -i examples/raw/args-oap.sarg.reads.tsv -o tmp
argnorm argsoap --mode reads -i examples/hamronized/args-oap.sarg.reads.tsv -o tmp --hamronized
argnorm argsoap --mode orfs -i examples/raw/args-oap.sarg.orfs.tsv -o tmp
argnorm argsoap --mode orfs -i examples/hamronized/args-oap.sarg.orfs.tsv -o tmp --hamronized

echo "Testing DeepARG normalization..."
argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o tmp
argnorm deeparg -i examples/hamronized/deeparg.deeparg.orfs.tsv -o tmp --hamronized

echo "Testing Abricate normalization..."
argnorm abricate --db ncbi -i examples/raw/abricate.ncbi.tsv -o tmp
argnorm abricate --db megares -i examples/raw/abricate.megares.tsv -o tmp
argnorm abricate --db argannot -i examples/raw/abricate.argannot.tsv -o tmp
# argnorm abricate --db resfinder -i examples/raw/abricate.resfinder.tsv -o tmp

argnorm abricate --db ncbi -i examples/hamronized/abricate.ncbi.tsv -o tmp --hamronized
argnorm abricate --db megares -i examples/hamronized/abricate.megares.tsv -o tmp --hamronized
argnorm abricate --db argannot -i examples/hamronized/abricate.argannot.tsv -o tmp --hamronized
argnorm abricate --db resfinder -i examples/hamronized/abricate.resfinder.tsv -o tmp --hamronized

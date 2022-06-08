#!/bin/bash


echo "Testing normalization of hamronized inputs."
argnorm deeparg -i examples/hamronized/deeparg.deeparg.orfs.tsv -o tmp
argnorm megares -i examples/hamronized/abricate.megares.tsv -o tmp
argnorm argannot -i examples/hamronized/abricate.argannot.tsv -o tmp
argnorm resfinder -i examples/hamronized/abricate.resfinder.tsv -o tmp
argnorm ncbi -i examples/hamronized/abricate.ncbi.tsv -o tmp

echo "Testing normalization of raw inputs."
argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o tmp --raw
argnorm megares -i examples/raw/abricate.megares.tsv -o tmp --raw
argnorm argannot -i examples/raw/abricate.argannot.tsv -o tmp --raw
argnorm ncbi -i examples/raw/abricate.ncbi.tsv -o tmp --raw

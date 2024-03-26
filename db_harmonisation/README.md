# Mapping databases to RGI

This is useful to _build_ the argNorm mapping databases, it is **not** intended for regular users of argNorm.

## Setting up conda environment

```bash
conda create --name rgi --channel conda-forge --channel bioconda --channel defaults rgi
conda activate rgi
conda install jug
```

## Running the code

This is a [Jug](https://jug.rtfd.io/) project, so you run `jug execute`

```bash
jug execute crude_db_harmonisation.py
```

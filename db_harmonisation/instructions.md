# Setting up conda environment

```
conda create --name rgi --channel conda-forge --channel bioconda --channel defaults rgi
```

```
conda activate rgi
```

```
pip install jug
```

# Running the code

```
jug execute crude_db_harmonisation.py
```
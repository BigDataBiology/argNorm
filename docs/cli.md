# Using argNorm as a command line tool

## Quick example

Here is a quick demo of running argNorm on the command line.

### Step 1: Install argNorm

Install argNorm and check installation
```
pip install argnorm
argnorm -h
```

`argnorm -h` or `argnorm --help` will display all the available options to run argNorm with.

```
> argnorm -h
usage: argnorm [-h] [--db {argannot,deeparg,megares,ncbi,resfinder,resfinderfg,sarg,groot-db,groot-core-db,groot-argannot,groot-resfinder,groot-card}]
               [-i INPUT] [--hamronized] [--hamronization_skip_unsupported_tool] [-o OUTPUT]
               {argsoap,abricate,deeparg,resfinder,amrfinderplus,groot,hamronization}

argNorm normalizes ARG annotation results from different tools and databases to the same ontology, namely ARO (Antibiotic Resistance Ontology).

positional arguments:
  {argsoap,abricate,deeparg,resfinder,amrfinderplus,groot,hamronization}
                        The tool you used to do ARG annotation.

options:
  -h, --help            show this help message and exit
  --db {argannot,deeparg,megares,ncbi,resfinder,resfinderfg,sarg,groot-db,groot-core-db,groot-argannot,groot-resfinder,groot-card}
                        The database you used to do ARG annotation.
  -i INPUT, --input INPUT
                        The annotation result you have
  --hamronized          Use hamronization as a tool instead
  --hamronization_skip_unsupported_tool
                        Skip rows with unsupported tools for hamronization outputs
  -o OUTPUT, --output OUTPUT
                        The file to save normalization results
```

### Step 2: Create working directory & download sample data

argNorm adds ARO mappings and drug categories as additional columns to the outputs of ARG annotation tools.

For this example, we will run argNorm on the ARG annotation output from the [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) tool (with the [ResFinder database](https://bitbucket.org/genomicepidemiology/resfinder_db/raw/8aad1d20603fbec937cdae55024568de6dbd609f/all.fsa)). 

Create a folder called `argNorm_tutorial` and store the downloaded data file in it. Navigate into the `argNorm_tutorial` folder.

```
mkdir argNorm_tutorial
cd argNorm_tutorial
```

Click [here](https://raw.githubusercontent.com/BigDataBiology/argNorm/main/examples/raw/resfinder.resfinder.orfs.tsv) to download the input data.

If you are on Linux:
```
wget https://raw.githubusercontent.com/BigDataBiology/argNorm/main/examples/raw/resfinder.resfinder.orfs.tsv
```

### Step 3: Running argNorm

Here is a basic outline of most argNorm commands:

```bash
argnorm [tool] -i [original_annotation.tsv] -o [argnorm_result.tsv] [--db]
```

Here, `tool` refers to the ARG annotation tool used (ResFinder in this case). `original_annotation.tsv` is the path to the input data and `argnorm_result.tsv` is the path to output file where the resulting table from argNorm will be stored. `--db` is the ARG databases used along with `tool` to perform annotation. ResFinder does not require a `--db` (argNorm will automatically load up the ResFinder database), however, `--db` is required for the ARG annotation tools `groot` and `abricate`.


To run argNorm on the input data, use this command in your terminal:

```
argnorm resfinder -i ./resfinder.resfinder.orfs.tsv -o ./resfinder.resfinder.orfs.normed.tsv
```

The argNorm result will be stored in the file `resfinder.resfinder.orf.normed.tsv`.

## Examples for each tool

Create a folder called `argnorm_cli_tutorial` and download the [sample data](https://github.com/BigDataBiology/argNorm/tree/main/examples) inside it.

```
mkdir argnorm_cli_tutorial
cd argnorm_cli_tutorial
```

Here is a basic outline of calling argNorm.

```bash
argnorm [tool] -i [original_annotation.tsv] -o [annotation_result_with_aro.tsv]
```

### ARGs-OAP

```bash
argnorm argsoap -i examples/raw/args-oap.sarg.reads.tsv -o outputs/raw/args-oap.sarg.reads.tsv
```

### DeepARG

```bash
argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o outputs/raw/deeparg.deeparg.orfs.tsv
```

### ABRicate

When using abricate, it is necessary to specify the database used:

```bash
argnorm abricate --db ncbi -i examples/raw/abricate.ncbi.tsv -o outputs/raw/abricate.ncbi.tsv
argnorm abricate --db megares -i examples/raw/abricate.megares.tsv -o outputs/raw/abricate.megarest.tsv
argnorm abricate --db argannot -i examples/raw/abricate.argannot.tsv -o outputs/raw/abricate.argannot.tsv
```

### ResFinder

```bash
argnorm resfinder -i examples/raw/resfinder.resfinder.orfs.tsv -o outputs/raw/resfinder.resfinder.orfs.tsv
argnorm resfinder -i examples/raw/resfinder.resfinder.reads.tsv -o outputs/raw/resfinder.resfinder.reads.tsv
```

### AMRFinderPlus
```bash
argnorm amrfinderplus -i examples/raw/amrfinderplus.ncbi.orfs.v3.10.30.tsv -o outputs/raw/amrfinderplus.ncbi.orfs.v3.10.30.tsv
argnorm amrfinderplus -i examples/raw/amrfinderplus.ncbi.orfs.v4.tsv -o outputs/raw/amrfinderplus.ncbi.orfs.v4.tsv
```

### GROOT
```bash
argnorm groot -i examples/raw/groot.argannot.tsv -o outputs/raw/groot.argannot.tsv --db groot-argannot
argnorm groot -i examples/raw/groot.resfinder.tsv -o outputs/raw/groot.resfinder.tsv --db groot-resfinder
argnorm groot -i examples/raw/groot.card.tsv -o outputs/raw/groot.card.tsv --db groot-card
argnorm groot -i examples/raw/groot.groot-db.tsv -o outputs/raw/groot.groot-db.tsv --db groot-db
argnorm groot -i examples/raw/groot.groot-core-db.tsv -o ouptuts/raw/groot.groot-core-db.tsv --db groot-core-db
```

### Hamronization

```bash
argnorm hamronization -i examples/hamronized/args-oap.sarg.reads.tsv -o outputs/hamronized/args-oap.sarg.reads.tsv
argnorm hamronization -i examples/hamronized/deeparg.deeparg.orfs.tsv -o outputs/hamronized/deeparg.deeparg.orfs.tsv
argnorm hamronization -i examples/hamronized/abricate.ncbi.tsv -o outputs/hamronized/abricate.ncbi.tsv 
argnorm hamronization -i examples/hamronized/abricate.megares.tsv -o outputs/hamronized/abricate.megares.tsv 
argnorm hamronization -i examples/hamronized/abricate.argannot.tsv -o outputs/hamronized/abricate.argannot.tsv 
argnorm hamronization -i examples/hamronized/abricate.resfinder.tsv -o outputs/hamronized/abricate.resfinder.tsv 
argnorm hamronization -i examples/hamronized/resfinder.resfinder.orfs.tsv -o outputs/hamronized/resfinder.resfinder.orfs.tsv 
argnorm hamronization -i examples/hamronized/resfinder.resfinder.reads.tsv -o outputs/hamronized/resfinder.resfinder.reads.tsv 
argnorm hamronization -i examples/hamronized/amrfinderplus.ncbi.orfs.tsv -o outputs/hamronized/amrfinderplus.ncbi.orfs.tsv
argnorm hamronization -i examples/hamronized/groot.argannot.tsv -o outputs/hamronized/groot.argannot.tsv 
argnorm hamronization -i examples/hamronized/groot.resfinder.tsv -o outputs/hamronized/groot.resfinder.tsv
argnorm hamronization -i examples/hamronized/groot.card.tsv -o outputs/hamronized/groot.card.tsv
argnorm hamronization -i examples/hamronized/groot.groot-db.tsv -o outputs/hamronized/groot.groot-db.tsv
argnorm hamronization -i examples/hamronized/groot.groot-core-db.tsv -o outputs/hamronized/groot.groot-core-db.tsv
argnorm hamronization -i examples/hamronized/combined_hamronization.tsv -o outputs/hamronized/combined_hamronization.tsv
argnorm hamronization -i examples/hamronized/combined_hamronization_full.tsv -o outputs/hamronized/combined_hamronization_full.tsv ----hamronization_skip_unsupported_tool
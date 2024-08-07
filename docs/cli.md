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
usage: argnorm [-h]
               [--db {sarg,ncbi,resfinder,deeparg,megares,argannot,resfinderfg}]
               [--hamronized] [-i INPUT] [-o OUTPUT]
               {argsoap,abricate,deeparg,resfinder,amrfinderplus}

argNorm normalizes ARG annotation results from different tools and databases to the same ontology, namely ARO (Antibiotic Resistance Ontology).

positional arguments:
  {argsoap,abricate,deeparg,resfinder,amrfinderplus}
                        The tool you used to do ARG annotation.

optional arguments:
  -h, --help            show this help message and exit
  --db {sarg,ncbi,resfinder,deeparg,megares,argannot,resfinderfg}
                        The database you used to do ARG annotation.
  --hamronized          Use this if the input is hamronized (processed using
                        the hAMRonization tool)
  -i INPUT, --input INPUT
                        The annotation result you have
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
argnorm [tool] -i [original_annotation.tsv] -o [argnorm_result.tsv] [--hamronized]
```

Here, `tool` refers to the ARG annotation tool used (ResFinder in this case). `original_annotation.tsv` is the path to the input data and `argnorm_result.tsv` is the path to output file where the resulting table from argNorm will be stored. `--hamronized` is an option to indicate if the input data is a result of using the [hAMRonization package](https://github.com/pha4ge/hAMRonization). In our example, the input data is not a result of using the hAMRonization package, and so the `--hamronized` option can be omitted.


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

argnorm argsoap -i examples/hamronized/args-oap.sarg.reads.tsv -o outputs/hamronized/args-oap.sarg.reads.tsv --hamronized
```

### DeepARG

```bash
argnorm deeparg -i examples/raw/deeparg.deeparg.orfs.tsv -o outputs/raw/deeparg.deeparg.orfs.tsv

argnorm deeparg -i examples/hamronized/deeparg.deeparg.orfs.tsv -o outputs/hamronized/deeparg.deeparg.orfs.tsv --hamronized
```

### ABRicate

When using abricate, it is necessary to specify the database used:

#### Hamronized
```bash
argnorm abricate --db ncbi -i examples/hamronized/abricate.ncbi.tsv -o outputs/hamronized/abricate.ncbi.tsv --hamronized
argnorm abricate --db megares -i examples/hamronized/abricate.megares.tsv -o outputs/hamronized/abricate.megares.tsv --hamronized
argnorm abricate --db argannot -i examples/hamronized/abricate.argannot.tsv -o outputs/hamronized/abricate.argannot.tsv --hamronized
argnorm abricate --db resfinder -i examples/hamronized/abricate.resfinder.tsv -o outputs/hamronized/abricate.resfinder.tsv --hamronized
```

#### Raw
```bash
argnorm abricate --db ncbi -i examples/raw/abricate.ncbi.tsv -o outputs/raw/abricate.ncbi.tsv
argnorm abricate --db megares -i examples/raw/abricate.megares.tsv -o outputs/raw/abricate.megarest.tsv
argnorm abricate --db argannot -i examples/raw/abricate.argannot.tsv -o outputs/raw/abricate.argannot.tsv
```

### ResFinder

#### Hamronized
```bash
argnorm resfinder -i examples/hamronized/resfinder.resfinder.orfs.tsv -o outputs/hamronized/resfinder.resfinder.orfs.tsv --hamronized
argnorm resfinder -i examples/hamronized/resfinder.resfinder.reads.tsv -o outputs/hamronized/resfinder.resfinder.reads.tsv --hamronized
```

#### Raw
```bash
argnorm resfinder -i examples/raw/resfinder.resfinder.orfs.tsv -o outputs/raw/resfinder.resfinder.orfs.tsv
argnorm resfinder -i examples/raw/resfinder.resfinder.reads.tsv -o outputs/raw/resfinder.resfinder.reads.tsv
```

### AMRFinderPlus
```bash
argnorm amrfinderplus -i examples/raw/amrfinderplus.ncbi.orfs.tsv -o outputs/raw/amrfinderplus.ncbi.orfs.tsv

argnorm amrfinderplus -i examples/hamronized/amrfinderplus.ncbi.orfs.tsv -o outputs/hamronized/amrfinderplus.ncbi.orfs.tsv
```

### GROOT
```bash
argnorm groot -i examples/raw/groot.argannot.tsv -o outputs/raw/groot.argannot.tsv --db groot-argannot
argnorm groot -i examples/raw/groot.resfinder.tsv -o outputs/raw/groot.resfinder.tsv --db groot-resfinder
argnorm groot -i examples/raw/groot.card.tsv -o outputs/raw/groot.card.tsv --db groot-card
argnorm groot -i examples/raw/groot.groot-db.tsv -o outputs/raw/groot.groot-db.tsv --db groot-db
argnorm groot -i examples/raw/groot.groot-core-db.tsv -o ouptuts/raw/groot.groot-core-db.tsv --db groot-core-db

argnorm groot -i examples/hamronized/groot.argannot.tsv -o outputs/hamronized/groot.argannot.tsv --db groot-argannot --hamronized
argnorm groot -i examples/hamronized/groot.resfinder.tsv -o outputs/hamronized/groot.resfinder.tsv --db groot-resfinder --hamronized
argnorm groot -i examples/hamronized/groot.card.tsv -o outputs/hamronized/groot.card.tsv --db groot-card --hamronized
argnorm groot -i examples/hamronized/groot.groot-db.tsv -o outputs/hamronized/groot.groot-db.tsv --db groot-db --hamronized
argnorm groot -i examples/hamronized/groot.groot-core-db.tsv -o outputs/hamronized/groot.groot-core-db.tsv --db groot-core-db --hamronized
```
# What's New

argNorm is under active development and outputs may change between versions. Particular versions of argNorm should, however, be reproducible as the databases are bundled with the tool.

## Version 1.1.0 - released on 24 June 2025

The way that NCBI mappings are handled has been changed, which fixes several issues, particularly in handling hAMROnized data. The [CHANGELOG](https://github.com/BigDataBiology/argNorm/blob/main/CHANGELOG.md) has internal details, but users should see a few changes in the output (for the better). In some cases, previous versions of argNorm returned too optimistic mappings.

## Version 1.0.0 - released on 19 April 2025

Adds citation information to the README and the CLI.

Otherwise, no major changes.

## Version 0.8.0 - released on 12 April 2025

*This is a test release for version 1.0, which will be released as soon as the manuscript is published.*

### User-facing changes

- Added columns to output with English names for AROs, drugs, and drug classes (`ARO_name`, `confers_resistance_to_names`, `resistance_to_drug_classes_names`, corresponding to `ARO`, `confers_resistance_to`, and `resistance_to_drug_classes` columns).
- Added `--version` argument

### Bug fixes
- Fixed reference gene preprocessing for abricate's resfinder option—previous version missed to include a part of the reference genes for abricate's resfinder option.

## Version 0.7.0 - 5 March 2025

### User-facing changes

#### argNorm version in output

The version of argNorm used for normalization is now added to the top of argNorm TSV outputs. For example:

```
# argNorm version: 0.7.0
input_sequence_id	input_file_name	gene_symbol	gene_name	reference_database_id ...
. REST OF ARG ANNOTATION OUTPUT TSV TABLE
.
.
```

**Note** ⚠️: this will break any previous scripts analyzing argnorm cli outputs before this update! Use the `skiprows=1` argument (or equivalent) when loading argnorm output dataframes to ignore the comment with the argnorm version as shown below for Pandas:

```python
import pandas as pd
df = pd.read_csv(<PATH TO ARGNORM OUTPUT>, sep='\t', skiprows=1)
```

#### argNorm mappings updated for CARD and ARO v4.0.0

CARD and ARO were updated to v4.0.0 (released 18 Dec 2024). This brought significant changes to both CARD and the ARO, with 1200+ new beta-lactamase genes being added.

> **NOTE 💡 : ARO mappings for many ARGs have changed. An ARO mapping from a previous version of argNorm might not be present in the latest ARO or can even point to a completely different gene.** 

#### Improved manual curation

We automatically found all cases where the mappings were returning a gene in a different gene class than the original one and manually curated them.

#### Added `Cut_Off` column to outputs
- A new column called `Cut_Off` is now present in all argNorm outputs.
- These cut off scores are procured from RGI outputs after running input databases through RGI, and these cut-off scores correspond to RGI's Discovery paradigm (learn more here: https://github.com/arpcard/rgi/blob/master/docs/rgi_main.rst).
- RGI's cut-off scores are of three types: `Loose`, `Strict`, and `Perfect`.
- argNorm also adds another score, `Manual`, to signify genes that are mapped to ARO terms manually without using RGI.

From RGI documentation (https://github.com/arpcard/rgi/blob/master/docs/rgi_main.rst):
> The RGI analyzes genome or proteome sequences under a Perfect, Strict, and Loose (a.k.a. Discovery) paradigm. The Perfect algorithm is most often applied to clinical surveillance as it detects perfect matches to the curated reference sequences in CARD. In contrast, the Strict algorithm detects previously unknown variants of known AMR genes, including secondary screen for key mutations, using detection models with CARD's curated similarity cut-offs to ensure the detected variant is likely a functional AMR gene. The Loose algorithm works outside of the detection model cut-offs to provide detection of new, emergent threats and more distant homologs of AMR genes, but will also catalog homologous sequences and spurious partial matches that may not have a role in AMR. Combined with phenotypic screening, the Loose algorithm allows researchers to hone in on new AMR genes.
>
> Within the Perfect, Strict, and Loose paradigm, RGI currently supports CARD's protein homolog models, protein variant models, protein over-expression models, and rRNA mutation models:
> 
> Protein Homolog Models (PHM) detect protein sequences based on their similarity to a curated reference sequence, using curated BLASTP bitscore cut-offs, for example NDM-1. Protein Homolog Models apply to all genes that confer resistance through their presence in an organism, such as the presence of a beta-lactamase gene on a plasmid. PHMs include a reference sequence and a bitscore cut-off for detection using BLASTP. A Perfect RGI match is 100% identical to the reference protein sequence along its entire length, a Strict RGI match is not identical but the bit-score of the matched sequence is greater than the curated BLASTP bit-score cutoff, Loose RGI matches have a bit-score less than the curated BLASTP bit-score cut-off.

#### Support amrfinderplus v4 alongisde amrfinderplus v3
- Both the latest version of amrfinderplus (v4.0.19) and the previous v3.10.30 are now supported.

#### Added `hamronization` as a subcommand to the CLI

A better UX for hamronization has been implemented. The `hamronization` subcommand now allows users to hamronize results from all supported tools and databases. The `hamronization` subcommand can be used as follows:

```bash
argnorm hamronization -i PATH_TO_INPUT -o PATH_TO_OUTPUT [--hamronization_skip_unsupported_tool]
```

Note that the `--hamronization_skip_unsupported_tool` flag can be used to skip rows with unsupported tools. A warning will be raised rather than an exception.

On the internal library, this corresponds to the use of the `HamronizationNormalizer` class.

#### Better error messages if output file is missing

Now, an explicit error message is raised if the output argument is not passed on the command line.

### Internal improvements

#### Update `confers_resistance_to()` to use `regulates`, `part_of`, and `participates_in` ARO relationships

Previously, argNorm used the `is_a` ARO relationship along with `confers_resistance_to_drug_class` and `confers_resistance_to_antibiotic` to map ARGs to the drugs they confer resistance to. While this worked well for most genes, some ARGs such as those coding for efflux pumps/proteins (e.g. `ARO:3003548`, `ARO:3000826`, `ARO:3003066`) were previously not mapped to any drugs. This is because none of their superclasses mapped to drugs/antibiotics via `confers_resistance_to_antibiotic` or `confers_resistance_to_drug_class`. However, these genes were related to other ARGs that did map to drugs via the `regulates`, `part_of`, or `participates_in` ARO relationships. argNorm now also utilizes these three relationships.

#### Updated to modern build systems

Now, using `uv` for managing virtual environments (including on the CI) and `pixi` for installing all the tools for creating the databases.

#### Added testing for Python 3.13

Python 3.13 is now officially supported and tested on the CI.

## 0.6.0 - 21 August 2024

### User-facing changes

#### GROOT support
- argNorm supports [GROOT v1.1.2](https://github.com/will-rowe/groot). Use the `groot` tool parameter with the `groot-db`, `groot-core-db`, `groot-argannot`, `groot-card`, and `groot-resfinder` `db` parameters on the command line. The `GrootNormalizer` is available for use in Python scripts.
- `__version__` attribute added to the package (accessible as `argnorm.__version__` or `argnorm.lib.__version__`)
- atomic writing is used for outputs (https://github.com/untitaker/python-atomicwrites/tree/master). When everything is alright, this should not make a difference, but it avoids certain categories of failures.

### funcscan integration
- argNorm has been included as an nf-core module: https://nf-co.re/modules/argnorm/
- argNorm will also be available on the funcscan pipeline: https://github.com/nf-core/funcscan/pull/410

### Internal changes
- SARG db link was updated as old link is down
- RGI outputs in `crude_db_harmonisation` are concatenated so frequencies of `perfect`, `strict`, and `loose` hits can be calculated from concatenated file 

## 0.5.0 - 25 June 2024

### User-facing changes

### Improved drug categorization
- `drugs_to_drug_classes()` also uses the 'has_part' ARO relationship now to get drug classes for antibiotic mixtures. In case of antibiotic mixtures, the drug classes of the drugs associated with 'has_part' are returned rather than 'antibiotic mixture' (ARO:3000707).
- 'antibiotic mixture' will not be reported as a drug class, rather the individual antibiotic classes making up the antibiotic mixture will be reported.

#### Improved curation
- **manual curation (argannot)**: `(Tet)tetH:EF460464:6286-7839:1554` was incorrectly annotated as ARO:3004797 which is a beta-lactamase due to a loose RGI hit. This was manually curated to ARO:3000175.
- **Improved curation**:
    - resfinder_curation: grdA_1_QJX10702 -> 3007380 & EstDL136_1_JN242251 -> 3000557
    - megares_curation: MEG_2865|Drugs|Phenicol|Chloramphenicol_hydrolase|ESTD -> 3000557

#### Bugfixes
- `confers_resistance_to()` now gets drugs information even if it is encoded at a higher level in the ARO. For example, OXA-19 previously only returned cephalosporin and penam, but now will also return oxacillin (from AMR gene family).
- `drugs_to_drug_classes()` now correctly only returns the immediate child of 'antibiotic molecule' as the drug class (this was previously not the case for certain corner cases).
- **inconsistent ARO versions** deeparg, megares, resfinderfg & sarg curation: ARO:3004445 -> ARO:3005440, this was due to a change in the ARO and the ARO number for the RSA2 gene changing, but the version of ARO bundled with argNorm was out of sync.

### Internal changes

- AROs were previously handled as integers in the `get_aro_mapping_table()` function and this posed challenges when ARO numbers such as 'ARO:0010004' (leading zeros leading to issues). To fix this, AROs are now treated as strings so leading zeros can be maintained.

## 0.4.0 - 10 June

### User-facing changes
- Command line tool accept database/tool names in case-independent way (by @sebastianLedzianowski)
- A few ARO mappings were missing in the manual curation and they have been added. See the [CHANGELOG](https://github.com/BigDataBiology/argNorm/blob/main/CHANGELOG.md) for more details.
- `lib.map_to_aro` returns `None` if there is no mapping (raises an exception if the name is missing)

### Internal changes

- Bundle a specific version of ARO with the package instead of downloading it from the internet (ensures reproducibility)

## 0.3.0 - 27 April 2024

### User-facing changes
- Improved Resfinder mappings (40 gene clusters and 9 reverse complements were manually curated)
- Updated ARG-ANNOT mappings (a total of 10 mappings changed)
- argNorm is now more usable as a library
- Remove warnings

### Internal changes
- Code has been refactored to be simpler


## 0.2.0 - 26 March 2024

#### ARO Mapping & Normalization

- Updated mappings and manual curation tables for latest RGI
- Hamronized ResFinderFG support
- Removed python syntax in output

#### Drug Categorization

- Improved drug categorization by using superclasses whenever direct drug categorization is not possible
- Added better column headings for drug categorization (confers_resistance_to and resistance_to_drug_class)

#### Internal Improvements: Testing

- Improved pytest testing
- Added integration tests

## 0.1.0 - 20 December, 2023

- Added hamronized support for AMRFinderPlus
- Fixed ARO:nan issue (added manually curated mapping tables and integrated it with normalizers)
- Added drug categorization feature and integrated it with normalizers
- Added AMRFinderPlusNormalizer, ResFinderNormalizer
- Added specific smoke tests for ARGSOAPNormalizer, DeepARGNormalizer, AbricateNormalizer, AMRFinderPlusNormalizer and ResFinderNormalizer

## 0.0.1 - 13 June, 2022

- First release
- Initial source code started
- Normalizers: added BaseNormalizer, ARGSOAPNormalizer, DeepARGNormalizer, AbricateNormalizer
- Testing: added basic ARO column test

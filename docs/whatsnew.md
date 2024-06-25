# What's New

ArgNorm is under active development and outputs may change between versions. Particular versions of argNorm should, however, be reproducible as the databases are bundled with the tool.

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

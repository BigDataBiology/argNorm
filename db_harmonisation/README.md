# db_harmonisation

This is useful to _build_ the argNorm mapping databases, it is **not** intended for regular users of argNorm.

# Mapping Databases to RGI

Applicable to all databases except MEGARes

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

## Constructing MEGARes Mappings

> `construct_megares_mapping.py` is integrated with `crude_db_harmonisation.py` so `crude_db_harmonisation.py` will get megares mappings as well.

MEGARes is composed of genes from:
- CARD
- NCBI
- AMRFinder
- ARG-ANNOT
- BacMet: http://bacmet.biomedicine.gu.se/, https://academic.oup.com/nar/article/42/D1/D737/1064363
- Lahey
- ResFinder
- PointFinder
- rRNACARD (same as CARD, only 2 entries have this)

```
rRNACARD,>rRNACARD|gb|NR_103957.1|+|0-2910|ARO:3004836|Ngon_23S_AZM,MEG_8254|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation,

rRNACARD,>rRNACARD|gb|AE017221.1|+|1534489-1537382|ARO:3005083|Tthe_23S_PLM,MEG_8443|Drugs|Pleuromutilin|Pleuromutilin-resistant_23S_rRNA_mutation|P23S|RequiresSNPConfirmation,
```

> BacMet is ignored as metal and biocide resistance is not supported by CARD.

Megares to external header mappings file has one to many mappings!
```
- ARG-ANNOT,(Bla)SHV-1:AF148850:6-866:861,MEG_6240|Drugs|betalactams|Class_A_betalactamases|SHV
- ARG-ANNOT,(Bla)SHV-1:AF148850:6-866:861,MEG_6241|Drugs|betalactams|Class_A_betalactamases|SHV
&
- ARG-ANNOT,(Bla)OXY1-5:AY077483:181-1056:876,MEG_5284|Drugs|betalactams|Class_A_betalactamases|OXY
- ARG-ANNOT,(Bla)OXY1-5:AY077483:181-1056:876,MEG_5285|Drugs|betalactams|Class_A_betalactamases|OXY
```

Megares to external header mappings file has CARD mappings which can be directly used. However, CARD mappings have been run through RGI as there are cases where the ARO mapping for a gene from CARD is not even in the ARO:

```
CARD,CARD|gb|HM560971.1|-|119071-119269|ARO:3004468|CrpP [Pseudomonas aeruginosa] ,MEG_2133|Drugs|Fluoroquinolones|Fluoroquinolone_resistance_phosphotransferase|CRPP,

Here, ARO:3004468 does not even exist in the ARO.
```

Resfinder & argannot mappings are directly looked up using argNorm mappings. CARD, NCBI, AMRFinder, Lahey & PointFinder are ran through RGI. Mappings for resfinder & argannot that can't be found on argNorm's mappings are also run through RGI.

```
ResFinder,sitABCD_1_AY598030,MEG_8700|Multi-compound|Biocide_and_metal_resistance|Biocide_and_metal_ABC_efflux_pumps|SITABCD,

The Megares to external header mappings file states this is in resfinder, however resfinder does not have biocide & metal resistance entries and sitABCD_1_AY598030 can't be found in resfinder's mappings.
```

If a sequence translates clean or its reverse complement translates clean, the sequence is considered a CDS and the translated amino acid sequences are run through RGI using `protein` mode. Everything else is run through RGI using `contig` mode.

The RGI outputs of CDSs & contigs are combined with resfinder & argannot mappings to generate megares mappings.

> Note: as BacMet is now ignored, there will be less megares mappings (compared to when all of megares was run through RGI using `contig` mode) as loose hits from RGI won't be associated with BacMet entries.
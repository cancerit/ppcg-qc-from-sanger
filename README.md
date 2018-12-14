# ppcg-qc-from-sanger

[![Docker Repository on Quay](https://quay.io/repository/wtsicgp/ppcg-qc-from-sanger/status "Docker Repository on Quay")](https://quay.io/repository/wtsicgp/ppcg-qc-from-sanger)

ppcg-qc-from-sanger is used to extract PPCG defined QC metrics from Sanger Variant Calling pipeline results.

## Inputs

### `-tb`/`--tumour_bas`: Tumour sample BAS file(s), `-nb`/`--normal_bas`: Normal sample BAS file(s)

  A BAS file is one of the outputs of Sanger mapping pileline (cgpmap), which contains metrics collected during mapping and it has `.bas` suffix.

  The option can take files and directories as input. The option needs to be repeated for each entry if there're multiple entries. When input is a directory, only files with '.bam.bas' extension at the directory ***root*** will be taken. Subfolders are ignored.

### `-rt`/`--variant_call_tar`: Vaiant calling result tar file(s)

  The result tar file of Sanger Variant Calling pipeline (cgpwgs Version 1.1.2), and/or directories containing those files. Again, only a directory ***root*** is searched. If a `tar.gz` file is found in input and it's not a valid Sanger Variant Calling pipeline output, the tool will exit non-zero.

  The option needs to be repeated for each entry if there're multiple entries.

  ppcg-qc-from-sanger uses valid tar files to figure out which samples' BAS files are required.

### `-o`/`--output_tar`: Output tar file

  The output file. It has to have a `.tar.gz` extension. If it exists when the ppcg-qc-from-sanger starts to run, it'll exit non-zero.

### `-mt`/`--metadata`: Metadata file(s) (optional)

  Ideally this will be the tsv files generated by `validate_sample_meta.pl` of [cgpNgsQc](https://github.com/cancerit/cgpNgsQc) when validating BAM's metadata, with two extra columns: `Sequencer` and `Sequencing year`,  but you can also generate this file yourself.

  In order to have a full output, your meta files should have the following columns: `Donor ID`, `Donor UUID`, `Sample ID`, `Sample UUID`, `Sequencer` and `Sequencing Year`. Either `Sample ID` or `Sample UUID` should match the sample name in its BAM file header, as ppcg-qc-from-sanger will search samples names extracted from BAS files (e.g. sample name in a BAM file header) in the two colums to map samples' metadata to their QC metrics. ppcg-qc-from-sanger will populate 'NA's in the columns if their coresponding metadata is not found.

  The option can also take directories as input, in which case all `.tsv` files at the directories' ***root*** will be taken as metadata. However if same `Sample_ID` or `Sample_UUID` is found more than once in the files, ppcg-qc-from-sanger will exit non-zero. The option needs to be repeated for each entry if there're multiple entries.

### `-cv`/`--count_variants`: Count variants (optional)

  This flag is default to false. If specified, ppcg-qc-from-sanger will also count number of variants called by Sanger Variant Calling pipeline.

### `-gs`/`--genome_size`: Genome size (optional)

  It's required for calculating sequencing depth. Default to the sum of length of GRCh37 chromosomes.

## Outputs

Output of ppcg-qc-from-sanger is a tar file containing the following files:

### `ppcg_sanger_metrics.txt`

A tsv file containing the following columns:

1. **Tumour sample name**: Tumour sample name found in the BAM file header.
1. **Tumour sample ID**: `Sample_ID` in the metadata, if provided.
1. **Tumour sample UUID**: `Sample_UUID` in the metadata, if provided.
1. **Tumour sequencing year**: `Sequencing_Year` in the metadata, if provided.
1. **Tumour sequencer**: The sequencer on which the sample was sequenced, if provided in the metadata input in a column named `Sequencer`.
1. **Tumour ReadGroup IDs**: ReadGroup (RG) ID found in the sample BAM.
1. **Tumour depth per RG**: Total mapped bases (number of mapped reads * read length) divided by the genome size per RG.
1. **Tumour total depth**: Total depth.
1. **Tumour fraction of mapped reads per RG**: Number of mapped reads divided by the number of total reads per RG.
1. **Tumour median fraction of mapped reads**: Median fraction of mapped reads for the sample.
1. **Tumour insert size per RG**:  Mean insert size per RG.
1. **Tumour median Insert size**: Median of mean insert size.
1. **Tumour insert size sd per RG**: Insert size standard deviation per RG.
1. **Tumour median insert size sd**: Median of the insert size standard deviations.
1. **Tumour r1 GC content per RG**: GC base content of all read 1 in read pairs per RG.
1. **Tumour median r1 GC content**: Median r1 GC base content for the sample.
1. **Tumour r2 GC content per RG**: GC base content of all read 2 in read pairs per RG.
1. **Tumour median r2 GC content**: Median r2 GC base content for the sample.
1. **Tumour fraction of duplicated reads per RG**: Number of duplicated reads divided by the number of total reads per RG.
1. **Tumour median fraction of duplicated reads**: Median fraction of duplicated reads for the sample.
1. **Tumour fraction of mis-matched pairs per RG**: 1 minus the fraction of mapped pairs in total read pairs per RG.
1. **Tumour median fraction of mis-matched pairs**: Median mismatch fraction.
1. **Tumour contamination per RG**: Cross individual contamination per RG.
1. **Tumour median contamination**: Median cross individual contamination for the sample.
1. **Tumour sex**: Deduced sex from genotypes of 4 SNPs on sex chromosomes.
1. **Tumour fraction of matched sex with Normal**: Fraction of the 4 SNPs on sex chromosomes that have matched genotypes in the normal sample.
1. **Tumour fraction of matched genotype with Normal**: Fraction of 92 autosome SNPs that have matched genotypes in the normal sample
1. **Normal contamination in Tumour**: Estimate of fraction of normal cells contaminating the tumour sample.
1. **Normal sample name**:
1. **Normal sample ID**:
1. **Normal sample UUID**:
1. **Normal sequencing year**:
1. **Normal sequencer**:
1. **Normal ReadGroup IDs**:
1. **Normal depth per RG**:
1. **Normal total depth**:
1. **Normal fraction of mapped reads per RG**:
1. **Normal median fraction of mapped reads**:
1. **Normal insert size per RG**:
1. **Normal median Insert size**:
1. **Normal insert size sd per RG**:
1. **Normal median insert size sd**:
1. **Normal r1 GC content per RG**:
1. **Normal median r1 GC content**:
1. **Normal r2 GC content per RG**:
1. **Normal median r2 GC content**:
1. **Normal fraction of duplicated reads per RG**:
1. **Normal median fraction of duplicated reads**:
1. **Normal fraction of mis-matched pairs per RG**:
1. **Normal median fraction of mis-matched pairs**:
1. **Normal contamination per RG**:
1. **Normal median contamination**:
1. **Donor ID**: `Donor_ID` in the metadata, if provided.
1. **Donor UUID**: `Donor_UUID` in the metadata, if provided.

Additional columns if uses `--count_variants` flag:

1. **Number of SNVs** (tumour only): the number of 'PASS'ed variants in the CAVEMAN flagged VCF output `/` all variants in the file
1. **SNVs per Chr (chromosome:[PASS, ALL])** (tumour only): JSON format string: keys are chromosome names, values are lists consisting of the number of 'PASS'ed SNVs and the number of all SNVs.
1. **Number of INDELs** (tumour only): the number of 'PASS'ed variants in the PINDEL flagged VCF output `/` all variants in the file
1. **INDELs per Chr (chromosome:[PASS, ALL])** (tumour only): JSON format string: keys are chromosome names, values are lists consisting of the number of 'PASS'ed INDELs and the number of all INDELs.
1. **Number of SVs** (tumour only): the number of variants with 'BAS' info tag in the BRASS VCF output, devided by 2 `/` all variants in the file devided by 2 *(in the VCF, one SV is recorded with two rows)*
1. **SVs per Chr (chromosome:[PASS, ALL])** (tumour only): JSON format string: keys are chromosome names, values are lists consisting of the number of variants with 'BAS' info tag in the BRASS VCF output, devided by 2 and all variants in the file devided by 2.
1. **Number of CNVs** (tumour only): number of variants which have different copy numbers between Normal and Tumour sample in the ASCAT VCF output

### `.tsv` genotyping files

Four files containing the genotyping information of the two samples (the tumour and the normal):

* **{tumour_sample}.full_gender.tsv** and **{normal_sample}.full_gender.tsv**

  Genotypes of 4 SNPs on sex chromosomes.

* **{tumour_sample}.full_genotype.tsv** and **{normal_sample}.full_genotype.tsv**

  Genotypes of 92 SNPs on autosomes.

## INSTALL

Installation is via `pip`.  Simply execute with the path to the packaged distribution:

```bash
pip install -pip install https://github.com/cancerit/ppcg-qc-from-sanger/archive/master.tar.gz
```

## Development environment

### Setup VirtualEnv

```bash
cd $PROJECTROOT
hash virtualenv || pip3 install virtualenv
virtualenv -p python3 env
source env/bin/activate
pip install -r requirements.txt
python setup.py develop # so bin scripts can find module

## If changed requirements please run:
pip freeze | grep -v `echo ${PWD##*/}` > requirements.txt
```

For testing/coverage (`./run_tests.sh`)

```bash
source env/bin/activate # if not already in env
pip install pytest
pip install pytest-cov
pip install pep8
pip install radon
gem install --user-install mdl
pip install requests_mock
```

Test that `mdl` is available, if not add the following to your path variable:

```bash
export PATH=$HOME/.gem/ruby/X.X.X/bin:$PATH
```

## LICENCE

Copyright (c) 2018 Genome Research Ltd.

Author: CancerIT <cgpit@sanger.ac.uk>

This file is part of ppcg-qc-from-sanger.

ppcg-qc-from-sanger is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

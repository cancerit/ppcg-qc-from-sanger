# ppcg-qc-from-sanger

The tools is used to extract PPCG defined QC metrics from Sanger variant calling pipeline results.

## Inputs

### tumour and normal sample bas files

  BAS file is one of the output of Sanger mapping pileline (cgpmap), which contains metrics collected during mapping and it has `.bas` suffix.

### vaiant calling result tar file

  This is the result tar file of Sanger variant calling pipeline (cgpwgs Version 1.1.2), it should have `.tar.gz` suffix if named correctly.

### Genome size

  It's required for calculating sequencing depth. Default to the sum of chromosome lengths in GRCh37.

## Outputs

### File: ppcg_metrics_from_cgpwgs.txt

A tsv file containing the following columns:

1. **Sample name**: sample name
1. **Sample Type**: `Tumour` or `Normal`
1. **Depth**: total mapped bases divided by the genome size
1. **Fraction of mapped reads**: number of mapped reads divided by the number of total reads in the sample
1. **Insert size**: mean of mean insert size of lanes
1. **Insert size sd**: mean of insert size standard division of lanes
1. **r1 GC content**: mean GC base content of all read 1 in read pairs
1. **r2 GC content**: mean GC base content of all read 2 in read pairs
1. **Fraction of duplicated reads**: number of mapped pairs divided by the number of total reads
1. **Fraction of mis-matched pair**: number of mismatched pairs divided by the number of total read pairs
1. **Contamination**: contamination (SNP based)
1. **Gender** (tumour only): deduced gender from genotypes of 4 SNPs on sex chromosomes
1. **Fraction of matched gender** (tumour only): fraction of the 4 SNPs on sex chromosomes that have matched genotypes in the normal sample
1. **Fraction of matched genotype** (tumour only): fraction of 92 autosome SNPs that have matched genotypes in the normal sample
1. **Normal Contamination** (tumour only): normal sample contamination in the tumour sample

### Folder: genotypes_from_cgpwgs

In this folder, there are four files containing the genotyping information of the two samples (the tumour and the normal).

* **{tumour_sample}.full_gender.tsv** and **{normal_sample}.full_gender.tsv**

  Genotypes of 4 SNPs on sex chromosomes.

* **{tumour_sample}.full_genotype.tsv** and **{normal_sample}.full_genotype.tsv**

  Genotypes of 92 SNPs on autosomes.

## INSTALL

Installation is via `pip`.  Simply execute with the path to the packaged distribution:

```bash
pip install --find-links=~/wheels ppcg-qc-from-sanger
```

## Development environment

### Setup VirtualEnv

```
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

```
source env/bin/activate # if not already in env
pip install pytest
pip install pytest-cov
pip install pep8
pip install radon
gem install --user-install mdl
pip install requests_mock
```

Test that `mdl` is available, if not add the following to your path variable:

```
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

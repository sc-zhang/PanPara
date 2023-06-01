## Introduction
PanPara is a tool for identifying paralog genes with pan-genome.

## Dependencies
### Software
* ncbi-blast+
* MCScanX

## Usage
```bash
usage: panpara.py [-h] -l LIST -s CDS -b BED [-d IDEN] [-c COVERAGE] -o OUTPUT [-t THREAD]

options:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  list file, each row contain one sample name, first row means reference
  -s CDS, --cds CDS     cds directory, must exists all samples end with ".cds", for example: "sample1.cds"
  -b BED, --bed BED     bed directory, must exists all samples end with ".bed" , for example: "sample1.bed"
  -d IDEN, --iden IDEN  identity threshold, default=0.8
  -c COVERAGE, --coverage COVERAGE
                        the threshold of alignment coverage, default=0.8
  -o OUTPUT, --output OUTPUT
                        output directory
  -t THREAD, --thread THREAD
                        threads for running some steps, default=6
```

## Output
**final.csv** contain all indentified paralog genes.
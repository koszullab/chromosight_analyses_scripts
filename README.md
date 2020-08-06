## Codes and functions for bioanalysis of Chromosight article ##

This page presents the different codes and functions developed to perform the bioanalyses described in the article **Computer vision for pattern detection in chromosome contact maps** by Cyril Matthey-Doret et al. The codes presented here should allow to reproduce the different plots and analyses from the main text and the supplementary data. 

### Table of contents

* [Dependencies](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#raw-data-extraction-and-alignment)
* [Building contacts map](https://github.com/koszullab/E.coli.analysis/blob/master/README.md#building-contacts-map)


Scripts and codes can be run on OS X and other Unix-based systems, and necessitate:
#### *Python (>=3.)*
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `SRA Tool kit ` / [sra_sdk/2.9.6](https://www.ncbi.nlm.nih.gov/books/NBK158900/)


## Raw data extraction and alignment
#### Data extraction

FASTQ files of the reads were deposited in the NCBI database under the GEO accession number GSE107301. A SRA executable called fastq-dump from SRA can be used to extract and split both reads of pair-end sequences: 
```bash
fasterq-dump --split-3 SRR1514669 -O .
```

#### Alignment

The aligment and construction of matrices have been performed as perviously described. 

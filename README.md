## Codes and functions for bioanalysis of Chromosight article ##

This page presents the different codes and functions developed to perform the bioanalyses described in the article **Computer vision for pattern detection in chromosome contact maps** by Cyril Matthey-Doret et al. The codes presented here should allow to reproduce the different plots and analyses from the main text and the supplementary data. 

Preprint can be found on https://www.biorxiv.org/content/10.1101/2020.03.08.981910v3.full

Docs of the algorithm available at https://chromosight.readthedocs.io

This github page is a companion page of the Chromosight algorithme page:
https://github.com/koszullab/chromosight

### Table of contents

* [Dependencies](https://github.com/koszullab/chromosight_codes_for_bioanalysis/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/koszullab/chromosight_codes_for_bioanalysis/blob/master/README.md#raw-data-extraction-and-alignment)
*


Scripts and codes can be run on OS X and other Unix-based systems, and necessitate:
#### *Python (>=3)*
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython
* pandas

#### External programs

* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `SRA Tool kit ` / [sra_sdk/2.9.6](https://www.ncbi.nlm.nih.gov/books/NBK158900/)


#### Data extraction

FASTQ files of the reads were deposited in the NCBI database under the GEO accession number GSE107301. A SRA executable called fastq-dump from SRA can be used to extract and split both reads of pair-end sequences: 
```bash
fasterq-dump --split-3 SRR1514669 -O .
```

#### Alignment 

The aligment and construction of matrices have been performed as perviously described. 



## Detection analysis

This document presents the command lines used for the detection analyses presented in the main text. 
Contact data as cool files can be dowloaded on zenodo [doi:  10.5281/zenodo.3742095](doi:  10.5281/zenodo.3742095)

version of Chromosight 1.1.2

## Fig2.a and b
### *Saccharomyces cerevisiae*

```chromosight detect --pattern=loops_small --min-dist 5000 --max-dist 200000  --perc-undetected=30 contacts2_AT198_AT199_2000_df.cool  out_1.1_AT198_AT199_pear0.5```

## Fig. 2.c
###  *Schizosaccharomyces pombe*
```chromosight detect --pattern=loops_small --pearson=0.4 --min-dist 5000 --max-dist 200000 contacts2_SRR5149256_2000_df.cool SRR5149256_Pombe_pear04```

## Fig. 2.d
### *Saccharomyces cerevisiae*
```chromosight detect --pattern=loops_small --min-dist 5000 --max-dist 200000 --perc-undetected=50 --perc-zero=10 contacts2_SRR8769554_HiC_SacCerW303_stop-alpha-factor-G1_2000_df.cool out_G1_SRR8769554_HiC_SacCerW303_pear05```

#### For comparison of groups of loops

```bash script_common_loops5.bh contacts2_AT198_AT199_2000_df.cool mitotic contacts2_SRR8769554_HiC_SacCerW303_stop-alpha-factor-G1_2000_df.cool G1```

## Fig. 3

### *Human*

```chromosight detect --pattern=loops_small  --min-dist=15000 --max-dist=2000000 contacts2_1_SRR6675327.cool out_contacts2_1_SRR6675327_pear05```

```chromosight detect --pattern=borders --pearson=0.4   contacts2_1_SRR6675327.cool out_contacts2_1_SRR6675327_borders_04```

```chromosight detect --pattern=hairpins --pearson=0.4  contacts2_1_SRR6675327.cool out_contacts2_1_SRR6675327_hairpins_04```


### *Bacillus subtilis*

```chromosight detect --pearson=0.3 --min-dist 12000 --max-dist 4300000 --perc-undetected=30 contacts2_df_2000_SRR2214080.cool out_bacillus_03```

### *Epstein Barr Virus (EBV)*

```chromosight detect --pattern=loops --min-dist=7000 --perc-undetected=30  --perc-zero=30 contacts2_SRR2312566_RAW_Epstein_500_500_df.cool  out_SRR2312566_epstein```









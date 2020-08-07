## Main results from the chromosight paper

This repo contains instructions to reproduce the different figures and results from our paper, **[Computer vision for pattern detection in chromosome contact maps](https://www.biorxiv.org/content/10.1101/2020.03.08.981910v3.full)** by Cyril Matthey-Doret et al., which showcases its signature program, [chromosight](https://github.com/koszullab/chromosight). See also the [official documentation of the program](https://chromosight.readthedocs.io), complete with demos and tutorials for your own use cases.


### Table of contents

* Requirements
* [Raw data extraction and alignment](https://github.com/koszullab/chromosight_codes_for_bioanalysis/blob/master/README.md#raw-data-extraction-and-alignment)

#### Requirements

##### Environment

Chromosight is designed to run on UNIX-like environments (e.g. Linux, OS X, Windows Subsystems for Linux, etc.) and has been tested on Ubuntu >= 14.10. It is entirely written in Python 3.6 with no planned support for earlier versions, including Python 2. The version used in the preprint is 1.1.2 and is installable from PyPI:

```sh
python3 -m pip install --user chromosight==1.1.2
```

##### Python packages

If for some reason you need to install chromosight manually, then the following packages are required (from requirements.txt):

* cooler
* docopt
* jsonschema
* matplotlib
* numpy
* scikit-learn
* scipy>= 1.3

##### External programs

In order to pull Hi-C data from the SRA and perform the alignment to generate the contact maps, you will need the following programs (directly available on the  Ubuntu repositories):

* [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [`sra-toolkit`](https://www.ncbi.nlm.nih.gov/books/NBK158900/) (`sratoolkit` on homebrew)

#### Data extraction

FASTQ reads are available under the GEO accession number GSE107301. They can be pulled and split pairwise with the following command: 

```bash
fasterq-dump --split-3 SRR1514669 -O .
```

## Detection analysis
 
Contact data as cool files can be dowloaded on zenodo [doi:10.5281/zenodo.3742095](https://zenodo.org/record/3742095)

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

### *Candida Albicans*

TO DO 


#### Computation and visualisation of Loop spectrum

From cool files, we first quantified the loop signals for different pairs of cohesin peaks using quantification mode of chromosight: 

```chromosight quantify --pattern=loops --perc-zero=100 --perc-undetected=100 pairs_intra_groupe_Gm12878Rad21V0416101UniPk.txt.reordered.1Mb contacts2_1_SRR6675327.cool contacts2_1_SRR6675327.cool.quantified```

We then use the python code [spectrum_cycle2_imple.py](https://github.com/koszullab/chromosight_analyses_scripts/blob/master/python_codes/spectrum_cycle2_imple.py) to computute the lowess signal from the scatter plot of the loop scores versus distances between 
peaks of cohesin and create the corresponding plot. 







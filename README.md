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
* [`bedtools`](https://github.com/arq5x/bedtools2) (`bedtools`)
* [`cooler`](https://github.com/mirnylab/cooler) (`cooler`)

#### Data extraction

FASTQ reads are available on SRA servers. They can be pulled and split pairwise with the following command: 

```bash
fasterq-dump --split-3 SRR1514669 -O .
```

## Detection analysis
 
Contact data as cool files can be dowloaded on zenodo [doi:10.5281/zenodo.3742095](https://zenodo.org/record/3742095)

## Fig2.
### *Saccharomyces cerevisiae*

```chromosight detect --pattern=loops_small --min-dist 5000 --max-dist 200000  --perc-undetected=30 SRR7706227_SRR7706226_hic_scer_mitotic_2kb.cool SRR7706227_SRR7706226_scer_mitotic_p0.5```

###  *Schizosaccharomyces pombe*
```chromosight detect --pattern=loops_small --pearson=0.4 --min-dist 5000 --max-dist 200000 SRR5149256_hic_spo_2kb.cool SRR5149256_spo_loops_p04```

### *Saccharomyces cerevisiae*
```chromosight detect --pattern=loops_small --min-dist 5000 --max-dist 200000 --perc-undetected=50 --perc-zero=10 SRR8769554_hic_scer_g1_2kb.cool SRR8769554_scer_g1_loops_p05```

#### For comparison of groups of loops

```bash script_common_loops5.bh SRR7706227_SRR7706226_hic_scer_mitotic_2kb.cool mitotic SRR8769554_hic_scer_g1_2kb.cool G1```

## Fig. 3

### *Human*

```chromosight detect --pattern=loops_small  --min-dist=15000 --max-dist=2000000 SRR6675327_hic_hsap_GM12878_10kb.cool SRR6675327_GM12878_loops_p05```

```chromosight detect --pattern=borders --pearson=0.4   SRR6675327_hic_hsap_GM12878_10kb.cool SRR6675327_GM12878_borders_p04```

```chromosight detect --pattern=hairpins --pearson=0.4  SRR6675327_hic_hsap_GM12878_10kb.cool SRR6675327_GM12878_hairpins_04```


### *Bacillus subtilis*

```chromosight detect --pearson=0.3 --min-dist 12000 --max-dist 4300000 --perc-undetected=30 SRR2214080_3cseq_bsub_2kb.cool out_bsub_loops_p03```

### *Epstein Barr Virus (EBV)*

```chromosight detect --pattern=loops --min-dist=7000 --perc-undetected=30  --perc-zero=30 SRR2312566_chiapet_ebv_500bp.cool SRR2312566_ebv_loops```

## Fig. 4

### *Candida Albicans*

```bash

sc_map="SRR7706226_hic_scer_mitotic.mcool::/resolutions/10000"
ca_map="SRR3381672_hic_calb_10kb.cool"

# Interactively pick 5-10 centros from S. cerevisiae
# Note the grep and sed commands capture the standard output (coords of clicks)
#and send them to a file
chromosight generate-config --n-mads 15 \
                            --click "$sc_map" \
                            --win-size 41 \
                            -e centromeres sc_centro \
  | grep "^x = " \
  | sed 's/[xy] \?= //g' \
  | sed 's/, /\t/' \
  > input_coords.tsv

# De novo centromere detection on C. albicans using kernel built from S. cerevisiae
chromosight detect -p0.7 \
                  --n-mads 15 \
                  -t12 \
                  --kernel-config sc_centro.json \
                  --inter "$ca_map" \
                  sc_to_ca_centro

# Plot both Hi-C maps. Scer with the click coordinates and Calb with the detected regions
python plot_sc_ca_maps.py "$sc_map" \
                          "$ca_map" \
                          input_coords.tsv \
                          sc_to_ca_centro.tsv

```

The script `plot_sc_ca_maps.py` is in the `python_codes` directory.

## Supplementary figures

### Computation and visualisation of Loop spectrum

To generate possible pairs from bed files you can use:

```bash
MINDIST=10000
MAXDIST=1000000
bedtools window -a rad21_hg38_top4000.bed \
                -b rad21_hg38_top4000.bed \
                -w $MAXDIST \
    | awk -vmd=$MINDIST '$1 == $4 && ($5 - $2) >= md {print}' \
    | sort -k1,1 -k2,2n -k4,4 -k5,5n \
    > rad21_hg38_top4000_comb_10kb_1mb.bed2d
```

From cool files, we work with 10 kb resolution:
```cooler coarsen 4DNFIMH3J7RW.mcool::/resolutions/10000 -o 4DNFIMH3J7RW.mcool.10000```

We then quantified the loop signals for different pairs of cohesin peaks using quantification mode of chromosight: 

```chromosight quantify --pattern=loops --perc-zero=100 --perc-undetected=100 rad21_hg38_comb_10kb_1mb.bed2d 4DNFIMH3J7RW.mcool.10000 4DNFIMH3J7RW.mcool.10000.quantified```

We then use the python code [spectrum_cycle2_imple.py](https://github.com/koszullab/chromosight_analyses_scripts/blob/master/python_codes/spectrum_cycle2_imple.py) to compute the lowess signal from the scatter plot of the loop scores versus distances between 
peaks of cohesin and create the corresponding plot. 



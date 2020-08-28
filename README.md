# Codes developped for the analyses presented in the Chromosight paper

This repository contains instructions to reproduce the different figures and results from our paper, **[Computer vision for pattern detection in chromosome contact maps](https://www.biorxiv.org/content/10.1101/2020.03.08.981910v3.full)** by Cyril Matthey-Doret et al., which showcases its signature program, [chromosight](https://github.com/koszullab/chromosight). See also the [official documentation of the program](https://chromosight.readthedocs.io), complete with demos and tutorials for your own use cases.

## Table of contents

* [Requirements](https://github.com/koszullab/chromosight_analyses_scripts/blob/master/README.md#requirements)
* [Raw data extraction and alignment](https://github.com/koszullab/chromosight_analyses_scripts/blob/master/README.md#data-extraction)
* [Generation of simulated data](https://github.com/koszullab/chromosight_analyses_scripts/blob/master/README.md#generation-of-simulated-data)
* [Analyses](https://github.com/koszullab/chromosight_analyses_scripts/blob/master/README.md#analyses)
* [Supplementary figures](https://github.com/koszullab/chromosight_analyses_scripts/blob/master/README.md#supplementary-figures)

## Requirements

### Environment

Chromosight is designed to run on UNIX-like environments (e.g. Linux, OS X, Windows Subsystems for Linux, etc.) and has been tested on Ubuntu >= 14.10. It is entirely written in Python 3.6 with no planned support for earlier versions, including Python 2. The version used in the preprint is 1.1.2 and is installable from PyPI or conda:

```sh
python3 -m pip install --user chromosight
```
or
```sh
conda install -c bioconda chromosight
```

### Python packages

If for some reason you need to install chromosight manually (without using pip or conda), then the following packages are required (from requirements.txt):

* cooler
* docopt
* jsonschema
* matplotlib
* numpy
* scikit-learn
* scipy>= 1.3

### External programs

In order to pull Hi-C data from the SRA and perform the alignment to generate the contact maps, you will need the following programs (directly available on the  Ubuntu repositories):

* [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [`sra-toolkit`](https://www.ncbi.nlm.nih.gov/books/NBK158900/) (`sratoolkit` on homebrew)
* [`bedtools`](https://github.com/arq5x/bedtools2) (`bedtools`)
* [`cooler`](https://github.com/mirnylab/cooler) (`cooler`)

## Data extraction

FASTQ reads are available on SRA servers. They can be pulled and split pairwise with the following command: 

```bash
fasterq-dump --split-3 SRR1514669 -O .
```

## Generation of simulated data
In order to carry out a test bench of our algorithm and also to be able to compare its performance with other methods, we generated simulated contact maps containing loops and domains. These contact maps are based on the statistical features of contact maps of S.Cerevisiae yeast. 
The principle is that the general characteristics of the experimental map should be found in the simulated maps. For this, the frequency of contact as a function of genomic distance, the number of chromosomal domains, the presence of loops must be as close as possible to the experimental data so that it is difficult for an eye to tell the difference. 
The algorithm test_simulated_maps8.py calculates these different signals from experimental and recreate data with similar properties. 
The loops and borders are positioned by multiplying by the corresponding pattern. A contact probability map is constructed from which the contacts are drawn respecting the number of reads of the experiment.


## Analyses

This section contains the command lines and scripts used to generate figures from the manuscript.
Contact data (cool files) and BED files used in the manuscript can be dowloaded on zenodo [doi:10.5281/zenodo.3742094](https://doi.org/10.5281/zenodo.3742094)

### Fig. 1

Benchmark input and output data, as well as the code used, are archived and documented in the [Zenodo entry](https://doi.org/10.5281/zenodo.3742094).

### Fig. 2

Generally speaking, the chromosight parameters used for loop detection in figure 2 are the following:

**Saccharomyces cerevisiae, M phase**

```bash
chromosight detect --pattern=loops_small \
                   --min-dist 5000 \
                   --max-dist 200000 \
                   --perc-undetected=10 \
                   SRR7706227_SRR7706226_hic_scer_mitotic_2kb.cool \
                   SRR7706227_SRR7706226_scer_mitotic_p05
```

### *Saccharomyces cerevisiae, G1 phase*
```bash
chromosight detect --pattern=loops_small \
                   --min-dist 5000 \
                   --max-dist 200000 \
                   --perc-undetected=50 \
                   --perc-zero=10 \
                   SRR8769554_hic_scer_g1_2kb.cool \
                   SRR8769554_scer_g1_loops_p05
```

**Schizosaccharomyces pombe**

```bash
chromosight detect --pattern=loops_small \
                   --pearson=0.4 \
                   --min-dist 5000 \
                   --max-dist 200000 \
                   SRR5149256_hic_spo_2kb.cool \
                   SRR5149256_spo_loops_p04
```

For comparison of the number of loops detected in G1 and M phase, both matrices have been subsampled to the same number of contacts. The whole process, along with the generation of the Venn diagram in figure 2c is streamlined in bash script `common_loops_fig2.sh`, which can be used as follows:

`bash detect_loops_venn2.sh SRR7706227_SRR7706226_hic_scer_mitotic_2kb.cool mitotic SRR8769554_hic_scer_g1_2kb.cool G1`

### Fig. 3

**Human**

Loop detection in human matrices were done with the following commands:

```bash
chromosight detect --pattern=loops_small \
                   --min-dist=15000 \
                   --max-dist=2000000 \
                   SRR6675327_hic_hsap_GM12878_10kb.cool \
                   SRR6675327_GM12878_loops_p05

chromosight detect --pattern=borders \
                   --pearson=0.4 \
                   SRR6675327_hic_hsap_GM12878_10kb.cool \
                   SRR6675327_GM12878_borders_p04

chromosight detect --pattern=hairpins \
                   --pearson=0.4 \
                   SRR6675327_hic_hsap_GM12878_10kb.cool \
                   SRR6675327_GM12878_hairpins_04
```


**Bacillus subtilis**

```bash
chromosight detect --pearson=0.3 \
                   --min-dist 12000 \
                   --max-dist 4300000 \
                   --perc-undetected=30 \
                   SRR2214080_3cseq_bsub_2kb.cool \
                   out_bsub_loops_p03
```

**Epstein Barr Virus (EBV)**

```bash
chromosight detect --pattern=loops \
                   --min-dist=7000 \
                   --perc-undetected=30  \
                   --perc-zero=30 \
                   SRR2312566_chiapet_ebv_500bp.cool \
                   SRR2312566_ebv_loops
```

### Fig. 4

Below are the commands required to interactively pick inter-telomeric coordinates from the _S. cerevisiae_ contact map, generate for custom kernel generation, ans use the resulting kernel for detection on the _C. albicans_ contact map. Both contact maps are available in the zenodo entry.

**Candida Albicans**

```bash

sc_map="SRR7706226_hic_scer_mitotic.mcool::/resolutions/10000"
ca_map="SRR3381672_hic_calb_10kb.cool"

# Interactively pick 5-10 centros from S. cerevisiae
# Note the grep and sed commands capture the standard output (coords of clicks)
# and send them to a file. They are only required to show each individual window on the figure
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

We then quantified the loop signals for different pairs of cohesin peaks using quantification mode of chromosight: 

```bash
chromosight quantify --pattern=loops \
                      --perc-zero=100 \
                      --perc-undetected=100 \
                      rad21_hg38_comb_10kb_1mb.bed2d \
                      4DNFIMH3J7RW.mcool::/resolutions/10000 4DNFIMH3J7RW_10kb_rad21_hg38_comb_10kb_1mb_loops
```

We then use the python code [spectrum_cycle2_imple.py](https://github.com/koszullab/chromosight_analyses_scripts/blob/master/python_codes/spectrum_cycle2_imple.py) to compute the lowess signal from the scatter plot of the loop scores versus distances between 
peaks of cohesin and create the corresponding plot. 

### Detection of loops in other protocols than Hi-C

```bash
chromosight detect --pattern=loops_small \
                   --threads=10 \
                   --min-dist=15000 \
                   --max-dist=2000000 \
                   4DNFI81RQ431.mcool::/resolutions/10000 \
                   out_4DNFI81RQ431_10kb_loops_small`
```

### Comparison of loop calls around CTCF across loop callers

Loops are detected on GSE63525_GM12878_insitu_primary, from Rao et al. 2014. The contact data and HiCCUPS calls were retrieved from the [GEO entry](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525).

> The contact data is only availabl in hic format, which can be converted to mcool using [hic2cool](https://github.com/4dn-dcic/hic2cool)

For other loop callers, loops were detected with default parameters at 10kb resolution. Output calls available in the `GSE63525_loop_calls` folder of the `processed_files.tar.gz` tarball of [the zenodo entry](https://doi.org/10.5281/zenodo.3742094):

```bash
COOL="GSE63525_GM12878_insitu_primary.mcool::/resolutions/10000"
hicDetectLoops --matrix "$COOL" --outFileName hicexplorer/hicexplorer_loops
cooltools compute-expected -p 4 -o cooltools/cooltools_expected.tsv "$COOL"
cooltools call-dots "$COOL" cooltools/cooltools_expected.tsv -o cooltools/cooltools_loops.tsv
chromosight detect "$COOL" chromosight/chromosight_loops_small_GSE63525
```
The output calls when then visualized using `python zoom_compare_calls.py`. The script is available in the `python_codes` directory.

### Comparison of loop calls genome wide across loop callers

> Note: all files required for this section are available in the `GSE63525_loop_calls` directory in the `processed_files.tar.gz` archive of the zenodo record.

The output from the 4 commands above were used. The script `zoom_compare_calls.py` also generates a coordinate file for each software named "<software_coords>". Those files all have the same format: two tab-separated columns corresponding to the row and column coordinates of loops detected with the software.
Thos coordinates files can be fed as input to the `common_patterns_upsetplot.py` script as follows:

> Note: This script requires the python package [upsetplot](https://upsetplot.readthedocs.io/en/latest/index.html), which can be installed using pip.

`python common_patterns_upsetplot.py chromosight/chromosight_coords cooltools/cooltools_coords hicexplorer/hicexplorer_coords hiccups/hiccups_coords`
This script finds the common loops found by every combination of softwares, allowing an arbitrary jitter for each pattern (+/- 1pixel by default).
It prints the table with loop counts for each combination, and display an upsetplot.

To compute the overlap between CTCF peaks and the loop calls from each software, we used the script `common_chip_patterns.py`. Note the first argument is used as the reference set (here, ChIP-seq):

```bash
python common_chip_patterns.py \
       ctcf_chipseq/ctcf_hg19_comb_10kb_2mb_coords \
       chromosight/chromosight_coords \
       hicexplorer/hicexplorer_coords \
       hiccups/hiccups_coords \
       cooltools/cooltools_coords
```



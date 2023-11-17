#!/bin/bash
# KournaK 
# script to compare groups of detected patterns in 2 conditions and generate a Venn diagram
# Ex:  bash common_loops_venn2.sh lib1.cool lib1 lib2.cool lib2

# Input parameters:

bank1_file="$1"
name_bank1="$2"
bank2_file="$3"
name_bank2="$4"

win_size=17

# Sub-sampling: 
N1=$(cooler info $bank1_file | grep sum | awk '{print int($NF)}')
N2=$(cooler info $bank2_file | grep sum | awk '{print int($NF)}')
Nmin=$((( $N1 <= $N2 )) && echo "$N1" || echo "$N2")
echo 'Number of reads for both banks'
echo $Nmin 
Nmin=$(($Nmin -1))
echo $Nmin 

# Chromosight de-novo detection on both matrices, subsampled to equal coverage:

chromosight detect --threads 8 \
                   --subsample=$Nmin \
                   --perc-undetected=50 \
                   --perc-zero=10 \
                   --min-dist 5000 \
                   --max-dist 200000 \
                   $bank1_file \
                   $name_bank1

chromosight detect --threads 8 \
                   --subsample=$Nmin \
                   --perc-undetected=50 \
                   --perc-zero=10 \
                   --min-dist 5000 \
                   --max-dist 200000 \
                   $bank2_file \
                   $name_bank2

# Computation of overlap of groups of detections: 

python common_loops_venn2.py "$name_bank1.tsv" "$name_bank1" "$name_bank2.tsv" "$name_bank2"

# Quantification in each sub-group and in each bank:
 
chromosight quantify "group1_detected_in_${name_bank1}.txt" \
                     --subsample=$Nmin \
                     --win-size=$win_size \
                     "$bank1_file" \
                     "out_in_${name_bank1}.group1.${name_bank1}"

chromosight quantify "group2_detected_in_${name_bank2}.txt" \
                     --subsample=$Nmin \
                     --win-size=$win_size \
                     "$bank1_file" \
                     "out_in_${name_bank1}.group2.${name_bank2}"

chromosight quantify "group12_detected_in_${name_bank1}_${name_bank2}.txt" \
                     --subsample=$Nmin \
                     --win-size=$win_size \
                     "$bank1_file" \
                     "out_in_${name_bank1}.group12"

chromosight quantify "group1_detected_in_${name_bank1}.txt" \
                     --subsample=$Nmin \
                     --win-size=$win_size \
                     "$bank2_file" \
                     "out_in_${name_bank2}.group1.${name_bank1}"

chromosight quantify "group2_detected_in_${name_bank2}.txt" \
                     --subsample=$Nmin \
                     --win-size=$win_size \
                     "$bank2_file" \
                     "out_in_${name_bank2}.group2.${name_bank2}"

chromosight quantify "group12_detected_in_${name_bank1}_${name_bank2}.txt" \
                     --subsample=$Nmin \
                     --win-size=$win_size \
                     "$bank2_file" \
                     "out_in_${name_bank2}.group12"

# merging of pdf plots: 
# Merge pileups of group1 into a single pdf
pdfjam "out_in_${name_bank1}.group1.${name_bank1}.pdf" "out_in_${name_bank2}.group1.${name_bank1}.pdf"  --nup 2x1 --landscape --outfile all_pileups_g1.pdf
# Same for group2
pdfjam "out_in_${name_bank1}.group2.${name_bank2}.pdf" "out_in_${name_bank2}.group2.${name_bank2}.pdf"  --nup 2x1 --landscape --outfile all_pileups_g2.pdf
# Same for intersection
pdfjam "out_in_${name_bank1}.group12.pdf"  "out_in_${name_bank2}.group12.pdf" --nup 2x1 --landscape --outfile all_pileups_g12.pdf
# Combine all pileups into a single pdf
pdfjam all_pileups_g1.pdf all_pileups_g12.pdf all_pileups_g2.pdf --nup 6x1 --landscape --outfile "all_pileups_${name_bank1}.${name_bank2}.pdf"
pdfjam "all_pileups_${name_bank1}.${name_bank2}.pdf" "venn_diagram_group_loops_${name_bank1}"_"${name_bank2}.pdf" --nup 1x2  --outfile "final_${name_bank1}.${name_bank2}.pdf"

echo "Computations and plots are finished! Just relax, go to Breguet."



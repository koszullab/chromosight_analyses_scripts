Python codes used in the analyses.

* `chromo_simul.py`: Generate simulated contact maps with known loop and border coordinates for a single chromosome
* `common_chip_patterns.py`: visualise the proportion of common 2d positions between a reference file and a set of files.
* `common_patterns_upsetplot.py`: Generate an upsetplot to visualize the proportion of common positions between multiple files containing 2D coordinates.
* `plot_sc_ca_maps.py`: Script used to plot S. cerevisiae and C. albicans contact maps along with the coordinates manually selected on Scer (using generate-config with the `--click` option) and coordinates detected on Calb.
* `score.py`: Given a file with detected coordinates, and another file with real coordinates, this script computes the precision, recall and F1 score based on the confusion matrix.
* `zoom_compare_calls.py`: Script used to vizualise loop calls from different softwares in the same region. 

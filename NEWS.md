# UMIcountR 0.5.0
* new functions to process complex molecular spike-in set
* add complex spike-in set barcode annotation (`molspike_barcodes_infos_fivep_final`)

# UMIcountR 0.2.2
Handle some unexpected situations.

# UMIcountR 0.2.1
Fix implementation of ngram by passing split UMI length correctly.

# UMIcountR 0.2.0

* Large change in the error correction utilities. Now builds on BK-Tree datastructure for massive increase in computational efficiency in UMI error corrections. Builds on pybktree package in python called via reticulate.

# UMIcountR 0.1.2

* Add option to subset the length of the spUMI when loading.
* Correct handling of gene name / contig when loading.
* Removed bug with hd1 of spUMI.
* Fix bug if ggrastr not installed.

# UMIcountR 0.1.1

* Add biocViews attribute for correct installation behavior.

# UMIcountR 0.1.0

* Initial version of the UMIcountR package with `extract_spike_dat`, `get_overrepresented_spikes`, `plot_spike_distances`, `return_corrected_umi` and `subsample_recompute` functions.

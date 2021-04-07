
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UMIcountR

<!-- badges: start -->

<!-- badges: end -->

## Installation

\#You can install the released version of UMIcountR from
[CRAN](https://CRAN.R-project.org) with:

\#`r #install.packages("UMIcountR") #`

You can install UMIcountR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cziegenhain/UMIcountR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(UMIcountR)
## basic example code
#load reads from the provided example bam file (Smart-seq3 data)
bam_path <- system.file("extdata", "Smartseq3.TTACCTGCCAGATTCG.bam", package = "UMIcountR", mustWork = TRUE)
spikedat <- extract_spike_dat(bam_path, match_seq_before_UMI = "GAGCCTGGGGGAACAGGTAGG", match_seq_after_UMI = "CTCGGAGGAGAAA")
#> [1] "Reading in data from bam file..."
#> [1] "Hamming correct spikeUMIs..."
```

After loading the data, we can see the data structure:

``` r
str(spikedat)
#> Classes 'data.table' and 'data.frame':   47727 obs. of  13 variables:
#>  $ contig      : Factor w/ 195 levels "1","10","11",..: 195 195 195 195 195 195 195 195 195 195 ...
#>  $ pos         : int  5641 5641 5641 5641 5641 5641 5641 5641 5641 5641 ...
#>  $ CIGAR       : chr  "53M" "53M" "53M" "53M" ...
#>  $ seq         : chr  "GAGCCTGGGGGAACAGGTAGGTAGTGTTGACTACTCGAGCTCGGAGGAGAAAA" "GAGCCTGGGGGAACAGGTAGGACTTGCGCGGTGAGCAAGCTCGGAGGAGAAAA" "GAGCCTGGGGGAACAGGTAGGTTCCAAAAGCAACTCGAGCTCGGAGGAGAAAA" "GAGCCTGGGGGAACAGGTAGGCTTCGTATATTCATTGAGCTCGGAGGAGAAAA" ...
#>  $ BC          : chr  "TTACCTGCCAGATTCG" "TTACCTGCCAGATTCG" "TTACCTGCCAGATTCG" "TTACCTGCCAGATTCG" ...
#>  $ QU          : chr  "EEEEEEEE" "EEEEEEEE" "EEEEEEEE" "EEEEEEEE" ...
#>  $ UX          : chr  "ACTGAGTG" "AGTGGACA" "AAAGGCCC" "AATCATAA" ...
#>  $ UB          : chr  "ACTGAGTG" "AGCGGACA" "AAAGTCCC" "AATCATGA" ...
#>  $ TSSseq      : chr  "GAGCCTGGGGGAACAGGTAGG" "GAGCCTGGGGGAACAGGTAGG" "GAGCCTGGGGGAACAGGTAGG" "GAGCCTGGGGGAACAGGTAGG" ...
#>  $ spikeUMI    : chr  "TAGTGTTGACTACTCGAG" "ACTTGCGCGGTGAGCAAG" "TTCCAAAAGCAACTCGAG" "CTTCGTATATTCATTGAG" ...
#>  $ seqAfterUMI : chr  "CTCGGAGGAGAAAA" "CTCGGAGGAGAAAA" "CTCGGAGGAGAAAA" "CTCGGAGGAGAAAA" ...
#>  $ spikeUMI_hd1: chr  "TAGTGTTGACTACTCGAG" "ACTTGCGCGGTGAGCAAG" "TTCCAAAAGCAACTCGAG" "CTTCGTATATTCATTGAG" ...
#>  $ spikeUMI_hd2: chr  "TAGTGTTGACTACTCGAG" "ACTTGCGCGGTGAGCAAG" "TTCCAAAAGCAACTCGAG" "CTTCGTATATTCATTGAG" ...
#>  - attr(*, ".internal.selfref")=<externalptr>
```

Next, we can run the filtering for overrespresented spUMIs:
<img src="man/figures/README-pressure-1.png" width="100%" />

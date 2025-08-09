# BHiCect2

<!-- badges: start -->

<!-- badges: end -->

The goal of BHiCect2 is to cluster HiC data as described in our [manuscript](). <!-- markdownlint-disable-line MD042 -->
Briefly, we decompose intra-chromosomal HiC data into nested clusters of chromosome regions across multiple resolutions
starting from the complete chromosome all the way to DNA-loops at the maximum resolution provided.

## Installation

You can install the current version of BHiCect2 from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("princeps091-binf/BHiCect2")
```

## Example

BHiCect2 offers a core function to cluster the input HiC data.
The expected input data to BHiCect2 is a list of dataframes containing the HiC data in a three columns format for the various resolution provided by the user.

```r
library(BHiCect2)
## basic example using chromosome 22 data from human GM12878 data (Rao et al. 2014)
data(chr_dat_l)

# manual entry for resolutions
res_set<-c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num<-c(1e6L,5e5L,1e5L,5e4L,1e4L,5e3L)
names(res_num)<-res_set

BHiCect_results<-BHiCect(res_set,res_num,chr_dat_l,4)

```

An overview that describes the main components of the package.
For more complex packages, this will point to vignettes for more details.
This is also a good place to describe how your package fits into the ecosystem of its target domain.

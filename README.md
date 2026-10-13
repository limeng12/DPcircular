# DPcircular: Circular Sequence Detection via DP

## Overview

`DPcircular` detects circular (repeat-unit) structures in DNA/RNA sequences using a modified Smith-Waterman local alignment algorithm. It compresses homopolymer runs, recursively splits confirmed circles, and builds consensus sequences for linear units.

### Key Features

-   dynamic programming
-   **Recursive circle detection** with multi-unit splitting
-   **Reverse-complement filtering** to remove false positives
-   Homopolymer compression (runs >= 6 bp)
-   **Consensus sequence** generation via MUSCLE when > 2 units are found
-   Dot plot visualization of self-alignments
-   High-performance **Rcpp/C++** implementation, O(n^2) time complexity.

## Installation

### Prerequisites

-   R (>= 4.0)
-   A C++11 compatible compiler (GCC, Clang, or Rtools for Windows)

### Install Dependencies

``` r
install.packages(c("Rcpp", "stringr", "stringi", "ggplot2",
                   "scattermore", "ggpointdensity"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "msa"))
```

### Install the Package

``` r

library(devtools)
install_github("limeng12/DPcircular")
```

## Quick Start

``` r

library(Biostrings)
library(DPcircular)

# Read sequences
readfile <-  system.file("extdata", "example_seq.txt", package = "DPcircular");

dna_seq <- readDNAStringSet(readfile)
dna_seq_char <- as.character(dna_seq)

# Detect circular structures
res <- get_circular_count(dna_seq_char[2])

res$is_circuler    # TRUE / FALSE
res$circle_num     # number of circular units
res$units          # linear units
res$consensus      # consensus sequence (if units > 2)

# Visualize self-alignment
dotplot_optimized(dna_seq_char[2], dna_seq_char[2], "xx")
```

## Key Functions

### `get_circular_count()`

Recursive detection of circular repeat units using modified Smith-Waterman self-alignment.

``` r
res <- get_circular_count(seq)
```

**Returns**

| Field | Description |
|-------|-------------|
| `is_circuler` | Logical; TRUE if a circle was detected |
| `circle_num` | Number of detected circular units |
| `units` | Character vector of linear units |
| `consensus` | Consensus sequence (units > 2), else NULL |
| `score` | Circular score (normal_score / path_length) |

### `dotplot_optimized()`

Dot plot visualization of self-alignment via sliding-window match rate.

``` r
dotplot_optimized(seq1, seq2, read_name,
                  window = 10, step = 1,
                  min_match = 0.9,
                  max_points = 100000,
                  use_rcpp = TRUE)
```

| Param | Default | Description |
|-------|---------|-------------|
| `seq1`, `seq2` | — | Sequences to compare |
| `read_name` | — | Axis label |
| `window` | `10` | Sliding window size |
| `step` | `1` | Window step |
| `min_match` | `0.9` | Minimum match rate to plot |

## Algorithm

1.  **Homopolymer compression**: runs >= 6 bp collapsed to one character with recorded frequency.
2.  **Compressed self-alignment**: only cells with `j < i - offset` are filled (`offset = n/2 - min(n/4, 2050)`).
3.  **Reverse-complement filter**: rejects candidates whose reverse-complement score is > 2x forward score.
4.  **Recursive splitting**: confirmed circles are cut at `start_j` and pushed back onto the stack.
5.  **Consensus**: when `units > 2`, MUSCLE alignment + `consensusString()` produces an IUPAC consensus.

## Troubleshooting

1.  **`scattermore` not found**: `install.packages("scattermore")`; otherwise falls back to pure R.
2.  **Compilation errors**: Verify a C++11 compiler is configured with R (`Rcpp::sourceCpp`).
3.  **No circles detected**: Check `path_len / n` and `start_j` thresholds.

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.


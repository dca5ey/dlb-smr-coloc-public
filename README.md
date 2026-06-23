# Cell-type-specific colocalization analysis: DLB GWAS × brain eQTLs

This repository contains the R code used to perform the Bayesian
colocalization analysis reported in Table 1 of:

> [Author]. Cell-type-specific integration of dementia with Lewy bodies
> GWAS and aging brain eQTLs. [Journal, year, DOI — add once published]

## What this code does

Tests for colocalization between genetic risk for dementia with Lewy
bodies (DLB) and cis-eQTLs across five candidate genes (APOE, BIN1,
GBA, SNCA, TMEM175) and six brain cell types, using the `coloc` R
package (Giambartolomei et al., 2014).

## Requirements

- R (version 4.6.0 or later recommended)
- The `renv` package, to restore the exact package versions used

## Setup

1. Clone this repository
2. In R, run:
```r
   install.packages("renv")
   renv::restore()
```
   This installs the exact package versions (including `coloc` v5.2.3)
   used in the original analysis.

## Data Setup (Required Before Running)

This repository does not include the underlying genetic data, in
accordance with the original data sources' access terms. You will
need to download and place the following files yourself:

**GWAS summary statistics**
- Source: GWAS Catalog, accession [GCST90001390](https://www.ebi.ac.uk/gwas/)
- Original publication: Chia et al., 2021, Nature Genetics
- Reformat to tab-delimited columns: `SNP  A1  A2  freq  b  se  p  n`
- Save as: `data/LBD_GWAS_SMR_formatted.txt`

**Cell-type-specific eQTL summary statistics**
- Source: AD Knowledge Portal, Synapse ID [syn52335807](https://adknowledgeportal.synapse.org/)
- Original publication: Fujita et al., 2024, Nature Genetics
- Requires a free Synapse account; some files may require a Data Use
  Certificate — see the AD Knowledge Portal for current access terms
- Place per-gene, per-cell-type `.esd` files (tab-delimited:
  `Chr SNP Bp A1 A2 Freq Beta se p`) in: `data/eqtl_esd/`
- Expected naming convention: `{CellType}_{Gene}.esd`
  (e.g., `Mic_SNCA.esd`)

## Running the Analysis

```r
source("run_analysis.R")
```

This will loop through all available gene × cell-type combinations,
run colocalization analysis on each, and save results to
`results/coloc_results.csv`. Expect 28 rows in the final output if all
data is correctly in place.

## A Note on SMR

The original study also performed summary-data-based Mendelian
randomization (SMR) and the HEIDI test, reported alongside the
colocalization results in Table 1. SMR is a standalone command-line
program (not an R package) and is not included in this repository.

- SMR software: available from the [Yang Lab](https://yanglab.westlake.edu.cn/software/smr/)
  (versions 1.03 and 1.4.1 were used; see Methods for details on a
  known ARM64 v1.4.1 compatibility issue with older-format BESD files)
- LD reference panel: 1000 Genomes Project Phase 3, European subset,
  pre-filtered via the `bigsnpr` R package (Privé et al., 2018)

SMR output is reported as-is in Table 1; this repository reproduces
the colocalization half of the analysis only.

## Output

Running the analysis reproduces the PP.H4 (posterior probability of a
shared causal variant) and PP.H3 (posterior probability of two
distinct causal variants) values reported in Table 1, including the
three combinations meeting the PP.H4 > 0.75 threshold for evidence of
colocalization: SNCA in microglia and astrocytes, and BIN1 in
microglia.

## License

MIT License - see LICENSE file for details.

Copyright (c) 2026 @dca5ey (github username)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ON ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Contact

Questions about this code can be directed to @dca5ey (github username)
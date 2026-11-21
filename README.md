# Haplotype Characterization

Run SNP quality checks, haplotype splitting, line assignments, and the Excel report from one R script.

## Run it in RStudio

Install the workbook package once.

```r
install.packages("openxlsx")
```

Open the repository folder in RStudio and run the following, using the path to your original SNP CSV.

```r
source("run_haplotype.R")
run_haplotype("snp_haplotype2_525.csv", "haplotype_results")
```

Run the same workflow from a terminal.

```sh
Rscript run_haplotype.R snp_haplotype2_525.csv haplotype_results
```

With omitted output folder, the script creates `haplotype_results` next to the input. Choose a new folder for a second run or pass `overwrite = TRUE` in R to replace results.

## Input and results

The CSV needs a unique `Line` or `germplasmName` for each sample and marker columns beginning with `SNP_chr`. Use A, C, G, and T for homozygous calls. The workflow recognizes slash or asterisk calls and `1` as unresolved or heterozygous, and treats `FAIL`, `99`, and blanks as missing. Original marker values stay in the line report. Include `general_grouping` if you want a cross-tabulation against your existing groups.

- **output.xlsx** contains SNP quality and frequency summaries, split markers, line assignments, representative haplotypes, and group allele similarity.
- **final_merged_with_specific2.csv** contains the grouping path and original marker calls for each line. A `0` denotes a depth where the line's group was never split. A `U` in the specific path denotes an unresolved or missing call at a split marker. Those lines follow the major branch, following the original convention.
- **filtered_commonsnps.csv** contains one representative row per group, with tied calls reported together.
- **grm.csv** contains the original haplotype workflow's group allele-similarity score.

## Change the split settings

The default is up to three levels of splitting, using the highest within-group minor allele frequency at each split. Markers used for splitting may have at most 10% missing and 10% unresolved or heterozygous calls. The workbook's filtered view uses 20% for the latter limit.

Choose a different rule or set thresholds in the same call.

```r
run_haplotype(
  "snp_haplotype2_525.csv", "haplotype_results_threshold",
  split_rule = "min_maf_fallback", min_maf = 0.35, fallback_rank = 5L
)
```

`min_maf` applies the MAF threshold and ranks eligible markers by the smallest count of known homozygous calls, corresponding to the original threshold option. `min_maf_fallback` also considers the top MAF markers when none meets the threshold. Splits use markers with exactly two observed homozygous alleles in the group. Ties use input marker order and alphabetical allele order.

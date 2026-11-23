# Haplotype Characterization

Run SNP quality checks, haplotype splitting, line assignments, and the Excel report from one R script.

## Use the Shiny app

Install the app packages once.

```r
install.packages(c("shiny", "openxlsx"))
```

Open the repository folder in RStudio and click **Run App** in `app.R`, or run `shiny::runApp(".")`.

Upload a **CSV** or **XLSX** containing the original SNP calls. The app runs the workflow as soon as the file uploads. If several Excel sheets have line IDs and SNP columns, the app uses the first and lets you switch sheets.

Use the left controls to set SNP filtering limits, maximum split depth, and the SNP selection rule. Percent fields use whole percentages such as 35 for a MAF of 0.35.

- **Preview SNP Filters** opens the SNP Filtering tab, with marker counts, exclusion reasons, and a QC CSV.
- **Run Haplotype Analysis** applies your edits and opens the Haplotype Results tab, with the groups and report downloads.

Download the finished Excel report or a ZIP containing the line assignments, representative haplotypes, and allele similarity CSVs. Store line IDs as text in Excel to preserve leading zeros.

### Review uncertain lines

Open **Group Review** to review lines with a `U` at a split SNP. The app compares their SNP calls with signatures built from lines in the current upload with resolved calls at their splits. It uses QC-passed markers and looks for consistent alleles that differ between groups. Matching homozygous calls give full support. A call such as `A/G` gives partial support to either matching allele and cannot separate two groups that expect `A` and `G`. Missing calls and calls such as `1` and `*/*` give no allele-level evidence.

The **Group Review Suggestions** worksheet in the Excel report shows the original group, candidates, supporting SNP counts, and the reason for a suggestion or a **Needs Review** result. The signature scores describe agreement with this upload's reference lines. Group labels are specific to the split tree from the current run. In Shiny, inspect a line's colored calls and group signatures, then save and download review decisions before running a new analysis.

At each SNP, a group signature requires calls from at least two resolved reference lines and at least half the references in that group, with at least 80% agreement on the allele. A suggestion requires calls at a quarter or more of the SNPs that distinguish the candidates, at least two SNPs favoring the best group, at least 60% weighted agreement, and a lead of at least 15 percentage points in the head-to-head comparisons. A line also stays under review if its original group lacks a usable signature.

Install `uwot` once to use the **UMAP** tab.

```r
install.packages("uwot")
```

Click **Build UMAP** to explore lines with comparable calls at half or more of the QC-passed markers, with a minimum of two markers. Colors show original groups, and outlines identify lines flagged for review. Click a point to inspect its SNP letters. Drag a box and click **Review Selected Lines** to see the flagged lines in that area together. Use visible clusters as candidates for review and check the underlying SNP evidence before assigning a group.

## Run the script directly

From the repository folder, you can also use the same workflow in RStudio with your SNP file path.

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

The CSV or Excel worksheet needs a unique `germplasmName` for each sample and marker columns beginning with `SNP_chr`. Older files can use `Line`. When both columns are present, the workflow uses `germplasmName`. The workbook and CSV exports use `germplasmName` for sample IDs. Use A, C, G, and T for homozygous calls. The workflow recognizes slash or asterisk calls and `1` as unresolved or heterozygous, and treats `FAIL`, `99`, and blanks as missing. Original marker values stay in the line report. Include `general_grouping` if you want a cross-tabulation against your existing groups.

If the input starts with `germplasmName`, the workflow adds any missing `specific_grouping` or `general_grouping` columns and fills them with `NoName`. When it adds a column, it places both grouping columns before `germplasmName` in the sample reports. The cross-tabulation uses `general_grouping` supplied in the input.

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

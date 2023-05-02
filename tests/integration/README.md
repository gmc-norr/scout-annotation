# Workflow test data

## VCF files

- `data/HD832_chr7_twist-solid-0.1.5-alpha.vcf`: a subset of an ensembled VCF file from an HD832 sample that has been run through the [Twist Solid pipeline](https://github.com/Genomic-Medicine-Sweden/Twist_Solid) (v0.1.5-alpha).

## Tests

There are a number of cases that get tested.

- Panel filtering but no SNV filtering with both empty and non-empty results
- No panel filtering but SNV filtering with both empty and non-empty results
- Both panel filtering and SNV filtering with both empty and non-empty results

sample  | panel                   | snv_filtering | track        | expected variants
--------|-------------------------|---------------|--------------|-------------------
sample1 | test_panel1             |               | cancer       | >0
sample2 | test_panel3             |               | rare_disease | 0
sample3 |                         | strict        | cancer       | 0
sample4 |                         | rare_disease  | rare_disease | >0
sample5 | test_panel1,test_panel2 | cancer        | cancer       | >0
sample6 | test_panel3             | rare_disease  | rare_disease | 0
sample7 | test_panel2             | strict        | rare_disease | 0

### Panels

- test_panel1: EGFR, BRAF
- test_panel2: MET, BRAF
- test_panel3: NEG (Non-Existing Gene)

### Filtering

- cancer: somatic variant filtering
- rare_disease: constitutional variant filtering
- strict: filtering so strict that nothing makes it through

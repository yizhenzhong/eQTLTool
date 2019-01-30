# Hierarchical correction method to call significant eQTLs

This method is described in [Huang et al](https://academic.oup.com/nar/article/46/22/e133/5090771).

There are three steps in the correction method as descrbed in the paper:
- p values of all cis SNPs are adjusted for multiple testing for each gene using Benjamini and Yekutieli (BY) method as the locally adjusted pvalues.
- The minimum BY-adjusted pvalue for each gene are corrected using Benjamini and Hochberg (BH) method as the globally adjusted p values (BY-BH pvalues).
- For a chosen threshold (for example, 0.05), find the largest BY-BH pvalues under the threshold and the corresponding BY-adjusted p value. This BY-adjusted p value will be used as the threshold to call signfiant eQTLs for each gene.

Usage

- 1st argument is the matrix eQTL output
- 2nd argument is the threshold (default is 0.05)

```bash
python get_BY_BH_correction.py \
matrixeqtl.result.out \
0.01 
```


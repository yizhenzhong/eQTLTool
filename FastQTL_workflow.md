# To run eQTL mapping with [FastQTL](http://fastqtl.sourceforge.net/)

Note: Use fastQTL.static

program should be able to call directly


## Input

All the input files should be put in the working_dir/data

- vcf file (vcf.gz and vcf.gz.tbi) of FastQTL format
- bed file (bed.gz and bed.gz.tbi) of FastQTL format
- covariate file (cov.txt) of FastQTL format

## Step1
Run permutation mode to get the number of eGenes

```bash
python eQTL_pipeline.py vcf.gz bed.gz cov.txt \
  --working_dir ./  \
  --analysis permute \
  --prefix eQTL \  
```


## Step2
Get number of eGene

what it does:
- FDR correction of permutation pvalue for each gene
- get the genes whose FDR-corrected pvalue < 0.05 (eGene) and write to working_dir/result/eQTL.sig (three columns are genename, nominal pvalue threshold for this gene, and permutation-pvalue)


```bash
python eQTL_pipeline.py vcf.gz bed.gz cov.txt \
  --working_dir ./  \
  --analysis eGene \
  --prefix eQTL  
```

## Step3

Nominal pass of eQTL mapping

```bash
python eQTL_pipeline.py vcf.gz bed.gz cov.txt \
  --working_dir ./  \
  --analysis nominal \
  --prefix eQTL   
```

## Step4

Call significant eQTLs for each gene

```bash
python eQTL_pipeline.py vcf.gz bed.gz cov.txt \
  --working_dir ./  \
  --analysis call \
  --prefix eQTL   
```

what it does
- cancatenate all nominal pass results to working_dir/result/eQTL_cis_full.txt.gz
- write nomimal associations for eGenes only to working_dir/result/eQTL_cis_full.sig.txt
- get eQTLs whose nominal pvalue is <= per-gene threshold in step 2 and write to working_dir/result/eQTL.significant_eQTLs.txt






<!---
To obtain a gene-specific significance threshold to call eQTL due to the varying allele frequency and LD structure, we used permutations together with empirical beta-distribution approximation approach implemented in QTLTools to model the null distribution of associations at each gene. Gene expression was permuted to obtain the null associations with cis SNPs while retaining the LD structure of genotype. The best association for each permutation across all associations was extracted to estimate the beta-distribution parameters with maximum likelihood estimation method. The nominal p-value of the best association for each gene is compared against the β-distribution to obtain the Beta distribution-adjusted empirical p-value. Beta distribution-adjusted empirical p-values from FastQTL were used to calculate q-values(30), and a false discovery rate (FDR) threshold of ≤5% was applied to identify genes with at least one genome-wide significant cis-eQTL (“eGenes”). QTLTools permutation mode was used with the setting “--permute 1000 10000”.  
To identify the list of all significant variant-gene pairs associated with eGenes, the β-distribution adjusted p-value of the gene whose q value is closest to the 0.05 FDR threshold was then used to calculate a nominal p-value threshold for each gene based on the beta distribution parameters estimated for each eGene. For each gene, variants with a nominal p-value below the gene-level threshold were considered significant and included in the final list of variant-gene pairs. 

-->



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










# Run Matrix eQTL

[Matrix eQTL](https://cran.r-project.org/web/packages/MatrixEQTL/index.html)
```bash
Rscript eqtl.matrixEL.r dosage_file \
  snp_location_file \
  expression_file \
  gene_location \
  covariate \
  output_file_name \
```

### Input:
####  Dosage file
- one header
- one column for SNP ID
- Tab separated


#### snp_location_file
- one header
- space separated

#### Gene expression file
- space separated
- one header and one column for gene id

#### expression position file
- space separated
- one header


#### Covariate file
- covariate * sample
- first column is the covariate ID
- space separated
- one header of sample ID


Note:

modify cisDist=1e6 to specify your cis region

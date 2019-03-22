### mashr workflow


"In this notebook, we did not apply the inference to all the gene-snp pairs. Rather we focused on the "top" gene-snp pairs as a demonstration. It should be straightforward to configure the Posterior computatoin step to work on all gene-snp pairs instead."


## Select strongest/10 random signals
write common genes list to file: common_gene.txt

[Run pipeline to select strong and random signals from eQTL mapping across multiple conditions](https://github.com/yizhenzhong/eQTLTool/blob/master/get_mashr_input.py)

Results from a subset “strong” tests. These tests were identified by taking the “top eQTL” in each gene based on univariate SNP-gene association tests. (Here, “top eQTL” for a given gene is defined as the SNP with the largest (univariate) Z statistic among all tissues). ##largest in single tissue or all tissue

```python
  lp = data.dump('pval') #pandas data.frame for one gene
  lp = lp[np.all(np.isfinite(lp), axis=1) & np.all(np.isfinite(shat), axis=1)]
  lp = -np.log10(lp)
  rowidx = np.where(data['rownames'] == lp.max(axis=1).idxmax())[0][0] #first find the max for each row and become a series and then find the index of the max in that column
```

## Estimate Data-driven covariance
[If you have only access to Z scores, you can set Bhat to the Z scores, and set Shat to be the matrix with all 1s].
produce covariance matrices based on the top 5 PCs of the strong signals. vi
```{r}
data   = mash_set_data(simdata$Bhat, simdata$Shat)
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))
``
## Apply Extreme Deconvolution
```{r}

U.ed = cov_ed(data, U.pca, subset=strong)
```

## Run mash r 
```{r}
U.c = cov_canonical(data)  
m   = mash(data, c(U.c,U.ed))
print(get_loglik(m),digits = 10)
```

Finally, the gene expression measurements in the GTEx study are correlated due to sample overlap (sometimes multiple measurements were obtained from the same individual). Therefore, we have also estimated a correlation matrix, which is stored in dat$vhat: #this is unclear how to calculate it


data/MatrixEQTLSumStats.Portable.ld2.Z.rds
```{R}
data = mash_set_data(simdata$Bhat, simdata$Shat) #contain beta and standard errors for each associations

```
#### Set up the covariance matrix

- canonical: 

```{r}
U.c = cov_canonical(data)  
print(names(U.c))
```
- data-driven
Compute a sparse factorization of the (centered) z-scores using the SFA software, with K = 5 factors, and save the factors in an .rds file. This will be used to construct the mixture-of-multivariate normals prior. This step is labeled sfa, and should only take a few minutes to run.
#### Fit the model

```{r}
m.c = mash(data, U.c)
```


I tested when fitting using only the z score, the results are simular
```{r}
simdata$Bnew <- matrix(1, nrow=nrow(simdata$Bhat), ncol=ncol(simdata$Bhat)) 
simdata$Zhat <- simdata$Bhat/simdata$Shat
data = mash_set_data(simdata$Zhat, simdata$Bnew)
data_Z = mash_set_data(simdata$Zhat, simdata$Bnew)
m.c.z = mash(data_Z, U.c)
head(get_lfsr(m.c.z))

         condition_1 condition_2 condition_3 condition_4 condition_5
effect_1   0.7561728   0.7706145   0.8053138   0.8111543   0.8189726
effect_2   0.7245512   0.6906096   0.7760496   0.6923021   0.7012937
effect_3   0.7387626   0.7514716   0.8105185   0.8390738   0.8501844
effect_4   0.7855852   0.8454591   0.8360694   0.8660196   0.8503558
effect_5   0.8044666   0.8371182   0.8808539   0.8786020   0.8758416
effect_6   0.7207353   0.6833960   0.7894523   0.7003147   0.7569664

print(length(get_significant_results(m.c.z)))
[1] 150
print(get_loglik(m.c.z))
[1] -16120.32

```

* This step must be preformed using all the tests (or a large random subset), because this is where mash learns that many tests are null and corrects for it.

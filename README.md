# eQTLTool

## This toolkit provides downstream analysis after eQTL mapping.

Including:
* matching null SNPs for eQTLs (--match)
* eQTL enrichment in functional annotations compared with a null set (--enrich)

### Required input/argument:
* eQTL file (--eqtl) : eQTL file including required columns (chr, pos, geneid, snpid), optional columns (beta, pvalue).
* type of analysis (--analysis): choose from [match, enrich].

### Optional input/argument:
* Ouput dir (--output_dir): default is current directory.
* Output prefix (--output_prefix): default is eqtl.
* file including features to match (--match_feature): required columns (snpid, MAF, minDIST, LDSC). The snpid column will be matched with eQTL file snpid.
* file including features to test for enrichment (--enrich_feature): required columns (snpid(1st), enrichment features). The snpid column will be matched with eQTL file snpid.
* index of features to test enrichment (--index). Choose from 1 to the total number of features in enrich_feature file).


### Analysis:

##### match: match eQTL to null set of SNPs by quantile bins
Cut the MAF, LDSC, and minDIST with 10 quantile bins. For each eQTL, find SNPs in the same bin for all three criterias. Sampled with replacement for 10000 times.

- output file: prefix+match.txt nrow is the number of eQTLs, ncol is 10000 (10000 null set) + 1 (eQTLs to be matched)
                prefix+count.txt how many eQTLs matched for each queried SNPs
- example: 

```python
python funcs.py --analysis match --eqtl hepatocytes.txt \
        --output_dir test \
        --output_prefix hepatocytes \
        --match_feature match_feature.txt
```
##### enrich: find the number of feature in true eQTLs and matched null set
- input: --enrich_feature, (matched null set from previous match analysis will be readed automatically, should use the same output_prefix)

- example:
```python
 python funcs.py --analysis enrich --eqtl hepatocytes.txt \
         --output_dir test \
         --output_prefix hepatocytes \
         --enrich_feature encode_features.txt
         --index 1
```



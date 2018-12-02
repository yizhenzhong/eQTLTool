#This toolkit provides downstream analysis after eQTL mapping.

Including:
* matching null SNPs for eQTLs
* eQTL enrichment in functional annotations compared with a null set

###Required input/argument:
* eQTL file (--eqtl) : eQTL files including required columns (chr, pos, geneid, snpid), optional columns (beta, pvalue)
* type of analysis (--analysis): choose from match and enrich

###Optional input/argument:
* Ouput dir (--output_dir): default is current directory
* Output prefix (--output_prefix): default is eqtl
* files including features to match (--match_feature): required columns (snpid (will be matched with eQTL table),
        MAF, minDIST, LDSC



###Analysis:

match: match eQTL to null set of SNPs by MAF within 5%, LDSC within 20% and minDIST within 20%
- output file: nrow is the number of eQTLs, ncol is 10000 + 1 (eQTLs to be matched)



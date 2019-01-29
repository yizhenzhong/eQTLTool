suppressMessages(library(qvalue))

args <- commandArgs(trailingOnly = TRUE)

ifile = args[1]
fdr = as.numeric(args[2]);

cat("Processing fastQTL concatenated output [", ifile, "] controlling for FDR =", fdr * 100, "%\n");

#Read data
D = read.table(ifile, hea=FALSE, stringsAsFactors=FALSE)
D = D[which(!is.na(D[, 10])),] #a second permutation p-value obtained via beta approximation
cat("  * Number of molecular phenotypes =" , nrow(D), "\n")
cat("  * Correlation between Beta approx. and Empirical p-values =", round(cor(D[, 10], D[, 11]), 4), "\n")

#Run qvalue on pvalues for best signals
Q = qvalue(D[, 11])
print(summary(Q))
print(Q$pi0)
D$qval = Q$qvalue
cat("  * Proportion of significant phenotypes =" , round((1 - Q$pi0) * 100, 2), "%\n")

#Determine significance threshold
set0 = D[which(D$qval <= fdr),] 
set1 = D[which(D$qval > fdr),]
pthreshold = (sort(set1$V10)[1] - sort(-1.0 * set0$V10)[1]) / 2 #threshold of fdr 0.05
cat("  * Corrected p-value threshold = ", pthreshold, "\n")

#Calculate nominal pvalue thresholds
D$nthresholds = qbeta(pthreshold, D$V3, D$V4, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#Write output
write.table(D[, c(1, 13, 12)], args[3], quote=FALSE, row.names=FALSE, col.names=FALSE) #write the fdr converted beta-pvalue and beta converted nominal pvalue

write.table(D[D$qval<0.05, c(1, 13, 12)], paste0(args[3],".sig"), quote=FALSE, row.names=FALSE, col.names=FALSE)

cat("Done\n")

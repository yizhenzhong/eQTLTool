#######################################

#eQTL estimation
#######################################
#library(LAMatrix)
library(MatrixEQTL)
#library(methods)
useModel = modelLINEAR

print("file1: dosage, file2: snp location, file3: expression: file4: gene location, file5: covariate, file6: output")
args = commandArgs(trailingOnly=TRUE)
SNP_file_name=args[1]
snps_Location_File_name=args[2]
expression_file=args[3]
gene_location_file_name=args[4]
covariates_file_name=args[5]

output_file_name_tra=tempfile()
output_File_name_cis=args[6]
pvOutputThreshold_cis=1
pvOutputThreshold_tra=0
errorCovariance=numeric()
cisDist=1e6

#########################################
#read covariate 
###########################################

cvrt=SlicedData$new()
table = read.table(covariates_file_name,header = T,check.names = F, stringsAsFactors = FALSE)
cvrt$CreateFromMatrix(as.matrix(table[,-1])) # expected the first column is the ID and from column 2 is the covariates
print(str(cvrt))
#########################################
#read gene expression
###########################################
gene=SlicedData$new();
gene$fileDelimiter=" "
gene$fileOmitCharacters="NA"
gene$fileSkipColumns=1
gene$fileSkipRows=1
gene$fileSliceSize=2000
gene$LoadFile(expression_file)

print(str(gene))
##############################################\
#read snp location and gene position
#################################################


snpspos=read.table(snps_Location_File_name,header=T,stringsAsFactors = FALSE)
genepos=read.table(gene_location_file_name,header =T,stringsAsFactors = FALSE)
head(snpspos)
head(genepos)
############################################
#read snp
############################################

snp=SlicedData$new()
snp$fileDelimiter="\t"
snp$fileOmitCharacters="NA"
snp$fileSkipColumns=1
snp$fileSkipRows=1
snp$fileSliceSize=2000
snp$LoadFile(SNP_file_name)

print(str(snp))
############################################
#read covariate
############################################

#cvrt=SlicedData$new()
#table = read.table(covaiates_file_name,header = T,check.names = F)
#cvrt$CreateFromMatrix(as.matrix(table[,-1]))

#cvrt$fileDelimiter=" "
#cvrt$fileOmitCharacters="NA"
#cvrt$fileSkipColumns=1
#cvrt$fileSkipRows=1
#cvrt$fileSliceSize=2000
#cvrt$LoadFile(covariates_file_name)


#write.table(snpspos,snps_Location_File_name,col.nTImes = FALSE,row.names = FALSE,quote = FALSE)
#############################################
#compute eQTL
#############################################
me=Matrix_eQTL_main(snps=snp,
                    gene=gene,
                    cvrt=cvrt,
                    output_file_name = output_file_name_tra,
                    pvOutputThreshold = pvOutputThreshold_tra,
                    useModel = useModel,
                    errorCovariance = errorCovariance,
                    verbose = TRUE,
                    output_file_name.cis = output_File_name_cis,
                    pvOutputThreshold.cis = pvOutputThreshold_cis,
                    snpspos = snpspos,
                    genepos = genepos,
                    cisDist = cisDist,
                    #noFDRsaveMemory = TRUE,
                    pvalue.hist = "qqplot",
                    min.pv.by.genesnp = FALSE,
                    noFDRsaveMemory = FALSE)


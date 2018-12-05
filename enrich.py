

from random import choices
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict


def match(featureTable, OUT_NAME, eqtl_object):
        outfile = OUT_NAME + "_match.txt"
        outcount = OUT_NAME + "_count.txt"
        
        #get 10 quantile bins, and update the data with bin id
        featureTable.update(pd.qcut(featureTable['MAF'], q=10, labels=range(1,11)))
        featureTable.update(pd.qcut(featureTable['minDIST'], q=10, labels=range(1,11)))
        featureTable.update(pd.qcut(featureTable['LDSC'], q=10, labels=range(1,11)))

        #split the data into eQTL and non eQTL table
        QTLfeature = featureTable[featureTable['snpid'].isin(set(
                 eqtl_object.snpid))].reset_index(drop=True)
        NonQTLfeature = featureTable[-featureTable['snpid'].isin(set(
                eqtl_object.snpid))].reset_index(drop=True)

        print("QTLtable shape:". QTLfeature.shape) 

        #match null SNPs for each eQTL 
        data = []
        count = []
        matched_snp = []
        for index, row in QTLfeature.iterrows():
                temp_snps = NonQTLfeature.loc[(NonQTLfeature['LDSC'] == row.LDSC) & (NonQTLfeature['MAF'] == row.MAF) & (NonQTLfeature['minDIST'] == row.minDIST)]['snpid']
                matched_snp.append(row['snpid'])
                count.append(temp_snps.count())
                try:
                        data.append(temp_snps.sample(10, replace=True ).tolist())
                except: 
                        print(row["snpid"]i, "is not matched to any SNPs!")
                        continue	
        df = pd.DataFrame(data, index = matched_snp)
        count_table = pd.DataFrame({'count':count}, index=matched_snp)
        df.to_csv(outfile, sep="\t", index=True, header=False)
        count_table.to_csv(outcount, sep="\t", index=True)
        return df


def featureDict(enrich_feature, feature):
        feature_dict = dict(zip(enrich_feature.snpid, enrich_feature[feature]))
        return feature_dict

def sumFeature(feature_dict, snplist):
        return sum([feature_dict[i] for i in snplist])

def enrich(ENRICH_FEATURE, OUTPUT_DIR, OUT_PREFIX, eqtl_object, index):
        '''
        EnrichFeature: name of the file name for the feature to test enrichment, header is the name of feature, first column is snpid.
        OUT_NAME: prefix of output file
        eqtl_object: eqtl class

        Return: count of eQTL and null set for each feature.

        '''
        
        enrich_feature = pd.read_table(ENRICH_FEATURE, sep="\t") 
        nullset = pd.read_table(OUTPUT_DIR + "/" + OUT_PREFIX + "_match.txt", sep='\t') #read the generated null set, snps, row is eQTL, column is matched SNPs
        features = list(enrich_feature.columns.values) #get features for enrichment
        feature = features[index]
        
        feature_dict = featureDict(enrich_feature, feature)
        result = [0]*10001
        for n, column in enumerate(nullset):
                snplist = nullset[column]
                result[n] = sumFeature(feature_dict, snplist)
        
        outfile = OUTPUT_DIR + "/enrich/" + OUT_PREFIX + feature + "_enrich.txt"
        df = pd.DataFrame({feature:result})
        df.to_csv(outfile, sep="\t", index=False)






























			

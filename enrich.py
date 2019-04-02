

import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict


def match(featureTable, OUT_NAME, eqtl_object, index=None):
        outfile = OUT_NAME + "_match_{}.txt".format(str(index))
        outcount = OUT_NAME + "_count_{}.txt".format(str(index))
        #get 10 quantile bins, and update the data with bin id
        featureTable.update(pd.qcut(featureTable['MAF'], q=10))
        featureTable.update(pd.qcut(featureTable['minDIST'], q=10, duplicates="drop"))
        featureTable.update(pd.qcut(featureTable['LDSC'], q=10))

        #split the data into eQTL and non eQTL table
        QTLfeature = featureTable[featureTable['snpid'].isin(set(
                 eqtl_object.snpid))].reset_index(drop=True)
        
        NonQTLfeature = featureTable[-featureTable['snpid'].isin(set(
                eqtl_object.snpid))].reset_index(drop=True)

        print("QTLtable shape:", QTLfeature.shape) 

        #match null SNPs for each eQTL 

        if index:
                QTLfeature = QTLfeature.iloc[int(index)*500:(int(index)+1)*500, ]
        data = pd.DataFrame(index=QTLfeature["snpid"],columns=range(10000)) 
        count = []
        print(QTLfeature.shape)
        for index, row in QTLfeature.iterrows():
                print(index)
                
                temp_snps = NonQTLfeature.loc[(NonQTLfeature['LDSC'] == row.LDSC) & (NonQTLfeature['MAF'] == row.MAF) & (NonQTLfeature['minDIST'] == row.minDIST)]['snpid']
                count.append(temp_snps.count())
                try:
                        data.loc[row["snpid"]]=temp_snps.sample(10000, replace=True ).tolist()
                except: 
                        print(row["snpid"], "is not matched to any SNPs!")
                        continue	
               
        count_table = pd.DataFrame({'count':count}, index=QTLfeature['snpid'])
        data.to_csv(outfile, sep="\t", index=True, header=False)
        count_table.to_csv(outcount, sep="\t", index=True)
        return data 


def featureDict(enrich_feature, feature):
        feature_dict = dict(zip(enrich_feature.snpid, enrich_feature[feature]))
        return feature_dict

def sumFeature(feature_dict, snplist):
        return sum([feature_dict[i] for i in snplist])

def enrich(ENRICH_FEATURE, OUTPUT_DIR, OUT_PREFIX, eqtl_object, match, index):
        '''
        EnrichFeature: name of the file name for the feature to test enrichment, header is the name of feature, first column is snpid.
        OUT_NAME: prefix of output file
        eqtl_object: eqtl class

        Return: count of eQTL and null set for each feature.

        '''
        
        enrich_feature = pd.read_table(ENRICH_FEATURE, sep="\t") 
        if match:
                nullset=pd.read_table(match, names=range(1001))
        else:                
                nullset = pd.read_table(OUTPUT_DIR + "/" + OUT_PREFIX + "_match.txt", sep='\t') #read the generated null set, snps, row is eQTL, column is matched SNPs
        print(nullset.iloc[:5,:5])
        features = list(enrich_feature.columns.values) #get features for enrichment
        feature = features[int(index)]
        print(feature) 
        feature_dict = featureDict(enrich_feature, feature)
        result = [0]*1001
        for n, column in enumerate(nullset):
                snplist = nullset[column]
                result[n] = sumFeature(feature_dict, snplist)
        
        outfile = OUTPUT_DIR + "/enrich/" + OUT_PREFIX + "_" + feature + "_enrich.txt"
        df = pd.DataFrame({feature:result})
        df.to_csv(outfile, sep="\t", index=False)






























			

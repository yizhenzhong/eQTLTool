

from random import choices
import pandas as pd
from collections import OrderedDict


def match(QTLfeature, NonQTLfeature, outfile):
        data = [[str(i) for i in range(10001)]]
        for i in range(QTLfeature.shape[0]):
                temp_maf = QTLfeature.at[i,"MAF"]
                temp_ldsc = QTLfeature.at[i,"LDSC"]
                temp_dist = QTLfeature.at[i,"minDIST"]

                try:
                        index1=  NonQTLfeature['MAF'].between(temp_maf*0.95, temp_maf*1.05) & NonQTLfeature['minDIST'].between(temp_dist*0.8, temp_dist*1.2)
                        if temp_ldsc > 0:
                                index2 = NonQTLfeature['LDSC'].between(temp_ldsc*(0.8), temp_ldsc*1.2) 
                        else:
                                index2 = NonQTLfeature['LDSC'].between(temp_ldsc*(1.2), temp_ldsc*0.8)

                        snpin =  NonQTLfeature[index1&index2]['snpid']
                        print(len(snpin))
                        data.append([QTLfeature.at[i,"snpid"]]+ choices(snpin.tolist(), k=10000))
                except: 
                        print(i, QTLfeature.at[i,"snpid"])
                        continue	
                
        df = pd.DataFrame(data[1:], columns=data[0])
        print(df.shape)
        df.to_csv(outfile, sep="\t", index=False)
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






























			

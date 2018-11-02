

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
			snpin = NonQTLfeature['snpid'][(NonQTLfeature['MAF'].between(temp_maf*0.95, temp_maf*1.05) & (NonQTLfeature['LDSC'].between(temp_ldsc*(0.5), temp_ldsc*-1.5) | NonQTLfeature['LDSC'].between(temp_ldsc*(-1.5), temp_ldsc*0.5)) & NonQTLfeature['minDIST'].between(temp_dist*0.5, temp_dist*1.5))]

			data.append([QTLfeature.at[i,"snpid"]]+ choices(snpin.tolist(), k=10000))
		except: 
			print(i, QTLfeature.at[i,"snpid"])
			continue	
		
	df = pd.DataFrame(data[1:], columns=data[0])
	print(df.shape)
	df.to_csv(outfile, sep="\t", index=False)
	return df




def enrich(ENRICH_FEATURE, OUT_NAME, eqtl_object):
	'''
	EnrichFeature: name of the file name for the feature to test enrichment, header is the name of feature, first column is snpid.
	OUT_NAME: prefix of output file
	eqtl_object: eqtl class

	Return: count of eQTL and null set for each feature.

	'''

	enrich_feature = pd.read_table(ENRICH_FEATURE, sep="\t") 
	nullset = pd.read_table(OUT_NAME + "_match.txt", sep='\t') #read the generated null set, snps, row is eQTL, column is matched SNPs
	features = list(enrich_feature.columns.values) #get features for enrichment
	
	res = OrderedDict()
	for feature in features[1:]:
		count_eQTL = sum(enrich_feature.loc[enrich_feature['snpid'].isin(set(eqtl_object.snpid)), feature]) #count how many eQTLs have this feaure
		res[feature] = [count_eQTL]

			
		count_null = list(map(lambda i: sum(enrich_feature.loc[enrich_feature['snpid'].isin(set(nullset.iloc[:,i])), feature]), range(1, 10)))
			 #found in ths null set, how many snps have this feature
		res[feature] = res[feature] + count_null #append the count to the dict

		dfEnrich = pd.DataFrame.from_dict(res)
		outfile = OUT_NAME + "_enrich.txt"
		dfEnrich.to_csv(outfile, sep="\t", index=False)


























			

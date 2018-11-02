import pandas as pd
import sys
import os
import random

def write_fgwas(fgwas_table, fgwas_anno, out, chrs):
    print("Writing chr " + str(chrs))
    out_fn = out + "fgwas_input_chr" + str(chrs) + ".txt.gz"
    anno = pd.read_table(fgwas_anno + "chr" + str(chrs) + ".annot.wdist.wcoding.gz", compression='gzip', sep=' ')
    anno = anno.rename(columns={'chr':'CHR', 'pos': 'POS'})
    fgwas_input = pd.merge(fgwas_table, anno, on=['CHR', 'POS'])
    fgwas_input.to_csv(out_fn, sep=' ', index=False, compression='gzip')
    print(fgwas_input.shape)
    return fgwas_input



def combine_fgwas(out):
    '''
    input: output dir
    
    '''

    fgwas_chr = [out + "fgwas_input_chr" + str(chrs) + ".txt.gz" for chrs in range(1,23)]
    fgwas_pd = [pd.read_table(f, compression='gzip', sep=' ') for f in fgwas_chr]
    result = pd.concat(fgwas_pd)

    result['F'] = [round(random.uniform(0,1),3) for i in range(result.shape[0])] #random because se is used
    result['N'] = [60 for i in range(result.shape[0])] #random because se is used
    result = result.rename(columns={'snpid': 'SNPID'})

    #keep the SNP with smaller pvalue, this step could be removed if the position is updated
    result = result.sort_values(['pvalue']) #sort by pvalue first
    result = result.drop_duplicates(subset=['CHR','POS'], keep='first') #keep the smaller pvalue

    #sort the table by chr and pos
    result = result.sort_values(['CHR', 'POS'])

    #get the segnumber
    seg = {}
    segnumber = []
    index = 1
    for i in result['geneid']:
    	if i in seg:
    		segnumber.append(seg[i])
    	else:
    		seg[i] = index
    		segnumber.append(seg[i])
    		index = index + 1
    result['SEGNUMBER'] = segnumber
    print(result.shape[0], "SNPs of interest")
    result = result.sort_values(['SEGNUMBER'])

    result["POS"] = range(result.shape[0]) #update the position because it is not used in the analysis
    result.to_csv(out + "fgwas_input.txt.gz", sep=' ', index=False, compression='gzip')



def run(fgwas, input_f, out, filter_anno=None):
    
    '''
    run fgwas with an optional filter file
    '''

    if not os.path.exists(out + "/fgwas_out"):
        os.system("mkdir " + out + "/fgwas_out")   
 
    annos = []
    if filter_anno:
            f_annos = open(filter_anno)
    else:
            f_annos = open("./utility/fgwas_453_anno.txt")
    for line in f_annos:
            annos.append(line.strip())
    total = len(annos)
            

    for i in range(len(annos)):
            
            print("Performing " + str(i+1) + "/" + str(total), " annotations...")                   
            script = fgwas + " -i "+ input_f + " -w " + annos[i] + " -fine -k 100 -o ./fgwas_out/" + annos[i]
            print(script)
            try:
                    os.system(script)
            except:
                    continue
               


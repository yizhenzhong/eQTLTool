import sys
import multiprocessing as mp
import pandas as pd
import pickle
import os
import re


def grep(g, eqtl_dir, condi, output_dir):
        print("grep {} from condition{}".format(g, condi))
        if not os.path.exists("{}{}_{}.txt".format( output_dir, g,condi)):
                script = "grep {} {}condition{}_10peer_covariate_dosage.out >{}{}_{}.txt".format(g, eqtl_dir,condi, output_dir, g,condi)
                os.system(script)


def read(f):
        df = pd.read_csv(f, sep='\t', names=["SNP" ,    "gene",    "beta",    "t-stat",  "p-value", "FDR"],index_col=None)
        return df

        
def merge(f):
        merged = reduce(lambda  left,right: pd.merge(left,right,on="SNP",how='outer'), f)
        return merged

def subset(merged, c):
        headers = [n for n, x in enumerate(list(merged.columns.values)) if re.match(c,x)]
        eqtl = merged["SNP"].map(str) + "_" + merged["gene"]
        merged = merged.iloc[:, headers]
        merged.insert(0, "eQTL", eqtl)
        merged.columns=["eQTL"]+["condition"+str(i) for i in range(1, 8)]

        return merged


def get_random(merged):
        return merged.sample(10)


def get_max(merged):
        return merged.loc[merged.iloc[:,range(1,8)].max(axis=1).idxmax()].to_frame().transpose()


def main():
        # Define an output queue
        f="/projects/b1047/zhong/Hepatocyte_project/data/expression/common_genes.txt"
        gene_list = [s.strip() for s in open(f)] #get all the genes
        res = {}
        res["random"] = {}
        res["max"] = {}

        eqtl_dir = "/projects/b1047/zhong/Hepatocyte_project/results/eQTL/"
        output_dir = "/projects/b1047/zhong/Hepatocyte_project/results/mashr/input/"
        index = int(sys.argv[1])
        for gene in gene_list[index*100: (index+1)*100]:#can do parallel
                print(gene)
                for condi in range(1, 8):
                        grep(gene, eqtl_dir, condi, output_dir)
                dfs = [read("{}{}_{}.txt".format(output_dir, gene,condi)) for condi in range(1,8)]
                merged = merge(dfs)
                if merged.size >0:
                        merged_t = subset(merged, "t-stat")
                        res["random"][gene] = get_random(merged_t)
                        res["max"][gene] = get_max(merged_t)


        pd.concat(res["random"].values()).to_csv("{}randon_{}.txt".format( output_dir, index), sep="\t", index=False)
        pd.concat(res["max"].values()).to_csv("{}max_{}.txt".format( output_dir, index), sep="\t", index=False   )


        
if __name__ == '__main__':
        main()

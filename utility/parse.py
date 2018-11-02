import glob
import pandas as pd


def parse_fgwas(out):

        '''parse the fgwas output files,
        output a table ranked by the likelihood

        input: output directory
        '''
        
        out_fn = open(out + "fgwas_individual_llk.txt", "w")
        out_fn.write("\t".join(["annotation", "CI_lo", "estimate","CI_hi", "llk"])+ "\n")
        files = glob.glob(out + "*params")
        for fn in files:
                f = open(fn)
                next(f)
                for line in f:
                        items = line.split()
                fn_llk = fn.replace("params", "llk")
         
                f = open(fn_llk)
                for line in f:
                        if "ln(lk)" in line:  
                                items.append(line.split()[1])
                out_fn.write("\t".join(items)+ "\n")
                
        out_fn.close()
        f = pd.read_table(out + "fgwas_individual_llk.txt")
        f = f.sort_values("llk", ascending=False)
        f.to_csv(out + "fgwas_individual_llk.txt", sep = "\t", index=False)


parse_fgwas('./fgwas_out_common/')

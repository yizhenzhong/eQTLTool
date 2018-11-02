import glob
import pandas as pd


def parse_fgwas(out):

        '''parse the fgwas output files,
        output a table ranked by the likelihood

        input: output directory
        '''
        
        out_fn = open(out + "fgwas_individual_llk.txt", "w")
        out_fn.write("\t".join(["annotation", "llk"])+ "\n")
        files = glob.glob(out + "fgwas_out/*llk")
        for fn in files:
                llk = []
                llk.append(fn.split("/")[-1])
                f = open(fn)
                for line in f:
                        if "ln(lk)" in line:
                                llk.append(line.split()[1])
                out_fn.write("\t".join(llk)+ "\n")
        out_fn.close()
        f = pd.read_table(out + "fgwas_individual_llk.txt")
        f = f.sort_values("llk", ascending=False)
        f.to_csv(out + "fgwas_individual_llk.txt", sep = "\t", index=False)


parse_fgwas('./')

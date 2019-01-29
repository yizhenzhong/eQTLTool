#!/usr/bin/python3
import subprocess
import argparse
import get_eqtl
import os
import subset

def cd(cd_path):
        saved_path = os.getcwd()
        os.chdir(cd_path)
        yield
        os.chdir(saved_path)


def runCis(vcf, bed, cov, dirs, prefix, analysis, samples):
        for i in range(1,30):
                print("cis-eQTL of chunck " + str(i)+ "/30")
                cmd = "cp submit.sh "+ dirs +"/src/" + prefix + "_cis_eQTL_"+ analysis + "_" + str(i)+".sh"
                
                os.system(cmd)
                cmd1 = "cd " + dirs +  "\n\n\n\n"

                
                #cmd2 = "QTLtools_1.1_Ubuntu12.04_x86_64 cis \\\n" \
                cmd2 = "fastQTL.static  \\\n" \
                + "--vcf " + dirs + "/data/" + vcf + " \\\n" \
                + "--bed " + dirs + "/data/" + bed + " \\\n" \
                + "--cov " + dirs + "/data/" + cov + " \\\n" \
                + "--chunk " + str(i) + " 30" + " \\\n" 

                print(cmd2)
                print(samples)
                
                if analysis == "nominal":
                        cmd2 = cmd2  \
                        + "--out " + dirs + "/result/" + prefix + "_cis_" + str(i) + "_30.txt" \
                        #+ "--nominal 1 \\\n"

                elif analysis == "permute":
                        cmd2 = cmd2 + "--permute 1000 \\\n" \
                        + "--out " + dirs + "/result/" + prefix + "_permute_" + str(i) + "_30.txt"

                with open(dirs +"/src/" + prefix + "_cis_eQTL_"+ analysis + "_" + str(i)+".sh" , "a") as myfile:
                        myfile.write(cmd1)
                        myfile.write(cmd2)
                myfile.close()
                os.system("chmod +x " + dirs +"/src/" + prefix + "_cis_eQTL_"+ analysis + "_"+str(i)+".sh")
                os.system("msub " + dirs + "/src/" + prefix + "_cis_eQTL_"+ analysis + "_" + str(i)+".sh")

if __name__=='__main__':
        parser = argparse.ArgumentParser(description='eQTL analysis with QTLTools')
        parser.add_argument('vcf', help='vcf file')
        parser.add_argument('bed', help='bed file')
        parser.add_argument('cov', help='covariate file')
        parser.add_argument('--working_dir', help='directory', default='./')
        parser.add_argument('--analysis',default=None, help='analysis')
        parser.add_argument('--prefix', default='eQTL', help='result file prefix')
        parser.add_argument('--sample', default=None, help='incuded samples')
        args = parser.parse_args()


        vcf = args.vcf
        bed = args.bed
        cov = args.cov
        prefix = args.prefix
        dirs = args.working_dir
        samples = args.sample 
       
        if samples and args.analysis == "subset": 
                '''if samples not empty,
                do the subset, put new input files in data with the given prefix
                ''' 
                print("write vcf")
                #subset.subsetVCF(dirs+"/data/"+vcf,  dirs+"/data/"+prefix+".vcf.gz", dirs+"/data/"+samples)
                
                #subset.index(dirs+"/data/"+prefix+".vcf.gz")
                print("write bed")
                #subset.subsetBed(dirs+"/data/"+bed, dirs+"/data/"+prefix+".bed", dirs+"/data/"+samples, 4)
                #subset.index( dirs+"/data/"+prefix+".bed")
                print("write covariate")
                subset.subsetBed(dirs+"/data/"+cov, dirs+"/data/"+prefix+".covariate.txt", dirs+"/data/"+samples, 1)
        if samples:                
                vcf = prefix+".vcf.gz"
                bed = prefix+".bed.gz"
                cov = prefix+".covariate.txt" #update new input files



        if args.analysis == "nominal":
                print("cis-eQTL mapping nomianl mode")
                runCis(vcf, bed, cov, dirs, prefix, "nominal", samples)


        elif args.analysis == "permute":
                print("cis-eQTL mapping permutation mode")
                runCis(vcf, bed, cov, dirs, prefix, "permute", samples)

        elif args.analysis == "eGene":
                cmd = "cat " + dirs + "/result/" + prefix + "_permute_*_30.txt | gzip -c > " +  dirs + "/result/" + prefix + "_permute_full.txt.gz"
                os.system(cmd)
                cmd = "Rscript calulateNominalPvalueThresholds.R " + dirs + "/result/" + prefix + "_permute_full.txt.gz" + " 0.05 " + dirs + "/result/" + prefix
                
                os.system(cmd)

        elif args.analysis == "call":
                print("call significant eQTLs")
                cmd = "cat " + dirs + "/result/" + prefix + "_cis_*_30.txt | gzip -c > "  + dirs + "/result/" + prefix + "_cis_full.txt.gz"

                print(cmd)
                #os.system(cmd)
                egene = get_eqtl.threshold(dirs+'/result/' + prefix + ".sig")
                #print(egene)
                cmd = "awk \'NR==FNR {a[$1];next} ($1 in a) {print $0}\' " \
                + dirs + "/result/"+prefix + ".sig " \
                + "<(zcat " + dirs + "/result/" +prefix + "_cis_full.txt.gz)" \
                + "> " + dirs + "/result/"+prefix + "_cis_full.sig.txt"
                print(cmd)
                #subprocess.call(cmd, shell=True)

                get_eqtl.call(egene, dirs+"/result/"+prefix + "_cis_full.sig.txt", dirs+"/result/"+prefix + ".significant_eQTLs.txt") 


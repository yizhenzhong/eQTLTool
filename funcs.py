#!/usr/bin/env python3

'''
(c) 2018 Yizhen Zhong

'''

import argparse
from collections import defaultdict
import datetime
import gzip
import numpy as np
import os
import sys
import subprocess
import pandas as pd
import random

import fgwas as fgwas
import utility as utility
import enrich as enrich
from collections import OrderedDict




class eqtls(object):

    '''

    Parameters
    -----------
    eqtl input file
    required columns: chr, pos, ref, alt, geneid, snpid, zscore, se


    Attributes
    -----------
    egenes : number of egenes


    Methods
    ----------
    fwgas(snp, gene, pvalue, zscore):
        run fwgas for functional enrichement
    coloc(snp, gene, pvalue, zsocre):
        run coloc for colocalization analysis
    overlapping(snp, gene, other):
        overlap eqtls with another dataset
    CaVEMaN(snp, gene, zscore):
        run CaVEMaN analysis for estimating causal probability
    '''

    def __init__(self, eqtl, OUT_NAME):
        self.eqtl =  eqtl
        self.outname = OUT_NAME
        self.table = pd.read_table(self.eqtl)
        self.chrs = self.table["chr"]
        self.pos = self.table["pos"]
        #self.ref = self.table["ref"]
        #self.alt = self.table["alt"]
        self.geneid = self.table["geneid"]
        self.snpid = self.table["snpid"]

        #self.se = self.table["se"]
        #self.zscore = self.table["zscore"]

        
        try:
            self.beta = self.table["beta"]
        except:
            self.beta = None #will replace with a function to calculate beta
        
        try:
            self.pvalue = self.table["pvalue"]
        except:
            self.pvalue = None #will replace with a function to calculate pvalue
    
    
    def egenes(self):
        egenes = list(set(list(self.geneid)))
        return egenes
    


    def prepare_fgwas(self, out_dir):
        out = out_dir + "fgwas_input.txt"
        self.fgwas = self.table
        self.fgwas['chr'] = 'chr' + self.fgwas['chr'].astype(str)        
        self.fgwas = self.table.rename(columns={'chr':'CHR', 'pos': 'POS', 'ref' : 'REF', 'alt': 'ALT',
                                               'zscore': 'Z', 'snpid': 'SNPID'})
        
        #s = write_fgwas(self.fgwas, "./1000-genomes-master/", out)
        self.fgwas.to_csv(out, sep=' ', index=False)
        return self.fgwas


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', action='store', dest="analysis", default=False, choices={"fgwas", "match", "enrich"}, help="perform eqtl functional analysis", required=True)
    parser.add_argument('--eqtl', action='store', dest='eqtl', default=None, help="eQTL input file", required=True)
    parser.add_argument('--output_dir', action='store', dest='output_dir', default=None, help ="Output directory")
    parser.add_argument('--output_prefix', action='store', dest='output_prefix', default=None, help="Output prefix")
    parser.add_argument('--fgwas_anno', action='store', dest='fgwas_anno', default=None, help="fgwas annotation directory")
    parser.add_argument('--fgwas_filter', action='store', dest='fgwas_filter', default=None, help="fgwas annotation filter file")
    parser.add_argument('--fgwas_exe', action='store', dest='fgwas_exe', default=None, help="fgwas exe path")
    parser.add_argument('--match_feature', action='store', dest='match_feature', default=None, help="features to match")
    parser.add_argument('--enrich_feature', action='store', dest='enrich_feature', default=None, help="features to match")
    
    args = parser.parse_args()

    #Parse arguments
    ANALYSIS = args.analysis 
    EQTL = args.eqtl 
    FGWAS_ANNO = args.fgwas_anno 
    FGWAS_FILTER = args.fgwas_filter 
    FGWAS_EXE = args.fgwas_exe
    MATCH_FEATURE = args.match_feature
    ENRICH_FEATURE = args.enrich_feature


    #Get output file prefix
    if not args.output_dir:
        OUTPUT_DIR = "./"
    else: 
        OUTPUT_DIR = args.output_dir
    if not args.output_prefix:
        OUTPUT_PREFIX = "eqtl"
    else: 
        OUTPUT_PREFIX = args.output_prefix

    OUT_NAME = OUTPUT_DIR + OUTPUT_PREFIX

    #Error check
    if args.analysis not in ['fgwas', 'match', 'enrich']:
        print("Error: Specified analysis is not supported")
    if not os.path.exists(OUTPUT_DIR):
        print("Error: Output directory does not exists")


    eqtl_object = eqtls(EQTL, OUT_NAME)

    if args.analysis == "match":
    	outfile = OUT_NAME + "_match.txt"
    	fatureTable = pd.read_table(FEATURE, header=None, names=["chr","snpid","A1","A2","MAF","count","minDIST","LDSC"], sep='\t')
    	print(fatureTable.head())
    	QTLfeature = fatureTable[fatureTable['snpid'].isin(set(
    		eqtl_object.snpid))].reset_index(drop=True)
    	print(QTLfeature.shape)
    	NonQTLfeature = fatureTable[-fatureTable['snpid'].isin(set(eqtl_object.snpid))].reset_index(drop=True)
    	print(NonQTLfeature.shape)
    	enrich.match(QTLfeature, NonQTLfeature, outfile)

    if args.analysis == "enrich":
    	enrich.enrich(ENRICH_FEATURE, OUT_NAME, eqtl_object)    
	    
	    
	    	

    '''
    if args.analysis == "fgwas":
    	input_f = OUT_NAME + "fgwas_input.txt.gz"
		if not os.path.exists(OUT_NAME + "_fgwas/"):
		    os.system("mkdir " + OUT_NAME + "_fgwas/")
		    os.system("mkdsir " + OUT_NAME + "_fgwas/input/")    
		#fgwas_input = eqtl_object.prepare_fgwas(OUT_NAME + "_fgwas/")
		#for i in range(1,23):
		#    fgwas.write_fgwas(fgwas_input, FGWAS_ANNO, OUT_NAME + "_fgwas/input/", i)
		fgwas.combine_fgwas(OUT_NAME + "_fgwas/input/")
		#fgwas.run(FGWAS_EXE, input_f, filter_anno=FGWAS_FILTER, OUT_NAME + "_fgwas/")
	'''






if __name__ == '__main__':
    main()

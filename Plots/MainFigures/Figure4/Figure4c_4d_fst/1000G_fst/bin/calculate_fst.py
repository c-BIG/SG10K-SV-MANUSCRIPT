#!/usr/bin/env python

import logging
import argparse
from cyvcf2 import VCF, Writer
import os
import sys

import allel
import numpy as np
import pandas as pd
from random import *
#random.seed(42)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_vcf", required=True, default=None,
                        help="Path to input 1000Genome VCF.")
    parser.add_argument("--annot", required=True, default=None,
                        help="Path to input 1000Genome annotations.")
    parser.add_argument("--seed", required=False, default=100,
                        help="Default seed for Numpy permut.")
    parser.add_argument("--out_csv", required=True, default=None,
                        help="Path to output csv.")
    parser.add_argument("--permut_count", required=False, default=2000,
                        help="Number of times to permutate for p-value calculation. Default: 2000.")
    parser.add_argument("--loglevel", required=False, default="INFO",
                        help="Set logging level to INFO (default), WARNING or DEBUG.")
    args = parser.parse_args()
        
    # checks
    if not os.path.exists(args.in_vcf):
        logging.error("VCF not found: %s" % args.in_vcf)
        sys.exit(1)
        
    return args

def set_logging(loglevel):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % loglevel)
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=numeric_level)

def load_data( vcf ):
    logging.info("Loading VCF...")
    VCF = allel.read_vcf(vcf)
    sampleList = VCF["samples"]
    sampleList.sort()
    return VCF, sampleList

def proc_annotation( annot_file , sampleList ):
    raw_annot = pd.read_table(annot_file)
    raw_annot.columns = [c.strip().lower().replace(' ', '_') for c in raw_annot.columns]
    annot = raw_annot[raw_annot.sample_name.isin(sampleList)]
    annot = annot.sort_values('sample_name')
    annot.reset_index(inplace=True)
    return annot
    
def calc_fst( vcf , annotation , permut_count ):
    values = ["EAS","SAS"]
    
    # Get index of 2 population
    index_list = annotation.index[annotation['superpopulation_code'].isin(values)].to_numpy()
    full_index = annotation.index.to_numpy()
    index_list_comp = full_index[list(set(full_index)-set(index_list))]
    
    df = pd.DataFrame(columns=['ID','chrom','start','fst','p_value' ])
        
    # loop through variants
    for num, var_id in enumerate(vcf['variants/ID']):
        
        # Get GenotypeArray of 2 pop
        GTarray = allel.GenotypeArray([vcf['calldata/GT'][num]])
        #GTarray = allel.GenotypeArray(vcf['calldata/GT'])[num]
        ac1 = GTarray.count_alleles(subpop=index_list)
        ac2 = GTarray.count_alleles(subpop=index_list_comp)

        fst = get_hudson_fst(ac1 , ac2 )
        fst_permutated_list = permutate_genotype( GTarray , annotation , permut_count )
        p_value = get_larger_than_fst( fst , fst_permutated_list ) / permut_count

        entry = pd.DataFrame.from_dict({
            "ID": [var_id] ,
            "chrom": [vcf['variants/CHROM'][num]] ,
            "start": [vcf['variants/POS'][num]] ,
            "fst": fst ,
            "p_value": p_value 
            # remove randomised fst_list output, "randomised_genotype_fst_list": [randomised_genotype_fst_list]
        })
        df = pd.concat([df, entry] , ignore_index=True)
    return df

def get_hudson_fst(ac_pop1, ac_pop2):
    num, den = allel.hudson_fst(ac_pop1, ac_pop2)
    fst = num / den
    return fst

def permutate_genotype( GTarray , annotation , permut_count ):
    randomised_fst_list = []
    values = ["EAS","SAS"]
    
    for i in range( permut_count ):
        annotation["randomised_code"] = np.random.permutation(annotation["superpopulation_code"].values)
        
        index_list = annotation.index[annotation['randomised_code'].isin(values)].to_numpy()
        full_index = annotation.index.to_numpy()
        index_list_comp = full_index[list(set(full_index)-set(index_list))]
    
        ac1 = GTarray.count_alleles(subpop=index_list)
        ac2 = GTarray.count_alleles(subpop=index_list_comp)
        fst = get_hudson_fst(ac1 , ac2 )
        randomised_fst_list.append(fst)

    return randomised_fst_list

def get_larger_than_fst( fst , randomised_fst_list ):
    larger_than_count = sum( i > fst for i in randomised_fst_list)
    return larger_than_count
    
def done(args):
    logging.info("DONE: %s" % args.out_csv)

if __name__ == "__main__":
    args = parse_args()
    set_logging(args.loglevel)

    # load vcf
    vcf, sampleList = load_data(args.in_vcf)
    annotation = proc_annotation(args.annot , sampleList) 
    
    if( not (annotation['sample_name'] == sampleList).all() ):
        print("annotation table and sample-list does not match, exiting...")
        sys.exit(200)
    
    # calculate fst from the observed genotypes
    dfs = calc_fst(vcf , annotation , int(args.permut_count) )
    dfs.to_csv(args.out_csv)
 

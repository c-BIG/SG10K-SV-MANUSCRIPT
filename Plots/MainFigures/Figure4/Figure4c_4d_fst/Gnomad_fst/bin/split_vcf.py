#!/usr/bin/env python

import argparse
import os
import logging
import subprocess
from cyvcf2 import VCF, Writer

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_vcf", required=True, default=None,
                        help="Path to input Gnomad VCF.")
    parser.add_argument("--core", required=False, default=1,
                        help="Number of core to run. Default: 2.")
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

def count_vcf_lines( vcf , core ):
    #result  = subprocess.run(['bcftools' , 'view' , '-H ' , vcf , ' | wc -l '], stdout=subprocess.PIPE , shell=True, text=True)
    result  = subprocess.check_output("bcftools view -H %s | wc -l " % vcf, shell=True, text=True).rstrip() 
    logging.info("%s lines in %s " % ( result , vcf )  )
    return result

def load_data( vcf , core ):
    logging.info("Loading VCF...")
    vcf = VCF( vcf , threads=int(core))
    return vcf
   
def get_header( vcf , core):
    header  = subprocess.check_output("bcftools view -H %s | wc -l " % vcf, shell=True, text=True).rstrip() 
    logging.info("%s lines in %s " % ( result , vcf )  )
    return vcfHeader
   
if __name__ == "__main__":
    args = parse_args()
    set_logging(args.loglevel)
    
    # bcftools count lines
    vcf_lines_count = count_vcf_lines(args.in_vcf , args.split_number )
    vcf_lines_count = get_header( args.in_vcf , core )
    
    # load vcf
    vcf = load_data(args.in_vcf , args.split_number)
    
    # split vcf into split_number files of equal chunks
    split_vcf(vcf , args.split_number , vcf_lines_count )

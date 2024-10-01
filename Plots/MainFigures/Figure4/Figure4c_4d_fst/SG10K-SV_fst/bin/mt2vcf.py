#!/usr/bin/env python3

import logging
import argparse
import allel
import numpy as np
import pandas as pd
import random
random.seed(42)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_mt", required=True, default=None,
                        help="Path to input hail matrix table.")                 
    parser.add_argument("--loglevel", required=False, default="INFO",
                        help="Set logging level to INFO (default), WARNING or DEBUG.")
    args = parser.parse_args()
    return args

def set_logging(loglevel):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % loglevel)
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=numeric_level)


def load_data(args):
    logging.info("Loading hail mt...")
    mt = hl.read_matrix_table(args.in_mt)

    logging.info("Saving temporary VCF...")
    mt = mt.key_rows_by("locus", "alleles")
    out_vcf = "SG10K-SV.tmp.vcf.bgz"
    #hl.export_vcf(mt, out_vcf , tabix = True)
    logging.info("Done: out_vcf" )

    logging.info("Loading ethnicity info...")
    ethnicity = np.array(mt.metadata.Self_Reported_Ethnicity.collect())
    np.save('ethnicity' , ethnicity , allow_pickle=True )
    ethnicity.tofile("ethnicity.txt" , sep="\n")
    return out_vcf, ethnicity
    
if __name__ == "__main__":
    args = parse_args()
    set_logging(args.loglevel)

    # load vcf
    vcf, ethnicity = load_data(args)

#!/usr/bin/env python3

import logging
import argparse
import allel
import numpy as np
import pandas as pd
import random
random.seed(42)
#import hail as hl
#hl.init(default_reference='GRCh38')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_vcf", required=True, default=None,
                        help="Path to input hail matrix table.")
    parser.add_argument("--ethnicity", required=True, default=None,
                        help="Path to input hail matrix table.")
    parser.add_argument("--out_csv", required=True, default=None,
                        help="Path to output csv.")
    parser.add_argument("--bootstrap", required=False, default=0,
                        help="Number of iterations when shuffling ethnicity labels. Default: 0.")
    parser.add_argument("--chr", required=False, default="all_chr",
                        help="Subset the analysis on a specific chromosome (e.g. chr22). Default: all_chr.")
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
    logging.info("Loading VCF into scikit allel...")
    sel_fields = ["variants/ID", "variants/CHROM", "variants/POS", "variants/END",
                "variants/SVLEN", "variants/SVTYPE", "calldata/GT"]
    vcf = allel.read_vcf(args.in_vcf, fields = sel_fields)

    logging.info("Loading ethnicity info...")
    #ethnicity = np.array(mt.metadata.Self_Reported_Ethnicity.collect())
    ethnic = np.load(args.ethnicity , allow_pickle=True )
    return vcf, ethnic


def calc_fst(vcf, ethnicity, args, iter="observed"):
    logging.info("Calculating fst - iteration: %s..." % iter)

    # obtain ac by ethnicity
    ac_chi, ac_ind, ac_mal = get_ac(vcf, ethnicity)

    # calculate hudson fst per site for each pair of populations
    fst_chi_ind = hudson_fst(ac_chi, ac_ind)
    fst_chi_mal = hudson_fst(ac_chi, ac_mal)
    fst_ind_mal = hudson_fst(ac_ind, ac_mal)

    # format output into a dataframe
    r = format_output(vcf, fst_chi_ind, fst_chi_mal, fst_ind_mal, iter)

    return r


def get_ac(vcf, ethnicity):
    gt_array = allel.GenotypeArray(vcf["calldata/GT"])

    chi = np.flatnonzero(ethnicity == "Chinese")
    ind = np.flatnonzero(ethnicity == "Indian")
    mal = np.flatnonzero(ethnicity == "Malay")

    ac_chi = gt_array.count_alleles(subpop = chi)
    ac_ind = gt_array.count_alleles(subpop = ind)
    ac_mal = gt_array.count_alleles(subpop = mal)

    return ac_chi, ac_ind, ac_mal


def hudson_fst(ac_pop1, ac_pop2):
    num, den = allel.hudson_fst(ac_pop1, ac_pop2)
    fst = num / den
    return fst


def format_output(vcf, fst_chi_ind, fst_chi_mal, fst_ind_mal, iter):
    r = pd.DataFrame({
        "rsid": vcf["variants/ID"],
        "svtype": vcf["variants/SVTYPE"],
        "chrom": vcf["variants/CHROM"],
        "pos": vcf["variants/POS"],
        "end": vcf["variants/END"],
        "svlen": vcf["variants/SVLEN"],
        "fst_chi_ind": fst_chi_ind,
        "fst_chi_mal": fst_chi_mal,
        "fst_ind_mal": fst_ind_mal
    })
    r["iter"] = iter

    # NOTE: subsequent data tyding steps implemented in R include (1) fix negative/missing values
    # and (2) calculate fst_max and keep track of which comparison it belongs to

    return r


def save_output(r, args):
    logging.info("Writting outputs to file...")
    r.to_csv(args.out_csv, index = False)


def done(args):
    logging.info("DONE: %s" % args.out_csv)


if __name__ == "__main__":
    args = parse_args()
    set_logging(args.loglevel)

    # load vcf
    vcf, ethnicity = load_data(args)

    # initialise list to store fst results
    fst_datasets = list()

    # calculate fst from the observed genotypes
    r = calc_fst(vcf, ethnicity, args)
    fst_datasets.append(r)

    # bootstrap
    i = 0
    while i < int(args.bootstrap):
        i += 1
        random_ethnicity = np.random.permutation(ethnicity)
        random_r = calc_fst(vcf, random_ethnicity, args, iter="bootstrap-%d" % i)
        fst_datasets.append(random_r)

    # save output
    result = pd.concat(fst_datasets)
    save_output(result, args)

    # done
    done(args)

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
                        help="Path to input Gnomad VCF.")
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
    vcf = VCF( vcf )

    return vcf

def calc_fst( vcf , permut_count ):
    logging.info("Calculating fst for Gnomad Asian (EAS) vs Grouped Non-Asian (AFR,AMR,EUR,OTH)..." )
    
    df = pd.DataFrame(columns=['gnomad_id','chrom','start','stop','svlen','svtype','SG10K_SV_rsid','REF_EAS_AC','ALT_EAS_AC','REF_NonEAS_AC','ALT_NonEAS_AC','fst','call_num' ,'p_value' ])
    #file = open( output ,"w")
    
    # assumed from description below
    # https://www.ncbi.nlm.nih.gov/sites/dbvarapp/studies/nstd166/
    gnomad_sample_size = 10847
    
    for variant in vcf:
        # Initialise values
        EAS_REF_AC = 0
        EAS_AC = 0
        NonEAS_REF_AC = 0
        NonEAS_AC = 0
        EAS_AF = 0
        call_num = 0
        p_value_genotype = 0
        fst_genotype = 0
        
        if variant.FILTER == 'MULTIALLELIC' :
            # EXAMPLE data line
            # CNV #
            # 10      135503000       gnomAD-SV_v2.1_MCNV_10_602      N       <CN=0>,<CN=1>,<CN=2>,<CN=3>,<CN=4>,<CN=5>,<CN=6>,<CN=7>,<CN=8>,<CN=9>,<CN=10>,<CN=11>,<CN=12>   999     MULTIALLELIC    
            # DBVARID=nsv4042429;END=135524350;SVTYPE=MCNV;SVLEN=21350;ALGORITHMS=depth;EVIDENCE=BAF,RD;PCRPLUS_DEPLETED;PROTEIN_CODING__INTERGENIC;PROTEIN_CODING__NEAREST_TSS=FRG2B;
            # AN=10307;AC=0,88,427,5579,2946,1002,215,33,10,5,0,1,1;AF=0,0.008538,0.041428,0.541283,0.285825,0.097215,0.02086,0.003202,0.00097,0.000485,0,9.7e-05,9.7e-05;
            # AFR_AN=4543;AFR_AC=0,18,136,2640,1237,409,85,10,4,4,0,0,0;AFR_AF=0,0.003962,0.029936,0.581114,0.272287,0.090029,0.01871,0.002201,0.00088,0.00088,0,0,0;
            # AMR_AN=735;AMR_AC=0,2,10,317,242,129,28,6,1,0,0,0,0;AMR_AF=0,0.002721,0.013605,0.431293,0.329252,0.17551,0.038095,0.008163,0.001361,0,0,0,0;
            # EAS_AN=1208;EAS_AC=0,48,124,648,284,83,18,0,2,0,0,1,0;EAS_AF=0,0.039735,0.102649,0.536424,0.235099,0.068709,0.014901,0,0.001656,0,0,0.000828,0;
            # EUR_AN=3726;EUR_AC=0,19,153,1926,1160,372,75,16,3,1,0,0,1;EUR_AF=0,0.005099,0.041063,0.516908,0.311326,0.099839,0.020129,0.004294,0.000805,0.000268,0,0,0.000268;
            # OTH_AN=95;OTH_AC=0,1,4,48,23,9,9,1,0,0,0,0,0;OTH_AF=0,0.010526,0.042105,0.505263,0.242105,0.094737,0.094737,0.010526,0,0,0,0,0;
            # SG10K-SV_rsid=SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,SG10K_SV_DUP_chr10_133765614_61,
            
            # Ways to get fst for CNV 
            # https://onlinelibrary.wiley.com/doi/10.1002/humu.21601
            
            # Assumes CNV AC/AN is in the order CN0,CN1,CN2,CN3 etc
            # Assumes CN2 (array[2]) is Hom_ref
            # All other CN status are alternate allele
            EAS_REF_AC = int(list(variant.INFO.get('EAS_AC'))[2])
            EAS_AC = int(variant.INFO.get('EAS_AN')) - EAS_REF_AC
            NonEAS_REF_AC = int(list(variant.INFO.get('AC'))[2] - EAS_REF_AC)
            NonEAS_AC = int(variant.INFO.get('AN')) - EAS_AC - EAS_REF_AC - NonEAS_AC
            
            EAS_AF = 1 - float(list(variant.INFO.get('EAS_AF'))[2])
            call_num =  int(variant.INFO.get('AN'))
                        
        else:
            # Non-CNV # 
            # 1       54665   gnomAD-SV_v2.1_INS_1_1  N       <INS>   248     PASS    
            # DBVARID=nsv4536437;END=54666;SVTYPE=INS;SVLEN=52;CHR2=1;POS2=54716;END2=54717;ALGORITHMS=manta;EVIDENCE=SR;HIGH_SR_BACKGROUND;PROTEIN_CODING__INTERGENIC;PROTEIN_CODING__NEAREST_TSS=OR4F5;
            # AN=21306;AC=2;AF=9.4e-05;N_BI_GENOS=10653;N_HOMREF=10651;N_HET=2;N_HOMALT=0;FREQ_HOMREF=0.999812;FREQ_HET=0.000187741;FREQ_HOMALT=0;
            # AFR_AN=9380;AFR_AC=1;AFR_AF=0.000107;AFR_N_BI_GENOS=4690;AFR_N_HOMREF=4689;AFR_N_HET=1;AFR_N_HOMALT=0;AFR_FREQ_HOMREF=0.999787;AFR_FREQ_HET=0.00021322;AFR_FREQ_HOMALT=0;
            # AMR_AN=1908;AMR_AC=0;AMR_AF=0;AMR_N_BI_GENOS=954;AMR_N_HOMREF=954;AMR_N_HET=0;AMR_N_HOMALT=0;AMR_FREQ_HOMREF=1;AMR_FREQ_HET=0;AMR_FREQ_HOMALT=0;
            # EAS_AN=2366;EAS_AC=0;EAS_AF=0;EAS_N_BI_GENOS=1183;EAS_N_HOMREF=1183;EAS_N_HET=0;EAS_N_HOMALT=0;EAS_FREQ_HOMREF=1;EAS_FREQ_HET=0;EAS_FREQ_HOMALT=0;
            # EUR_AN=7462;EUR_AC=1;EUR_AF=0.000134;EUR_N_BI_GENOS=3731;EUR_N_HOMREF=3730;EUR_N_HET=1;EUR_N_HOMALT=0;EUR_FREQ_HOMREF=0.999732;EUR_FREQ_HET=0.000268025;EUR_FREQ_HOMALT=0;
            # OTH_AN=190;OTH_AC=0;OTH_AF=0;OTH_N_BI_GENOS=95;OTH_N_HOMREF=95;OTH_N_HET=0;OTH_N_HOMALT=0;OTH_FREQ_HOMREF=1;OTH_FREQ_HET=0;OTH_FREQ_HOMALT=0;
            # POPMAX_AF=0.000134;
            # SG10K-SV_rsid=SG10K_SV_DUP_chr1_54720_54,
            EAS_AC = int(variant.INFO.get('EAS_AC'))
            EAS_REF_AC = int(variant.INFO.get('EAS_AN')) - EAS_AC
            NonEAS_AC = int(variant.INFO.get('AFR_AC')) + int(variant.INFO.get('AMR_AC')) + int(variant.INFO.get('EUR_AC')) + int(variant.INFO.get('OTH_AC'))
            NonEAS_REF_AC = int(variant.INFO.get('AN')) - EAS_AC - EAS_REF_AC - NonEAS_AC
            
            EAS_AF = variant.INFO.get('EAS_AF')
            call_num =  int(variant.INFO.get('N_BI_GENOS'))
            
            # Genotype infor
            EAS_HOMREF = int(variant.INFO.get('EAS_N_HOMREF'))
            EAS_HET = int(variant.INFO.get('EAS_N_HET'))
            EAS_HOMALT = int(variant.INFO.get('EAS_N_HOMALT'))
            NonEAS_HOMREF = int(variant.INFO.get('N_HOMREF')) - EAS_HOMREF
            NonEAS_HET = int(variant.INFO.get('N_HET')) - EAS_HET
            NonEAS_HOMALT = int(variant.INFO.get('N_HOMALT')) - EAS_HOMALT
            
            genotype_EAS = []
            genotype_NonEAS = []
            # PERUMTATE ON GENOTYPE LIST
            # EAS genotype list
            gt_EAS_HOMREF = [[0, 0]] * EAS_HOMREF
            gt_EAS_HET = [[0, 1]] * EAS_HET
            gt_EAS_HOMALT = [[1, 1]] * EAS_HOMALT    
            genotype_EAS.extend(gt_EAS_HOMREF)
            genotype_EAS.extend(gt_EAS_HET)
            genotype_EAS.extend(gt_EAS_HOMALT)
            
            # NonEAS genotype list
            gt_NonEAS_HOMREF = [[0, 0]] * NonEAS_HOMREF
            gt_NonEAS_HET = [[0, 1]] * NonEAS_HET
            gt_NonEAS_HOMALT = [[1, 1]] * NonEAS_HOMALT    
            genotype_NonEAS.extend(gt_NonEAS_HOMREF)
            genotype_NonEAS.extend(gt_NonEAS_HET)
            genotype_NonEAS.extend(gt_NonEAS_HOMALT)
            
            # Calculate fst from genotype list
            g_EAS = allel.GenotypeArray([ genotype_EAS ])
            ac_EAS = g_EAS.count_alleles()
            EAS_call_num = EAS_HOMREF + EAS_HET + EAS_HOMALT
            g_NonEAS = allel.GenotypeArray([ genotype_NonEAS ])
            ac_NonEAS = g_NonEAS.count_alleles()
            num, den =  allel.hudson_fst(ac_EAS , ac_NonEAS )
            fst_genotype = float(num/den)
            #print("ORIGINAL_AC\t"+str(ac_EAS)+"\t"+str(ac_NonEAS))
            
            # print(str(variant.ID)+"\t"+str(fst_genotype)+"\t"+ str(ac_EAS) +"\t"+ str(ac_NonEAS))
            combined_gt_list = genotype_EAS.copy()
            combined_gt_list.extend(genotype_NonEAS)
            randomised_genotype_fst_list = permutate_genotype( combined_gt_list , EAS_call_num , permut_count )
            p_value_genotype = get_larger_than_fst( fst_genotype , randomised_genotype_fst_list ) / permut_count
            # print(str(fst_genotype) + "\t" + str(p_value_genotype) + "\t"+ str(randomised_genotype_fst_list))
            
        # Assumes only biallelic
        EAS_ac_array = allel.AlleleCountsArray([[EAS_REF_AC,EAS_AC]])
        Non_EAS_ac_array = allel.AlleleCountsArray([[NonEAS_REF_AC,NonEAS_AC]])
        num, den =  allel.hudson_fst(EAS_ac_array , Non_EAS_ac_array )
        fst = float(num/den)
        # Get call_rate from N_HOMREF/HET/HOMALT
        
        NonEAS_AF = NonEAS_AC / (NonEAS_AC+NonEAS_REF_AC)
        #randomised_fst_list = permutate_ac( EAS_REF_AC , EAS_AC , NonEAS_REF_AC , NonEAS_AC , permut_count )
        #p_value = get_larger_than_fst( fst , randomised_fst_list ) / permut_count
                
        entry = pd.DataFrame.from_dict({
            "gnomad_id": [variant.ID] ,
            "chrom": [variant.CHROM] ,
            "start": [variant.POS] ,
            "stop": [variant.end] ,
            "svlen": [variant.INFO.get('SVLEN')] ,
            "svtype": [variant.INFO.get('SVTYPE')] ,
            "SG10K_SV_rsid": [variant.INFO.get('SG10K-SV_rsid')] ,
            "REF_EAS_AC": [EAS_REF_AC] ,
            "ALT_EAS_AC": [EAS_AC] ,
            "EAS_AF": [EAS_AF] ,
            "REF_NonEAS_AC": [NonEAS_REF_AC] ,
            "ALT_NonEAS_AC": [NonEAS_AC] ,
            "NonEAS_AF": [NonEAS_AF] ,
            "fst": [fst] ,
            "fst_gt": [fst_genotype] ,
            "call_num": [call_num] ,
            "p_value": [p_value_genotype] 
            # remove randomised fst_list output, "randomised_genotype_fst_list": [randomised_genotype_fst_list]
        })
        df = pd.concat([df, entry] , ignore_index=True)
        
        #variant_info = variant.ID + "\t" + variant.CHROM  + "\t" + str(variant.POS) + "\t" + str(variant.end) + "\t" + str(variant.INFO.get('SVLEN')) + "\t" + variant.INFO.get('SVTYPE') + "\t" + variant.INFO.get('SG10K-SV_rsid') + "\t" + str(EAS_REF_AC) + "\t" + str(EAS_AC) + "\t" + str(NonEAS_REF_AC) + "\t" + str(NonEAS_AC) 
        #out_line = variant_info + "\t" + str(fst) + "\t" + str(call_num)+"\n"
        #file.write(out_line)
    #file.close()
    return df

def permutate_ac( EAS_REF_AC , EAS_AC , NonEAS_REF_AC , NonEAS_AC , permut_count ):
    total_EAS = int(EAS_REF_AC) + int(EAS_AC)
    total_NonEAS = int(NonEAS_REF_AC) + int(NonEAS_AC)
    #print("%s\t%s\n" % (total_EAS,total_NonEAS) )
    randomised_fst_list = []
    for i in range( permut_count ):
        # Set EAS_AC as random number from 1 to total_EAS
        rand_EAS_AC = randint(1, total_EAS )
        # Set EAS_REF_AC as total_EAS - rand_EAS_AC
        rand_EAS_REF_AC = total_EAS - rand_EAS_AC
        
        # Set NonEAS_AC as random number from 1 to total_NonEAS
        rand_NonEAS_AC = randint(1, total_NonEAS )
        # Set NonEAS_REF_AC as total_NonEAS - rand_EAS_AC
        rand_NonEAS_REF_AC = total_NonEAS - rand_NonEAS_AC
        
        # Get fst for randomised data
        rand_EAS_ac_array = allel.AlleleCountsArray([[rand_EAS_REF_AC,rand_EAS_AC]])
        rand_Non_EAS_ac_array = allel.AlleleCountsArray([[rand_NonEAS_REF_AC,rand_NonEAS_AC]])
        num, den =  allel.hudson_fst(rand_EAS_ac_array , rand_Non_EAS_ac_array )
        fst = float(num/den)
        randomised_fst_list.append(fst)
    return randomised_fst_list

def permutate_genotype( combined_gt_list , EAS_call_num , permut_count ):
    randomised_fst_list = []

    for i in range( permut_count ):
        shuffle(combined_gt_list)
        # splice to EAS and NonEAS list
        rand_EAS_g = combined_gt_list[0:EAS_call_num]
        rand_NonEAS_g = combined_gt_list[EAS_call_num:]
        # fst cal
        rand_EAS_g = allel.GenotypeArray([ rand_EAS_g ])
        rand_EAS_ac = rand_EAS_g.count_alleles()
        rand_NonEAS_g = allel.GenotypeArray([ rand_NonEAS_g ])
        rand_NonEAS_ac = rand_NonEAS_g.count_alleles()
        num, den =  allel.hudson_fst(rand_EAS_ac , rand_NonEAS_ac )
        fst_genotype = float(num/den)
        randomised_fst_list.append(fst_genotype)
        #print("rand_fst\t" + str(fst_genotype))
        #print("rand_EAS_ac\t" + str(rand_EAS_ac) + "rand_NonEAS_ac\t" + str(rand_NonEAS_ac))
    return randomised_fst_list

def get_larger_than_fst( fst , randomised_fst_list ):
    larger_than_count = sum( i > fst for i in randomised_fst_list)
    return larger_than_count

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
    vcf = load_data(args.in_vcf)

    # calculate fst from the observed genotypes
    dfs = calc_fst(vcf , int(args.permut_count) )
    dfs.to_csv(args.out_csv)
 

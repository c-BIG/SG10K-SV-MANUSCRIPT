#!/usr/bin/env python

import fileinput
import json
import collections
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--intron_bed", required=True, default=None ,
                        help="Input intronic bed")
    parser.add_argument("--intergenic_bed", required=True, default=None  ,
                        help="Input intergenic bed")
    args = parser.parse_args()
    
    return args


def summarise_line(line , feature , id_dict ):
    srs = pd.Series(line.strip().split("\t"))
    
    # Generate SG10K-SV ID + feature key-pair and check existence in dictionary 
    key = srs[3] + "/" + feature
    svtype = srs[4]
    basepairs = 0
    
    if( svtype == 'INS' ):
        basepairs = 1
    elif( feature == 'intron' ):
        basepairs = srs[11]
    elif( feature == 'intergenic' ):
        basepairs = int(srs[2]) - int(srs[1])
    
    if key in id_dict:
        # return empty
        return ( None , 0 , id_dict )
    else:
        id_dict[key] = 1
        return ( key , basepairs , id_dict )

if __name__ == "__main__":
    args = parse_args()
    line_tuples = []
    
    # Empty dictionary for checking SV IDs are counted once per event only
    SV_ID_Dict = {}
    
    intron_bases = 0
    intergenic_bases = 0
    
    f_intron = open( args.intron_bed , 'r')
    for line in f_intron:
        feature = 'intron'
        key, basepairs, SV_ID_Dict = summarise_line(line , feature , SV_ID_Dict )
        if key:
            line_tuples.append(feature)
            intron_bases += int(basepairs)

    f_intergenic = open( args.intergenic_bed , 'r')
    for line in f_intergenic:
        feature = 'intergenic'
        key, basepairs, SV_ID_Dict = summarise_line(line , feature , SV_ID_Dict )
        if key:
            line_tuples.append(feature)
            intergenic_bases += int(basepairs)

    value_counts = collections.Counter(line_tuples)
    summary = []
    value_dict = {
        "intron" : value_counts['intron'],
        "intergenic" : value_counts['intergenic'],
        "intron_bases": intron_bases,
        "intergenic_bases": intergenic_bases
    }
    summary.append(value_dict)

    print(json.dumps(summary, indent="\t"))
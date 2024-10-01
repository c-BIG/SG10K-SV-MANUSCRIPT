#!/usr/bin/env python

import fileinput
import json
import collections
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_bed", required=True, default=None ,
                        help="Input bed file")
    parser.add_argument("--interest_region", required=True, default='gene'  ,
                        help="Type of region-of-interest 1) ccre 2) intron 3) GENCODE GTF features eg. exon, gene etc ")
    args = parser.parse_args()
    
    return args


def summarise_line(line , interest_region , id_dict ):
    # line example =>   chr start stop SG10K-SV-ID svtype (from original SG10K-SV) - 5 column
    # a) Gene example      chr start stop UCSC_basic_ID 0 strand subtype (from processed bed file from UCSC) - 7 column
    # b) UCSC cCRE bed file       16 column
            
    srs = pd.Series(line.strip().split("\t"))
    
    svtype = srs[4]
    subtype = ''
    if ( interest_region == "ccre" ):
        # srs[20] = cCRE eg.CTCF, enhP etc
        subtype = srs[17]
    elif ( interest_region == "gene"):
        subtype = srs[11]
    else:
        subtype = interest_region

    # Generate SG10K-SV ID + svtype key-pair and check existence in dictionary 
    key = srs[3] + "/" + subtype
    if key in id_dict:
        # return empty
        return ( None , None , id_dict )
    else:
        id_dict[key] = 1
        return ( svtype , subtype , id_dict )

if __name__ == "__main__":
    args = parse_args()
    line_tuples = []
    
    # Empty dictionary for checking SV IDs are counted once per event only
    SV_ID_Dict = {}
    
    f = open( args.input_bed , 'r')
    for line in f:
        svtype = ""
        svtype, subtype, SV_ID_Dict = summarise_line(line , args.interest_region , SV_ID_Dict )
        if svtype:
            count_key = svtype + "/" + subtype
            line_tuples.append(count_key)

    value_counts = collections.Counter(line_tuples)
    summary = []
    for value, count in value_counts.items():
        svtype, subtype = value.split("/")
        value_dict = {
            "svtype": svtype,
            args.interest_region: subtype,
            "count": count
        }
        summary.append(value_dict)

    print(json.dumps(summary, indent="\t"))
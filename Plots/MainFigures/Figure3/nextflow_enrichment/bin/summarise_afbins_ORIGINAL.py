#!/usr/bin/env python

import sys
import argparse
import fileinput
import json
import collections
import pandas as pd

def summarise_line(line , grp_col , col , id_dict ):
    # line example =>   chr start stop SG10K-SV-ID svtype (from original SG10K-SV) - 5 column
    # a) Gene example      chr start stop UCSC_basic_ID 0 strand subtype (from processed bed file from UCSC) - 7 column
    # b) UCSC cCRE bed file       16 column
            
    srs = pd.Series(line.strip().split("\t"))
    
    svtype = srs[grp_col]
    subtype = srs[col]

    # Generate SG10K-SV ID + svtype key-pair and check existence in dictionary 
    key = srs[3] + "/" + subtype
    if key in id_dict:
        # return empty
        return ( None , None , id_dict )
    else:
        id_dict[key] = 1
        return ( svtype , subtype , id_dict )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--column', help='Which column to work on', required=True)
    parser.add_argument('--grouping_column', help='Which column to GROUP on', default = 4, nargs='?', const=1, type=int )
    parser.add_argument('--input_file', help='Input can be "-" for stdin', required=True)
    args = parser.parse_args()
    
    # Empty dictionary for checking SV IDs are counted once per event only
    SV_ID_Dict = {}
    line_tuples = []

    for line in sys.stdin:
        svtype = ""
        subtype = ""
        svtype, subtype, SV_ID_Dict = summarise_line(line, int(args.grouping_column) , int(args.column) , SV_ID_Dict )
        if svtype:
            count_key = svtype + "/" + subtype
            line_tuples.append(count_key)
            
    value_counts = collections.Counter(line_tuples)
    summary = []
    for value, count in value_counts.items():
        svtype, ccre = value.split("/")
        value_dict = {
            "svtype": svtype,
            "feature": ccre,
            "count": count
        }
        summary.append(value_dict)

    print(json.dumps(summary, indent="\t"))
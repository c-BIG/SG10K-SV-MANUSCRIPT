#!/bin/bash

IN_VCF=$1

bcftools query -f'%CHROM\t%POS0\t%SVSIZE\t%ID\t%SVTYPE\n' ${IN_VCF} | awk '{ if ( $3 < 1 ) {printf "%s\t%s\t%i\t%s\t%s\n",$1,$2,$2+1,$4,$5} else printf "%s\t%s\t%i\t%s\t%s\n",$1,$2,$2+$3,$4,$5}' 

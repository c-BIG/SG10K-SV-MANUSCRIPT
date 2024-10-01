#!/bin/bash

# Download data to scratch first and change directory
cd /mnt/volume1/SG10K-SV-NOTEBOOK/noncoding

conda activate sv_rng

nextflow run main.nf -profile local -resume --iterations 10000 \
--ccre true \
--bed '/mnt/volume1/SG10K-SV-NOTEBOOK/noncoding/resources/encodeCcreCombined.bed.gz' \
--results_dir 'result_r1.3_10k_ccre' --intermediate_dir 'intermediate_r1.3_10k_ccre' \
> log.10k_r1.3_ccre 2> logError.10k_r1.3_ccre
SUCCESS=$?

if [ $SUCCESS -eq 0 ]; then
echo OK
sleep 10s
rm -rf work

nextflow run main.nf -profile local -resume --iterations 10000 \
--ccre false \
--bed '/mnt/volume1/SG10K-SV-NOTEBOOK/noncoding/resources/gencode.v40.GENIC_COMPONENT.bed.gz' \
--results_dir 'result_r1.3_10k_gene' --intermediate_dir 'intermediate_r1.3_10k_gene' \
> log.10k_r1.3_gene 2> logError.10k_r1.3_gene
SUCCESS_PROC2=$?
	
	if [ $SUCCESS_PROC2 -eq 0 ]; then
	echo OK
	sleep 10s
	rm -rf work
	fi
fi

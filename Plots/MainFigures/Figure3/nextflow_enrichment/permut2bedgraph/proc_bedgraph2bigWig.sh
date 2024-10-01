#!/bin/bash

cd /mnt/volume1/SG10K-SV-NOTEBOOK/noncoding

THREAD=16
MEM='112G'
FOLDER="/mnt/volume1/SG10K-SV-NOTEBOOK/noncoding/result_r1.3_10k_ccre"
BIN="/mnt/volume1/SG10K-SV-NOTEBOOK/noncoding/bin"

AFBIN="common"

echo -e `date +"%F\t%T\t"`"[START] Cat  ${AFBIN}";
cat ${FOLDER}/shuffled_bed/${AFBIN}/*.bed > ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed
echo -e `date +"%F\t%T\t"`"[STOP] Cat ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Sort ${AFBIN}";
${BIN}/bedops_bin/sort-bed --max-mem ${MEM} --tmpdir $PWD/tmp ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed > ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed
bgzip --threads ${THREAD} ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed 
echo -e `date +"%F\t%T\t"`"[STOP] Sort ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Make bedgraph ${AFBIN}";
bedtools genomecov -bg -i ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed.gz -g resources/hg38.chrom.sizes \
>> ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph
echo -e `date +"%F\t%T\t"`"[STOP] Make bedgraph ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Make bigWig ${AFBIN}";
${BIN}/bedGraphToBigWig ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph \
resources/hg38.chrom.sizes \
${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bw
if [ $? -eq 0 ]; then
    echo "Make bigWig SUCCEEDED, doing cleanup"
	rm ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed
	rm ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph
	rm ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed.gz
else
    echo FAIL
	exit 3 
fi

AFBIN="rare"
echo -e `date +"%F\t%T\t"`"[START] Cat  ${AFBIN}";
cat ${FOLDER}/shuffled_bed/${AFBIN}/*.bed > ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed
echo -e `date +"%F\t%T\t"`"[STOP] Cat ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Sort ${AFBIN}";
${BIN}/bedops_bin/sort-bed --max-mem ${MEM} --tmpdir $PWD/tmp ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed > ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed
bgzip --threads ${THREAD} ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed 
echo -e `date +"%F\t%T\t"`"[STOP] Sort ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Make bedgraph ${AFBIN}";
bedtools genomecov -bg -i ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed.gz -g resources/hg38.chrom.sizes \
>> ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph
echo -e `date +"%F\t%T\t"`"[STOP] Make bedgraph ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Make bigWig ${AFBIN}";
${BIN}/bedGraphToBigWig ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph \
resources/hg38.chrom.sizes \
${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bw
if [ $? -eq 0 ]; then
    echo "Make bigWig SUCCEEDED, doing cleanup"
	rm ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed
	rm ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph
	rm ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed.gz
else
    echo FAIL
	exit 3 
fi

AFBIN="ultrarare"
echo -e `date +"%F\t%T\t"`"[START] Cat  ${AFBIN}";
cat ${FOLDER}/shuffled_bed/${AFBIN}/*.bed > ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed
echo -e `date +"%F\t%T\t"`"[STOP] Cat ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Sort ${AFBIN}";
${BIN}/bedops_bin/sort-bed --max-mem ${MEM} --tmpdir $PWD/tmp ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed > ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed
bgzip --threads ${THREAD} ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed 
echo -e `date +"%F\t%T\t"`"[STOP] Sort ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Make bedgraph ${AFBIN}";
bedtools genomecov -bg -i ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed.gz -g resources/hg38.chrom.sizes \
>> ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph
echo -e `date +"%F\t%T\t"`"[STOP] Make bedgraph ${AFBIN}";

echo -e `date +"%F\t%T\t"`"[START] Make bigWig ${AFBIN}";
${BIN}/bedGraphToBigWig ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph \
resources/hg38.chrom.sizes \
${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bw
if [ $? -eq 0 ]; then
    echo "Make bigWig SUCCEEDED, doing cleanup"
	rm ${FOLDER}/shuffled_bed/${AFBIN}.cat.bed
	rm ${FOLDER}/shuffled_bed/${AFBIN}.sorted.shuffled.bedgraph
	rm ${FOLDER}/shuffled_bed/${AFBIN}.sorted.bed.gz
else
    echo FAIL
	exit 3 
fi
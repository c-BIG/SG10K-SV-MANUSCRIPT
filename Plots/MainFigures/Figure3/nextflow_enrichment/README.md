# Enrichment analysis

1. Go to directory
```
cd SG10K-SV-NOTEBOOK/noncoding
```

2. Get exon file and ccre file from UCSC browser
    
**resources/UCSC_for_SV_enrichment/**

    a. gencode.v40.basic.exon.bed.gz
    b. encodeCcreCombined.bed.gz
    c. gencode.v40.GENIC_COMPONENT.bed.gz
    d. gatk-bundle_assembly38.genome.txt

3. Run nextflow script

a. CCRE
```
nextflow run main.nf -profile local -resume --iterations 10000 \
--ccre true --bed encodeCcreCombined.bed.gz \
--exons_bed gencode.v40.basic.exon.bed.gz \
--results_dir result  --intermediate_dir intermediate_results
```
b. Genic

```
nextflow run main_ALL_v2.nf -profile local -resume --iterations 10000 \
--ccre false \
--bed gencode.v40.GENIC_COMPONENT.bed.gz \
--results_dir result_gene --intermediate_dir intermediate_results_gene
```

4. Plot in Jupyter Notebook
```SV_noncoding_enrichment_GENCODE_V40.ipynb```

```SV_GENIC_enrichment_GENCODE_V40.ipynb```

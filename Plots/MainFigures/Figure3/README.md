# Figure 3a and 3b - Enrichment analysis

1. Go to directory
```
cd nextflow_enrichment
```

2. Get files of relevant version from UCSC browser

    a. gencode.v40.basic.exon.bed.gz
    b. encodeCcreCombined.bed.gz
    c. gencode.v40.GENIC_COMPONENT.bed.gz

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
nextflow run main.nf -profile local -resume --iterations 10000 \
--ccre false \
--bed gencode.v40.GENIC_COMPONENT.bed.gz \
--results_dir result_gene --intermediate_dir intermediate_results_gene
```

4. Plot in Jupyter Notebook
```Figure3a.ipynb```

```Figure3b.ipynb```
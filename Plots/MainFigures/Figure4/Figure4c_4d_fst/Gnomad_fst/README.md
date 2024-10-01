# Fixation index calculation of gnomAD-SV events

1. Go to directory
```
cd SG10K-SV-NOTEBOOK/Manuscript_Figures/figure_4c_fst/Gnomad_fst
```

2. Download Gnomad-SV **ACCESSIONED** data ( NOTE: hg19) with needed tags to calcalate Fst \

```
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd166/gnomad_v2.1_sv.sites.accessioned.vcf.gz
```

3. Run Nextflow to calculate Fst 
```
nextflow run main.nf --resume
```

4. Massage output
```
cd results/fst_csv
cat *.csv > ../gnomad_fst.csv

```
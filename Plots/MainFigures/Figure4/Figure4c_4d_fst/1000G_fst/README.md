# Fixation index calculation of gnomAD-SV events

1. Go to directory
```
cd SG10K-SV-NOTEBOOK/Manuscript_Figures/figure_4c_4d_fst/1000G_fst
```

2. Download 1000G VCF and sample grouping

```
mkdir data
cd data
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz.tbi

# Download annotation of samples from webpage below
wget https://www.internationalgenome.org/data-portal/data-collection/phase-3

3. Run Nextflow to calculate Fst 
```
nextflow run main.nf --resume
```

4. Massage output
```
cd results
cat fst_csv/*.csv > fst_1000G.csv
```
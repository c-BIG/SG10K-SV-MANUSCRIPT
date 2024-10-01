# Fixation index calculation of SG10K-SV events

1. Go to directory
```
cd SG10K-SV-MANUSCRIPT/Plots/MainFigures/Figure4c_4d_fst/SG10K-SV_fst
```

2. Download SG10K-SV matrixTable for ethnic group output (ethnicity.txt.npy)
```
mkdir data
cd data
# matrixTable for ethnic group output
aws s3 sync s3://<bucket>/release1.4/SG10K-SV-Release-1.4.mt/ . 
python ../bin/mt2vcf.py --in_mt SG10K-SV-Release-1.4.mt

# get output file "ethnicity.txt.npy"
```

2. Run Nextflow to calculate Fst 
```
nextflow run main.nf --resume
```
3. Massage output
```
cd results/fst_csv
echo "rsid,svtype,chrom,pos,end,svlen,fst_chi_ind,fst_chi_mal,fst_ind_mal,iter" \
> fst_hudson_r1.4_plus_bootstrap_1000.csv 
cat *.vcf.gz.csv | sed 1d >> fst_hudson_r1.4_plus_bootstrap_1000.csv

gzip fst_hudson_r1.4_plus_bootstrap_1000.csv

```
4. Move gzipped csv to top level for plotting
```
mv fst_hudson_r1.4_plus_bootstrap_1000.csv.gz ../../../fst_hudson_r1.4_plus_bootstrap_1000.csv.gz
```


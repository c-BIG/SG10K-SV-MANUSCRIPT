# SG10K-SV Figure 4c and 4d

1. Run Fst calculation of SG10K-SV events \
SG10K-SV-MANUSCRIPT/Plots/MainFigures/Figure4c_4d_fst/SG10K-SV_fst

2. Obtain required input file : fst_hudson_r1.3_plus_bootstrap_1000.csv.gz

3. Run bcftools to obtain "n_events_by_pop.csv" for Fst plots

`echo "rsid,svtype,af_alt_chinese,af_alt_indian,af_alt_malay\n" > n_events_by_pop.csv`

`bcftools query -f '%ID\t%SVTYPE\t%AF_CHINESE\t%AF_INDIAN\t%AF_MALAY\n' SG10K-SV-Release-1.4-HighConfidenceSV-withMetadata.vcf.bgz >> n_events_by_pop_r1.4.csv `

4. Run R script to generate Figure 4c and 4d \
a. Figure4c.R
b. Figure4d.R

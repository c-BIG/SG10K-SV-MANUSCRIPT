# read the file
ld_table <- read.table("SG10K-SV-SG10K-health-AllChr.filtered.final.ld", header=TRUE, sep="\t")

ld_table_filtered <- ld_table[which(ld_table$SNP_A != ld_table$SNP_B),]

ld_snp_sv <- ld_table_filtered[-grep("SG10K_SV", ld_table_filtered$SNP_B),]

# need to find the max LD per SV
library(dplyr)

maxld <- ld_snp_sv %>%
  group_by(SNP_A) %>%
  summarise(max = max(R2, na.rm=TRUE))

maxld$label <- "ALL"

ld_dup_snp <- maxld[grep("SG10K_SV_DUP", maxld$SNP_A),]
ld_del_snp <- maxld[grep("SG10K_SV_DEL", maxld$SNP_A),]
ld_ins_snp <- maxld[grep("SG10K_SV_INS", maxld$SNP_A),]

ld_dup_snp$label <- "DUP" 
ld_del_snp$label <- "DEL"
ld_ins_snp$label <- "INS"

library(ggplot2)
bigtable <- rbind(maxld, ld_dup_snp, ld_del_snp, ld_ins_snp)
ggplot(bigtable, aes(x=label, y=max, fill=label)) + geom_violin(width = 1) + 
  geom_boxplot(width=0.02, outlier.shape = NA) + 
  scale_fill_manual(values=c("#999999", '#4DBBD5FF','#00A087FF', '#DE8F44')) +
  theme_bw()

length(unique(ld_snp_sv$SNP_A))
#[1] 5697

6772-5697
#[1] 1075
1075/6772
#[1] 0.1587419

ld_r0.8<-ld_snp_sv[which(ld_snp_sv$R2>=0.8),]
dim(ld_r0.8)

length(unique(ld_r0.8$SNP_A))
#[1] 3909
length(unique(ld_r0.8$SNP_B))
#[1] 164992


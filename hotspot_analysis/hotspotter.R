library(GenomicRanges)
library(primatR)

bedfile <- read.table("SG10K-SV-Release-1.4-HighConfidenceSV-withMetadata.variantsonly.bed", header=FALSE, sep="\t")
head(bedfile)

res <-GRanges(seqnames = bedfile$V1,
            ranges = IRanges(start = bedfile$V2,
                             end = bedfile$V3,
                             names = bedfile$V5),
            type = bedfile$V4)

seqlengths(res) <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)

hotspots <- hotspotter(res, bw=200000,num.trial=10000, pval = 5e-03)

#export to bed
hotspot_df <- data.frame(seqnames=seqnames(hotspots), starts=start(hotspots)-1, ends=end(hotspots), pvalue=elementMetadata(hotspots)$pvalue, numevents=elementMetadata(hotspots)$num.events)
write.table(hotspot_df, "SG10K-SV-Release-1.4-HighConfidenceSV_hotspots_hotspotter_pvalue5e-3_10000iterations.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

bedgnomad <- read.table("nstd166.GRCh38.variant_region.bed", header=FALSE, sep="\t")
head(bedgnomad)

chromosomes_to_extract <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
finalbedgnomad <- bedgnomad[bedgnomad$V1 %in% chromosomes_to_extract,]

resgnomad <-GRanges(seqnames = finalbedgnomad$V1,
              ranges = IRanges(start = finalbedgnomad$V2,
                               end = finalbedgnomad$V3,
                               names = finalbedgnomad$V5),
              type = finalbedgnomad$V4)


seqlengths(resgnomad) <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)

hotspots_gnomad <- hotspotter(resgnomad, bw=200000,num.trial=10000, pval = 5e-03)

#export to bed
hotspot_gnomad_df <- data.frame(seqnames=seqnames(hotspots_gnomad), starts=start(hotspots_gnomad)-1, ends=end(hotspots_gnomad), pvalue=elementMetadata(hotspots_gnomad)$pvalue, numevents=elementMetadata(hotspots_gnomad)$num.events)
write.table(hotspot_gnomad_df, "gnomad-hotspots_hotspotter_pvalue5e-3_10000iterations.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

# Eichler long read data
bedeichler <- read.table("variants_freeze4_sv_insdel_alt.bed", header=FALSE, sep="\t")
head(bedeichler)

chromosomes_to_extract <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
chrom_to_extract <- paste("chr", chromosomes_to_extract,sep="")
finalbedeichler <- bedeichler[bedeichler$V1 %in% chrom_to_extract,]
resichler <-GRanges(seqnames = finalbedeichler$V1,
                    ranges = IRanges(start = finalbedeichler$V2,
                                     end = finalbedeichler$V3,
                                     names = finalbedeichler$V5),
                    type = finalbedeichler$V4)

resichler <- sortSeqlevels(resichler)
seqlengths(resichler) <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)

hotspots_1000glr<- hotspotter(resichler, bw=200000,num.trial=10000, pval = 5e-03)

#export to bed
hotspots_1000glr_df <- data.frame(seqnames=seqnames(hotspots_1000glr), starts=start(hotspots_1000glr)-1, ends=end(hotspots_1000glr), pvalue=elementMetadata(hotspots_1000glr)$pvalue, numevents=elementMetadata(hotspots_1000glr)$num.events)
write.table(hotspots_1000glr_df, "1000G_longread_hotspotter_pvalue5e-3_10000iterations.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)


# check the number of bps affected by hotspots
hotspot_df$size <- hotspot_df$ends - hotspot_df$starts
sum_of_hotspots_bp <- sum(hotspot_df$size)
sum_of_hotspots_bp
#211448993

# Check what number did gnomAD has
hotspot_gnomad_df$size <- hotspot_gnomad_df$ends - hotspot_gnomad_df$starts
sum(hotspot_gnomad_df$size)
# 379873456

# Check Eichler's paper
eichler_hotspots <- read.table("eichler_1000G.bed", sep="\t", header=FALSE)
eichler_hotspots$size <- eichler_hotspots$V3 - eichler_hotspots$V2
sum(eichler_hotspots$size)
# 284399556

# Annotate the regions to check if they overlap with the ends of the chromosome

# run the script containing all the functions
source("~/Functions.R")

chrom <- read.csv("hg38.chrom.sizes", sep="\t", header=FALSE)
head(chrom)
colnames(chrom) <- c("Chr", "Size")
chrom$mbp <- chrom$Size-5000000
head(chrom)

hotspot_df_chrsize <- add_chromosome_size(chrom, hotspot_df)
head(hotspot_df_chrsize)

### Load the centromeres coodrinates
cytobands <- read.csv("GRCH38_cytoband_coordinates.tsv", sep="\t", header=TRUE)
head(cytobands)
centromeres <- cytobands[cytobands$gieStain %in% "acen",]
centromere_df <- get_centromere_df(centromeres)
hotspot_df_chrsize_centro <-add_centromere_coordinates(centromeredf = centromere_df, hotspots = hotspot_df_chrsize)
head(hotspot_df_chrsize_centro)

library(dplyr)

# Annotate if the bins are in centromeres
hotspot_df_chrsize_centro_1 <- hotspot_df_chrsize_centro %>%  
  mutate(in_centromeres = case_when(ends < centromerestart | starts > centromereend ~ "not in centromere",                                    
                                    TRUE ~ "in centromere"  )
  )

hotspot_df_chrsize_centro_1_2 <- hotspot_df_chrsize_centro_1 %>%  
  mutate(in_start_end = case_when(ends < 5000000 | starts < 5000000 ~ "in 5mbp of chr start",                                  
                                  ends > chrsize5mbp | starts > chrsize5mbp ~ "in 5mbp of chr end",                                  
                                  TRUE ~ "not in"  )
  )


write.table(hotspot_df_chrsize_centro_1_2, "SG10K-SV-Release-1.4-HighConfidenceSV_hotspots_hotspotter_pvalue5e-3_10000iterations_with_annotations.bed", row.names = FALSE,  sep="\t")

# Find the number of intervals in chr ends
hotspots_in_chr_ends <- subset(hotspot_df_chrsize_centro_1_2, in_start_end != "not in")
hotspots_in_centro <- subset( hotspot_df_chrsize_centro_1_2,in_centromeres == "in centromere")
# get regions not in centromeres
hotspots_not_in_centro <- subset( hotspot_df_chrsize_centro_1_2,in_centromeres == "not in centromere")
hotspots_not_in_centro_chrends <- subset( hotspots_not_in_centro, in_start_end == "not in")
write.table(hotspots_not_in_centro_chrends, "SG10K-SV-Release-1.4-HighConfidenceSV_hotspots_hotspotter_pvalue5e-3_10000iterations_with_annotations_filtered_for_washU.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)



















# read csv file
truvari_csv <- read.table("truvari_summary_v4.csv", sep="\t", header=TRUE)
head(truvari_csv)

# extract 
gt2 <- truvari_csv[which(truvari_csv$svtype== "3SVTYPE"  & truvari_csv$Raw == "GT2" ) ,]
gt2

smoove <- truvari_csv[which(truvari_csv$svtype== "3SVTYPE"  & truvari_csv$tools == "SMOOVE" ) ,]
head(smoove)

combine_gt2_smoove <- rbind(gt2, smoove)
write.table(combine_gt2_smoove, "Truvari_supplmentarytable.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names=TRUE)

# extract manta, 15x
manta_15x_gt2 <- truvari_csv[which(truvari_csv$svtype== "3SVTYPE" & truvari_csv$coverage == "15x" & truvari_csv$Raw == "GT2" & truvari_csv$tools == "MANTA" ),]
manta_15x_gt2

# extract manta, 30x
manta_30x_gt2 <- truvari_csv[which(truvari_csv$svtype== "3SVTYPE" & truvari_csv$coverage == "30x" & truvari_csv$Raw == "GT2" & truvari_csv$tools == "MANTA"),]
manta_30x_gt2

sampleid <- manta_15x_gt2$sampleID 
sampleid

delta_df <- data.frame(matrix(NA, nrow = 34, ncol = 3))
colnames(delta_df) <- c("Sampleid", "delta_FP", "delta_FN")
delta_df

count = 0
for(i in sampleid)
{
  print(i)
  count <- count + 1
  manta_30x_gt2_FP <- manta_30x_gt2[which(manta_30x_gt2$sampleID == i), "FP"]
  manta_30x_gt2_FN <- manta_30x_gt2[which(manta_30x_gt2$sampleID == i), "FN"]
  
  manta_15x_gt2_FP <- manta_15x_gt2[which(manta_15x_gt2$sampleID == i), "FP"]
  manta_15x_gt2_FN <- manta_15x_gt2[which(manta_15x_gt2$sampleID == i), "FN"]
 
  delta_FP <- (manta_15x_gt2_FP/manta_30x_gt2_FP)*100
  delta_FN <- (1-(manta_30x_gt2_FN/manta_15x_gt2_FN))*100
  
  delta_df[count, "Sampleid"] <- i
  delta_df[count, "delta_FP"] <- delta_FP
  delta_df[count, "delta_FN"] <- delta_FN
  
}

delta_df
average_delta_FP <- mean(delta_df$delta_FP)
average_delta_FN <- mean(delta_df$delta_FN) 
average_delta_FP
average_delta_FN

# check del for 10 samples
manta_15x_gt2 <- truvari_csv[ truvari_csv$svtype== "DEL" & truvari_csv$coverage == "15x" & truvari_csv$Raw == "GT2" & truvari_csv$tools == "MANTA" ,]
s<-c("HG00512", "HG00514", "HG00864","HG01596","HG02492","HG03009","HG03683","HG03732","NA18939","NA20847")
manta_15x_gt2[which(manta_15x_gt2$sampleID %in% s),]
# Check ins for 10 samples
manta_15x_gt2 <- truvari_csv[ truvari_csv$svtype== "INS" & truvari_csv$coverage == "15x" & truvari_csv$Raw == "GT2" & truvari_csv$tools == "MANTA" ,]
manta_15x_gt2[which(manta_15x_gt2$sampleID %in% s),]


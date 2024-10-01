# Updated to run on Release 1.4

### 1. Prep work with bcftools

SG10K-SV release 1.4
```
# Download/Obtain SG10K-SV VCF
aws s3 cp s3://<bucket>/release1.4/SG10K-SV-Release-1.4.vcf.bgz . 
aws s3 cp s3://<bucket>/release1.4/SG10K-SV-Release-1.4.vcf.bgz.tbi .

# Get ethnic AC (Fig 4B) and all sample AF
echo -e "ID,SVTYPE,AF,AC_CHINESE,AC_MALAY,AC_INDIAN" > n_events_by_pop_r1.4.csv
bcftools query -f "%ID,%SVTYPE,%AF,%AC_CHINESE,%AC_MALAY,%AC_INDIAN" SG10K-SV-Release-1.4.vcf.bgz >> n_events_by_pop_r1.4.csv

```

### 2. R Command on RStudio
```
library("tidyverse")
library("janitor")
library("ggplot2")
library("repr")
library("ggsci")
library("hrbrthemes")
library("UpSetR")
library(svglite)

# local work dir
setwd("working_dir")

in_csv = "n_events_by_pop_r1.4.csv"
data = read_csv(in_csv, col_types = cols()) %>% 
clean_names()

# Plot - upsetr - all events
fdata = data %>% 
mutate(sg_chinese = ifelse(ac_chinese >= 1, 1, 0)) %>% 
mutate(sg_indian = ifelse(ac_indian  >= 1, 1, 0)) %>% 
mutate(sg_malay = ifelse(ac_malay  >= 1, 1, 0)) %>% 

mutate(sg_chinese = replace_na(sg_chinese, 0)) %>% 
mutate(sg_indian = replace_na(sg_indian, 0)) %>% 
mutate(sg_malay = replace_na(sg_malay, 0)) %>% 

mutate(n_categories = sg_chinese + sg_indian + sg_malay) %>% 

# discard events seen in others
filter(n_categories > 0)

## sanity checks
fdata %>% head

fdata %>% 
group_by(n_categories) %>% 
summarise(n = n()) %>% 
adorn_totals()

# n_categories     n
#            0   128
#            1 38019
#            2  9538
#            3 25350
#        Total 73035

ffdata = fdata

# PLOT
svglite("Fig_4b.svg")

options(repr.plot.width=8, repr.plot.height=5, repr.plot.res=300, repr.plot.quality=100)

ffdata %>% 
select(id, sg_chinese, sg_indian, sg_malay) %>% 
data.frame() %>% 
upset(order.by = "freq", empty.intersections = "on")

dev.off()

```



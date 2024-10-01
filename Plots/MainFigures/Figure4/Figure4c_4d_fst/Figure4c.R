# R command
library("tidyverse")
library("janitor")
library("ggplot2")
library("repr")
library("ggsci")
library("hrbrthemes")
library("UpSetR")
library("superheat")
setwd("work_dir")

# Working on SG10K-SV Release 1.4
data = read_csv('fst_hudson_r1.4_plus_bootstrap_1000.csv.gz', col_types = cols()) %>%
clean_names %>%
bind_rows()

# Cleanup data
data = data %>%
filter(!grepl('rsid', data$rsid))

# Set col to numeric 
data = data  %>% 
mutate_at( c("pos","end", "fst_chi_ind","fst_chi_mal" , "fst_ind_mal") , as.numeric )

# Function
col_max = function(x) {
    xs = x[c("fst_chi_ind", "fst_chi_mal", "fst_ind_mal")]
    r = names(xs[which.max(xs)])
    return(r)
}

calc_pvals = function(x) {
    
    # Using RSID directly
    # this_id = unique(x$id) %>% unlist()
    this_rsid = unique(x$rsid) %>% unlist()
    this_svtype = unique(x$svtype) %>% unlist()
    this_chrom = unique(x$chrom) %>% unlist()
    this_pos = unique(x$pos) %>% unlist()
    this_end = unique(x$end) %>% unlist()
    this_svlen = unique(x$svlen) %>% unlist()

    
    # observed
    o = x %>% 
        filter(iter == "observed")
    obs_fst_max = o$fst_max
    obs_fst_max_comparison = o$fst_max_comparison

    # bootstrap
    b = x %>% 
        filter(iter != "observed")

    n_iter = length(b$iter)
    
    n_higher_fst = b %>% 
        filter(fst_max >= obs_fst_max) %>% 
        count() %>% unlist()

    prop_higher_fst = n_higher_fst / n_iter

    # result
    r = data.frame(
        rsid = this_rsid,
        svtype = this_svtype,
        chrom = this_chrom,
        pos = this_pos,
        end = this_end,
        svlen = this_svlen,
        obs_fst_max = obs_fst_max,
        obs_fst_max_comparison = obs_fst_max_comparison,
        n_higher_fst = n_higher_fst,
        n_iter = n_iter,
        p_value = prop_higher_fst
    )

    return(r)
    
}

############
# Check output
# dimension of data
dim(data)
# Quick view
glimpse(data)
###########

# Data Cleanup
data = data %>% 
# fix negative/missing values
mutate(fst_chi_ind = ifelse(fst_chi_ind < 0 | is.na(fst_chi_ind), 0, fst_chi_ind)) %>% 
mutate(fst_chi_mal = ifelse(fst_chi_mal < 0 | is.na(fst_chi_mal), 0, fst_chi_mal)) %>% 
mutate(fst_ind_mal = ifelse(fst_ind_mal < 0 | is.na(fst_ind_mal), 0, fst_ind_mal)) %>% 
# calculate FST_max and keep track of which comparison it belongs to
mutate(fst_max = pmax(fst_chi_ind, fst_chi_mal, fst_ind_mal)) %>% 
mutate(fst_max_comparison = apply(., 1, col_max))


# p-value calculation
N_PERMUTATIONS = 1000

pvals = data %>% 
group_by(rsid) %>% 
group_map(~ calc_pvals(.x), .keep = TRUE) %>% 
bind_rows() %>% 
mutate(p_value = ifelse(p_value == 0, 1/N_PERMUTATIONS, p_value))
     
# Get FDR
pvals = pvals %>% 
mutate(fdr = p.adjust(p_value, method = "fdr"))

# Significant event
MEAN_FST = mean(pvals$obs_fst_max)

# Load functional annotations + call rate
# Get SVTK2 annotatin file
f = "svtk_annotation_for_fst_computation.tsv"
svtk = read_tsv(file = f, col_types = cols()) %>% 
clean_names() %>% 
mutate(functional = ifelse(svtk_annot2 != "[]", TRUE, FALSE)) %>% 
select(rsid, call_rate, svtk_annot2, functional)

# Check
svtk %>% head

# merge
data_annot = pvals %>% 
merge(svtk, by = "rsid", all.x = TRUE) %>% 
mutate(category = case_when(
  (significant == TRUE & functional == TRUE) ~ 'significant and functional',
  (significant == TRUE & functional == FALSE) ~ 'significant',
  significant == FALSE ~ 'not significant',
  TRUE ~ as.character(significant)
))

data_annot %>% head

######## PLOT
######## Figure 4c
options(repr.plot.width=12, repr.plot.height=6, repr.plot.res=300, repr.plot.quality=100)

pdf("fig_4c.pdf")

data_annot %>% 
filter(category != "significant and functional") %>% 
ggplot(aes(x = call_rate, y = obs_fst_max, color = category)) +
geom_point(alpha = 0.25, size = 0.9) +
geom_point(data = data_annot %>% filter(category == "significant and functional")) +
scale_color_jama() +
ylim(0, 1) +
labs(y= "Observed Max Fst", x = "Call Rate") +
theme_bw() +
theme(
    legend.position = c(.85, .90),
    legend.box.background = element_rect(color="grey", size=0.75),
    legend.background = element_rect(fill="white")
)

dev.off()

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
# Get SVTK2 annotation file
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

######## PLOT
######## Figure 4d
#  Filter significant SV
sig_sv = data_annot %>% 
filter(significant == TRUE & functional == TRUE ) 

dim(sig_sv)

# Import n_events_by_pop.csv
af_csv = "n_events_by_pop_r1.4.csv"
af = read_csv(af_csv, col_types = cols()) 

af <- af %>% 
    mutate(af_alt_indian = replace_na(af_alt_indian, 0)) %>%
	mutate_at( c("af_alt_chinese","af_alt_indian", "af_alt_malay") , as.numeric )

# Merge sig_sv with data
sig_sv = sig_sv %>% 
left_join( af %>% select(-svtype), by = "rsid")

# Check
sig_sv %>% 
group_by(svtype) %>% 
count() %>% 
adorn_totals()

# filter
h_data = sig_sv %>% 
select(rsid , chrom, pos, end, svtype, af_alt_chinese, af_alt_indian, af_alt_malay, obs_fst_max, obs_fst_max_comparison) %>% 
unique() %>% 
mutate(filter_col = pmax(af_alt_chinese, af_alt_indian, af_alt_malay)) %>% 
filter(filter_col > 0.05)

# Check
h_data %>% 
group_by(svtype) %>% 
count() %>% 
adorn_totals()
     
## PLOT heatmap

mh_data = h_data %>% 
select(af_alt_chinese, af_alt_indian, af_alt_malay) %>% 
as.matrix()

# Note spaces at the end of row names to include margins in the heatmap
colnames(mh_data) = c("CHIESE_AF", "INDIAN_AF", "MALAY_AF")
rownames(mh_data) = h_data$rsid

###

pdf("Fig_4d.pdf")

options(repr.plot.width=10, repr.plot.height=40, repr.plot.res=1200, repr.plot.quality=100)
superheat(
	## base matrix
	mh_data,
	## left labels
	left.label.size = 1.5,
	left.label.text.alignment = "right",
	left.label.text.size = 1.5,
	left.label.col = "white",
	## bottom labels
	bottom.label.size = 0.02,
	bottom.label.text.size = 2,
	bottom.label.col = "white",
	## colors & legend
	heat.pal = c("#8ecae6", "#023047", "#ca6702", "#ffb703"),
	heat.pal.values = c(0, 0.1, 0.4, 1),
	heat.lim = c(0, 1),
	legend = TRUE,
	## row grouping
	scale = FALSE,
	pretty.order.rows = TRUE,
	# cell labels
	round(as.matrix(mh_data), 2),
	X.text.size = 1.8,
	X.text.col = "white",
	## barplot with fst_max
	yr = h_data$obs_fst_max,
	yr.plot.type = "bar",
	yr.plot.size = 0.66,
	yr.axis.name = "fst_max",
	yr.axis.size = 5,
	yr.axis.name.size = 5,
	yr.obs.col = rep("turquoise4", nrow(h_data)),
)

dev.off()

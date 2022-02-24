# FIND INDEPENDENT SIGNALS IN SUMMARY STATISTICS FROM MATTIA - 20220223

import sys
import pandas as pd
PAD = 1500000
data = pd.read_csv(sys.argv[1], sep='\t')
data = data.sort_values('all_inv_var_meta_p').reset_index(drop=True)
data['lead'] = False
leads = pd.DataFrame(columns=list(data.columns))
i=0
while len(data) > 0:
    data.loc[0, 'lead'] = True
    leads.loc[i] = data.loc[0]
    i = i+1
    chr = data.loc[0, '#CHR']
    pos = data.loc[0, 'POS']
    data = data[(data['#CHR'] != chr) | (data['POS'] < pos-PAD) | (data['POS'] > pos+PAD)].reset_index(drop=True)
    
leads.to_csv(sys.argv[1] + '.leads.txt', na_rep='NA', sep='\t', index=False)

##################################################
########## POLYGENIC SCORES - 20220215
##################################################

# PRS DATA
cp /finngen/library-red/finngen_R8/prs_1.0/data/*.sscore ./prs

# LOAD PRS DATA
# FIX FILES THAT BEGIN WITH NUMBERS IN TERMINAL
#mv finngen_R8_20171017_MW_eGFR_overall_EA_nstud42.dbgap_v2.txt.sscore finngen_R8_MW_eGFR_overall_EA_nstud42.dbgap_v2.txt.sscore
prs_metadata = fread(file = paste0(path, "prs_metadata_R8.csv"), header = TRUE, stringsAsFactors = FALSE)
prs_metadata$basename = str_replace_all(prs_metadata$basename, "[-+]", "_")
prs_file_list = prs_metadata$scorename
prs_list = prs_metadata$basename
for (i in 1:length(prs_file_list)) {
	prs_temp = fread(file = paste0(path, "prs/", prs_file_list[i]), header = TRUE, stringsAsFactors = FALSE)
	prs_temp = prs_temp[,c("IID", "SCORE1_AVG")]
	colnames(prs_temp) = c("finngenid", prs_list[i])
	total = merge(total, prs_temp, by = "finngenid", all.x = TRUE)
	print(prs_list[i])
}
dim(total) # 388465

##################################################
########## DO ASSOCIATION ANALYSIS - 20220215
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

total_new = fread(file = paste0(path, "healthcare_cost_20211222.txt"), header = TRUE, stringsAsFactors = FALSE)

prs_metadata = fread(file = paste0(path, "prs_metadata_R8.csv"), header = TRUE, stringsAsFactors = FALSE)
prs_metadata$basename = str_replace_all(prs_metadata$basename, "[-+]", "_")
prs_list = prs_metadata$basename

# MODEL DRUG COST ~ PRS - 20220215
mod_res = NULL
for (i in 1:length(prs_list)) {
	print(prs_list[i])
	temp = total_new[which(!is.na(total_new[[prs_list[i]]])),]
	# temp = total_new[which(!is.na(total_new[[prs_list[i]]]) & total_new$rx_statin == 1 & total_new$fu_end_age >= 60),]
	temp$total_cost = log(temp$total_cost + 1)
	#temp$total_cost = log((temp$total_cost / temp$length_of_fu) + 1)
	temp$prs = scale(temp[[prs_list[i]]])
	mod = lm(formula = total_cost ~ prs + birth_year + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = temp)
	coef = coef(summary(mod))["prs", "Estimate"]
	if (coef < 0) {
		temp$prs = -scale(temp[[prs_list[i]]])
		mod = lm(formula = total_cost ~ prs + birth_year + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = temp)
		coef = coef(summary(mod))["prs", "Estimate"]
	}
	ste = coef(summary(mod))["prs", "Std. Error"]
	pval = coef(summary(mod))["prs", "Pr(>|t|)"]
	mod_res = rbind(mod_res, c(prs_list[i], coef, ste, pval))
}
mod_res = as.data.frame(mod_res)
colnames(mod_res) = c("prs", "coef", "ste", "pval")
mod_res$coef = as.numeric(mod_res$coef)
mod_res$ste = as.numeric(mod_res$ste)
mod_res$pval = as.numeric(mod_res$pval)
mod_res$beta = (exp(mod_res$coef)-1)*100
mod_res$beta_lower = (exp(mod_res$coef - 1.96 * mod_res$ste)-1)*100
mod_res$beta_upper = (exp(mod_res$coef + 1.96 * mod_res$ste)-1)*100
mod_res = merge(prs_metadata, mod_res, by.x = "basename", by.y = "prs", all = TRUE)

fwrite(mod_res, paste0(path, "association_analysis_20220131.txt"), row.names = FALSE, col.names = TRUE)
mod_res = fread(file = paste0(path, "association_analysis_20220131.txt"), header = TRUE, stringsAsFactors = FALSE)

# PLOT MODEL RESULTS - 20220215
mod_temp = mod_res[which(mod_res$category == "Measurement" & mod_res$include == 1),]
ggplot() +
	geom_point(data = mod_temp, mapping = aes(x = reorder(basename, -beta), y = beta, color = ifelse(pval < 0.05/nrow(mod_temp), "Yes", "No")), size = 1) +
	coord_flip () +
	geom_errorbar(data = mod_temp, mapping = aes(x = reorder(basename, -beta), ymin = beta_lower, ymax = beta_upper, color = ifelse(pval < 0.05/nrow(mod_temp), "Yes", "No")), width = 0.1, size = 1) +
	scale_x_discrete(limits = mod_temp$basename[order(-mod_temp$beta)], labels = mod_temp$phenotype_label[order(-mod_temp$beta)]) + 
	scale_y_continuous(breaks = seq(floor(min(mod_temp$beta_lower)), ceiling(max(mod_temp$beta_upper)), 1), limits = c(floor(min(mod_temp$beta_lower)), ceiling(max(mod_temp$beta_upper)))) + 
	labs(x = "", y = "Percent Change in Healthcare Cost", color = "Significant?") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# MODEL DRUG COST ~ PRS BY AGE INTERVAL - 20220215
age_low = seq(0, 80, 20)
age_high = seq(20, 100, 20)
mod_res = NULL
for (j in 1:length(age_low)) {
	print(paste0(age_low[j], " to ", age_high[j]))
	for (i in 1:length(prs_list)) {
		print(prs_list[i])		
		temp = long_new %>% filter(event_age >= age_res$age_low[j] & event_age < age_res$age_high[j]) %>% group_by(finngenid) %>% summarise(sum_cost = sum(cost, na.rm = TRUE))
		temp = merge(total_new, temp, by = "finngenid")
		temp = temp[which(!is.na(temp[[prs_list[i]]])),]
		temp$person_time = ifelse(temp$fu_end_age < age_res$age_high[j], temp$fu_end_age - age_res$age_low[j], 5)
		temp$total_cost = log(temp$total_cost + 1)
		temp$prs = scale(temp[[prs_list[i]]])
		mod = lm(formula = total_cost ~ prs + person_time + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = temp)
		coef = coef(summary(mod))["prs", "Estimate"]
		if (coef < 0) {
			temp$prs = -scale(temp[[prs_list[i]]])
			mod = lm(formula = total_cost ~ prs + person_time + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = temp)
			coef = coef(summary(mod))["prs", "Estimate"]
		}
		ste = coef(summary(mod))["prs", "Std. Error"]
		pval = coef(summary(mod))["prs", "Pr(>|t|)"]
		mod_res = rbind(mod_res, c(age_low[j], age_high[j], prs_list[i], coef, ste, pval))
	}
}
mod_res = as.data.frame(mod_res)
colnames(mod_res) = c("low", "high", "prs", "coef", "ste", "pval")
mod_res$low = as.numeric(mod_res$low)
mod_res$high = as.numeric(mod_res$high)
mod_res$coef = as.numeric(mod_res$coef)
mod_res$ste = as.numeric(mod_res$ste)
mod_res$pval = as.numeric(mod_res$pval)
mod_res$beta = (exp(mod_res$coef)-1)*100
mod_res$beta_lower = (exp(mod_res$coef - 1.96 * mod_res$ste)-1)*100
mod_res$beta_upper = (exp(mod_res$coef + 1.96 * mod_res$ste)-1)*100
mod_res = merge(prs_metadata, mod_res, by.x = "basename", by.y = "prs", all = TRUE)

# PLOT MODEL RESULTS - 20220215
mod_temp = mod_res[which(mod_res$category == "Behavior" & mod_res$include == 1),]
ggplot() +
	geom_line(data = mod_temp, mapping = aes(x = low, y = beta, color = basename), size = 1) +
	geom_point(data = mod_temp, mapping = aes(x = low, y = beta, shape = ifelse(pval < 0.05/nrow(mod_res), "Yes", "No")), size = 3) +
	geom_errorbar(data = mod_temp, mapping = aes(x = low, ymin = beta_lower, ymax = beta_upper, color = basename), width = 0.1, size = 1) +
	scale_x_continuous(labels = age_low, breaks = age_low) +
	scale_shape_manual(values = c(1, 16)) +
	scale_color_discrete(limits = mod_temp$basename, labels = mod_temp$phenotype_label) + 
  labs(x = "Age Interval", y = "Percent Change in Healthcare Cost", color = "Phenotype", shape = "Significant?") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# MODEL DRUG COST ~ PRS BY SEX - 20220215
mod_res = NULL
for (i in 1:length(prs_list)) {
	temp = total_new[which(!is.na(total_new[[prs_list[i]]])),]
	temp$total_cost = log(temp$total_cost + 1)
	temp$prs = scale(temp[[prs_list[i]]])
	mod = lm(formula = total_cost ~ prs * sex + birth_year + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = temp)
	coef = coef(summary(mod))["prs:sexmale", "Estimate"]
	ste = coef(summary(mod))["prs:sexmale", "Std. Error"]
	pval = coef(summary(mod))["prs:sexmale", "Pr(>|t|)"]
	mod_res = rbind(mod_res, c(prs_list[i], coef, ste, pval))
}
mod_res = as.data.frame(mod_res)
colnames(mod_res) = c("prs", "coef", "ste", "pval")
mod_res$coef = as.numeric(mod_res$coef)
mod_res$ste = as.numeric(mod_res$ste)
mod_res$pval = as.numeric(mod_res$pval)
mod_res$beta = (exp(mod_res$coef)-1)*100
mod_res$beta_lower = (exp(mod_res$coef - 1.96 * mod_res$ste)-1)*100
mod_res$beta_upper = (exp(mod_res$coef + 1.96 * mod_res$ste)-1)*100
mod_res = merge(prs_metadata, mod_res, by.x = "basename", by.y = "prs", all = TRUE)
mod_res[which(mod_res$pval < 0.05 / nrow(mod_res)),c("phenotype_label", "coef", "ste", "pval")]

# DEFINE IMPORTANT PHENOTYPES - 20220215

pheno_list = c("cad.add.160614.website.txt", 
	"UKB_ICBPmeta750k_SBPsummaryResults.txt",
	"SavageJansen_2018_intelligence_metaanalysis_formatted.txt",
	"MTAG_CP.to10K.txt",
	"PGC3_SCZ_wave3_public_without_frequencies.v1.tsv",
	"PGC_UKB_depression_genome_wide.txt",
	"GSCAN_CigarettesPerDay.txt",
	"sumstats_neuroticism_ctg_format_formatted.txt",
	"MTAG_EA.to10K.txt",
	"lifegen_phase2_bothpl_alldr_2017_09_18.tsv",
	"GSEM.GWAS.EXTERNALIZING.SHARE.v20191014.txt",
	"Mahajan.NatGenet2018b.T2D.European.txt",
	"AD_sumstats_Jansenetal.txt",
	"GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt",
	"GSCAN_DrinksPerWeek.txt",
	"biomarkers_30890_both_sexes_irnt.tsv",
	"Meta_analysis_Locke_et_al_UKBiobank_2018_bmi.txt",
	"continuous_LDLC_both_sexes_medadj_irnt.tsv",
	"chronic_pain_bgen.stats")

# PLOT DISTRIBUTION OF HEALTHCARE COSTS BY PRS PERCENTILE - 20220215

temp = total_new %>% mutate_at(pheno_list, ntile, 10) %>% pivot_longer(pheno_list, names_to = "prs", values_to = "decile") %>% group_by(prs, decile, sex) %>% summarise(cost = mean(total_cost, na.rm = TRUE), 
	confint_lower = t.test(total_cost)$conf.int[1],
	confint_upper = t.test(total_cost)$conf.int[2])
temp = as.data.frame(temp)

ggplot() +
	geom_point(data = temp, mapping = aes(x = decile, y = cost, color = sex), size = 3) +
	geom_line(data = temp, mapping = aes(x = decile, y = cost, color = sex), size = 1) +
	#geom_errorbar(data = temp, mapping = aes(x = decile, ymin = confint_lower, ymax = confint_upper, color = sex), width = 0.1, size = 1) +
	scale_x_discrete(limits = seq(1, 10, 1), labels = seq(1, 10, 1)) + 
	labs(x = "PRS Decile", y = "Median Healthcare Cost", color = "Sex") + 
	#facet_wrap(~prs) + 
	facet_wrap(~prs, scales = "free") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

##################################################
########## PLOT DISTRIBUTION OF HEALTHCARE COSTS BY AGE INTERVAL - 20220215
##################################################

age_res = data.frame(age_low = seq(0, 95, 5),
	age_high = seq(5, 100, 5),
	female_time = NA,
	female_cost = NA,
	male_time = NA,
	male_cost = NA,
	stringsAsFactors = FALSE)

# PLOT DISTRIBUTION OF HEALTHCARE COSTS BY PRS PERCENTILE - 20220215
#id_prs = unique(total$finngenid[which(total$cad.add.160614.website > quantile(total$cad.add.160614.website, 0.90, na.rm = TRUE))])

# PLOT DISTRIBUTION OF HEALTHCARE COSTS BY HEALTHCARE COSTS PERCENTILE - 20220215
#id_cost = unique(total_new$finngenid[which(total_new$total_cost > quantile(total_new$total_cost, 0.99, na.rm = TRUE))])

calculate_person_time_cost = function(data, column) {
	for (i in 1:nrow(age_res)) {
		print(paste0(age_res$age_low[i], " to ", age_res$age_high[i]))
		# Get subset of individuals followed up during ith age interval
		temp = endpoint_new %>% filter(fu_start_age < age_res$age_high[i] & fu_end_age >= age_res$age_low[i])
		# Change from sex to PRS percentile
		#temp$sex = ifelse(temp$finngenid %in% id_prs, "female", "male")
		# Change from sex to healthcare cost percentile
		#temp$sex = ifelse(temp$finngenid %in% id_cost, "female", "male")
		# Calculate person-time observed for each individual
		temp$person_time = ifelse(temp$fu_end_age < age_res$age_high[i], temp$fu_end_age - age_res$age_low[i], 5)
		# Calculate total person-time
		sum_female_time = temp %>% filter(sex == "female") %>% summarise(sum_female_time = sum(person_time, na.rm = TRUE)) %>% pull(sum_female_time)
		sum_male_time = temp %>% filter(sex == "male") %>% summarise(sum_male_time = sum(person_time, na.rm = TRUE)) %>% pull(sum_male_time)
		# Calculate total cost
		sum_female_cost = data %>% filter(finngenid %in% temp$finngenid[which(temp$sex == "female")] & event_age >= age_res$age_low[i] & event_age < age_res$age_high[i]) %>% summarise(sum_cost = sum(!!as.name(column), na.rm = TRUE)) %>% pull(sum_cost)
		sum_male_cost = data %>% filter(finngenid %in% temp$finngenid[which(temp$sex == "male")] & event_age >= age_res$age_low[i] & event_age < age_res$age_high[i]) %>% summarise(sum_cost = sum(!!as.name(column), na.rm = TRUE)) %>% pull(sum_cost)
		# Set follow-up time
		age_res$female_n[i] = length(unique(temp$finngenid[which(temp$person_time > 0 & temp$sex == "female")]))
		age_res$male_n[i] = length(unique(temp$finngenid[which(temp$person_time > 0 & temp$sex == "male")]))
		age_res$female_time[i] = sum_female_time
		age_res$female_cost[i] = sum_female_cost
		age_res$male_time[i] = sum_male_time
		age_res$male_cost[i] = sum_male_cost
	}
	age_res$person_n = age_res$female_n + age_res$male_n
	age_res$person_time = age_res$female_time + age_res$male_time
	age_res$person_cost = age_res$female_cost + age_res$male_cost
	age_res$log_female_cost = log(age_res$female_cost + 1)
	age_res$log_male_cost = log(age_res$male_cost + 1)
	age_res$log_person_cost = log(age_res$person_cost + 1)
	age_res$cost_per_female_time = age_res$female_cost / age_res$female_time
	age_res$cost_per_male_time = age_res$male_cost / age_res$male_time
	age_res$cost_per_person_time = age_res$person_cost / age_res$person_time
	age_res$log_cost_per_female_time = age_res$log_female_cost / age_res$female_time
	age_res$log_cost_per_male_time = age_res$log_male_cost / age_res$male_time
	age_res$log_cost_per_person_time = age_res$log_person_cost / age_res$person_time
	return(age_res)
}

drug_age_res = calculate_person_time_cost(drug_temp, "drug_cost")
hilmo_age_res = calculate_person_time_cost(hilmo_temp, "cost")
avohilmo_age_res = calculate_person_time_cost(avohilmo_temp, "cost")
total_age_res = calculate_person_time_cost(long_new, "cost")

age_res = drug_age_res
age_res = hilmo_age_res
age_res = avohilmo_age_res
age_res = total_age_res

ggplot() +
	geom_line(data = age_res, mapping = aes(x = age_low, y = cost_per_person_time), color = "black", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = cost_per_female_time), color = "red", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = cost_per_male_time), color = "blue", size = 1) +
	labs(title = "Total", x = "Age", y = "Cost per Person-Time") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_line(data = age_res, mapping = aes(x = age_low, y = log_cost_per_person_time), color = "black", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = log_cost_per_female_time), color = "red", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = log_cost_per_male_time), color = "blue", size = 1) +
	labs(title = "Total", x = "Age", y = "ln(Cost) per Person-Time") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_line(data = age_res, mapping = aes(x = age_low, y = person_n), color = "black", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = female_n), color = "red", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = male_n), color = "blue", size = 1) +
	labs(title = "Total", x = "Age", y = "Number of Persons Contributing Person-Time") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# PLOT DISTRIBUTION OF SUBTYPE PERCENT OF HEALTHCARE COSTS BY AGE INTERVAL - 20220215
comb_age_res = as.data.frame(rbind(cbind(drug_age_res$age_low, drug_age_res$age_high, drug_age_res$female_cost, drug_age_res$male_cost, drug_age_res$person_cost, "Drug"),
	cbind(hilmo_age_res$age_low, hilmo_age_res$age_high, hilmo_age_res$female_cost, hilmo_age_res$male_cost, hilmo_age_res$person_cost, "HILMO"),
	cbind(avohilmo_age_res$age_low, avohilmo_age_res$age_high, avohilmo_age_res$female_cost, avohilmo_age_res$male_cost, avohilmo_age_res$person_cost, "AvoHILMO")))
colnames(comb_age_res) = c("age_low", "age_high", "female_cost", "male_cost", "person_cost", "class")
comb_age_res$age_low = as.numeric(comb_age_res$age_low)
comb_age_res$age_high = as.numeric(comb_age_res$age_high)
comb_age_res$female_cost = as.numeric(comb_age_res$female_cost)
comb_age_res$male_cost = as.numeric(comb_age_res$male_cost)
comb_age_res$person_cost = as.numeric(comb_age_res$person_cost)

ggplot() + 
	geom_bar(data = comb_age_res[which(comb_age_res$class != "Total"),], mapping = aes(x = age_low, y = person_cost, fill = class), position = "fill", stat = "identity") + 
	scale_x_discrete(limits = seq(0,100,10), labels = seq(0,100,10)) + 
	labs(x = "Age", y = "Percent of Total Cost", fill = "") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# PLOT DISTRIBUTION OF SUBTYPE PERCENT OF HEALTHCARE COSTS BY AGE INTERVAL AND SEX - 20220215
comb_age_sex_res = as.data.frame(rbind(cbind(drug_age_res$age_low, drug_age_res$age_high, drug_age_res$female_cost, total_age_res$female_cost, "Drug", "F"),
	cbind(drug_age_res$age_low, drug_age_res$age_high, drug_age_res$male_cost, total_age_res$male_cost, "Drug", "M"),
	cbind(drug_age_res$age_low, drug_age_res$age_high, drug_age_res$person_cost, total_age_res$person_cost, "Drug", "X"),
	cbind(hilmo_age_res$age_low, hilmo_age_res$age_high, hilmo_age_res$female_cost, total_age_res$female_cost, "HILMO", "F"),
	cbind(hilmo_age_res$age_low, hilmo_age_res$age_high, hilmo_age_res$male_cost, total_age_res$male_cost, "HILMO", "M"),
	cbind(hilmo_age_res$age_low, hilmo_age_res$age_high, hilmo_age_res$person_cost, total_age_res$person_cost, "HILMO", "X"),
	cbind(avohilmo_age_res$age_low, avohilmo_age_res$age_high, avohilmo_age_res$female_cost, total_age_res$female_cost, "AvoHILMO", "F"),
	cbind(avohilmo_age_res$age_low, avohilmo_age_res$age_high, avohilmo_age_res$male_cost, total_age_res$male_cost, "AvoHILMO", "M"),
	cbind(avohilmo_age_res$age_low, avohilmo_age_res$age_high, avohilmo_age_res$person_cost, total_age_res$person_cost, "AvoHILMO", "X")))
colnames(comb_age_sex_res) = c("age_low", "age_high", "cost", "total_cost", "class", "sex")
comb_age_sex_res$age_low = as.numeric(comb_age_sex_res$age_low)
comb_age_sex_res$age_high = as.numeric(comb_age_sex_res$age_high)
comb_age_sex_res$cost = as.numeric(comb_age_sex_res$cost)
comb_age_sex_res$total_cost = as.numeric(comb_age_sex_res$total_cost)
comb_age_sex_res$percent = comb_age_sex_res$cost/comb_age_sex_res$total_cost

ggplot() + 
	geom_bar(data = comb_age_sex_res[which(comb_age_sex_res$sex != "X"),], mapping = aes(x = sex, y = percent, fill = class), position = "stack", stat = "identity") + 
	facet_wrap(~age_low, nrow = 1) +
	labs(x = "", y = "Percent of Total Cost", fill = "") + 
	theme(panel.spacing = unit(0, "lines"), axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# CHECK N VUMBERS - 20220214
length(unique(long_all$finngenid)) # 355907
length(unique(long_new$finngenid)) # 352651
length(unique(drug_all$finngenid)) # 354734
length(unique(drug_new$finngenid)) # 354734
length(unique(drug_temp$finngenid)) # 351593
dim(drug_cost) # 351593
length(unique(hilmo_all$finngenid)) # 352632
length(unique(hilmo_new$finngenid)) # 352632
length(unique(hilmo_temp$finngenid)) # 347065
dim(hilmo_cost) # 347065
length(unique(avohilmo_all$finngenid)) # 323016

# CHECK NA VALUES - 20220214
if(!require(Hmisc)) {install.packages("Hmisc"); library(Hmisc)}
length(unique(total_new$finngenid[which(total_new$finngenid %nin% drug_cost$finngenid & total_new$finngenid %nin% hilmo_cost$finngenid & total_new$finngenid %nin% avohilmo_cost$finngenid)]))
length(total_new$finngenid[which(is.na(total_new$drug_cost) & is.na(total_new$hilmo_cost) & is.na(total_new$avohilmo_cost))])
length(unique(total_new$finngenid[which(is.na(total_new$drug_cost))]))
length(unique(total_new$finngenid[which(is.na(total_new$hilmo_cost))]))
length(unique(total_new$finngenid[which(is.na(total_new$avohilmo_cost))]))
length(unique(total_new$finngenid[which(is.na(total_new$drug_cost) & is.na(total_new$hilmo_cost) & is.na(total_new$avohilmo_cost))]))
length(unique(total_new$finngenid[which(!is.na(total_new$drug_cost) & !is.na(total_new$hilmo_cost) & !is.na(total_new$avohilmo_cost))]))
length(unique(total_new$finngenid[which(!is.na(total_new$drug_cost) & is.na(total_new$hilmo_cost) & is.na(total_new$avohilmo_cost))]))
length(unique(total_new$finngenid[which(is.na(total_new$drug_cost) & !is.na(total_new$hilmo_cost) & is.na(total_new$avohilmo_cost))]))
length(unique(total_new$finngenid[which(is.na(total_new$drug_cost) & is.na(total_new$hilmo_cost) & !is.na(total_new$avohilmo_cost))]))
length(unique(total_new$finngenid[which(!is.na(total_new$drug_cost) & !is.na(total_new$hilmo_cost) & is.na(total_new$avohilmo_cost))]))
length(unique(total_new$finngenid[which(!is.na(total_new$drug_cost) & is.na(total_new$hilmo_cost) & !is.na(total_new$avohilmo_cost))]))
length(unique(total_new$finngenid[which(is.na(total_new$drug_cost) & !is.na(total_new$hilmo_cost) & !is.na(total_new$avohilmo_cost))]))

# PREDICT MISSING AVOHILMO COSTS - 20220214
if(!require(quantreg)) {install.packages("quantreg"); library(quantreg)}
avohilmo_pred = merge(endpoint_new, avohilmo_cost, by = "finngenid")
avohilmo_pred$avohilmo_cost_per_person_day = avohilmo_pred$avohilmo_cost / avohilmo_pred$person_day_since_2011
avohilmo_pred$avohilmo_cost_per_person_day[which(avohilmo_pred$avohilmo_cost_per_person_day == Inf)] = NA
avohilmo_mod = rq(avohilmo_cost_per_person_day ~ birth_year + sex, data = avohilmo_pred)
avohilmo_pred$predicted_avohilmo_cost_per_person_day = predict(object = avohilmo_mod, newdata = avohilmo_pred)
avohilmo_pred$predicted_avohilmo_cost_per_person_day[which(avohilmo_pred$predicted_avohilmo_cost_per_person_day < 0)] = 0
avohilmo_pred$predicted_avohilmo_cost = avohilmo_pred$predicted_avohilmo_cost * avohilmo_pred$person_day_from_1998_to_2011
avohilmo_pred$total_avohilmo_cost = avohilmo_pred$avohilmo_cost + avohilmo_pred$predicted_avohilmo_cost
avohilmo_cost = avohilmo_pred[,c("finngenid", "total_avohilmo_cost")]
colnames(avohilmo_cost) = c("finngenid", "avohilmo_cost")

# COMBINE DRUG, HILMO, AND AVOHILMO INTO LONG DATAFRAME - 20220211
drug_temp = drug_new[which(drug_new$approx_event_year >= 1998),]
drug_temp$type = "kela"
hilmo_temp = hilmo_new[which(hilmo_new$approx_event_year >= 1998),]
hilmo_temp$type = "hilmo"
avohilmo_temp = avohilmo_new[which(avohilmo_new$approx_event_year >= 2011),]
avohilmo_temp$type = "avohilmo"
long_new = rbindlist(list(setnames(drug_temp[,c("finngenid", "event_age", "approx_event_year", "drug_cost", "type")], c("finngenid", "event_age", "approx_event_year", "cost", "type")),
	hilmo_temp[,c("finngenid", "event_age", "approx_event_year", "cost", "type")],
	avohilmo_temp[,c("finngenid", "event_age", "approx_event_year", "cost", "type")]))
fwrite(long_new, "/finngen/red/jiwoo/longitudinal_data_for_ag.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)

# TEST AVOHILMO COSTS - 20220211
par(mfrow=c(1,5))
hist(avohilmo_temp$predicted_avohilmo_cost / avohilmo_temp$person_day_from_1998_to_2011 * 365.25, main = "1998-2011", xlab = "Cost per Person-Year", breaks = 50)
hist(avohilmo_temp$avohilmo_cost / avohilmo_temp$person_day_since_2011 * 365.25, main = "2011-2021", xlab = "Cost per Person-Year", breaks = 1000)
hist(avohilmo_temp$total_avohilmo_cost / avohilmo_temp$person_day_since_1998 * 365.25, main = "1998-2021", xlab = "Cost per Person-Year", breaks = 1000)
hist(avohilmo_temp$avohilmo_cost / avohilmo_temp$person_day_since_2011 * 365.25, main = "2011-2021 ZOOMED IN", xlab = "Cost per Person-Year", xlim = c(0, 1000), breaks = 1000)
hist(avohilmo_temp$total_avohilmo_cost / avohilmo_temp$person_day_since_1998 * 365.25, main = "1998-2021 ZOOMED IN", xlab = "Cost per Person-Year", xlim = c(0, 1000), breaks = 100)

summary(avohilmo_temp$predicted_avohilmo_cost / avohilmo_temp$person_day_from_1998_to_2011 * 365.25)
summary(avohilmo_temp$avohilmo_cost / avohilmo_temp$person_day_since_2011 * 365.25)
summary(avohilmo_temp$total_avohilmo_cost / avohilmo_temp$person_day_since_1998 * 365.25)

# EXCLUDE PRIVATE PRACTICE INDIVIDUALS - 20220208
id_private = total_new$finngenid[which((total_new$avohilmo_cost == 0 | is.na(total_new$avohilmo_cost)) & (total_new$hilmo_cost != 0 | !is.na(total_new$hilmo_cost) | total_new$drug_cost != 0 | !is.na(total_new$drug_cost)))]
total_new = total_new[!which(total_new$finngenid %in% id_private),]
dim(total_new)
# 236101

# ESTIMATION OF LIFETIME COSTS - 20220207

id1 = 1:100
sex = sample(c(0, 1), 100, replace = TRUE)
drug_percentile = sample(0:100, 100, replace = TRUE)/100
hilmo_percentile = sample(0:100, 100, replace = TRUE)/100
avohilmo_percentile = sample(0:100, 100, replace = TRUE)/100

id2 = sample(1:100, 1000, replace = TRUE)
age = sample(seq(0, 100, 5), 1000, replace = TRUE) # here we are assuming that there is no 10 year censoring
drug_cost = rnorm(1000, 100, 10)
hilmo_cost = rnorm(1000, 100, 10)
avohilmo_cost = rnorm(1000, 5000, 1000)

df1 = data.frame("id" = id1, "sex" = sex, "drug_percentile" = drug_percentile, "hilmo_percentile" = hilmo_percentile, "avohilmo_percentile" = avohilmo_percentile)
df2 = data.frame("id" = id2, "age" = age, "drug_cost" = drug_cost, "hilmo_cost" = hilmo_cost, "avohilmo_cost" = avohilmo_cost)

df = merge(df1, df2, by = "id")

res = NULL
age_interval = seq(0, 100, 5)
for (i in 1:length(age_interval)) {
	temp = df[which(df$age >= age_interval[i] & df$age < age_interval[i] + 5),]
	drug_mod = lm(drug_cost ~ sex + drug_percentile + hilmo_percentile + avohilmo_percentile, data = temp)
	hilmo_mod = lm(hilmo_cost ~ sex + drug_percentile + hilmo_percentile + avohilmo_percentile, data = temp)
	avohilmo_mod = lm(avohilmo_cost ~ sex + drug_percentile + hilmo_percentile + avohilmo_percentile, data = temp)
	drug_pred = predict(drug_mod, newdata = df1)
	hilmo_pred = predict(hilmo_mod, newdata = df1)
	avohilmo_pred = predict(avohilmo_mod, newdata = df1)
	cost_pred = drug_pred + hilmo_pred + avohilmo_pred
	res = cbind(res, cost_pred)
}
res = as.data.frame(res)
colnames(res) = paste0("age", age_interval)
res$lifetime_cost = rowSums(res)

df$lifetime_cost = res$lifetime_cost

##################################################
########## EXTRACT VARIANTS - 20220207
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

sbp_variants_all = fread(file = paste0(path, "gwas/sbp_variants_all.txt"), header = FALSE, stringsAsFactors = FALSE)
sbp_variants_sig = fread(file = paste0(path, "gwas/sbp_variants_sig.txt"), header = FALSE, stringsAsFactors = FALSE)
colnames(sbp_variants_all) = c("rsid")
colnames(sbp_variants_sig) = c("rsid")

#rsid = fread(file = paste0(path, "MR/finngen_rsids.tsv"), sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#colnames(rsid) = c("rsid","chr_pos_hg38", "SNP")
rsid_cpra = fread(file = paste0(path, "gwas/finngen_rsids_cpra.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)

dim(sbp_variants_all) # 145
sbp_variants_all = merge(sbp_variants_all, rsid_cpra, by = "rsid", all.x = TRUE)
dim(sbp_variants_all) # 148

dim(sbp_variants_sig) # 37
sbp_variants_sig = merge(sbp_variants_sig, rsid_cpra, by = "rsid", all.x = TRUE)
dim(sbp_variants_sig) # 37

sbp_variants_all$allele = substr(sbp_variants_all$SNP, nchar(sbp_variants_all$SNP), nchar(sbp_variants_all$SNP))
fwrite(list(sbp_variants_all$SNP), paste0(path, "gwas/sbp_variants_all.csv"), col.names = FALSE)
fwrite(sbp_variants_all[,c("SNP", "allele")], paste0(path, "gwas/sbp_variants_all_alt.csv"), col.names = FALSE, quote = FALSE, sep = "\t")

sbp_variants_sig$allele = substr(sbp_variants_sig$SNP, nchar(sbp_variants_sig$SNP), nchar(sbp_variants_sig$SNP))
fwrite(list(sbp_variants_sig$SNP), paste0(path, "gwas/sbp_variants_sig.csv"), col.names = FALSE)
fwrite(sbp_variants_sig[,c("SNP", "allele")], paste0(path, "gwas/sbp_variants_sig_alt.csv"), col.names = FALSE, quote = FALSE, sep = "\t")

plink --bfile /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8 \
--extract sbp_variants_all.csv \
--recode-allele sbp_variants_all_alt.csv \
--recode A --threads 16 --out sbp_variants_all

plink --bfile /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8 \
--extract sbp_variants_sig.csv \
--recode-allele sbp_variants_sig_alt.csv \
--recode A --threads 16 --out sbp_variants_sig

# DO MENDELIAN RANDOMIZATION - 20220118

#res_single = mr_singlesnp(dat)
#fwrite(list(res_single$SNP[-which(grepl("all", res_single$SNP, ignore.case = TRUE))]), paste0(path, "sbp_variants_all.txt"), col.names = FALSE, row.names = FALSE)
#fwrite(list(res_single$SNP[-which(grepl("all", res_single$SNP, ignore.case = TRUE) | res_single$p >= 0.05)]), paste0(path, "sbp_variants_sig.txt"), col.names = FALSE, row.names = FALSE)

# MODEL DRUG COST ~ PRS BY AGE INTERVAL - 20211222
age_low = seq(0, 100, 20)
age_high = seq(20, 120, 20)
mod_res = NULL
for (j in 1:length(age_low)) {
	print(paste0(age_low[j], " to ", age_high[j]))
	for (i in 1:length(prs_list)) {
		print(prs_list[i])		
		temp = total_new[which(!is.na(total_new[[prs_list[i]]]) & total_new$fu_end_age < age_high[j] & total_new$fu_end_age >= age_low[j]),]
		temp$total_cost = log(temp$total_cost + 1)
		temp$prs = scale(temp[[prs_list[i]]])
		mod = lm(formula = total_cost ~ prs + bl_age + length_of_fu + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = temp)
		coef = coef(summary(mod))["prs", "Estimate"]
		if (coef < 0) {
			temp$prs = -scale(temp[[prs_list[i]]])
			mod = lm(formula = total_cost ~ prs + bl_age + length_of_fu + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = temp)
			coef = coef(summary(mod))["prs", "Estimate"]
		}
		ste = coef(summary(mod))["prs", "Std. Error"]
		pval = coef(summary(mod))["prs", "Pr(>|t|)"]
		mod_res = rbind(mod_res, c(age_low[j], age_high[j], nrow(temp), prs_list[i], coef, ste, pval))
	}
}
mod_res = as.data.frame(mod_res)
colnames(mod_res) = c("low", "high", "n", "prs", "coef", "ste", "pval")
mod_res$low = as.numeric(mod_res$low)
mod_res$high = as.numeric(mod_res$high)
mod_res$n = as.numeric(mod_res$n)
mod_res$coef = as.numeric(mod_res$coef)
mod_res$ste = as.numeric(mod_res$ste)
mod_res$pval = as.numeric(mod_res$pval)
mod_res$beta = (exp(mod_res$coef)-1)*100
mod_res$beta_lower = (exp(mod_res$coef - 1.96 * mod_res$ste)-1)*100
mod_res$beta_upper = (exp(mod_res$coef + 1.96 * mod_res$ste)-1)*100
mod_res = merge(prs_metadata, mod_res, by.x = "basename", by.y = "prs", all = TRUE)

# PLOT MODEL RESULTS
mod_temp = mod_res
mod_temp = mod_res[which(mod_res$basename == pheno),]
mod_temp = mod_res[which(mod_res$basename %in% pheno_list & mod_res$include == 1),]
mod_temp = mod_res[which(mod_res$category == "Behavior" & mod_res$include == 1),]

ggplot() +
	geom_line(data = mod_temp, mapping = aes(x = low, y = beta, color = basename), size = 1) +
	geom_point(data = mod_temp, mapping = aes(x = low, y = beta, shape = ifelse(pval < 0.05/nrow(mod_res), "Yes", "No")), size = 3) +
	geom_errorbar(data = mod_temp, mapping = aes(x = low, ymin = beta_lower, ymax = beta_upper, color = basename), width = 0.1, size = 1) +
	scale_x_continuous(labels = age_low, breaks = age_low) +
	scale_shape_manual(values = c(1, 16)) +
	scale_color_discrete(limits = mod_temp$basename, labels = mod_temp$phenotype_label) + 
  labs(x = "Age Interval", y = "Percent Change in Healthcare Cost", color = "Phenotype", shape = "Significant?") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_line(data = mod_res, mapping = aes(x = low, y = n), color = "black", size = 1) +
	labs(x = "Age Interval", y = "Number of Individuals per Age Interval") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# DO QUALITY CONTROL FOR PEOPLE WHO DO NOT HAVE AN EVENT WITHIN FIVE YEARS OF FOLLOW UP - 20211209
#long_new = rbind(drug_new[])
#long_new = long_new %>% group_by(finngenid) %>% summarise(last_event_age = max(event_age))
#long_new = as.data.frame(long_new)
#long_new = merge(endpoint_all, long_new, by = "finngenid", all.x = TRUE)
#long_new$length_of_last = long_new$fu_end_age - long_new$last_event_age

# DO MENDELIAN RANDOMIZATION - 20211203

#ao_df[(grepl("stroke", ao_df$trait, ignore.case = TRUE) & grepl("neale", ao_df$author, ignore.case = TRUE)), c("id", "trait", "author", "sex", "population", "unit", "category")] # ukb-d-C_STROKE

#ao_df[(grepl("coronary heart disease", ao_df$trait, ignore.case = TRUE)), c("id", "trait", "author", "sex", "population", "unit", "category")] # ieu-a-7

#ao_df[(grepl("atrial fibrillation", ao_df$trait, ignore.case = TRUE)), c("id", "trait", "author", "sex", "population", "unit", "category")] # ukb-a-536

#ao_df[(grepl("type 2 diabetes", ao_df$trait, ignore.case = TRUE) & grepl("mahajan", ao_df$author, ignore.case = TRUE)), c("id", "trait", "author", "sex", "population", "unit", "category")] # ebi-a-GCST007517

#ao_df[(grepl("depression", ao_df$trait, ignore.case = TRUE) & grepl("neale", ao_df$author, ignore.case = TRUE)), c("id", "trait", "author", "sex", "population", "unit", "category")] # ukb-d-F5_DEPRESSIO

##################################################
########## DOWNLOAD FROM SANDBOX - 20211202
##################################################

gcloud auth login
gsutil cp gs://... .

##################################################
########## GENERATE PRS - 20211202
##################################################

plink --bfile mydata --score myprofile.raw

plink \
    --bfile EUR.QC \
    --score Height.QC.Transformed 3 4 12 header \
    --q-score-range range_list SNP.pvalue \
    --extract EUR.valid.snp \
    --out EUR

cd ~/jiwoo/gwas_cost_repository

cp /finngen/green/P15151_model.txt .

plink \
	--bed /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8.bed \
	--bim /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8.bim \
	--fam /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8.fam \
	--score P15151_model.txt 1 5 6 header 

#need cov (EUR.cov)
#need post QCed sum stat (Height.QC.gz)
#need phenotype of the samples (EUR.height)

# MENDELIAN RANDOMIZATION - 20211202

# Get instruments
exposure_dat <- extract_instruments(outcomes = "ukb-d-30720_raw")
exposure_dat_test <- extract_instruments(outcomes = "ieu-a-2")

# Get effects of instruments on outcome
cost = fread(file = paste0(path, "MR/total_cost.pheweb"), header = TRUE, stringsAsFactors = FALSE)
cost_log = fread(file = paste0(path, "MR/log_total_cost.pheweb"), header = TRUE, stringsAsFactors = FALSE)
colnames(cost) = c("chr", "pos", "ref", "alt", "pval", "beta", "se", "maf", "n_hom_cases", "n_het_cases", "n_hom_controls", "n_het_controls")
colnames(cost_log) = c("chr", "pos", "ref", "alt", "pval", "beta", "se", "maf", "n_hom_cases", "n_het_cases", "n_hom_controls", "n_het_controls")
rsid_cpra = fread(file = paste0(path, "MR/finngen_rsids_cpra.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rsid_cpra = rsid_cpra[which(rsid_cpra$rsid %in% exposure_dat_test$SNP),]
rsid_cpra$chr = vapply(strsplit(rsid_cpra$chr_pos_hg38, "_"), `[`, 1, FUN.VALUE=character(1))
rsid_cpra$pos = as.numeric(vapply(strsplit(rsid_cpra$chr_pos_hg38, "_"), `[`, 2, FUN.VALUE=character(1)))
rsid_cpra = rsid_cpra[!duplicated(rsid_cpra$rsid),]
cost = merge(rsid_cpra, cost, by = c("chr", "pos"))
cost_log = merge(rsid_cpra, cost_log, by = c("chr", "pos"))
fwrite(cost, paste0(path, "MR/total_cost_test.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
fwrite(cost_log, paste0(path, "MR/log_total_cost_test.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
# 196 SNPs
# 178 SNPs
# 174 SNPs

# 79 SNPs 
# 77 SNPs 

outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,
	filename = paste0(path, "MR/total_cost.csv"),
	sep = ",",
	snp_col = "rsid",
	beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  pval_col = "pval"
)
outcome_log_dat <- read_outcome_data(snps = exposure_dat$SNP,
	filename = paste0(path, "MR/log_total_cost.csv"),
	sep = ",",
	snp_col = "rsid",
	beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  pval_col = "pval"
)
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)
dat_log <- harmonise_data(exposure_dat, outcome_log_dat)
# Perform MR
res <- mr(dat) 
res_log <- mr(dat_log)
# Plot MR
mr_scatter_plot(res, dat) 
mr_scatter_plot(res_log, dat_log) 

outcome_dat_test <- read_outcome_data(snps = exposure_dat_test$SNP,
	filename = paste0(path, "MR/total_cost_test.csv"),
	sep = ",",
	snp_col = "rsid",
	beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  pval_col = "pval"
)
outcome_log_dat_test <- read_outcome_data(snps = exposure_dat_test$SNP,
	filename = paste0(path, "MR/log_total_cost_test.csv"),
	sep = ",",
	snp_col = "rsid",
	beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  pval_col = "pval"
)

# Harmonise the exposure and outcome data
dat_test <- harmonise_data(exposure_dat_test, outcome_dat_test)
dat_log_test <- harmonise_data(exposure_dat_test, outcome_log_dat_test)
# Perform MR
res_test <- mr(dat_test) 
res_log_test <- mr(dat_log_test)
# Plot MR
mr_scatter_plot(res_test, dat_test) 
mr_scatter_plot(res_log_test, dat_log_test) 

ao <- available_outcomes()
ao_df <- as.data.frame(ao)

# LOAD PRS DATA - 20211130

#prs_file_list = list.files(path = paste0(path, "prs/"), pattern = "\\.sscore$")
#prs_list = str_replace_all(str_remove(str_remove(prs_file_list, fixed("finngen_R8_")), fixed(".sscore")), "[-+]", "_")

##################################################
########## DO ASSOCIATION ANALYSIS - 20211130
##################################################

ggplot() +
	geom_point(data = mod_temp, mapping = aes(x = reorder(basename, -abs_coef), y = abs_coef, color = ifelse(pval < 0.05/nrow(mod_temp), "Yes", "No")), size = 1) +
	coord_flip () +
	geom_errorbar(data = mod_temp, mapping = aes(x = reorder(basename, -abs_coef), ymin = abs(confint_lower), ymax = abs(confint_upper), color = ifelse(pval < 0.05/nrow(mod_temp), "Yes", "No")), width = 0.1, size = 1) +
	scale_x_discrete(limits = mod_temp$basename[order(-mod_temp$abs_coef)], labels = mod_temp$phenotype_label[order(-mod_temp$abs_coef)]) + 
	labs(x = "", y = "Estimate", color = "Significant?") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

##################################################
########## CLEAN UK BIOBANK SUMMARY STATISTICS - 20211129
##################################################

cp /finngen/green/jiwoo/... .
gunzip ...

pval_threshold = 5 * 10 ^ -8
pheno = "biomarkers-30710-both_sexes-irnt"
website = "https://www.ncbi.nlm.nih.gov/gene/1401"
chr_num = 1
start_pos = 159712289
end_pos = 159714589

ukbiobank_sumstats = fread(paste0(path, "mr/", pheno, "_ukbiobank_sumstats.tsv"), header = TRUE)
ukbiobank_sumstats$chr[which(ukbiobank_sumstats$chr == "X")] = 23
ukbiobank_sumstats$chr = as.numeric(ukbiobank_sumstats$chr)
dim(ukbiobank_sumstats) 
# 28987534
ukbiobank_sumstats_new = ukbiobank_sumstats[which(ukbiobank_sumstats$chr == chr_num & ukbiobank_sumstats$pos >= start_pos & ukbiobank_sumstats$pos <= end_pos),]
#ukbiobank_sumstats_new = ukbiobank_sumstats[which((ukbiobank_sumstats$chr == chr_num & ukbiobank_sumstats$pos >= start_pos & ukbiobank_sumstats$pos <= end_pos) | (ukbiobank_sumstats$pval_meta < pval_threshold)),]
dim(ukbiobank_sumstats_new)
# 38
ukbiobank_sumstats_new$chr_pos_ref_alt = paste0("chr", ukbiobank_sumstats_new$chr, "_", ukbiobank_sumstats_new$pos, "_", ukbiobank_sumstats_new$ref, "_", ukbiobank_sumstats_new$alt)
write.table(ukbiobank_sumstats_new$chr_pos_ref_alt, paste0(path, "mr/", pheno, "_ukbiobank_snps_hg37.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

##################################################
########## CALCULATES FINNGEN BETAS - 20211129
##################################################

plink --bfile /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8 --extract biomarkers-30710-both_sexes-irnt_ukbiobank_snps_hg38.txt --recodeA --out biomarkers-30710-both_sexes-irnt_finngen_matrix

total_mini = total_new[,c("finngenid", "sex", "birth_year", "death_year", "death_age", "death", "bl_age", "bl_year", "fu_start_age", "fu_end_age", "length_of_fu", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "total_cost")]
pheno = "biomarkers-30710-both_sexes-irnt"
geno_mat = fread(paste0(path, "mr/", pheno, "_finngen_matrix.raw"), header = TRUE, stringsAsFactors = FALSE)
geno_mat = geno_mat[,-c("FID", "PAT", "MAT", "SEX", "PHENOTYPE")]
colnames(geno_mat)[1] = "finngenid"

geno_mat_new = merge(total_mini, geno_mat, by = "finngenid", all.x = TRUE)
snp_list = colnames(geno_mat_new)[grep("chr", colnames(geno_mat_new))]

finngen_sumstats = NULL
for (i in 1:length(snp_list)) {
	print(snp_list[i])
	geno_mat_temp = geno_mat_new
	geno_mat_temp$age_sex = geno_mat_temp$bl_age * ifelse(geno_mat_temp$sex == "female", 1, 0)
	geno_mat_temp$age_squared = geno_mat_temp$bl_age ^ 2
	geno_mat_temp$age_squared_sex = geno_mat_temp$age_squared * ifelse(geno_mat_temp$sex == "female", 1, 0)
	geno_mat_temp$log_total_cost = log(geno_mat_temp$total_cost + 1)
	form_mod = as.formula(sprintf("total_cost ~ %s + bl_age + sex + age_sex + bl_age^2 + age_squared_sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", snp_list[i]))
	form_mod_log = as.formula(sprintf("log_total_cost ~ %s + bl_age + sex + age_sex + age_squared + age_squared_sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", snp_list[i]))
	glm_mod = lm(formula = form_mod, data = geno_mat_temp)
	glm_mod_log = lm(formula = form_mod_log, data = geno_mat_temp)
	coef = coef(summary(glm_mod))[snp_list[i], "Estimate"]
	std = coef(summary(glm_mod))[snp_list[i], "Std. Error"]
	pval = coef(summary(glm_mod))[snp_list[i], "Pr(>|t|)"]
	coef_log = coef(summary(glm_mod_log))[snp_list[i], "Estimate"]
	std_log = coef(summary(glm_mod_log))[snp_list[i], "Std. Error"]
	pval_log = coef(summary(glm_mod_log))[snp_list[i], "Pr(>|t|)"]
	finngen_sumstats = rbind(finngen_sumstats, c(snp_list[i], coef, std, pval, coef_log, std_log, pval_log))
}
finngen_sumstats = as.data.frame(finngen_sumstats)
colnames(finngen_sumstats) = c("snp", "coef_finngen", "std_finngen", "pval_finngen", "coef_log_finngen", "std_log_finngen", "pval_log_finngen")
finngen_sumstats$coef_finngen = as.numeric(finngen_sumstats$coef_finngen)
finngen_sumstats$std_finngen = as.numeric(finngen_sumstats$std_finngen)
finngen_sumstats$pval_finngen = as.numeric(finngen_sumstats$pval_finngen)
finngen_sumstats$coef_log_finngen = as.numeric(finngen_sumstats$coef_log_finngen)
finngen_sumstats$std_log_finngen = as.numeric(finngen_sumstats$std_log_finngen)
finngen_sumstats$pval_log_finngen = as.numeric(finngen_sumstats$pval_log_finngen)
finngen_sumstats$chr_pos_ref_alt = substring(finngen_sumstats$snp, 1, nchar(finngen_sumstats$snp) - 6)

hg37 = fread(paste0(path, "mr/", pheno, "_ukbiobank_snps_hg37.txt"), header = FALSE, stringsAsFactors = FALSE)
hg38 = fread(paste0(path, "mr/", pheno, "_ukbiobank_snps_hg38.txt"), header = FALSE, stringsAsFactors = FALSE)
hg37_to_hg38 = cbind(hg37, hg38)
colnames(hg37_to_hg38) = c("hg37", "hg38")

ukbiobank_sumstats_mini = merge(ukbiobank_sumstats_new, hg37_to_hg38, by.x = "chr_pos_ref_alt", by.y = "hg37", all = TRUE)
ukbiobank_sumstats_mini$chr_pos_ref_alt = substring(ukbiobank_sumstats_mini$hg38, 1, nchar(ukbiobank_sumstats_mini$hg38) - 4)
ukbiobank_sumstats_mini = ukbiobank_sumstats_mini[,c("chr_pos_ref_alt", "beta_EUR", "se_EUR", "pval_EUR")]
colnames(ukbiobank_sumstats_mini) = c("chr_pos_ref_alt", "beta_ukbiobank", "std_ukbiobank", "pval_ukbiobank")
sumstats_all = merge(finngen_sumstats, ukbiobank_sumstats_mini, by = "chr_pos_ref_alt", all.x = TRUE)

ggplot() +
	geom_point(data = sumstats_all, mapping = aes(x = coef_log_finngen, y = beta_ukbiobank), size = 3) +
	#scale_x_discrete(limits = seq(1, 10, 1), labels = seq(1, 10, 1)) + 
	labs(x = "Beta for log(Healthcare Cost)", y = "Beta for C-Reactive Protein") + 
	#facet_wrap(~prs) + 
	#facet_wrap(~prs, scales = "free") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

##################################################
########## MAKE SNP LIST FOR FINNGEN FROM UK BIOBANK - 20211125
##################################################

if(!require(data.table)) {install.packages("data.table"); library(data.table)}

path = "/home/ivm/jiwoo/healthcare_cost_repository/"
pheno = "biomarkers-30710-both_sexes-irnt"
full_snp_list = fread(paste0(path, "mr/full_snp_list.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(full_snp_list) = c("chr", "pos", "rsid")
full_snp_list$chr_pos = paste0(substring(full_snp_list$chr, 4), ":", full_snp_list$pos)
full_snp_list_new = full_snp_list[,c("chr_pos", "rsid")]
write.table(full_snp_list_new, paste0(path, "mr/full_snp_list_new.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

total = fread(paste0(path, "mr/", pheno, "_ukbiobank_sumstats.tsv"), header = TRUE, stringsAsFactors = FALSE)
my_snp_list = total[,c("chr", "pos", "ref", "alt")]
my_snp_list$chr_pos = paste0("chr", my_snp_list$chr, "_", my_snp_list$pos, "_", my_snp_list$ref, "_", my_snp_list$alt)
write.table(my_snp_list$chr_pos, paste0(path, "mr/", pheno, "_ukbiobank_snp_list.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

#plink --update-map full_snp_list_new.txt biomarkers-30710-both_sexes-irnt_ukbiobank_snp_list.bim --make-bed --out out

#my_snp_list = total_new[,c("chr", "pos", "ref", "alt")]
#my_snp_list$end = my_snp_list$pos
#my_snp_list$allele = paste0(my_snp_list$ref, "/", my_snp_list$alt)
#my_snp_list$strand = "+"
#my_snp_list = my_snp_list[,c("chr", "pos", "end", "allele", "strand")]
#colnames(my_snp_list) = c("chrom", "start", "end", "allele", "strand")
#write.table(my_snp_list, paste0(path, "my_snp_list.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

##################################################
########## MODEL LIFETIME COST - 20211119
##################################################

df <- data.frame(age_group=1:10)
df <- df %>%
  mutate(five_year_costs=1000 + age_group*100 + rnorm(n=10, mean=0, sd=100),
         p_death = seq(0.1, 1, 0.1))


#simulate costs for 1 person
simulate_costs <- function(data) {
  #everyone starts with first age group costs
  costs <- data$five_year_costs[1]
  for (i in 1:nrow(data)) {
    #if person survives to next age group, add costs from that age group
    if (runif(1) >  data$p_death[i] & i+1<nrow(data)) {
      costs <- costs + data$five_year_costs[i+1]
    } else { #else the person dies, so return the costs accumulated so far
      return(costs)
    }
  }
  return(costs)
}

simulate_costs(df)

#run the simulation 10000 times for 1 person
sim_costs <- replicate(n=10000, expr=simulate_costs(df))
summary(sim_costs) #mean approximates the expected lifetime costs
hist(sim_costs)

#calculate probabilities analytically
df$prob_dying_at_agegroup <- NA
for (i in 1:nrow(df)) {
  
  if (i==1) {
    #probability of dying at first age group
    df$prob_dying_at_agegroup[i] <- df$p_death[i]
  } else if (i>1) {
    #probability of dying at i'th age group, conditional on having survived to that point
      df$prob_dying_at_agegroup[i] <- prod(1-df$p_death[1:(i-1)]) * df$p_death[i]
  }
   
}

#calculate cumulative costs
df$cumcost <- cumsum(df$five_year_costs)

#expected costs
df$expected_cost_from_dying <- df$cumcost * df$prob_dying
sum(df$expected_cost_from_dying) #this should be same as the mean of the simulated costs

# MODEL DRUG COST ~ PRS BY AGE INTERVAL

ggplot() +
	geom_point(data = mod_res[which(mod_res$prs == "Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted"),], mapping = aes(x = low, y = coef, color = ifelse(pval < 0.05/nrow(mod_res), "Yes", "No")), size = 3) +
	geom_errorbar(data = mod_res[which(mod_res$prs == "Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted"),], mapping = aes(x = low, ymin = coef - std, ymax = coef + std, color = ifelse(pval < 0.05/nrow(mod_res), "Yes", "No")), width = 0.1, size = 1) +
	scale_x_continuous(labels = age_low, breaks = age_low) +
	labs(title = "Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted", x = "Age Interval", y = "Estimate", color = "Significant?") + 
	theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))


ggplot() +
	geom_point(data = mod_res[which(mod_res$prs == "PGC_UKB_depression_genome_wide"),], mapping = aes(x = low, y = coef, color = ifelse(pval < 0.05/nrow(mod_res), "Yes", "No")), size = 3) +
	geom_errorbar(data = mod_res[which(mod_res$prs == "PGC_UKB_depression_genome_wide"),], mapping = aes(x = low, ymin = coef - std, ymax = coef + std, color = ifelse(pval < 0.05/nrow(mod_res), "Yes", "No")), width = 0.1, size = 1) +
	scale_x_continuous(labels = age_low, breaks = age_low) +
	labs(title = "PGC_UKB_depression_genome_wide", x = "Age Interval", y = "Estimate", color = "Significant?") + 
	theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))


ggplot() +
	geom_point(data = mod_res[which(mod_res$prs == "UKB_ICBPmeta750k_SBPsummaryResults"),], mapping = aes(x = low, y = coef, color = ifelse(pval < 0.05/nrow(mod_res), "Yes", "No")), size = 3) +
	geom_errorbar(data = mod_res[which(mod_res$prs == "UKB_ICBPmeta750k_SBPsummaryResults"),], mapping = aes(x = low, ymin = coef - std, ymax = coef + std, color = ifelse(pval < 0.05/nrow(mod_res), "Yes", "No")), width = 0.1, size = 1) +
	scale_x_continuous(labels = age_low, breaks = age_low) +
	labs(title = "UKB_ICBPmeta750k_SBPsummaryResults", x = "Age Interval", y = "Estimate", color = "Significant?") + 
	theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))


# PLOT DISTRIBUTION OF DRUG COST BY PRS PERCENTILE
temp = total_new %>% mutate(prs = ntile(!!as.name(trait), 10)) %>% group_by(prs, sex) %>% summarise(cost = mean(total_cost, na.rm = TRUE), 
	confint_lower = t.test(total_cost)$conf.int[1],
	confint_upper = t.test(total_cost)$conf.int[2])

# PLOT DISTRIBUTION OF DRUG COST BY PRS PERCENTILE - 20211108
# cad.add.160614.website
# Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted
# PGC_UKB_depression_genome_wide
# UKB_ICBPmeta750k_SBPsummaryResults

plot(total$UKB_ICBPmeta750k_SBPsummaryResults, total$DRUG_COST, xlab = "SBP PRS", ylab = "Drug Cost")
plot(total$UKB_ICBPmeta750k_SBPsummaryResults, total$LOG_DRUG_COST), xlab = "SBP PRS", ylab = "log10(Drug Cost)")

temp = total %>% mutate(prs = ntile(Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted, 20)) %>% group_by(prs) %>% summarise(cost = median(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted", ylab = "Median Cost")
temp = total %>% mutate(prs = ntile(PGC_UKB_depression_genome_wide, 20)) %>% group_by(prs) %>% summarise(cost = median(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "PGC_UKB_depression_genome_wide", ylab = "Median Cost")
temp = total %>% mutate(prs = ntile(UKB_ICBPmeta750k_SBPsummaryResults, 20)) %>% group_by(prs) %>% summarise(cost = median(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "UKB_ICBPmeta750k_SBPsummaryResults", ylab = "Median Cost")

temp = total %>% mutate(prs = ntile(Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted, 20)) %>% group_by(prs) %>% summarise(cost = median(log10(DRUG_COST), na.rm = TRUE))
plot(temp, xlab = "Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted", ylab = "Median log10(Cost)")
temp = total %>% mutate(prs = ntile(PGC_UKB_depression_genome_wide, 20)) %>% group_by(prs) %>% summarise(cost = median(log10(DRUG_COST), na.rm = TRUE))
plot(temp, xlab = "PGC_UKB_depression_genome_wide", ylab = "Median log10(Cost)")
temp = total %>% mutate(prs = ntile(UKB_ICBPmeta750k_SBPsummaryResults, 20)) %>% group_by(prs) %>% summarise(cost = median(log10(DRUG_COST), na.rm = TRUE))
plot(temp, xlab = "UKB_ICBPmeta750k_SBPsummaryResults", ylab = "Median log10(Cost)")

temp = total %>% mutate(prs = ntile(cad.add.160614.website, 100)) %>% group_by(prs) %>% summarise(cost = mean(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "cad.add.160614.website", ylab = "Cost")
temp = total %>% mutate(prs = ntile(Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted, 100)) %>% group_by(prs) %>% summarise(cost = mean(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted", ylab = "Cost")
temp = total %>% mutate(prs = ntile(PGC_UKB_depression_genome_wide, 100)) %>% group_by(prs) %>% summarise(cost = mean(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "PGC_UKB_depression_genome_wide", ylab = "Cost")
temp = total %>% mutate(prs = ntile(UKB_ICBPmeta750k_SBPsummaryResults, 100)) %>% group_by(prs) %>% summarise(cost = mean(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "UKB_ICBPmeta750k_SBPsummaryResults", ylab = "Cost")

temp = total %>% mutate(prs = ntile(cad.add.160614.website, 20)) %>% group_by(prs) %>% summarise(cost = mean(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "cad.add.160614.website", ylab = "Mean Cost")
temp = total %>% mutate(prs = ntile(Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted, 20)) %>% group_by(prs) %>% summarise(cost = mean(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted", ylab = "Mean Cost")
temp = total %>% mutate(prs = ntile(PGC_UKB_depression_genome_wide, 20)) %>% group_by(prs) %>% summarise(cost = mean(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "PGC_UKB_depression_genome_wide", ylab = "Mean Cost")
temp = total %>% mutate(prs = ntile(UKB_ICBPmeta750k_SBPsummaryResults, 20)) %>% group_by(prs) %>% summarise(cost = mean(DRUG_COST, na.rm = TRUE))
plot(temp, xlab = "UKB_ICBPmeta750k_SBPsummaryResults", ylab = "Mean Cost")

# PLOT DISTRIBUTION OF HEALTHCARE COSTS BY AGE INTERVAL - 20211108
ggplot() +
	geom_line(data = age_res, mapping = aes(x = age_low, y = person_cost), color = "black", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = female_cost), color = "red", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = male_cost), color = "blue", size = 1) +
	labs(title = "Total", x = "Age", y = "Cost") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_line(data = age_res, mapping = aes(x = age_low, y = log_person_cost), color = "black", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = log_female_cost), color = "red", size = 1) +
	geom_line(data = age_res, mapping = aes(x = age_low, y = log_male_cost), color = "blue", size = 1) +
	labs(title = "Total", x = "Age", y = "ln(Cost)") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# MINIMUM DATA - 20211108
cp /finngen/library-red/finngen_R8/phenotype_2.0/data/finngen_R8_minimum.txt.gz .
gunzip finngen_R8_minimum.txt.gz

head -1 finngen_R8_minimum.txt | tr '\t' '\n' | cat -n | grep "FINNGENID"
head -1 finngen_R8_minimum.txt | tr '\t' '\n' | cat -n | grep "movedabroad"

awk '{print $1, $15}' finngen_R8_minimum.txt > finngen_R8_minimum_clean.txt

# FIGURE OUT AVOHILMO SPIKES - 20211104
temp = total_new[which(log(total_new$avohilmo_cost + 1) < 6),]
ggplot() +
	geom_histogram(data = temp, mapping = aes(x = avohilmo_cost, fill = "AvoHILMO"), alpha = 0.5, bins = 100) +
	labs(x = "Cost", y = "Frequency", fill = "") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))
tab = as.data.frame(table(temp$avohilmo_cost))
tab = tab[order(-tab$Freq),]
top = as.numeric(as.character(tab$Var1[1:3]))
id_spike = avohilmo_cost$finngenid[which(avohilmo_cost$avohilmo_cost == top[2])]
df_spike = avohilmo_new_2011[which(avohilmo_new_2011$finngenid %in% id_spike),]
df_spike = df_spike[which(!is.na(df_spike$cost)),]

# EXCLUDE IMMIGRANT INDIVIDUALS - 20211104
min_all = fread(file = paste0(path, "finngen_R8_minimum.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(min_all) = c("finngenid", "movedabroad")
id_immigrant = total_new$finngenid[which(total_new$death == 0 & (as.Date(total_new$bl_year, format = "%Y") + total_new$length_of_fu))]
summary(as.Date(total_new$bl_year, format = "%Y") + years(total_new$length_of_fu))

# PLOT DISTRIBUTION OF DRUG COST - 20211104
ggplot() +
	geom_histogram(data = total, mapping = aes(x = DRUG_COST), color = "black", fill = "gray", bins = 100) +
	labs(x = "Drug Cost", y = "Frequency") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_histogram(data = total, mapping = aes(x = LOG_DRUG_COST), color = "black", fill = "gray", bins = 100) +
	labs(x = "log10(Drug Cost)", y = "Frequency") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_histogram(data = total, mapping = aes(x = LOG_DRUG_COST, fill = SEX), color = "black", alpha = 0.5, bins = 100) +
	labs(x = "log10(Drug Cost)", y = "Frequency", fill = "Sex") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_histogram(data = total, mapping = aes(x = ASINH_DRUG_COST), color = "black", fill = "gray", bins = 100) +
	labs(x = "asinh(Drug Cost)", y = "Frequency") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_histogram(data = total, mapping = aes(x = ASINH_DRUG_COST, fill = SEX), color = "black", alpha = 0.5, bins = 100) +
	labs(x = "asinh(Drug Cost)", y = "Frequency", fill = "Sex") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# REQUEST FOR FILE DOWNLOAD - 20211103
drug_table = as.data.frame(table(drug_new$vnr[is.na(drug_new$drug_cost)]))
drug_table = drug_table[order(-drug_table$Freq),]
drug_table = drug_table[which(drug_table$Freq > 5),]
# 17647 -> 12641
avohilmo_table = as.data.frame(avohilmo_new %>% filter(is.na(avohilmo_cost)) %>% group_by(profession_label, service_type_label, contact_type_label) %>% summarise(n = n()))
avohilmo_table = avohilmo_table[order(-avohilmo_table$n),]
avohilmo_table = avohilmo_table[which(avohilmo_table$n > 5),]
# 96 -> 88
hilmo_table = as.data.frame(hilmo_new %>% filter(is.na(hilmo_cost)) %>% group_by(service_type, specialty_label, hospital_type) %>% summarise(n = n()))
hilmo_table = hilmo_table[order(-hilmo_table$n),]
hilmo_table = hilmo_table[which(hilmo_table$n > 5),]
# 29 -> 25

write.csv(drug_table, "drug_missing_data_summary.csv", row.names = FALSE)
write.csv(avohilmo_table, "avohilmo_missing_data_summary.csv", row.names = FALSE)
write.csv(hilmo_table, "hilmo_missing_data_summary.csv", row.names = FALSE)

temp = avohilmo_new[which(is.na(avohilmo_new$profession_label)),]
temp = as.data.frame(table(temp$profession))
temp = temp[order(-temp$Freq),]
temp = temp[which(temp$Freq > 5),]
# 133 -> 115
write.csv(avohilmo_profession_summary, "avohilmo_profession_summary.csv", row.names = FALSE)

Files:
drug_missing_data_summary.csv
avohilmo_missing_data_summary.csv
hilmo_missing_data_summary.csv

These file contains the number of rows of missing data for each type of event. 
Knowing the number of rows of missing data will be helpful in imputing missing data. 
All events have n > 5. 

File:
avohilmo_profession_summary.csv

This file contains the number of rows for each type of event. 
All events have n > 5. 

The number of rows of missing data were counted, based on column type. 

# LOAD DRUG DATA - 20211029
drug_all = fread(file = paste0(path, "finngen_R8_detailed_longitudinal_purch.txt"), header = TRUE, stringsAsFactors = FALSE)
drug_all = drug_all[,c("FINNGENID", "EVENT_AGE", "APPROX_EVENT_DAY", "CODE3")]
colnames(drug_all)[4] = "VNR"
drug_all$APPROX_EVENT_YEAR = as.integer(format(as.Date(drug_all$APPROX_EVENT_DAY), format = "%Y"))
drug_new = merge(drug_all, vnr_all, by.x = c("VNR", "APPROX_EVENT_YEAR"), by.y = c("VNR", "YEAR"), all.x = TRUE)
dim(drug_new) # 77491626
length(unique(drug_new$FINNGENID)) # 354734

# ENDPOINT LONGITUDINAL DATA - 20211029
cp /finngen/library-red/finngen_R8/phenotype_3.0/data/finngen_R8_endpoint_longitudinal.txt.gz .
gunzip finngen_R8_endpoint_longitudinal.txt.gz

# DETAILED LONGITUDINAL DATA - 20211029
cp /finngen/library-red/finngen_R8/phenotype_3.0/data/finngen_R8_detailed_longitudinal.txt.gz .
gunzip finngen_R8_detailed_longitudinal.txt.gz

head -1 finngen_R8_detailed_longitudinal.txt | tr '\t' '\n' | cat -n | grep "SOURCE"
awk -F, '$2 == "PURCH"' finngen_R8_detailed_longitudinal.txt > finngen_R8_detailed_longitudinal_purch.txt
awk '{if ($2 == "PURCH" || $2 == "SOURCE"){print $0}}' finngen_R8_detailed_longitudinal.txt > finngen_R8_detailed_longitudinal_purch.txt

# CALCULATE DRUG COST - 20210920
fid_list = unique(drug_all$FINNGENID)
cost_res = NULL
for (i in 1:length(fid_list)) {
	cost = sum(drug_new[which(drug_new$FINNGENID == fid_list[i]), "cost"])
	cost_res = rbind(cost_res, c(fid_list[i], cost))
	if (i %% 100 == 0) {
		print(i)
	}
}
cost_res = as.data.frame(cost_res)
colnames(cost_res) = c("FINNGENID", "DRUG_COST")
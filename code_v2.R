##################################################
########## DOWNLOAD DATA
##################################################

cd ~/jiwoo/healthcare_cost_repository

# ENDPOINT DATA
cp /finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_endpoint_1.0.txt.gz .
gunzip finngen_R9_endpoint_1.0.txt.gz

head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "FINNGENID"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "BL_AGE"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "BL_YEAR"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "FU_END_AGE"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "SEX"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "DEATH"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "DEATH_AGE"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "DEATH_YEAR"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "STATIN"
head -1 finngen_R9_endpoint_1.0.txt | tr '\t' '\n' | cat -n | grep "I9_CHD"
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $5654, $5655, $5656, $15614, $15615, $15616}' finngen_R9_endpoint_1.0.txt > finngen_R9_endpoint_1.0_clean.txt

##################################################
########## LOAD DATA
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

# DEFINE FUNCTIONS
convert_camel_case <- function(x) {
  string_vec <- x
  # substitute underscores into spaces
  string_vec <- gsub(" ", "_", string_vec)
  # remove apostrophes
  string_vec <- gsub("\\'", "", string_vec)
  # remove commas
  string_vec <- gsub(",", "", string_vec)
  # remove dashes
  string_vec <- gsub("\\-", "_", string_vec)
  # remove \
  string_vec <- gsub("/", "_", string_vec)
  # convert lowercase
  string_vec <- tolower(string_vec)
  return(string_vec)
}

start_year = 1995
start_year_hilmo_kela = 1998
start_year_avohilmo = 2011
end_year = 2020

# LOAD ENDPOINT DATA
endpoint_all = fread(file = paste0(path, "finngen_R9_endpoint_1.0_clean.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(endpoint_all) = tolower(colnames(endpoint_all))
endpoint_all$birth_year = endpoint_all$bl_year - endpoint_all$bl_age
endpoint_all$fu_end_year = endpoint_all$birth_year + endpoint_all$fu_end_age
endpoint_all$fu_end_year = ifelse(endpoint_all$fu_end_year >= end_year, end_year, endpoint_all$fu_end_year)
endpoint_all$person_day_in_hilmo_kela = (endpoint_all$fu_end_year - start_year_hilmo_kela) * 365.25
endpoint_all$person_day_in_avohilmo = (endpoint_all$fu_end_year - start_year_avohilmo) * 365.25
endpoint_all$person_day_in_avohilmo[which(endpoint_all$person_day_in_avohilmo < 0)] = 0
endpoint_all = endpoint_all[,c("finngenid", "sex", "birth_year", "death_year", "death_age", "death", "bl_age", "bl_year", "fu_end_age", "fu_end_year", "person_day_in_hilmo_kela", "person_day_in_avohilmo", "rx_statin", "rx_statin_age", "rx_statin_year", "i9_chd", "i9_chd_age", "i9_chd_year")]
endpoint_new = endpoint_all[which(endpoint_all$person_day_in_hilmo_kela > 0),]
dim(endpoint_all) # 392423
dim(endpoint_new) # 388465

# LOAD PCA DATA
pca_all = fread(file = "/finngen/library-red/finngen_R9/pca_1.0/data/finngen_R9.eigenvec.txt", header = TRUE, stringsAsFactors = FALSE)
pca_new = pca_all[,-"#FID"]
colnames(pca_new)[which(colnames(pca_new) == "IID")] = "finngenid"
dim(pca_new) # 377498
total = merge(endpoint_new, pca_new, by = "finngenid")
dim(total) # 373460

# LOAD KELA DICTIONARIES
kela_dict = fread(file = paste0(path, "avg_vnro_price.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(kela_dict) = c("year", "vnr", "cost")
kela_dict$cost[which(kela_dict$cost == "Inf")] = NA

# LOAD HILMO DICTIONARIES
hilmo_inpat_dict = fread(file = paste0(path, "hilmo_inpatient_unit_costs_imputed.tsv"))
hilmo_inpat_dict$specialty <- as.character(hilmo_inpat_dict$specialty)
hilmo_outpat_dict = fread(file = paste0(path, "hilmo_outpatient_unit_costs_imputed.tsv"))
hilmo_outpat_dict$specialty <- as.character(hilmo_outpat_dict$specialty)

# LOAD AVOHILMO DICTIONARIES
avohilmo_dict = fread(file = paste0(path, "avohilmo_unit_costs_long.txt"), header = TRUE, stringsAsFactors = FALSE)
avohilmo_dict = avohilmo_dict[,-"contact_type"]
contact_dict = fread(file = paste0(path, "avohilmo_contact_type_labels.txt"), header = TRUE, stringsAsFactors = FALSE)
profession_dict = fread(file = paste0(path, "avohilmo_profession_labels.txt"), header = TRUE, stringsAsFactors = FALSE)
profession_dict$profession = as.character(profession_dict$profession)
service_dict  = fread(file = paste0(path, "avohilmo_service_type_labels.txt"), header = TRUE, stringsAsFactors = FALSE)

# LOAD LONGITUDINAL DATA
long_all = fread(file = paste0(path, "finngen_R9_service_sector_detailed_longitudinal.txt"), header = TRUE, stringsAsFactors = FALSE)
long_all$APPROX_EVENT_YEAR = as.integer(format(as.Date(long_all$APPROX_EVENT_DAY), format = "%Y"))
colnames(long_all) = tolower(colnames(long_all))
long_all = long_all %>% filter(approx_event_year >= start_year_hilmo_kela, approx_event_year < end_year)

##################################################
########## MERGE DATA
##################################################

# MERGE KELA DATA
kela_all = long_all[which(long_all$source == "PURCH"), c("finngenid", "code3", "event_age", "approx_event_year")]
colnames(kela_all)[which(colnames(kela_all) == "code3")] = "vnr"
kela_new = merge(kela_all, kela_dict, by.x = c("vnr", "approx_event_year"), by.y = c("vnr", "year"), all.x = TRUE)

# MERGE HILMO DATA
hilmo_all = long_all %>% 
	filter(source %in% c("INPAT", "OUTPAT")) %>%
	select(finngenid, source, code4, code5, code6, code7, event_age, approx_event_year, index) %>%
	distinct()
colnames(hilmo_all) = c("finngenid", "source", "length_of_stay", "service_type", "specialty", "hospital_type", "event_age", "approx_event_year", "index")
hilmo_all$specialty = substr(hilmo_all$specialty, 0, 2)
hilmo_all = hilmo_all %>% mutate(service_type = case_when(source == "OUTPAT" & service_type == "91" ~ "emergency",
												source == "OUTPAT" ~ "outpatient",
												source == "INPAT" ~ "ward"))
hilmo_all = hilmo_all %>% mutate(hospital_type = case_when(hospital_type == "Central Hospital" ~ "central",
												hospital_type == "Other Hospital" ~ "other",
												hospital_type == "University Hospital" ~ "university"))

# ESTIMATE INPATIENT HILMO COSTS
hilmo_inpat = hilmo_all %>% filter(source == "INPAT")
hilmo_inpat = left_join(hilmo_inpat, as.data.frame(hilmo_inpat_dict)[,c("service_type", "specialty", "hospital_type", "specialty_label", "costs_per_day", "fixed_costs")], by = c("service_type", "specialty", "hospital_type"))
hilmo_inpat = hilmo_inpat %>% filter(!is.na(specialty_label) & !is.na(length_of_stay))
hilmo_inpat$length_of_stay = ifelse(hilmo_inpat$length_of_stay == 0, 1, hilmo_inpat$length_of_stay)
hilmo_inpat$costs_per_day = ifelse(hilmo_inpat$length_of_stay > 100, 250, hilmo_inpat$costs_per_day)
hilmo_inpat = hilmo_inpat %>%
	mutate(cost = fixed_costs + length_of_stay * costs_per_day) %>%
	select(finngenid, source, service_type, specialty, hospital_type, specialty_label, event_age, approx_event_year, cost)

# ESTIMATE OUTPATIENT HILMO COSTS
hilmo_outpat = hilmo_all %>% filter(source == "OUTPAT")
hilmo_outpat$specialty = ifelse(hilmo_outpat$specialty %in% hilmo_outpat_dict$specialty, hilmo_outpat$specialty, 10)
hilmo_outpat = left_join(hilmo_outpat, as.data.frame(hilmo_outpat_dict)[,c("service_type", "specialty", "hospital_type", "specialty_label", "cost")], by = c("service_type", "specialty", "hospital_type"))
hilmo_outpat = hilmo_outpat %>% filter(!is.na(cost)) %>%
	select(finngenid, source, service_type, specialty, hospital_type, specialty_label, event_age, approx_event_year, cost)

# MERGE HILMO DATA 
hilmo_new = rbind(hilmo_inpat, hilmo_outpat)

# MERGE AVOHILMO DATA
avohilmo_all = long_all[which(long_all$source == "PRIM_OUT"), c("finngenid", "code5", "code6", "code7", "event_age", "approx_event_year", "index")]
colnames(avohilmo_all) = c("finngenid", "contact_type", "service_type", "profession", "event_age", "approx_event_year", "index")
avohilmo_all = avohilmo_all %>% distinct() %>% select(-index)
avohilmo_new = merge(avohilmo_all, profession_dict, by = "profession", all.x = TRUE)
avohilmo_new = merge(avohilmo_new, service_dict, by = "service_type", all.x = TRUE)
avohilmo_new = merge(avohilmo_new, contact_dict, by = "contact_type", all.x = TRUE)
avohilmo_new = merge(avohilmo_new, avohilmo_dict, by = c("profession_label", "service_type_label", "contact_type_label"), all.x = TRUE)
avohilmo_new$profession_label = ifelse(is.na(avohilmo_new$profession_label), "other", avohilmo_new$profession_label)

# COMBINE DRUG, HILMO, AND AVOHILMO FOR SENSITIVITY GWAS
#temp = total
#temp$fu_end_age_minus_five = temp$fu_end_age - 5
#temp = temp[,c("finngenid", "fu_end_age_minus_five")]
#
#kela_temp = merge(kela_new, temp, by = "finngenid")
#hilmo_temp = merge(hilmo_temp, temp, by = "finngenid")
#avohilmo_temp = merge(avohilmo_temp, temp, by = "finngenid")
#
#kela_temp = kela_temp[which(kela_temp$event_age > kela_temp$fu_end_age_minus_five), -ncol(kela_temp), with = FALSE]
#hilmo_temp = hilmo_temp[which(hilmo_temp$event_age > hilmo_temp$fu_end_age_minus_five), -ncol(hilmo_temp), with = FALSE]
#avohilmo_temp = avohilmo_temp[which(avohilmo_temp$event_age > avohilmo_temp$fu_end_age_minus_five), -ncol(avohilmo_temp), with = FALSE]
#
#kela_temp = kela_new %>% filter(approx_event_year >= start_year_hilmo_kela, approx_event_year < end_year) %>% as.data.frame()
#hilmo_temp = hilmo_new %>% filter(approx_event_year >= start_year_hilmo_kela, approx_event_year < end_year) %>% as.data.frame()
#avohilmo_temp = avohilmo_new %>% filter(approx_event_year >= start_year_hilmo_kela, approx_event_year < end_year) %>% as.data.frame()
#
#long_new = rbindlist(list(setnames(kela_temp[,c("finngenid", "event_age", "cost")], c("finngenid", "event_age", "cost")),
#	hilmo_temp[,c("finngenid", "event_age", "cost")],
#	avohilmo_temp[,c("finngenid", "event_age", "cost")]))

# COMBINE DRUG, HILMO, AND AVOHILMO INTO LONG DATAFRAME
kela_temp = kela_new %>% filter(approx_event_year >= start_year_hilmo_kela, approx_event_year < end_year) %>% as.data.frame()
hilmo_temp = hilmo_new %>% filter(approx_event_year >= start_year_hilmo_kela, approx_event_year < end_year) %>% as.data.frame()
avohilmo_temp = avohilmo_new %>% filter(approx_event_year >= start_year_hilmo_kela, approx_event_year < end_year) %>% as.data.frame()
long_new = rbindlist(list(setnames(kela_temp[,c("finngenid", "event_age", "cost")], c("finngenid", "event_age", "cost")),
	hilmo_temp[,c("finngenid", "event_age", "cost")],
	avohilmo_temp[,c("finngenid", "event_age", "cost")]))

##################################################
########## CALCULATE HEALTHCARE COSTS
##################################################

# CALCULATE HEALTHCARE COSTS
kela_cost = kela_temp %>% group_by(finngenid) %>% summarise(kela_cost = sum(cost, na.rm = TRUE)) %>% as.data.frame()
hilmo_cost = hilmo_temp %>% group_by(finngenid) %>% summarise(hilmo_cost = sum(cost, na.rm = TRUE)) %>% as.data.frame()
avohilmo_cost = avohilmo_temp %>% group_by(finngenid) %>% summarise(avohilmo_cost = sum(cost, na.rm = TRUE)) %>% as.data.frame()

# ADJUST PERSON DAY FOR SENSITIVITY GWAS
total$person_day_in_hilmo_kela = ifelse(total$fu_end_year - 5 > start_year_hilmo_kela, 5 * 365.25, (total$fu_end_year - start_year_hilmo_kela) * 365.25)
total$person_day_in_avohilmo = ifelse(total$fu_end_year - 5 > start_year_avohilmo, 5 * 365.25, (total$fu_end_year - start_year_avohilmo) * 365.25)
total$person_day_in_avohilmo[which(total$person_day_in_avohilmo < 0)] = 0

# MERGE HEALTHCARE COSTS
total_new = merge(total, kela_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, hilmo_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, avohilmo_cost, by = "finngenid", all.x = TRUE)
dim(total_new) # 373460
total_new = total_new[!which(is.na(total_new$kela_cost) & is.na(total_new$hilmo_cost) & is.na(total_new$avohilmo_cost)),]
dim(total_new) # 373166

total_new$kela_cost_per_person_year = (total_new$kela_cost / total_new$person_day_in_hilmo_kela) * 365.25
total_new$hilmo_cost_per_person_year = (total_new$hilmo_cost / total_new$person_day_in_hilmo_kela) * 365.25
total_new$avohilmo_cost_per_person_year = (total_new$avohilmo_cost / total_new$person_day_in_avohilmo) * 365.25
total_new$avohilmo_cost_per_person_year[which(total_new$avohilmo_cost_per_person_year == Inf)] = 0
total_new$kela_cost_per_person_year[which(is.na(total_new$kela_cost_per_person_year))] = 0
total_new$hilmo_cost_per_person_year[which(is.na(total_new$hilmo_cost_per_person_year))] = 0
total_new$avohilmo_cost_per_person_year[which(((!is.na(total_new$kela_cost) | !is.na(total_new$hilmo_cost)) & is.na(total_new$avohilmo_cost)) | total_new$avohilmo_cost_per_person_year == Inf | is.na(total_new$avohilmo_cost_per_person_year))] = 
	median(total_new$avohilmo_cost_per_person_year, na.rm = TRUE)

total_new$log_kela_cost_per_person_year = log(total_new$kela_cost_per_person_year + 1)
total_new$log_hilmo_cost_per_person_year = log(total_new$hilmo_cost_per_person_year + 1)
total_new$log_avohilmo_cost_per_person_year = log(total_new$avohilmo_cost_per_person_year + 1)

total_new$total_cost = rowSums(total_new[,c("kela_cost", "hilmo_cost", "avohilmo_cost")], na.rm = TRUE)
total_new$total_cost_per_person_year = rowSums(total_new[,c("kela_cost_per_person_year", "hilmo_cost_per_person_year", "avohilmo_cost_per_person_year")], na.rm = TRUE)
total_new$log_total_cost_per_person_year = log(total_new$total_cost_per_person_year + 1)
fwrite(total_new, paste0(path, "healthcare_cost_five_year_2023123.txt"), row.names = FALSE, col.names = TRUE)

##################################################
########## FINALIZE HEALTHCARE COSTS
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}

# FINALIZE DATAFRAME
total_new = fread(file = paste0(path, "healthcare_cost_full_year_2023123.txt"), header = TRUE, stringsAsFactors = FALSE)
#fwrite(total_new, "/finngen/red/jiwoo/healthcare_cost.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# PLOT THE DIFFERENCE IN EUROS AT DIFFERENT COST PERCENTILES
#seq = seq(0, 0.99, by = 0.01)
#quantiles = quantile(total_new$total_cost_per_person_year, probs = seq)
#differences = quantiles * 0.2278
#data = as.data.frame(cbind(seq, differences))
#data$color = ifelse(data$seq == 0.5, "red", "black")
#plot(data$seq, data$differences, xlab = "Percentile", ylab = "Estimated Absolute Change in Euros", col = data$col, pch = 20)
#text(0.50, 1000, "298.99 euro increase per\n1 SD increase in WC for individuals\nin 50th percentile of WC", col = "red")

# MAKE GWAS DATAFRAME
cov_all = fread(file = "/finngen/library-red/finngen_R9/analysis_covariates/R9_COV_V1.FID.txt.gz", header = TRUE, stringsAsFactors = FALSE)
temp = total_new[,c("finngenid", "birth_year", "log_kela_cost_per_person_year", "log_hilmo_cost_per_person_year", "log_avohilmo_cost_per_person_year", "log_total_cost_per_person_year")]
cov_new = merge(cov_all, temp, by.x = "FID", by.y = "finngenid")
#cov_new$fu_end_age_squared = cov_new$fu_end_age ^ 2
cov_new$log_total_cost_per_female_year = ifelse(cov_new$SEX == "female", cov_new$log_total_cost_per_person_year, NA)
cov_new$log_total_cost_per_male_year = ifelse(cov_new$SEX == "male", cov_new$log_total_cost_per_person_year, NA)
cov_new$log_total_cost_per_person_year_0_30 = ifelse(cov_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP < 30, cov_new$log_total_cost_per_person_year, NA)
cov_new$log_total_cost_per_person_year_30_60 = ifelse(cov_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP >= 30 & cov_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP < 60, cov_new$log_total_cost_per_person_year, NA)
cov_new$log_total_cost_per_person_year_60_90 = ifelse(cov_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP >= 60, cov_new$log_total_cost_per_person_year, NA)
cov_new = cov_new %>% mutate_all(na_if," ")

# FORMAT PHENO FILE FOR SENSITIVITY GWAS
#cov_temp = cov_new[,c(1:3,7:119,140,142:151,168)]
#cov_temp = cov_temp[which(!is.na(cov_temp$log_total_cost_per_person_year)),]
#fwrite(cov_temp, paste0("/finngen/red/jiwoo/regenie_healthcare_cost/pheno/log_cost_five_age_agesq.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#cov_temp = fread(file = paste0("/finngen/red/jiwoo/regenie_healthcare_cost/pheno/log_cost_five_age_agesq.txt"), header = TRUE, stringsAsFactors = FALSE)
#fwrite(cov_temp, paste0("/finngen/red/jiwoo/regenie_healthcare_cost/pheno/log_cost_five_age_agesq.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
#fwrite(cov_new, "/finngen/red/jiwoo/regenie_healthcare_cost/pheno_five_years_2023123.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
#cov_new = fread("/finngen/red/jiwoo/regenie_healthcare_cost/pheno_five_years_2023123.txt", header = TRUE, stringsAsFactors = FALSE)
#fwrite(cov_new, "/finngen/red/jiwoo/regenie_healthcare_cost/pheno_five_years_2023123.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

phenolist = fread("/finngen/red/jiwoo/regenie_healthcare_cost/phenolist.txt", header = FALSE)
colnames(phenolist) = c("phenotype")
for (i in 1:length(phenolist$phenotype)) {
	print(phenolist$phenotype[i])
	write(phenolist$phenotype[i], paste0("/finngen/red/jiwoo/regenie_healthcare_cost/phenolist/", phenolist$phenotype[i], ".txt"))
	if (grepl("male", phenolist$phenotype[i])) {
		cov_temp = cov_new[,c(1:2,7:12,14:119,142:151,164,169)]
		cov_temp[[phenolist$phenotype[i]]] = cov_new[[phenolist$phenotype[i]]]
	} else {
		cov_temp = cov_new[,c(1:2,7:12,14:119,140,142:151,164,169)]
		cov_temp[[phenolist$phenotype[i]]] = cov_new[[phenolist$phenotype[i]]]
	}
	cov_temp = cov_temp[which(!is.na(cov_temp[[phenolist$phenotype[i]]]))]
	fwrite(cov_temp, paste0("/finngen/red/jiwoo/regenie_healthcare_cost/pheno/", phenolist$phenotype[i], ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	temp = fread(file = paste0("/finngen/red/jiwoo/regenie_healthcare_cost/pheno/", phenolist$phenotype[i], ".txt"), header = TRUE, stringsAsFactors = FALSE)
	fwrite(temp, paste0("/finngen/red/jiwoo/regenie_healthcare_cost/pheno/", phenolist$phenotype[i], ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
}

finngen-cli rw -w ./regenie.wdl -i ./json/log_total_cost_per_person_year_new.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_kela_cost_per_person_year.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_hilmo_cost_per_person_year.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_avohilmo_cost_per_person_year.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_total_cost_per_female_year.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_total_cost_per_male_year.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_total_cost_per_person_year_0_30.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_total_cost_per_person_year_30_60_new.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_total_cost_per_person_year_60_90.json -d ./regenie_sub_r8mem.zip

# SENSITIVITY GWAS
finngen-cli rw -w ./regenie.wdl -i ./json/log_cost_five_age_agesq.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_cost_five_age_agesq_birth.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_cost_full_age_agesq.json -d ./regenie_sub_r8mem.zip
finngen-cli rw -w ./regenie.wdl -i ./json/log_cost_full_age_agesq_birth.json -d ./regenie_sub_r8mem.zip

##################################################
########## PLOT DISTRIBUTION OF HEALTHCARE COSTS
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require(gridExtra)) {install.packages("gridExtra"); library(gridExtra)}

total_new = fread(file = "/finngen/red/jiwoo/regenie_healthcare_cost/pheno.txt", header = TRUE, stringsAsFactors = FALSE)

# PLOT DISTRIBUTION OF HEALTHCARE COSTS
options(scipen=999)

a = ggplot() +
	geom_histogram(data = total_new, mapping = aes(x = log_total_cost_per_person_year), fill = "black", alpha = 0.25) +
	geom_density(data = total_new, mapping = aes(x = log_total_cost_per_person_year, y = ..density..*nrow(total_new)*0.5), fill = "black", alpha = 0.5, adjust = 5) +
	scale_x_continuous(breaks = log(c(1, seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000), seq(10000, 100000, 10000), seq(100000, 1000000, 100000), seq(1000000, 10000000, 1000000))),
		labels = c(1, rep("", 9), 10, rep("", 9), 100, rep("", 9), 1000, rep("", 9), 10000, rep("", 9), 100000, rep("", 9), 1000000, rep("", 10))) +
	scale_y_continuous(limits = c(0, 70000), breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)) +
	labs(x = "", y = "", fill = "", tag = "A") + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 15, angle = 0, vjust = 0, color = "black"), axis.title.x = element_text(size = 15, color = "black"), axis.text.y = element_text(size = 15, color = "black"), axis.title.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

b = ggplot() +
	geom_histogram(data = total_new, mapping = aes(x = log_avohilmo_cost_per_person_year, fill = "Primary Care"), alpha = 0.25) +
	geom_density(data = total_new, mapping = aes(x = log_avohilmo_cost_per_person_year, y = ..density..*nrow(total_new)*0.5, fill = "Primary Care"), alpha = 0.5, adjust = 5) +
	geom_histogram(data = total_new, mapping = aes(x = log_hilmo_cost_per_person_year, fill = "Secondary Care"), alpha = 0.25) +
	geom_density(data = total_new, mapping = aes(x = log_hilmo_cost_per_person_year, y = ..density..*nrow(total_new)*0.5, fill = "Secondary Care"), alpha = 0.5, adjust = 5) +
	geom_histogram(data = total_new, mapping = aes(x = log_kela_cost_per_person_year, fill = "Medication"), alpha = 0.25) +
	geom_density(data = total_new, mapping = aes(x = log_kela_cost_per_person_year, y = ..density..*nrow(total_new)*0.5, fill = "Medication"), alpha = 0.5, adjust = 5) +
	scale_x_continuous(breaks = log(c(1, seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000), seq(10000, 100000, 10000), seq(100000, 1000000, 100000), seq(1000000, 10000000, 1000000))),
		labels = c(1, rep("", 9), 10, rep("", 9), 100, rep("", 9), 1000, rep("", 9), 10000, rep("", 9), 100000, rep("", 9), 1000000, rep("", 10))) +
	scale_y_continuous(limits = c(0, 70000), breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)) +
	scale_fill_manual(limits = c("Primary Care", "Secondary Care", "Medication"), values = c("skyblue", "dodgerblue", "navyblue")) +
	labs(x = "", y = "", fill = "", tag = "B") + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 15, angle = 0, vjust = 0, color = "black"), axis.title.x = element_text(size = 15, color = "black"), axis.text.y = element_text(size = 15, color = "black"), axis.title.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

c = ggplot() +
	geom_histogram(data = total_new[which(total_new$SEX == "female"),], mapping = aes(x = log_total_cost_per_person_year, fill = "Female"), alpha = 0.25) +
	geom_density(data = total_new[which(total_new$SEX == "female"),], mapping = aes(x = log_total_cost_per_person_year, y = ..density..*length(unique(total_new$FID[which(total_new$SEX == "female")]))*0.5, fill = "Female"), alpha = 0.5, adjust = 5) +
	geom_histogram(data = total_new[which(total_new$SEX == "male"),], mapping = aes(x = log_total_cost_per_person_year, fill = "Male"), alpha = 0.25) +
	geom_density(data = total_new[which(total_new$SEX == "male"),], mapping = aes(x = log_total_cost_per_person_year, y = ..density..*length(unique(total_new$FID[which(total_new$SEX == "male")]))*0.5, fill = "Male"), alpha = 0.5, adjust = 5) +
	scale_x_continuous(breaks = log(c(1, seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000), seq(10000, 100000, 10000), seq(100000, 1000000, 100000), seq(1000000, 10000000, 1000000))),
		labels = c(1, rep("", 9), 10, rep("", 9), 100, rep("", 9), 1000, rep("", 9), 10000, rep("", 9), 100000, rep("", 9), 1000000, rep("", 10))) +
	scale_y_continuous(limits = c(0, 70000), breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)) +
	scale_fill_manual(limits = c("Female", "Male"), values = c("dark green", "yellowgreen")) +
	labs(x = "", y = "", fill = "", tag = "C") + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 15, angle = 0, vjust = 0, color = "black"), axis.title.x = element_text(size = 15, color = "black"), axis.text.y = element_text(size = 15, color = "black"), axis.title.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_density(data = total_new[which(total_new$SEX == "female"),], mapping = aes(x = log_total_cost_per_person_year, y = ..density.., fill = "Female"), alpha = 0.5, adjust = 5) +
	geom_density(data = total_new[which(total_new$SEX == "male"),], mapping = aes(x = log_total_cost_per_person_year, y = ..density.., fill = "Male"), alpha = 0.5, adjust = 5)

d = ggplot() +
	geom_histogram(data = total_new[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP > 60)], mapping = aes(x = log_total_cost_per_person_year, fill = ">60 years old"), alpha = 0.5) +
	geom_density(data = total_new[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP > 60)], mapping = aes(x = log_total_cost_per_person_year, y = ..density..*length(unique(total_new$FID[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP > 60)]))*0.5, fill = ">60 years old"), alpha = 0.25, adjust = 5) +
	geom_histogram(data = total_new[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP < 60 & total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP >= 30)], mapping = aes(x = log_total_cost_per_person_year, fill = "30-60 years old"), alpha = 0.5) +
	geom_density(data = total_new[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP < 60 & total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP >= 30)], mapping = aes(x = log_total_cost_per_person_year, y = ..density..*length(unique(total_new$FID[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP < 60 & total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP >= 30)]))*0.5, fill = "30-60 years old"), alpha = 0.25, adjust = 5) +
	geom_histogram(data = total_new[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP < 30)], mapping = aes(x = log_total_cost_per_person_year, fill = "<30 years old"), alpha = 0.5) +
	geom_density(data = total_new[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP < 30)], mapping = aes(x = log_total_cost_per_person_year, y = ..density..*length(unique(total_new$FID[which(total_new$AGE_AT_DEATH_OR_END_OF_FOLLOWUP < 30)]))*0.5, fill = "<30 years old"), alpha = 0.25, adjust = 5) +
	scale_x_continuous(breaks = log(c(1, seq(1, 10, 1), seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000), seq(10000, 100000, 10000), seq(100000, 1000000, 100000), seq(1000000, 10000000, 1000000))),
		labels = c(1, rep("", 9), 10, rep("", 9), 100, rep("", 9), 1000, rep("", 9), 10000, rep("", 9), 100000, rep("", 9), 1000000, rep("", 10))) +
	scale_y_continuous(limits = c(0, 70000), breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)) +
	scale_fill_manual(limits = c("<30 years old", "30-60 years old", ">60 years old"), values = c("coral", "red", "red4")) +
	labs(x = "", y = "", fill = "", tag = "D") + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 15, angle = 0, vjust = 0, color = "black"), axis.title.x = element_text(size = 15, color = "black"), axis.text.y = element_text(size = 15, color = "black"), axis.title.y = element_text(size = 15, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

grid.arrange(a, b, c, d, left = textGrob("Number of Individuals", rot = 90, gp = gpar(fontsize = 25)), bottom = textGrob("Annual Healthcare Cost in Euros", gp = gpar(fontsize = 25)))

##################################################
########## CALCLUATE PRS FOR UK BIOBANK
##################################################

finngen-cli request-workflow -w /finngen/red/jiwoo/prs_healthcare_cost/prs.wdl -i /finngen/red/jiwoo/prs_healthcare_cost/prs.json

##################################################
########## HARMONIZE FINLAND
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20230123"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

sumstats = c(
	"log_avohilmo_cost_per_person_year",
	"log_hilmo_cost_per_person_year",
	"log_kela_cost_per_person_year",
	"log_total_cost_per_person_year",
	"log_total_cost_per_female_year",
	"log_total_cost_per_male_year",
	"log_total_cost_per_person_year_0_30",
	"log_total_cost_per_person_year_30_60",
	"log_total_cost_per_person_year_60_90")

sumstats = c("log_cost_five_age_agesq", "log_cost_five_age_agesq_birth", "log_cost_full_age_agesq", "log_cost_full_age_agesq_birth")

options(scipen = 999)

for (i in 1:length(sumstats)) {
	sumstat = fread(file = paste0(path, "gwas/", today, "/", sumstats[i]), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
	colnames(sumstat)[1] = "chr"
	sumstat$cpra = paste0("chr", sumstat$chr, "_", sumstat$pos, "_", sumstat$ref, "_", sumstat$alt)
	sumstat = sumstat[,c("chr", "pos", "cpra", "ref", "alt", "pval", "mlogp", "beta", "sebeta", "af_alt")]
	fwrite(sumstat, paste0(path, "gwas/", today, "/", sumstats[i], ".txt.gz"), sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE, na = NA)
	print(sumstats[i])
}

awk 'function nn(x){split(x,xx,"");str=n[xx[1]];for(i=2;i<=length(xx);i++){str=str""n[xx[i]]};return str}BEGIN{n["A"]=1;n["C"]=2;n["G"]=3;n["T"]=4}NR==FNR{split($3,a,"_");split(a[1],b,"hr");id=b[2]":"a[2]":"a[3]":"a[4];if(nn(a[4])<nn(a[3])){id=b[2]":"a[2]":"a[4]":"a[3]};cpra[id]++;next}($2 in cpra){print $1" "$2}' <(zcat log_cost_full_age_agesq.txt.gz) <(zcat /medpop/esp2/jiwoolee/cost_rep/dbsnp_151_b38_rsid_cpra.txt.gz) > dbsnp_151_b38_rsid_cpra.txt

awk 'BEGIN{OFS=" "}NR==FNR{rsid[$2]=$1;next}FNR==1{print;next}{split($3,a,"_");split(a[1],b,"hr");cpra=b[2]":"a[2]":"a[3]":"a[4];cpar=b[2]":"a[2]":"a[4]":"a[3]}(cpra in rsid){$3=rsid[cpra];print;next}(cpar in rsid){$3=rsid[cpar];print;next}{print}' dbsnp_151_b38_rsid_cpra.txt <(zcat log_cost_full_age_agesq.txt.gz) | gzip --best > log_cost_full_age_agesq_final.txt.gz

##################################################
########## HARMONIZE NETHERLANDS
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/
R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

options(scipen = 999)

total_neth = fread(file = paste0(path, "gwas/", today, "/netherlands_total.fastGWA"), header = TRUE, stringsAsFactors = FALSE)
total_neth = total_neth[!which(total_neth$CHR == "CHR"),]
total_neth$CHR = paste0("chr", total_neth$CHR)
total_neth$POSPOS = as.numeric(total_neth$POS) + 1
total_neth$POSPOS = as.character(total_neth$POSPOS)
total_neth$hg37 = paste0(total_neth$CHR, ":", total_neth$POS, "_", total_neth$POSPOS)
fwrite(total_neth[,c("CHR", "POS", "POSPOS")], paste0(path, "netherlands_hg37.bed"), sep = "\t", col.names = FALSE, row.names = FALSE)

./liftOver netherlands_hg37.bed hg19ToHg38.over.chain netherlands_hg38.bed unMapped

hg37 = fread(file = paste0(path, "netherlands_hg37.bed"), header = FALSE, stringsAsFactors = FALSE)
hg38 = fread(file = paste0(path, "netherlands_hg38.bed"), header = FALSE, stringsAsFactors = FALSE)
unmapped = fread(file = paste0(path, "unMapped"), fill = TRUE, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
unmapped = unmapped[!grepl("#Deleted", unmapped$V1),]
fwrite(unmapped, paste0(path, "unMapped_cleaned"))
unmapped_cleaned = fread(paste0(path, "unMapped_cleaned"), sep = "\t")

hg37$V4 = paste0(hg37$V1, ":", hg37$V2, "_", hg37$V3)
hg38$V4 = paste0(hg38$V1, ":", hg38$V2, "_", hg38$V3)
unmapped_cleaned$V4 = paste0(unmapped_cleaned$V1, ":", unmapped_cleaned$V2, "_", unmapped_cleaned$V3)

hg37 = hg37[!which(hg37$V4 %in% unmapped_cleaned$V4),]
hg37_hg38 = cbind(hg37[,"V4"], hg38)
colnames(hg37_hg38) = c("hg37", "hg38_chr", "hg38_start", "hg38_end", "hg38")

total_new = merge(total_neth, hg37_hg38, by = "hg37")
total_temp = total_new[,c("hg38_chr", "hg38_start", "hg38", "A2", "A1", "AF1")]
total_temp$hg38_chr = substring(total_temp$hg38_chr, 4)
temp = strsplit(total_temp$hg38, "_")
total_temp$hg38 = sapply(temp, head, 1) 
total_temp$hg38 = paste0(total_temp$hg38, "_", total_temp$A2, "_", total_temp$A1)
colnames(total_temp) = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ")
total_temp$ID = paste0("chr", total_temp$CHROM, "_", total_temp$GENPOS, "_", total_temp$ALLELE0, "_", total_temp$ALLELE1)
fwrite(total_temp, paste0(path, "netherlands_hg38.txt"), sep = " ")

awk 'function nn(x){split(x,xx,"");str=n[xx[1]];for(i=2;i<=length(xx);i++){str=str""n[xx[i]]};return str}BEGIN{n["A"]=1;n["C"]=2;n["G"]=3;n["T"]=4}NR==FNR{split($3,a,"_");split(a[1],b,"hr");id=b[2]":"a[2]":"a[3]":"a[4];if(nn(a[4])<nn(a[3])){id=b[2]":"a[2]":"a[4]":"a[3]};cpra[id]++;next}($2 in cpra){print $1" "$2}' <(zcat log_avohilmo_cost_per_person_year.gz) <(zcat /medpop/esp2/jiwoolee/cost_rep/dbsnp_151_b38_rsid_cpra.txt.gz) > dbsnp_151_b38_rsid_cpra_FinnGen_GWAS_only.txt

awk 'BEGIN{OFS=" "}NR==FNR{rsid[$2]=$1;next}FNR==1{print;next}{split($3,a,"_");split(a[1],b,"hr");cpra=b[2]":"a[2]":"a[3]":"a[4];cpar=b[2]":"a[2]":"a[4]":"a[3]}(cpra in rsid){$3=rsid[cpra];print;next}(cpar in rsid){$3=rsid[cpar];print;next}{print}' dbsnp_151_b38_rsid_cpra_FinnGen_GWAS_only.txt <(zcat log_avohilmo_cost_per_person_year.gz) | gzip --best > log_avohilmo_cost_per_person_year_final.gz

total_neth = fread(file = paste0(path, "gwas/", today, "/netherlands_total.fastGWA"), header = TRUE, stringsAsFactors = FALSE)
pharmacy = fread(file = paste0(path, "gwas/", today, "/netherlands_pharmacy.fastGWA"), header = TRUE, stringsAsFactors = FALSE)
mental = fread(file = paste0(path, "gwas/", today, "/netherlands_mental_health.fastGWA"), header = TRUE, stringsAsFactors = FALSE)
hospital = fread(file = paste0(path, "gwas/", today, "/netherlands_hospital.fastGWA"), header = TRUE, stringsAsFactors = FALSE)
gp = fread(file = paste0(path, "gwas/", today, "/netherlands_gp.fastGWA"), header = TRUE, stringsAsFactors = FALSE)

hg37_hg38$CHROM = substr(hg37_hg38$hg38_chr, 4, nchar(hg37_hg38$hg38_chr))
hg37_hg38$hg38_start = as.character(hg37_hg38$hg38_start)

rsid = fread(file = paste0(path, "/netherlands_hg38_rsid.txt.gz"), header = TRUE)
rsid = rsid[,c("CHROM", "GENPOS", "ID")]
rsid$GENPOS = as.character(rsid$GENPOS)

my_map = merge(hg37_hg38, rsid, by.x = c("CHROM", "hg38_start"), by.y = c("CHROM", "GENPOS"))
temp = strsplit(my_map$hg37, ":")
my_map$hg37_chr = sapply(temp, head, 1) 
my_map$hg37_chr = substring(my_map$hg37_chr, 4)
temp = strsplit(my_map$hg37, "_")
my_map$hg37_pos = sapply(temp, tail, 1) 
my_map$hg37_pos = as.character(as.numeric(my_map$hg37_pos)-1)
my_new_map = my_map[,c("hg37_chr", "hg37_pos", "ID")]

neth_new = merge(total_neth, my_new_map, by.x = c("CHR", "POS"), by.y = c("hg37_chr", "hg37_pos"))
pharmacy_new = merge(pharmacy, my_new_map, by.x = c("CHR", "POS"), by.y = c("hg37_chr", "hg37_pos"))
mental_new = merge(mental, my_new_map, by.x = c("CHR", "POS"), by.y = c("hg37_chr", "hg37_pos"))
hospital_new = merge(hospital, my_new_map, by.x = c("CHR", "POS"), by.y = c("hg37_chr", "hg37_pos"))
gp_new = merge(gp, my_new_map, by.x = c("CHR", "POS"), by.y = c("hg37_chr", "hg37_pos"))
neth_new = neth_new[,-"SNP"]
pharmacy_new = pharmacy_new[,-"SNP"]
mental_new = mental_new[,-"SNP"]
hospital_new = hospital_new[,-"SNP"]
gp_new = gp_new[,-"SNP"]
neth_new$beta = as.numeric(neth_new$beta)
pharmacy_new$beta = as.numeric(pharmacy_new$beta)
mental_new$beta = as.numeric(mental_new$beta)
hospital_new$beta = as.numeric(hospital_new$beta)
gp_new$beta = as.numeric(gp_new$beta)
fwrite(neth_new, paste0(path, "gwas/", today, "/netherlands_total.txt"), sep = " ")
fwrite(pharmacy_new, paste0(path, "gwas/", today, "/netherlands_pharmacy.txt"), sep = " ")
fwrite(mental_new, paste0(path, "gwas/", today, "/netherlands_mental.txt"), sep = " ")
fwrite(hospital_new, paste0(path, "gwas/", today, "/netherlands_hospital.txt"), sep = " ")
fwrite(gp_new, paste0(path, "gwas/", today, "/netherlands_gp.txt"), sep = " ")

##################################################
########## MUNGE SUMMARY STATISTICS
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/ldsc
source activate ldsc

python munge_sumstats.py \
--sumstats /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/log_total_cost_per_person_year_final.txt.gz \
--N 373160 \
--snp cpra \
--a1 alt \
--a2 ref \
--p pval \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/log_total_cost_per_person_year_final \
--merge-alleles ./eur_w_ld_chr/w_hm3.snplist

python munge_sumstats.py \
--sumstats /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/ukbb_secondary_care_costs_imputed.txt \
--N 307048 \
--snp SNP \
--a1 ALLELE1 \
--a2 ALLELE0 \
--p P_LINREG \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/ukbb_secondary_care_costs_imputed \
--merge-alleles ./eur_w_ld_chr/w_hm3.snplist

python munge_sumstats.py \
--sumstats /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/netherlands_total.txt \
--snp ID \
--a1 A1 \
--a2 A2 \
--p p \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/netherlands_total \
--merge-alleles ./eur_w_ld_chr/w_hm3.snplist

#git checkout b02f2a6 # revert to a previous version because of value error of concantenation
python munge_sumstats.py \
--sumstats /medpop/esp2/jiwoolee/cost_rep/gwas/20230123/log_cost_full_age_agesq_birth_final.txt.gz \
--N 373160 \
--snp cpra \
--a1 alt \
--a2 ref \
--p pval \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20230123/log_cost_full_age_agesq_birth_final \
--merge-alleles ./eur_w_ld_chr/w_hm3.snplist

##################################################
########## CALCULATE GENETIC CORRELATION
##################################################

python ldsc.py \
--rg /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/log_avohilmo_cost_per_person_year_final.sumstats.gz,/medpop/esp2/jiwoolee/cost_rep/gwas/20220505/netherlands_gp.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/primary

python ldsc.py \
--rg /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/log_hilmo_cost_per_person_year_final.sumstats.gz,/medpop/esp2/jiwoolee/cost_rep/gwas/20220505/test.sumstats.gz,/medpop/esp2/jiwoolee/cost_rep/gwas/20220505/netherlands_hospital.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/test

python ldsc.py \
--rg /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/log_kela_cost_per_person_year_final.sumstats.gz,/medpop/esp2/jiwoolee/cost_rep/gwas/20220505/netherlands_pharmacy.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/medication

python ldsc.py \
--rg /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/log_total_cost_per_person_year_final.sumstats.gz,/medpop/esp2/jiwoolee/cost_rep/gwas/20220505/netherlands_total.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20220505/total





python ldsc.py \
--rg /medpop/esp2/jiwoolee/cost_rep/gwas/20230123/log_cost_five_age_agesq.txt.gz,/medpop/esp2/jiwoolee/cost_rep/gwas/20230123/log_cost_five_age_agesq_birth.txt.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out /medpop/esp2/jiwoolee/cost_rep/gwas/20230123/fives

##################################################
########## FORMAT PHENOTYPE FILE
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(TwoSampleMR)) {install.packages("TwoSampleMR"); library(TwoSampleMR)}
setDTthreads(0)

# List available GWAS
#ao <- available_outcomes()
#ao_df <- as.data.frame(ao)

phenos = fread(paste0(path, "phenotypes.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
exposure_id = phenos$gwas_id

glgc_hdl = fread(paste0(path, "gwas/", today, "/nonFinnish_only_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results"), header = TRUE, stringsAsFactors = FALSE)
glgc_hdl = format_data(glgc_hdl, type = "exposure", snp_col = "rsID", beta_col = "EFFECT_SIZE", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF", eaf_col = "POOLED_ALT_AF", pval_col = "pvalue")
glgc_hdl_new = glgc_hdl[which(glgc_hdl$pval.exposure < 5e-8),]
glgc_hdl_new = clump_data(glgc_hdl_new)
glgc_ldl = fread(paste0(path, "gwas/", today, "/nonFinnish_only_LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results"), header = TRUE, stringsAsFactors = FALSE)
glgc_ldl = format_data(glgc_ldl, type = "exposure", snp_col = "rsID", beta_col = "EFFECT_SIZE", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF", eaf_col = "POOLED_ALT_AF", pval_col = "pvalue")
glgc_ldl_new = glgc_ldl[which(glgc_ldl$pval.exposure < 5e-8),]
glgc_ldl_new = clump_data(glgc_ldl_new)
glgc_tg = fread(paste0(path, "gwas/", today, "/nonFinnish_only_logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results"), header = TRUE, stringsAsFactors = FALSE)
glgc_tg = format_data(glgc_tg, type = "exposure", snp_col = "rsID", beta_col = "EFFECT_SIZE", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF", eaf_col = "POOLED_ALT_AF", pval_col = "pvalue")
glgc_tg_new = glgc_tg[which(glgc_tg$pval.exposure < 5e-8),]
glgc_tg_new = clump_data(glgc_tg_new)
glgc = c("ieu-a-299", "ieu-a-300", "ieu-a-302")

f_res = NULL
for (i in 1:length(exposure_id)) {
	print(exposure_id[i])
	if (exposure_id[i] == glgc[1]) {
		exposure_dat = glgc_hdl_new		
	} else if (exposure_id[i] == glgc[2]) {
		exposure_dat = glgc_ldl_new
	} else if (exposure_id[i] == glgc[3]) {
		exposure_dat = glgc_tg_new
	} else {
		exposure_dat = extract_instruments(outcomes = exposure_id[i])
	}
	exposure_dat$beta.exposure = exposure_dat$beta.exposure * phenos$transformation[i]
	exposure_dat$se.exposure = exposure_dat$se.exposure * phenos$transformation[i]
	n = as.numeric(phenos[i, "sample_size"])
	k = nrow(exposure_dat)
	
	exposure_dat$maf = ifelse(exposure_dat$eaf.exposure < 0.50, exposure_dat$eaf.exposure, 1 - exposure_dat$eaf.exposure)
	exposure_dat$rsq = 2*exposure_dat$beta.exposure^2*exposure_dat$maf*(1-exposure_dat$maf)/(2*exposure_dat$beta.exposure^2*exposure_dat$maf*(1-exposure_dat$maf) + exposure_dat$se.exposure^2*2*n*exposure_dat$maf*(1-exposure_dat$maf))

	rsq = sum(exposure_dat$rsq, na.rm = TRUE)
	f = ((n - k - 1) / k) * (rsq/(1-rsq))
	res = c(phenos[i, "my_id"], n, k, rsq, f)
	f_res = rbind(f_res, res)
}
f_res = as.data.frame(f_res)
colnames(f_res) = c("my_id", "n", "k", "rsq", "f")
f_temp = f_res[,c("my_id", "k", "rsq", "f")]
colnames(f_temp) = c("my_id", "n_var", "rsq", "f_statistic")
f_temp$my_id = as.character(f_temp$my_id)
phenos_new = merge(phenos, f_temp, by = "my_id")
fwrite(phenos_new, paste0(path, "/phenotypes_new.csv"), col.names = TRUE, row.names = FALSE, quote = FALSE)

# sumstat = fread(paste0(path, "gwas/", today, "/log_total_cost_per_person_year_final.txt.gz"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)

##################################################
########## DO MENDELIAN RANDOMIZATION	FOR FINLAND
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(TwoSampleMR)) {install.packages("TwoSampleMR"); library(TwoSampleMR)}

setDTthreads(0)

# List available GWAS
#ao <- available_outcomes()
#ao_df <- as.data.frame(ao)

phenos = fread(paste0(path, "phenotypes.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
exposure_id = phenos$my_id

outcome_id = c(
	"log_avohilmo_cost_per_person_year",
	"log_hilmo_cost_per_person_year",
	"log_kela_cost_per_person_year",
	"log_total_cost_per_person_year",
	"log_total_cost_per_female_year",
	"log_total_cost_per_male_year",
	"log_total_cost_per_person_year_0_30",
	"log_total_cost_per_person_year_30_60",
	"log_total_cost_per_person_year_60_90"
)

glgc_hdl = fread(paste0(path, "gwas/", today, "/nonFinnish_only_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results"), header = TRUE, stringsAsFactors = FALSE)
glgc_hdl = format_data(glgc_hdl, type = "exposure", snp_col = "rsID", beta_col = "EFFECT_SIZE", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF", eaf_col = "POOLED_ALT_AF", pval_col = "pvalue")
glgc_hdl_new = glgc_hdl[which(glgc_hdl$pval.exposure < 5e-8),]
glgc_hdl_new = clump_data(glgc_hdl_new)
glgc_ldl = fread(paste0(path, "gwas/", today, "/nonFinnish_only_LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results"), header = TRUE, stringsAsFactors = FALSE)
glgc_ldl = format_data(glgc_ldl, type = "exposure", snp_col = "rsID", beta_col = "EFFECT_SIZE", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF", eaf_col = "POOLED_ALT_AF", pval_col = "pvalue")
glgc_ldl_new = glgc_ldl[which(glgc_ldl$pval.exposure < 5e-8),]
glgc_ldl_new = clump_data(glgc_ldl_new)
glgc_tg = fread(paste0(path, "gwas/", today, "/nonFinnish_only_logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results"), header = TRUE, stringsAsFactors = FALSE)
glgc_tg = format_data(glgc_tg, type = "exposure", snp_col = "rsID", beta_col = "EFFECT_SIZE", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF", eaf_col = "POOLED_ALT_AF", pval_col = "pvalue")
glgc_tg_new = glgc_tg[which(glgc_tg$pval.exposure < 5e-8),]
glgc_tg_new = clump_data(glgc_tg_new)
glgc = c("ieu-a-299", "ieu-a-300", "ieu-a-302")

mr_res = NULL
for (j in 1:length(outcome_id)) {
	print(outcome_id[j])
	sumstat = fread(paste0(path, "gwas/", today, "/", outcome_id[j], "_final.txt.gz"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)

	for (i in 1:length(exposure_id)) {
		print(exposure_id[i])
		# Get instruments
		#print("Loading exposure data...")
		#if (phenos$gwas_id[i] == glgc[1]) {
		#	next
		#	exposure_dat = glgc_hdl_new		
		#	exposure_dat$id.exposure = glgc[1]
		#} else if (phenos$gwas_id[i] == glgc[2]) {
		#	next
		#	exposure_dat = glgc_ldl_new
		#	exposure_dat$id.exposure = glgc[2]
		#} else if (phenos$gwas_id[i] == glgc[3]) {
		#	next
		#	exposure_dat = glgc_tg_new
		#	exposure_dat$id.exposure = glgc[3]
		#} else {
			exposure_dat = extract_instruments(outcomes = phenos$gwas_id[i])
		#}
		exposure_dat$beta.exposure = exposure_dat$beta.exposure * phenos$transformation[i]
		exposure_dat$se.exposure = exposure_dat$se.exposure * phenos$transformation[i]
		# Get effects of instruments on outcome
		#print("Loading outcome data...")
		outcome_temp = sumstat[which(sumstat$cpra %in% exposure_dat$SNP),]
		fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
		outcome_dat = read_outcome_data(snps = NULL,
			filename = paste0(path, "mr/", today, "/temp.csv"),
			sep = ",",
			snp_col = "cpra",
			beta_col = "beta",
		  se_col = "sebeta",
		  effect_allele_col = "alt",
		  other_allele_col = "ref",
		  eaf_col = "af_alt",
		  pval_col = "pval"
		)
		# Harmonise the exposure and outcome data
		#print("Harmonizing data...")
		dat = harmonise_data(exposure_dat, outcome_dat)
		# Perform MR
		#print("Performing MR...")
		res = mr(dat)
		if (outcome_id[j] == "log_total_cost_per_person_year") {
			print("Plotting...")
			scatter = mr_scatter_plot(res, dat)
			#ggsave(paste0(path, "mr/", today, "/", phenos$my_id[i], "_scatter.jpg"), scatter[[1]], device = "jpg", type = "cairo", height = 10, width = 10)
			#dev.off()
		}
		res$id.exposure = exposure_id[i]
		res$exposure = phenos$phenotype[i]
		res$outcome = outcome_id[j]
		mr_res = rbind(mr_res, res)
	}
	rm(sumstat)
}

mr_temp = mr_res
mr_temp = as.data.frame(mr_temp)
mr_temp$b = as.numeric(mr_temp$b)
mr_temp$se = as.numeric(mr_temp$se)
mr_temp$pval = as.numeric(mr_temp$pval)
mr_temp$percent = round((exp(mr_temp$b) - 1) * 100, 2)
mr_temp$confint_lower = round((exp(mr_temp$b - 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint_upper = round((exp(mr_temp$b + 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint = paste0("[", mr_temp$confint_lower, ", ", mr_temp$confint_upper, "]")
mr_temp = mr_temp[order(mr_temp$pval),]
fwrite(mr_temp, paste0(path, "mr/", today, "/fi_mr_res_total_new.csv"))

##################################################
########## DO MENDELIAN RANDOMIZATION	FOR UNITED KINGDOM AND NETHERLANDS
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(TwoSampleMR)) {install.packages("TwoSampleMR"); library(TwoSampleMR)}

setDTthreads(0)

# List available GWAS
#ao <- available_outcomes()
#ao_df <- as.data.frame(ao)

phenos = fread(paste0(path, "phenotypes.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
phenos = phenos[which(phenos$my_id %in% c("ieu-a-835_irnt","ieu-a-835_raw","ukb-a-360_irnt","ukb-a-360_raw","ukb-a-382_irnt","ukb-a-382_raw","ieu-a-61_irnt","ieu-a-61_raw","ieu-b-38_irnt","ieu-b_raw")),]
exposure_id = phenos$gwas_id

sumstat = fread(paste0(path, "gwas/", today, "/netherlands_total.txt"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
mr_res = NULL
for (i in 1:length(exposure_id)) {
	print(exposure_id[i])
	# Get instruments
	print("Loading exposure data...")
	exposure_dat = extract_instruments(outcomes = exposure_id[i])
	exposure_dat$beta.exposure = exposure_dat$beta.exposure * phenos$transformation[i]
	exposure_dat$se.exposure = exposure_dat$se.exposure * phenos$transformation[i]
	# Get effects of instruments on outcome
	print("Loading outcome data...")
	outcome_temp = sumstat[which(sumstat$ID %in% exposure_dat$SNP),]
	fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
	outcome_dat = read_outcome_data(snps = NULL,
		filename = paste0(path, "mr/", today, "/temp.csv"),
		sep = ",",
		snp_col = "ID",
		beta_col = "beta",
	  se_col = "se",
	  effect_allele_col = "A1",
	  other_allele_col = "A2",
	  eaf_col = "AF1",
	  pval_col = "p"
	)
	# Harmonise the exposure and outcome data
	print("Harmonizing data...")
	dat = harmonise_data(exposure_dat, outcome_dat)
	# Perform MR
	print("Performing MR...")
	res = mr(dat)
	if (phenos$transformation[i] != 1) {
		res$id.exposure = paste0(res$id.exposure, "_raw")
	}
	res$exposure = phenos$phenotype[i]
	res$outcome = "netherlands"
	mr_res = rbind(mr_res, res)
}
mr_temp = mr_res
mr_temp = as.data.frame(mr_temp)
mr_temp$b = as.numeric(mr_temp$b)
mr_temp$se = as.numeric(mr_temp$se)
mr_temp$pval = as.numeric(mr_temp$pval)
mr_temp$percent = round((exp(mr_temp$b) - 1) * 100, 2)
mr_temp$confint_lower = round((exp(mr_temp$b - 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint_upper = round((exp(mr_temp$b + 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint = paste0("[", mr_temp$confint_lower, ", ", mr_temp$confint_upper, "]")
mr_temp = mr_temp[order(mr_temp$pval),]
fwrite(mr_temp, paste0(path, "mr/", today, "/nl_mr_res_total.csv"))


sumstat = fread(paste0(path, "gwas/", today, "/ukbb_secondary_care_costs_imputed.txt"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
mr_res = NULL
for (i in 1:length(exposure_id)) {
	print(exposure_id[i])
	# Get instruments
	print("Loading exposure data...")
	exposure_dat = extract_instruments(outcomes = exposure_id[i])
	exposure_dat$beta.exposure = exposure_dat$beta.exposure * phenos$transformation[i]
	exposure_dat$se.exposure = exposure_dat$se.exposure * phenos$transformation[i]
	# Get effects of instruments on outcome
	print("Loading outcome data...")
	outcome_temp = sumstat[which(sumstat$SNP %in% exposure_dat$SNP),]
	fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
	outcome_dat = read_outcome_data(snps = NULL,
		filename = paste0(path, "mr/", today, "/temp.csv"),
		sep = ",",
		snp_col = "SNP",
		beta_col = "BETA",
	  se_col = "SE",
	  effect_allele_col = "ALLELE1",
	  other_allele_col = "ALLELE0",
	  eaf_col = "A1FREQ",
	  pval_col = "P_LINREG"
	)
	# Harmonise the exposure and outcome data
	print("Harmonizing data...")
	dat = harmonise_data(exposure_dat, outcome_dat)
	# Perform MR
	print("Performing MR...")
	res = mr(dat)
	if (phenos$transformation[i] != 1) {
		res$id.exposure = paste0(res$id.exposure, "_raw")
	}
	res$exposure = phenos$phenotype[i]
	res$outcome = "ukbb_secondary"
	mr_res = rbind(mr_res, res)
}
mr_temp = mr_res
mr_temp = as.data.frame(mr_temp)
mr_temp$b = as.numeric(mr_temp$b)
mr_temp$se = as.numeric(mr_temp$se)
mr_temp$pval = as.numeric(mr_temp$pval)
mr_temp$percent = round((exp(mr_temp$b) - 1) * 100, 2)
mr_temp$confint_lower = round((exp(mr_temp$b - 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint_upper = round((exp(mr_temp$b + 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint = paste0("[", mr_temp$confint_lower, ", ", mr_temp$confint_upper, "]")
mr_temp = mr_temp[order(mr_temp$pval),]
fwrite(mr_temp, paste0(path, "mr/", today, "/uk_mr_res_total.csv"))
	
##################################################
########## DO MULTIVARIABLE MENDELIAN RANDOMIZATION	
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(TwoSampleMR)) {install.packages("TwoSampleMR"); library(TwoSampleMR)}

setDTthreads(0)

sumstat = fread(paste0(path, "gwas/", today, "/log_total_cost_per_person_year_final.txt.gz"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)

mr_res = NULL
exposure_list = c("ukb-a-360", "ieu-a-835", "ukb-a-382")
exposure_verbose_list = c("Systolic blood pressure", "Body mass index", "Waist circumference")
mediator_list = c("ukb-a-473", "ukb-a-534", "ukb-a-306", "ukb-a-543", "ukb-d-C_STROKE", "ukb-b-18009")
mediator_verbose_list = c("Back pain", "Chronic ischemic heart disease", "Type 2 diabetes", "Chronic obstructive pulmonary disease", "Stroke", "Blood pressure medication")
for (i in 1:length(exposure_list)) {
	print(exposure_verbose_list[i])
	exposure_dat = extract_instruments(outcomes = exposure_list[i])
	outcome_temp = sumstat[which(sumstat$cpra %in% exposure_dat$SNP),]
	fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
	outcome_dat = read_outcome_data(snps = NULL, filename = paste0(path, "mr/", today, "/temp.csv"), sep = ",", snp_col = "cpra", beta_col = "beta", se_col = "sebeta", effect_allele_col = "alt", other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval")
	dat = harmonise_data(exposure_dat, outcome_dat)
	res = mr(dat, method_list = c("mr_ivw")) 
	res = res[,c("id.exposure", "exposure", "id.outcome", "outcome", "nsnp", "b", "se", "pval", "method")]
	res$method = exposure_verbose_list[i]
	mr_res = rbind(mr_res, res)
	for (j in 1:length(mediator_list)) {
		exposure_id = c(exposure_list[i], mediator_list[j])
		exposure_dat = mv_extract_exposures(exposure_id)
		outcome_temp = sumstat[which(sumstat$cpra %in% exposure_dat$SNP),]
		fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
		outcome_dat = read_outcome_data(snps = NULL, filename = paste0(path, "mr/", today, "/temp.csv"), sep = ",", snp_col = "cpra", beta_col = "beta", se_col = "sebeta", effect_allele_col = "alt", other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval")
		dat = mv_harmonise_data(exposure_dat, outcome_dat)
		res = as.data.frame(mv_multiple(dat))
		res$type = paste0(exposure_verbose_list[i], " (", mediator_verbose_list[j], ")")
		colnames(res) = c("id.exposure", "exposure", "id.outcome", "outcome", "nsnp", "b", "se", "pval", "method")
		mr_res = rbind(mr_res, res)
	}
	exposure_id = c(exposure_list[i], mediator_list)
	exposure_dat = mv_extract_exposures(exposure_id)
	outcome_temp = sumstat[which(sumstat$cpra %in% exposure_dat$SNP),]
	fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
	outcome_dat = read_outcome_data(snps = NULL, filename = paste0(path, "mr/", today, "/temp.csv"), sep = ",", snp_col = "cpra", beta_col = "beta", se_col = "sebeta", effect_allele_col = "alt", other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval")
	dat = mv_harmonise_data(exposure_dat, outcome_dat)
	res = as.data.frame(mv_multiple(dat))
	res$type = paste0(exposure_verbose_list[i], " (All)")
	colnames(res) = c("id.exposure", "exposure", "id.outcome", "outcome", "nsnp", "b", "se", "pval", "method")
	mr_res = rbind(mr_res, res)
}
mr_temp = mr_res
mr_temp = as.data.frame(mr_temp)
mr_temp$b = as.numeric(mr_temp$b)
mr_temp$se = as.numeric(mr_temp$se)
mr_temp$pval = as.numeric(mr_temp$pval)
mr_temp$percent = round((exp(mr_temp$b) - 1) * 100, 2)
mr_temp$confint_lower = round((exp(mr_temp$b - 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint_upper = round((exp(mr_temp$b + 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint = paste0("[", mr_temp$confint_lower, ", ", mr_temp$confint_upper, "]")
fwrite(mr_temp, paste0(path, "mr/", today, "/mvmr_res_total.csv"))
	
##################################################
########## PLOT MENDELIAN RANDOMIZATION	
##################################################

path = "C:/Jiwoo_Lee/Healthcare_Cost_Research_2021/"
today = "20220505"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}
if(!require(ggplot2)) {install_github("ggplot2"); library(ggplot2)}
if(!require(gridExtra)) {install_github("gridExtra"); library(gridExtra)}

phenos = fread(paste0(path, "phenotypes.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
phenos = phenos[!which(phenos$notes == "non_ukbb" | phenos$notes == "statin_adjusted"),]

phenotype_raw = phenos$my_id[grepl("raw", phenos$my_id)]
phenotype_norm = phenos$my_id[grepl("irnt", phenos$my_id)]

total = fread(paste0(path, "MR/", today, "/fi_mr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE)
total_new = merge(phenos, total, by.x = "my_id", by.y = "id.exposure")

threshold = 15 # 17

total_sig = total_new[which(total_new$my_id %in% phenotype_norm & total_new$pval < 0.05/threshold & total_new$outcome == "log_total_cost_per_person_year"),]

# MAKE FOREST PLOT FROM ONE OUTCOME AND MANY EXPOSURES
temp = total_new[which(total_new$my_id %in% phenotype_norm & total_new$outcome == "log_total_cost_per_person_year" & total_new$method == "Inverse variance weighted"),]
temp$star = ifelse(temp$pval < 0.05/nrow(temp), "*", "")
temp$pound = ifelse(temp$my_id %in% names(table(total_sig$my_id)[which(table(total_sig$my_id) >= 3)]), "#", "")
temp$label = paste0(temp$percent, " ", temp$confint, " ", temp$star, " ", temp$pound)

#temp = temp[which(temp$phenotype %in% c("Serum creatinine", "Serum cystatin C", "Creatinine", "Cystatin C")),]

ggplot() +
	geom_hline(mapping = aes(yintercept = 0), size = 1, lty = "dashed") +
	geom_point(data = temp, mapping = aes(x = reorder(phenotype, -percent), y = percent, color = ifelse(pval < 0.05/nrow(temp), "Yes", "No")), size = 3, pch = 15) +
	coord_flip () +
	geom_errorbar(data = temp, mapping = aes(x = reorder(phenotype, -percent), ymin = confint_lower, ymax = confint_upper, color = ifelse(pval < 0.05/nrow(temp), "Yes", "No")), width = 0, size = 1, alpha = 0.5) +
	geom_text(data = temp, mapping = aes(x = reorder(phenotype, -percent), y = percent, label = label), vjust = -1, size = 5) +
	scale_color_manual(values = c("gray", "black")) +
	scale_x_discrete(limits = temp$phenotype[order(-temp$percent)], labels = temp$phenotype[order(-temp$percent)]) +
	scale_y_continuous(breaks = seq(-5, 25, 5)) + 
	labs(x = "", y = "Percent Change in Healthcare Cost per\n1 SD Increase in Risk Factor", color = "") + 
	guides(color = FALSE) +
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 25, angle = 45, vjust = 0.5, color = "black"), axis.title.x = element_text(size = 25, color = "black"), axis.text.y = element_text(size = 25, color = "black"), axis.title.y = element_text(size = 25, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))
ggsave(paste0(path, "MR/", today, "/main_mr_res_new.jpg"), device = "jpg", type = "cairo", height = 10, width = 16, family = "Arial")

# MAKE FOREST PLOT FROM ONE EXPOSURE AND MANY OUTCOMES (CARE AND SEX SENSITIVITY ANALYSIS)
key = total_new[which(total_new$my_id %in% phenotype_norm & total_new$outcome == "log_total_cost_per_person_year" & total_new$method == "Inverse variance weighted" & total_new$pval < 0.05/threshold),]
key = key[order(key$percent),]
key$number = 1:nrow(key)
key = key[,c("exposure", "number")]

outcome_care = c("log_avohilmo_cost_per_person_year", "log_hilmo_cost_per_person_year", "log_kela_cost_per_person_year")
care = total_new[which(total_new$outcome %in% outcome_care & total_new$my_id %in% phenotype_norm & total_new$method == "Inverse variance weighted"),]
care = merge(care, key, by = "exposure")
#care = merge(phenos, care, by.x = "my_id", by.y = "exposure")
care$number = 10*care$number
care$sub_number = ifelse(care$outcome == "log_avohilmo_cost_per_person_year", care$number - 2,
	ifelse(care$outcome == "log_hilmo_cost_per_person_year", care$number + 2, care$number))
care$label = ifelse(care$sub_number %% 10 == 0, care$phenotype, " ")
care = as.data.frame(care)

res = NULL
for (i in 1:length(unique(care$phenotype))) {
	care_temp = care[which(care$phenotype == unique(care$phenotype)[i]),]
	b_kela = care_temp$b[which(care_temp$outcome == "log_kela_cost_per_person_year")]
	b_hilmo = care_temp$b[which(care_temp$outcome == "log_hilmo_cost_per_person_year")]
	b_avohilmo = care_temp$b[which(care_temp$outcome == "log_avohilmo_cost_per_person_year")]
	se_kela = care_temp$se[which(care_temp$outcome == "log_kela_cost_per_person_year")]
	se_hilmo = care_temp$se[which(care_temp$outcome == "log_hilmo_cost_per_person_year")]
	se_avohilmo = care_temp$se[which(care_temp$outcome == "log_avohilmo_cost_per_person_year")]
	z_hilmo_kela = (b_hilmo - b_kela)/sqrt((se_hilmo)^2+(se_kela)^2)
	z_avohilmo_kela = (b_avohilmo - b_kela)/sqrt((se_avohilmo)^2+(se_kela)^2)
	z_hilmo_avohilmo = (b_hilmo - b_avohilmo)/sqrt((se_hilmo)^2+(se_avohilmo)^2)
	res = rbind(res, c(unique(care$phenotype)[i], b_kela, b_hilmo, b_avohilmo, se_kela, se_hilmo, se_avohilmo, z_hilmo_kela, z_avohilmo_kela, z_hilmo_avohilmo, pnorm(z_hilmo_kela), pnorm(z_avohilmo_kela), pnorm(z_hilmo_avohilmo)))
}
res = as.data.frame(res)
colnames(res) = c("phenotype", "b_kela", "b_hilmo", "b_avohilmo", "se_kela", "se_hilmo", "se_avohilmo", "z_hilmo_kela", "z_avohilmo_kela", "z_hilmo_avohilmo", "p_hilmo_kela", "p_avohilmo_kela", "p_hilmo_avohilmo")
res$b_kela = as.numeric(res$b_kela)
res$b_hilmo = as.numeric(res$b_hilmo)
res$b_avohilmo = as.numeric(res$b_avohilmo)
res$se_kela = as.numeric(res$se_kela)
res$se_hilmo = as.numeric(res$se_hilmo)
res$se_avohilmo = as.numeric(res$se_avohilmo)
res$z_hilmo_kela = as.numeric(res$z_hilmo_kela)
res$z_avohilmo_kela = as.numeric(res$z_avohilmo_kela)
res$z_hilmo_avohilmo = as.numeric(res$z_hilmo_avohilmo)
res$p_hilmo_kela = as.numeric(res$p_hilmo_kela)
res$p_avohilmo_kela = as.numeric(res$p_avohilmo_kela)
res$p_hilmo_avohilmo = as.numeric(res$p_hilmo_avohilmo)
care$star_hilmo_kela = ifelse(care$label %in% res$phenotype[which(res$p_hilmo_kela < 0.05 /nrow(res))], "*", "")
care$star_avohilmo_kela = ifelse(care$label %in% res$phenotype[which(res$p_avohilmo_kela < 0.05 /nrow(res))], "*", "")
care$star_hilmo_avohilmo = ifelse(care$label %in% res$phenotype[which(res$p_hilmo_avohilmo < 0.05 /nrow(res))], "*", "")

A = ggplot() +
	geom_hline(mapping = aes(yintercept = 0), size = 1, lty = "dashed") +
	geom_point(data = care, mapping = aes(x = -sub_number, y = percent, color = outcome, alpha = ifelse(pval < 0.05/threshold, 1, 0)), size = 3, pch = 15) +
	coord_flip() +
	geom_errorbar(data = care, mapping = aes(x = -sub_number, ymin = confint_lower, ymax = confint_upper, color = outcome, alpha = ifelse(pval < 0.05/threshold, 0.7, 0)), width = 0, size = 1) +
	geom_segment(size = 1, mapping = aes(x = -(care$sub_number[which(care$star_hilmo_kela == "*")]+0.25), xend = -(care$sub_number[which(care$star_hilmo_kela == "*")]+1.75), y = care$confint_upper[which(care$star_hilmo_kela == "*")] + 3, yend = care$confint_upper[which(care$star_hilmo_kela == "*")] + 3)) +
	geom_text(mapping = aes(x = -(care$sub_number[which(care$star_hilmo_kela == "*")]+1.5), y = care$confint_upper[which(care$star_hilmo_kela == "*")] + 5, label = care$star_hilmo_kela[which(care$star_hilmo_kela == "*")]), size = 10) +
	geom_segment(size = 1, mapping = aes(x = -(care$sub_number[which(care$star_avohilmo_kela == "*")]-0.25), xend = -(care$sub_number[which(care$star_avohilmo_kela == "*")]-1.75), y = care$confint_upper[which(care$star_avohilmo_kela == "*")] + 3, yend = care$confint_upper[which(care$star_avohilmo_kela == "*")] + 3)) +
	geom_text(mapping = aes(x = -(care$sub_number[which(care$star_avohilmo_kela == "*")]-0.5), y = care$confint_upper[which(care$star_avohilmo_kela == "*")] + 5, label = care$star_avohilmo_kela[which(care$star_avohilmo_kela == "*")]), size = 10) +
	scale_x_continuous(breaks = -care$sub_number, labels = care$label) +
	scale_color_manual(limits = c("log_avohilmo_cost_per_person_year", "log_hilmo_cost_per_person_year", "log_kela_cost_per_person_year"), labels = c("Primary Care", "Secondary Care", "Medication"), values = c("skyblue", "dodgerblue", "navyblue")) +
	scale_alpha_continuous(range = c(0.3, 1)) +
	labs(x = "", y = "% Change in Healthcare Costs per\n1 SD Increase in Risk Factor", color = "", tag = "A") +
	guides(alpha = FALSE) + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 25, angle = 45, vjust = 0.5, color = "black"), axis.title.x = element_text(size = 25, color = "black"), axis.text.y = element_text(size = 25, color = "black"), axis.title.y = element_text(size = 25, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))
#ggsave(paste0(path, "MR/", today, "/", "sensitivity_care.jpg"), device = "jpg", type = "cairo", height = 10, width = 10, family = "Arial")

outcome_sex = c("log_total_cost_per_female_year", "log_total_cost_per_male_year")
sex = total_new[which(total_new$outcome %in% outcome_sex & total_new$my_id %in% phenotype_norm & total_new$method == "Inverse variance weighted"),]
sex = merge(sex, key, by = "exposure")
#sex = merge(phenos, sex[which(sex$exposure %in% phenotype_norm)], by.x = "my_id", by.y = "exposure")
sex$number = 2*sex$number
sex$sub_number = ifelse(sex$outcome == "log_total_cost_per_female_year", sex$number - 0.5, sex$number)
sex$label = ifelse(sex$sub_number %% 1 == 0, sex$phenotype, " ")
sex = as.data.frame(sex)

res = NULL
for (i in 1:length(unique(sex$phenotype))) {
	sex_temp = sex[which(sex$phenotype == unique(sex$phenotype)[i]),]
	b_female = sex_temp$b[which(sex_temp$outcome == "log_total_cost_per_female_year")]
	b_male = sex_temp$b[which(sex_temp$outcome == "log_total_cost_per_male_year")]
	se_female = sex_temp$se[which(sex_temp$outcome == "log_total_cost_per_female_year")]
	se_male = sex_temp$se[which(sex_temp$outcome == "log_total_cost_per_male_year")]
	z = (b_female - b_male)/sqrt((se_female)^2+(se_male)^2)
	res = rbind(res, c(unique(sex$phenotype)[i], b_female, b_male, se_female, se_male, z, pnorm(z)))
}
res = as.data.frame(res)
colnames(res) = c("phenotype", "b_female", "b_male", "se_female", "se_male", "z", "p")
res$b_female = as.numeric(res$b_female)
res$b_male = as.numeric(res$b_male)
res$se_female = as.numeric(res$se_female)
res$se_male = as.numeric(res$se_male)
res$z = as.numeric(res$z)
res$p = as.numeric(res$p)
sex$star = ifelse(sex$label %in% res$phenotype[which(res$p < 0.05/nrow(res))], "*", "")

B = ggplot() +
	geom_hline(mapping = aes(yintercept = 0), size = 1, lty = "dashed") +
	geom_point(data = sex, mapping = aes(x = -sub_number, y = percent, color = outcome, alpha = ifelse(pval < 0.05/threshold, 1, 0)), size = 3, pch = 15) +
	coord_flip() +
	geom_errorbar(data = sex, mapping = aes(x = -sub_number, ymin = confint_lower, ymax = confint_upper, color = outcome, alpha = ifelse(pval < 0.05/threshold, 0.7, 0)), width = 0, size = 1) +
	geom_segment(size = 1, mapping = aes(x = -(sex$sub_number[which(sex$star == "*")]-0.5), xend = -sex$sub_number[which(sex$star == "*")], y = sex$confint_upper[which(sex$star == "*")] + 1, yend = sex$confint_upper[which(sex$star == "*")] + 1)) +
	geom_text(mapping = aes(x = -(sex$sub_number[which(sex$star == "*")]-0.05) , y = sex$confint_upper[which(sex$star == "*")] + 3, label = sex$star[which(sex$star == "*")]), size = 10) +
	scale_x_continuous(breaks = -sex$sub_number, labels = sex$label) +
	scale_color_manual(limits = c("log_total_cost_per_female_year", "log_total_cost_per_male_year"), labels = c("Female", "Male"), values = c("dark green", "yellowgreen")) +
	scale_alpha_continuous(range = c(0.3, 1)) + 
	labs(x = "", y = "% Change in Healthcare Costs per\n1 SD Increase in Risk Factor", color = "", tag = "B") +
	guides(alpha = FALSE) + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 25, angle = 45, vjust = 0.5, color = "black"), axis.title.x = element_text(size = 25, color = "black"), axis.text.y = element_text(size = 25, color = "black"), axis.title.y = element_text(size = 25, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))
#ggsave(paste0(path, "MR/", today, "/", "sensitivity_sex.jpg"), device = "jpg", type = "cairo", height = 10, width = 10, family = "Arial")

outcome_age = c("log_total_cost_per_person_year_0_30", "log_total_cost_per_person_year_30_60", "log_total_cost_per_person_year_60_90")
age = total_new[which(total_new$outcome %in% outcome_age & total_new$my_id %in% phenotype_norm & total_new$method == "Inverse variance weighted"),]
age = merge(age, key, by = "exposure")
age$number = 10*age$number
age$sub_number = ifelse(age$outcome == "log_total_cost_per_person_year_0_30", age$number - 2,
	ifelse(age$outcome == "log_total_cost_per_person_year_60_90", age$number + 2, age$number))
age$label = ifelse(age$sub_number %% 10 == 0, age$phenotype, " ")
age = as.data.frame(age)

res = NULL
for (i in 1:length(unique(age$phenotype))) {
	age_temp = age[which(age$phenotype == unique(age$phenotype)[i]),]
	b_young = age_temp$b[which(age_temp$outcome == "log_total_cost_per_person_year_0_30")]
	b_mid = age_temp$b[which(age_temp$outcome == "log_total_cost_per_person_year_30_60")]
	b_old = age_temp$b[which(age_temp$outcome == "log_total_cost_per_person_year_60_90")]
	se_young = age_temp$se[which(age_temp$outcome == "log_total_cost_per_person_year_0_30")]
	se_mid = age_temp$se[which(age_temp$outcome == "log_total_cost_per_person_year_30_60")]
	se_old = age_temp$se[which(age_temp$outcome == "log_total_cost_per_person_year_60_90")]
	z_young_mid = (b_young - b_mid)/sqrt((se_young)^2+(se_mid)^2)
	z_mid_old = (b_mid - b_old)/sqrt((se_mid)^2+(se_old)^2)
	z_young_old = (b_young - b_old)/sqrt((se_young)^2+(se_old)^2)
	res = rbind(res, c(unique(age$phenotype)[i], b_young, b_mid, b_old, se_young, se_mid, se_old, z_young_mid, z_mid_old, z_young_old, pnorm(z_young_mid), pnorm(z_mid_old), pnorm(z_young_old)))
}
res = as.data.frame(res)
colnames(res) = c("phenotype", "b_young", "b_mid", "b_old", "se_young", "se_mid", "se_old", "z_young_mid", "z_mid_old", "z_young_old", "p_yound_mid", "p_mid_old", "p_young_old")
res$b_young = as.numeric(res$b_young)
res$b_mid = as.numeric(res$b_mid)
res$b_old = as.numeric(res$b_old)
res$se_young = as.numeric(res$se_young)
res$se_mid = as.numeric(res$se_mid)
res$se_old = as.numeric(res$se_old)
res$z_young_mid = as.numeric(res$z_young_mid)
res$z_mid_old = as.numeric(res$z_mid_old)
res$z_young_old = as.numeric(res$z_young_old)
res$p_yound_mid = as.numeric(res$p_yound_mid)
res$p_mid_old = as.numeric(res$p_mid_old)
res$p_young_old = as.numeric(res$p_young_old)
age$star_young_mid = ifelse(age$label %in% res$phenotype[which(res$p_young_mid < 0.05 /nrow(res))], "*", "")
age$star_mid_old = ifelse(age$label %in% res$phenotype[which(res$p_mid_old < 0.05 /nrow(res))], "*", "")
age$star_young_old = ifelse(age$label %in% res$phenotype[which(res$p_young_old < 0.05 /nrow(res))], "*", "")

C = ggplot() +
	geom_hline(mapping = aes(yintercept = 0), size = 1, lty = "dashed") +
	geom_point(data = age, mapping = aes(x = -sub_number, y = percent, color = outcome, alpha = ifelse(pval < 0.05/threshold, 1, 0)), size = 3, pch = 15) +
	coord_flip() +
	geom_errorbar(data = age, mapping = aes(x = -sub_number, ymin = confint_lower, ymax = confint_upper, color = outcome, alpha = ifelse(pval < 0.05/threshold, 0.7, 0)), width = 0, size = 1) +
	geom_segment(size = 1, mapping = aes(x = -(age$sub_number[which(age$star_mid_old == "*")]+0.25), xend = -(age$sub_number[which(age$star_mid_old == "*")]+1.75), y = age$confint_upper[which(age$star_mid_old == "*")] + 13, yend = age$confint_upper[which(age$star_mid_old == "*")] + 13)) +
	geom_text(mapping = aes(x = -(age$sub_number[which(age$star_mid_old == "*")]+1.75), y = age$confint_upper[which(age$star_mid_old == "*")] + 15, label = age$star_mid_old[which(age$star_mid_old == "*")]), size = 10) +
	geom_segment(size = 1, mapping = aes(x = -(age$sub_number[which(age$star_young_mid == "*")]-0.25), xend = -(age$sub_number[which(age$star_young_mid == "*")]-1.75), y = age$confint_upper[which(age$star_young_mid == "*")] + 13, yend = age$confint_upper[which(age$star_young_mid == "*")] + 13)) +
	geom_text(mapping = aes(x = -(age$sub_number[which(age$star_young_mid == "*")]-0.25), y = age$confint_upper[which(age$star_young_mid == "*")] + 15, label = age$star_young_mid[which(age$star_young_mid == "*")]), size = 10) +
	scale_x_continuous(breaks = -age$sub_number, labels = age$label) +
	scale_color_manual(limits = c("log_total_cost_per_person_year_0_30", "log_total_cost_per_person_year_30_60", "log_total_cost_per_person_year_60_90"), labels = c("<30 years old", "30-60 years old", ">60 years old"), values = c("coral", "red", "red4")) +
	scale_alpha_continuous(range = c(0.3, 1)) +
	labs(x = "", y = "% Change in Healthcare Costs per\n1 SD Increase in Risk Factor", color = "", tag = "C") +
	guides(alpha = FALSE) + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 25, angle = 45, vjust = 0.5, color = "black"), axis.title.x = element_text(size = 25, color = "black"), axis.text.y = element_text(size = 25, color = "black"), axis.title.y = element_text(size = 25, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))
#ggsave(paste0(path, "MR/", today, "/", "sensitivity_age.jpg"), device = "jpg", type = "cairo", height = 10, width = 10, family = "Arial")

G = arrangeGrob(A, B, C, nrow = 1)
ggsave(paste0(path, "MR/", today, "/sensitivity.jpg"), G, device = "jpg", type = "cairo", height = 10, width = 30, family = "Arial")

##################################################
########## PLOT MULTIVARIABLE MENDELIAN RANDOMIZATION	
##################################################

path = "C:/Jiwoo Lee/Healthcare_Cost_Research_2021/"
today = "20220505"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}
if(!require(ggplot2)) {install_github("ggplot2"); library(ggplot2)}
if(!require(BSDA)) {install_github("BSDA"); library(BSDA)}

exposure_list = c("ukb-a-360", "ieu-a-835", "ukb-a-382")

mvmr = fread(paste0(path, "MR/", today, "/mvmr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE)
mvmr = mvmr[which(mvmr$id.exposure %in% exposure_list),]

mvmr$number = c(seq(1, 2, length.out = 8), 
	seq(3, 4, length.out = 8),
	seq(5, 6, length.out = 8))
mvmr$color = as.factor(c(rep(seq(1, 2, length.out = 8), 3)))

threshold = nrow(mvmr)

temp = mvmr[which(mvmr$id.exposure == "ukb-a-360"),]
temp = temp[!c(nrow(temp)-1, nrow(temp)),]
c = -(temp$number[1]+temp$number[2])/2
C = ggplot() +
	geom_hline(mapping = aes(yintercept = 0), size = 1, lty = "dashed") +
	geom_vline(mapping = aes(xintercept = c), size = 1, lty = "dashed") + 
	geom_errorbar(data = temp, mapping = aes(x = -number, ymin = confint_lower, ymax = confint_upper, color = color), width = 0, size = 1) +
	geom_point(data = temp, mapping = aes(x = -number, y = percent, color = color, alpha = ifelse(pval < 0.05/threshold, 1, 0)), size = 3, pch = 15) +
	coord_flip() +
	scale_x_continuous(breaks = -temp$number, labels = c("", "", "", "", "", "")) +
	#scale_color_manual(limits = c("log_total_cost_per_person_year", "log_total_cost_per_person_year_cad_exclusion", "log_total_cost_per_person_year_statin_exclusion", "log_total_cost_per_person_year_htn_med_exclusion"), labels = c("All individuals", "Individuals without CAD", "Individuals without statin medication", "Individuals without HTN medication"), values = c("purple", "mediumpurple", "darkslateblue", "black")) +
	scale_alpha_continuous(range = c(0.3, 1)) +
	labs(x = "", y = "% Change in Healthcare Cost per\n1 SD Increase in SBP", color = "", tag = "C") +
	guides(alpha = FALSE, color = FALSE) + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 25, angle = 45, vjust = 0.5, color = "black"), axis.title.x = element_text(size = 25, color = "black"), axis.text.y = element_text(size = 25, color = "black"), axis.title.y = element_text(size = 25, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))
#ggsave(paste0(path, "MR/", today, "/", "mvmr_sbp.jpg"), device = "jpg", type = "cairo", height = 10, width = 10, family = "Arial")

temp = mvmr[which(mvmr$id.exposure == "ieu-a-835"),]
temp = temp[!c(nrow(temp)-1, nrow(temp)),]
b = -(temp$number[1]+temp$number[2])/2
B = ggplot() +
	geom_hline(mapping = aes(yintercept = 0), size = 1, lty = "dashed") +
	geom_vline(mapping = aes(xintercept = b), size = 1, lty = "dashed") + 
	geom_errorbar(data = temp, mapping = aes(x = -number, ymin = confint_lower, ymax = confint_upper, color = color), width = 0, size = 1) +
	geom_point(data = temp, mapping = aes(x = -number, y = percent, color = color, alpha = ifelse(pval < 0.05/threshold, 1, 0)), size = 3, pch = 15) +
	coord_flip() +
	scale_x_continuous(breaks = -temp$number, labels = c("", "", "", "", "", "")) +
	#scale_color_manual(limits = c("log_total_cost_per_person_year", "log_total_cost_per_person_year_cad_exclusion", "log_total_cost_per_person_year_statin_exclusion", "log_total_cost_per_person_year_htn_med_exclusion"), labels = c("All individuals", "Individuals without CAD", "Individuals without statin medication", "Individuals without HTN medication"), values = c("purple", "mediumpurple", "darkslateblue", "black")) +
	scale_alpha_continuous(range = c(0.3, 1)) +
	labs(x = "", y = "Percent Change in Healthcare Cost per\n1 SD Increase in BMI", color = "", tag = "B") +
	guides(alpha = FALSE, color = FALSE) + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 25, angle = 45, vjust = 0.5, color = "black"), axis.title.x = element_text(size = 25, color = "black"), axis.text.y = element_text(size = 25, color = "black"), axis.title.y = element_text(size = 25, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))
#ggsave(paste0(path, "MR/", today, "/", "mvmr_bmi.jpg"), device = "jpg", type = "cairo", height = 10, width = 10, family = "Arial")

temp = mvmr[which(mvmr$id.exposure == "ukb-a-382"),]
temp = temp[!c(nrow(temp)-1, nrow(temp)),]
a = -(temp$number[1]+temp$number[2])/2
A = ggplot() +
	geom_hline(mapping = aes(yintercept = 0), size = 1, lty = "dashed") +
	geom_vline(mapping = aes(xintercept = a), size = 1, lty = "dashed") + 
	geom_errorbar(data = temp, mapping = aes(x = -number, ymin = confint_lower, ymax = confint_upper, color = color), width = 0, size = 1) +
	geom_point(data = temp, mapping = aes(x = -number, y = percent, color = color, alpha = ifelse(pval < 0.05/threshold, 1, 0)), size = 3, pch = 15) +
	coord_flip() +
	scale_x_continuous(breaks = -temp$number, labels = c("MAIN EFFECT", "Back pain", "Chronic ischemic\nheart disease", "Type 2 diabetes", "Chronic obstructive\npulmonary disease", "Stroke")) +
	#scale_color_manual(limits = c("log_total_cost_per_person_year", "log_total_cost_per_person_year_cad_exclusion", "log_total_cost_per_person_year_statin_exclusion", "log_total_cost_per_person_year_htn_med_exclusion"), labels = c("All individuals", "Individuals without CAD", "Individuals without statin medication", "Individuals without HTN medication"), values = c("purple", "mediumpurple", "darkslateblue", "black")) +
	scale_alpha_continuous(range = c(0.3, 1)) +
	labs(x = "Effect after adjustment for:", y = "Percent Change in Healthcare Cost per\n1 SD Increase in WC", color = "", tag = "A") +
	guides(alpha = FALSE, color = FALSE) + 
	theme(plot.tag = element_text(size = 25), legend.position = "bottom", axis.text.x = element_text(size = 25, angle = 45, vjust = 0.5, color = "black"), axis.title.x = element_text(size = 25, color = "black"), axis.text.y = element_text(size = 25, color = "black"), axis.title.y = element_text(size = 25, color = "black"), legend.text = element_text(size = 15), legend.title = element_text(size = 25), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))
#ggsave(paste0(path, "MR/", today, "/", "mvmr_wc.jpg"), device = "jpg", type = "cairo", height = 10, width = 10, family = "Arial")

total = fread(paste0(path, "MR/", today, "/fi_mr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE)
total[which(total$outcome == "log_total_cost_per_person_year" & total$exposure %in% c("Waist circumference", "Body mass index", "Systolic blood pressure") & total$method == "Inverse variance weighted"),]
total[which(total$outcome == "log_total_cost_per_person_year" & total$method == "Inverse variance weighted"),]
total[which(total$id.exposure == "ukb-a-382_irnt" & total$outcome == "log_total_cost_per_person_year"),]
total[which(total$id.exposure == "ieu-a-835_irnt" & total$outcome == "log_total_cost_per_person_year"),]
total[which(total$id.exposure == "ukb-a-360_irnt" & total$outcome == "log_total_cost_per_person_year"),]

total = fread(paste0(path, "MR/", today, "/uk_mr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE)
total[which(total$exposure %in% c("Waist circumference", "Body mass index", "Systolic blood pressure") & total$method == "Inverse variance weighted"),]

total = fread(paste0(path, "MR/", today, "/nl_mr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE)
total[which(total$exposure %in% c("Waist circumference", "Body mass index", "Systolic blood pressure") & total$method == "Inverse variance weighted"),]

G = arrangeGrob(A, B, C, nrow = 1)
ggsave(paste0(path, "MR/", today, "/mvmr_new.jpg"), G, device = "jpg", type = "cairo", height = 10, width = 30)

##################################################
########## ANALYZE POLYGENIC RISK SCORE
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}

total_new = fread(file = paste0(path, "healthcare_cost_20220505.txt"), header = TRUE, stringsAsFactors = FALSE)

gp = fread(file = paste0(path, "prs/finngen_R9/finngen_R9_netherlands_gp.fastGWA.sscore"), header = TRUE, stringsAsFactors = FALSE)
gp_new = gp[,c("IID", "SCORE1_AVG")]
colnames(gp_new) = c("finngenid", "neth_pri_prs")
hospital = fread(file = paste0(path, "prs/finngen_R9/finngen_R9_netherlands_hospital.fastGWA.sscore"), header = TRUE, stringsAsFactors = FALSE)
hospital_new = hospital[,c("IID", "SCORE1_AVG")]
colnames(hospital_new) = c("finngenid", "neth_sec_prs")
pharmacy = fread(file = paste0(path, "prs/finngen_R9/finngen_R9_netherlands_pharmacy.fastGWA.sscore"), header = TRUE, stringsAsFactors = FALSE)
pharmacy_new = pharmacy[,c("IID", "SCORE1_AVG")]
colnames(pharmacy_new) = c("finngenid", "neth_drug_prs")
neth = fread(file = paste0(path, "prs/finngen_R9/finngen_R9_netherlands_total.fastGWA.sscore"), header = TRUE, stringsAsFactors = FALSE)
neth_new = neth[,c("IID", "SCORE1_AVG")]
colnames(neth_new) = c("finngenid", "neth_total_prs")
ukbb = fread(file = paste0(path, "prs/finngen_R9/ukbb_secondary_care_costs_prs_cs.sscore"), header = TRUE, stringsAsFactors = FALSE)
ukbb_new = ukbb[,c("IID", "SCORE1_AVG")]
colnames(ukbb_new) = c("finngenid", "ukbb_sec_prs")

total_new = merge(total_new, gp_new, by = "finngenid")
total_new = merge(total_new, hospital_new, by = "finngenid")
total_new = merge(total_new, pharmacy_new, by = "finngenid")
total_new = merge(total_new, neth_new, by = "finngenid")
total_new = merge(total_new, ukbb_new, by = "finngenid")

#total_temp = total_new
#total_temp = total_temp %>% mutate(gp_decile = ntile(gp_prs, 10))
#total_temp = total_temp %>% mutate(hospital_decile = ntile(hospital_prs, 10))
#total_temp = total_temp %>% mutate(pharmacy_decile = ntile(pharmacy_prs, 10))
#total_temp = total_temp %>% mutate(total_decile = ntile(total_prs, 10))

ukbb_sec_mod = lm(hilmo_cost ~ scale(ukbb_sec_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
ukbb_sec_mod_year = lm(hilmo_cost_per_person_year ~ scale(ukbb_sec_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
ukbb_sec_mod_log = lm(log_hilmo_cost_per_person_year ~ scale(ukbb_sec_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)

neth_sec_mod = lm(hilmo_cost ~ scale(neth_sec_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
neth_sec_mod_year = lm(hilmo_cost_per_person_year ~ scale(neth_sec_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
neth_sec_mod_log = lm(log_hilmo_cost_per_person_year ~ scale(neth_sec_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)

neth_total_mod = lm(total_cost ~ scale(neth_total_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
neth_total_mod_year = lm(total_cost_per_person_year ~ scale(neth_total_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
neth_total_mod_log = lm(log_total_cost_per_person_year ~ scale(neth_total_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)

neth_pri_mod = lm(avohilmo_cost ~ scale(neth_pri_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
neth_pri_mod_year = lm(avohilmo_cost_per_person_year ~ scale(neth_pri_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
neth_pri_mod_log = lm(log_avohilmo_cost_per_person_year ~ scale(neth_pri_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)

neth_drug_mod = lm(kela_cost ~ scale(neth_drug_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
neth_drug_mod_year = lm(kela_cost_per_person_year ~ scale(neth_drug_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)
neth_drug_mod_log = lm(log_kela_cost_per_person_year ~ scale(neth_drug_prs) + sex + birth_year + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_new)

res_total = as.data.frame(rbind(c("UK Secondary", coef(summary(ukbb_sec_mod))[2,]),
	c("NL Total", coef(summary(neth_total_mod))[2,]),
	c("NL Primary", coef(summary(neth_pri_mod))[2,]),
	c("NL Secondary", coef(summary(neth_sec_mod))[2,]),
	c("NL Drug", coef(summary(neth_drug_mod))[2,])
	))
colnames(res_total) = c("model", "b", "se", "t", "p")
res_total$b = as.numeric(res_total$b)
res_total$se = as.numeric(res_total$se)
res_total$p = as.numeric(res_total$p)

ggplot() +
	coord_flip() +
	geom_point(data = res_total, mapping = aes(x = model, y = b), size = 3) +
	geom_errorbar(data = res_total, mapping = aes(x = model, ymin = b-1.96*se, ymax = b+1.96*se), size = 1, width = 0.1) +
	geom_label(data = res_total, mapping = aes(x = model, y = b, label = paste("P = ", signif(p, digits = 3))), nudge_x = 0.25, size = 5) +
	labs(x = "PRS", y = "Increase in Total Healthcare Costs Associated with 1 SD Increase in PRS") +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

res_year = as.data.frame(rbind(c("UK Secondary", coef(summary(ukbb_sec_mod_year))[2,]),
	c("NL Total", coef(summary(neth_total_mod_year))[2,]),
	c("NL Primary", coef(summary(neth_pri_mod_year))[2,]),
	c("NL Secondary", coef(summary(neth_sec_mod_year))[2,]),
	c("NL Drug", coef(summary(neth_drug_mod_year))[2,])
	))
colnames(res_year) = c("model", "b", "se", "t", "p")
res_year$b = as.numeric(res_year$b)
res_year$se = as.numeric(res_year$se)
res_year$p = as.numeric(res_year$p)

ggplot() +
	coord_flip() +
	geom_point(data = res_year, mapping = aes(x = model, y = b), size = 3) +
	geom_errorbar(data = res_year, mapping = aes(x = model, ymin = b-1.96*se, ymax = b+1.96*se), size = 1, width = 0.1) +
	geom_label(data = res_year, mapping = aes(x = model, y = b, label = paste("P = ", signif(p, digits = 3))), nudge_x = 0.25, size = 5) +
	labs(x = "PRS", y = "Increase in Annual Healthcare Costs Associated with 1 SD Increase in PRS") +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

res_log = as.data.frame(rbind(c("UK Secondary", coef(summary(ukbb_sec_mod_log))[2,]),
	c("NL Total", coef(summary(neth_total_mod_log))[2,]),
	c("NL Primary", coef(summary(neth_pri_mod_log))[2,]),
	c("NL Secondary", coef(summary(neth_sec_mod_log))[2,]),
	c("NL Drug", coef(summary(neth_drug_mod_log))[2,])
	))
colnames(res_log) = c("model", "b", "se", "t", "p")
res_log$b = as.numeric(res_log$b)
res_log$se = as.numeric(res_log$se)
res_log$p = as.numeric(res_log$p)

res_log$beta = (exp(res_log$b)-1)*100
res_log$beta_lower = (exp(res_log$b - 1.96 * res_log$se)-1)*100
res_log$beta_upper = (exp(res_log$b + 1.96 * res_log$se)-1)*100

ggplot() +
	coord_flip() +
	geom_point(data = res_log, mapping = aes(x = model, y = beta), size = 3) +
	geom_errorbar(data = res_log, mapping = aes(x = model, ymin = beta_lower, ymax = beta_upper), size = 1, width = 0.1) +
	geom_label(data = res_log, mapping = aes(x = model, y = beta, label = paste("P = ", signif(p, digits = 3))), nudge_x = 0.25, size = 5) +
	labs(x = "Model", y = "Percent Change in Annual Healthcare Costs Associated with 1 SD Increase in PRS") +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

#ggplot() + 
#	geom_violin(data = total_temp, mapping = aes(x = as.factor(total_decile), y = log_total_cost_per_person_year))

#cor.test(total_new$log_kela_cost_per_person_year, total_new$pharmacy_prs)

##################################################
########## DO LINEAR REGRESSION WITH WEIGHTS
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}

bmi = extract_instruments(outcomes = "ieu-a-835")
sbp = extract_instruments(outcomes = "ukb-a-360")
wc = extract_instruments(outcomes = "ukb-a-382")
total = rbind(rbind(bmi, sbp), wc)

sumstat = fread(file = paste0(path, "gwas/20220505/log_total_cost_per_person_year_final.txt"), header = TRUE, stringsAsFactors = FALSE)
sumstat_new = sumstat[which(sumstat$rsid %in% total$SNP),]
sumstat_new$cpra = paste0("chr", sumstat_new$chr, "_", sumstat_new$pos, "_", sumstat_new$ref, "_", sumstat_new$alt)
fwrite(list(sumstat_new$cpra), paste0(path, "bmi_sbp_wc_cpra_list.txt"), col.names = FALSE, row.names = FALSE)
fwrite(sumstat_new[,c("cpra", "alt")], paste0(path, "bmi_sbp_wc_cpra_alt_list.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")

plink --bfile /finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9 \
--extract bmi_sbp_wc_cpra_list.txt \
--recode-allele bmi_sbp_wc_cpra_alt_list.txt \
--recode A --threads 16 --out bmi_sbp_wc

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}

sumstat = fread(file = paste0(path, "log_total_cost_per_person_year"), header = TRUE, stringsAsFactors = FALSE)
colnames(sumstat)[1] = "chr"
sumstat$cpra = paste0("chr", sumstat$chr, "_", sumstat$pos, "_", sumstat$ref, "_", sumstat$alt)

total = fread(file = paste0(path, "pheno.txt"), header = TRUE, stringsAsFactors = FALSE)
variants = fread(file = paste0(path, "bmi_sbp_wc.raw"))
variants = variants[,-c("FID", "SEX")]
variant_list = colnames(variants)[5:length(colnames(variants))]
weights = fread(file = "/finngen/red/yllza/raked.weights.csv", header = TRUE, stringsAsFactors = FALSE)
weights = weights[,-"V1"]
total_new = merge(total, variants, by.x = "FID", by.y = "IID")
total_new = merge(total_new, weights, by.x = "FID", by.y = "FINNGENID")
total_temp = total_new

mod_res = NULL
for (i in 1:length(variant_list)) {
	#temp = cbind(variants$IID, variants[[variant_list[i]]])
	#colnames(temp) = c("finngenid", "variant")
	#total_new = merge(total, temp, by = "finngenid")
	total_temp$variant = as.numeric(as.character(total_temp[[variant_list[i]]]))
	mod_unweight = lm(formula = log_total_cost_per_person_year ~ variant + birth_year + birth_year_squared + SEX_IMPUTED + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_temp)
	coef_unweight = coef(summary(mod_unweight))["variant", "Estimate"]
	ste_unweight = coef(summary(mod_unweight))["variant", "Std. Error"]
	pval_unweight = coef(summary(mod_unweight))["variant", "Pr(>|t|)"]
	mod_weight = lm(formula = log_total_cost_per_person_year ~ variant + birth_year + birth_year_squared + SEX_IMPUTED + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_temp, weights = total_temp$raked.weight)
	coef_weight = coef(summary(mod_weight))["variant", "Estimate"]
	ste_weight = coef(summary(mod_weight))["variant", "Std. Error"]
	pval_weight = coef(summary(mod_weight))["variant", "Pr(>|t|)"]
	mod_res = rbind(mod_res, c(variant_list[i], coef_unweight, ste_unweight, pval_unweight, coef_weight, ste_weight, pval_weight))
	print(i)
}
mod_res = as.data.frame(mod_res)
colnames(mod_res) = c("variant", "beta_lr_unweight", "se_lr_unweight", "pval_lr_unweight", "beta_lr_weight", "se_lr_weight", "pval_lr_weight")
mod_res$beta_lr_unweight = as.numeric(mod_res$beta_lr_unweight)
mod_res$se_lr_unweight = as.numeric(mod_res$se_lr_unweight)
mod_res$pval_lr_unweight = as.numeric(mod_res$pval_lr_unweight)
mod_res$beta_lr_weight = as.numeric(mod_res$beta_lr_weight)
mod_res$se_lr_weight = as.numeric(mod_res$se_lr_weight)
mod_res$pval_lr_weight = as.numeric(mod_res$pval_lr_weight)
sumstat$variant = paste0(sumstat$cpra, "_", sumstat$alt)
comb = merge(mod_res, sumstat[,c("variant", "beta", "sebeta", "pval")], by = "variant")
fwrite(comb, paste0(path, "weighted_linear_regression_results.txt"), sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

par(mfrow=c(2,3))
plot(comb$beta, comb$beta_lr_unweight, xlab = "GWAS BETA", ylab = "UNWEIGHTED REGRESSION BETA")
abline(a = 0, b = 1)
text(0.005, 0.02, paste0("R^2: ", round(cor(comb$beta, comb$beta_lr_unweight)^2, 5)), col = "red")

plot(comb$beta, comb$beta_lr_weight, xlab = "GWAS BETA", ylab = "WEIGHTED REGRESSION BETA")
abline(a = 0, b = 1)
text(0.005, 0.02, paste0("R^2: ", round(cor(comb$beta, comb$beta_lr_weight)^2, 5)), col = "red")

plot(comb$beta_lr_unweight, comb$beta_lr_weight, xlab = "UNWEIGHTED REGRESSION BETA", ylab = "WEIGHTED REGRESSION BETA")
abline(a = 0, b = 1)
text(0.005, 0.02, paste0("R^2: ", round(cor(comb$beta_lr_unweight, comb$beta_lr_weight)^2, 5)), col = "red")

plot(comb$sebeta, comb$se_lr_unweight, xlab = "GWAS SE", ylab = "UNWEIGHTED REGRESSION SE")
abline(a = 0, b = 1)
text(0.015, 0.01, paste0("R^2: ", round(cor(comb$sebeta, comb$se_lr_unweight)^2, 5)), col = "red")

plot(comb$sebeta, comb$se_lr_weight, xlab = "GWAS SE", ylab = "WEIGHTED REGRESSION SE")
abline(a = 0, b = 1)
text(0.015, 0.01, paste0("R^2: ", round(cor(comb$sebeta, comb$se_lr_weight)^2, 5)), col = "red")

plot(comb$se_lr_unweight, comb$se_lr_weight, xlab = "UNWEIGHTED REGRESSION SE", ylab = "WEIGHTED REGRESSION SE")
abline(a = 0, b = 1)
text(0.015, 0.01, paste0("R^2: ", round(cor(comb$se_lr_unweight, comb$se_lr_weight)^2, 5)), col = "red")

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/
R

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"
total = fread(file = paste0(path, "gwas/20220505/weighted_linear_regression_results.txt"), header = TRUE, stringsAsFactors = FALSE)

total$chr = sapply(temp, "[", 1)
total$pos = sapply(temp, "[", 2)
total$ref = sapply(temp, "[", 3)
total$alt = sapply(temp, "[", 4)
total$cpra = paste0(total$chr, "_", total$pos, "_", total$ref, "_", total$alt)
total$chr = substring(total$chr, 4)
total = total[,c("chr", "pos", "cpra", "ref", "alt", "variant", 
	"beta_lr_unweight", "se_lr_unweight", "pval_lr_unweight",
	"beta_lr_weight", "se_lr_weight", "pval_lr_weight",
	"beta", "sebeta", "pval")]
fwrite(total, paste0(path, "gwas/20220505/weighted_linear_regression_results.txt.gz"), sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE, na = NA)

awk 'function nn(x){split(x,xx,"");str=n[xx[1]];for(i=2;i<=length(xx);i++){str=str""n[xx[i]]};return str}BEGIN{n["A"]=1;n["C"]=2;n["G"]=3;n["T"]=4}NR==FNR{split($3,a,"_");split(a[1],b,"hr");id=b[2]":"a[2]":"a[3]":"a[4];if(nn(a[4])<nn(a[3])){id=b[2]":"a[2]":"a[4]":"a[3]};cpra[id]++;next}($2 in cpra){print $1" "$2}' <(zcat weighted_linear_regression_results.txt.gz) <(zcat /medpop/esp2/jiwoolee/cost_rep/dbsnp_151_b38_rsid_cpra.txt.gz) > dbsnp_151_b38_rsid_cpra.txt

awk 'BEGIN{OFS=" "}NR==FNR{rsid[$2]=$1;next}FNR==1{print;next}{split($3,a,"_");split(a[1],b,"hr");cpra=b[2]":"a[2]":"a[3]":"a[4];cpar=b[2]":"a[2]":"a[4]":"a[3]}(cpra in rsid){$3=rsid[cpra];print;next}(cpar in rsid){$3=rsid[cpar];print;next}{print}' dbsnp_151_b38_rsid_cpra.txt <(zcat weighted_linear_regression_results.txt.gz) | gzip --best > weighted_linear_regression_results_final.txt.gz

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}

sumstat = fread(file = paste0(path, "gwas/20220505/weighted_linear_regression_results_final.txt.gz"), header = TRUE, stringsAsFactors = FALSE)

setDTthreads(0)

# List available GWAS
#ao <- available_outcomes()
#ao_df <- as.data.frame(ao)

phenos = fread(paste0(path, "phenotypes.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
phenos = phenos[which(phenos$my_id %in% c("ieu-a-835_irnt","ieu-a-835_raw","ukb-a-360_irnt","ukb-a-360_raw","ukb-a-382_irnt","ukb-a-382_raw")),]
exposure_id = phenos$gwas_id

mr_res = NULL
for (i in 1:length(exposure_id)) {
	print(exposure_id[i])
	# Get instruments
	print("Loading exposure data...")
	exposure_dat = extract_instruments(outcomes = exposure_id[i])
	exposure_dat$beta.exposure = exposure_dat$beta.exposure * phenos$transformation[i]
	exposure_dat$se.exposure = exposure_dat$se.exposure * phenos$transformation[i]
	# Get effects of instruments on outcome
	print("Loading outcome data...")
	outcome_temp = sumstat[which(sumstat$cpra %in% exposure_dat$SNP),]
	fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
	outcome_dat = read_outcome_data(snps = NULL,
		filename = paste0(path, "mr/", today, "/temp.csv"),
		sep = ",",
		snp_col = "cpra",
		beta_col = "beta_lr_weight",
	  se_col = "se_lr_weight",
	  effect_allele_col = "alt",
	  other_allele_col = "ref",
	  pval_col = "pval_lr_weight"
	)
	# Harmonise the exposure and outcome data
	print("Harmonizing data...")
	dat = harmonise_data(exposure_dat, outcome_dat)
	# Perform MR
	print("Performing MR...")
	res = mr(dat)
	if (phenos$transformation[i] != 1) {
		res$id.exposure = paste0(res$id.exposure, "_raw")
	}
	res$exposure = phenos$phenotype[i]
	res$outcome = "weight"
	mr_res = rbind(mr_res, res)
}
mr_temp = mr_res
mr_temp = as.data.frame(mr_temp)
mr_temp$b = as.numeric(mr_temp$b)
mr_temp$se = as.numeric(mr_temp$se)
mr_temp$pval = as.numeric(mr_temp$pval)
mr_temp$percent = round((exp(mr_temp$b) - 1) * 100, 2)
mr_temp$confint_lower = round((exp(mr_temp$b - 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint_upper = round((exp(mr_temp$b + 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint = paste0("[", mr_temp$confint_lower, ", ", mr_temp$confint_upper, "]")
mr_temp = mr_temp[order(mr_temp$pval),]
fwrite(mr_temp, paste0(path, "mr/", today, "/weight_mr_res_total.csv"))

##################################################
########## FIGURE OUT WHAT IS WRONG WITH LDL CHOLESTEROL
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}

ldl = extract_instruments(outcomes = "ieu-b-110")
sumstat = fread(file = paste0(path, "gwas/20220505/log_total_cost_per_person_year_final.txt.gz"), header = TRUE, stringsAsFactors = FALSE)
sumstat_new = sumstat[which(sumstat$cpra %in% ldl$SNP),]
sumstat_new$snp = sumstat_new$cpra
sumstat_new$cpra = paste0("chr", sumstat_new$chr, "_", sumstat_new$pos, "_", sumstat_new$ref, "_", sumstat_new$alt)
fwrite(list(sumstat_new$cpra), paste0(path, "ldl_cpra_list.txt"), col.names = FALSE, row.names = FALSE)
fwrite(sumstat_new[,c("cpra", "alt")], paste0(path, "ldl_cpra_alt_list.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")

plink --bfile /finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9 \
--extract ldl_cpra_list.txt \
--recode-allele ldl_cpra_alt_list.txt \
--recode A --threads 16 --out ldl

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require(gridExtra)) {install.packages("gridExtra"); library(gridExtra)}

total = fread(file = paste0(path, "healthcare_cost_20220505.txt"), header = TRUE, stringsAsFactors = FALSE)
variants = fread(file = paste0(path, "gwas/ldl.raw"))
variants = variants[,-c("FID", "SEX")]
variant_list = colnames(variants)[5:length(colnames(variants))]
total_new = merge(total, variants, by.x = "finngenid", by.y = "IID")
total_temp = total_new
total_temp$birth_year_squared = total_temp$birth_year ^ 2

mod_res = NULL
for (i in 1:length(variant_list)) {
	total_temp$variant = as.numeric(as.character(total_temp[[variant_list[i]]]))
	mod_unweight = lm(formula = log_total_cost_per_person_year ~ variant + birth_year + birth_year_squared + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_temp[which(total_temp$rx_statin == 0),])
	coef_unweight = coef(summary(mod_unweight))["variant", "Estimate"]
	ste_unweight = coef(summary(mod_unweight))["variant", "Std. Error"]
	pval_unweight = coef(summary(mod_unweight))["variant", "Pr(>|t|)"]
	mod_weight = lm(formula = log_total_cost_per_person_year ~ variant + birth_year + birth_year_squared + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_temp[which(total_temp$rx_statin == 1),])
	coef_weight = coef(summary(mod_weight))["variant", "Estimate"]
	ste_weight = coef(summary(mod_weight))["variant", "Std. Error"]
	pval_weight = coef(summary(mod_weight))["variant", "Pr(>|t|)"]
	mod = lm(formula = log_total_cost_per_person_year ~ variant + birth_year + birth_year_squared + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = total_temp)
	coef = coef(summary(mod))["variant", "Estimate"]
	ste = coef(summary(mod))["variant", "Std. Error"]
	pval = coef(summary(mod))["variant", "Pr(>|t|)"]
	mod_res = rbind(mod_res, c(variant_list[i], coef_unweight, ste_unweight, pval_unweight, coef_weight, ste_weight, pval_weight, coef, ste, pval))
	print(i)
}
mod_res = as.data.frame(mod_res)
colnames(mod_res) = c("variant", "beta_no_statin", "ste_no_statin", "pval_no_statin", "beta_yes_statin", "ste_yes_statin", "pval_yes_statin", "beta", "ste", "pval")
mod_res$beta_no_statin = as.numeric(mod_res$beta_no_statin)
mod_res$ste_no_statin = as.numeric(mod_res$ste_no_statin)
mod_res$pval_no_statin = as.numeric(mod_res$pval_no_statin)
mod_res$beta_yes_statin = as.numeric(mod_res$beta_yes_statin)
mod_res$ste_yes_statin = as.numeric(mod_res$ste_yes_statin)
mod_res$pval_yes_statin = as.numeric(mod_res$pval_yes_statin)
mod_res$beta = as.numeric(mod_res$beta)
mod_res$ste = as.numeric(mod_res$ste)
mod_res$pval = as.numeric(mod_res$pval)
fwrite(mod_res, paste0(path, "statin_stratified_linear_regression_results.txt"), sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}

key = fread(file = paste0(path, "gwas/20220505/log_total_cost_per_person_year_final.txt.gz"), header = TRUE, stringsAsFactors = FALSE)
sumstat = fread(file = paste0(path, "gwas/20220505/statin_stratified_linear_regression_results.txt"), header = TRUE, stringsAsFactors = FALSE)
temp = strsplit(sumstat$variant, "_")
sumstat$chr = sapply(temp, "[", 1)
sumstat$pos = as.numeric(sapply(temp, "[", 2))
sumstat$ref = sapply(temp, "[", 3)
sumstat$alt = sapply(temp, "[", 4)
sumstat$chr = as.numeric(substring(sumstat$chr, 4))

sumstat = merge(key[,c("chr", "pos", "cpra")], sumstat, by = c("chr", "pos"))

exposure_dat = extract_instruments(outcomes = "ieu-b-110")
#exposure_dat = extract_instruments(outcomes = "ukb-b-1668")

outcome_temp = sumstat[which(sumstat$cpra %in% exposure_dat$SNP),]
fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
outcome_dat = read_outcome_data(snps = NULL,
	filename = paste0(path, "mr/", today, "/temp.csv"),
	sep = ",",
	snp_col = "cpra",
	beta_col = "beta_no_statin",
  se_col = "ste_no_statin",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval_no_statin"
)
# Harmonise the exposure and outcome data
print("Harmonizing data...")
dat = harmonise_data(exposure_dat, outcome_dat)
# Perform MR
print("Performing MR...")
res = mr(dat)

yes_statin = res
no_statin = res
all_statin = res

all_statin$category = "All Statin"
no_statin$category = "No Statin"
yes_statin$category = "Yes Statin"

total = rbind(all_statin, rbind(no_statin, yes_statin))
total$percent = (exp(total$b)-1)*100
total$percent_lower = (exp(total$b-1.96*total$se)-1)*100
total$percent_upper = (exp(total$b+1.96*total$se)-1)*100
total_new = total[which(total$method == "Inverse variance weighted"),]
fwrite(total_new, "ldl_sensitivity_analysis.txt")

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
#gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}

#sumstat = fread(paste0(path, "gwas/", today, "/log_total_cost_per_person_year_final.txt.gz"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
#sumstat = fread(paste0("/medpop/esp2/jiwoolee/met_rep/gwas/ukbb_cad_sumstat_clean.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
sumstat = fread(paste0(path, "gwas/", today, "/finngen_R9_I9_CHD.gz"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)

exposure_dat = extract_instruments(outcomes = "ieu-b-110")
#exposure_dat[which(exposure_dat$SNP == "rs472495"),]
#exposure_dat = extract_instruments(outcomes = "ieu-a-300")
#exposure_dat[which(exposure_dat$SNP == "rs472495"),]

outcome_temp = sumstat[which(sumstat$rsids %in% exposure_dat$SNP),]
fwrite(outcome_temp, paste0(path, "mr/", today, "/temp.csv"))
outcome_dat = read_outcome_data(snps = NULL,
	filename = paste0(path, "mr/", today, "/temp.csv"),
	sep = ",",
	snp_col = "rsids",
	beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval"
)
# Harmonise the exposure and outcome data
print("Harmonizing data...")
dat = harmonise_data(exposure_dat, outcome_dat)
# Perform MR
print("Performing MR...")
res = mr(dat)

##################################################
########## EXTRACT VARIANTS TO MAKE LDL ALLELE SCORE
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

glgc = fread(file = paste0(path, "gwas/glgc_ldl_rsid_list.txt"), header = FALSE, stringsAsFactors = FALSE)
ukbb = fread(file = paste0(path, "gwas/ukbb_ldl_rsid_list.txt"), header = FALSE, stringsAsFactors = FALSE)
colnames(glgc) = c("rsid")
colnames(ukbb) = c("rsid")

sumstat = fread(paste0(path, "gwas/", today, "/", "log_total_cost_per_person_year", ".txt"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)

glgc_new = sumstat[which(sumstat$rsid %in% glgc$rsid),]
ukbb_new = sumstat[which(sumstat$rsid %in% ukbb$rsid),]
glgc_new$cpra = paste0("chr", glgc_new$chr_pos_hg38, "_", glgc_new$ref, "_", glgc_new$alt)
ukbb_new$cpra = paste0("chr", ukbb_new$chr_pos_hg38, "_", ukbb_new$ref, "_", ukbb_new$alt)
ukbb_temp = ukbb_new %>% distinct() 
fwrite(glgc_new[,c("cpra", "alt", "beta")], paste0(path, "gwas/", today, "/glgc_ldl_sumstat.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
fwrite(ukbb_temp[,c("cpra", "alt", "beta")], paste0(path, "gwas/", today, "/ukbb_ldl_sumstat.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")

plink --bfile /finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9 \
--score glgc_ldl_sumstat.txt 1 2 3 center --out glgc_ldl

plink --bfile /finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9 \
--score ukbb_ldl_sumstat.txt 1 2 3 center --out ukbb_ldl

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"
today = "20220505"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}

#sumstat = fread(file = paste0(path, "gwas/", today, "/log_total_cost_per_person_year.txt"), sep = " ", header = TRUE, stringsAsFactors = FALSE)
#sumstat$cpra = paste0("chr", sumstat$chr, "_", sumstat$pos, "_", sumstat$ref, "_", sumstat$alt)

total = fread(file = paste0(path, "healthcare_cost_20220505.txt"), header = TRUE, stringsAsFactors = FALSE)
glgc_ldl = fread(paste0(path, "gwas/glgc_ldl.profile"))
ukbb_ldl = fread(paste0(path, "gwas/ukbb_ldl.profile"))

glgc_ldl = glgc_ldl[,c("IID", "SCORE")]
colnames(glgc_ldl) = c("finngenid", "glgc_ldl_score")

ukbb_ldl = ukbb_ldl[,c("IID", "SCORE")]
colnames(ukbb_ldl) = c("finngenid", "ukbb_ldl_score")

total_new = merge(total, glgc_ldl, by = "finngenid")
total_new = merge(total_new, ukbb_ldl, by = "finngenid")

glgc_mod = lm(log_total_cost_per_person_year~scale(glgc_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new)
glgc_mod_nostat = lm(log_total_cost_per_person_year~scale(glgc_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new[which(total_new$rx_statin == 0),])
glgc_mod_stat = lm(log_total_cost_per_person_year~scale(glgc_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new[which(total_new$rx_statin == 1),])

ukbb_mod = lm(log_total_cost_per_person_year~scale(ukbb_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new)
ukbb_mod_nostat = lm(log_total_cost_per_person_year~scale(ukbb_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new[which(total_new$rx_statin == 0),])
ukbb_mod_stat = lm(log_total_cost_per_person_year~scale(ukbb_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new[which(total_new$rx_statin == 1),])

coef(summary(glgc_mod))[1:5,] # all individuals, sumstat adjusted for statin
coef(summary(glgc_mod_stat))[1:5,] # statin individuals, sumstat adjusted for statin
coef(summary(glgc_mod_nostat))[1:5,] # no statin individuals, sumstat adjusted for statin

coef(summary(ukbb_mod))[1:5,] # all individuals, sumstat NOT adjusted for statin
coef(summary(ukbb_mod_stat))[1:5,] # statin individuals, sumstat NOT adjusted for statin
coef(summary(ukbb_mod_nostat))[1:5,] # no statin individuals, sumstat NOT adjusted for statin

glgc_chd = glm(i9_chd~scale(glgc_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new, family = "binomial")
glgc_chd_nostat = glm(i9_chd~scale(glgc_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new[which(total_new$rx_statin == 0),], family = "binomial")
glgc_chd_stat = glm(i9_chd~scale(glgc_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new[which(total_new$rx_statin == 1),], family = "binomial")

ukbb_chd = glm(i9_chd~scale(ukbb_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new, family = "binomial")
ukbb_chd_nostat = glm(i9_chd~scale(ukbb_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new[which(total_new$rx_statin == 0),], family = "binomial")
ukbb_chd_stat = glm(i9_chd~scale(ukbb_ldl_score)+birth_year+birth_year^2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = total_new[which(total_new$rx_statin == 1),], family = "binomial")

coef(summary(glgc_chd))[1:5,] # all individuals, sumstat adjusted for statin
coef(summary(glgc_chd_stat))[1:5,] # statin individuals, sumstat adjusted for statin
coef(summary(glgc_chd_nostat))[1:5,] # no statin individuals, sumstat adjusted for statin

coef(summary(ukbb_chd))[1:5,] # all individuals, sumstat NOT adjusted for statin
coef(summary(ukbb_chd_stat))[1:5,] # statin individuals, sumstat NOT adjusted for statin
coef(summary(ukbb_chd_nostat))[1:5,] # no statin individuals, sumstat NOT adjusted for statin

res_cost = as.data.frame(rbind(c("Cost~LDL in ALL (statin-adjusted sumstats)", coef(summary(glgc_mod))[2,]),
	c("Cost~LDL in statin users (statin-adjusted sumstats)", coef(summary(glgc_mod_stat))[2,]),
	c("Cost~LDL in non-statin users (statin-adjusted sumstats)", coef(summary(glgc_mod_nostat))[2,]),
	c("Cost~LDL in ALL (non-statin-adjusted sumstats)", coef(summary(ukbb_mod))[2,]),
	c("Cost~LDL in statin users (non-statin-adjusted sumstats)", coef(summary(ukbb_mod_stat))[2,]),
	c("Cost~LDL in non-statin users (non-statin-adjusted sumstats)", coef(summary(ukbb_mod_nostat))[2,])))
colnames(res_cost) = c("model", "b", "se", "t", "p")
res_cost$my_order = 1:nrow(res_cost)
res_cost$b = as.numeric(res_cost$b)
res_cost$se = as.numeric(res_cost$se)
res_cost$p = as.numeric(res_cost$p)

res_cost$beta = (exp(res_cost$b)-1)*100
res_cost$beta_lower = (exp(res_cost$b - 1.96 * res_cost$se)-1)*100
res_cost$beta_upper = (exp(res_cost$b + 1.96 * res_cost$se)-1)*100

ggplot() +
	coord_flip() +
	geom_point(data = res_cost, mapping = aes(x = my_order, y = beta, color = ifelse(res_cost$my_order <= 3, "red", "blue")), size = 3) +
	geom_errorbar(data = res_cost, mapping = aes(x = my_order, ymin = beta_lower, ymax = beta_upper, color = ifelse(res_cost$my_order <= 3, "red", "blue"))) +
	geom_label(data = res_cost, mapping = aes(x = my_order + 0.5, y = beta, label = paste("P = ", p))) +
	scale_x_continuous(breaks = res_cost$my_order, labels = res_cost$model) +
	labs(x = "Model", y = "Percent Change in Healthcare Costs") +
	guides(color = FALSE) +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

res_chd = as.data.frame(rbind(c("CHD~LDL in ALL (statin-adjusted sumstats)", coef(summary(glgc_chd))[2,]),
	c("CHD~LDL in statin users (statin-adjusted sumstats)", coef(summary(glgc_chd_stat))[2,]),
	c("CHD~LDL in non-statin users (statin-adjusted sumstats)", coef(summary(glgc_chd_nostat))[2,]),
	c("CHD~LDL in ALL (non-statin-adjusted sumstats)", coef(summary(ukbb_chd))[2,]),
	c("CHD~LDL in statin users (non-statin-adjusted sumstats)", coef(summary(ukbb_chd_stat))[2,]),
	c("CHD~LDL in non-statin users (non-statin-adjusted sumstats)", coef(summary(ukbb_chd_nostat))[2,])))
colnames(res_chd) = c("model", "b", "se", "t", "p")
res_chd$my_order = 1:nrow(res_chd)
res_chd$b = as.numeric(res_chd$b)
res_chd$se = as.numeric(res_chd$se)
res_chd$p = as.numeric(res_chd$p)
res_chd$beta = exp(res_chd$b)
res_chd$beta_lower = exp(res_chd$b - 1.96*res_chd$se)
res_chd$beta_upper = exp(res_chd$b + 1.96*res_chd$se)

ggplot() +
	coord_flip() +
	geom_point(data = res_chd, mapping = aes(x = my_order, y = beta, color = ifelse(res_chd$my_order <= 3, "red", "blue")), size = 3) +
	geom_errorbar(data = res_chd, mapping = aes(x = my_order, ymin = beta_lower, ymax = beta_upper, color = ifelse(res_chd$my_order <= 3, "red", "blue"))) +
	geom_label(data = res_chd, mapping = aes(x = my_order + 0.5, y = beta, label = paste("P = ", p))) +
	scale_x_continuous(breaks = res_chd$my_order, labels = res_chd$model) +
	labs(x = "Model", y = "Odds Ratio") +
	guides(color = FALSE) +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

temp = rbind(coef(summary(ukbb_mod))["scale(ukbb_ldl_score)", c("Estimate", "Std. Error", "Pr(>|t|)")],
	rbind(coef(summary(ukbb_mod_stat))["scale(ukbb_ldl_score)", c("Estimate", "Std. Error", "Pr(>|t|)")],
		coef(summary(ukbb_mod_nostat))["scale(ukbb_ldl_score)", c("Estimate", "Std. Error", "Pr(>|t|)")]))
temp = as.data.frame(temp)
temp$category = c("Regression of Cost~LDL PRS (All Statin)", "Regression of Cost~LDL PRS (Yes Statin)", "Regression of Cost~LDL PRS (No Statin)")
temp$percent = 100*(exp(temp$Estimate)-1)
temp$percent_lower = 100*(exp(temp$Estimate-1.96*temp[,"Std. Error"])-1)
temp$percent_upper = 100*(exp(temp$Estimate+1.96*temp[,"Std. Error"])-1)
temp$pval = temp[,"Pr(>|t|)"]
temp = temp[,c("category", "percent", "percent_lower", "percent_upper", "pval")]

total = fread("~/jiwoo/healthcare_cost_repository/ldl_sensitivity_analysis.txt")
total = total[,c("category", "percent", "percent_lower", "percent_upper", "pval")]
total$category = c("MR of Cost~LDL (All Statin)", "MR of Cost~LDL (Yes Statin)", "MR of Cost~LDL (No Statin)")

total_new = rbind(total, temp)
ggplot() +
	coord_flip() +
	geom_point(data = total_new, mapping = aes(x = category, y = percent), size = 3) +
	geom_errorbar(data = total_new, mapping = aes(x = category, ymin = percent_lower, ymax = percent_upper)) +
	geom_label(data = total_new, mapping = aes(x = category, y = percent, label = paste("P = ", pval))) +
	#scale_x_continuous(breaks = total$my_order, labels = total$model) +
	labs(x = "", y = "Percent Change in Healthcare Costs") +
	guides(color = FALSE) +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

##################################################
########## PLOT DIFFERENT ANALYSES TOGETHER
##################################################

path = "C:/Jiwoo Lee/Healthcare_Cost_Research_2021/"
today = "20220505"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}
if(!require(ggplot2)) {install_github("ggplot2"); library(ggplot2)}
if(!require(gridExtra)) {install_github("gridExtra"); library(gridExtra)}

exposures = c("ukb-a-382", "ukb-a-360", "ieu-a-835")
fg = fread(paste0(path, "MR/20220505/mr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
fg = fg[which(fg$id.exposure %in% exposures & fg$outcome == "log_total_cost_per_person_year" & fg$method == "Inverse variance weighted")]
fg$outcome = "Finland"
fg$percent = fg$percent/100*1312.53
fg$confint_lower = fg$confint_lower/100*1312.53
fg$confint_upper = fg$confint_upper/100*1312.53
nl = fread(paste0(path, "MR/20220505/nl_mr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
nl = nl[which(nl$id.exposure %in% exposures & nl$method == "Inverse variance weighted")]
nl$outcome = "Netherlands"
nl$percent = nl$percent/100*1312.53
nl$confint_lower = nl$confint_lower/100*1623.8492
nl$confint_upper = nl$confint_upper/100*1623.8492
wt = fread(paste0(path, "MR/20220505/weight_mr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
wt = wt[which(wt$id.exposure %in% exposures & wt$method == "Inverse variance weighted")]
wt$outcome = "Finland (Weighted)"
wt$percent = wt$percent/100*1312.53
wt$confint_lower = wt$confint_lower/100*1312.53
wt$confint_upper = wt$confint_upper/100*1312.53
uk = fread(paste0(path, "MR/20220505/uk_mr_res_total.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)
uk = uk[which(uk$id.exposure %in% exposures & uk$method == "Inverse variance weighted")]
uk$outcome = "United Kingdom"
uk$percent = uk$b * 1.1921
uk$confint_lower = (uk$b - uk$se) * 1.1921
uk$confint_upper = (uk$b + uk$se) * 1.1921

total = rbind(fg, rbind(nl, rbind(uk, wt)))
total$exposure[1:3] = c("Waist circumference", "Body mass index", "Systolic blood pressure")
total = total[order(total$exposure, total$outcome),]
total$number = c(1:4, 6:9, 11:14)
total$confint_lower = as.numeric(total$confint_lower)
total$exposure = ifelse(total$number %in% c(2, 7, 12), total$exposure, "")

ggplot() +
	geom_point(data = total, mapping = aes(x = -number, y = percent, color = outcome), size = 3) +
	geom_errorbar(data = total, mapping = aes(x = -number, ymin = confint_lower, ymax = confint_upper, color = outcome)) +
	coord_flip() + 
	scale_x_continuous(breaks = -total$number-0.5, labels = total$exposure) +
	#scale_y_continuous(breaks = seq() +
	scale_color_discrete(name = NULL) +
	labs(x = "", y = "Annual Change in Healthcare Cost in Euros per 1 SD Increase in Risk Factor") +
	theme(axis.ticks.y = element_blank(), legend.position = "bottom", axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
ggsave(paste0(path, "MR/", today, "/", "all_comparison.jpg"), device = "jpg", type = "cairo", height = 10, width = 10)

##################################################
########## DO REVISION ANALYSIS
##################################################

c(264.5807, 483.7110, 713.6994, 981.9374, 1312.5304)

snps = fread(paste0(path, "SNPlist.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, fill = TRUE)

##################################################
########## QC SUMMARY STATISTICS
##################################################

ish -l h_rt=12:00:00 -l h_vmem=50g
use R-4.1 
use Anaconda3
use Google-Cloud-SDK
gcloud auth login --no-launch-browser

cd /medpop/esp2/jiwoolee/

R

path = "/medpop/esp2/jiwoolee/cost_rep/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(ggrepel)) {install.packages("ggrepel"); library(ggplot2)}
if(!require(TwoSampleMR)) {install.packages("TwoSampleMR"); library(TwoSampleMR)}

 snps$LOG10P > (-log(5e-8)

snps = fread(paste0(path, "gwas/20220505/log_total_cost_per_person_year.regenie"))
snps_new = snps[which(snps$A1FREQ > 0.001 & snps$INFO > 0.8),c("CHROM", "GENPOS", "INFO")]
colnames(snps_new) = c("chr", "pos", "info")

total = fread(paste0(path, "gwas/20220505/log_total_cost_per_person_year_final.txt.gz"))

total_new = merge(total, snps_new, by = c("chr", "pos"))

total_formatted = format_data(total_new,
	type = "outcome",
	snp_col = "cpra",
	beta_col = "beta",
	se_col = "sebeta",
	effect_allele_col = "alt",
	other_allele_col = "ref",
	eaf_col = "af_alt",
	pval_col = "pval")
colnames(total_formatted) = c("chr", "pos", "rsid", "ref", "alt", "pval", "beta", "se", "eaf", "info", "outcome", "keep", "pval_origin", "id")

total_clumped = NULL
chromosomes = unique(total_formatted$chr)
for (i in 1:length(chromosomes)) {
	print(chromosomes[i])
	temp = total_formatted[which(total_formatted$chr == chromosomes[i]),]
	final = ieugwasr::ld_clump(temp)
	total_clumped = rbind(total_clumped, final)
}

backup = total_clumped

total_final = total_clumped[which(total_clumped$pval < 5e-8),]

fwrite(total_clumped, paste0(path, "gwas_cost_clumped.csv"))
fwrite(total_final, paste0(path, "gwas_cost_clumped_significant.csv"))





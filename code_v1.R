##################################################
########## DOWNLOAD DATA
##################################################

cd ~/jiwoo/healthcare_cost_repository

# PCA DATA
cp /finngen/library-red/finngen_R9/pca_1.0/data/finngen_R9.eigenvec.txt .

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

# LONGITUDINAL DATA
cp /finngen/pipeline/DATA_IMPORT_28102021/Ganna_data_1.0/data/finngen_R8_Ganna_detailed_longitudinal.txt.gz .
gunzip finngen_R8_Ganna_detailed_longitudinal.txt.gz

# COVARIATE DATA
cp /finngen/library-red/finngen_R9/analysis_covariates/R9_COV_PHENO_V1.txt.gz .
gunzip R9_COV_PHENO_V1.txt.gz

# GITHUB DATA
cp /finngen/green/jiwoo/* .

# SAKARI CODE
cp /finngen/red/sakari/simulate_lifetime_costs_bootstrap_sandbox.R ./code/

# SAKARI RESULTS
cp /finngen/red/jiwoo/sakari_folder/* .

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

# LOAD ENDPOINT DATA
endpoint_all = fread(file = paste0(path, "finngen_R8_endpoint_clean.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(endpoint_all) = tolower(colnames(endpoint_all))
endpoint_all$length_of_fu = endpoint_all$fu_end_age - endpoint_all$bl_age
endpoint_all$birth_year = endpoint_all$bl_year - endpoint_all$bl_age
endpoint_all$fu_start_age = ifelse(1995 - endpoint_all$birth_year < 0, 0, 1995 - endpoint_all$birth_year)
endpoint_all$fu_end_year = endpoint_all$birth_year + endpoint_all$fu_end_age
endpoint_all$person_day_since_1998 = (endpoint_all$fu_end_year - 1998) * 365.25
endpoint_all$person_day_since_2011 = (endpoint_all$fu_end_year - 2011) * 365.25
endpoint_all$person_day_since_2011[which(endpoint_all$person_day_since_2011 < 0)] = 0
endpoint_all = endpoint_all[,c("finngenid", "sex", "birth_year", "death_year", "death_age", "death", "bl_age", "bl_year", "fu_start_age", "fu_end_year", "fu_end_age", "length_of_fu", "person_day_since_1998", "person_day_since_2011", "rx_statin", "rx_statin_age", "rx_statin_year", "i9_chd", "i9_chd_age", "i9_chd_year")]
endpoint_new = endpoint_all[which(endpoint_all$person_day_since_1998 > 0),]
dim(endpoint_all) # 356077
dim(endpoint_new) # 352374

# LOAD PCA DATA
pca_all = fread(file = paste0(path, "finngen_R8.eigenvec.txt"), header = TRUE, stringsAsFactors = FALSE)
pca_new = pca_all[,-"#FID"]
colnames(pca_new)[which(colnames(pca_new) == "IID")] = "finngenid"
dim(pca_new) # 342448
length(unique(pca_new$finngenid)) # 342448
total = merge(endpoint_new, pca_new, by = "finngenid", all.x = TRUE)
dim(total) # 352374

# LOAD DRUG DICTIONARIES
drug_dict = fread(file = paste0(path, "avg_vnro_price.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(drug_dict) = c("year", "vnr", "drug_cost")
drug_dict$drug_cost[which(drug_dict$drug_cost == "Inf")] = NA

# LOAD HILMO DICTIONARIES
hilmo_dict = fread(file = paste0(path, "hilmo_unit_costs_long.txt"), header = TRUE, stringsAsFactors = FALSE)
ward_days = fread(file = paste0(path, "ward_days.txt"), header = TRUE, stringsAsFactors = FALSE)
hilmo_dict = merge(hilmo_dict, ward_days, by = c("year", "hospital_type", "specialty"))
hilmo_dict$specialty = as.character(hilmo_dict$specialty)
hilmo_dict = hilmo_dict[which(hilmo_dict$year == 2017),]

# LOAD AVOHILMO DICTIONARIES
avohilmo_dict = fread(file = paste0(path, "avohilmo_unit_costs_long.txt"), header = TRUE, stringsAsFactors = FALSE)
avohilmo_dict = avohilmo_dict[,-"contact_type"]
contact_dict = fread(file = paste0(path, "avohilmo_contact_type_labels.txt"), header = TRUE, stringsAsFactors = FALSE)
profession_dict = fread(file = paste0(path, "avohilmo_profession_labels.txt"), header = TRUE, stringsAsFactors = FALSE)
profession_dict$profession = as.character(profession_dict$profession)
service_dict  = fread(file = paste0(path, "avohilmo_service_type_labels.txt"), header = TRUE, stringsAsFactors = FALSE)

# LOAD LONGITUDINAL DATA
long_all = fread(file = paste0(path, "finngen_R8_Ganna_detailed_longitudinal.txt"), header = TRUE, stringsAsFactors = FALSE)
long_all$APPROX_EVENT_YEAR = as.integer(format(as.Date(long_all$APPROX_EVENT_DAY), format = "%Y"))
colnames(long_all) = tolower(colnames(long_all))

##################################################
########## MERGE DATA
##################################################

# MERGE DRUG DATA
drug_all = long_all[which(long_all$source == "PURCH"), c("finngenid", "code3", "event_age", "approx_event_year")]
colnames(drug_all)[which(colnames(drug_all) == "code3")] = "vnr"
drug_new = merge(drug_all, drug_dict, by.x = c("vnr", "approx_event_year"), by.y = c("vnr", "year"), all.x = TRUE)

# MERGE HILMO DATA
hilmo_all = long_all[which(long_all$source == "INPAT" | long_all$source == "OUTPAT" | long_all$source == "OPER_IN" | long_all$source == "OPER_OUT"), c("finngenid", "source", "code5", "code6", "code7", "event_age", "approx_event_year")]
colnames(hilmo_all) = c("finngenid", "source", "service_type", "specialty", "hospital_type", "event_age", "approx_event_year")
hilmo_all = hilmo_all %>% mutate(service_type = case_when(source == "OUTPAT" & service_type == "91" ~ "emergency",
												source == "OUTPAT" ~ "outpatient",
												source == "INPAT" ~ "ward"))
hilmo_all = hilmo_all %>% mutate(hospital_type = case_when(hospital_type == "Central Hospital" ~ "central",
												hospital_type == "Other Hospital" ~ "other",
												hospital_type == "University Hospital" ~ "university"))
hilmo_new = merge(hilmo_all, hilmo_dict, by = c("service_type", "specialty", "hospital_type"), all.x = TRUE)
hilmo_new$hospital_type = ifelse(is.na(hilmo_new$hospital_type), "all", hilmo_new$hospital_type)

# MERGE AVOHILMO DATA
avohilmo_all = long_all[which(long_all$source == "PRIM_OUT"), c("finngenid", "code5", "code6", "code7", "event_age", "approx_event_year")]
colnames(avohilmo_all) = c("finngenid", "contact_type", "service_type", "profession", "event_age", "approx_event_year")
avohilmo_new = merge(avohilmo_all, profession_dict, by = "profession", all.x = TRUE)
avohilmo_new = merge(avohilmo_new, service_dict, by = "service_type", all.x = TRUE)
avohilmo_new = merge(avohilmo_new, contact_dict, by = "contact_type", all.x = TRUE)
avohilmo_new = merge(avohilmo_new, avohilmo_dict, by = c("profession_label", "service_type_label", "contact_type_label"), all.x = TRUE)
avohilmo_new$profession_label = ifelse(is.na(avohilmo_new$profession_label), "other", avohilmo_new$profession_label)

# COMBINE DRUG, HILMO, AND AVOHILMO INTO LONG DATAFRAME
drug_temp = drug_new[which(drug_new$approx_event_year >= 1998),]
hilmo_temp = hilmo_new[which(hilmo_new$approx_event_year >= 1998),]
avohilmo_temp = avohilmo_new[which(avohilmo_new$approx_event_year >= 2011),]
long_new = rbindlist(list(setnames(drug_temp[,c("finngenid", "event_age", "drug_cost")], c("finngenid", "event_age", "cost")),
	hilmo_temp[,c("finngenid", "event_age", "cost")],
	avohilmo_temp[,c("finngenid", "event_age", "cost")]))

##################################################
########## CALCULATE HEALTHCARE COSTS
##################################################

# CALCULATE DRUG COSTS
drug_cost = drug_temp %>% group_by(finngenid) %>% summarise(drug_cost = sum(drug_cost, na.rm = TRUE))
drug_cost = as.data.frame(drug_cost)

# CALCULATE SUBTYPE COSTS
calculate_subtype_cost = function(data, source, subtype) {
	res = data %>% group_by(finngenid, !!as.name(subtype)) %>% summarise(cost = sum(cost, na.rm = TRUE))
	res = as.data.frame(res)
	res[,subtype] = paste0(source, "_", subtype, "_", convert_camel_case(res[,subtype]), "_cost")
	res = res %>% spread(!!as.name(subtype), cost)
	return(res)
}

# CALCULATE HILMO COSTS
hilmo_service_type_cost = calculate_subtype_cost(hilmo_temp, "hilmo", "service_type")
hilmo_specialty_label_cost = calculate_subtype_cost(hilmo_temp, "hilmo", "specialty_label")
hilmo_hospital_type_cost = calculate_subtype_cost(hilmo_temp, "hilmo", "hospital_type")

#CALCULATE AVOHILMO COSTS
avohilmo_profession_label_cost = calculate_subtype_cost(avohilmo_temp, "avohilmo", "profession_label")
avohilmo_service_type_label_cost = calculate_subtype_cost(avohilmo_temp, "avohilmo", "service_type_label")
avohilmo_contact_type_label_cost = calculate_subtype_cost(avohilmo_temp, "avohilmo", "contact_type_label")

# CHECK HILMO COSTS 
hilmo_service_type = as.data.frame(cbind(finngenid = hilmo_service_type_cost$finngenid, hilmo_service_type_cost = rowSums(hilmo_service_type_cost[grep('^hilmo', colnames(hilmo_service_type_cost))], na.rm = TRUE)))
hilmo_specialty_label = as.data.frame(cbind(finngenid = hilmo_specialty_label_cost$finngenid, hilmo_specialty_label_cost = rowSums(hilmo_specialty_label_cost[grep('^hilmo', colnames(hilmo_specialty_label_cost))], na.rm = TRUE)))
hilmo_hospital_type = as.data.frame(cbind(finngenid = hilmo_hospital_type_cost$finngenid, hilmo_hospital_type_cost = rowSums(hilmo_hospital_type_cost[grep('^hilmo', colnames(hilmo_hospital_type_cost))], na.rm = TRUE)))
hilmo_cost = merge(hilmo_service_type, hilmo_specialty_label, by = "finngenid", all = TRUE)
hilmo_cost = merge(hilmo_cost, hilmo_hospital_type, by = "finngenid", all = TRUE)
hilmo_cost$hilmo_service_type_cost = as.numeric(hilmo_cost$hilmo_service_type_cost)
hilmo_cost$hilmo_specialty_label_cost = as.numeric(hilmo_cost$hilmo_specialty_label_cost)
hilmo_cost$hilmo_hospital_type_cost = as.numeric(hilmo_cost$hilmo_hospital_type_cost)
hilmo_cost$hilmo_cost = rowMeans(hilmo_cost[grep('^hilmo', colnames(hilmo_cost))], na.rm = TRUE)
length(which(hilmo_cost$hilmo_service_type_cost - hilmo_cost$hilmo_specialty_label_cost > 1))
length(which(hilmo_cost$hilmo_service_type_cost - hilmo_cost$hilmo_hospital_type_cost > 1))
length(which(hilmo_cost$hilmo_specialty_label_cost - hilmo_cost$hilmo_hospital_type_cost > 1))
hilmo_cost = hilmo_cost[,c("finngenid", "hilmo_cost")]

# CHECK AVOHILMO COSTS
avohilmo_profession_label = as.data.frame(cbind(finngenid = avohilmo_profession_label_cost$finngenid, avohilmo_profession_label_cost = rowSums(avohilmo_profession_label_cost[grep('^avohilmo', colnames(avohilmo_profession_label_cost))], na.rm = TRUE)))
avohilmo_service_type_label = as.data.frame(cbind(finngenid = avohilmo_service_type_label_cost$finngenid, avohilmo_service_type_label_cost = rowSums(avohilmo_service_type_label_cost[grep('^avohilmo', colnames(avohilmo_service_type_label_cost))], na.rm = TRUE)))
avohilmo_contact_type_label = as.data.frame(cbind(finngenid = avohilmo_contact_type_label_cost$finngenid, avohilmo_contact_type_label_cost = rowSums(avohilmo_contact_type_label_cost[grep('^avohilmo', colnames(avohilmo_contact_type_label_cost))], na.rm = TRUE)))
avohilmo_cost = merge(avohilmo_profession_label, avohilmo_service_type_label, by = "finngenid", all = TRUE)
avohilmo_cost = merge(avohilmo_cost, avohilmo_contact_type_label, by = "finngenid", all = TRUE)
avohilmo_cost$avohilmo_profession_label_cost = as.numeric(avohilmo_cost$avohilmo_profession_label_cost)
avohilmo_cost$avohilmo_service_type_label_cost = as.numeric(avohilmo_cost$avohilmo_service_type_label_cost)
avohilmo_cost$avohilmo_contact_type_label_cost = as.numeric(avohilmo_cost$avohilmo_contact_type_label_cost)
avohilmo_cost$avohilmo_cost = rowMeans(avohilmo_cost[grep('^avohilmo', colnames(avohilmo_cost))], na.rm = TRUE)
length(which(avohilmo_cost$avohilmo_profession_label_cost - avohilmo_cost$avohilmo_service_type_label_cost > 1))
length(which(avohilmo_cost$avohilmo_profession_label_cost - avohilmo_cost$avohilmo_contact_type_label_cost > 1))
length(which(avohilmo_cost$avohilmo_service_type_label_cost - avohilmo_cost$avohilmo_contact_type_label_cost > 1))
avohilmo_cost = avohilmo_cost[,c("finngenid", "avohilmo_cost")]

# MERGE HEALTHCARE COSTS
total_new = merge(total, drug_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, hilmo_service_type_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, hilmo_specialty_label_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, hilmo_hospital_type_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, avohilmo_profession_label_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, avohilmo_service_type_label_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, avohilmo_contact_type_label_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, hilmo_cost, by = "finngenid", all.x = TRUE)
total_new = merge(total_new, avohilmo_cost, by = "finngenid", all.x = TRUE)

dim(total_new) # 352374
total_new = total_new[!which(is.na(total_new$drug_cost) & is.na(total_new$hilmo_cost) & is.na(total_new$avohilmo_cost)),]
dim(total_new) # 352080

total_new$drug_cost_per_person_year = (total_new$drug_cost / total_new$person_day_since_1998) * 365.25
total_new$hilmo_cost_per_person_year = (total_new$hilmo_cost / total_new$person_day_since_1998) * 365.25
total_new$avohilmo_cost_per_person_year = (total_new$avohilmo_cost / total_new$person_day_since_2011) * 365.25

total_new$drug_cost_per_person_year[which(is.na(total_new$drug_cost_per_person_year))] = 0
total_new$hilmo_cost_per_person_year[which(is.na(total_new$hilmo_cost_per_person_year))] = 0
total_new$avohilmo_cost_per_person_year[which(((!is.na(total_new$drug_cost) | !is.na(total_new$hilmo_cost)) & is.na(total_new$avohilmo_cost)) | total_new$avohilmo_cost_per_person_year == Inf | is.na(total_new$avohilmo_cost_per_person_year))] = 
	median(total_new$avohilmo_cost_per_person_year, na.rm = TRUE)

total_new$log_drug_cost_per_person_year = log(total_new$drug_cost_per_person_year + 1)
total_new$log_hilmo_cost_per_person_year = log(total_new$hilmo_cost_per_person_year + 1)
total_new$log_avohilmo_cost_per_person_year = log(total_new$avohilmo_cost_per_person_year + 1)

total_new$total_cost = rowSums(total_new[,c("drug_cost", "hilmo_cost", "avohilmo_cost")], na.rm = TRUE)
total_new$total_cost_per_person_year = rowSums(total_new[,c("drug_cost_per_person_year", "hilmo_cost_per_person_year", "avohilmo_cost_per_person_year")], na.rm = TRUE)
total_new$log_total_cost_per_person_year = log(total_new$total_cost_per_person_year + 1)

##################################################
########## FINALIZE HEALTHCARE COSTS
##################################################

# FINALIZE DATAFRAME
fwrite(total_new, paste0(path, "healthcare_cost_20220214.txt"), row.names = FALSE, col.names = TRUE)

# ADD VARIANTS TO SAKARI DATAFRAME
variants = fread(paste0(path, "gwas/bmi_variants_sig.raw"))
variants = variants[,-c("FID", "PAT", "MAT", "SEX", "PHENOTYPE")]
colnames(variants)[which(colnames(variants) == "IID")] = "finngenid"

scores = fread(paste0(path, "gwas/plink.profile"))
scores = scores[,-c("FID", "PHENO", "CNT")]
colnames(scores) = c("finngenid", "unweighted_allele_score", "weighted_allele_score")

# MAKE SAKARI DATAFRAME WITH DRUG, HILMO, AND AVOHILMO >= 2011
drug_temp = drug_new[which(drug_new$approx_event_year >= 2011),]
hilmo_temp = hilmo_new[which(hilmo_new$approx_event_year >= 2011),]
avohilmo_temp = avohilmo_new[which(avohilmo_new$approx_event_year >= 2011),]
long_new = rbindlist(list(setnames(drug_temp[,c("finngenid", "event_age", "drug_cost")], c("finngenid", "event_age", "cost")),
	hilmo_temp[,c("finngenid", "event_age", "cost")],
	avohilmo_temp[,c("finngenid", "event_age", "cost")]))
# MAKE SAKARI DATAFRAME WITH DRUG AND HILMO >= 1998
drug_temp = drug_new[which(drug_new$approx_event_year >= 1998),]
hilmo_temp = hilmo_new[which(hilmo_new$approx_event_year >= 1998),]
long_new = rbindlist(list(setnames(drug_temp[,c("finngenid", "event_age", "drug_cost")], c("finngenid", "event_age", "cost")),
	hilmo_temp[,c("finngenid", "event_age", "cost")]))
# MAKE SAKARI DATAFRAME WITH AVOHILMO >= 2011
avohilmo_temp = avohilmo_new[which(avohilmo_new$approx_event_year >= 2011),]
long_new = avohilmo_new[,c("finngenid", "event_age", "cost")]
# MAKE SAKARI DATAFRAME
sakari = long_new
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
sakari$age_interval = round_any(sakari$event_age, 5, f = floor)
sakari$age_interval = ifelse(sakari$age_interval > 80, 80, sakari$age_interval)
detach(package:plyr)
sakari = sakari %>% group_by(finngenid, age_interval) %>% summarise(cost = sum(cost, na.rm = TRUE)) 
sakari = as.data.frame(sakari)
sakari_new = merge(total, sakari, by = "finngenid", all.x = TRUE)
sakari_new$person_time = ifelse((sakari_new$age_interval == 80) | (sakari_new$fu_end_age < sakari_new$age_interval + 5), sakari_new$fu_end_age - sakari_new$age_interval, 5)
dim(sakari_new)
sakari_new = sakari_new[which(sakari_new$person_time >= 0.5),]
dim(sakari_new)
sakari_new$cost_per_person_time = sakari_new$cost / sakari_new$person_time
sakari_new = sakari_new %>% mutate_at(all_of(prs_list), ntile, 10)
sakari_final = sakari_new[,c(1, 160:163, 2, 38:159)]
colnames(sakari_final)[2:6] = c("agegroup", "cost", "py", "cost_per_py", "female")
sakari_final = merge(sakari_final, variants, by = "finngenid")
sakari_final = merge(sakari_final, scores, by = "finngenid")
fwrite(sakari_final, "/finngen/red/jiwoo/bmi_cost_per_py_avohilmo.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)

# MAKE GWAS DATAFRAME
cov_all = fread(file = paste0(path, "finngen_R8_cov_1.0.txt"), header = TRUE, stringsAsFactors = FALSE)
temp = total_new[,c("finngenid", "i9_chd", "birth_year", "drug_cost_per_person_year", "hilmo_cost_per_person_year", "avohilmo_cost_per_person_year", "total_cost_per_person_year",
	"log_drug_cost_per_person_year", "log_hilmo_cost_per_person_year", "log_avohilmo_cost_per_person_year", "log_total_cost_per_person_year")]
cov_new = merge(cov_all, temp, by.x = "FINNGENID", by.y = "finngenid")
cov_new$birth_year_squared = cov_new$birth_year ^ 2
cov_new$total_cost_per_female_year = ifelse(cov_new$SEX == "female", cov_new$total_cost_per_person_year, NA)
cov_new$total_cost_per_male_year = ifelse(cov_new$SEX == "male", cov_new$total_cost_per_person_year, NA)
cov_new$log_total_cost_per_female_year = ifelse(cov_new$SEX == "female", cov_new$log_total_cost_per_person_year, NA)
cov_new$log_total_cost_per_male_year = ifelse(cov_new$SEX == "male", cov_new$log_total_cost_per_person_year, NA)
cov_new$total_cost_per_person_year_cad_exclusion = ifelse(cov_new$i9_chd == 0, cov_new$total_cost_per_person_year, NA)
cov_new$log_total_cost_per_person_year_cad_exclusion = ifelse(cov_new$i9_chd == 0, cov_new$log_total_cost_per_person_year, NA)
fwrite(cov_new, "/finngen/red/jiwoo/saige_healthcare_cost/healthcare_cost_20220214.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

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

total_new = fread(file = paste0(path, "healthcare_cost_20220214.txt"), header = TRUE, stringsAsFactors = FALSE)

# PLOT DISTRIBUTION OF HEALTHCARE COSTS
options(scipen = 999)
ggplot() +
	geom_histogram(data = total_new, mapping = aes(x = total_cost_per_person_year, fill = "Total"), alpha = 0.5, bins = 100) +
	geom_histogram(data = total_new, mapping = aes(x = drug_cost_per_person_year, fill = "Kela"), alpha = 0.5, bins = 100) +
	geom_histogram(data = total_new, mapping = aes(x = hilmo_cost_per_person_year, fill = "HILMO"), alpha = 0.5, bins = 100) +
	geom_histogram(data = total_new, mapping = aes(x = avohilmo_cost_per_person_year, fill = "AvoHILMO"), alpha = 0.5, bins = 100) +
	labs(x = "Annual Healthcare Cost", y = "Frequency", fill = "") + 
	xlim(0, 10000) +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_histogram(data = total_new, mapping = aes(x = log_total_cost_per_person_year, fill = "Total"), alpha = 0.5, bins = 100) +
	geom_histogram(data = total_new, mapping = aes(x = log_drug_cost_per_person_year, fill = "Kela"), alpha = 0.5, bins = 100) +
	geom_histogram(data = total_new, mapping = aes(x = log_hilmo_cost_per_person_year, fill = "HILMO"), alpha = 0.5, bins = 100) +
	geom_histogram(data = total_new, mapping = aes(x = log_avohilmo_cost_per_person_year, fill = "AvoHILMO"), alpha = 0.5, bins = 100) +
	labs(x = "log(Annual Healthcare Cost)", y = "Frequency", fill = "") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

##################################################
########## CONVERT CPRA TO RSID
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

sumstat_date = "20220217"
sumstat_name = "log_total_cost_per_person_year"
sumstat = fread(file = paste0(path, "gwas/", sumstat_date, "/", sumstat_name, ".pheweb"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
colnames(sumstat) = c("chr","pos","ref","alt","pval","beta","se","maf","n_hom_cases","n_het_cases","n_hom_controls","n_het_controls")
sumstat$chr_pos_hg38 = paste0(sumstat$chr, "_", sumstat$pos)
#rsid = fread(file = paste0(path, "MR/finngen_rsids.tsv"), sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#colnames(rsid) = c("rsid","chr_pos_hg38", "SNP")
rsid_cpra = fread(file = paste0(path, "gwas/finngen_rsids_cpra.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
sumstat_new = merge(rsid_cpra, sumstat, by = "chr_pos_hg38", all.y = TRUE)
sumstat_new = sumstat_new[,-3]
fwrite(sumstat_new, paste0(path, "gwas/", sumstat_date, "/", sumstat_name, "_", sumstat_date, ".txt"), sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)


finngen-cli request-workflow --wdl ./saige.wdl --dependencies ./saige_sub.wdl.zip --input ./saige.json

##################################################
########## CALCULATE HERITABILITY 
##################################################

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 w_hm3.snplist.bz2
tar -jxvf eur_w_ld_chr.tar.bz2

# go to ubuntu
# conda activate ldsc
# cd ldsc

python munge_sumstats.py \
--sumstats /mnt/c/Users/jiwoo/Desktop/gwas/20220222/total_cost_per_person_year_20220217.txt \
--N 352080 \
--snp rsid \
--a1 ref \
--a2 alt \
--p pval \
--a1-inc \
--out /mnt/c/Users/jiwoo/Desktop/gwas/20220222/total_cost_per_person_year_20220217 \
--merge-alleles /mnt/c/Users/jiwoo/Desktop/gwas/w_hm3.snplist

python ldsc.py \
--h2 /mnt/c/Users/jiwoo/Desktop/gwas/20220222/total_cost_per_person_year_20220217.sumstats.gz \
--ref-ld-chr /mnt/c/Users/jiwoo/Desktop/gwas/eur_w_ld_chr/ \
--w-ld-chr /mnt/c/Users/jiwoo/Desktop/gwas/eur_w_ld_chr/ \
--out /mnt/c/Users/jiwoo/Desktop/gwas/20220222/total_cost_per_person_year_20220217_h2

#uniedit munge sumstats

##################################################
########## DO MENDELIAN RANDOMIZATION	
##################################################

path = "C:/Jiwoo Lee/Healthcare_Cost_Research_2021/"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}
if(!require(ggplot2)) {install_github("ggplot2"); library(ggplot2)}

setDTthreads(0)

# List available GWAS
#ao <- available_outcomes()
#ao_df <- as.data.frame(ao)

exposure_id = "ukb-a-360"
exposure_id = c(
"ieu-a-299",
"ieu-a-300",
"ieu-a-302",
"ieu-a-835",
"ieu-b-25",
"ieu-b-73",
"ukb-a-382",
"ukb-d-30600_irnt",
"ukb-d-30620_irnt",
"ukb-d-30700_irnt",
"ukb-d-30710_irnt",
"ukb-d-30720_irnt",
"ukb-d-30740_irnt",
"ukb-d-30890_irnt",
"ukb-a-360"
)	

sumstat_name = "total_cost"
sumstat = fread(paste0("C:/Users/jiwoo/Desktop/gwas/", sumstat_name, ".txt"), header = TRUE, stringsAsFactors = FALSE, verbose = TRUE)

mr_res = NULL
mr_res_mini = NULL
for (i in 1:length(exposure_id)) {
	print(exposure_id[i])
	# Get instruments
	print("Loading exposure data...")
	exposure_dat = extract_instruments(outcomes = exposure_id[i])
	# Get effects of instruments on outcome
	print("Loading outcome data...")
	outcome_temp = sumstat[which(sumstat$rsid %in% exposure_dat$SNP),]
	fwrite(outcome_temp, paste0(path, "MR/temp.csv"))
	outcome_dat = read_outcome_data(snps = NULL,
		filename = paste0(path, "MR/temp.csv"),
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
	print("Harmonizing data...")
	dat = harmonise_data(exposure_dat, outcome_dat)
	#dat_mini = dat[order(dat$pval.exposure),]
	#dat_mini = dat_mini[1:10,]
	# Perform MR
	print("Performing MR...")
	res = mr(dat, method_list = c("mr_ivw")) 
	#res = mr(dat) 
	#res_mini = mr(dat_mini, method_list = c("mr_ivw"))
	mr_res = rbind(mr_res, res)
	#mr_res_mini = rbind(mr_res_mini, res_mini)
}
mr_res = as.data.frame(mr_res)
mr_res$b = as.numeric(mr_res$b)
mr_res$se = as.numeric(mr_res$se)
mr_res$pval = as.numeric(mr_res$pval)
mr_temp = mr_res
mr_temp$percent = round((exp(mr_temp$b) - 1) * 100, 2)
mr_temp$confint_lower = round((exp(mr_temp$b - 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint_upper = round((exp(mr_temp$b + 1.96 * mr_temp$se) - 1) * 100, 2)
mr_temp$confint = paste0("[", mr_temp$confint_lower, ", ", mr_temp$confint_upper, "]")
mr_temp = mr_temp[order(mr_temp$pval),]
fwrite(mr_temp, paste0(path, "MR/", sumstat_name, "_irnt_20220204.csv"))

#mr_res_mini = as.data.frame(mr_res_mini)
#mr_res_mini$b = as.numeric(mr_res_mini$b)
#mr_res_mini$se = as.numeric(mr_res_mini$se)
#mr_res_mini$pval = as.numeric(mr_res_mini$pval)
#mr_temp = mr_res_mini
#mr_temp$percent = round((exp(mr_temp$b) - 1) * 100, 2)
#mr_temp$confint_lower = round((exp(mr_temp$b - 1.96 * mr_temp$se) - 1) * 100, 2)
#mr_temp$confint_upper = round((exp(mr_temp$b + 1.96 * mr_temp$se) - 1) * 100, 2)
#mr_temp$confint = paste0("[", mr_temp$confint_lower, ", ", mr_temp$confint_upper, "]")
#mr_temp = mr_temp[order(mr_temp$pval),]
#fwrite(mr_temp, paste0(path, "MR/", sumstat_name, "_irnt_mini.csv"))

my_scatter_plot <- function(mr_results, dat, x_label, y_label){
	dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
	requireNamespace("ggplot2", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
	{
		d <- plyr::mutate(d)
		if(nrow(d) < 2 | sum(d$mr_keep) == 0)
		{
			return(blank_plot("Insufficient number of SNPs"))
		}
		d <- subset(d, mr_keep)
		index <- d$beta.exposure < 0
		d$beta.exposure[index] <- d$beta.exposure[index] * -1
		d$beta.outcome[index] <- d$beta.outcome[index] * -1
		mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
		mrres$a <- 0
		if("MR Egger" %in% mrres$method)
		{
			temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
			mrres$a[mrres$method == "MR Egger"] <- temp$b_i
		}

		if("MR Egger (bootstrap)" %in% mrres$method)
		{
			temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
			mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
		}
		#g = gridExtra::tableGrob(mr_results[,c("method", "b", "se", "pval")], theme = ttheme_default(base_size = 7))
		ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
				ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
				ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
				ggplot2::geom_point() +
				ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
				ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
				ggplot2::labs(colour="MR Test", x=paste0("SNP Effect on ", x_label), y=paste0("SNP Effect on ", y_label)) +
				ggplot2::theme(legend.position="top", legend.direction="vertical", panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black")) + 
				ggplot2::guides(colour=ggplot2::guide_legend(ncol=2)) #+
				#ggplot2::annotation_custom(g, xmin=max(d$beta.exposure)*0.5, ymin=max(d$beta.outcome)*0.75)
	})
	mrres}
my_scatter_plot(res, dat, "Systolic Blood Pressure", "Log(Healthcare Cost)")

ggplot() +
	geom_point(data = res, mapping = aes(x = reorder(method, -b), y = b), size = 3) +
	coord_flip () +
	geom_errorbar(data = res, mapping = aes(x = reorder(method, -b), ymin = b - 1.96 * se, ymax = b + 1.96 * se), width = 0.1, size = 1) +
	scale_x_discrete(limits = res$method[order(-res$b)], labels = res$method[order(-res$b)]) + 
	#scale_y_continuous(breaks = seq(floor(min(mod_temp$beta_lower)), ceiling(max(mod_temp$beta_upper)), 1), limits = c(floor(min(mod_temp$beta_lower)), ceiling(max(mod_temp$beta_upper)))) + 
	labs(x = "", y = "Beta", color = "Significant?") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

##################################################
########## PLOT MENDELIAN RANDOMIZATION	
##################################################

path = "C:/Jiwoo Lee/Healthcare_Cost_Research_2021/"

if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(TwoSampleMR)) {install_github("TwoSampleMR"); library(TwoSampleMR)}
if(!require(ggplot2)) {install_github("ggplot2"); library(ggplot2)}

# MAKE FOREST PLOT

total = fread(paste0(path, "MR/total_cost_irnt_20220204.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE)
mod_temp = total
ggplot() +
	geom_point(data = mod_temp, mapping = aes(x = reorder(exposure, -b), y = b, color = ifelse(pval < 0.05, "Yes", "No")), size = 3) +
	coord_flip () +
	geom_errorbar(data = mod_temp, mapping = aes(x = reorder(exposure, -b), ymin = b - 1.96 * se, ymax = b + 1.96 * se, color = ifelse(pval < 0.05, "Yes", "No")), width = 0.5, size = 2) +
	scale_x_discrete(limits = mod_temp$exposure[order(-mod_temp$b)], labels = mod_temp$exposure[order(-mod_temp$b)]) + 
	#scale_y_continuous(breaks = seq(floor(min(mod_temp$beta_lower)), ceiling(max(mod_temp$beta_upper)), 1), limits = c(floor(min(mod_temp$beta_lower)), ceiling(max(mod_temp$beta_upper)))) + 
	labs(x = "", y = "Beta (IVW)", color = "Significant?") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

log_total = fread(paste0(path, "MR/log_total_cost_irnt_20220204.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE)
mod_temp = log_total
ggplot() +
	geom_point(data = mod_temp, mapping = aes(x = reorder(exposure, -percent), y = percent, color = ifelse(pval < 0.05, "Yes", "No")), size = 3) +
	coord_flip () +
	geom_errorbar(data = mod_temp, mapping = aes(x = reorder(exposure, -percent), ymin = confint_lower, ymax = confint_upper, color = ifelse(pval < 0.05, "Yes", "No")), width = 0.5, size = 2) +
	scale_x_discrete(limits = mod_temp$exposure[order(-mod_temp$percent)], labels = mod_temp$exposure[order(-mod_temp$percent)]) + 
	#scale_y_continuous(breaks = seq(floor(min(mod_temp$beta_lower)), ceiling(max(mod_temp$beta_upper)), 1), limits = c(floor(min(mod_temp$beta_lower)), ceiling(max(mod_temp$beta_upper)))) + 
	labs(x = "", y = "Percent Change in Healthcare Cost (IVW)", color = "Significant?") + 
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

# MAKE FOREST PLOT WITH CAD EXCLUSION

biomarker = "ieu-a-300"
file_list = c("total_cost", "log_total_cost", "total_cost_20220105", "log_total_cost_20220105",
	"drug_cost", "log_drug_cost", "drug_cost_20220105", "log_drug_cost_20220105",
	"hilmo_cost", "log_hilmo_cost", "hilmo_cost_20220105", "log_hilmo_cost_20220105",
	"avohilmo_cost", "log_avohilmo_cost", "avohilmo_cost_20220105", "log_avohilmo_cost_20220105")
cost_list = c("Total", "Total", "Total", "Total",
	"Drug", "Drug", "Drug", "Drug",
	"HILMO", "HILMO", "HILMO", "HILMO",
	"AvoHILMO", "AvoHILMO", "AvoHILMO", "AvoHILMO")
log_list = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1) 
cad_list = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)

final = NULL
for(i in 1:length(file_list)) {
	sumstat = fread(paste0(path, "MR/", file_list[i], "_irnt.csv"), header = TRUE, stringsAsFactors = FALSE, verbose = FALSE)
	final = rbind(final, cbind(sumstat[which(sumstat$id.exposure == biomarker)], "file" = file_list[i], "type" = cost_list[i], "log" = log_list[i], "cad" = cad_list[i]))
	print(cost_list[i])
}
final$percent = as.numeric(final$percent)
final$confint_lower = as.numeric(final$confint_lower)
final$confint_upper = as.numeric(final$confint_upper)
final$label = ifelse(final$cad == 0, paste0(final$type, ", CAD included"), paste0(final$type, ", CAD excluded"))
final_temp = final[which(final$method == "Inverse variance weighted"),]

ggplot() +
	geom_line(data = final_temp[which(final_temp$log == 1),], mapping = aes(x = label, y = percent, color = cad), size = 1) +
	geom_point(data = final_temp[which(final_temp$log == 1),], mapping = aes(x = label, y = percent, color = cad), size = 3) +
	geom_errorbar(data = final_temp[which(final_temp$log == 1),], mapping = aes(x = label, ymin = confint_lower, ymax = confint_upper, color = cad), width = 0.1, size = 1) +
	coord_flip() + 
  labs(title = final_temp$exposure[1], x = "", y = "Percent Change in Healthcare Cost", color = "") + 
  guides(color = FALSE) +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

ggplot() +
	geom_line(data = final_temp[which(final_temp$log == 0),], mapping = aes(x = label, y = b, color = cad), size = 1) +
	geom_point(data = final_temp[which(final_temp$log == 0),], mapping = aes(x = label, y = b, color = cad), size = 3) +
	geom_errorbar(data = final_temp[which(final_temp$log == 0),], mapping = aes(x = label, ymin = b - 1.96*se, ymax = b + 1.96*se, color = cad), width = 0.1, size = 1) +
	coord_flip() + 
  labs(title = final_temp$exposure[1], x = "", y = "Beta +/- 95% CI", color = "") + 
  guides(color = FALSE) +
	theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15), panel.grid.major = element_line(size = 0.1, color = "gray"), panel.grid.minor = element_line(size = 0.1, color = "gray"), panel.background = element_blank(), strip.background = element_rect(fill = "gray95"), axis.line = element_line(colour = "black"))

##################################################
########## STANDARDIZE RSID AND CPRA AND GENOME BUILD
##################################################

cd ~/jiwoo/healthcare_cost_repository

R

path = "/home/ivm/jiwoo/healthcare_cost_repository/"

# LOAD PACKAGES
if(!require(data.table)) {install.packages("data.table"); library(data.table)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

variants = fread(file = paste0(path, "gwas/bmi_variants_sig.txt"), header = FALSE, stringsAsFactors = FALSE)
colnames(variants) = c("rsid")

#rsid = fread(file = paste0(path, "MR/finngen_rsids.tsv"), sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#colnames(rsid) = c("rsid","chr_pos_hg38", "SNP")
rsid_cpra = fread(file = paste0(path, "gwas/finngen_rsids_cpra.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)

dim(variants) # 37
variants = merge(variants, rsid_cpra, by = "rsid", all.x = TRUE)
dim(variants) # 37

variants$allele = substr(variants$SNP, nchar(variants$SNP), nchar(variants$SNP))
fwrite(list(variants$SNP), paste0(path, "gwas/bmi_variants_sig.csv"), col.names = FALSE)
fwrite(variants[,c("SNP", "allele")], paste0(path, "gwas/bmi_variants_sig_alt.csv"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

weights = fread(file = paste0(path, "gwas/bmi_rsid_weights_sig.txt"), header = FALSE, stringsAsFactors = FALSE)
colnames(weights) = c("rsid", "allele", "weight")
weights = merge(weights, rsid_cpra, by = "rsid", all.x = TRUE)
fwrite(weights[,c("SNP", "allele", "weight")], paste0(path, "gwas/bmi_cpra_weights_sig.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

##################################################
########## EXTRACT VARIANTS
##################################################

plink --bfile /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8 \
--extract bmi_variants_sig.csv \
--recode-allele bmi_variants_sig_alt.csv \
--recode A --threads 16 --out bmi_variants_sig

##################################################
########## CALCULATE ALLELE FREQUENCIES
##################################################

plink --bfile /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8 \
--extract sbp_variants_all.csv \
--freq 

##################################################
########## CALCULATE ALLELE SCORES
##################################################

plink --bfile /finngen/library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8 \
--score bmi_cpra_weights_sig.txt center
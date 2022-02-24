setwd("/home/jsjukara/costs/")
library(data.table)
library(tidyverse)
library(lubridate)

# load HILMO unit costs table scraped from the THL report and do some wrangling
unit_costs_hilmo <- read_delim("unit_costs/unit_costs_hilmo_21_4_2020.csv", delim=";")[2:7]
names(unit_costs_hilmo) <- tolower(names(unit_costs_hilmo))
names(unit_costs_hilmo)[c(4,5)] <- c("hospital_type", "pala_class")
unit_costs_hilmo <- unit_costs_hilmo[!is.na(unit_costs_hilmo$ea),]
unit_costs_hilmo$unit_cost <- as.numeric(unit_costs_hilmo$unit_cost)
unit_costs_hilmo$ea <- as.character(unit_costs_hilmo$ea)

startyear <- 1972
endyear <- 2019

# load custom more detailed longitudinal phenotype file for FinnGen and do some wrangling
d <- fread("zcat finngen_R5_v3_custom_detailed_longitudinal_allcodes_REFINERY_ONLY.gz")
dnames <- names(d)
d$year <- decimal_date(as.Date(d$APPROX_EVENT_DAY))
d <- d[year >= startyear]
d <- d[year < endyear]


d$rowindex <- 1:nrow(d)
dhilmo <- d[SOURCE %in% c("INPAT", "OUTPAT")]
dprim <- d[SOURCE == "PRIM_OUT"]

# modify specialty
dhilmo$CODE6 <- as.integer(gsub("[^0-9.-]", "", dhilmo$CODE6))
dhilmo$year <- decimal_date(as.Date(dhilmo$APPROX_EVENT_DAY))
dhilmo <- dhilmo[year >= startyear]
dhilmo <- dhilmo[CATEGORY=="M"]
names(dhilmo) <- tolower(names(dhilmo))
names(dhilmo)[c(5:11)] <- c("pdgo", "pdge", "atc", "hospdays", "pala", "ea", "hospital_type")
# make specialties that are missing from unit costs table "-1"
missing_eas <- unique(dhilmo$ea)[!(unique(dhilmo$ea)%in%unique(unit_costs_hilmo$ea))]
dhilmo$ea[dhilmo$ea %in% missing_eas] <- "-1"
dhilmo$hospdays <- dhilmo$hospdays+1


# impute some missing values for the unit costs 
unit_costs_hilmo$thl_days[is.na(unit_costs_hilmo$thl_days)] <- 4.9
day_ward_costs <- unit_costs_hilmo %>% filter(pala_class=="ward")
day_ward_costs$thl_days[day_ward_costs$ea==15] <- mean(dhilmo$hospdays[dhilmo$ea==15 & dhilmo$source=="INPAT"])

day_ward_costs$hilmo_days_mean <- NA
day_ward_costs$hilmo_days_median <- NA
for (i in 1:nrow(day_ward_costs)) {
    ea <- day_ward_costs$ea[i]
    day_ward_costs$hilmo_days_mean[i] <- mean(dhilmo$hospdays[dhilmo$source=="INPAT" & dhilmo$ea==ea & dhilmo$year>=1998])
    day_ward_costs$hilmo_days_median[i] <- median(dhilmo$hospdays[dhilmo$source=="INPAT" & dhilmo$ea==ea & dhilmo$year>=1998])
}

# assuming 408€ per incremental ward episode day, estimate costs per episode based on a model of:
# costs = fixed costs + 408 * ward_days
day_ward_costs$cost_fixed <- day_ward_costs$unit_cost - 408 * day_ward_costs$thl_days
day_ward_costs$cost_fixed <- ifelse(day_ward_costs$cost_fixed<0, 0, day_ward_costs$cost_fixed)
day_ward_costs$cost_per_thlday <- day_ward_costs$unit_cost/day_ward_costs$thl_days
day_ward_costs$cost_per_day_wofixed <- (day_ward_costs$unit_cost-day_ward_costs$cost_fixed)/day_ward_costs$thl_days


pala_outpatient <- c("83", "92", "93", "94") #päiväsairaanhoito, ajanvaraus, ajanvaraus(uusintak.), konsultaatio
pala_ward <- c("1", "2", "5", "6", "7", "8")

dhilmo <- dhilmo %>%
    mutate(pala_class = case_when(pala == "91" & source == "OUTPAT" ~ "emergency",
                                 source == "OUTPAT" ~ "outpatient",
                                 source == "INPAT" ~ "ward"))

# join unit costs
dhilmo <- left_join(dhilmo, unit_costs_hilmo[unit_costs_hilmo$hospital_type=="all",c("ea", "pala_class", "unit_cost")], by=c("ea", "pala_class"))
# join fixed costs and thl_days
dhilmo <- left_join(dhilmo, day_ward_costs[,c("ea", "pala_class", "cost_fixed", "thl_days")],by=c("ea", "pala_class"))

# adjust ward episode unit cost based on the model of fixed costs + costs/day
dhilmo$unit_cost_w_fixed <- ifelse(dhilmo$pala_class == "ward",
                                   dhilmo$cost_fixed + (dhilmo$unit_cost - dhilmo$cost_fixed) * (dhilmo$hospdays / dhilmo$thl_days),
                                   dhilmo$unit_cost)

# for some specialties, instead use specific costs/day estimates from the report
dhilmo <- dhilmo %>%
    mutate(unit_cost_w_fixed = case_when(ea == "98" & pala_class == "ward" & hospdays<=90 ~hospdays*234,
                                ea == "98" & pala_class == "ward" & hospdays>90 ~ 90*234 + (hospdays-90)*191,
                                ea == "70" & pala_class == "ward" & hospdays<=90 ~hospdays*408,
                                ea == "70" & pala_class == "ward" & hospdays>90 ~ 90*408 + (hospdays-90)*408*0.816,
                                ea == "74" & pala_class == "ward" & hospdays<=90 ~hospdays*568,
                                ea == "74" & pala_class == "ward" & hospdays>90 ~ 90*568 + (hospdays-90)*568*0.816,
                                ea == "75" & pala_class == "ward" & hospdays<=90 ~hospdays*638,
                                ea == "75" & pala_class == "ward" & hospdays>90 ~ 90*638 + (hospdays-90)*638*0.816,
                                TRUE ~ unit_cost_w_fixed),
          unit_cost = case_when(ea == "98" & pala_class == "ward" & hospdays<=90 ~hospdays*234,
                                ea == "98" & pala_class == "ward" & hospdays>90 ~ 90*234 + (hospdays-90)*191,
                                ea == "70" & pala_class == "ward" & hospdays<=90 ~hospdays*408,
                                ea == "70" & pala_class == "ward" & hospdays>90 ~ 90*408 + (hospdays-90)*408*0.816,
                                ea == "74" & pala_class == "ward" & hospdays<=90 ~hospdays*568,
                                ea == "74" & pala_class == "ward" & hospdays>90 ~ 90*568 + (hospdays-90)*568*0.816,
                                ea == "75" & pala_class == "ward" & hospdays<=90 ~hospdays*638,
                                ea == "75" & pala_class == "ward" & hospdays>90 ~ 90*638 + (hospdays-90)*638*0.816,
                                TRUE ~ unit_cost))


#convert 2011 euro estimates to 2018 euros using price index of public expenditure
dhilmo$unit_cost <- dhilmo$unit_cost/0.9341038

names(dhilmo)[1] <- "FINNGENID"

# load additional phenotype file to get year of birth and death
 
endbig <- fread("/home/jsjukara/finngen_R5_V2_endpoint.txt", data.table = FALSE)

keepcols_endbig <- c("FINNGENID", "BL_AGE", "BL_YEAR", "SEX",
                     'DEATH','DEATH_AGE','DEATH_YEAR')

fg_df <- endbig[,keepcols_endbig]


fg_df$YEAR_OF_BIRTH <- fg_df$BL_YEAR - fg_df$BL_AGE

# make end of follow-up be death or 2019 if did not die during follow-up
fg_df$FU_END <- ifelse(is.na(fg_df$DEATH_YEAR), 2019, fg_df$DEATH_YEAR)

## Calculate individual sums of costs for overall data
# calculate sum of costs by id for inpatient and outpatient data from HILMO after 1998
indiv_costs_post1998 <- dhilmo %>%
    filter(year >= 1998) %>%
    group_by(FINNGENID) %>%
    summarize(HILMO_COSTS_POST1998=sum(unit_cost))

# calculate sum of costs by id for inpatient data only from HILMO for whole follow-up
indiv_costs_inpatient <- dhilmo %>%
    filter(source == "INPAT") %>%
    group_by(FINNGENID) %>%
    summarize(HILMO_COSTS_INPATIENT=sum(unit_cost))

# create data frame containing all ids
indiv_costs <- data.frame(FINNGENID = unique(fg_df$FINNGENID))

# join costs
indiv_costs <- left_join(indiv_costs, indiv_costs_post1998)

# make post1998 costs for IDs that died before 1998 NA, make NA costs for others 0.
ids_pre1998 <- fg_df$FINNGENID[fg_df$DEATH_YEAR<=1998]
indiv_costs$HILMO_COSTS_POST1998 <- ifelse(indiv_costs$FINNGENID %in% ids_pre1998, NA,
                                          ifelse(is.na(indiv_costs$HILMO_COSTS_POST1998), 0, indiv_costs$HILMO_COSTS_POST1998))


indiv_costs <- left_join(indiv_costs, indiv_costs_inpatient)
                  
# make inpatient costs for IDs that died before 1972 NA, make NA costs for others 0.
ids_pre1972 <- fg_df$FINNGENID[fg_df$DEATH_YEAR<=1972]
indiv_costs$HILMO_COSTS_INPATIENT <- ifelse(indiv_costs$FINNGENID %in% ids_pre1972, NA,
                                          ifelse(is.na(indiv_costs$HILMO_COSTS_INPATIENT), 0, indiv_costs$HILMO_COSTS_INPATIENT))

# save costs
fwrite(indiv_costs, file="hilmo_costs/data/finngen_R5_overall_hilmo_costs.txt.gz")


## Calculate sums of costs by age group
# join year of birth to hilmo and calculate age at event
dhilmo <- left_join(dhilmo, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH")], by="FINNGENID")

dhilmo$age_at_event <- dhilmo$year - dhilmo$YEAR_OF_BIRTH



dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          source == "INPAT") %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c("FINNGENID", sort(names(dhilmo_age[2:10])))]

# make data frame with all FINNGENID's, even those without any costs in the data
dhilmo_age_full <- data.frame(FINNGENID = unique(fg_df$FINNGENID))
dhilmo_age_full <- left_join(dhilmo_age_full, dhilmo_age, by="FINNGENID")
dhilmo_age_full <- left_join(dhilmo_age_full, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

dhilmo_age <- dhilmo_age_full

age_low <- seq(0, 80, 10)

# convert temporarily to matrix for fast assignment (at least 1000x faster)
dhilmo_agem <- as.matrix(dhilmo_age)

# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:nrow(dhilmo_age)) {  
    age_start <- 1972 - dhilmo_age$YEAR_OF_BIRTH[i]  # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2019, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2019
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]

    # set NA costs to 0 if follow-up started before upper bound of age interval
    # and follow-up ended after the lower bound of age interval
    indices <- which(age_start < age_low + 10 & age_end > age_low & is.na(dhilmo_agem[i,2:10])) + 1
    dhilmo_agem[i,indices] <- 0
    
    # set costs to NA for age interval if follow-up started after upper bound of age interval
    indices <- which(age_start > age_low + 10) + 1
    dhilmo_agem[i,indices] <- NA
}

# convert back to data frame
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)
fwrite(dhilmo_age[,1:10], file="hilmo_costs/data/finngen_R5_hilmo_inpatient_costs_by_age.txt.gz")



dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          year >= 1998) %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c("FINNGENID", sort(names(dhilmo_age[2:10])))]

# make data frame with all FINNGENID's, even those without any costs in the data
dhilmo_age_full <- data.frame(FINNGENID = unique(fg_df$FINNGENID))
dhilmo_age_full <- left_join(dhilmo_age_full, dhilmo_age, by="FINNGENID")
dhilmo_age_full <- left_join(dhilmo_age_full, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

dhilmo_age <- dhilmo_age_full

# convert temporarily to matrix for fast assignment (at least 1000x faster)
dhilmo_agem <- as.matrix(dhilmo_age)

# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:nrow(dhilmo_age)) {  
    age_start <- 1998 - dhilmo_age$YEAR_OF_BIRTH[i]  # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2019, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2019
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]

    # set NA costs to 0 if follow-up started before upper bound of age interval
    # and follow-up ended after the lower bound of age interval
    indices <- which(age_start < age_low + 10 & age_end > age_low & is.na(dhilmo_agem[i,2:10])) + 1
    dhilmo_agem[i,indices] <- 0
    
    # set costs to NA for age interval if follow-up started after upper bound of age interval
    indices <- which(age_start > age_low + 10) + 1
    dhilmo_agem[i,indices] <- NA
}
# convert back to data frame
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)
head(dhilmo_age)
fwrite(dhilmo_age[,1:10], file="hilmo_costs/data/finngen_R5_hilmo_costs_by_age_post1998.txt.gz")

d1 <- fread("zcat hilmo_costs/finngen_R5_overall_hilmo_costs.txt.gz", data.table=FALSE)
d2 <- fread("zcat hilmo_costs/finngen_R5_hilmo_inpatient_costs_by_age.txt.gz", data.table=FALSE)
d3 <- fread("zcat hilmo_costs/finngen_R5_hilmo_costs_by_age_post1998.txt.gz", data.table=FALSE)

head(d1)
head(d2)
head(d3)

str(d1)
str(d2)
str(d3)

hist(dhilmo$year)

sum(is.na(d2[,2:10]), na.rm=T)
sum(is.na(d3[,2:10]), na.rm=T)

dhilmo %>% filter(FINNGENID == "FG2224RB9B")

# load additional phenotype file to get year of birth and death
 
endbig <- fread("/home/jsjukara/finngen_R5_V2_endpoint.txt", data.table = FALSE)

keepcols_endbig <- c("FINNGENID", "BL_AGE", "BL_YEAR", "SEX",
                     'DEATH','DEATH_AGE','DEATH_YEAR')

fg_df <- endbig[,keepcols_endbig]


fg_df$YEAR_OF_BIRTH <- fg_df$BL_YEAR - fg_df$BL_AGE

# make end of follow-up be death or 2017 if did not die during follow-up
fg_df$FU_END <- ifelse(is.na(fg_df$DEATH_YEAR), 2017, fg_df$DEATH_YEAR)

## Calculate individual sums of costs for overall data
# calculate sum of costs by id for inpatient and outpatient data from HILMO after 1998
indiv_costs_post1998 <- dhilmo %>%
    filter(year >= 1998) %>%
    group_by(FINNGENID) %>%
    summarize(HILMO_COSTS_POST1998=sum(unit_cost))

# calculate sum of costs by id for inpatient data only from HILMO for whole follow-up
indiv_costs_inpatient <- dhilmo %>%
    filter(source == "INPAT") %>%
    group_by(FINNGENID) %>%
    summarize(HILMO_COSTS_INPATIENT=sum(unit_cost))

# create data frame containing all ids
indiv_costs <- data.frame(FINNGENID = unique(fg_df$FINNGENID))

# join costs
indiv_costs <- left_join(indiv_costs, indiv_costs_post1998)

# make post1998 costs for IDs that died before 1998 NA, make NA costs for others 0.
ids_pre1998 <- fg_df$FINNGENID[fg_df$DEATH_YEAR<=1998]
indiv_costs$HILMO_COSTS_POST1998 <- ifelse(indiv_costs$FINNGENID %in% ids_pre1998, NA,
                                          ifelse(is.na(indiv_costs$HILMO_COSTS_POST1998), 0, indiv_costs$HILMO_COSTS_POST1998))


indiv_costs <- left_join(indiv_costs, indiv_costs_inpatient)
                  
# make inpatient costs for IDs that died before 1972 NA, make NA costs for others 0.
ids_pre1972 <- fg_df$FINNGENID[fg_df$DEATH_YEAR<=1972]
indiv_costs$HILMO_COSTS_INPATIENT <- ifelse(indiv_costs$FINNGENID %in% ids_pre1972, NA,
                                          ifelse(is.na(indiv_costs$HILMO_COSTS_INPATIENT), 0, indiv_costs$HILMO_COSTS_INPATIENT))

# save costs
fwrite(indiv_costs, file="hilmo_costs/finngen_R5_overall_hilmo_costs.txt.gz")


## Calculate sums of costs by age group
# join year of birth to hilmo and calculate age at event
dhilmo <- left_join(dhilmo, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH")], by="FINNGENID")

dhilmo$age_at_event <- dhilmo$year - dhilmo$YEAR_OF_BIRTH



dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          source == "INPAT") %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c("FINNGENID", sort(names(dhilmo_age[2:10])))]

# make data frame with all FINNGENID's, even those without any costs in the data
dhilmo_age_full <- data.frame(FINNGENID = unique(fg_df$FINNGENID))
dhilmo_age_full <- left_join(dhilmo_age_full, dhilmo_age, by="FINNGENID")
dhilmo_age_full <- left_join(dhilmo_age_full, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

dhilmo_age <- dhilmo_age_full

age_low <- seq(0, 80, 10)

# convert temporarily to matrix for fast assignment (at least 1000x faster)
dhilmo_agem <- as.matrix(dhilmo_age)

# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:nrow(dhilmo_age)) {  
    age_start <- 1972 - dhilmo_age$YEAR_OF_BIRTH[i]  # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2017, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2017
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]

    # set NA costs to 0 if follow-up started before upper bound of age interval
    # and follow-up ended after the lower bound of age interval
    indices <- which(age_start < age_low + 10 & age_end > age_low & is.na(dhilmo_agem[i,2:10])) + 1
    dhilmo_agem[i,indices] <- 0
    
    # set costs to NA for age interval if follow-up started after upper bound of age interval
    indices <- which(age_start > age_low + 10) + 1
    dhilmo_agem[i,indices] <- NA
}

# convert back to data frame
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)
fwrite(dhilmo_age, file="hilmo_costs/finngen_R5_hilmo_inpatient_costs_by_age.txt.gz")



dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          year >= 1998) %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c("FINNGENID", sort(names(dhilmo_age[2:10])))]

# make data frame with all FINNGENID's, even those without any costs in the data
dhilmo_age_full <- data.frame(FINNGENID = unique(fg_df$FINNGENID))
dhilmo_age_full <- left_join(dhilmo_age_full, dhilmo_age, by="FINNGENID")
dhilmo_age_full <- left_join(dhilmo_age_full, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

dhilmo_age <- dhilmo_age_full

age_low <- seq(0, 80, 10)

# convert temporarily to matrix for fast assignment (at least 1000x faster)
dhilmo_agem <- as.matrix(dhilmo_age)

# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:nrow(dhilmo_age)) {  
    age_start <- 1972 - dhilmo_age$YEAR_OF_BIRTH[i]  # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2017, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2017
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]

    # set NA costs to 0 if follow-up started before upper bound of age interval
    # and follow-up ended after the lower bound of age interval
    indices <- which(age_start < age_low + 10 & age_end > age_low & is.na(dhilmo_agem[i,2:10])) + 1
    dhilmo_agem[i,indices] <- 0
    
    # set costs to NA for age interval if follow-up started after upper bound of age interval
    indices <- which(age_start > age_low + 10) + 1
    dhilmo_agem[i,indices] <- NA
}
# convert back to data frame
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)
head(dhilmo_age)
fwrite(dhilmo_age, file="hilmo_costs/finngen_R5_hilmo_costs_by_age_post1998.txt.gz")

# load additional phenotype file to get year of birth and death
 
endbig <- fread("/home/jsjukara/finngen_R5_V2_endpoint.txt", data.table = FALSE)

keepcols_endbig <- c("FINNGENID", "BL_AGE", "BL_YEAR", "SEX",
                     'DEATH','DEATH_AGE','DEATH_YEAR')

fg_df <- endbig[,keepcols_endbig]


fg_df$YEAR_OF_BIRTH <- fg_df$BL_YEAR - fg_df$BL_AGE

# make end of follow-up be death or 2017 if did not die during follow-up
fg_df$FU_END <- ifelse(is.na(fg_df$DEATH_YEAR), 2017, fg_df$DEATH_YEAR)

## Calculate individual sums of costs for overall data
# calculate sum of costs by id for inpatient and outpatient data from HILMO after 1998
indiv_costs_post1998 <- dhilmo %>%
    filter(year >= 1998) %>%
    group_by(FINNGENID) %>%
    summarize(HILMO_COSTS_POST1998=sum(unit_cost))

# calculate sum of costs by id for inpatient data only from HILMO for whole follow-up
indiv_costs_inpatient <- dhilmo %>%
    filter(source == "INPAT") %>%
    group_by(FINNGENID) %>%
    summarize(HILMO_COSTS_INPATIENT=sum(unit_cost))

# create data frame containing all ids
indiv_costs <- data.frame(FINNGENID = unique(fg_df$FINNGENID))

# join costs
indiv_costs <- left_join(indiv_costs, indiv_costs_post1998)

# make post1998 costs for IDs that died before 1998 NA, make NA costs for others 0.
ids_pre1998 <- fg_df$FINNGENID[fg_df$DEATH_YEAR<=1998]
indiv_costs$HILMO_COSTS_POST1998 <- ifelse(indiv_costs$FINNGENID %in% ids_pre1998, NA,
                                          ifelse(is.na(indiv_costs$HILMO_COSTS_POST1998), 0, indiv_costs$HILMO_COSTS_POST1998))


indiv_costs <- left_join(indiv_costs, indiv_costs_inpatient)
                  
# make inpatient costs for IDs that died before 1972 NA, make NA costs for others 0.
ids_pre1972 <- fg_df$FINNGENID[fg_df$DEATH_YEAR<=1972]
indiv_costs$HILMO_COSTS_INPATIENT <- ifelse(indiv_costs$FINNGENID %in% ids_pre1972, NA,
                                          ifelse(is.na(indiv_costs$HILMO_COSTS_INPATIENT), 0, indiv_costs$HILMO_COSTS_INPATIENT))

# save costs
fwrite(indiv_costs, file="hilmo_costs/finngen_R5_overall_hilmo_costs.txt.gz")


## Calculate sums of costs by age group
# join year of birth to hilmo and calculate age at event
dhilmo <- left_join(dhilmo, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH")], by="FINNGENID")

dhilmo$age_at_event <- dhilmo$year - dhilmo$YEAR_OF_BIRTH




dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          source == "INPAT") %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c("FINNGENID", sort(names(dhilmo_age[2:10])))]

# make data frame with all FINNGENID's, even those without any costs in the data
dhilmo_age_full <- data.frame(FINNGENID = unique(fg_df$FINNGENID))
dhilmo_age_full <- left_join(dhilmo_age_full, dhilmo_age, by="FINNGENID")
dhilmo_age_full <- left_join(dhilmo_age_full, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

dhilmo_age <- dhilmo_age_full

age_low <- seq(0, 80, 10)

# convert temporarily to matrix for fast assignment (at least 1000x faster)
dhilmo_agem <- as.matrix(dhilmo_age)

# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:nrow(dhilmo_age)) {  
    age_start <- 1972 - dhilmo_age$YEAR_OF_BIRTH[i]  # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2017, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2017
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]

    # set NA costs to 0 if follow-up started before upper bound of age interval
    # and follow-up ended after the lower bound of age interval
    indices <- which(age_start < age_low + 10 & age_end > age_low & is.na(dhilmo_agem[i,2:10])) + 1
    dhilmo_agem[i,indices] <- 0
    
    # set costs to NA for age interval if follow-up started after upper bound of age interval
    indices <- which(age_start > age_low + 10) + 1
    dhilmo_agem[i,indices] <- NA
}

# convert back to data frame
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)
#fwrite(dhilmo_age, file="hilmo_costs/finngen_R5_hilmo_inpatient_costs_by_age.txt.gz")


dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          year >= 1998) %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c("FINNGENID", sort(names(dhilmo_age[2:10])))]

# make data frame with all FINNGENID's, even those without any costs in the data
dhilmo_age_full <- data.frame(FINNGENID = unique(fg_df$FINNGENID))
dhilmo_age_full <- left_join(dhilmo_age_full, dhilmo_age, by="FINNGENID")
dhilmo_age_full <- left_join(dhilmo_age_full, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

dhilmo_age <- dhilmo_age_full

age_low <- seq(0, 80, 10)

# convert temporarily to matrix for fast assignment (at least 1000x faster)
dhilmo_agem <- as.matrix(dhilmo_age)

# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:nrow(dhilmo_age)) {  
    age_start <- 1972 - dhilmo_age$YEAR_OF_BIRTH[i]  # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2017, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2017
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]

    # set NA costs to 0 if follow-up started before upper bound of age interval
    # and follow-up ended after the lower bound of age interval
    indices <- which(age_start < age_low + 10 & age_end > age_low & is.na(dhilmo_agem[i,2:10])) + 1
    dhilmo_agem[i,indices] <- 0
    
    # set costs to NA for age interval if follow-up started after upper bound of age interval
    indices <- which(age_start > age_low + 10) + 1
    dhilmo_agem[i,indices] <- NA
}
# convert back to data frame
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)
head(dhilmo_age)
#fwrite(dhilmo_age, file="hilmo_costs/finngen_R5_hilmo_costs_by_age_post1998.txt.gz")

d1 <- fread()
d2
d3

dhilmo_age$eof <- ifelse(is.na(dhilmo_age$DEATH_YEAR), 2017, dhilmo_age$DEATH_YEAR)
dhilmo_age$age_eof <- as.numeric(dhilmo_age$eof) - as.numeric(dhilmo_age$YEAR_OF_BIRTH)
dhilmo_age$age_sof <- 1972 - as.numeric(dhilmo_age$YEAR_OF_BIRTH)
head(dhilmo_age)
tail(dhilmo_age)


dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               post1998 = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          year>1998) %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c(1,10,9,8,3,4,2,5,6,7)]
#names(dhilmo_age)[2:10] <- paste(names(dhilmo_age)[2:10], "_HILMO_POST1998", sep="")
dhilmo_age_post1998 <- data.frame(FINNGENID = unique(endbig$FINNGENID))
dhilmo_age_post1998 <- left_join(dhilmo_age_post1998, dhilmo_age, by="FINNGENID")
head(dhilmo_age_post1998)

dhilmo_age <- left_join(dhilmo_age, fg_df[,c("FINNGENID", "YEAR_OF_BIRTH", "DEATH_YEAR")])

# define helper vector that contains lower bounds of the age intervals
age_low <- seq(0,80,10)
tic()
# loop to make costs that are NA to 0 if patient was followed-up during that age interval
for (i in 1:100) {
    id <- dhilmo_age$FINNGENID[i]
    
    # calculate age at start of follow-up
    yob <- fg_df$YEAR_OF_BIRTH[fg_df$FINNGENID == id]
    age_start <- yob - 1972 # can be negative
    
    # calculate age at end of follow-up/death
    yod <- fg_df$DEATH_YEAR[fg_df$FINNGENID == id]
    yod <- ifelse(is.na(yod), 2017, yod) # if no death recorded, assume complete follow-up until 2017
    age_end <- yod - yob
    
    # loop over the age categories for the individual
    # setting NA costs to 0 if follow-up started before lower bound of age interval
    # and follow-up ended after the lower bound of age interval
    for (j in 2:10) {
        if (age_start < age_low[j-1] &
            age_end > age_low[j-1] &
            is.na(dhilmo_age[i,j])) {
            dhilmo_age[i,j] <- 0
        }
    }
}
toc()

age_low <- seq(0,80,10)
library(tictoc)
tic()
# loop to make costs that are NA to 0 if patienincit was followed-up during that age interval
for (i in 1:100) {  
    
    # calculate age at start of follow-up
    age_start <- dhilmo_age$YEAR_OF_BIRTH[i] - 1972 # can be negative
    
    # calculate age at end of follow-up/death
    yod <- ifelse(is.na(dhilmo_age$DEATH_YEAR[i]), 2017, dhilmo_age$DEATH_YEAR[i]) # if no death recorded, assume complete follow-up until 2017
    age_end <- yod - dhilmo_age$YEAR_OF_BIRTH[i]
    
    # loop over the age categories for the individual
    # setting NA costs to 0 if follow-up started before lower bound of age interval
    # and follow-up ended after the lower bound of age interval
    checkindices <- which(age_start < age_low & age_end > age_low & is.na(dhilmo_age[i,2:10])) + 1
    
    dhilmo_age[i,checkindices] <- 0
    
}
toc()



head(dhilmo_agem)
dhilmo_age <- as.data.frame(dhilmo_agem, stringsAsFactors=FALSE)

str(dhilmo_age)

head(dhilmo_age)

((nrow(dhilmo_age)/100)*0.026)/(60*60)

names(fg_df)

dhilmo_age_post1998[rowSums(dhilmo_age_post1998[,2:10], na.rm=TRUE)==0,]

dhilmo_age <- dhilmo %>% mutate(age_group = case_when(age_at_event<10 ~ "0_10",
                                                 age_at_event<20 ~ "10_20",
                                                 age_at_event<30 ~ "20_30",
                                                 age_at_event<40 ~ "30_40",
                                                 age_at_event<50 ~ "40_50",
                                                 age_at_event<60 ~ "50_60",
                                                 age_at_event<70 ~ "60_70",
                                                 age_at_event<80 ~ "70_80",
                                                 age_at_event<90 ~ "80_90"),
                               inpat = case_when(year>1998 ~ 1,
                                                   year<=1998 ~ 0)) %>%
    filter(!is.na(age_group),
          source=="INPAT") %>% 
    group_by(FINNGENID, age_group) %>%
    summarize(HILMO_COSTS=sum(unit_cost))

dhilmo_age <- dhilmo_age %>%
    pivot_wider(names_from="age_group", values_from="HILMO_COSTS")

dhilmo_age <- dhilmo_age[,c(1,10,9,2,3,4,5,6,7,8)]
#names(dhilmo_age)[2:10] <- paste(names(dhilmo_age)[2:10], "_HILMO_inpat", sep="")
dhilmo_age_inpat <- data.frame(FINNGENID = unique(endbig$FINNGENID))
dhilmo_age_inpat <- left_join(dhilmo_age_inpat, dhilmo_age, by="FINNGENID")
head(dhilmo_age_inpat)

head(dhilmo_age_HILMO %>% summarize(across(contains("_"), ~sum(.x, na.rm = TRUE))))



summary(dhilmo$age_group)

fwrite(indiv_costs, file="hilmo_costs/finngen_R5_hilmo_costs.txt.gz")
fwrite(indiv_costs_post1998, file="hilmo_costs/finngen_R5_hilmo_costs_post1998.txt.gz")

dat1 <- fread("zcat hilmo_costs/data/finngen_R5_hilmo_costs_post1998.txt.gz", data.table=FALSE)
dat <- fread("zcat hilmo_costs/data/finngen_R5_hilmo_costs.txt.gz", data.table=FALSE)

head(dhilmo)

summary(dat)
summary(dat1)

library(lattice)
densityplot(log(dat$HILMO_COSTS))

library(lattice)
densityplot(log(dat1$HILMO_COSTS))

dhilmo %>% filter(unit_cost_w_fixed>10^5 & ea != "70" & year>=1998) %>% arrange(desc(unit_cost_w_fixed))

table(dhilmo[dhilmo$unit_cost_w_fixed>3136940 & !(dhilmo$ea %in% c("70", "74", "75", "98")),]$ea)

sum(dhilmo$unit_cost_w_fixed)/sum(dhilmo$unit_cost)

unit_costs_hilmo

nrow(dhilmo[is.na(dhilmo$unit_cost),])

dhilmo$unit_cost <- dhilmo$unit_cost_w_fixed

temp <- dhilmo %>%
    filter(source=="OUTPAT") %>%
    group_by(finngenid, source) %>%
    summarize(n=n()) %>%
    arrange(desc(n))

ids_highoutpat <- temp$finngenid[1:50]

temp <- dhilmo %>%
    filter(source=="INPAT") %>%
    group_by(finngenid, source) %>%
    summarize(n=n()) %>%
    arrange(desc(n))

ids_highinpat <- temp$finngenid[1:50]

nrow(dhilmo[dhilmo$finngenid %in% ids_highoutpat,])


table(dhilmo[dhilmo$finngenid %in% ids_highoutpat,"pdgo"])[order(table(dhilmo[dhilmo$finngenid %in% ids_highoutpat,"pdgo"]))]




table(dhilmo[dhilmo$finngenid %in% ids_highinpat,"pdgo"])[order(table(dhilmo[dhilmo$finngenid %in% ids_highinpat,"pdgo"]))]




dhilmo[dhilmo$hospdays>90 & !(dhilmo$ea %in% c("70", "74", "75", "98")),]

sum(dhilmo$unit_cost[dhilmo$source=="INPAT"],na.rm=T)

dhilmo[which.max(dhilmo$hospdays),]

options(repr.plot.width=10, repr.plot.height=5)
dhilmo %>%
    mutate(year = as.integer(floor(year)), missing_cost=is.na(unit_cost)) %>%
    group_by(missing_cost, year) %>%
    summarize(n=n()) %>%
    ggplot(aes(x=year, y=n, colour=missing_cost)) +
    geom_line() +
    geom_point()

options(repr.plot.width=15, repr.plot.height=5)
dhilmo %>%
    mutate(year = as.integer(floor(year)), missing_cost=is.na(unit_cost)) %>%
    filter(year<2019) %>%
    group_by(year, pala_class, missing_cost) %>%
    summarize(n=n()) %>%
    ggplot(aes(x=year, y=n)) +
    facet_wrap(~pala_class, scales="free") +
    geom_line() +
    geom_point() +
    #scale_x_continuous(breaks=seq(1972,2018,2)) +
    geom_vline(xintercept=c(1998,2018)) +
    theme_bw()

year_costs <- dhilmo %>%
    filter(year<2019, year >= 1998) %>%
    mutate(year_int = floor(year)) %>%
    group_by(year_int, pala_class) %>%
    summarize(sum_cost=sum(unit_cost)) %>%
    arrange(pala_class, year_int) %>%
    group_by(pala_class) %>%
    mutate(cumsum=cumsum(sum_cost))

options(repr.plot.width=9, repr.plot.height=7)
year_costs %>%
    ggplot(aes(x=year_int, y=cumsum, fill=pala_class)) +
    geom_area() +
    theme_bw() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0), breaks=seq(1998,2018,2))

options(repr.plot.width=9, repr.plot.height=7)
year_costs %>%
    ggplot(aes(x=year_int, y=sum_cost, colour=pala_class)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0), breaks=seq(1998,2018,2))

ea_costs <- dhilmo %>%
    filter(year<2019, year >= 1998) %>%
    mutate(year_int = floor(year)) %>%
    group_by(year_int, pala_class) %>%
    summarize(sum_cost=sum(unit_cost)) %>%
    arrange(pala_class, year_int) %>%
    group_by(pala_class) %>%
    mutate(cumsum=cumsum(sum_cost))


ea_costs <- dhilmo %>%
    filter(year<2019, year >= 1998) %>%
    mutate(year_int = floor(year)) %>%
    group_by(ea) %>%
    summarize(sum_cost=round(sum(unit_cost)/10^6,1)) %>%
    arrange(desc(sum_cost))

left_join(ea_costs, unit_costs_hilmo[unit_costs_hilmo$pala_class=="ward",c("label_ea", "ea")])

indiv_costs <- dhilmo %>%
    filter(year >= 1998) %>%
    group_by(finngenid) %>%
    summarize(sum_cost=sum(unit_cost)) %>%
    arrange(desc(sum_cost))

head(indiv_costs)

summary(indiv_costs$sum_cost)

library(lattice)
densityplot(log(indiv_costs$sum_cost, base=10))

densityplot(log(sim, base=10))

qqnorm(log(indiv_costs$sum_cost, base=10))
qqline(log(indiv_costs$sum_cost, base=10))

mean(log(indiv_costs$sum_cost))
sd(log(indiv_costs$sum_cost))

n_individuals <- 217672
costs <- exp(rnorm(mean=9.555, sd=1.385, n=n_individuals))

nrow(indiv_costs)

boxplot(log(indiv_costs$sum_cost, base=10))

boxplot(indiv_costs$sum_cost)

write.csv(indiv_costs, "costs_by_finngenid_30062020.csv", row.names=FALSE)

unit_costs_avo <- read.csv("unit_costs/AvoHILMO unit costs from pdf not tidy.csv", sep=",", as.is=TRUE)
unit_costs_avo <- unit_costs_avo %>% pivot_longer(names_to="yhteystapa", values_to="unit_cost", c("R10", "R50", "R55", "R60", "R20", "R70"))
head(unit_costs_avo)
names(unit_costs_avo) <- tolower(names(unit_costs_avo))
unique(unit_costs_avo$toimintayksikkö)
unique(unit_costs_avo$ammatti)

table(dprim$CATEGORY)

names(dprim)

unique(unit_costs_avo$ammatti)

# modify specialty
dprim$year <- decimal_date(as.Date(dprim$APPROX_EVENT_DAY))
names(dprim) <- tolower(names(dprim))
names(dprim)[c(9:11)] <- c("yhteystapa", "palvelumuoto", "ammattiluokitus")

head(dprim)

# ammattiryhmäs from http://www.stat.fi/meta/luokitukset/ammatti/001-2001/222.html
dprim <- dprim %>%
    mutate(ammatti=case_when(substr(ammattiluokitus, 0, 4) == "2221" ~ "Lääkäri",
                            substr(ammattiluokitus, 0, 3) == "323" ~ "Sairaanhoitaja",
                            TRUE ~ "NA"),
          toimintayksikkö=case_when(palvelumuoto %in% c("T11", "T21", "T22") ~"Vastaanotto",
                                   TRUE ~ "NA"))

table(dprim$ammatti)
table(dprim$toimintayksikkö)

names(dprim)

dprim <- left_join(dprim, unit_costs_avo, by=c("toimintayksikkö", "ammatti", "yhteystapa"))

sum(dprim$unit_cost, na.rm=T)





options(repr.plot.width=6, repr.plot.height=5)
dprim %>%
    mutate(year_int = as.integer(floor(year))) %>%
    filter(year_int<=2018) %>%
    group_by(year_int) %>%
    summarize(n=n()) %>%
    ggplot(aes(x=year_int, y=n)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks=seq(2010,2018)) +
    labs(x="Year",y = "N of episodes") +
    theme_bw()

names(d)[1] <- "finngenid"

d <- left_join(d, dhilmo[,c("rowindex","unit_cost")], by="rowindex")
d <- left_join(d, dprim[,c("rowindex","unit_cost")], by="rowindex")

d$unit_cost <- ifelse(!is.na(d$unit_cost.x), d$unit_cost.x,
                     ifelse(!is.na(d$unit_cost.y), d$unit_cost.y, NA))

summary(d$unit_cost)

names(d)[1] <- "FINNGENID"

names(d)

d <- d[,c(dnames, "unit_cost")]

summary(d$unit_cost)

3136946/17003150 

head(temp)

temp <- d %>%
    filter(!is.na(unit_cost)) %>%
    mutate(year=as.integer(floor(decimal_date(as.Date(APPROX_EVENT_DAY))))) %>%
    group_by(FINNGENID, SOURCE, year) %>%
    summarize(costs=sum(unit_cost))
write.csv(temp, "costs_by_id_source_year_23_4_2020.csv", row.names=FALSE)

summary(temp$costs)

unique(temp$SOURCE)

rm(d)

#load small and big endpoint file
end <- fread("zcat ../R5_COV_PHENO_V1.txt.gz")
end <- as.data.frame(end)

keepcols_end <- c("FINNGENID", "BL_AGE", "BL_YEAR", "SEX", "AGE_AT_DEATH_OR_NOW", "regionofbirth", "cohort")
keepcols_endbig <- c("FINNGENID","FU_END_AGE", 'DEATH','DEATH_AGE','DEATH_YEAR')
end <- end[,keepcols_end]

endbig <- fread("../finngen_R5_V2_endpoint.txt")
endbig <- as.data.frame(endbig)

endbig <- endbig[,keepcols_endbig]

#keepcols_endbig <- keepcols_endbig[(keepcols_endbig %in% names(endbig))]
fg_df <- endbig
fg_df <- left_join(fg_df, end, by="FINNGENID")
fg_df <- fg_df[!is.na(fg_df$AGE_AT_DEATH_OR_NOW),]

fg_df$YEAR_OF_BIRTH <- fg_df$BL_YEAR - fg_df$BL_AGE
fg_df$FU_START_AGE <- 1972 - fg_df$YEAR_OF_BIRTH
fg_df$FEMALE <- ifelse(fg_df$SEX=="female", 1, 0)

save.image(file = "costs_workspace.RData")

workspace.size <- function() {
  ws <- sum(sapply(ls(envir=globalenv()), function(x)object.size(get(x))))
  class(ws) <- "object_size"
  ws
}

workspace.size()



setwd("/home/jsjukara/")
load("./costs/costs_workspace.RData")
library(data.table)
library(tidyverse)
library(lubridate)
#load phenotype data
d <- fread("chr19_44908684_T_C.tsv")
d <- t(d)
rsid <- d[1]
vec <- d[4:length(d)]
vec <- as.integer(vec)
df <- data.frame(FINNGENID = rownames(d)[4:length(d)], rsid = vec)
names(df)[2] <- rsid
rsids <- names(df[2:ncol(df)])

total_costs <- temp %>%
    group_by(FINNGENID) %>%
    summarize(costs=sum(costs))

total_costs <- left_join(total_costs, df, by="FINNGENID")

fit <- lm(log(costs, base=10) ~ chr19_44908684_T_C)
summary(fit)

unit_costs_hilmo


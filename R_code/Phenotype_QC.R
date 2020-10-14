###
# Quality control on phenotype data, analysis groups definition, removing of outliers
###


#load data and perform QC
library(data.table)
local_path = ""


gen_and_phen = readRDS(paste(local_path, ".rds", sep = "")) # table of phenotypes and polygenic score
#rename columns with short names
col_descriprion = fread('.csv') # UKB phenotype columns description
setnames(gen_and_phen, old = col_descriprion[UDI %in% colnames(gen_and_phen)]$UDI,
         new=col_descriprion[UDI %in% colnames(gen_and_phen)]$short_description)
gen_and_phen[,V38:=NULL]

# add binary column about Gbatch 
gen_and_phen[, Centre:=as.factor(Centre)]
gen_and_phen[, Sex:=as.factor(Sex)]
gen_and_phen[, Gbatch:=as.factor(Gbatch)]

# year of visit
gen_and_phen$year_visit = format(as.Date(gen_and_phen$Date_of_measurement),'%Y')

# - select only independent individuals with known height, sex, YOB
gen_and_phen = gen_and_phen[! is.na(Height)][! is.na(Sex)][!is.na(Year_of_birth)][!is.na(Income)]
gen_and_phen$log_height = log10(gen_and_phen$Height)


ethnic_discription = fread('.tsv', select = c(1,2)) # description of 'ethnicity' field
ethnic_stats = gen_and_phen[,.N,by=`Ethnic background`]
ethnic_stats = merge(ethnic_stats, ethnic_discription, 
                     by.x = "Ethnic background", by.y="coding")

place_of_birth_nonUK = fread('.tsv', select = c(1,2)) # description of 'Country of birth non UK' field
non_UK_stats = gen_and_phen[,.N,by=Country_of_birth_non_UK]
non_UK_stats = merge(non_UK_stats, place_of_birth_nonUK, 
                     by.x = "Country_of_birth_non_UK", by.y="coding")

place_of_birth_UK = fread('.tsv') # description of 'Country of birth non UK' field
UK_stats = gen_and_phen[,.N,by=Country_of_birth_UK]
UK_stats = merge(UK_stats, place_of_birth_UK, 
                     by.x = "Country_of_birth_UK", by.y="coding")

cross_table = gen_and_phen[,.N, by=c("Country_of_birth_non_UK", "Country_of_birth_UK", 
                                     "Ethnic background", "Genetic ethnic grouping")
                           ][N>0][order(N, decreasing = T)]
cross_table$Country_of_birth_non_UK = place_of_birth_nonUK[match(
  cross_table$Country_of_birth_non_UK, coding)]$meaning
cross_table$Country_of_birth_UK = place_of_birth_UK[match(
  cross_table$Country_of_birth_UK, coding)]$meaning
cross_table$`Ethnic background` = ethnic_discription[match(
  cross_table$`Ethnic background`, coding)]$meaning

gen_and_phen$country_of_birth = place_of_birth_UK[match(
  gen_and_phen$Country_of_birth_UK, coding)]$meaning
gen_and_phen[country_of_birth=="Elsewhere"]$country_of_birth = place_of_birth_nonUK[match(
  gen_and_phen[country_of_birth=="Elsewhere"]$Country_of_birth_non_UK, coding)]$meaning

gen_and_phen[, country_of_birth:=as.factor(country_of_birth)]


cross_table2 = gen_and_phen[,.N, by=c("country_of_birth",
                                     "Ethnic background", "Genetic ethnic grouping")
                           ][N>400][order(N, decreasing = T)]
cross_table2$`Ethnic background` = ethnic_discription[match(
  cross_table2$`Ethnic background`, coding)]$meaning

# generate custom analysis groups
gen_and_phen[, analysis_group := as.character(NA)]

gen_and_phen[`Ethnic background`==1001][
  `Genetic ethnic grouping`==1][country_of_birth=="England"]$analysis_group = "English"

gen_and_phen[`Ethnic background`==1001][
  `Genetic ethnic grouping`==1][country_of_birth=="Scotland"]$analysis_group = "Scottish"

gen_and_phen[`Ethnic background`==1001][
  `Genetic ethnic grouping`==1][country_of_birth=="Wales"]$analysis_group = "Welsh"

gen_and_phen[`Ethnic background`==1001][
  is.na(`Genetic ethnic grouping`)]$analysis_group = "Other_British"

gen_and_phen[`Ethnic background`==1002][
  country_of_birth %in% c("England", "Republic of Ireland")]$analysis_group = "Irish"

gen_and_phen[`Ethnic background`==1003]$analysis_group = "Other_white"

gen_and_phen = gen_and_phen[!is.na(analysis_group)]


# QC - remove outliers
number_of_sd  = 4.75
gen_and_phen[, outlier := abs(Height - mean(Height)) > number_of_sd * sd(Height), 
                            by = list(analysis_group, Sex)]
gen_and_phen = gen_and_phen[outlier==FALSE]

sitting_to_standing = 0.75
gen_and_phen = gen_and_phen[! is.na(Sitting_height)]
gen_and_phen[, outlier := Sitting_height / Height > sitting_to_standing, 
             by = list(analysis_group, Sex)]

gen_and_phen = gen_and_phen[outlier==FALSE]



saveRDS(gen_and_phen, ".rds")

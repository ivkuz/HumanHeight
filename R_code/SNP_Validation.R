###
# SNP validation
###

library(data.table)

gen_and_phen = fread(".tsv") # table with genotypes and phenotypes
local_path = ""


gen_and_phen[, Centre:=as.factor(Centre)]
gen_and_phen[, Sex:=as.factor(Sex)]
gen_and_phen[, Gbatch:=as.factor(Gbatch)]
gen_and_phen[, UKBiLEVEAX:=as.factor(UKBiLEVEAX)]
gen_and_phen[, norm_height:=scale(Height), by = analysis_group]


mfi_table = fread(paste(local_path, '.tsv', sep = "")) # table with summary statistics from GIANT for top SNPs
mfi_table = mfi_table[rs_id_O_E %in% colnames(gen_and_phen)]

#normalize genotypes
gen_and_phen = gen_and_phen[analysis_group=="English"]
cols <- mfi_table$rs_id_O_E
gen_and_phen[, (cols) := lapply(.SD, scale), .SDcols=cols]


model = "norm_height ~ Sex + Year_of_birth + Centre + Gbatch + PC1 + PC2 + PC3  + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"

get_snp_stat_from_ukb = function(rs_id_O_E){
  summary_lm = summary(lm(as.formula(paste(model, rs_id_O_E, sep = '+')), 
                          data = gen_and_phen[analysis_group=="English"]))$coefficients
  return(summary_lm[rs_id_O_E,])
}

validation_matrix = t(sapply(mfi_table$rs_id_O_E, get_snp_stat_from_ukb))
validation_table = as.data.table(validation_matrix, keep.rownames = T)
setnames(validation_table, old = "rn", new="rs_id_O_E")

write.table(validation_table, file=paste(local_path, '.tsv', sep = ""),
            sep="\t", quote = F, row.names = F)

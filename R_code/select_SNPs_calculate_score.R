###
# Selection SNPS for PGHS and PGHS calculation
###

#working with genetic scores
current.directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current.directory)
library(data.table)
library(limma)
local_path = ""


mfi_table = fread(paste(local_path, 'imp_and_beta_giant.tsv', sep = ""))
giant_table <- fread(paste(local_path,'ng.3097-S2.csv',sep = ""), select=c(1:10,17))
mfi_table$pvalue_giant = giant_table[match(mfi_table$rs_id, SNP)]$`P-value`
mfi_table[, EAF_diff := abs(EAF_GIANT- EAF_UKBB)]
rm(giant_table)


# draw venn diargram
table_for_venn = data.table(mfi_table$rs_id)
table_for_venn$`p>5*10-8` = mfi_table$pvalue_giant >5*10**(-8)
table_for_venn$`p<5*10-10` = mfi_table$pvalue_giant <5*10**(-10)
table_for_venn$`IQ>0.9` = mfi_table$imputation_quality > 0.9
table_for_venn$`EAF_DIF<0.05` = with(mfi_table, abs(EAF_GIANT- EAF_UKBB)) < 0.05


table_for_venn[,V1:=NULL]
venn_counts = vennCounts(table_for_venn)
vennDiagram(venn_counts)

#select good SNPS
rs_ids_to_calculate_score = mfi_table[imputation_quality > 0.9][pvalue_giant < 5*10**(-10)][EAF_diff < 0.05]$rs_id_O_E


genotype_matrix = readRDS(paste(local_path, "genotype_matrix.rds", sep=""))
betas = mfi_table[match(rs_ids_to_calculate_score, rs_id_O_E)]$Beta

genetic_scores = genotype_matrix[, rs_ids_to_calculate_score] %*% betas
genetic_scores = as.data.table(genetic_scores, keep.rownames = T)
genetic_scores[, rn := as.numeric(rn)]
setnames(genetic_scores, old = "V1", new = "score")
write.table(genetic_scores,  paste(local_path, "genetic_scores", sep = ""), 
            quote = F, row.names = F, sep = '\t')


processed_phenotypes = fread(paste(local_path, 'processed_phenotypes.tsv', sep = ""))
# add genetic scores to the phenotype table
gen_and_phen = merge(processed_phenotypes, genetic_scores, by.x="ID", by.y='rn')

# center and scale risk score
gen_and_phen[, norm_score := (score - mean(score)) / sd(score), by = analysis_group]


# merge genotypes to processed phenotypes
genotype_matrix = genotype_matrix[rownames(genotype_matrix) %in% processed_phenotypes$ID, ]
genotype_matrix = as.data.table(genotype_matrix, keep.rownames = T)
genotype_matrix[, rn := as.numeric(rn)]
gen_and_phen = merge(gen_and_phen, genotype_matrix , by.x="ID", by.y='rn')


col_descriprion = fread('../Tables/UKB_phenotype_columns_description - ML.csv')
setnames(gen_and_phen, old = col_descriprion[UDI %in% colnames(gen_and_phen)]$UDI,
         new=col_descriprion[UDI %in% colnames(gen_and_phen)]$short_description)

columns_for_ML = c(col_descriprion$short_description, rs_ids_to_calculate_score)



write.table(gen_and_phen[,..columns_for_ML], file = "/mnt/ukbb/storage/projects/UKB-nonAdditive-41601/shared/ML_dataset.tsv", 
            sep="\t", quote = F, row.names = F)

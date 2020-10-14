###
# PGHS validation analysis
###

library(data.table)
local_path = "/mnt/ukbb/home/sslavskii/height/"

mfi_table = fread(paste(local_path, 'imp_and_beta_giant.tsv', sep = ""))
giant_table <- fread(paste(local_path,'ng.3097-S2.csv',sep = ""), select=c(1:10,17))
mfi_table$pvalue_giant = giant_table[match(mfi_table$rs_id, SNP)]$`P-value`
mfi_table$SE_giant = giant_table[match(mfi_table$rs_id, SNP)]$SE
mfi_table[, EAF_diff := abs(EAF_GIANT- EAF_UKBB)]
rm(giant_table)


validation_table = fread(paste(local_path, 'validation_table.tsv', sep = ""))
setnames(validation_table, old=c("Estimate", "Pr(>|t|)"), new = c("Beta_UKB", "Pvalue_UKB"))


validation_table = merge(validation_table, mfi_table, by = "rs_id_O_E")
validation_table[,y:= -sign(Beta_UKB)*log10(Pvalue_UKB)]
validation_table[,x:= -sign(Beta_GIANT)*log10(pvalue_giant)]
validation_table[, pass_qc :=  pvalue_giant < 5*10**(-10) & EAF_diff < 0.05 & imputation_quality > 0.9]



validatio_threshold = 0.05/419


pdf("../Figures/GIANT_UKB.pdf")
with(validation_table[abs(x)<12], plot(x, y, 
                                       xlab = "sign(Z)*logP GIANT", ylab = "sign(Z)*logP UKB", 
                                       col = as.factor(Pvalue_UKB< validatio_threshold))) 
abline(v = 0)
abline(v = c(-log10(5*10**-8), log10(5*10**-8)), col="blue")
abline(v = c(-log10(5*10**-10), log10(5*10**-10)), col="darkgreen")
abline(h = c(-log10(5*10**-8), log10(5*10**-8)), col="blue")
abline(h = c(-log10(5*10**-10), log10(5*10**-10)), col="darkgreen")
legend("bottom", bty="n", legend=paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))
legend("bottomright", bty="n", legend=paste("slope=", format(summary(fit)$coefficients["x",1], digits=4)))
dev.off()

number_of_snps = nrow(validation_table[pvalue_giant < 5*10**-8])
nrow(validation_table[Pvalue_UKB < 0.05/number_of_snps ])/419
replication_rate_8 = nrow(validation_table[pvalue_giant < 5*10**-8][Pvalue_UKB < 0.05/number_of_snps ])/number_of_snps
replication_rate_10 = nrow(validation_table[pvalue_giant < 5*10**-10][Pvalue_UKB < 0.05/number_of_snps ])/nrow(validation_table[pvalue_giant < 5*10**-10])
dev.off()





pdf("../Figures/GIANT_validation.pdf")
layout(matrix(c(1,1,2,2,0,3,3,4), ncol = 4, byrow = TRUE))

par(mar = c(5,5,1,0))
with(validation_table, plot(x=Beta_GIANT / SE_giant, y=`t value`, xlab = "Z, GIANT", ylab = "Z, UKB", bty="n", 
                            xlim = c(-25,25), ylim = c(-25,25),
                            col = as.factor(pass_qc), pch = ((Pvalue_UKB> validatio_threshold) +1 ), axes=F  )) 
axis(1, col='gray44')
axis(2, col='gray44')
abline(h=0, col='gray44')
abline(v = 0, col='gray44')
mtext('A', side = 3, at = -19, line = -2)
abline(fit <- lm(`t value` ~ I(Beta_GIANT / SE_giant), data = validation_table), col='gray44')
legend("bottom", bty="n", legend=paste("R2=", format(summary(fit)$adj.r.squared, digits=4)))
legend("bottomright", bty="n", legend=paste("slope=", format(summary(fit)$coefficients["I(Beta_GIANT/SE_giant)",1], digits=4)))

par(mar = c(5,5,1,1))
with(validation_table[abs(log10(pvalue_giant)) < 12], plot(x, y, xlab = "sign(Z)*log10(P), GIANT", ylab = "sign(Z)*log10(P), UKB", 
                                                           col = as.factor(pass_qc), pch = ((Pvalue_UKB> validatio_threshold) +1 ),
                                                           bty="n", xlim = c(-12,12), ylim = c(-21,21), axes=F)) 
axis(1, col='gray44')
axis(2, col='gray44')
abline(h=0, col='gray44')
abline(v = 0, col='gray44')
abline(v = c(-log10(5*10**-8), log10(5*10**-8)), col="blue")
abline(h = c(-log10(5*10**-8), log10(5*10**-8)), col="blue")
mtext('B', side = 3, at = -11, line = -2)


par(mar = c(5,5,1,0))
with(validation_table, plot(x= EAF_GIANT, y=EAF_UKBB, xlab = "EAF, GIANT", ylab = "EAF, UKB", bty="n", xlim = c(0,1), ylim = c(0,1),
                            col = as.factor(pass_qc), pch = ((Pvalue_UKB> validatio_threshold) +1 ), axes = F )) 
axis(1, col='gray44')
axis(2, col='gray44')
abline(a = 0.05,b = 1, col='gray44')
abline(a = -0.05,b = 1, col='gray44')
mtext('C', side = 3, at = 0.05, line = -2)

par(mar = c(5,0,1,0))
plot(1, type = "n", axes=F, xlab="", ylab="")
legend(x=0.6, y = 1.2,
       col = c("red", col="black"), 
       legend = c("Used for PGHS contruction",
                  "Excluded from PGHS"),
       bty='n',
       pch=1,
       title.adj = 0.2)

legend(x=0.6, y = 1.4,
       legend = c("Passed replication",
                  "Did not pass replication"),
       bty='n',
       pch=c(1,2),
       title.adj = 0.2)

dev.off()

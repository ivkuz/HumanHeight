###
# Calculate residuals for height and log-height
###

library(data.table)

rp <- fread("") # Table with PGHS and phenotypes
addreg <- lm(Height ~ Sex + norm_score + predF3, data = rp)
multreg <- lm(log10(Height) ~ Sex + norm_score + predF3, data = rp)
residuals_add <- rp$Height - predict(addreg)
residuals_log <- log10(rp$Height) - predict(multreg)
residuals_mult <- rp$Height - 10^(predict(multreg))
resid <- data.table(add = residuals_add, log = residuals_log, 
                        mult = residuals_mult)
save(resid, file = "residuals.RData")


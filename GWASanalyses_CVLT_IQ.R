# RNA GWAS
# swedish data

# load data and add neurocog variable
setwd("/mnt/nfs/swe_gwas/ABZ/RNA_GWAS")
load("swedenclean.rdata")
CVLT = read.table("GWAS-CVLT.txt",header=T)
IQ = read.table("GWAS-IQ.txt",header=T)
swedenclean$CVLT = CVLT[match(swedenclean$StudyID,CVLT$IID),3]
swedenclean <- swedenclean[order(swedenclean$CVLT,na.last=F),]
swedenclean$IQ = IQ[match(swedenclean$StudyID,IQ$IID),3]
swedenclean[1:30,18564:18567]

# GWAS 
# run mixed regression for each marker
# save summary statistics

library(lme4)
RNA = names(swedenclean)[7:18565]

# linear mixed model; predictor = expression 
# covar = age, sex; random = family ID
models <- lapply(RNA, function(x) {
  lmer(substitute(CVLT ~ i + IQ + Age + Sex + (1|Family), 
                  list(i = as.name(x))), data = swedenclean)
})

# summary statistics
model_stats <- lapply(models, function(x) {
  coef(summary(x))
})

# save results, age, sex separately
results = NULL
age = NULL
sex = NULL
for (i in 1:18559) {
  temp_results <- as.data.frame(model_stats[[i]])
  results <- rbind(results,temp_results[2,])
  age <- rbind(age,temp_results[3,])
  sex <- rbind(sex,temp_results[4,])
}

# compute degrees of freedom and p-values
colnames(results)[2:3] <- c("SE","tvalue")
require(pbkrtest)
results$df.KR <- sapply(models,function(x) {
  get_Lb_ddf(x, fixef(x))
})
results$p.KR <- 2 * (1 - pt(abs(results$tvalue),results$df.KR))

rownames(age) = RNA
rownames(sex) = RNA
colnames(age)[2:3] <- c("SE","tvalue")
colnames(sex)[2:3] <- c("SE","tvalue")
age$p.KR <- 2 * (1 - pt(abs(age$tvalue),results$df.KR))
sex$p.KR <- 2 * (1 - pt(abs(sex$tvalue),results$df.KR))

results <- results[order(results$p.KR),]
age <- age[order(age$p.KR),]
sex <- sex[order(sex$p.KR),]

results[1:25,]

# write results
write.table(results,"RNA-GWAS-CVLT_IQ_results.txt",sep="\t",
            col.names=T,row.names=T,quote=F)
write.table(age,"RNA-GWAS-CVLT_IQ_age-effects.txt",sep="\t",
            col.names=T,row.names=T,quote=F)
write.table(sex,"RNA-GWAS-CVLT_IQ_sex-effects.txt",sep="\t",
            col.names=T,row.names=T,quote=F)


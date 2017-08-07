# RNA GWAS
# swedish data
# exclude outlier (ID = 2008601953)

# load data and add neurocog variable
setwd("/mnt/nfs/swe_gwas/ABZ/RNA_GWAS")
load("swedenclean.rdata")
CVLT = read.table("GWAS-CVLT.txt",header=T)
swedenclean$CVLT = CVLT[match(swedenclean$StudyID,CVLT$IID),3]
swedenclean <- swedenclean[order(swedenclean$CVLT,na.last=F),]
swedenclean <- swedenclean[c(1:191,193),] # exclude outlier
swedenclean[180:192,18564:18566]

results_orig <- read.table("RNA-GWAS-CVLT_nooutlier_results.txt",header=T)
results_orig <- results_orig[order(results_orig$p.KR),]
hits <- rownames(results_orig)[1:76]

hitsdf <- swedenclean[,colnames(swedenclean)%in%hits]
hitsdf <- cbind(swedenclean[,c(1:4,18566)],hitsdf)

# GWAS 
# run mixed regression for each marker
# save summary statistics

library(nlme)
#RNA = names(swedenclean)[7:18565]
#RNAtest = RNA[1:2]
#RNAtest10 = RNA[1:10]

# linear mixed model; predictor = expression 
# covar = age, sex; random = family ID
models <- lapply(hits, function(x) {
  lme(eval(substitute(CVLT ~ i + Age + Sex, list(i = as.name(x)))), random= ~1|Family, data = swedenclean)
})

# summary statistics
model_stats <- lapply(models, function(x) summary(x)$tTable)

# save results
results = NULL
for (i in 1:76) {
  temp_results <- as.data.frame(model_stats[[i]])
  results <- rbind(results,temp_results[2,])
}

# compute degrees of freedom and p-values
colnames(results)[4:5] <- c("tvalue","pvalue")
results <- results[order(results$pvalue),]
results$p.bon <- p.adjust(results$pvalue,method="bonferroni",n=18559)
results$p.FDR <- p.adjust(results$pvalue,method="BH",n=18559)
results$marker <- rownames(results)

# write results
write.table(results,"RNA-GWAS-CVLT_nooutlier_results_nlme.txt",sep="\t",
            col.names=T,row.names=T,quote=F)



# C-C gene expression
# swedish twins
# check for up- and down-regulation
# case-control (N = 117) - MDD included in control

# load data
setwd("~/ABZ/Yale/RNA_SwedishTwins/CVLT_GWAS_CorrectNaming")
load("~/ABZ/Yale/RNA_SwedishTwins/swedenclean.rdata")
dx = read.table("~/ABZ/Yale/RNA_SwedishTwins/swedishNP_PCA.txt",header=T)
swedenclean$dx_new <- dx[match(swedenclean$StudyID,dx$IID),7]
swedenclean$dx_new <- as.factor(swedenclean$dx_new)

# case-control (N = 117)
genecheck <- subset(swedenclean,dx_new==1|dx_new==7|dx_new==5|dx_new==6,
                    select=c(1:2,18566,7:18565))
genecheck$dx_group <- genecheck$dx_new
genecheck$dx_group <- droplevels(genecheck$dx_group)
levels(genecheck$dx_group) = c("SZ","HC","HC","HC")
genecheck$dx_group <- droplevels(genecheck$dx_group)

# z-score based on control mean, SD
zs_sweden <- genecheck[order(genecheck$dx_group),]
zs_sweden <- zs_sweden[,c(1:3,18563,4:18562)]
for (i in 18564:37122) {
  meanhc <- mean(zs_sweden[37:117,i-18559])
  sdhc <- sd(zs_sweden[37:117,i-18559])
  zs_sweden[,i] <- (zs_sweden[,i-18559]-meanhc)/sdhc
  names(zs_sweden)[i] <- paste(names(zs_sweden)[i-18559],"z",sep="_")
}

colnames(zs_sweden)[1:2] = c("IID","FID")

# just patients
zs.swe.pt <- zs_sweden[1:36,]

# melt data
library(reshape)
zs.melt <- melt(zs.swe.pt[,c(1:2,18564:37122)],id=c("IID","FID"))
colnames(zs.melt)[3:4] = c("gene","expression")
zs.table <- aggregate(zs.melt,by="gene",mean)
zs.table <- zs.table[order(zs.table$expression),]
#zs.se <- summarySE(zs.melt,measurevar="expression",groupvars="gene")
#zs.se <- zs.se[order(zs.se$expression)]

# graph z-scores of patients vs. controls
library(ggplot2)
ggplot(zs.table, aes(y=expression,x=gene)) + 
  geom_bar() +
  geom_hline(yintercept=0) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=12),
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black")) +
  ylab("Patient Expression Levels") +
  xlab("18,559 Genes (Genome-Wide)")



################################################################### 
# testing
zs_sweden_test <- genecheck[order(genecheck$dx_group),]
zs_sweden_test <- zs_sweden_test[,c(1:3,18563,4:18562)]
for (i in 18564:18600) {
  meanhc <- mean(zs_sweden_test[37:117,i-18559])
  sdhc <- sd(zs_sweden_test[37:117,i-18559])
  zs_sweden_test[,i] <- (zs_sweden_test[,i-18559]-meanhc)/sdhc
  names(zs_sweden_test)[i] <- paste(names(zs_sweden_test)[i-18559],"z",sep="_")
}

colnames(zs_sweden)[1:2] = c("IID","FID")

library(reshape)
genes <- genecheck[,3:18562]
genes <- melt(genes,id="dx_new")
dxmeans <- cast(genes, variable~dx_new, mean)
colnames(dxmeans) <- c("gene","scz","hc")
dxmeans$cc_diff <- dxmeans$scz - dxmeans$hc
dxmeans <- dxmeans[order(dxmeans$cc_diff),]
# 783 genes under-expressed in SZ; 217 over-expressed

###################################################################
# summarySE function from 
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

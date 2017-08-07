# bubble graphs for posters
# RNA GWAS - CNS 2016

# graph the functional annotation clusters
# cluster x enrichment score 
# size of bubble is max number of genes per term

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/RNA_GWAS/"
setwd(workdir)

libs <- c("dplyr", "ggplot2","stringr")
invisible(lapply(libs, require, character.only = TRUE))

# load in data
func.ann <- read.table("func-ann_RNA.txt", header=T, sep=",")

# graph
ggplot(func.ann, aes(x=cluster, y=enrichscore, label=cluster)) + 
  geom_point(position="jitter",aes(colour=maxgenes, size=maxgenes)) +
  geom_text(size=5,aes(label=function(x) str_wrap(x,5))) +
  scale_colour_gradientn(colours=topo.colors(5),guide="none") + 
  scale_size_continuous(range=c(15,40),name="Genes in \nCluster",
                        breaks=c(3,9,15)) + 
  scale_x_discrete(expand=c(.15,0)) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_text(size=13, color="black"),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  xlab("Cluster") +
  ylab("Enrichment Score") +
  expand_limits(y=c(-.2,3))
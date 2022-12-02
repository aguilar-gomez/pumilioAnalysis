library(ggplot2)
library(tidyverse)
library(ggrepel)

setwd("~/Documents/Opumilio/Ohana/selscanv4")
annotation_genome <- read.delim("~/Documents/Opumilio/Rescaffolded/Rescaffolded_annoTranscriptome.bed", header=FALSE)
annotation_resca<-annotation_genome[c(1,2,3,12)]
colnames(annotation_resca)<-c("scaffold","start","end","genes")
positions <- read.delim("pumilio.beagle.pos", header=FALSE)
position<-paste0(positions$V1,":",positions$V2)


k=0

filename<-paste0("scan.pumiliok",k,".txt")
scan.pumilio <- read.delim(filename)
rownames(scan.pumilio)<-position
scan.pumilio$scaffold<-positions$V1
scan.pumilio$pos<-positions$V2

small<-scan.pumilio[scan.pumilio$lle.ratio>0,]

#K1 SC #K2 HP #K3 CM #K4 TB #K5 CL #K6 AG #K7 PP
outquantile<-.9995
#cuts<-c(20,30,35,40,20,65,55)
#Find outliers
cutoff<-quantile(small$lle.ratio,outquantile)
#cutoff<-cuts[k+1]
hits<-scan.pumilio[scan.pumilio$lle.ratio>cutoff,]
anno<-merge(hits,annotation_resca)

annotated<-anno[(anno$pos>=anno$start) & (anno$pos<=anno$end),]
annotated<-separate(annotated, genes, into = c("first", "rest"), sep = "\\s",
               extra = "merge")

#There is anotehr outlier SNP in the scaffold
hitsPerWindow<-annotated%>% group_by(gene=tolower(first))%>%filter(n() > 1)

anno_reduced<-hitsPerWindow %>% group_by(gene=tolower(first)) %>% top_n(1, lle.ratio)
anno_reduced<-anno_reduced %>% group_by(gene) %>% top_n(1, start)
anno_reduced<-anno_reduced %>% group_by(gene) %>% top_n(1, end)
anno_reduced<-anno_reduced %>% group_by(gene) %>% top_n(1, pos)
anno_reduced<-unique(anno_reduced[c("scaffold","pos","lle.ratio","gene")])

#colors=c('#008000','#FF0000','#7FCC12','#F0E442', '#0000FF','#FFA500','#8B2121')
colors<-c("#7FCC12", '#8B2121','#008000','#F0E442', 
'#0000FF','#FF0000','#FFA500')




half_data<-scan.pumilio[scan.pumilio$lle.ratio>quantile(scan.pumilio$lle.ratio,.5),]
half_data$SNP<-1:nrow(half_data)

half_anno<-merge(half_data,anno_reduced,all.x = T)
ggplot(data=half_anno,aes(x=SNP,y=lle.ratio,label=gene))+geom_point(color=colors[k+1])+
  ggtitle(paste0("selection on k=",k+1))+theme_classic(base_size=20)+
  geom_label_repel(size=5,max.overlaps = 60,fontface = 'bold') +ylim(c(0,NA))
ggsave(filename=paste0("k",k+1,"selscan_anno_2SNPS.png"),width=10,height=4)



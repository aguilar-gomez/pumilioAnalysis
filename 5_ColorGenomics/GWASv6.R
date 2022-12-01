setwd("~/Documents/Opumilio/gemma/version6/")
library(ggplot2)
library(tidyverse)
library(ggrepel)

simpleQQPlot<-function (observedPValues,min=0,max=1,color) {
  x=-log10(sort(runif(n=length(observedPValues),min,max)))
  y=-log10(sort(observedPValues))
  all1 = c(x,y)
  range = c(min(all1), max(all1))
  plot(x,y,
       xlab = "-log10(expectedPValues)",
       ylab = "-log10(observedPValues)",
       xlim=range, ylim=range,col=color
  )
  abline(0, 1, col = "red")
}


#Reference annotation  
colors<-c("black","gray","black",
          "purple","blue","#008000",
          "gold","red","black","magenta",
          "hotpink","orange","turquoise","magenta","magenta")
names<-c("class1","brightness","S1_Ultraviolet",
         "S1V","blue","green","yellow",
         "red","blackprop","vgg16_k6",
         "pc1_sizeSpots","pc2_spotedness","sex","vgg16_k3","vgg16_k4")
annotation_genome <- read.delim("~/Documents/Opumilio/Rescaffolded/Rescaffolded_annoTranscriptome.bed", header=FALSE)
annotation_resca<-annotation_genome[c(1,2,3,12)]
colnames(annotation_resca)<-c("chr","start","end","genes")

#5,6,8
k<-9
pum.cov10.assoc <- read.delim(paste0("n",k,"__rescaffold.assoc.txt"))
#pum.cov10.assoc <- read.delim("./sex_nopcs_rescaffold.assoc.txt")
#pum.cov10.assoc <- read.delim("./k4_rescaffold.assoc.txt")
#pum.cov10.assoc<-pum.cov10.assoc[is.na(pum.cov10.assoc$p_wald)==F,]
color<-colors[k]

#Unannotated plot

valid<-na.omit(pum.cov10.assoc)
small<-pum.cov10.assoc[pum.cov10.assoc$p_wald<0.1,]
cutoff<-quantile(small$p_wald,na.rm = T,.0001)
hits<-pum.cov10.assoc[pum.cov10.assoc$p_wald<cutoff,]

#Annotate

anno<-merge(hits,annotation_resca)

annotated<-anno[(anno$ps>=anno$start) & (anno$ps<=anno$end),]
annotated<-separate(annotated, genes, into = c("first", "rest"), sep = "\\s",
                    extra = "merge")
annotated$first[annotated$first=="BRMS1L"]="BRMS1"
annotated$first[annotated$first=="PPP4R1L"]="PPP4R1"
#There is anotehr outlier SNP in the scaffold
hitsPerWindow<-annotated%>% group_by(gene=tolower(first))%>%filter(n() > 1)

anno_reduced<-hitsPerWindow %>% group_by(gene=tolower(first)) %>% top_n(1, -log(p_wald))
anno_reduced<-anno_reduced %>% group_by(gene) %>% top_n(1, start)
anno_reduced<-anno_reduced %>% group_by(gene) %>% top_n(1, end)
anno_reduced<-anno_reduced %>% group_by(gene) %>% top_n(1, ps)
anno_reduced<-unique(anno_reduced[c("chr","ps","af","beta","p_wald","gene")])

anno_small<-merge(small,anno_reduced,all.x = T)

#Sort
chr.size<-data.frame(sort(table(anno_small$chr),decreasing = TRUE))
colnames(chr.size)<-c("chr","size")
all<-merge(chr.size,anno_small,by="chr")
sortedscan<-all[order(all[2],decreasing = TRUE),]
sortedscan$SNPs<-c(1:length(sortedscan$chr))
n_of_chr<-length(unique(small$chr))

sortedscan$gene[sortedscan$gene==""]=NA
sortedscan$colorme[is.na(sortedscan$gene)]<-"0"

chrcolor<-sortedscan$chr[!is.na(sortedscan$gene)]
#sortedscan$colorme[!is.na(sortedscan$gene1)]
sortedscan$colorme[sortedscan$chr%in%chrcolor]<-"1"



h<-ggplot(sortedscan, aes(x=SNPs,y=-log10(p_wald),color=colorme,label=gene,alpha=colorme))+
  geom_point(size=1)+theme_classic(base_size=20)+xlab("Scaffolds")+ 
  scale_color_manual(values=rep(c("gray48",color),n_of_chr))

h+theme(plot.title = element_text(face="italic"),legend.position = "none")+
  geom_label_repel(size=5,max.overlaps = 60,fontface = 'bold')

ggsave(paste0("version6",names[k],".png"), width = 10, height = 5)


write.table(anno_reduced,
            paste0(names[k],".genename"),
            quote = F,sep = "\t",row.names = F, col.names = F)



png(paste0("qqplot_pval",names[k],".png"), width=600, height=600)
simpleQQPlot(pum.cov10.assoc$p_wald,min=0,max=1,color)
dev.off()


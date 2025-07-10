setwd("~/Documents/Berkeley/Pumilio/PCAv11")
library(ggplot2)
library(ggrepel)
library(RcppCNPy)
library("grid")


ALL<- as.matrix(read.table("pca_pumilio_v11_newversion.cov")) 
newlabels<- read.delim("sorted_labels_RA", header=FALSE)
colnames(newlabels)<-c("Individual","Random","Location")
pumilio <- read.delim("pumilio.fam", header=FALSE)
colnames(pumilio)<-c("code","Population","bamname","h","o","l","a")
#Labels rescaffold
pumilio.res <- read.table("pumilio.r.bamlist", quote="\"", comment.char="")
pumilio.res$Individual<-gsub("_S.*bam","",pumilio.res$V1)

pumilio$Individual<-gsub("_S.*bam","",pumilio$bamname)
pumilio$Individual[pumilio$Individual=="TB03"]<-"PP37"
reslabels<-merge(pumilio.res,pumilio)



# Eigenvalues
eig <- eigen(ALL, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

 # Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Population <- factor(newlabels$Location)
reslabels$Population[reslabels$Population=="PP"]<-"PO"
reslabels$Population[reslabels$Population=="CMR"]<-"BMR"
reslabels$Population[reslabels$Population=="CMY"]<-"BMY"
reslabels$Population[reslabels$Population=="AG"]<-"SK"

reslabels$Population[reslabels$Individual=="PP34"]<-"PD"
reslabels$Population[reslabels$Individual=="PP35"]<-"PD"
reslabels$Population[reslabels$Individual=="PP36"]<-"PD"


PC$Population <- reslabels$Population
PC$Population <- factor(PC$Population,
                        levels=c("CL", "BMR", "BMY",
                                 "HP","TB","PO","PD",
                                 "DBB","DBI","DBR",
                                 "RA","SK","SC"))
comp<-c(1,2)
comp<-c(3,4)
title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")
xlabel = paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)",sep="")
ylabel = paste("PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="")
x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")
rownames(PC)<-as.character(reslabels$Individual)
PC$ind<-reslabels$Individual
cbPalette <- c("#7FCC12", "#efe133", "#A38C0E" ,
               "orange","red","#008000","#00C000",
               "blue", "purple","magenta4",
               "dodgerblue4","black","brown")
set.seed(42)  # Set seed for reproducibility
PC <- PC[sample(nrow(PC)), ]

# Create the plot
g <- ggplot() +
  geom_point(data = PC, aes_string(x = x_axis, y = y_axis, color = "Population"), alpha = 0.8, size = 4) +
  scale_colour_manual(values = cbPalette) +
  theme_bw(base_size = 25) +
  xlab(xlabel) +
  ylab(ylabel) +
  theme(  text = element_text(family = "Avenir"),
    panel.border = element_rect(colour = "black", fill = NA, size = 3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

g
ggsave(filename=paste0(x_axis,y_axis,"_PCAngsd_pumilio_allpops_2025.png"),width=7.5,height=9)









# Shuffle the rows of the data for random point order
set.seed(42)  # Optional: For reproducibility
PC <- PC[sample(nrow(PC)), ]

# Create the Cartesian coordinate plot
g <- ggplot() +
  geom_point(data = PC, aes_string(x = x_axis, y = y_axis, color = "Population"), alpha = 0.8, size = 3) +
  geom_hline(yintercept = 0, color = "black", size = 1) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, color = "black", size = 1) +  # Vertical line at x=0
  scale_colour_manual(values = cbPalette) +
  theme_minimal(base_size = 25) +  # Use a minimal theme
  xlab(xlabel) +
  ylab(ylabel) +
  theme(
    panel.border = element_blank(),  # Remove panel border
    panel.grid = element_blank(),    # Remove grid lines
    axis.line = element_blank(),     # Remove default axis lines
    axis.text = element_text(size = 15),  # Adjust axis text size
    legend.position = "top",         # Place legend at the top
    legend.title = element_blank()   # Remove legend title
  )

g

ggsave(filename=paste0(x_axis,y_axis,"_PCAngsd_pumilio_allpops_v11.png"),width=7.2,height=8)










#print(g,angle = 45)
ggsave(filename=paste0(x_axis,y_axis,"_PCAngsd_pumilio_allpops_v8.png"),width=7.2,height=8)


g<-ggplot(data=PC, aes_string(x=x_axis, y=y_axis, color="Population",label="ind")) + geom_point(alpha = .8,size=3)+
  scale_colour_manual(values=cbPalette)+geom_text()+
  theme_bw(base_size = 25) +xlab(xlabel)+ylab(ylabel)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top",legend.title=element_blank())+
  geom_label_repel(size=5,max.overlaps = 60,fontface = 'bold')

g
ggsave(filename=paste0(x_axis,y_axis,"_PCAngsd_pumilio_allpops_rescafffold_labels.png"),width=7.2,height=8)



plot(1:15,cumsum(eig$val)[1:15],xlab = "Number of Components",ylab="Variance (%)",main="PCA Phrynocephalus")
plot(1:30,cumsum(eig$val)[1:30],xlab = "Number of Components",ylab="Variance (%)",main="PCA Phrynocephalus")
PC$Population <- factor(newlabels$Sex)
  comp<-c(3,4)
  title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")
  
  x_axis = paste("PC",comp[1],sep="")
  y_axis = paste("PC",comp[2],sep="")
  rownames(PC)<-labels$Individual
  
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
  ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Population",shape="Population"),size=4,alpha = .8)+
    ggtitle(title) + scale_colour_manual(values=cbPalette)+ theme_bw(base_size = 20) 
  
  
#Add intercept
PC$cov1<-1
gemma_rescaffold_cov<-PC[c("cov1",paste0("PC",1:10))]
gemma_rescaffold_cov<-PC[c("cov1",paste0("PC",1:20))]
gemma_rescaffold_cov<-PC[c("cov1",paste0("PC",1:30))]
gemma_rescaffold_cov<-PC[c("cov1",paste0("PC",1:40))]
write.table(gemma_rescaffold_cov,"PCA_loadings_rescaffold_v4_PC40",append = FALSE,quote=FALSE, sep="\t",col.names = FALSE,row.names = FALSE,)

library(rgl)
library(dplyr)

comp<-c(1,2,3)
xlabel = paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)",sep="")
ylabel = paste("PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="")
zlabel = paste("PC",comp[3]," (",signif(eig$val[comp[3]], digits=3)*100,"%)",sep="")


# Adding column based on other column:
PC<-PC %>%mutate(pintar = case_when(
    (Population=="CL") ~ "#7FCC12",
    (Population=="CMR") ~ "yellow",
    (Population=="CMY") ~ "#A38C0E",
    (Population=="HP") ~ "orange",
    (Population=="TB") ~ "red",
    (Population=="PP") ~ "#008000",
    (Population=="DBB") ~ "blue",
    (Population=="DBI") ~  "purple",
    (Population=="DBR") ~ "magenta4",
    (Population=="RA") ~ "dodgerblue4",
    (Population=="AG") ~ "black",
    (Population=="SC") ~ "brown"
  ))

open3d()
plot3d( 
  x=PC$PC1, y=PC$PC2, z=PC$PC3, 
  col = PC$pintar, 
  radius = .05,
  size="10",
  type="p",
  xlab=xlabel, ylab=ylabel, zlab=zlabel)



grid3d(side = "x")
grid3d(side = "y")
grid3d(side = "z")
rgl.postscript("Pca123grid.eps","eps")
rgl.postscript("Pca123.pdf","pdf")



clustfile<-reslabels[c("V1","V1","Population")]

clustfile<-reslabels[c("V1","V1","code")]
clustfile$code[clustfile$code=="AG"]<-"SK"
clustfile$code[clustfile$code=="PP"]<-"PO"
clustfile$code[clustfile$V1=="PP34_S265.r.bam"]<-"PD"
clustfile$code[clustfile$V1=="PP35_S269.r.bam"]<-"PD"
clustfile$code[clustfile$V1=="PP36_S254.r.bam"]<-"PD"

write.table(clustfile,"pumilio.clust3",append = FALSE,quote=FALSE, sep="\t",col.names = FALSE,row.names = FALSE)


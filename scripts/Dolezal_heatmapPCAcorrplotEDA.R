#---- libraries & settings -------

library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)#needed for rotate_x_text(45)

library(RColorBrewer)
library(svglite)
library(stringr) #needed for grep colnames syntax & str_split
library(factoextra) #PCA
library(ggrepel)# needed for geom_label_repel
library(pheatmap)#Heatmap
library(corrplot)
library(variancePartition)
library(BiocParallel)
library(lme4)
library(lmerTest)
library(emmeans)

getRversion()
sessionInfo()
setwd("your//path//")
getwd()
path="your//path//"

date="jan2025"
dir.create(paste0("path","//results_",date,"//"))
resultsDir<-paste0("path","//results_",date,"//")
resultsDir


#------- colour palettes -------------
#display.brewer.all(colorblindFriendly = T)
#display.brewer.pal(12, name = "Paired")

palette.bin<-brewer.pal(n=11,name="RdBu")[c(1,3,5,7)]
palette.bin
"#67001F" "#D6604D" "#FDDBC7" "#D1E5F0"


#xdata has flies in rows, bodylength measures in columns, 
# in size info,numbered (categorical) and sieve size (numeric) in columns # bin.cols

# xdata columns: flyID,thorax.r1,thorax.r2,wing.r1,wing.r2,bin,sieve.size

# bin.cols: 6:7

metric.variables<-c("thorax.r1", 
                    "thorax.r2",
                    "wing.r1",
                    "wing.r2")    

#-------- violin & boxplots plots for numeric data -------
for (i in metric.variables)
{
  print(i)
  p<-ggplot2::ggplot(data=xdata,
                     aes(x=bin,
                         y=xdata[,i])) +
    geom_boxplot() +
    geom_violin() +
    geom_point(pch = 19,
               size=2, 
               position = position_jitterdodge())+
    scale_fill_manual(values=palette.bin)+
    
    theme_classic()  +
    theme(axis.text = element_text(size = 10, angle=0, hjust=1),
          text = element_text(size = 10))+
    labs(x="",y=i,title="")
  print(p)
  ggsave(file=paste0(resultsDir,"ViolinBoxplot_",i,".svg"),
         height=8,
         width=8)
}


#------------ correlation plot ------------
#example https://www.nature.com/articles/s41598-023-49151-9#Sec2 Figure 5

M<-cor(xdata[,-c(1,bin.cols)],use="pairwise.complete.obs", method="pearson")
M

cormatrix<-cor.mtest(xdata[,-c(1)]) #perform cor.test for a matrix
fdrs<-matrix(p.adjust(cormatrix$p,method="fdr"),byrow = T,
             nrow = 13, ncol = 13)
colnames(fdrs)<-metric.variables
rownames(fdrs)<-metric.variables

cormatrix$p.fdr<-fdrs
cormatrix$pearson<-M
cormatrix$pearson
cormatrix$lowCI
cormatrix$uppCI
cormatrix$p.fdr
cormatrix$p

svg(paste(resultsDir_corrplots,"Corrplot.svg",sep=""),
    width=8, height=9)
corrplot(M,p.mat=cormatrix$fdr,
         type="lower",
         tl.col = "black",
         tl.srt = 45,
         method="circle",
         number.cex = 0.8,
         sig.level=0.05,
         insig="blank",
         order="hclust",
         hclust.method="ward.D2",
         diag = FALSE)$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2))
dev.off()





#--------- pairs plot function -------------------
pairsplot<-function(dataset,phenotypestring,textsize,plotsize,phenogroup,colorby,colorpalette,cortestmethod)
{
  print(phenotypestring)
  dat<-dataset
  
  upper.panel<-function(x,y)
  {
    points(x,y, pch = 19, col = colorpalette[dat[,colorby]])
  }
  panel.cor <- function(x, y)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y,
                   method=cortestmethod,
                   use="pairwise.complete.obs"), digits=2)
    txt <- paste0("R = ", r)
    text(0.5, 0.5, txt, cex=textsize)
  }
  svg(paste(resultsDir,"PairsPlot_",phenogroup,".svg",sep=""),
      width=plotsize, height=plotsize)
  pairs(as.formula(phenotypestring),
        data=dataset, 
        lower.panel = panel.cor,
        upper.panel = upper.panel)
  dev.off()
}

vars<-colnames(xdata[-c(1)])
vars

varselection <- paste0('"', paste(vars , collapse='+'), '"')
print(varselection, quote=F)
varselection

pairsplot(xdata,
          "~ thorax.r1+thorax.r2+wing.r1+wing.r2",
          1.2,30,
          "all.raw.measures",
          "bin",
          "palette.bin",
          "spearman")




#---------- heatmap -----
#example https://www.nature.com/articles/s41598-023-49151-9#Sec2 Figure 4
#optional annotation for rows. 
anno_rows<-data.frame(phenotype =factor(c(rep("thorax",2),
                                          rep("wing",2))))
anno_rows

colnames(xdata_residuals)

rownames(anno_rows) = c(
  "thorax.r1",
  "thorax.r2",
  "wing.r1",
  "wing.r2")

#define the variables to be used for column annotation
anno_cols<-metadata[,c("bin","person"),drop=F]# drop=F, if you specify only one variable for annotation
rownames(anno_cols)<-metadata$flyID # variable to be printed at the end of each column

#specify the colours to be used for annotation, must be a list
anno_colours<-list(
  person=c(Siraj="#A6CEE3",
           Misa="#B2DF8A"),
  bin = c("1"=palette.bin[1],
          "2"=palette.bin[2],
          "3"=palette.bin[3],
          "4"=palette.bin[4]))


pheatmapfunction<-function(dataset,scaling,distance,title,plot_width,plot_height)
{
  exportfilename<-paste0(resultsDir,"Heatmap.",title,".",scaling,".",distance,".svg")
  print(exportfilename)
  svg(exportfilename, width = plot_width, height =plot_height)
  par(mar=c(5.1,10.1,4.1,2.1))
  #  bottom, left, top and right margins
  pheatmap_sorted_raw<-pheatmap::pheatmap(#na.omit(dataset[,-c(1:2)]),
    dataset[,-c(1,bin_cols)],
    annotation_col = anno_cols, 
    annotation_colors = anno_colours,
    scale=scaling,
    na_col="black",#color for missing data
    fontsize=10,
    cutree_cols = 13,
    cutree_rows = 3,
    show_rownames = F,
    main="",
    #main=paste(title,scaling,distance,sep=" "),
    clustering_distance_rows = distance,
    clustering_distance_cols=distance,
    clustering_method="ward.D2",
    cluster_cols = T,
    cluster_rows = T,
    cellheight=15)

  print(pheatmap_sorted_raw)
  dev.off()
}

#distance can be "correlation" or "euclidean"
pheatmapfunction(xdata,"row","correlation","time.human.contact",13,12)

#--------------------  PCA -------------------------------
# https://www.sthda.com/english/wiki/wiki.php?id_contents=7891
# https://rpkgs.datanovia.com/factoextra/
#example https://www.nature.com/articles/s41598-023-49151-9#Sec2 Figure 3


PCA <- function(dataset,columns_numeric,columns_meta,axis1,axis2,colorby,colorvec,whichpheno,type)
{
  #dat<<-dataset[complete.cases(dataset),]
  dat<<-dataset
  print(colorby)
  print(colorvec)
  meta<<-dat[,columns_meta]
  pcadata<<-dat[,columns_numeric]
  pca_result <<- prcomp((pcadata),
                        center = TRUE,
                        scale. = TRUE) 
  # #--------Individual Plot PCs-----------
  indplot<-fviz_pca_ind(pca_result,
                        axes = c(axis1, axis2),
                        #title = "",
                        title=paste(whichpheno,type,sep=" "),
                        #subtitle = ,
                        #legend.title = legend_title,
                        pointsize = 5,
                        pointshape = 21,
                        addEllipses=F,
                        label="ind",
                        #label="none",
                        #legend="none",
                        fill.ind =meta[[colorby]],
                        mean.point = FALSE, #do not plot group mean
                        palette = colorvec,#vector of colors
                        repel = TRUE,
                        max.overlaps=50)
  print(indplot +
          theme(axis.title.x = element_text(size=14),
                axis.title.y = element_text(size=14),
                axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12)) +
          scale_shape_manual(values=c(0,1,2,6)))
  
  ggsave(paste(resultsDir,"IndPlot_PC",axis1,axis2,"_","colorby_",colorby,"_",whichpheno,"_",type,
               "_withlabels.svg",sep=""),width = 10, height = 10)
  
  
  biplot<-fviz_pca_biplot(pca_result, 
                          axes = c(axis1,axis2),
                          #title=paste(whichpheno,type,sep=" "),
                          title = "",
                          subtitle = "",
                          geom.ind = "point",
                          fill.ind =meta[[colorby]],
                          mean.point = FALSE,
                          pointshape = 21,
                          pointsize = 5,
                          labelsize = 5,
                          palette = colorvec,
                          addEllipses = F,
                          # Variables
                          col.var = "black",
                          arrowsize=1.1,
                          legend.title = "",
                          legend="",
                          repel=TRUE)
  print(biplot +
          theme(axis.title.x = element_text(size=14),
                axis.title.y = element_text(size=14),
                axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12)) +
          scale_shape_manual(values=c(0,1,2,6))
  )
  
  ggsave(paste0(resultsDir,"Biplot_PC",axis1,axis2,".","colorby.",colorby,".",whichpheno,".",type,
               ".withoutlabels.svg"),width = 8, height = 8)
  
}
PCA(xdata,c(2:5),
    c(6:7),1,2,
    "bin",palette.bin,
    "all",
    "raw")




# ----------hypothesis testing via LM --------------
variables<-colnames(xdata[,-c(1,bin.cols)])
# xdata has flies in rows variables in columns, first column is flyID 

for (i in variables)
{
  print(i)
  pseudosize<-xdata[,i]
  lm.result<-lm(pseudosize~bin, data=xdata)
  print(summary(lm.result))
  print(emmeans::emmeans(lm.result,pairwise~bin))
  
  #print(emmeans::emmeans(lmer_res,pairwise~time.human.contact))
  xdata<-wide.meta.NA[,1:25]
  xdata$dummy<-dummy
  sel=!is.na(xdata[,i])# vector of T and F for non missing rows
  #print(sel)
  xdata$predicted[sel]<- predict(lm.result)
  p<-ggplot2::ggplot(data=xdata,
                     aes(x=bin,
                         y=predicted,
                         color=bin))+
    geom_point(size=4)+
    facet_wrap(~bin,ncol=4)+
    scale_colour_manual(values=palette.bin)+
    theme(legend.position="none") +
    #add the raw data as separate layer
    geom_point(data=xdata,
               aes(x=bin,
                   y=pseudosize,
                   color=bin),
               size=2,
               shape=8,
               inherit.aes = FALSE,
               position = position_dodge(width= 3))+
    facet_wrap(~bin,ncol=4)+
    scale_colour_manual(values=palette.bin) +
    ylab("observed and predicted") +
    xlab("sieve") +
    theme(legend.title=element_blank(),
          legend.text=element_text(size=22)) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    theme_classic(base_size = 22) +
    theme(legend.title=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.text.x=element_blank() )+
    ggtitle(i)
  print(p)
}

#---------- citations -----
getRversion()
citation("dplyr")

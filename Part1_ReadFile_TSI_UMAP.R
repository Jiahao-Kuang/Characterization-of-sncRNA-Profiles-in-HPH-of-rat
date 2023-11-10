#TSI tissue specific index 组织特异性参数：
#range from 0 to 1 ;https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4856985/
# single tissues (TSI > 0.85), housekeeping miRNAs (TSI < 0.5). 
if (F) {
  #To evaluate the variability of expression patterns, 
  #we calculated a tissue specificity index (TSI) for 
  # each miRNA analogously to the TSI ‘tau’ for mRNAs 
  # originally developed by Yanai et al. (7). This 
  # specificity index is a quantitative, graded scalar 
  # measure for the specificity of expression of a miRNA 
  # with respect to different organs. The values range from 0 to 1, 
  # with scores close to 0 represent miRNAs expressed in many or all tissues (i.e. housekeepers)
  # and scores close to 1 miRNAs expressed in only one specific tissue (i.e. tissue-specific miRNAs).
}#more information

TSI <- function(x,n){
  max.exp <- rowMaxs(as.matrix(x))
  min.exp <- rowMins(as.matrix(x))
  x.j.i <- (x-min.exp)/(max.exp-min.exp)
  x.j.i <- 1-x.j.i
  done <- na.omit(rowSums(x.j.i)/(n-1))
  return(done)
}

##最大最小归一化
mmscale <- function(x,n){
  max.exp <- rowMaxs(as.matrix(x))
  min.exp <- rowMins(as.matrix(x))
  x.j.i <- (x-min.exp)/(max.exp-min.exp)
}

# ##统计每种RNA占比
# readxl::read_xlsx()


#####
library(DESeq2)
#rs PART
setwd("C:/Users/colet/Desktop/re snc atlas/quantified/")

library(tidyverse)
library(readxl)
library(writexl)
library(readxl)
#读取sample id
sample <- read_xlsx("sampleID.xlsx")
#read.table(file = "sampleid.txt",header = TRUE)
#sort files
if (T) {
  ##read files as a large list
  listys=list.files(pattern = "ys")
  listrs=list.files(pattern = "rs")
  
  ##rename every samples with readable names
  ys.in = map(listys, ~ read.table(.,header = F))
  rs.in = map(listrs, ~ read.table(.,header = F))
  
  aa = 1
  i = 1
  
  for (i in 1:48) {
    ys.in[[i]]$V1 <- gsub(ys.in[[i]]$V1,pattern = "rno*",replacement = "Rny")
  }
  
  
  ys.in.b=data.frame(id=ys.in[[1]][,1])
  
  
  for (i in 1:48) {
    spn <- sample$readable_name[sample$sample_name==gsub("-",".", paste(strsplit(listys[i],"_")[[1]][1],strsplit(listys[i],"_")[[1]][2],sep = "_"))]
    colnames(ys.in[[i]]) <- c("id",spn) 
    names(ys.in[i]) <- spn
    ys.in.b <- merge(ys.in.b,ys.in[[i]],by="id")
  }
  ys.merge <- data.frame(ys.in.b[,-1],row.names = ys.in.b[,1])
  
  ys.df=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  i=1
  for (i in 1:8) {
    ys.df[[i]] <- ys.merge[,c(i,i+8,i+16,i+24,i+32,i+40)]
  }
  
  ##ysRNA处理成对照与处理组分开的表格
  ys.df.nor=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  ys.df.hyp=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  ys.df.nor.all=data.frame(rep(NA,nrow(ys.df[[1]])))
  ys.df.hyp.all=data.frame(rep(NA,nrow(ys.df[[1]])))
  ys.df.all=data.frame(rep(NA,nrow(ys.df[[1]])))
  for (i in 1:8) {
    ys.df.nor[[i]] <- ys.df[[i]][,c(1,2,3)]
    ys.df.hyp[[i]] <- ys.df[[i]][,c(4,5,6)]
    ys.df.nor.all <- cbind(ys.df.nor.all,ys.df.nor[[i]])
    ys.df.hyp.all <- cbind(ys.df.hyp.all,ys.df.hyp[[i]])
    ys.df.all <- cbind(ys.df.all,ys.df[[i]])
  }
  ys.df.all <- ys.df.all[,-1]
  ys.df.nor.all <- ys.df.nor.all[,-1]
  ys.df.hyp.all <- ys.df.hyp.all[,-1]
  
  ##火山图
  { 
    rny.df.all <- cbind(ys.df.nor.all,ys.df.hyp.all)
    
    df.col.1 <- data.frame(id=colnames(rny.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                        "Liver","Lung","Intestines",
                                                                        "Heart","Spleen"),each=3)),2)),
                           condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
    i=7
    a=1
    p.rny.vol <- NULL  
    for (i in c(seq(1,22,3))) {
      df.col.tem <- df.col.1[c(i,i+1,i+2,i+24,i+25,i+26),]
      
      df.col.tem$condition <- relevel(as.factor(df.col.tem$condition), ref = "Normoxia")
      
      rny.df.all.tem <- rny.df.all[,c(i,i+1,i+2,i+24,i+25,i+26)]
      dds.1 <- DESeqDataSetFromMatrix(countData = data.frame(rny.df.all.tem),
                                      colData = df.col.tem%>%dplyr::select(-tissue), design = ~condition)
      res.tem <- na.omit(data.frame(results(DESeq(dds.1),pAdjustMethod = "BH")))
      
      res.vol.1 <- data.frame(res.tem)%>%filter(baseMean>1)
      
      res.vol.1 <- rownames_to_column(res.vol.1,var = "id")
      
      ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
             ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                    "Both")) -> res.vol.1$exin
      ifelse(res.vol.1$pvalue < .05,1.1,1) -> res.vol.1$size
      
      res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(pvalue)>1)->id.vol
      
      library(ggrepel)
      ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(pvalue)))+
        geom_point(aes(size=-log10(pvalue), color=exin))+
        scale_size_continuous(range = c(0,2.5),
                              # breaks = seq(0,10,2)
        )+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(pvalue)"))+
        geom_hline(yintercept = log10(10),linetype=5)+
        geom_vline(xintercept = .5,color="#f0a1a8")+
        geom_vline(xintercept = -.5,color="#4994c4")+
        theme_minimal()+
        geom_text_repel(data = id.vol%>%filter(-log10(pvalue)>1.5),aes(x=log2FoldChange,y=-log10(pvalue),label = id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
        labs(title = as.character(unique(df.col.tem$tissue)))+
        theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
        # theme(legend.position = legend.justification=c(0,0))+
        scale_y_continuous(limits = c(0,8))+
        scale_x_continuous(limits = c(-2,2))->vol.rny
      vol.rny
      
      ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
        geom_point(aes(size=-log10(pvalue), color=exin))+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
        geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
        theme_light()+
        geom_hline(yintercept = .5,color="#f0a1a8")+
        geom_hline(yintercept = -.5,color="#4994c4")+
        scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
        theme(legend.position = "none",
              # axis.text = element_text(size = 10),
              # axis.title = element_text(size=10),
              # axis.line = element_line(color = "black"),
              # plot.background = element_rect(color="black",size=1),
              panel.background = element_rect(color = "white"),
              panel.grid = element_blank())->vol.2.rny
      vol.2.rny
      
      vol.rny+
        # theme(legend.position = "none")+
        annotation_custom(
          grob = ggplotGrob(vol.2.rny),
          xmin = -2.15,
          xmax = -.10,
          ymin = 4.6,
          ymax = 8.5
        )->uni.vol.rny
      uni.vol.rny
      
      p.rny.vol[[a]] <- uni.vol.rny
      a <- a+1
    }
    library(patchwork)
    
    p.rny.vol[[1]]+p.rny.vol[[2]]+p.rny.vol[[3]]+p.rny.vol[[4]]+
      p.rny.vol[[5]]+p.rny.vol[[6]]+p.rny.vol[[7]]+p.rny.vol[[8]]->a
    
    library(cowplot)
    plot_grid(p.rny.vol[[1]],p.rny.vol[[2]],p.rny.vol[[3]],p.rny.vol[[4]],
              p.rny.vol[[5]],p.rny.vol[[6]],p.rny.vol[[7]],p.rny.vol[[8]],
              ncol = 4,
              label_size = 20,
              labels = "AUTO")->a
    
    # ggsave(a,filename = "rny.vol.patch.pdf",path = "C:\\Users\\colet\\Desktop\\re snc atlas",units = "cm",width = 60,height = 32)
    
    ##火山图
    
    
    dds<-DESeqDataSetFromMatrix(countData = data.frame(rny.df.all), colData = df.col.1, design = ~tissue+condition)
    
    
    
    # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
    dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
    res = results(dds, pAdjustMethod = "BH")
    res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
    
    
    
    
    
    res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
    library(data.table)
    res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
    ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
           ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                  "Both")) -> res.vol.1$exin
    ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
    
    res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(padj)>1)->id.vol
    
    library(ggrepel)
    ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
      geom_point(aes(size=-log10(padj), color=exin))+
      scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
      scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
      geom_hline(yintercept = log10(10),linetype=5)+
      geom_vline(xintercept = .5,color="#f0a1a8")+
      geom_vline(xintercept = -.5,color="#4994c4")+
      theme_minimal()+
      geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "rnyNA")+
      theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
      # theme(legend.position = legend.justification=c(0,0))+
      scale_y_continuous(limits = c(0,13.5))+
      scale_x_continuous(limits = c(-2,2))->vol.rny
    vol.rny
    
    ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
      geom_point(aes(size=-log10(padj), color=exin))+
      scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
      geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
      theme_light()+
      geom_hline(yintercept = .5,color="#f0a1a8")+
      geom_hline(yintercept = -.5,color="#4994c4")+
      scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
      theme(legend.position = "none",
            # axis.text = element_text(size = 10),
            # axis.title = element_text(size=10),
            # axis.line = element_line(color = "black"),
            # plot.background = element_rect(color="black",size=1),
            panel.background = element_rect(color = "white"),
            panel.grid = element_blank())->vol.2.rny
    vol.2.rny
    
    vol.rny+
      # theme(legend.position = "none")+
      annotation_custom(
        grob = ggplotGrob(vol.2.rny),
        xmin = -2.3,
        xmax = -.10,
        ymin = 7,
        ymax = 14.25
      )->uni.vol.rny
    
    uni.vol.rny
  }##火山图
  
  
  
  ##ys差异分析（对照对对照、处理对处理；组织别）
  if (T) {
    ##0.1all找常氧缺氧比对
    if (T) {
      ##0.1all找常氧缺氧比对
      df.col.a <- data.frame(id=colnames(ys.df.all),tissue=c(rep(c("Thymus","Kidney","Brain",
                                                                   "Liver","Lung","Intestines",
                                                                   "Heart","Spleen"),each=6)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=3),8)))
      
      dds<-DESeqDataSetFromMatrix(countData = ys.df.all, colData = df.col.a,design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
     
      ##火山图
      { 
      res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
      library(data.table)
      res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
      ifelse(res.vol.1$log2FoldChange > 1,"Hypoxia",
             ifelse(res.vol.1$log2FoldChange< -1,"Normoxia",
                    "Both")) -> res.vol.1$exin
      ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
      
      res.vol.1%>%filter(abs(log2FoldChange)>1&-log10(padj)>3)->id.vol
      
      library(ggrepel)
      ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
        geom_point(aes(size=-log10(padj), color=exin))+
        scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
        geom_hline(yintercept = log10(10),linetype=5)+
        geom_vline(xintercept = 1,color="#f0a1a8")+
        geom_vline(xintercept = -1,color="#4994c4")+
        theme_minimal()+
        geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id))+
        labs(title = "Rny")+
        theme(legend.justification=c(1,0), legend.position=c(1,0.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
        # theme(legend.position = legend.justification=c(0,0))+
        scale_y_continuous(limits = c(0,13.5))+
        scale_x_continuous(limits = c(-3,3))->vol.rny
      vol.rny
      
      ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
        geom_point(aes(size=-log10(padj), color=exin))+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
        geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
        theme_light()+
        scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
        theme(legend.position = "none",
              # axis.text = element_text(size = 10),
              # axis.title = element_text(size=10),
              # axis.line = element_line(color = "black"),
              # plot.background = element_rect(color="black",size=1),
              panel.background = element_rect(color = "white"),
              panel.grid = element_blank())->vol.2.rny
      vol.2.rny
      
      vol.rny+
        # theme(legend.position = "none")+
        annotation_custom(
          grob = ggplotGrob(vol.2.rny),
          xmin = -3.3,
          xmax = -.10,
          ymin = 7,
          ymax = 14.25
        )->uni.vol.rny
      
      uni.vol.rny
    }##火山图
      
      
      ysdf.sig <- ys.df.all[rownames(ys.df.all)%in%rownames(res.df[abs(res.df$log2FoldChange)>1,]),]
      ysdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        ysdf.sig.mean.1  <-  rowMeans(matrix(ysdf.sig[,c(i:i+2)]))
        ysdf.sig.mean <- cbind(ysdf.sig.mean,ysdf.sig.mean.1)
      } 
      colnames(ysdf.sig.mean) <- colnames(ysdf.sig)[seq(1,48,3)]
      # TSI(ysdf.sig)
      ysdf.sig$id <- rownames(ysdf.sig)
      
      library(data.table)
      ysdf.sig.long <- melt(ysdf.sig,id="id")
      name.idx <- strsplit(as.character(ysdf.sig.long$variable),"_")
      for (i in 1:nrow(ysdf.sig.long)) {
        ysdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        ysdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      ysdf.sig.long$condition[ysdf.sig.long$condition=="NOR"] <- "Normoxia"
      ysdf.sig.long$condition[ysdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      ysdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))
      
      ##两个柱子挨在一起的图
      if (T) {
        library(ggplot2)
        ysdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(x=id,y=log10(value)))+
          geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=condition))+
          stat_summary(fun = mean,
                       geom = "errorbar",
                       fun.max = function(x) mean(x) + sd(x),
                       fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
                       aes(group=condition,color=condition),
                       size=.5,
                       width=.25,
          )+
          # geom_curve()
          stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
          facet_wrap(~tissue,ncol=1,strip.position = "right")+
          labs(y="log10(CPM)",x=NULL,title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "2barEXP.rny.pdf",units = "cm",width = 16,height = 20,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子挨在一起的图
      
      ##两个柱子堆积成1
      if (T) {
        library(ggplot2)
        ysdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "fillbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (T) {
        library(ggplot2)
        ysdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "stackbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
      
      # # Create a tibble for LRT results
      # res_LRT_tb <- res %>%
      #   data.frame() %>%
      #   rownames_to_column(var="gene") %>% 
      #   as_tibble()
      # # Subset to return genes with padj < 0.05
      # sigLRT_genes <- res_LRT_tb %>% 
      #   filter(padj < .05)
      # 
      # # Get number of significant genes
      # nrow(sigLRT_genes)
      # 
      # # Subset results for faster cluster finding (for classroom demo purposes)
      # clustering_sig_genes <- sigLRT_genes %>%
      #   arrange(padj) %>%
      #   head(n=1000)
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
      # library(DEGreport)
      # # install.packages('C:/Users/colet/Desktop/lasso2_1.2-22.tar.gz', repo=NULL)
      # clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
      # 
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # 
      # BiocManager::install("DEGreport")
    }
    
    ##0.1all找组织比对(第一个是点图)
    if (T) {
      ##先算TSI
      df.col = data.frame(id = colnames(ys.df.nor.all),dex = c(rep("thymus",3),
                                                               rep("kidney",3),
                                                               rep("brain",3),
                                                               rep("liver",3),
                                                               rep("lung",3),
                                                               rep("intestines",3),
                                                               rep("heart",3),
                                                               rep("spleen",3)
      ))
      dds<-DESeqDataSetFromMatrix(countData = ys.df.nor.all, colData = df.col,design = ~dex)
      ## Prefiltering
      filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
      dds <- dds[!filt,]
      
      ## Perform DESeq2()
      dds = DESeq(dds)
      res = results(dds, pAdjustMethod = "BH")
      
      ##determine tissue specific(calculate TSI)
      # res.tsi <-data.frame(res)
      # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
      # res.sig <- na.omit(data.frame(res.sig))
      # res.sig.ys <- data.frame(merge(data.frame(id=rownames(ys.df.hyp.all),ys.df.nor.all),data.frame(id=rownames(res.sig))
      #                                ,by="id"),row.names = 1)
      
      tpm <- function(x){
        tpmres <- (x[,6:ncol(x)]/x[,5])/colSums(x[,6:ncol(x)]/x[,5])*1e6
        return(tpmres)
      }
      
      library(edgeR)
      # apply variance stabilizing transformation
      # vsted <- cpm(counts(dds))
      v = vst(dds, blind=FALSE)
      vsted = assay(v)
      
      ##calculate TSI
      tsi.ys.nor <- rep(0,nrow(vsted))
      for (i in c(1,4,7,10,13,16,19,21)) {
        tsi.ys.nor <- data.frame(tsi.ys.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
      }
      tsi.ys.nor <- tsi.ys.nor[,-1]
      colnames(tsi.ys.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
      avgexp.ys.nor <- tsi.ys.nor
      tsi.ys.nor <- data.frame(TSI=TSI(tsi.ys.nor,8))
      
      tsi.sig.id <- rownames_to_column(tsi.ys.nor,var="id")%>%filter(TSI>.85)
      
      ##0.1all找组织比对
      df.col.a <- data.frame(id=colnames(ys.df.all),tissue=c(rep(c("Thymus","Kidney","Brain",
                                                                   "Liver","Lung","Intestines",
                                                                   "Heart","Spleen"),each=6)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=3),8)))
      
      dds<-DESeqDataSetFromMatrix(countData = ys.df.all, colData = df.col.a,design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~condition, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
      
      # res.df.1 <- res.df
      # res.df.1$id <- rownames(res.df.1)
      
      ysdf.sig <- ys.df.all[rownames(ys.df.all)%in%rownames(res.df[abs(res.df$log2FoldChange)>1,]),]
      ysdf.sig <- merge(rownames_to_column(ysdf.sig,"id"),tsi.sig.id,by="id")
      ysdf.sig <- column_to_rownames(ysdf.sig,"id")
      # rownames(ysdf.sig) <- ysdf.sig$id
      ysdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        ysdf.sig.mean.1  <-  rowMeans(matrix(ysdf.sig[,c(i:i+2)]))
        ysdf.sig.mean <- cbind(ysdf.sig.mean,ysdf.sig.mean.1)
      } 
      colnames(ysdf.sig.mean) <- colnames(ysdf.sig)[seq(1,48,3)]
      # TSI(ysdf.sig)
      ysdf.sig$id <- rownames(ysdf.sig)
      
      library(data.table)
      ysdf.sig.long <- melt(ysdf.sig,id="id")
      name.idx <- strsplit(as.character(ysdf.sig.long$variable),"_")
      for (i in 1:nrow(ysdf.sig.long)) {
        ysdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        ysdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      ysdf.sig.long$condition[ysdf.sig.long$condition=="NOR"] <- "Normoxia"
      ysdf.sig.long$condition[ysdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      ysdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))->p.dat
      
      
      # dcast(data = p.dat,id~condition+tissue+value_mean+sd+se) 
      ##两个柱子挨在一起的图
      if (T) {
        library(ggplot2)
        # ysdf.sig.long%>%filter(value>0&condition=="Normoxia")%>%
        # p.dat%>%filter(value>0&condition=="Normoxia")%>%
        rny.dot.dat <- data.frame(merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id"),RNA="Rny")
        merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id")%>%
          filter(condition=="Normoxia")%>%
          ggplot(aes(y=id,x=tissue))+
          geom_point(aes(size=log10(value_mean),color=padj))+
          scale_size_continuous(range = c(-1,8),breaks = seq(.5,3,.5))+
          # scale_size_continuous(range = c(0,6),breaks = c(1,2,3))+
          # geom_point(fun = "mean", stat = "summary",aes(group=tissue,fill=condition,color=condition,size=log(value)),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=id,fill=condition,color=condition,size=log(value)))+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
        # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
        # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
        # facet_wrap(~tissue,nrow=1,strip.position = "top")+
        labs(y=NULL,x=NULL,title = "Rny")+
          # scale_color_gradient(low = "blue",high = "red")+
          viridis::scale_color_viridis()+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(size=guide_legend(reverse = F,title = "log10(CPM)"),
                 # color=guide_legend(title = "p.adj")
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(),
                axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
          )->p.rny.tsi.dot
        p.rny.tsi.dot
        # ggsave(filename = "TIS.dot.EXP.rny.pdf",units = "cm",width = 16,height = 8,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子挨在一起的图
      
      ##两个柱子堆积成1
      if (F) {
        library(ggplot2)
        ysdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.fillbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (F) {
        library(ggplot2)
        ysdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.stackbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
    }
    
  }##新加入新分析

  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(ys.df.nor.all),dex = c(rep("thymus",3),
                                                           rep("kidney",3),
                                                           rep("brain",3),
                                                           rep("liver",3),
                                                           rep("lung",3),
                                                           rep("intestines",3),
                                                           rep("heart",3),
                                                           rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = ys.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "BH")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.ys <- data.frame(merge(data.frame(id=rownames(ys.df.hyp.all),ys.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  tpm <- function(x){
    tpmres <- (x[,6:ncol(x)]/x[,5])/colSums(x[,6:ncol(x)]/x[,5])*1e6
    return(tpmres)
  }
  
  library(edgeR)
  # apply variance stabilizing transformation
  # vsted <- cpm(counts(dds))
  v = vst(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.ys.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.ys.nor <- data.frame(tsi.ys.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.ys.nor <- tsi.ys.nor[,-1]
  colnames(tsi.ys.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.ys.nor <- tsi.ys.nor
  tsi.ys.nor <- data.frame(tsi.ys.nor=TSI(tsi.ys.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.ys.nor
  p.pca.ys.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(ys.df.hyp.all),dex = c(rep("thymus",3),
                                                             rep("kidney",3),
                                                             rep("brain",3),
                                                             rep("liver",3),
                                                             rep("lung",3),
                                                             rep("intestines",3),
                                                             rep("heart",3),
                                                             rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = ys.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = vst(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.ys.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.ys.hyp <- data.frame(tsi.ys.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.ys.hyp <- tsi.ys.hyp[,-1]
  colnames(tsi.ys.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.ys.hyp <- tsi.ys.hyp
  tsi.ys.hyp <- data.frame(tsi.ys.hyp=TSI(tsi.ys.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.ys.hyp
  p.pca.ys.hyp
  
  p.pca.ys.nor-p.pca.ys.hyp -> p.pca.ys.nor.hyp.unite
  
  ##determine tissue specific ysRNA
  tsi.ys.nor  <- data.frame(id=rownames(tsi.ys.nor),tsi.ys.nor,row.names = NULL)
  # tsi.ys.nor.spc <- tsi.ys.nor[tsi.ys.nor$tsi.ys.nor>0.85,]
  tsi.ys.hyp  <- data.frame(id=rownames(tsi.ys.hyp),tsi.ys.hyp,row.names = NULL)
  # tsi.ys.hyp.spc <- tsi.ys.hyp[tsi.ys.hyp$tsi.ys.hyp>0.85,]
  
  ####出现了一种两种状况TSI变化在阈值线前后，尝试将TSI相减2022.10.05）
  tsi.ys.nor.1 <- data.frame(tsi.ys.nor,condition="Control")
  tsi.ys.hyp.1 <- data.frame(tsi.ys.hyp,condition="Hypoxia")
  tsi.ys.gap <- merge(tsi.ys.nor.1,tsi.ys.hyp.1,by="id")
  tsi.ys.gap$gap <- abs(tsi.ys.gap$tsi.ys.nor-tsi.ys.gap$tsi.ys.hyp)
  tsi.ys.gap.1 <- tsi.ys.gap[order(tsi.ys.gap$gap,decreasing = T),]
  
  tsi.ys.nor.spc.1 <- tsi.ys.gap.1[tsi.ys.gap.1$tsi.ys.nor>0.85&tsi.ys.gap.1$gap>=.1,]
  tsi.ys.nor.spc <- data.frame(tsi.ys.nor.spc.1[,c(1,2,6)])
  tsi.ys.hyp.spc.1 <- tsi.ys.gap.1[tsi.ys.gap.1$tsi.ys.hyp>0.85&tsi.ys.gap.1$gap>=.1,]
  tsi.ys.hyp.spc <- data.frame(tsi.ys.hyp.spc.1[,c(1,4,6)])
  ########################
  
  ##determine expression of tissue specific ysRNA
  library(data.table)
  library(ggplot2)
  ##ys对照组
  tsi.ys.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.ys.nor),avgexp.ys.nor),tsi.ys.nor.spc,by="id")
  tsi.ys.nor.spc.long <- melt(tsi.ys.nor.spc.long.1,id=c("id","tsi.ys.nor","gap"))
  colnames(tsi.ys.nor.spc.long) <- c("ID","TSI","Gap","Tissue","VST")
  tsi.ys.nor.spc.long <- tsi.ys.nor.spc.long[order(tsi.ys.nor.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.ys.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.ys.nor.spc 
  p.ys.nor.spc 
  
  ##ys处理组
  tsi.ys.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.ys.hyp),avgexp.ys.hyp),tsi.ys.hyp.spc,by="id")
  tsi.ys.hyp.spc.long <- melt(tsi.ys.hyp.spc.long.1,id=c("id","tsi.ys.hyp","gap"))
  colnames(tsi.ys.hyp.spc.long) <- c("ID","TSI","Gap","Tissue","VST")
  tsi.ys.hyp.spc.long <- tsi.ys.hyp.spc.long[order(tsi.ys.hyp.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.ys.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.ys.hyp.spc 
  p.ys.hyp.spc 
  
  tsi.ys.nor.spc.long <- data.frame(tsi.ys.nor.spc.long,condition="Control")
  tsi.ys.hyp.spc.long <- data.frame(tsi.ys.hyp.spc.long,condition="Hypoxia")
  
  tsi.ys.unite.spc.long <- rbind(tsi.ys.nor.spc.long,tsi.ys.hyp.spc.long)
  
  
  ###重新找出不显著的画灰点(2022.09.25后新加)
  tsi.ys.nor.nonspc <- data.frame(tsi.ys.nor[tsi.ys.nor$id%in%tsi.ys.unite.spc.long$ID,], condition = "Control")
  tsi.ys.nor.nonspc.long.1 <- merge(data.frame(id=tsi.ys.nor$id[tsi.ys.nor$id%in%tsi.ys.unite.spc.long$ID],
                                               avgexp.ys.nor[rownames(avgexp.ys.nor)%in%tsi.ys.unite.spc.long$ID,]),
                                    tsi.ys.nor.nonspc,by="id")
  tsi.ys.nor.nonspc.long <- melt(tsi.ys.nor.nonspc.long.1,id=c("id","condition","tsi.ys.nor"))
  
  tsi.ys.hyp.nonspc <- data.frame(tsi.ys.hyp[tsi.ys.hyp$id%in%tsi.ys.unite.spc.long$ID,], condition = "Hypoxia")
  tsi.ys.hyp.nonspc.long.1 <- merge(data.frame(id=tsi.ys.hyp$id[tsi.ys.hyp$id%in%tsi.ys.unite.spc.long$ID],
                                               avgexp.ys.hyp[rownames(avgexp.ys.hyp)%in%tsi.ys.unite.spc.long$ID,]),
                                    tsi.ys.hyp.nonspc,by="id")
  tsi.ys.hyp.nonspc.long <- melt(tsi.ys.hyp.nonspc.long.1,id=c("id","condition","tsi.ys.hyp"))
  
  merge(tsi.ys.hyp.nonspc.long,data.frame(ID=tsi.ys.unite.spc.long$ID,Gap=tsi.ys.unite.spc.long$Gap),by="ID")
  
  
  colnames(tsi.ys.nor.nonspc.long) <- c("ID","condition","TSI","Tissue","VST")
  colnames(tsi.ys.hyp.nonspc.long) <- c("ID","condition","TSI","Tissue","VST")
  
  tem <- merge(tsi.ys.nor.nonspc.long,tsi.ys.hyp.nonspc.long,by="ID")
  abs(tem$TSI.x-tem$TSI.y)[order(abs(tem$TSI.x-tem$TSI.y),decreasing = T)]
  
  tsi.ys.uni.nonspc.long <- rbind(tsi.ys.nor.nonspc.long, tsi.ys.hyp.nonspc.long)
 #######上为2022.09.25新加  
  # tsi.ys.hyp.nonspc <- data.frame(tsi.ys.hyp[tsi.ys.hyp$id%in%tsi.ys.unite.spc.long$ID,], condition = "Hypoxia")
  # tsi.ys.uni.nonspc <- rbind(tsi.ys.hyp.nonspc,tsi.ys.hyp.nonspc)
  tsi.ys.uni.nonspc.long
  
  
  ggplot()+
    geom_point(data = tsi.ys.uni.nonspc.long, aes(x=Tissue,y=ID,size=VST),color="grey50",fill="black",alpha=.5)+
    geom_point(data=tsi.ys.unite.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    geom_raster(data = data.frame(tsi.ys.unite.spc.long[,-6],condition="Hypoxia"),aes(x=9,y=ID,fill=Gap))+
    # geom_vline(xintercept = c(7.5,8.5))+
    facet_wrap(~condition)+
    labs(x="",y="")+
    scale_size_continuous(range = c(-1,6))+
    # , "grey", 
    scale_fill_gradient(high="#2fa1dd", low="white")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
    viridis::scale_color_viridis()->p.ys.unite.spc
   p.ys.unite.spc
   
   
   library(cowplot)
   plot_grid(p.ys.unite.spc, t.1,ncol = 1)
  # ggsave(p.ys.unite.spc,filename = "p.ys.unite.spc.pdf",units = "cm",width = 16,height = 50)
  
}##ysRNA

if (T) {
  ##read files as a large list
  setwd("C:/Users/colet/Desktop/re snc atlas/quantified/")
  listpi=list.files(pattern = "pi")
  
  ##rename every samples with readable names
  pi.in = map(listpi, ~ read.table(.,header = F))
  
  aa = 1
  i = 1
  pi.in.b=data.frame(id=pi.in[[1]][,2])
  
  for (i in 1:48) {
    spn <- sample$readable_name[sample$sample_name==gsub("-",".", paste(strsplit(listpi[i],"_")[[1]][1],strsplit(listpi[i],"_")[[1]][2],sep = "_"))]
    colnames(pi.in[[i]]) <- c(spn,"id") 
    names(pi.in[i]) <- spn
  }
  pi.allid = NULL
  for (i in 1:length(pi.in)) {
    pi.allid <- pi.in[[i]]$id
    pi.allid <- unique(pi.allid)
  }
  pi.in
  pi.in[[49]] <- as.data.frame(data.frame(id=pi.allid))
  pi.merge <- data.frame(reduce(pi.in,right_join),row.names = 2)
  for (i in 1:ncol(pi.merge)) {
    pi.merge[,i][is.na(pi.merge[,i])] <- 0
  }
  
  pi.df=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  i=1
  for (i in 1:8) {
    pi.df[[i]] <- pi.merge[,c(i,i+8,i+16,i+24,i+32,i+40)]
  }
  
  ##piRNA处理成对照与处理组分开的表格
  pi.df.nor=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  pi.df.hyp=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  pi.df.nor.all=data.frame(rep(NA,nrow(pi.df[[1]])))
  pi.df.hyp.all=data.frame(rep(NA,nrow(pi.df[[1]])))
  for (i in 1:8) {
    pi.df.nor[[i]] <- pi.df[[i]][,c(1,2,3)]
    pi.df.hyp[[i]] <- pi.df[[i]][,c(4,5,6)]
    pi.df.nor.all <- cbind(pi.df.nor.all,pi.df.nor[[i]])
    pi.df.hyp.all <- cbind(pi.df.hyp.all,pi.df.hyp[[i]])
  }
  
  pi.df.nor.all <- pi.df.nor.all[,-1]
  pi.df.hyp.all <- pi.df.hyp.all[,-1]
  pir.df.all <- cbind(pi.df.nor.all,pi.df.hyp.all)
  
  
  ##火山图
  { 
    df.col.1 <- data.frame(id=colnames(pir.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                        "Liver","Lung","Intestines",
                                                                        "Heart","Spleen"),each=3)),2)),
                           condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
    i=7
    a=1
    p.pir.vol <- NULL  
    for (i in c(seq(1,22,3))) {
      df.col.tem <- df.col.1[c(i,i+1,i+2,i+24,i+25,i+26),]
      
      df.col.tem$condition <- relevel(as.factor(df.col.tem$condition), ref = "Normoxia")
      
      pir.df.all.tem <- pir.df.all[,c(i,i+1,i+2,i+24,i+25,i+26)]
      dds.1 <- DESeqDataSetFromMatrix(countData = data.frame(pir.df.all.tem),
                                      colData = df.col.tem%>%dplyr::select(-tissue), design = ~condition)
      res.tem <- na.omit(data.frame(results(DESeq(dds.1),pAdjustMethod = "BH")))
      
      res.vol.1 <- data.frame(res.tem)%>%filter(baseMean>1)
      
      res.vol.1 <- rownames_to_column(res.vol.1,var = "id")
      
      ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
             ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                    "Both")) -> res.vol.1$exin
      ifelse(res.vol.1$pvalue < .05,1.1,1) -> res.vol.1$size
      
      res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(pvalue)>1)->id.vol
      
      library(ggrepel)
      ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(pvalue)))+
        geom_point(aes(size=-log10(pvalue), color=exin))+
        scale_size_continuous(range = c(0,2.5),
                              # breaks = seq(0,10,2)
        )+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(pvalue)"))+
        geom_hline(yintercept = log10(10),linetype=5)+
        geom_vline(xintercept = .5,color="#f0a1a8")+
        geom_vline(xintercept = -.5,color="#4994c4")+
        theme_minimal()+
        geom_text_repel(data = id.vol%>%filter(-log10(pvalue)>1.5),aes(x=log2FoldChange,y=-log10(pvalue),label = id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
        labs(title = as.character(unique(df.col.tem$tissue)))+
        theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
        # theme(legend.position = legend.justification=c(0,0))+
        scale_y_continuous(limits = c(0,8))+
        scale_x_continuous(limits = c(-2,2))->vol.pir
      vol.pir
      
      ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
        geom_point(aes(size=-log10(pvalue), color=exin))+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
        geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
        theme_light()+
        geom_hline(yintercept = .5,color="#f0a1a8")+
        geom_hline(yintercept = -.5,color="#4994c4")+
        scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
        theme(legend.position = "none",
              # axis.text = element_text(size = 10),
              # axis.title = element_text(size=10),
              # axis.line = element_line(color = "black"),
              # plot.background = element_rect(color="black",size=1),
              panel.background = element_rect(color = "white"),
              panel.grid = element_blank())->vol.2.pir
      vol.2.pir
      
      vol.pir+
        # theme(legend.position = "none")+
        annotation_custom(
          grob = ggplotGrob(vol.2.pir),
          xmin = -2.15,
          xmax = -.10,
          ymin = 4.6,
          ymax = 8.5
        )->uni.vol.pir
      uni.vol.pir
      
      p.pir.vol[[a]] <- uni.vol.pir
      a <- a+1
    }
    library(patchwork)
    
    p.pir.vol[[1]]+p.pir.vol[[2]]+p.pir.vol[[3]]+p.pir.vol[[4]]+
      p.pir.vol[[5]]+p.pir.vol[[6]]+p.pir.vol[[7]]+p.pir.vol[[8]]->a
    
    library(cowplot)
    plot_grid(p.pir.vol[[1]],p.pir.vol[[2]],p.pir.vol[[3]],p.pir.vol[[4]],
              p.pir.vol[[5]],p.pir.vol[[6]],p.pir.vol[[7]],p.pir.vol[[8]],
              ncol = 4,
              label_size = 20,
              labels = "AUTO")->a
    
    # ggsave(a,filename = "pir.vol.patch.pdf",path = "C:\\Users\\colet\\Desktop\\re snc atlas",units = "cm",width = 60,height = 32)
    
    ##火山图
    
    
    dds<-DESeqDataSetFromMatrix(countData = data.frame(pir.df.all), colData = df.col.1, design = ~tissue+condition)
    
    
    
    # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
    dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
    res = results(dds, pAdjustMethod = "BH")
    res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
    
    
    
    
    
    res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
    library(data.table)
    res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
    ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
           ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                  "Both")) -> res.vol.1$exin
    ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
    
    res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(padj)>1)->id.vol
    
    library(ggrepel)
    ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
      geom_point(aes(size=-log10(padj), color=exin))+
      scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
      scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
      geom_hline(yintercept = log10(10),linetype=5)+
      geom_vline(xintercept = .5,color="#f0a1a8")+
      geom_vline(xintercept = -.5,color="#4994c4")+
      theme_minimal()+
      geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "pirNA")+
      theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
      # theme(legend.position = legend.justification=c(0,0))+
      scale_y_continuous(limits = c(0,13.5))+
      scale_x_continuous(limits = c(-2,2))->vol.pir
    vol.pir
    
    ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
      geom_point(aes(size=-log10(padj), color=exin))+
      scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
      geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
      theme_light()+
      geom_hline(yintercept = .5,color="#f0a1a8")+
      geom_hline(yintercept = -.5,color="#4994c4")+
      scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
      theme(legend.position = "none",
            # axis.text = element_text(size = 10),
            # axis.title = element_text(size=10),
            # axis.line = element_line(color = "black"),
            # plot.background = element_rect(color="black",size=1),
            panel.background = element_rect(color = "white"),
            panel.grid = element_blank())->vol.2.pir
    vol.2.pir
    
    vol.pir+
      # theme(legend.position = "none")+
      annotation_custom(
        grob = ggplotGrob(vol.2.pir),
        xmin = -2.3,
        xmax = -.10,
        ymin = 7,
        ymax = 14.25
      )->uni.vol.pir
    
    uni.vol.pir
  }##火山图
  
  
  ##pi差异分析（对照对对照、处理对处理；组织别）
  
  if (T) {
    ##0.1all找常氧缺氧比对
    if (T) {
      ##0.1all找常氧缺氧比对
      df.col.a <- data.frame(id=colnames(pir.df.all),tissue=c(rep(c("Thymus","Kidney","Brain",
                                                                   "Liver","Lung","Intestines",
                                                                   "Heart","Spleen"),each=6)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=3),8)))
      
      dds<-DESeqDataSetFromMatrix(countData = pir.df.all, colData = df.col.a,design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
      
      ##火山图
      { 
        res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
        library(data.table)
        res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
        ifelse(res.vol.1$log2FoldChange > 1,"Hypoxia",
               ifelse(res.vol.1$log2FoldChange< -1,"Normoxia",
                      "Both")) -> res.vol.1$exin
        ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
        
        res.vol.1%>%filter(abs(log2FoldChange)>1&-log10(padj)>3)->id.vol
        
        library(ggrepel)
        ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
          geom_point(aes(size=-log10(padj), color=exin))+
          scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
          scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
          guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
          geom_hline(yintercept = log10(10),linetype=5)+
          geom_vline(xintercept = 1,color="#f0a1a8")+
          geom_vline(xintercept = -1,color="#4994c4")+
          theme_minimal()+
          geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id))+
          labs(title = "pir")+
          theme(legend.justification=c(1,0), legend.position=c(1,0.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
          # theme(legend.position = legend.justification=c(0,0))+
          scale_y_continuous(limits = c(0,13.5))+
          scale_x_continuous(limits = c(-3,3))->vol.pir
        vol.pir
        
        ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
          geom_point(aes(size=-log10(padj), color=exin))+
          scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
          guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
          geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
          theme_light()+
          scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
          theme(legend.position = "none",
                # axis.text = element_text(size = 10),
                # axis.title = element_text(size=10),
                # axis.line = element_line(color = "black"),
                # plot.background = element_rect(color="black",size=1),
                panel.background = element_rect(color = "white"),
                panel.grid = element_blank())->vol.2.pir
        vol.2.pir
        
        vol.pir+
          # theme(legend.position = "none")+
          annotation_custom(
            grob = ggplotGrob(vol.2.pir),
            xmin = -3.3,
            xmax = -.10,
            ymin = 7,
            ymax = 14.25
          )->uni.vol.pir
        
        uni.vol.pir
      }##火山图
      
      
      pirdf.sig <- pir.df.all[rownames(pir.df.all)%in%rownames(res.df[abs(res.df$log2FoldChange)>1,]),]
      pirdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        pirdf.sig.mean.1  <-  rowMeans(matrix(pirdf.sig[,c(i:i+2)]))
        pirdf.sig.mean <- cbind(pirdf.sig.mean,pirdf.sig.mean.1)
      } 
      
      colnames(pirdf.sig.mean) <- colnames(pirdf.sig)[seq(1,48,3)]
      
      # TSI(pirdf.sig)
      pirdf.sig$id <- rownames(pirdf.sig)
      
      library(data.table)
      pirdf.sig.long <- melt(pirdf.sig,id=c("id","TSI"))
      name.idx <- strsplit(as.character(pirdf.sig.long$variable),"_")
      for (i in 1:nrow(pirdf.sig.long)) {
        pirdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        pirdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      pirdf.sig.long$condition[pirdf.sig.long$condition=="NOR"] <- "Normoxia"
      pirdf.sig.long$condition[pirdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      pirdf.sig.long%>%dplyr::select(-c(variable,TSI))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))
      
      ##两个柱子挨在一起的图
      if (T) {
        library(ggplot2)
        pirdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(x=id,y=log10(value)))+
          geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=condition))+
          stat_summary(fun = mean,
                       geom = "errorbar",
                       fun.max = function(x) mean(x) + sd(x),
                       fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
                       aes(group=condition,color=condition),
                       size=.5,
                       width=.25,
          )+
          # geom_curve()
          stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
          facet_wrap(~tissue,ncol=1,strip.position = "right")+
          labs(y="log10(CPM)",x=NULL,title = "piRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "2barEXP.pir.pdf",units = "cm",width = 16,height = 20,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子挨在一起的图
      
      ##两个柱子堆积成1
      if (T) {
        library(ggplot2)
        pirdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "piRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "fillbarEXP.pir.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (T) {
        library(ggplot2)
        pirdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "piRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "stackbarEXP.pir.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
      
      # # Create a tibble for LRT results
      # res_LRT_tb <- res %>%
      #   data.frame() %>%
      #   rownames_to_column(var="gene") %>% 
      #   as_tibble()
      # # Subset to return genes with padj < 0.05
      # sigLRT_genes <- res_LRT_tb %>% 
      #   filter(padj < .05)
      # 
      # # Get number of significant genes
      # nrow(sigLRT_genes)
      # 
      # # Subset results for faster cluster finding (for classroom demo purposes)
      # clustering_sig_genes <- sigLRT_genes %>%
      #   arrange(padj) %>%
      #   head(n=1000)
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
      # library(DEGreport)
      # # install.packages('C:/Users/colet/Desktop/lasso2_1.2-22.tar.gz', repo=NULL)
      # clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
      # 
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # 
      # BiocManager::install("DEGreport")
    }
    
    ##0.1all找组织比对(第一个是点图)
    if (T) {
      ##先算TSI
      df.col = data.frame(id = colnames(pi.df.nor.all),dex = c(rep("thymus",3),
                                                               rep("kidney",3),
                                                               rep("brain",3),
                                                               rep("liver",3),
                                                               rep("lung",3),
                                                               rep("intestines",3),
                                                               rep("heart",3),
                                                               rep("spleen",3)
      ))
      dds<-DESeqDataSetFromMatrix(countData = pi.df.nor.all, colData = df.col,design = ~dex)
      ## Prefiltering
      filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
      dds <- dds[!filt,]
      
      ## Perform DESeq2()
      dds = DESeq(dds)
      res = results(dds, pAdjustMethod = "BH")
      
      ##determine tissue specific(calculate TSI)
      # res.tsi <-data.frame(res)
      # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
      # res.sig <- na.omit(data.frame(res.sig))
      # res.sig.pir <- data.frame(merge(data.frame(id=rownames(pir.df.hyp.all),pir.df.nor.all),data.frame(id=rownames(res.sig))
      #                                ,by="id"),row.names = 1)
      
      tpm <- function(x){
        tpmres <- (x[,6:ncol(x)]/x[,5])/colSums(x[,6:ncol(x)]/x[,5])*1e6
        return(tpmres)
      }
      
      library(edgeR)
      # apply variance stabilizing transformation
      # vsted <- cpm(counts(dds))
      v = varianceStabilizingTransformation(dds, blind=FALSE)
      vsted = assay(v)
      
      ##calculate TSI
      tsi.pir.nor <- rep(0,nrow(vsted))
      for (i in c(1,4,7,10,13,16,19,21)) {
        tsi.pir.nor <- data.frame(tsi.pir.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
      }
      tsi.pir.nor <- tsi.pir.nor[,-1]
      colnames(tsi.pir.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
      avgexp.pir.nor <- tsi.pir.nor
      tsi.pir.nor <- data.frame(TSI=TSI(tsi.pir.nor,8))
      
      tsi.sig.id <- rownames_to_column(tsi.pir.nor,var="id")%>%filter(TSI>.85)
      
      ##0.1all找组织比对
      df.col.a <- data.frame(id=colnames(pir.df.all),tissue=c(rep(c("Thymus","Kidney","Brain",
                                                                   "Liver","Lung","Intestines",
                                                                   "Heart","Spleen"),each=6)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=3),8)))
      
      dds<-DESeqDataSetFromMatrix(countData = pir.df.all, colData = df.col.a,design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~condition, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
      
      # res.df.1 <- res.df
      # res.df.1$id <- rownames(res.df.1)
      
      pirdf.sig <- pir.df.all[rownames(pir.df.all)%in%rownames(res.df[abs(res.df$log2FoldChange)>1,]),]
      pirdf.sig <- merge(rownames_to_column(pirdf.sig,"id"),tsi.sig.id,by="id")
      pirdf.sig <- column_to_rownames(pirdf.sig,"id")
      # rownames(pirdf.sig) <- pirdf.sig$id
      pirdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        pirdf.sig.mean.1  <-  rowMeans(matrix(pirdf.sig[,c(i:i+2)]))
        pirdf.sig.mean <- cbind(pirdf.sig.mean,pirdf.sig.mean.1)
      } 
      colnames(pirdf.sig.mean) <- colnames(pirdf.sig)[seq(1,48,3)]
      # TSI(pirdf.sig)
      pirdf.sig$id <- rownames(pirdf.sig)
      
      library(data.table)
      pirdf.sig.long <- melt(pirdf.sig,id="id")
      name.idx <- strsplit(as.character(pirdf.sig.long$variable),"_")
      for (i in 1:nrow(pirdf.sig.long)) {
        pirdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        pirdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      pirdf.sig.long$condition[pirdf.sig.long$condition=="NOR"] <- "Normoxia"
      pirdf.sig.long$condition[pirdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      pirdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))->p.dat
      
      
      # dcast(data = p.dat,id~condition+tissue+value_mean+sd+se) 
      ##点图（全部）
      if (T) {
        library(ggplot2)
        # pirdf.sig.long%>%filter(value>0&condition=="Normoxia")%>%
        # p.dat%>%filter(value>0&condition=="Normoxia")%>%
        pir.dot.dat <- data.frame(merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id"),RNA="pir")
        merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id")%>%
          filter(condition=="Normoxia")->pir.dot.dat.1
          
          pir.dot.dat.1%>%filter(value_mean>20)%>%head(30)->pir.dot.dat.2
          pir.dot.dat.2[pir.dot.dat.2$id!="piR-rno-1040",]%>%
          ggplot(aes(y=id,x=tissue))+
          geom_point(aes(size=log10(value_mean),color=padj))+
          scale_size_continuous(range = c(-1,8),breaks = seq(.5,3,.5))+
          # scale_size_continuous(range = c(0,6),breaks = c(1,2,3))+
          # geom_point(fun = "mean", stat = "summary",aes(group=tissue,fill=condition,color=condition,size=log(value)),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=id,fill=condition,color=condition,size=log(value)))+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
        # geom_curve()
        # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
        # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
        # facet_wrap(~tissue,nrow=1,strip.position = "top")+
        labs(y=NULL,x=NULL,title = "piRNA")+
          # scale_color_gradient(low = "blue",high = "red")+
          viridis::scale_color_viridis()+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(size=guide_legend(reverse = F,title = "log10(CPM)"),
                 # color=guide_legend(title = "p.adj")
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(),
                axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
          )->p.pir.tsi.dot
        p.pir.tsi.dot
        # ggsave(filename = "TIS.dot.EXP.pir.pdf",units = "cm",width = 16,height = 8,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##点图（全部）
      
      ##点图（rny+pir）
      if (T) {
        library(ggplot2)
        # pirdf.sig.long%>%filter(value>0&condition=="Normoxia")%>%
        # p.dat%>%filter(value>0&condition=="Normoxia")%>%
        pir.dot.dat <- data.frame(merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id"),RNA="pir")
        merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id")%>%
          filter(condition=="Normoxia")->pir.dot.dat.1
        
        pir.dot.dat.1%>%filter(value_mean>20)%>%head(30)->pir.dot.dat.2
        pir.dot.dat.2[pir.dot.dat.2$id!="piR-rno-1040",]->pir.dot.dat.3
          
        rbind(rny.dot.dat[,-13],pir.dot.dat.3)%>%filter(value_mean>15)->other.dot.dat
        other.dot.dat%>%
        ggplot(aes(y=id,x=tissue))+
          geom_point(aes(size=log10(value_mean),color=padj))+
          scale_size_continuous(range = c(-1,8),breaks = seq(.5,3,.5))+
          # scale_size_continuous(range = c(0,6),breaks = c(1,2,3))+
          # geom_point(fun = "mean", stat = "summary",aes(group=tissue,fill=condition,color=condition,size=log(value)),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=id,fill=condition,color=condition,size=log(value)))+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
        # geom_curve()
        # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
        # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
        # facet_wrap(~tissue,nrow=1,strip.position = "top")+
        labs(y=NULL,x=NULL,title = "other RNA")+
          # scale_color_gradient(low = "blue",high = "red")+
          viridis::scale_color_viridis()+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(size=guide_legend(reverse = F,title = "log10(CPM)"),
                 # color=guide_legend(title = "p.adj")
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(),
                axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
          )->p.other.tsi.dot
        p.other.tsi.dot
        # ggsave(filename = "TIS.dot.EXP.pir.pdf",units = "cm",width = 16,height = 8,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##点图（rny+pir）
      
      
      ##两个柱子堆积成1
      if (F) {
        library(ggplot2)
        pirdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "pir")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.fillbarEXP.pir.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (F) {
        library(ggplot2)
        pirdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "piRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.stackbarEXP.pir.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
    }
    
  }##新加入新分析
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(pi.df.nor.all),dex = c(rep("thymus",3),
                                                           rep("kidney",3),
                                                           rep("brain",3),
                                                           rep("liver",3),
                                                           rep("lung",3),
                                                           rep("intestines",3),
                                                           rep("heart",3),
                                                           rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = pi.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.pi <- data.frame(merge(data.frame(id=rownames(pi.df.hyp.all),pi.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  ## apply variance stabilizing transformation
  #v = vst(dds, blind=FALSE)
  v = varianceStabilizingTransformation(dds,blind = FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.pi.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.pi.nor <- data.frame(tsi.pi.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.pi.nor <- tsi.pi.nor[,-1]
  colnames(tsi.pi.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.pi.nor <- tsi.pi.nor
  tsi.pi.nor <- data.frame(tsi.pi.nor=TSI(tsi.pi.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.pi.nor
  p.pca.pi.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(pi.df.hyp.all),dex = c(rep("thymus",3),
                                                             rep("kidney",3),
                                                             rep("brain",3),
                                                             rep("liver",3),
                                                             rep("lung",3),
                                                             rep("intestines",3),
                                                             rep("heart",3),
                                                             rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = pi.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  #v = vst(dds, blind=FALSE)
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.pi.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.pi.hyp <- data.frame(tsi.pi.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.pi.hyp <- tsi.pi.hyp[,-1]
  colnames(tsi.pi.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.pi.hyp <- tsi.pi.hyp
  tsi.pi.hyp <- data.frame(tsi.pi.hyp=TSI(tsi.pi.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.pi.hyp
  p.pca.pi.hyp
  
  p.pca.pi.nor-p.pca.pi.hyp -> p.pca.pi.nor.hyp.unite
  
  ##determine tissue specific piRNA
  tsi.pi.nor  <- data.frame(id=rownames(tsi.pi.nor),tsi.pi.nor,row.names = NULL)
  tsi.pi.nor.spc <- tsi.pi.nor[tsi.pi.nor$tsi.pi.nor>0.85,]
  tsi.pi.hyp  <- data.frame(id=rownames(tsi.pi.hyp),tsi.pi.hyp,row.names = NULL)
  tsi.pi.hyp.spc <- tsi.pi.hyp[tsi.pi.hyp$tsi.pi.hyp>0.85,]
  
  ##determine expression of tissue specific piRNA
  library(data.table)
  library(ggplot2)
  ##pi对照组
  tsi.pi.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.pi.nor),avgexp.pi.nor),tsi.pi.nor.spc,by="id")
  tsi.pi.nor.spc.long <- melt(tsi.pi.nor.spc.long.1,id=c("id","tsi.pi.nor"))
  colnames(tsi.pi.nor.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.pi.nor.spc.long <- tsi.pi.nor.spc.long[order(tsi.pi.nor.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.pi.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.pi.nor.spc 
  p.pi.nor.spc 
  
  ##pi处理组
  tsi.pi.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.pi.hyp),avgexp.pi.hyp),tsi.pi.hyp.spc,by="id")
  tsi.pi.hyp.spc.long <- melt(tsi.pi.hyp.spc.long.1,id=c("id","tsi.pi.hyp"))
  colnames(tsi.pi.hyp.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.pi.hyp.spc.long <- tsi.pi.hyp.spc.long[order(tsi.pi.hyp.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.pi.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.pi.hyp.spc 
  p.pi.hyp.spc 
  
  tsi.pi.nor.spc.long = data.frame(tsi.pi.nor.spc.long,condition="Control")
  tsi.pi.hyp.spc.long = data.frame(tsi.pi.hyp.spc.long,condition="Hypoxia")
  tsi.pi.combine.spc.long <- rbind(tsi.pi.nor.spc.long,tsi.pi.hyp.spc.long)
  
  # tsi.pi.combine.spc.long <- rbind(data.frame(tsi.pi.nor.spc.long,Contition="Control"),data.frame(tsi.pi.hyp.spc.long,Conditoin="Hypoxia"))
  ggplot()+
    geom_point(data=tsi.pi.combine.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    scale_size_continuous(range = c(-1,6))+
    labs(x="",y="")+
    facet_wrap(~condition)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
    viridis::scale_color_viridis()->p.pi.combine.spc 
  p.pi.combine.spc
  # ggsave(filename = "p.pi.combine.spc.pdf",p.pi.combine.spc ,units = "cm",width = 16,height = 16)
  
}##piRNA

if (T) {
  ##read files as a large list
  listrs=list.files(pattern = "rs")
  
  ##rename every samples with readable names
  rs.in = map(listrs, ~ read.table(.,header = F))
  
  aa = 1
  i = 1
  rs.in.b=data.frame(id=rs.in[[1]][,1])
  
  for (i in 1:48) {
    spn <- sample$readable_name[sample$sample_name==gsub("-",".", paste(strsplit(listrs[i],"_")[[1]][1],strsplit(listrs[i],"_")[[1]][2],sep = "_"))]
    colnames(rs.in[[i]]) <- c("id",spn) 
    names(rs.in[i]) <- spn
    rs.in.b <- merge(rs.in.b,rs.in[[i]],by="id")
  }
  rs.merge <- data.frame(rs.in.b[,-1],row.names = rs.in.b[,1])
  
  rs.df=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  i=1
  for (i in 1:8) {
    rs.df[[i]] <- rs.merge[,c(i,i+8,i+16,i+24,i+32,i+40)]
  }
  
  ##rsRNA处理成对照与处理组分开的表格
  rs.df.nor=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  rs.df.hyp=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  rs.df.nor.all=data.frame(rep(NA,nrow(rs.df[[1]])))
  rs.df.hyp.all=data.frame(rep(NA,nrow(rs.df[[1]])))
  for (i in 1:8) {
    rs.df.nor[[i]] <- rs.df[[i]][,c(1,2,3)]
    rs.df.hyp[[i]] <- rs.df[[i]][,c(4,5,6)]
    rs.df.nor.all <- cbind(rs.df.nor.all,rs.df.nor[[i]])
    rs.df.hyp.all <- cbind(rs.df.hyp.all,rs.df.hyp[[i]])
  }
  
  rs.df.nor.all <- rs.df.nor.all[,-1]
  rs.df.hyp.all <- rs.df.hyp.all[,-1]
  ##rs差异分析（对照对对照、处理对处理；组织别）
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(rs.df.nor.all),dex = c(rep("thymus",3),
                                                           rep("kidney",3),
                                                           rep("brain",3),
                                                           rep("liver",3),
                                                           rep("lung",3),
                                                           rep("intestines",3),
                                                           rep("heart",3),
                                                           rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = rs.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.rs <- data.frame(merge(data.frame(id=rownames(rs.df.hyp.all),rs.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.rs.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.rs.nor <- data.frame(tsi.rs.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.rs.nor <- tsi.rs.nor[,-1]
  colnames(tsi.rs.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.rs.nor <- tsi.rs.nor
  tsi.rs.nor <- data.frame(tsi.rs.nor=TSI(tsi.rs.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.rs.nor
  p.pca.rs.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(rs.df.hyp.all),dex = c(rep("thymus",3),
                                                             rep("kidney",3),
                                                             rep("brain",3),
                                                             rep("liver",3),
                                                             rep("lung",3),
                                                             rep("intestines",3),
                                                             rep("heart",3),
                                                             rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = rs.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.rs.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.rs.hyp <- data.frame(tsi.rs.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.rs.hyp <- tsi.rs.hyp[,-1]
  colnames(tsi.rs.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.rs.hyp <- tsi.rs.hyp
  tsi.rs.hyp <- data.frame(tsi.rs.hyp=TSI(tsi.rs.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.rs.hyp
  p.pca.rs.hyp
  
  p.pca.rs.nor-p.pca.rs.hyp -> p.pca.rs.nor.hyp.unite
  
  ##determine tissue specific rsRNA
  tsi.rs.nor  <- data.frame(id=rownames(tsi.rs.nor),tsi.rs.nor,row.names = NULL)
  tsi.rs.nor.spc <- tsi.rs.nor[tsi.rs.nor$tsi.rs.nor>0.85,]
  tsi.rs.hyp  <- data.frame(id=rownames(tsi.rs.hyp),tsi.rs.hyp,row.names = NULL)
  tsi.rs.hyp.spc <- tsi.rs.hyp[tsi.rs.hyp$tsi.rs.hyp>0.85,]
  
  ##determine expression of tissue specific rsRNA
  library(data.table)
  library(ggplot2)
  ##rs对照组
  tsi.rs.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.rs.nor),avgexp.rs.nor),tsi.rs.nor.spc,by="id")
  tsi.rs.nor.spc.long <- melt(tsi.rs.nor.spc.long.1,id=c("id","tsi.rs.nor"))
  colnames(tsi.rs.nor.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.rs.nor.spc.long <- tsi.rs.nor.spc.long[order(tsi.rs.nor.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.rs.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.rs.nor.spc 
  p.rs.nor.spc 
  
  ##rs处理组
  tsi.rs.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.rs.hyp),avgexp.rs.hyp),tsi.rs.hyp.spc,by="id")
  tsi.rs.hyp.spc.long <- melt(tsi.rs.hyp.spc.long.1,id=c("id","tsi.rs.hyp"))
  colnames(tsi.rs.hyp.spc.long) <- c("ID","TSI","Tissue","VST")
  scale_size_continuous(range = c(-1,6))+
  tsi.rs.hyp.spc.long <- tsi.rs.hyp.spc.long[order(tsi.rs.hyp.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.rs.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.rs.hyp.spc 
  p.rs.hyp.spc 
  
  tsi.rs.nor.spc.long <- data.frame(tsi.rs.nor.spc.long,condition="Control")
  tsi.rs.hyp.spc.long <- data.frame(tsi.rs.hyp.spc.long,condition="Hypoxia")
  
  tsi.rs.combine.spc.long <- rbind(tsi.rs.nor.spc.long,tsi.rs.hyp.spc.long)
  
  ggplot()+
    geom_point(data=tsi.rs.combine.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    labs(x="",y="")+
    facet_wrap(~condition)+
    scale_size_continuous(range = c(-1,6))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))+
    viridis::scale_color_viridis()->p.rs.combine.spc 
  p.rs.combine.spc
  
  # ggsave(p.rs.combine.spc,filename = "p.rs.combine.spc.pdf",units = "cm",width = 16,height = 25)
}##rsRNA

if (T) {
  setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  mir.in <- read.table("miR.quantified",header = T)
  
  #去除mir
  mir.in <- mir.in[-grep("mir",mir.in$Geneid),]
  
  mir.in.c=data.frame(mir.in,row.names = 1)
  mir.in.b <- data.frame(mir.in[,7:ncol(mir.in)],row.names = mir.in$Geneid)
  
  i=1
  for (i in 1:ncol(mir.in.b)) {
    colnames(mir.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(mir.in.b),"_")[[i]][1],strsplit(colnames(mir.in.b),"_")[[i]][2],sep = "_")]
  }
  
  #去掉rno
  # paste(strsplit2(rownames(mir.df),"-")[,2],strsplit2(rownames(mir.df),"-")[,3],strsplit2(rownames(mir.df),"-")[,4],sep = "-")
  # library(clusterProfiler)
  # library(org.Rn.eg.db)
  # bitr(geneID = "mir-878",fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Rn.eg.db)
  
  mir.df=mir.in.b
  all.mir.df=mir.in.b
  # write.csv(mir.df,file = "mir.df.csv")
  ##miRNA处理成对照与处理组分开的表格
  mir.df.all <- data.frame(rep(0,nrow(mir.df)))
  mir.df.nor.all = data.frame(rep(0,nrow(mir.df)))
  mir.df.hyp.all = data.frame(rep(0,nrow(mir.df)))
  for (i in 1:8) {
    mir.df.nor.all <- cbind(mir.df.nor.all,mir.df[,c(i+8,i,i+16)])
  }
  i=9
  for (i in 25:32) {
    mir.df.hyp.all <- cbind(mir.df.hyp.all,mir.df[,c(i+8,i,i+16)])
  }
  mir.df.nor.all <- mir.df.nor.all[,-1]
  mir.df.hyp.all <- mir.df.hyp.all[,-1]
  mir.df.all <- cbind(mir.df.nor.all,mir.df.hyp.all)

  #找基因簇
  {
    
    mir.in.1 <- mir.in
    # mir.in.1$Chr <- strsplit2(mir.in.1$Chr,";")[,1]
    # mir.in.1$Start <- strsplit2(mir.in.1$Start,";")[,1]
    # mir.in.1$End <- strsplit2(mir.in.1$End,";")[,1]
    mir.pos <- mir.in.1[,c(1:4)]
    colnames(mir.pos)[1] <- "id"
    id.vol.lung <- id.vol
    id.vol.lung.pos <- merge(mir.pos,id.vol.lung,by="id")
    id.vol.lung.pos[order(id.vol.lung.pos$Chr),]
    #(lung)rno-miR-127-3p,rno-miR-136-5p,rno-miR-3543,	rno-miR-434-3p,	rno-miR-154-5p,rno-miR-541-5p

    # install.packages("gggenes")
    library(ggplot2)
    library(gggenes)
    data(example_genes)
    head(example_genes)
    
    lung <- c("rno-miR-127-3p","rno-miR-136-5p","rno-miR-3543","rno-miR-434-3p","rno-miR-154-5p","rno-miR-541-5p")
    # gene.dat <- id.vol.lung.pos[grep(id.vol.lung.pos$Chr,pattern = "chr6"),]
    
    res.vol.1%>%filter(exin!="Both")->id.brain
    id.brain$id
    merge(mir.pos,id.brain,by="id")->id.brain.1
    brain <- c("rno-miR-379-5p","rno-miR-411-3p","rno-miR-299a-3p","rno-miR-329-3p","rno-miR-494-3p","rno-miR-1193-3p")
    gene.dat <- id.vol.lung.pos[id.vol.lung.pos$id%in%lung,]
    gene.dat <- id.brain.1[id.brain.1$id%in%brain,]
    setwd("C:\\Users\\colet\\Desktop\\re snc atlas")
    write.csv(x = gene.dat,file = "genedat.lung.csv")
    gene.dat.1 <- read.csv(file = "genedat.lung.csv")
    
    ggplot(gene.dat, aes(xmin = Start, xmax = End, y = Chr, fill = id)) + 
      geom_gene_arrow() +  
      scale_fill_brewer(palette = "Set3") +  
      # geom_vline(xintercept = 0)+
      # geom_hline(yintercept = 0)+
      # theme_void()+
      theme(
          # legend.position = "none",
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.line = element_blank())->cluster.down
    cluster.down
    
    tem.brain <- merge(rownames_to_column(mir.df,"id"),data.frame(id=id.brain[,1]),by="id")
    tem.brain.1 <- merge(gene.dat[,c(1:4)],tem.brain,by="id")
    tem.brain.2 <- melt(tem.brain.1,id=c("id","Chr","Start","End"))%>%filter(id%in%brain)
    tem.brain.2$tissue <- str_to_title(strsplit2(tem.brain.2$variable,"_")[,1])
    tem.brain.2$condition <- str_to_title(strsplit2(tem.brain.2$variable,"_")[,2])
    tem.brain.2%>%filter(tissue=="Brain")->tem.brain.3
    
    tem.brain.3$Start <- as.numeric(tem.brain.3$Start)
    library(ggbreak)
    tem.brain.3%>%
    ggplot(aes(x=Start,y=log10(value),color=condition,fill=condition))+
      geom_boxplot(position = "dodge")+
      # geom_point(position_dodge(width = 2))+
      scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
      scale_color_manual(values = c("#f0a1a8","#4994c4"))+
      # guides(fill=guide_legend(reverse = T),size=guide_legend(title = "log10(pvalue)"))+
      theme_classic()+
      # geom_text_repel(data = id.vol%>%filter(-log10(pvalue)>1.5),aes(x=log2FoldChange,y=-log10(pvalue),label = id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "Brain", y="log10(TPM)")+
      theme(legend.justification=c(1,0),
            legend.position=c(1,.5),
            plot.title = element_text(size= 12, face = 2,hjust = .5),
            legend.title = element_text(face = 2,size = 10),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(color="grey"),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank())->cluster.up
    cluster.up
    ggsave(plot = cluster.up,filename = "cluster.up.pdf",units = "cm",width = 10,height = 5,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
    
    cluster.down+
      annotation_custom(
        grob = ggplotGrob(cluster.up),
        xmin = -1.25,
        xmax = 13,
        ymin = 1,
        ymax = 1.5
      )->uni.vol.mir
    uni.vol.mir
      
      # geom_point()
    geom_point(width = .65,position = "dodge")+
      # stat_summary(fun = "mean", geom = "point",aes(group=condition))+
      stat_summary(fun = mean,
                   geom = "errorbar",
                   fun.max = function(x) mean(x) + sd(x),
                   fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
                   aes(group=condition,color=condition),
                   size=.5,
                   width=.25,
      )
    akak
    
    

    
    
    
    
    ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
      geom_gene_arrow() +
      facet_wrap(~molecule,ncol=1)+
      theme_genes()
      facet_wrap(~ molecule, scales = "free", ncol = 1) +
      scale_fill_brewer(palette = "Set3")
    
    
  } #找基因簇
  
  
  ##mir差异分析（对照对对照、处理对处理；组织别）
  
  if (T) {

    ##0.1all找常氧缺氧比对
    if (T) {
      
      ##0.1all找常氧缺氧比对
      # df.col.a <- data.frame(id=colnames(mir.df),tissue=c(rep(c("Thymus","Kidney","Brain",
      #                                                              "Liver","Lung","Intestines",
      #                                                              "Heart","Spleen"),6)),
      #                        condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),1)))
      df.col.1 <- data.frame(id=colnames(mir.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                          "Liver","Lung","Intestines",
                                                                          "Heart","Spleen"),each=3)),2)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
      dds<-DESeqDataSetFromMatrix(countData = data.frame(mir.df.all), colData = df.col.1, design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)

            
      ##火山图
      { 
        df.col.1 <- data.frame(id=colnames(mir.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                            "Liver","Lung","Intestines",
                                                                            "Heart","Spleen"),each=3)),2)),
                               condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
        i=7
        a=3
        p.mir.vol <- NULL  
        for (i in c(seq(1,22,3))) {
          df.col.tem <- df.col.1[c(i,i+1,i+2,i+24,i+25,i+26),]
          
          df.col.tem$condition <- relevel(as.factor(df.col.tem$condition), ref = "Normoxia")
          
          mir.df.all.tem <- mir.df.all[,c(i,i+1,i+2,i+24,i+25,i+26)]
          dds.1 <- DESeqDataSetFromMatrix(countData = data.frame(mir.df.all.tem),
                                        colData = df.col.tem%>%dplyr::select(-tissue), design = ~condition)
          res.tem <- na.omit(data.frame(results(DESeq(dds.1),pAdjustMethod = "BH")))
          
          res.vol.1 <- data.frame(res.tem)%>%filter(baseMean>1)
        
          res.vol.1 <- rownames_to_column(res.vol.1,var = "id")
          
          ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
                 ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                        "Both")) -> res.vol.1$exin
          ifelse(res.vol.1$pvalue < .05,1.1,1) -> res.vol.1$size
          
          res.vol.1%>%filter(abs(log2FoldChange)>.5&pvalue <.05)->id.vol
          
          library(ggrepel)
          ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(pvalue)))+
            geom_point(aes(size=-log10(pvalue), color=exin))+
            scale_size_continuous(range = c(0,2.5),
                                  # breaks = seq(0,10,2)
                                  )+
            scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
            guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(pvalue)"))+
            geom_hline(yintercept = log10(10),linetype=5)+
            geom_vline(xintercept = .5,color="#f0a1a8")+
            geom_vline(xintercept = -.5,color="#4994c4")+
            theme_minimal()+
            geom_text_repel(data = id.vol%>%filter(-log10(pvalue)>1.5),aes(x=log2FoldChange,y=-log10(pvalue),label = id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
            labs(title = as.character(unique(df.col.tem$tissue)))+
            theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
            # theme(legend.position = legend.justification=c(0,0))+
            scale_y_continuous(limits = c(0,8))+
            scale_x_continuous(limits = c(-2,2))->vol.mir
          vol.mir
          
          ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
            geom_point(aes(size=-log10(pvalue), color=exin))+
            scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
            guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
            geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
            theme_light()+
            geom_hline(yintercept = .5,color="#f0a1a8")+
            geom_hline(yintercept = -.5,color="#4994c4")+
            scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
            theme(legend.position = "none",
                  # axis.text = element_text(size = 10),
                  # axis.title = element_text(size=10),
                  # axis.line = element_line(color = "black"),
                  # plot.background = element_rect(color="black",size=1),
                  panel.background = element_rect(color = "white"),
                  panel.grid = element_blank())->vol.2.mir
          vol.2.mir
          
          vol.mir+
            annotation_custom(
              grob = ggplotGrob(vol.2.mir),
              xmin = -2.15,
              xmax = -.10,
              ymin = 4.6,
              ymax = 8.5
            )->uni.vol.mir
          uni.vol.mir
          
          p.mir.vol[[a]] <- uni.vol.mir
          a <- a+1
          }
        library(patchwork)
        
        p.mir.vol[[1]]+p.mir.vol[[2]]+p.mir.vol[[3]]+p.mir.vol[[4]]+
          p.mir.vol[[5]]+p.mir.vol[[6]]+p.mir.vol[[7]]+p.mir.vol[[8]]->a
        
        library(cowplot)
        plot_grid(p.mir.vol[[1]],p.mir.vol[[2]],p.mir.vol[[3]],p.mir.vol[[4]],
                  p.mir.vol[[5]],p.mir.vol[[6]],p.mir.vol[[7]],p.mir.vol[[8]],
                  ncol = 4,
                  label_size = 20,
                  labels = "AUTO")->a
        
        # ggsave(a,filename = "mir.vol.patch.pdf",path = "C:\\Users\\colet\\Desktop\\re snc atlas",units = "cm",width = 60,height = 32)
        
##火山图
        
        
        dds<-DESeqDataSetFromMatrix(countData = data.frame(mir.df.all), colData = df.col.1, design = ~tissue+condition)
        
        
        
        # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
        dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
        res = results(dds, pAdjustMethod = "BH")
        res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
        
        
        
        
        
        res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
        library(data.table)
        res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
        ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
               ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                      "Both")) -> res.vol.1$exin
        ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
        
        res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(padj)>1)->id.vol
        
        library(ggrepel)
        ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
          geom_point(aes(size=-log10(padj), color=exin))+
          scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
          scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
          guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
          geom_hline(yintercept = log10(10),linetype=5)+
          geom_vline(xintercept = .5,color="#f0a1a8")+
          geom_vline(xintercept = -.5,color="#4994c4")+
          theme_minimal()+
          geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
          labs(title = "miRNA")+
          theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
          # theme(legend.position = legend.justification=c(0,0))+
          scale_y_continuous(limits = c(0,13.5))+
          scale_x_continuous(limits = c(-2,2))->vol.mir
        vol.mir
        
        ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
          geom_point(aes(size=-log10(padj), color=exin))+
          scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
          guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
          geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
          theme_light()+
          geom_hline(yintercept = .5,color="#f0a1a8")+
          geom_hline(yintercept = -.5,color="#4994c4")+
          scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
          theme(legend.position = "none",
                # axis.text = element_text(size = 10),
                # axis.title = element_text(size=10),
                # axis.line = element_line(color = "black"),
                # plot.background = element_rect(color="black",size=1),
                panel.background = element_rect(color = "white"),
                panel.grid = element_blank())->vol.2.mir
        vol.2.mir
        
        vol.mir+
          # theme(legend.position = "none")+
          annotation_custom(
            grob = ggplotGrob(vol.2.mir),
            xmin = -2.3,
            xmax = -.10,
            ymin = 7,
            ymax = 14.25
          )->uni.vol.mir
        
        uni.vol.mir
      }##火山图
      
      
      
      mirdf.sig <- mir.df[rownames(mir.df)%in%rownames(res.df[abs(res.df$log2FoldChange)>0,]),]
      mirdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        mirdf.sig.mean.1  <-  rowMeans(matrix(mirdf.sig[,c(i:i+2)]))
        mirdf.sig.mean <- cbind(mirdf.sig.mean,mirdf.sig.mean.1)
      } 
      colnames(mirdf.sig.mean) <- colnames(mirdf.sig)[seq(1,48,3)]
      # TSI(mirdf.sig)
      mirdf.sig$id <- rownames(mirdf.sig)
      
      library(data.table)
      mirdf.sig.long <- melt(mirdf.sig,id="id")
      name.idx <- strsplit(as.character(mirdf.sig.long$variable),"_")
      for (i in 1:nrow(mirdf.sig.long)) {
        mirdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        mirdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      mirdf.sig.long$condition[mirdf.sig.long$condition=="NOR"] <- "Normoxia"
      mirdf.sig.long$condition[mirdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      mirdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))
      
      ##两个柱子挨在一起的图
      if (T) {
        library(ggplot2)
        # mirdf.sig.long%>%filter(value>0)%>%filter(!grepl('let', id))%>%
        mirdf.sig.long%>%filter(value>0&!grepl('let|3068|3084|23|26|30|342|21-',id))%>%
          ggplot(aes(x=id,y=log2(value)))+
          geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=condition))+
          stat_summary(fun = mean,
                       geom = "errorbar",
                       fun.max = function(x) mean(x) + sd(x),
                       fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
                       aes(group=condition,color=condition),
                       size=.5,
                       width=.25,
          )+
          # geom_vline(xintercept = .5)+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
          facet_wrap(~tissue,ncol=1,strip.position = "right")+
          labs(y="log2(CPM)",x=NULL,title = "miRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "2barEXP.mir.pdf",units = "cm",width = 16,height = 20,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子挨在一起的图
      
      ##两个柱子堆积成1
      if (T) {
        library(ggplot2)
        # mirdf.sig.long%>%filter(value>0)%>%
        mirdf.sig.long%>%filter(value>0&!grepl('let|3068|3084|23|26|30|342|29|21-|339|664|10a',id))%>%
          ggplot(aes(y=id,x=log10(value)))+
          
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "miRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          geom_vline(xintercept = .5,linetype = 5)+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "fillbarEXP.mir.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (T) {
        library(ggplot2)
        # mirdf.sig.long%>%filter(value>0)%>%
        mirdf.sig.long%>%filter(value>0&!grepl('let|3068|3084|23|26|30|342|21-',id))%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "miRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )

        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "stackbarEXP.mir.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
      
      # # Create a tibble for LRT results
      # res_LRT_tb <- res %>%
      #   data.frame() %>%
      #   rownames_to_column(var="gene") %>% 
      #   as_tibble()
      # # Subset to return genes with padj < 0.05
      # sigLRT_genes <- res_LRT_tb %>% 
      #   filter(padj < .05)
      # 
      # # Get number of significant genes
      # nrow(sigLRT_genes)
      # 
      # # Subset results for faster cluster finding (for classroom demo purposes)
      # clustering_sig_genes <- sigLRT_genes %>%
      #   arrange(padj) %>%
      #   head(n=1000)
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
      # library(DEGreport)
      # # install.packages('C:/Users/colet/Desktop/lasso2_1.2-22.tar.gz', repo=NULL)
      # clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
      # 
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # 
      # BiocManager::install("DEGreport")
    }
    
    ##0.1all找组织比对(第一个是点图)
    if (T) {
      ##先算TSI
      df.col = data.frame(id = colnames(mir.df.nor.all),dex = c(rep("thymus",3),
                                                               rep("kidney",3),
                                                               rep("brain",3),
                                                               rep("liver",3),
                                                               rep("lung",3),
                                                               rep("intestines",3),
                                                               rep("heart",3),
                                                               rep("spleen",3)
      ))
      dds<-DESeqDataSetFromMatrix(countData = mir.df.nor.all, colData = df.col,design = ~dex)
      ## Prefiltering
      filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
      dds <- dds[!filt,]
      
      ## Perform DESeq2()
      dds = DESeq(dds)
      res = results(dds, pAdjustMethod = "BH")
      
      ##determine tissue specific(calculate TSI)
      # res.tsi <-data.frame(res)
      # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
      # res.sig <- na.omit(data.frame(res.sig))
      # res.sig.mir <- data.frame(merge(data.frame(id=rownames(mir.df.hyp.all),mir.df.nor.all),data.frame(id=rownames(res.sig))
      #                                ,by="id"),row.names = 1)
      
      tpm <- function(x){
        tpmres <- (x[,6:ncol(x)]/x[,5])/colSums(x[,6:ncol(x)]/x[,5])*1e6
        return(tpmres)
      }
      
      library(edgeR)
      # apply variance stabilizing transformation
      # vsted <- cpm(counts(dds))
      v = varianceStabilizingTransformation(dds, blind=FALSE)
      vsted = assay(v)
      
      ##calculate TSI
      tsi.mir.nor <- rep(0,nrow(vsted))
      for (i in c(1,4,7,10,13,16,19,21)) {
        tsi.mir.nor <- data.frame(tsi.mir.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
      }
      tsi.mir.nor <- tsi.mir.nor[,-1]
      colnames(tsi.mir.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
      avgexp.mir.nor <- tsi.mir.nor
      tsi.mir.nor <- data.frame(TSI=TSI(tsi.mir.nor,8))
      
      tsi.sig.id <- rownames_to_column(tsi.mir.nor,var="id")%>%filter(TSI>.85)
      
      ##0.1all找组织比对
      # df.col.a <- data.frame(id=colnames(mir.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
      #                                                              "Liver","Lung","Intestines",
      #                                                              "Heart","Spleen"),each=3)),2)),
      #                        condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),2)))
      df.col.a <- data.frame(id=colnames(mir.df),tissue=c(rep(c("Thymus","Kidney","Brain",
                                                                "Liver","Lung","Intestines",
                                                                "Heart","Spleen"),6)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),1)))
      dds<-DESeqDataSetFromMatrix(countData = mir.df, colData = df.col.a,design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~condition, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
      
      # res.df.1 <- res.df
      # res.df.1$id <- rownames(res.df.1)
      
      mirdf.sig <- mir.df[rownames(mir.df)%in%rownames(res.df[abs(res.df$log2FoldChange)>0,]),]
      mirdf.sig <- merge(rownames_to_column(mirdf.sig,"id"),tsi.sig.id,by="id")
      mirdf.sig <- column_to_rownames(mirdf.sig,"id")
      # rownames(mirdf.sig) <- mirdf.sig$id
      mirdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        mirdf.sig.mean.1  <-  rowMeans(matrix(mirdf.sig[,c(i:i+2)]))
        mirdf.sig.mean <- cbind(mirdf.sig.mean,mirdf.sig.mean.1)
      } 
      colnames(mirdf.sig.mean) <- colnames(mirdf.sig)[seq(1,48,3)]
      # TSI(mirdf.sig)
      mirdf.sig$id <- rownames(mirdf.sig)
      
      library(data.table)
      mirdf.sig.long <- melt(mirdf.sig,id=c("id","TSI"))
      name.idx <- strsplit(as.character(mirdf.sig.long$variable),"_")
      for (i in 1:nrow(mirdf.sig.long)) {
        mirdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        mirdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      mirdf.sig.long$condition[mirdf.sig.long$condition=="NOR"] <- "Normoxia"
      mirdf.sig.long$condition[mirdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      mirdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))->p.dat
      
      
      mir.dot.dat <- data.frame(merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id"),RNA="miRNA")
      rbind(rny.dot.dat,mir.dot.dat)
      # dcast(data = p.dat,id~condition+tissue+value_mean+sd+se) 
      ##点图，常氧下的组织特异性（全部）
      if (T) {
        library(ggplot2)
        # mirdf.sig.long%>%filter(value>0&condition=="Normoxia")%>%
        # p.dat%>%filter(value>0&condition=="Normoxia")%>%
        p.dat.log  <- merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id")%>%
                        filter(condition=="Normoxia")%>%
        # log10(p.dat.log$value_mean)
        # rbind(rny.dot.dat,mir.dot.dat)%>%
        #   filter(condition=="Normoxia")%>%
          ggplot(aes(y=id,x=tissue))+
          geom_point(aes(size=log(value_mean, 10),color=padj))+
          scale_size_continuous(range = c(-1,8),breaks = seq(.5,3,.5))+
          # facet_grid(~RNA)+
          # facet_wrap(~RNA,ncol = 1)+.
          # geom_point(fun = "mean", stat = "summary",aes(group=tissue,fill=condition,color=condition,size=log(value)),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=id,fill=condition,color=condition,size=log(value)))+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
        # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
        # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
        # facet_wrap(~tissue,nrow=1,strip.position = "top")+
        labs(y=NULL,x=NULL,title = "miRNA")+
          # scale_color_gradient(low = "blue",high = "red")+
          viridis::scale_color_viridis()+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(size=guide_legend(reverse = F,title = "log10(CPM)"),
                 # color=guide_legend(title = "p.adj")
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(),
                axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
          )->p.mir.tsi.dot
        p.mir.tsi.dot
        
        options(download.file.method = 'libcurl')
        options(url.method='libcurl')
        # BiocManager::install("lemon")
        library(lemon)
        # cowplot::plot_grid()
        
        # unique(rny.dot.dat$id)
        
        
        #####
        #拼点图 TSI
        patchwork::wrap_plots(guides = "collect",p.rny.tsi.dot, p.mir.tsi.dot)
        grid_arrange_shared_legend(, p.mir.tsi.dot, p.pir.tsi.dot,nrow = 3,ncol = 1,position='right')
        
        p1=p.rny.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(rny.dot.dat$id))/1.5,"cm"))
        p2=p.mir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(mir.dot.dat$id))/1.5,"cm"))
        p3=p.pir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(pir.dot.dat$id))/1.5,"cm"))
        p.trf.tsi.dot
        
        
        p <- grid_arrange_shared_legend(p1,p2,nrow = 2,ncol = 1,position='right')
        
        ggsave(plot = p,filename = "test.pdf",units = "cm",width = 60,height = 500,path = "C:\\Users\\colet\\Desktop\\re snc atlas",limitsize = F)
        p <- patchwork::wrap_plots(guides = "collect",p1, p2,p3,ncol = 1)
        #拼点图 TSI
        #####      
        # ggsave(filename = "TIS.dot.EXP.mir.pdf",units = "cm",width = 16,height = 55,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##点图，常氧下的组织特异性
      
      ##点图，常氧下的组织特异性（部分）
      if (T) {
        library(ggplot2)
        p.dat.log  <- merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id")%>%
          filter(condition=="Normoxia")
        
        
        p.dat.log[order(p.dat.log$padj),]%>%filter(value_mean>150&padj < .001)->p.dat.log.1
        strsplit2(p.dat.log.1$id,"-.p")[,1]->p.dat.log.1$id
        p.dat.log.1%>%filter(id!="rno-miR-132"&id!="rno-miR-3546")%>%
          ggplot(aes(y=id,x=tissue))+
          geom_point(aes(size=log(value_mean, 10),color=padj))+
          scale_size_continuous(range = c(-1,8),breaks = seq(.5,3,.5))+
        labs(y=NULL,x=NULL,title = "miRNA")+
          viridis::scale_color_viridis()+
          theme_gray()+
          guides(size=guide_legend(reverse = F,title = "log10(CPM)"),
                 # color=guide_legend(title = "p.adj")
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(),
                axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
          )->p.mir.tsi.dot.selected
        p.mir.tsi.dot.selected
        
        options(download.file.method = 'libcurl')
        options(url.method='libcurl')
        # BiocManager::install("lemon")
        library(lemon)
        # cowplot::plot_grid()
        
        # unique(rny.dot.dat$id)
        
        
        #####
        #拼点图 TSI
        patchwork::wrap_plots(guides = "collect",p.rny.tsi.dot, p.mir.tsi.dot)
        grid_arrange_shared_legend(, p.mir.tsi.dot, p.pir.tsi.dot,nrow = 3,ncol = 1,position='right')
        
        p1=p.rny.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(rny.dot.dat$id))/1.5,"cm"))
        p2=p.mir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(mir.dot.dat$id))/1.5,"cm"))
        p3=p.pir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(pir.dot.dat$id))/1.5,"cm"))
        p.trf.tsi.dot
        
        
        p <- grid_arrange_shared_legend(p1,p2,nrow = 2,ncol = 1,position='right')
        
        ggsave(plot = p,filename = "test.pdf",units = "cm",width = 60,height = 500,path = "C:\\Users\\colet\\Desktop\\re snc atlas",limitsize = F)
        p <- patchwork::wrap_plots(guides = "collect",p1, p2,p3,ncol = 1)
        #拼点图 TSI
        #####      
        # ggsave(filename = "TIS.dot.EXP.mir.pdf",units = "cm",width = 16,height = 55,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##点图，常氧下的组织特异性
      
      ##两个柱子堆积成1
      if (F) {
        library(ggplot2)
        mirdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.fillbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (F) {
        library(ggplot2)
        mirdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.stackbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
    }
    
  }##新加入新分析
  
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(mir.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = mir.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.mir <- data.frame(merge(data.frame(id=rownames(mir.df.hyp.all),mir.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  vsted.nor = assay(v)
  
  ##calculate TSI
  tsi.mir.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.mir.nor <- data.frame(tsi.mir.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.mir.nor <- tsi.mir.nor[,-1]
  colnames(tsi.mir.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.mir.nor <- tsi.mir.nor
  tsi.mir.nor <- data.frame(tsi.mir.nor=TSI(tsi.mir.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.mir.nor
  p.pca.mir.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(mir.df.hyp.all),dex = c(rep("thymus",3),
                                                              rep("kidney",3),
                                                              rep("brain",3),
                                                              rep("liver",3),
                                                              rep("lung",3),
                                                              rep("intestines",3),
                                                              rep("heart",3),
                                                              rep("spleen",3)
  ))
  
  
  t.mir.df.hyp.all <- data.frame(t(mir.df.hyp.all))
  t.mir.df.hyp.all$tissue <- df.col.2$dex
  
  
  dds<-DESeqDataSetFromMatrix(countData = mir.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  
  # ##挑选显著的id
  # ids <- rownames(res[data.frame(res)$pvalue<0.05,])
  # mir.df.hyp.scale <- data.frame(scale(mir.df.hyp.all))
  # mir.df.hyp.scale.1 <- mir.df.hyp.scale[rownames(mir.df.hyp.scale)%in%ids,]
  
  # apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  vsted.hyp = assay(v)
  
  
  # vsted=mir.df.hyp.scale.1
  
  ##calculate TSI
  tsi.mir.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.mir.hyp <- data.frame(tsi.mir.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.mir.hyp <- tsi.mir.hyp[,-1]
  colnames(tsi.mir.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.mir.hyp <- tsi.mir.hyp
  tsi.mir.hyp <- data.frame(tsi.mir.hyp=TSI(tsi.mir.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.mir.hyp
  p.pca.mir.hyp
  
  #  p.pca.mir.nor|p.pca.mir.hyp -> p.pca.mir.nor.hyp.unite
  
  ##determine tissue specific mirRNA
  tsi.mir.nor  <- data.frame(id=rownames(tsi.mir.nor),tsi.mir.nor,row.names = NULL)
  tsi.mir.nor.spc <- tsi.mir.nor[tsi.mir.nor$tsi.mir.nor>0,]
  tsi.mir.hyp  <- data.frame(id=rownames(tsi.mir.hyp),tsi.mir.hyp,row.names = NULL)
  tsi.mir.hyp.spc <- tsi.mir.hyp[tsi.mir.hyp$tsi.mir.hyp>0,]
  
  ##在均值数据框中确定reads大于10的
  # avgexp.mir.nor <- avgexp.mir.nor[rowSums(avgexp.mir.nor)>10,]
  # avgexp.mir.hyp <- avgexp.mir.hyp[rowSums(avgexp.mir.hyp)>10,]
  ##determine expression of tissue specific mirRNA
  library(data.table)
  library(ggplot2)
  ##mir对照组
  tsi.mir.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.mir.nor),avgexp.mir.nor),tsi.mir.nor.spc,by="id")
  # write.csv(tsi.mir.nor.spc.long.1,"tsi.mir.nor.spc.long.1.csv",row.names = F)
  ##取前100
  # tsi.mir.nor.long.2 <- head(tsi.mir.nor.spc.long.1[order(tsi.mir.nor.spc.long.1$tsi.mir.nor,decreasing = T),],)
  
  tsi.mir.nor.spc.long <- melt(tsi.mir.nor.spc.long.1,id=c("id","tsi.mir.nor"))
  colnames(tsi.mir.nor.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.mir.nor.spc.long <- tsi.mir.nor.spc.long[order(tsi.mir.nor.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.mir.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    scale_size_continuous(range=c(-1,6))+
    viridis::scale_color_viridis()->p.mir.nor.spc 
  p.mir.nor.spc 
  
  ##mir处理组
  tsi.mir.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.mir.hyp),avgexp.mir.hyp),tsi.mir.hyp.spc,by="id")
  
  #取前100
  # tsi.mir.hyp.spc.long.1 <- head(tsi.mir.hyp.spc.long.1[order(tsi.mir.hyp.spc.long.1$tsi.mir.hyp,decreasing = T),],)
  tsi.mir.hyp.spc.long.1
  
  
  tsi.mir.hyp.spc.long <- melt(tsi.mir.hyp.spc.long.1,id=c("id","tsi.mir.hyp"))
  colnames(tsi.mir.hyp.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.mir.hyp.spc.long <- tsi.mir.hyp.spc.long[order(tsi.mir.hyp.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.mir.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    scale_size_continuous(range=c(-1,6))+
    viridis::scale_color_viridis()->p.mir.hyp.spc 
  p.mir.hyp.spc 
  
  tsi.mir.nor.spc.long <- data.frame(tsi.mir.nor.spc.long,condition="Control")
  tsi.mir.hyp.spc.long <- data.frame(tsi.mir.hyp.spc.long,condition="Hypoxia")
  tsi.mir.combine.spc.long <- rbind(tsi.mir.nor.spc.long,tsi.mir.hyp.spc.long)
  #####
  ###重新找出不显著的画灰点(2022.09.25后新加)
  tsi.mir.nor.nonspc <- data.frame(tsi.mir.nor[tsi.mir.nor$id%in%tsi.mir.combine.spc.long$ID,], condition = "Control")
  tsi.mir.nor.nonspc.long.1 <- merge(data.frame(id=tsi.mir.nor$id[tsi.mir.nor$id%in%tsi.mir.combine.spc.long$ID],
                                               avgexp.mir.nor[rownames(avgexp.mir.nor)%in%tsi.mir.combine.spc.long$ID,]),
                                    tsi.mir.nor.nonspc,by="id")
  tsi.mir.nor.nonspc.long <- melt(tsi.mir.nor.nonspc.long.1,id=c("id","condition","tsi.mir.nor"))
  
  tsi.mir.hyp.nonspc <- data.frame(tsi.mir.hyp[tsi.mir.hyp$id%in%tsi.mir.combine.spc.long$ID,], condition = "Hypoxia")
  tsi.mir.hyp.nonspc.long.1 <- merge(data.frame(id=tsi.mir.hyp$id[tsi.mir.hyp$id%in%tsi.mir.combine.spc.long$ID],
                                               avgexp.mir.hyp[rownames(avgexp.mir.hyp)%in%tsi.mir.combine.spc.long$ID,]),
                                    tsi.mir.hyp.nonspc,by="id")
  tsi.mir.hyp.nonspc.long <- melt(tsi.mir.hyp.nonspc.long.1,id=c("id","condition","tsi.mir.hyp"))
  
  merge(tsi.mir.hyp.nonspc.long,data.frame(ID=tsi.mir.combine.spc.long$ID,Gap=tsi.mir.combine.spc.long$Gap),by="ID")
  
  ##准备表达数据
  id.list <- data.frame(id=unique(tsi.mir.nor.nonspc.long.1$id,tsi.mir.hyp.nonspc.long.1$id))
  
  tsi.mir.nor.nonspc.long.3=tsi.mir.nor.nonspc.long.1
  tsi.mir.hyp.nonspc.long.3=tsi.mir.hyp.nonspc.long.1
  
  colnames(tsi.mir.nor.nonspc.long.3)[-1] <- paste0(colnames(tsi.mir.nor.nonspc.long.3)[-1],"_Normoxia")
  colnames(tsi.mir.hyp.nonspc.long.3)[-1] <- paste0(colnames(tsi.mir.hyp.nonspc.long.3)[-1],"_Hypoxia")
  
  #计算基因的最大减最小
  expr.gap.nor <- rowMax(as.matrix(tsi.mir.nor.nonspc.long.3[,-c(1,10,11)]))- rowMin(as.matrix(tsi.mir.nor.nonspc.long.3[,-c(1,10,11)]))
  expr.sd.nor <- rowSds(as.matrix(tsi.mir.nor.nonspc.long.3[,-c(1,10,11)]))
  tsi.mir.nor.nonspc.long.3 <- cbind(tsi.mir.nor.nonspc.long.3,expr.gap.nor,expr.sd.nor)
  
  expr.gap.hyp <- rowMax(as.matrix(tsi.mir.hyp.nonspc.long.3[,-c(1,10,11)]))- rowMin(as.matrix(tsi.mir.hyp.nonspc.long.3[,-c(1,10,11)]))
  expr.sd.hyp <- rowSds(as.matrix(tsi.mir.hyp.nonspc.long.3[,-c(1,10,11)]))
  tsi.mir.hyp.nonspc.long.3 <- cbind(tsi.mir.hyp.nonspc.long.3,expr.gap.hyp,expr.sd.hyp)
  
  
  mir.exprset <- left_join(left_join(id.list,tsi.mir.nor.nonspc.long.3,by="id"),
                           tsi.mir.hyp.nonspc.long.3)
  write.csv(mir.exprset,file = "mir.ExprSet.csv",quote = FALSE, row.names = FALSE)
  
  
  mir.avg.expr.dat <- merge(tsi.mir.nor.nonspc.long.1,tsi.mir.hyp.nonspc.long.1,by="id")
  colnames(mir.avg.expr.dat)[1:10] <- gsub(pattern = ".x",replacement = "_Normoxia",colnames(mir.avg.expr.dat)[1:10])
  colnames(mir.avg.expr.dat)[12:21] <- gsub(pattern = ".y",replacement = "_Hypoxia",colnames(mir.avg.expr.dat)[12:21])
  
  colnames(tsi.mir.nor.nonspc.long) <- c("ID","condition","TSI","Tissue","VST")
  colnames(tsi.mir.hyp.nonspc.long) <- c("ID","condition","TSI","Tissue","VST")
  
  tem <- merge(tsi.mir.nor.nonspc.long,tsi.mir.hyp.nonspc.long,by="ID")
  abs(tem$TSI.x-tem$TSI.y)[order(abs(tem$TSI.x-tem$TSI.y),decreasing = T)]
  
  tsi.mir.uni.nonspc.long <- rbind(tsi.mir.nor.nonspc.long, tsi.mir.hyp.nonspc.long)
  ##### 
  ggplot()+
    geom_point(data = tsi.mir.uni.nonspc.long, aes(x=Tissue,y=ID,size=VST),color="grey50",fill="black",alpha=.5)+
    geom_point(data=tsi.mir.combine.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    # geom_point(data=tsi.ys.unite.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    facet_wrap(~condition)+
    scale_size_continuous(range=c(-1,6))+
    labs(x="",y="")+
    viridis::scale_color_viridis()+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))->p.mir.combine.spc 
  p.mir.combine.spc
  
  setwd("C:/Users/colet/Desktop/")
  # ggsave(filename = "p.mir.combine.spc.over.all.pdf",p.mir.combine.spc,units = "cm",width = 16,height = 40,limitsize = F)
  
  tsi.mir.nor.spc.long.3  <- tsi.mir.nor.spc.long.1
   colnames(tsi.mir.nor.spc.long.3)[-1] <- paste0(colnames(tsi.mir.nor.spc.long.3)[-1],"_control")
   tsi.mir.hyp.spc.long.3  <- tsi.mir.hyp.spc.long.1
   colnames(tsi.mir.hyp.spc.long.3)[-1] <- paste0(colnames(tsi.mir.hyp.spc.long.3)[-1],"_hypoxia")
   
  id  <- data.frame(id=unique(c(tsi.mir.nor.spc.long.3$id,tsi.mir.hyp.spc.long.3$id)))
  
  nunite.tsi.mir <- data.frame(left_join(left_join(id,tsi.mir.nor.spc.long.3),tsi.mir.hyp.spc.long.3))
  for (i in 1:ncol(nunite.tsi.mir)){
    nunite.tsi.mir[,i][is.na(nunite.tsi.mir[,i])] <- 0
  }
  nunite.tsi.mir$abs.gap <- abs(nunite.tsi.mir$tsi.mir.nor_control-nunite.tsi.mir$tsi.mir.hyp_hypoxia)
  nunite.tsi.mir$tsi.gap <- nunite.tsi.mir$tsi.mir.nor_control-nunite.tsi.mir$tsi.mir.hyp_hypoxia
  nunite.tsi.mir[order(nunite.tsi.mir$abs.gap,decreasing = T),]
  # write.csv(nunite.tsi.mir,file = "nunite.tsi.mir.csv")
  
  
  ####################3p/5p特异性
  ##5p/3p特异性
  i=1
  ##获取手臂名称和miR名称
  #1.对照组
  mir.df.vst.nor <- data.frame(vsted.nor)
  for (i in 1:nrow(mir.df.vst.nor)) {
    mir.df.vst.nor$arm[i] <- ifelse(
      strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]])]=="3p","3p",
      ifelse(strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]])]=="5p","5p","no"
      )
    )
    mir.df.vst.nor$mirname[i] <- paste(
      strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.nor)[3],"-")[[1]])-1)][1],
      strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.nor)[3],"-")[[1]])-1)][2],
      strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.nor)[3],"-")[[1]]))][3],
      sep = "-")
  }
  
  mir.df.vst.nor <- na.omit(mir.df.vst.nor)
  
  
  mir.df.vst.nor.3p <- mir.df.vst.nor[mir.df.vst.nor$arm=="3p",]
  mir.df.vst.nor.5p <- mir.df.vst.nor[mir.df.vst.nor$arm=="5p",]
  
  
  mir.df.vst.nor[,c(1:24)] <- log2(mir.df.vst.nor[,c(1:24)])
  
  library(data.table)
  mir.df.vst.nor.long <- melt(mir.df.vst.nor,id=c("mirname","arm"))
  colnames(mir.df.vst.nor.long) <- c("miRNA","Arm","Sample","Reads") 
  
  ##
  for (i in 1:nrow(mir.df.vst.nor.long)) {
    mir.df.vst.nor.long$Tissue[i] <- strsplit(as.character(mir.df.vst.nor.long$Sample)[i],"_")[[1]][1]
    mir.df.vst.nor.long$Condition[i] <- strsplit(as.character(mir.df.vst.nor.long$Sample)[i],"_")[[1]][2]
  }
  
  # BiocManager::install("introdataviz")
  # BiocManager::install("CBNplot")
  library(tidyverse)
  # library(introdataviz)
  library(ggunchained)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(ggprism)
  
  mir.df.vst.nor.long.rmarmless <- mir.df.vst.nor.long[mir.df.vst.nor.long$Arm!="no",]
  ggplot(data = mir.df.vst.nor.long.rmarmless%>%filter(Condition=="NOR"),aes(x=Tissue,y=Reads,fill=Arm))+
    geom_split_violin()+
    facet_wrap(~miRNA,ncol = 4)+
    stat_summary(fun.data = "mean_sd", position = position_dodge(0.15), geom = "errorbar", width=.1)+
    stat_summary(fun = "mean", geom = "point", position = position_dodge(0.15), show.legend = F)+
    stat_compare_means(aes(group = mir.df.vst.nor.long.rmarmless$Arm),label = "p.signif", method = "t.test")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1 ))-> p.mir.arm.nor
  #  ggsave("violin.arm.nor.pdf",units = "cm", width = 100,height = 600,limitsize = F)
  
  ##2、处理组
  mir.df.vst.hyp <- data.frame(vsted.hyp)
  for (i in 1:nrow(mir.df.vst.hyp)) {
    mir.df.vst.hyp$arm[i] <- ifelse(
      strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]])]=="3p","3p",
      ifelse(strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]])]=="5p","5p","no"
      )
    )
    mir.df.vst.hyp$mirname[i] <- paste(
      strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.hyp)[3],"-")[[1]])-1)][1],
      strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.hyp)[3],"-")[[1]])-1)][2],
      strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.hyp)[3],"-")[[1]]))][3],
      sep = "-")
  }
  
  mir.df.vst.hyp <- na.omit(mir.df.vst.hyp)
  
  
  mir.df.vst.hyp.3p <- mir.df.vst.hyp[mir.df.vst.hyp$arm=="3p",]
  mir.df.vst.hyp.5p <- mir.df.vst.hyp[mir.df.vst.hyp$arm=="5p",]
  
  
  mir.df.vst.hyp[,c(1:24)] <- log2(mir.df.vst.hyp[,c(1:24)])
  
  library(data.table)
  mir.df.vst.hyp.long <- melt(mir.df.vst.hyp,id=c("mirname","arm"))
  colnames(mir.df.vst.hyp.long) <- c("miRNA","Arm","Sample","Reads") 
  
  ##
  for (i in 1:nrow(mir.df.vst.hyp.long)) {
    mir.df.vst.hyp.long$Tissue[i] <- strsplit(as.character(mir.df.vst.hyp.long$Sample)[i],"_")[[1]][1]
    mir.df.vst.hyp.long$Condition[i] <- strsplit(as.character(mir.df.vst.hyp.long$Sample)[i],"_")[[1]][2]
  }
  
  # BiocManager::install("introdataviz")
  # BiocManager::install("CBNplot")
  library(tidyverse)
  # library(introdataviz)
  library(ggunchained)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(ggprism)
  
  mir.df.vst.hyp.long.rmarmless <- mir.df.vst.hyp.long[mir.df.vst.hyp.long$Arm!="no",]
  ggplot(data = mir.df.vst.hyp.long.rmarmless,aes(x=Tissue,y=Reads,fill=Arm))+
    geom_split_violin()+
    facet_wrap(~miRNA,ncol = 4)+
    stat_summary(fun.data = "mean_sd", position = position_dodge(0.15), geom = "errorbar", width=.1)+
    stat_summary(fun = "mean", geom = "point", position = position_dodge(0.15), show.legend = F)+
    stat_compare_means(aes(group = mir.df.vst.hyp.long.rmarmless$Arm),label = "p.signif", method = "t.test")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1 ))-> p.mir.arm.hyp
  #ggsave("violin.arm.hyp.pdf",units = "cm", width = 100,height = 600,limitsize = F)
  
  ##3、缺氧后是否出现3p、5p转换
  vsted <- right_join(data.frame(id=rownames(vsted.nor),vsted.nor),data.frame(id=rownames(vsted.hyp),vsted.hyp))
  vsted <- data.frame(vsted,row.names = 1)
  for (i in 1:ncol(vsted)) {
    vsted[,i][is.na(vsted[,i])]=0
  }
  
  mir.df.vst.all <- data.frame(vsted)
  for (i in 1:nrow(mir.df.vst.all)) {
    mir.df.vst.all$arm[i] <- ifelse(
      strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.all)[i],"-")[[1]])]=="3p","3p",
      ifelse(strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.all)[i],"-")[[1]])]=="5p","5p","no"
      )
    )
    mir.df.vst.all$mirname[i] <- paste(
      strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.all)[3],"-")[[1]])-1)][1],
      strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.all)[3],"-")[[1]])-1)][2],
      strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.all)[3],"-")[[1]]))][3],
      sep = "-")
  }
  
  mir.df.vst.all <- na.omit(mir.df.vst.all)
  
  
  mir.df.vst.all.3p <- mir.df.vst.all[mir.df.vst.all$arm=="3p",]
  mir.df.vst.all.5p <- mir.df.vst.all[mir.df.vst.all$arm=="5p",]
  
  
  mir.df.vst.all[,c(1:48)] <- log2(mir.df.vst.all[,c(1:48)])
  
  library(data.table)
  mir.df.vst.all.long <- melt(mir.df.vst.all,id=c("mirname","arm"))
  colnames(mir.df.vst.all.long) <- c("miRNA","Arm","Sample","Reads") 
  
  ##
  for (i in 1:nrow(mir.df.vst.all.long)) {
    mir.df.vst.all.long$Tissue[i] <- strsplit(as.character(mir.df.vst.all.long$Sample)[i],"_")[[1]][1]
    mir.df.vst.all.long$Condition[i] <- strsplit(as.character(mir.df.vst.all.long$Sample)[i],"_")[[1]][2]
  }
  
  # BiocManager::install("introdataviz")
  # BiocManager::install("CBNplot")
  library(tidyverse)
  # library(introdataviz)
  library(ggunchained)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(ggprism)
##未筛选  
  mir.df.vst.all.long.rmarmless <- mir.df.vst.all.long[mir.df.vst.all.long$Arm!="no",]
  ggplot(data = mir.df.vst.all.long.rmarmless,aes(x=Condition,y=Reads,fill=Arm))+
    geom_split_violin()+
    facet_wrap(Tissue~miRNA,ncol = 4)+
    stat_summary(fun.data = "mean_sd", position = position_dodge(0.15), geom = "errorbar", width=.1)+
    stat_summary(fun = "mean", geom = "point", position = position_dodge(0.15), show.legend = F)+
    stat_compare_means(aes(group = mir.df.vst.all.long.rmarmless$Arm),label = "p.signif", method = "t.test")+
    #  theme_prism()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1 ))-> p.mir.arm.all

##读取人工筛选后的列表 
  setwd("C:/Users/colet/Desktop/re snc atlas/")
  select.arm <- read.table("violin.arm.all.txt")
  select.arm.1 <- data.frame(miRNA=select.arm[seq(2,nrow(select.arm),2),],Tissue=select.arm[seq(1,nrow(select.arm),2),])
##匹配筛选后列表到未筛选列表  
  select.arm.matched <- mir.df.vst.all.long.rmarmless%>%filter(mir.df.vst.all.long.rmarmless$miRNA%in%select.arm.1$miRNA)
  select.arm.matched$Reads[select.arm.matched$Reads=="-Inf"] <- 0
  any(select.arm.matched$Reads=="-Inf")
  select.arm.matched$Reads[select.arm.matched$Reads=="-Inf"] <- 0
  select.arm.matched$Condition <- ifelse(select.arm.matched$Condition=="NOR","Control","Hypoxia")
  
  select.arm.matched$group <- paste(select.arm.matched$Tissue,select.arm.matched$Condition,sep = "-")
  shadow <- data.frame(x=c("heart-Control","heart-Hypoxia","kidney-Control","kidney-Hypoxia",
                           "lung-Control","lung-Hypoxia","thymus-Control","thymus-Hypoxia"),
                       y=4)
  
  
  ggplot(data = select.arm.matched,aes(x=group,y=Reads,fill=Arm))+
    #geom_split_violin()+
    # introdataviz::geom_split_violin(alpha = .4)+
     geom_split_violin(trim= F,scale = "area",color="white")+
     # geom_col(data = shadow,aes(x=x,y=y,color="grey50",fill="black"),alpha=.3)+
    geom_vline(xintercept=c(2.5,4.5,6.5,8.5,10.5,12.5,14.5), linetype=5,size=.5)+
    ylim(2,4)+
    facet_wrap(~miRNA,ncol=1)+
    scale_fill_manual(values = c("#2fa1dd", "#f87669"))+
    # stat_summary(fun.data = "mean_sd", position = position_dodge(0.15), geom = "errorbar", width=.1)+
    # stat_summary(fun = "mean", geom = "point", position = position_dodge(0.15), show.legend = F)+
    stat_compare_means(aes(group = select.arm.matched$Arm),label = "p.signif", method = "t.test")+
    # theme_bw()+
     theme_prism()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1 ))->p.select.arm
  # p.select.arm
  # ggsave(plot = p.select.arm,"p.select.arm.pdf",units = "cm", width = 16,height = 450,limitsize = F)
  
  ###循环出图
  mir.list <- unique(select.arm.matched$miRNA)
  p.arm <- vector("list", length = length(mir.list))
  p.arm <- list()
  i=1
  
  # install.packages("see")
  library(see)
  
  library(ggplot2)
  library(tidyverse)
  library(gghalves)
  for (i in 1:length(mir.list)) {
    select.arm.matched%>%filter(miRNA==as.character(mir.list[i]))%>%
    ggplot(aes(x=group,y=Reads,fill=Arm))+
      # geom_split_violin(trim= T,scale = "area",color="grey50")+
       geom_split_violin(trim= T,scale = "area")+
       geom_boxplot(width=0.75)+
      # geom_point(aes(color=Arm),size=.75)+
      # geom_violindot(size_dots = 1) +
      # geom_violindot()+
      # geom_split_violin(trim= F,scale = "area")+
      geom_vline(xintercept=c(2.5,4.5,6.5,8.5,10.5,12.5,14.5), linetype=5,size=.5,color="grey50")+
      ylim(1,4)+
      labs(title = as.character(mir.list[i]),x="")+
      scale_fill_manual(values = c("#2fa1dd", "#f87669"))+
      scale_color_manual(values = c("#2fa1dd", "#f87669"))+
      # scale_color_manual(values = c("black"))+
      stat_compare_means(aes(group = Arm),label = "p.signif", method = "t.test")+
      theme_light()+
      # theme_prism()+
      # theme(axis.title=element_text(size=12))+
      # theme(axis.text.x = element_text(angle = -90, hjust = 0))->p.arm[[i]]
       theme(axis.text.x = element_text(angle = 45, hjust = 1 ))->p.arm[[i]]
  }

  # install.packages("patchwork")
  library(patchwork)
  patchwork::wrap_plots(p.arm,ncol = 1,byrow = TRUE, plot_layout(guides = "collect"))->p.done.arm
  # p.done.arm
  # ggsave(plot = p.done.arm,"p.select.arm.pdf",units = "cm", width = 16,height = 900,limitsize = F)
  
  ##再次选择后的
  setwd("C:/Users/colet/Desktop/re snc atlas/")
  select.arm <- read.table("selectarm3.txt")
  select.arm.1 <- select.arm 
  
  ###挑选50个特别的miRNA 手臂特异性
  mir.list <- unique(select.arm.matched$miRNA)
  p.arm <- list()
  for (i in 1:length(mir.list)) {
    select.arm.matched%>%filter(miRNA==as.character(mir.list[i])&Reads>0)%>%
      ggplot(aes(x=group,y=Reads,fill=Arm))+
       # geom_jitter(aes(color=Arm),size=1)+
      # geom_split_violin(trim= T,scale = "area",color="grey50")+
      stat_summary(aes(col = Reads), fun.data = 'mean_sd', geom = "errorbar", colour = "black", width = .4,
                   position = position_dodge( .7)) +
      # geom_split_violin(trim= T,scale = "area",width=1,position = position_dodge(width = 1))+
      # geom_split_violin(trim = F,position = position_dodge(width = 1), scale = 'width',alpha=0.8,width=.2)+
       geom_boxplot(width=0.75)+
      # geom_point(aes(color=Arm),size=.75)+
      # geom_violin(trim = T,position = position_dodge(width = 1), scale = 'width',alpha=0.8,width=1) +
      # geom_violindot()+
      # geom_split_violin(trim= F,scale = "area")+
      # stat_summary(fun.data = "mean_sd", position = position_dodge(0.15), geom = "errorbar", width=.1)+
      geom_vline(xintercept=c(2.5,4.5,6.5,8.5,10.5,12.5,14.5), linetype=5,size=.5,color="grey50")+
      # ylim(0,5)+
      labs(title = as.character(mir.list[i]),x="")+
      scale_fill_manual(values = c("#87CFEA", "#FFC0CB"))+
      scale_color_manual(values = c("#87CFEA", "#FFC0CB"))+
      # scale_color_manual(values = c("#2fa1dd", "#f87669"))+
      # scale_color_manual(values = c("black"))+
      # stat_compare_means(aes(group = Arm),label = "p.signif", method = "t.test")+
      theme_light()+
      # theme_prism()+
      theme(axis.title=element_text(size=8))+
      # theme(axis.text.x = element_text(angle = -90, hjust = 0))->p.arm[[i]]
      theme(axis.text.x = element_text(angle = 45, hjust = 1 ))->p.arm[[i]]
  }
  # p.arm[[5]]
  # install.packages("patchwork")
  library(patchwork)
  patchwork::wrap_plots(p.arm,ncol = 3,byrow = TRUE, plot_layout(guides = "collect"))->p.done.arm.2
  # p.done.arm
  # ggsave(plot = p.done.arm.2,"p.select.arm.select.pdf",units = "cm", width = 48,height = 250,limitsize = F)
  ##miRswitch
#####################################
##最终挑选后的具有明显变化的手臂
  final.select  <- read.csv("selectarm.final.CSV")
  mir.list <- final.select%>%filter(!(Type%in%c("3p","5p")))
  mir.list  <- mir.list$miRNA
  mir.list  <- mir.list[mir.list!="rno-miR-3585"]
   p.arm <- list()
   # for (i in 1:length(mir.list)) {
   #   select.arm.matched%>%filter(miRNA==as.character(mir.list[i])&Reads>0)%>%
       select.arm.matched%>%filter(miRNA%in%as.character(mir.list))%>%
       ggplot(aes(x=group,y=Reads,fill=Arm))+
       # geom_jitter(aes(color=Arm),size=.75)+
       stat_summary(aes(col = Reads), fun.data = 'mean_sd', geom = "errorbar", colour = "black", width = .4,
                    position = position_dodge( .7)) +
       geom_boxplot(width=0.75)+
       facet_wrap(~miRNA,ncol=2)+
       geom_vline(xintercept=c(2.5,4.5,6.5,8.5,10.5,12.5,14.5), linetype=5,size=.5,color="grey50")+
        ylim(2.2,3)+
       labs(x="",y="")+
       scale_fill_manual(values = c("#87CFEA", "#FFC0CB"))+
       scale_color_manual(values = c("#87CFEA", "#FFC0CB"))+
       theme_light()+
       # stat_compare_means(aes(group = Arm),label = "p.signif", method = "t.test")+
       # theme(axis.title=element_text(size=3))+
       theme(axis.text.x = element_text(angle = 45, hjust = 1 ))->p.arm
   # }
   # p.arm[[1]]
   # library(patchwork)
   # patchwork::wrap_plots(p.arm,ncol = 1,byrow = TRUE, plot_layout(guides = "collect"))->p.done.arm.2
   # ggsave(plot = p.arm,"p.select.arm.select.all.final.pdf",units = "cm", width =28,height = 28,limitsize = F)
   ##画3p5p丰度雨云图    
   library(gghalves)
   
   mir.df.vst.all.long.rmarmless$Reads[mir.df.vst.all.long.rmarmless$Reads=="-Inf"]=0
   mir.df.vst.all.long.rmarmless <- mir.df.vst.all.long.rmarmless[mir.df.vst.all.long.rmarmless$Reads!=0,]
   fun.0.1 = function(x){(x-min(mir.df.vst.all.long.rmarmless$Reads))/(max(mir.df.vst.all.long.rmarmless$Reads)-min(mir.df.vst.all.long.rmarmless$Reads))}
   mir.df.vst.all.long.rmarmless$Reads <- fun.0.1(mir.df.vst.all.long.rmarmless$Reads)
   
   ggplot(mir.df.vst.all.long.rmarmless, aes(x = Reads, y = Tissue)) +
     theme_bw() +
     facet_grid(~Condition)+
     geom_point(aes(color = Arm), 
                position = position_jitterdodge(jitter.width = .2, jitter.height = .1,dodge.width = 0.6), 
                size = 3, shape = 20, alpha = 0.6) +
     geom_boxplot( aes(fill = Arm),position = position_nudge(x = -0.2), 
                   width = .1, alpha = .3, outlier.shape = NA) +
     xlim(0,1)+
     labs(x="Abundance")+
     geom_half_violin(aes(fill = Arm, color = Arm) ,position = position_nudge(x = -.25), 
                      adjust = 1, trim = T, alpha = .5, colour = NA, side = "l")+
     theme(axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           axis.text.x = element_text(size = 12,angle = 45,hjust = 1),
           axis.text.y = element_text(size = 12),
           strip.text.x = element_text(size = 12),
           legend.position = "right",
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()
     ) +
     ## 通过 scale 系列函数自定义标度
     scale_fill_brewer(palette = "Dark2")+
     scale_color_brewer(palette = "Dark2")->p.abundance.3p5p
   p.abundance.3p5p
   # setwd("C:/Users/colet/Desktop/")
   # ggsave(p.abundance.3p5p,filename = "p.abundance.3p5p.pdf",units = "cm",width = 16,height = 16)

#####################################
   
   # mir.in.b <- data.frame(vsted.nor)
   for (i in 1:nrow(mir.in.b)) {
     mir.in.b$arm[i] <- ifelse(
       strsplit(rownames(mir.in.b)[i],"-")[[1]][length(strsplit(rownames(mir.in.b)[3],"-")[[1]])]=="3p","3p",
       ifelse(strsplit(rownames(mir.in.b)[i],"-")[[1]][length(strsplit(rownames(mir.in.b)[3],"-")[[1]])]=="5p","5p","no"
       )
     )
     mir.in.b$mirname[i] <- paste(
       strsplit(rownames(mir.in.b)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.in.b)[3],"-")[[1]])-1)][1],
       strsplit(rownames(mir.in.b)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.in.b)[3],"-")[[1]])-1)][2],
       strsplit(rownames(mir.in.b)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.in.b)[3],"-")[[1]])-1)][3],
       sep = "-")
     
     
     
     # Ubiquitous ncRNA Transcripts
    mir.df  <- data.frame(all.mir.df[rowSums(all.mir.df)>300&rowSds(as.matrix(all.mir.df))<7.940,])
    mir.df.1 <- data.frame(id=rownames(mir.df),mir.df,row.names = NULL)
    long.mir.df.nor.all.1 <- melt(mir.df.1,id="id")
    
    for (i in 1:nrow(long.mir.df.nor.all.1)) {
      long.mir.df.nor.all.1$tissue[i] <- strsplit(as.character(long.mir.df.nor.all.1$variable[i]),"_")[[1]][1]
      long.mir.df.nor.all.1$condition[i] <- strsplit(as.character(long.mir.df.nor.all.1$variable[i]),"_")[[1]][2]
    }
    long.mir.df.nor.all.1$condition <- ifelse(long.mir.df.nor.all.1$condition=="HYP","Hypoxia","Control")
    

    ggplot(data = long.mir.df.nor.all.1,aes(x=tissue,y=value,color=condition))+
      geom_boxplot(width=0.75)+
      facet_wrap(~id,ncol=2)+
      # geom_vline(xintercept=c(2.5,4.5,6.5,8.5,10.5,12.5,14.5), linetype=5,size=.5,color="grey50")+
      # ylim(2.2,3)+
      labs(x="",y="Normalized Counts")+
      scale_color_aaas()+
      # scale_fill_manual(values = c("#87CFEA", "#FFC0CB"))+
      # scale_color_manual(values = c("#87CFEA", "#FFC0CB"))+
      theme_light()+
      # stat_compare_means(aes(group = Arm),label = "p.signif", method = "t.test")+
      # theme(axis.title=element_text(size=3))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1 ))->p.ubiquitous.mir
    p.ubiquitous.mir
      # ggsave(p.ubiquitous.mir,filename = "p.ubiquitous.mir.pdf",units = "cm",width = 16, height = 12)    
    
    
   }
   
   
   
  ##拼图和导出
  # library(cowplot)
  # p.pca.mir.nor.prism <- p.pca.mir.nor+theme_prism()+guides(x = "prism_offset_minor", y = "prism_offset_minor")#+theme(legend.position = "none")
  # p.pca.mir.hyp.prism <- p.pca.mir.hyp+theme_prism()+guides(x = "prism_offset_minor", y = "prism_offset_minor")#+theme(legend.position = "none")
  # ggsave(plot = p.pca.mir.nor.prism,"p.pca.mir.nor.prism.pdf",units = "cm", width = 16,height = 16,limitsize = F)
  # ggsave(plot = p.mir.nor.spc+theme(axis.text.x = element_text(angle=45, hjust = 1)),"p.mir.nor.spc.pdf",units = "cm", width = 12,height = 125,limitsize = F)
  # ggsave(plot = p.mir.hyp.spc+theme(axis.text.x = element_text(angle=45, hjust = 1)),"p.mir.hyp.spc.pdf",units = "cm", width = 12,height = 125,limitsize = F)
  # ggsave(plot = p.pca.mir.hyp,"p.pca.mir.hyp.prism.pdf",units = "cm", width = 16,height = 16,limitsize = F)
  # plot_grid(p.pca.mir.nor.prism,p.pca.mir.hyp.prism)
  # ggsave(plot = p.mir.arm.nor,"violin.arm.nor.pdf",units = "cm", width = 75,height = 600,limitsize = F)
  # ggsave(plot = p.mir.arm.hyp,"violin.arm.hyp.pdf",units = "cm", width = 75,height = 600,limitsize = F)
}##miRNA


if (T) {
  setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  snr.in <- read.table("snR.quantified",header = T)
  
  snr.in.b <- data.frame(snr.in[,7:ncol(snr.in)],row.names = snr.in$Geneid)
  
  i=1
  for (i in 1:ncol(snr.in.b)) {
    colnames(snr.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(snr.in.b),"_")[[i]][1],strsplit(colnames(snr.in.b),"_")[[i]][2],sep = "_")]
  }
  
  snr.df <- snr.in.b
  ##snrNA处理成对照与处理组分开的表格
  snr.df.all <- data.frame(rep(0,nrow(snr.df)))
  snr.df.nor.all = data.frame(rep(0,nrow(snr.df)))
  snr.df.hyp.all = data.frame(rep(0,nrow(snr.df)))
  for (i in 1:8) {
    snr.df.nor.all <- cbind(snr.df.nor.all,snr.df[,c(i+8,i,i+16)])
  }
  
  for (i in 25:32) {
    snr.df.hyp.all <- cbind(snr.df.hyp.all,snr.df[,c(i+8,i,i+16)])
  }
  snr.df.nor.all <- snr.df.nor.all[,-1]
  snr.df.hyp.all <- snr.df.hyp.all[,-1]
  
  ##snr差异分析（对照对对照、处理对处理；组织别）
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(snr.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = snr.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.snr <- data.frame(merge(data.frame(id=rownames(snr.df.hyp.all),snr.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.snr.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.snr.nor <- data.frame(tsi.snr.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.snr.nor <- tsi.snr.nor[,-1]
  colnames(tsi.snr.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.snr.nor <- tsi.snr.nor
  tsi.snr.nor <- data.frame(tsi.snr.nor=TSI(tsi.snr.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.snr.nor
  p.pca.snr.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(snr.df.hyp.all),dex = c(rep("thymus",3),
                                                              rep("kidney",3),
                                                              rep("brain",3),
                                                              rep("liver",3),
                                                              rep("lung",3),
                                                              rep("intestines",3),
                                                              rep("heart",3),
                                                              rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = snr.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.snr.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.snr.hyp <- data.frame(tsi.snr.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.snr.hyp <- tsi.snr.hyp[,-1]
  colnames(tsi.snr.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.snr.hyp <- tsi.snr.hyp
  tsi.snr.hyp <- data.frame(tsi.snr.hyp=TSI(tsi.snr.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.snr.hyp
  p.pca.snr.hyp
  
  # p.pca.snr.nor|p.pca.snr.hyp -> p.pca.snr.nor.hyp.unite
  
  ##determine tissue specific snrRNA
  tsi.snr.nor  <- data.frame(id=rownames(tsi.snr.nor),tsi.snr.nor,row.names = NULL)
  tsi.snr.nor.spc <- tsi.snr.nor[tsi.snr.nor$tsi.snr.nor>0.85,]
  tsi.snr.hyp  <- data.frame(id=rownames(tsi.snr.hyp),tsi.snr.hyp,row.names = NULL)
  tsi.snr.hyp.spc <- tsi.snr.hyp[tsi.snr.hyp$tsi.snr.hyp>0.85,]
  
  ##determine expression of tissue specific snrRNA
  library(data.table)
  library(ggplot2)
  ##snr对照组
  tsi.snr.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.snr.nor),avgexp.snr.nor),tsi.snr.nor.spc,by="id")
  tsi.snr.nor.spc.long <- melt(tsi.snr.nor.spc.long.1,id=c("id","tsi.snr.nor"))
  colnames(tsi.snr.nor.spc.long) <- c("ID","TSI","Tissue","Exp")
  tsi.snr.nor.spc.long <- tsi.snr.nor.spc.long[order(tsi.snr.nor.spc.long$TSI,decreasing = T),]
  
  # tsi.snr.nor.spc.long$ID = "U6-201"
  
  ggplot()+
    geom_point(data=tsi.snr.nor.spc.long,aes(x=Tissue,y=ID,size=Exp,color=TSI))+
    viridis::scale_color_viridis()->p.snr.nor.spc 
  p.snr.nor.spc
  
  ##snr处理组
  tsi.snr.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.snr.hyp),avgexp.snr.hyp),tsi.snr.hyp.spc,by="id")
  tsi.snr.hyp.spc.long <- melt(tsi.snr.hyp.spc.long.1,id=c("id","tsi.snr.hyp"))
  colnames(tsi.snr.hyp.spc.long) <- c("ID","TSI","Tissue","Exp")
  tsi.snr.hyp.spc.long <- tsi.snr.hyp.spc.long[order(tsi.snr.hyp.spc.long$TSI,decreasing = T),]
  
  # tsi.snr.hyp.spc.long$ID = "U6-201"
  
  ggplot()+
    geom_point(data=tsi.snr.hyp.spc.long,aes(x=Tissue,y=ID,size=Exp,color=TSI))+
    viridis::scale_color_viridis()->p.snr.hyp.spc 
  p.snr.hyp.spc 
  
  rbind(data.frame(tsi.snr.nor.spc.long,condition="Control",id="U6-201",VST=tsi.snr.nor.spc.long$Exp),data.frame(tsi.snr.hyp.spc.long,condition="Hypoxia",id="U6-201",VST=0))%>%
  ggplot()+
    geom_point(aes(x=Tissue,y=id,size=VST,color=TSI))+
    labs(x="",y="")+
    viridis::scale_color_viridis()+
    scale_size_continuous(range=c(-1,6))+
    facet_wrap(~condition)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))->p.snr.spc.combine
  # ggsave(p.snr.nor.spc,filename = "p.snr.nor.spc.pdf",units = "cm",width = 16,height = 5)
  # ggsave(p.snr.hyp.spc,filename = "p.snr.hyp.spc.pdf",units = "cm",width = 16,height = 5)
  # ggsave(p.snr.spc.combine,filename = "p.snr.spc.combine.pdf",units = "cm",width = 16,height = 5)
  
}##snRNA

if (T) {
  setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  sno.in <- read.table("snoR.quantified",header = T)
  
  sno.in.b <- data.frame(sno.in[,7:ncol(sno.in)],row.names = sno.in$Geneid)
  
  i=1
  for (i in 1:ncol(sno.in.b)) {
    colnames(sno.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(sno.in.b),"_")[[i]][1],strsplit(colnames(sno.in.b),"_")[[i]][2],sep = "_")]
  }
  
  sno.df <- sno.in.b
  ##snoNA处理成对照与处理组分开的表格
  sno.df.all <- data.frame(rep(0,nrow(sno.df)))
  sno.df.nor.all = data.frame(rep(0,nrow(sno.df)))
  sno.df.hyp.all = data.frame(rep(0,nrow(sno.df)))
  for (i in 1:8) {
    sno.df.nor.all <- cbind(sno.df.nor.all,sno.df[,c(i+8,i,i+16)])
  }
  
  for (i in 25:32) {
    sno.df.hyp.all <- cbind(sno.df.hyp.all,sno.df[,c(i+8,i,i+16)])
  }
  sno.df.nor.all <- sno.df.nor.all[,-1]
  sno.df.hyp.all <- sno.df.hyp.all[,-1]
  
  ##sno差异分析（对照对对照、处理对处理；组织别）
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(sno.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = sno.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.sno <- data.frame(merge(data.frame(id=rownames(sno.df.hyp.all),sno.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.sno.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.sno.nor <- data.frame(tsi.sno.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.sno.nor <- tsi.sno.nor[,-1]
  colnames(tsi.sno.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.sno.nor <- tsi.sno.nor
  tsi.sno.nor <- data.frame(tsi.sno.nor=TSI(tsi.sno.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.sno.nor
  p.pca.sno.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(sno.df.hyp.all),dex = c(rep("thymus",3),
                                                              rep("kidney",3),
                                                              rep("brain",3),
                                                              rep("liver",3),
                                                              rep("lung",3),
                                                              rep("intestines",3),
                                                              rep("heart",3),
                                                              rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = sno.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.sno.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.sno.hyp <- data.frame(tsi.sno.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.sno.hyp <- tsi.sno.hyp[,-1]
  colnames(tsi.sno.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.sno.hyp <- tsi.sno.hyp
  tsi.sno.hyp <- data.frame(tsi.sno.hyp=TSI(tsi.sno.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.sno.hyp
  p.pca.sno.hyp
  
  #  p.pca.sno.nor|p.pca.sno.hyp -> p.pca.sno.nor.hyp.unite
  
  ##determine tissue specific snoRNA
  tsi.sno.nor  <- data.frame(id=rownames(tsi.sno.nor),tsi.sno.nor,row.names = NULL)
  tsi.sno.nor.spc <- tsi.sno.nor[tsi.sno.nor$tsi.sno.nor>0.85,]
  tsi.sno.hyp  <- data.frame(id=rownames(tsi.sno.hyp),tsi.sno.hyp,row.names = NULL)
  tsi.sno.hyp.spc <- tsi.sno.hyp[tsi.sno.hyp$tsi.sno.hyp>0.85,]
  
  ##determine expression of tissue specific snoRNA
  library(data.table)
  library(ggplot2)
  ##sno对照组
  tsi.sno.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.sno.nor),avgexp.sno.nor),tsi.sno.nor.spc,by="id")
  tsi.sno.nor.spc.long <- melt(tsi.sno.nor.spc.long.1,id=c("id","tsi.sno.nor"))
  colnames(tsi.sno.nor.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.sno.nor.spc.long <- tsi.sno.nor.spc.long[order(tsi.sno.nor.spc.long$TSI,decreasing = T),]
  
  ##不知道为什么转换不了，手动查询后替换
  tsi.sno.nor.spc.long$ID <- gsub("ENSRNOG00000057766","SNORA70",tsi.sno.nor.spc.long$ID)
  
  # library(org.Rn.eg.db)
  # clusterProfiler::bitr(unique(tsi.sno.nor.spc.long$ID),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Rn.eg.db")
  
  ggplot()+
    geom_point(data=tsi.sno.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.sno.nor.spc 
  p.sno.nor.spc 
  
  ##sno处理组
  tsi.sno.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.sno.hyp),avgexp.sno.hyp),tsi.sno.hyp.spc,by="id")
  tsi.sno.hyp.spc.long <- melt(tsi.sno.hyp.spc.long.1,id=c("id","tsi.sno.hyp"))
  colnames(tsi.sno.hyp.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.sno.hyp.spc.long <- tsi.sno.hyp.spc.long[order(tsi.sno.hyp.spc.long$TSI,decreasing = T),]
  
  tsi.sno.hyp.spc.long$ID <- gsub("ENSRNOG00000057766","SNORA70",tsi.sno.hyp.spc.long$ID)
  # library(org.Rn.eg.db)
  # clusterProfiler::bitr(unique(tsi.sno.hyp.spc.long$ID),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Rn.eg.db")
  ggplot()+
    geom_point(data=tsi.sno.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.sno.hyp.spc 
  p.sno.hyp.spc 
  
  tsi.sno.nor.spc.long = data.frame(tsi.sno.nor.spc.long,condition="Control")
  tsi.sno.hyp.spc.long = data.frame(tsi.sno.hyp.spc.long,condition="Hypoxia")
  tsi.sno.long.combine <- rbind(tsi.sno.nor.spc.long,data.frame(tsi.sno.hyp.spc.long[tsi.sno.hyp.spc.long$ID!="SNORA70",-1],ID="SNORA70"))
  ggplot()+
    geom_point(data=tsi.sno.long.combine,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    facet_wrap(~condition,ncol = 2)+
    labs(x="",y="")+
    scale_size_continuous(range = c(-1,6))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
    viridis::scale_color_viridis()->p.sno.combine.spc 
  p.sno.combine.spc 
  # ggsave(filename = "p.sno.combine.spc.pdf",p.sno.combine.spc,units = "cm",width = 16,height = 5)
}##snoRNA

if (T) {
  setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  sca.in <- read.table("sca.quantified",header = T)
  
  sca.in.b <- data.frame(sno.in[,7:ncol(sno.in)],row.names = sno.in$Geneid)
}##scaRNA

if (T) {
  setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  lnc.in <- read.table("lnc.quantified",header = T)
  
  lnc.in.b <- data.frame(lnc.in[,7:ncol(lnc.in)],row.names = lnc.in$Geneid)
  
  lnc.in.name = NULL
   # for (i in 1:nrow(lnc.in.b)) {
   #   lnc.in.name$id[i] <- strsplit(rownames(lnc.in.b),"\\.")[[i]][1]
   # }
  
  
  
  # lnc.in.c <- lnc.in.b[grep("LOC",rownames(lnc.in.b)),]
  # lnc.in.b = lnc.in.c
  
  i=1
  for (i in 1:ncol(lnc.in.b)) {
    colnames(lnc.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(lnc.in.b),"_")[[i]][1],strsplit(colnames(lnc.in.b),"_")[[i]][2],sep = "_")]
  }
  
  lnc.df <- lnc.in.b
  ##lncRNA处理成对照与处理组分开的表格
  lnc.df.all <- data.frame(rep(0,nrow(lnc.df)))
  lnc.df.nor.all = data.frame(rep(0,nrow(lnc.df)))
  lnc.df.hyp.all = data.frame(rep(0,nrow(lnc.df)))
  for (i in 1:8) {
    lnc.df.nor.all <- cbind(lnc.df.nor.all,lnc.df[,c(i+8,i,i+16)])
  }
  
  for (i in 25:32) {
    lnc.df.hyp.all <- cbind(lnc.df.hyp.all,lnc.df[,c(i+8,i,i+16)])
  }
  lnc.df.nor.all <- lnc.df.nor.all[,-1]
  lnc.df.hyp.all <- lnc.df.hyp.all[,-1]
  
  ##火山图
  { 
    lnc.df.all <- cbind(lnc.df.nor.all,lnc.df.hyp.all)
    
    df.col.1 <- data.frame(id=colnames(lnc.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                        "Liver","Lung","Intestines",
                                                                        "Heart","Spleen"),each=3)),2)),
                           condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
    i=7
    a=1
    p.lnc.vol <- NULL  
    for (i in c(seq(1,22,3))) {
      df.col.tem <- df.col.1[c(i,i+1,i+2,i+24,i+25,i+26),]
      
      df.col.tem$condition <- relevel(as.factor(df.col.tem$condition), ref = "Normoxia")
      
      lnc.df.all.tem <- lnc.df.all[,c(i,i+1,i+2,i+24,i+25,i+26)]
      dds.1 <- DESeqDataSetFromMatrix(countData = data.frame(lnc.df.all.tem),
                                      colData = df.col.tem%>%dplyr::select(-tissue), design = ~condition)
      res.tem <- na.omit(data.frame(results(DESeq(dds.1),pAdjustMethod = "BH")))
      
      res.vol.1 <- data.frame(res.tem)%>%filter(baseMean>1)
      
      res.vol.1 <- rownames_to_column(res.vol.1,var = "id")
      
      ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
             ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                    "Both")) -> res.vol.1$exin
      ifelse(res.vol.1$pvalue < .05,1.1,1) -> res.vol.1$size
      
      res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(pvalue)>1)->id.vol
      
      library(ggrepel)
      ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(pvalue)))+
        geom_point(aes(size=-log10(pvalue), color=exin))+
        scale_size_continuous(range = c(0,2.5),
                              # breaks = seq(0,10,2)
        )+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(pvalue)"))+
        geom_hline(yintercept = log10(10),linetype=5)+
        geom_vline(xintercept = .5,color="#f0a1a8")+
        geom_vline(xintercept = -.5,color="#4994c4")+
        theme_minimal()+
        geom_text_repel(data = id.vol%>%filter(-log10(pvalue)>1.5),aes(x=log2FoldChange,y=-log10(pvalue),label = id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
        labs(title = as.character(unique(df.col.tem$tissue)))+
        theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
        # theme(legend.position = legend.justification=c(0,0))+
        scale_y_continuous(limits = c(0,8))+
        scale_x_continuous(limits = c(-2,2))->vol.lnc
      vol.lnc
      
      ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
        geom_point(aes(size=-log10(pvalue), color=exin))+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
        geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
        theme_light()+
        geom_hline(yintercept = .5,color="#f0a1a8")+
        geom_hline(yintercept = -.5,color="#4994c4")+
        scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
        theme(legend.position = "none",
              # axis.text = element_text(size = 10),
              # axis.title = element_text(size=10),
              # axis.line = element_line(color = "black"),
              # plot.background = element_rect(color="black",size=1),
              panel.background = element_rect(color = "white"),
              panel.grid = element_blank())->vol.2.lnc
      vol.2.lnc
      
      vol.lnc+
        # theme(legend.position = "none")+
        annotation_custom(
          grob = ggplotGrob(vol.2.lnc),
          xmin = -2.15,
          xmax = -.10,
          ymin = 4.6,
          ymax = 8.5
        )->uni.vol.lnc
      uni.vol.lnc
      
      p.lnc.vol[[a]] <- uni.vol.lnc
      a <- a+1
    }
    library(patchwork)
    
    p.lnc.vol[[1]]+p.lnc.vol[[2]]+p.lnc.vol[[3]]+p.lnc.vol[[4]]+
      p.lnc.vol[[5]]+p.lnc.vol[[6]]+p.lnc.vol[[7]]+p.lnc.vol[[8]]->a
    
    library(cowplot)
    plot_grid(p.lnc.vol[[1]],p.lnc.vol[[2]],p.lnc.vol[[3]],p.lnc.vol[[4]],
              p.lnc.vol[[5]],p.lnc.vol[[6]],p.lnc.vol[[7]],p.lnc.vol[[8]],
              ncol = 4,
              label_size = 20,
              labels = "AUTO")->a
    
    # ggsave(a,filename = "lnc.vol.patch.pdf",path = "C:\\Users\\colet\\Desktop\\re snc atlas",units = "cm",width = 60,height = 32)
    
    ##火山图
    
    
    dds<-DESeqDataSetFromMatrix(countData = data.frame(lnc.df.all), colData = df.col.1, design = ~tissue+condition)
    
    
    
    # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
    dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
    res = results(dds, pAdjustMethod = "BH")
    res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
    
    
    
    
    
    res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
    library(data.table)
    res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
    ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
           ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                  "Both")) -> res.vol.1$exin
    ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
    
    res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(padj)>1)->id.vol
    
    library(ggrepel)
    ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
      geom_point(aes(size=-log10(padj), color=exin))+
      scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
      scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
      geom_hline(yintercept = log10(10),linetype=5)+
      geom_vline(xintercept = .5,color="#f0a1a8")+
      geom_vline(xintercept = -.5,color="#4994c4")+
      theme_minimal()+
      geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "lncRNA")+
      theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
      # theme(legend.position = legend.justification=c(0,0))+
      scale_y_continuous(limits = c(0,13.5))+
      scale_x_continuous(limits = c(-2,2))->vol.lnc
    vol.lnc
    
    ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
      geom_point(aes(size=-log10(padj), color=exin))+
      scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
      geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
      theme_light()+
      geom_hline(yintercept = .5,color="#f0a1a8")+
      geom_hline(yintercept = -.5,color="#4994c4")+
      scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
      theme(legend.position = "none",
            # axis.text = element_text(size = 10),
            # axis.title = element_text(size=10),
            # axis.line = element_line(color = "black"),
            # plot.background = element_rect(color="black",size=1),
            panel.background = element_rect(color = "white"),
            panel.grid = element_blank())->vol.2.lnc
    vol.2.lnc
    
    vol.lnc+
      # theme(legend.position = "none")+
      annotation_custom(
        grob = ggplotGrob(vol.2.lnc),
        xmin = -2.3,
        xmax = -.10,
        ymin = 7,
        ymax = 14.25
      )->uni.vol.lnc
    
    uni.vol.lnc
  }##火山图
  
  
  
  
  ##lnc差异分析（对照对对照、处理对处理；组织别）
  
  if (T) {
    
    ##0.1all找常氧缺氧比对
    if (T) {
      
      ##0.1all找常氧缺氧比对
      # df.col.a <- data.frame(id=colnames(lnc.df),tissue=c(rep(c("Thymus","Kidney","Brain",
      #                                                              "Liver","Lung","Intestines",
      #                                                              "Heart","Spleen"),6)),
      #                        condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),1)))
      df.col.1 <- data.frame(id=colnames(lnc.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                          "Liver","Lung","Intestines",
                                                                          "Heart","Spleen"),each=3)),2)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
      dds<-DESeqDataSetFromMatrix(countData = data.frame(lnc.df.all), colData = df.col.1, design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
      
      
      ##火山图
      { 
        df.col.1 <- data.frame(id=colnames(lnc.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                            "Liver","Lung","Intestines",
                                                                            "Heart","Spleen"),each=3)),2)),
                               condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
        i=7
        a=1
        p.lnc.vol <- NULL  
        for (i in c(seq(1,22,3))) {
          df.col.tem <- df.col.1[c(i,i+1,i+2,i+24,i+25,i+26),]
          
          df.col.tem$condition <- relevel(as.factor(df.col.tem$condition), ref = "Normoxia")
          
          lnc.df.all.tem <- lnc.df.all[,c(i,i+1,i+2,i+24,i+25,i+26)]
          dds.1 <- DESeqDataSetFromMatrix(countData = data.frame(lnc.df.all.tem),
                                          colData = df.col.tem[,-2], design = ~condition)
          res.tem <- na.omit(data.frame(results(DESeq(dds.1),pAdjustMethod = "BH")))
          
          res.vol.1 <- data.frame(res.tem)%>%filter(baseMean>1)
          
          res.vol.1 <- rownames_to_column(res.vol.1,var = "id")
          
          ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
                 ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                        "Both")) -> res.vol.1$exin
          ifelse(res.vol.1$pvalue < .05,1.1,1) -> res.vol.1$size
          
          res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(pvalue)>1)->id.vol
          
          library(ggrepel)
          ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(pvalue)))+
            geom_point(aes(size=-log10(pvalue), color=exin))+
            scale_size_continuous(range = c(0,2.5),
                                  # breaks = seq(0,10,2)
            )+
            scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
            guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(pvalue)"))+
            geom_hline(yintercept = log10(10),linetype=5)+
            geom_vline(xintercept = .5,color="#f0a1a8")+
            geom_vline(xintercept = -.5,color="#4994c4")+
            theme_minimal()+
            geom_text_repel(data = id.vol%>%filter(-log10(pvalue)>1.5),aes(x=log2FoldChange,y=-log10(pvalue),label = id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
            labs(title = as.character(unique(df.col.tem$tissue)))+
            theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
            # theme(legend.position = legend.justification=c(0,0))+
            scale_y_continuous(limits = c(0,8))+
            scale_x_continuous(limits = c(-2,2))->vol.lnc
          vol.lnc
          
          ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
            geom_point(aes(size=-log10(pvalue), color=exin))+
            scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
            guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
            geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
            theme_light()+
            geom_hline(yintercept = .5,color="#f0a1a8")+
            geom_hline(yintercept = -.5,color="#4994c4")+
            scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
            theme(legend.position = "none",
                  # axis.text = element_text(size = 10),
                  # axis.title = element_text(size=10),
                  # axis.line = element_line(color = "black"),
                  # plot.background = element_rect(color="black",size=1),
                  panel.background = element_rect(color = "white"),
                  panel.grid = element_blank())->vol.2.lnc
          vol.2.lnc
          
          vol.lnc+
            # theme(legend.position = "none")+
            annotation_custom(
              grob = ggplotGrob(vol.2.lnc),
              xmin = -2.15,
              xmax = -.10,
              ymin = 4.6,
              ymax = 8.5
            )->uni.vol.lnc
          uni.vol.lnc
          
          p.lnc.vol[[a]] <- uni.vol.lnc
          a <- a+1
        }
        library(patchwork)
        
        p.lnc.vol[[1]]+p.lnc.vol[[2]]+p.lnc.vol[[3]]+p.lnc.vol[[4]]+
          p.lnc.vol[[5]]+p.lnc.vol[[6]]+p.lnc.vol[[7]]+p.lnc.vol[[8]]->a
        
        library(cowplot)
        plot_grid(p.lnc.vol[[1]],p.lnc.vol[[2]],p.lnc.vol[[3]],p.lnc.vol[[4]],
                  p.lnc.vol[[5]],p.lnc.vol[[6]],p.lnc.vol[[7]],p.lnc.vol[[8]],
                  ncol = 4,
                  label_size = 20,
                  labels = "AUTO")->a
        
        # ggsave(a,filename = "lnc.vol.patch.pdf",path = "C:\\Users\\colet\\Desktop\\re snc atlas",units = "cm",width = 60,height = 32)
        
        ##火山图
        
        
        dds<-DESeqDataSetFromMatrix(countData = data.frame(lnc.df.all), colData = df.col.1, design = ~tissue+condition)
        
        
        
        # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
        dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
        res = results(dds, pAdjustMethod = "BH")
        res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
        
        
        
        
        
        res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
        library(data.table)
        res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
        ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
               ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                      "Both")) -> res.vol.1$exin
        ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
        
        res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(padj)>1)->id.vol
        
        library(ggrepel)
        ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
          geom_point(aes(size=-log10(padj), color=exin))+
          scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
          scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
          guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
          geom_hline(yintercept = log10(10),linetype=5)+
          geom_vline(xintercept = .5,color="#f0a1a8")+
          geom_vline(xintercept = -.5,color="#4994c4")+
          theme_minimal()+
          geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
          labs(title = "lncRNA")+
          theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
          # theme(legend.position = legend.justification=c(0,0))+
          scale_y_continuous(limits = c(0,13.5))+
          scale_x_continuous(limits = c(-2,2))->vol.lnc
        vol.lnc
        
        ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
          geom_point(aes(size=-log10(padj), color=exin))+
          scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
          guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
          geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
          theme_light()+
          geom_hline(yintercept = .5,color="#f0a1a8")+
          geom_hline(yintercept = -.5,color="#4994c4")+
          scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
          theme(legend.position = "none",
                # axis.text = element_text(size = 10),
                # axis.title = element_text(size=10),
                # axis.line = element_line(color = "black"),
                # plot.background = element_rect(color="black",size=1),
                panel.background = element_rect(color = "white"),
                panel.grid = element_blank())->vol.2.lnc
        vol.2.lnc
        
        vol.lnc+
          # theme(legend.position = "none")+
          annotation_custom(
            grob = ggplotGrob(vol.2.lnc),
            xmin = -2.3,
            xmax = -.10,
            ymin = 7,
            ymax = 14.25
          )->uni.vol.lnc
        
        uni.vol.lnc
      }##火山图
      
      
      
      lncdf.sig <- lnc.df[rownames(lnc.df)%in%rownames(res.df[abs(res.df$log2FoldChange)>2,]),]
      lncdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        lncdf.sig.mean.1  <-  rowMeans(matrix(lncdf.sig[,c(i:i+2)]))
        lncdf.sig.mean <- cbind(lncdf.sig.mean,lncdf.sig.mean.1)
      } 
      colnames(lncdf.sig.mean) <- colnames(lncdf.sig)[seq(1,48,3)]
      # TSI(lncdf.sig)
      lncdf.sig$id <- rownames(lncdf.sig)
      
      library(data.table)
      lncdf.sig.long <- melt(lncdf.sig,id="id")
      name.idx <- strsplit(as.character(lncdf.sig.long$variable),"_")
      for (i in 1:nrow(lncdf.sig.long)) {
        lncdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        lncdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      lncdf.sig.long$condition[lncdf.sig.long$condition=="NOR"] <- "Normoxia"
      lncdf.sig.long$condition[lncdf.sig.long$condition=="HYP"] <- "Hypoxia"
      library(dplyr)
      
      lncdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))
      
      ##两个柱子挨在一起的图
      if (T) {
        library(ggplot2)
        # lncdf.sig.long%>%filter(value>0)%>%filter(!grepl('let', id))%>%
        lncdf.sig.long%>%filter(value>0&grepl('E',id))%>%
          ggplot(aes(x=id,y=log2(value)))+
          geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=condition))+
          stat_summary(fun = mean,
                       geom = "errorbar",
                       fun.max = function(x) mean(x) + sd(x),
                       fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
                       aes(group=condition,color=condition),
                       size=.5,
                       width=.25,
          )+
          # geom_vline(xintercept = .5)+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
          facet_wrap(~tissue,ncol=1,strip.position = "right")+
          labs(y="log2(CPM)",x=NULL,title = "lncRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "2barEXP.lnc.pdf",units = "cm",width = 16,height = 20,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子挨在一起的图
      
      ##两个柱子堆积成1
      if (T) {
        library(ggplot2)
        # lncdf.sig.long%>%filter(value>0)%>%
        lncdf.sig.long%>%filter(value>0&grepl('E',id))%>%
          ggplot(aes(y=id,x=log10(value)))+
          
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "lncRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          geom_vline(xintercept = .5,linetype = 5)+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "fillbarEXP.lnc.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (T) {
        library(ggplot2)
        # lncdf.sig.long%>%filter(value>0)%>%
        lncdf.sig.long%>%filter(value>20&grepl('E',id))%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "lncRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "stackbarEXP.lnc.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
      
      # # Create a tibble for LRT results
      # res_LRT_tb <- res %>%
      #   data.frame() %>%
      #   rownames_to_column(var="gene") %>% 
      #   as_tibble()
      # # Subset to return genes with padj < 0.05
      # sigLRT_genes <- res_LRT_tb %>% 
      #   filter(padj < .05)
      # 
      # # Get number of significant genes
      # nrow(sigLRT_genes)
      # 
      # # Subset results for faster cluster finding (for classroom demo purposes)
      # clustering_sig_genes <- sigLRT_genes %>%
      #   arrange(padj) %>%
      #   head(n=1000)
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
      # library(DEGreport)
      # # install.packages('C:/Users/colet/Desktop/lasso2_1.2-22.tar.gz', repo=NULL)
      # clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
      # 
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # 
      # BiocManager::install("DEGreport")
    }
    
    ##0.1all找组织比对(第一个是点图)
    if (T) {
      ##先算TSI
      df.col = data.frame(id = colnames(lnc.df.nor.all),dex = c(rep("thymus",3),
                                                                rep("kidney",3),
                                                                rep("brain",3),
                                                                rep("liver",3),
                                                                rep("lung",3),
                                                                rep("intestines",3),
                                                                rep("heart",3),
                                                                rep("spleen",3)
      ))
      dds<-DESeqDataSetFromMatrix(countData = lnc.df.nor.all, colData = df.col,design = ~dex)
      ## Prefiltering
      filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
      dds <- dds[!filt,]
      
      ## Perform DESeq2()
      dds = DESeq(dds)
      res = results(dds, pAdjustMethod = "BH")
      
      ##determine tissue specific(calculate TSI)
      # res.tsi <-data.frame(res)
      # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
      # res.sig <- na.omit(data.frame(res.sig))
      # res.sig.lnc <- data.frame(merge(data.frame(id=rownames(lnc.df.hyp.all),lnc.df.nor.all),data.frame(id=rownames(res.sig))
      #                                ,by="id"),row.names = 1)
      
      tpm <- function(x){
        tpmres <- (x[,6:ncol(x)]/x[,5])/colSums(x[,6:ncol(x)]/x[,5])*1e6
        return(tpmres)
      }
      
      library(edgeR)
      # apply variance stabilizing transformation
      # vsted <- cpm(counts(dds))
      v = varianceStabilizingTransformation(dds, blind=FALSE)
      vsted = assay(v)
      
      ##calculate TSI
      tsi.lnc.nor <- rep(0,nrow(vsted))
      for (i in c(1,4,7,10,13,16,19,21)) {
        tsi.lnc.nor <- data.frame(tsi.lnc.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
      }
      tsi.lnc.nor <- tsi.lnc.nor[,-1]
      colnames(tsi.lnc.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
      avgexp.lnc.nor <- tsi.lnc.nor
      tsi.lnc.nor <- data.frame(TSI=TSI(tsi.lnc.nor,8))
      
      tsi.sig.id <- rownames_to_column(tsi.lnc.nor,var="id")%>%filter(TSI>.85)
      
      ##0.1all找组织比对
      # df.col.a <- data.frame(id=colnames(lnc.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
      #                                                              "Liver","Lung","Intestines",
      #                                                              "Heart","Spleen"),each=3)),2)),
      #                        condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),2)))
      df.col.a <- data.frame(id=colnames(lnc.df),tissue=c(rep(c("Thymus","Kidney","Brain",
                                                                "Liver","Lung","Intestines",
                                                                "Heart","Spleen"),6)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),1)))
      dds<-DESeqDataSetFromMatrix(countData = lnc.df, colData = df.col.a,design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~condition, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
      
      # res.df.1 <- res.df
      # res.df.1$id <- rownames(res.df.1)
      
      lncdf.sig <- lnc.df[rownames(lnc.df)%in%rownames(res.df[abs(res.df$log2FoldChange)>0,]),]
      lncdf.sig <- merge(rownames_to_column(lncdf.sig,"id"),tsi.sig.id,by="id")
      lncdf.sig <- column_to_rownames(lncdf.sig,"id")
      # rownames(lncdf.sig) <- lncdf.sig$id
      lncdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        lncdf.sig.mean.1  <-  rowMeans(matrix(lncdf.sig[,c(i:i+2)]))
        lncdf.sig.mean <- cbind(lncdf.sig.mean,lncdf.sig.mean.1)
      } 
      colnames(lncdf.sig.mean) <- colnames(lncdf.sig)[seq(1,48,3)]
      # TSI(lncdf.sig)
      lncdf.sig$id <- rownames(lncdf.sig)
      
      library(data.table)
      lncdf.sig.long <- melt(lncdf.sig,id=c("id","TSI"))
      name.idx <- strsplit(as.character(lncdf.sig.long$variable),"_")
      for (i in 1:nrow(lncdf.sig.long)) {
        lncdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        lncdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      lncdf.sig.long$condition[lncdf.sig.long$condition=="NOR"] <- "Normoxia"
      lncdf.sig.long$condition[lncdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      lncdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))->p.dat
      lnc.dot.dat <- data.frame(merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id"),RNA="lncRNA")
      rbind(rny.dot.dat,lnc.dot.dat)
      # dcast(data = p.dat,id~condition+tissue+value_mean+sd+se) 
      ##点图（全部）
      if (T) {
        library(ggplot2)
        # lncdf.sig.long%>%filter(value>0&condition=="Normoxia")%>%
        # p.dat%>%filter(value>0&condition=="Normoxia")%>%
        p.dat.log  <- merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id")%>%
          filter(condition=="Normoxia")%>%
          # log10(p.dat.log$value_mean)
          # rbind(rny.dot.dat,lnc.dot.dat)%>%
          #   filter(condition=="Normoxia")%>%
          ggplot(aes(y=id,x=tissue))+
          geom_point(aes(size=log(value_mean, 10),color=padj))+
          scale_size_continuous(range = c(-1,8),breaks = seq(.5,3,.5))+
          # facet_grid(~RNA)+
          # facet_wrap(~RNA,ncol = 1)+.
          # geom_point(fun = "mean", stat = "summary",aes(group=tissue,fill=condition,color=condition,size=log(value)),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=id,fill=condition,color=condition,size=log(value)))+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
        # )+
        # geom_curve()
        # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
        # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
        # facet_wrap(~tissue,nrow=1,strip.position = "top")+
        labs(y=NULL,x=NULL,title = "lncRNA")+
          # scale_color_gradient(low = "blue",high = "red")+
          viridis::scale_color_viridis()+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(size=guide_legend(reverse = F,title = "log10(CPM)"),
                 # color=guide_legend(title = "p.adj")
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(),
                axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
          )->p.lnc.tsi.dot
        p.lnc.tsi.dot
        
        options(download.file.method = 'libcurl')
        options(url.method='libcurl')
        # BiocManager::install("lemon")
        library(lemon)
        cowplot::plot_grid()
        
        unique(rny.dot.dat$id)
        
        
        #####
        #拼点图 TSI
        patchwork::wrap_plots(guides = "collect",p.rny.tsi.dot, p.lnc.tsi.dot)
        grid_arrange_shared_legend(, p.lnc.tsi.dot, p.pir.tsi.dot,nrow = 3,ncol = 1,position='right')
        
        p1=p.rny.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(rny.dot.dat$id))/1.5,"cm"))
        p2=p.mir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(mir.dot.dat$id))/1.5,"cm"))
        p3=p.pir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(pir.dot.dat$id))/1.5,"cm"))
        p4=p.lnc.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(lnc.dot.dat$id))/1.5,"cm"))
        
        
        p.trf.tsi.dot
        
        
        p <- grid_arrange_shared_legend(p1,p2,nrow = 2,ncol = 1,position='right')
        
        # ggsave(plot = p,filename = "test.pdf",units = "cm",width = 60,height = 500,path = "C:\\Users\\colet\\Desktop\\re snc atlas",limitsize = F)
        p <- patchwork::wrap_plots(guides = "collect",p1, p2,p3,p4,ncol = 1)
        #拼点图 TSI
        #####      
        # ggsave(filename = "TIS.dot.EXP.lnc.pdf",units = "cm",width = 16,height = 55,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##点图（全部）
      
      ##点图（部分）
      if (T) {
        library(ggplot2)
        # lncdf.sig.long%>%filter(value>0&condition=="Normoxia")%>%
        # p.dat%>%filter(value>0&condition=="Normoxia")%>%
        p.dat.log  <- merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id")%>%
          filter(condition=="Normoxia")
        
        p.dat.log[grep(pattern = "ENSRAT|LOC",p.dat.log$id),]%>%filter(value_mean>20)->p.dat.log.1
        
        strsplit2(x = p.dat.log.1$id,split = "..Sep")[,1]->p.dat.log.1$id
          # log10(p.dat.log$value_mean)
          # rbind(rny.dot.dat,lnc.dot.dat)%>%
          #   filter(condition=="Normoxia")%>%
         
        p.dat.log.1%>%
         ggplot(aes(y=id,x=tissue))+
          geom_point(aes(size=log(value_mean, 10),color=padj))+
          scale_size_continuous(range = c(-1,8),breaks = seq(.5,3,.5))+
          # facet_grid(~RNA)+
          # facet_wrap(~RNA,ncol = 1)+.
          # geom_point(fun = "mean", stat = "summary",aes(group=tissue,fill=condition,color=condition,size=log(value)),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=id,fill=condition,color=condition,size=log(value)))+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
        # )+
        # geom_curve()
        # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
        # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
        # facet_wrap(~tissue,nrow=1,strip.position = "top")+
        labs(y=NULL,x=NULL,title = "lncRNA")+
          # scale_color_gradient(low = "blue",high = "red")+
          viridis::scale_color_viridis()+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(size=guide_legend(reverse = F,title = "log10(CPM)"),
                 # color=guide_legend(title = "p.adj")
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(),
                axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
          )->p.lnc.tsi.dot
        p.lnc.tsi.dot
        
        options(download.file.method = 'libcurl')
        options(url.method='libcurl')
        # BiocManager::install("lemon")
        library(lemon)
        cowplot::plot_grid()
        
        unique(rny.dot.dat$id)
        
        
        #####
        #拼点图 TSI
        patchwork::wrap_plots(guides = "collect",p.rny.tsi.dot, p.lnc.tsi.dot)
        grid_arrange_shared_legend(, p.lnc.tsi.dot, p.pir.tsi.dot,nrow = 3,ncol = 1,position='right')
        
        # p1=p.rny.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(rny.dot.dat$id))/1.5,"cm"))
        # p2=p.mir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(mir.dot.dat$id))/1.5,"cm"))
        # p3=p.pir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(pir.dot.dat$id))/1.5,"cm"))
        # p4=p.lnc.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(lnc.dot.dat$id))/1.5,"cm"))
        
        # p1=p.rny.tsi.dot%>%egg::set_panel_size(width = unit(7,"cm"),height = unit(length(unique(rny.dot.dat$id))/1.5,"cm"))
        p2=p.mir.tsi.dot.selected%>%egg::set_panel_size(width = unit(6,"cm"),height = unit(length(unique(mir.dot.dat$id))/4.2,"cm"))
        # p3=p.pir.tsi.dot%>%egg::set_panel_size(width = unit(7,"cm"),height = unit(length(unique(pir.dot.dat$id))/1.5,"cm"))
        p4=p.lnc.tsi.dot%>%egg::set_panel_size(width = unit(6,"cm"),height = unit(length(unique(p.dat.log.1$id))/4.2,"cm"))
        p5=p.other.tsi.dot%>%egg::set_panel_size(width = unit(6,"cm"),height = unit(length(unique(other.dot.dat$id))/2.4,"cm"))
        
        p.other.tsi.dot
        p.trf.tsi.dot
        
        
        p <- grid_arrange_shared_legend(p1,p2,nrow = 2,ncol = 1,position='right')
        
        # ggsave(plot = p,filename = "test.pdf",units = "cm",width = 60,height = 500,path = "C:\\Users\\colet\\Desktop\\re snc atlas",limitsize = F)
        p <- patchwork::wrap_plots(guides = "collect",p2,p4,p5,nrow = 1)
        p
        #拼点图 TSI
        #####      
        # ggsave(filename = "TIS.dot.EXP.lnc.pdf",units = "cm",width = 16,height = 55,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##点图（全部）
      
      
      
      
      ##两个柱子堆积成1
      if (F) {
        library(ggplot2)
        lncdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.fillbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (F) {
        library(ggplot2)
        lncdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.stackbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
    }
    
  }##新加入新分析
  
  
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(lnc.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = lnc.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.lnc <- data.frame(merge(data.frame(id=rownames(lnc.df.hyp.all),lnc.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.lnc.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.lnc.nor <- data.frame(tsi.lnc.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.lnc.nor <- tsi.lnc.nor[,-1]
  colnames(tsi.lnc.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.lnc.nor <- tsi.lnc.nor
  tsi.lnc.nor <- data.frame(tsi.lnc.nor=TSI(tsi.lnc.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.lnc.nor
  p.pca.lnc.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(lnc.df.hyp.all),dex = c(rep("thymus",3),
                                                              rep("kidney",3),
                                                              rep("brain",3),
                                                              rep("liver",3),
                                                              rep("lung",3),
                                                              rep("intestines",3),
                                                              rep("heart",3),
                                                              rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = lnc.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.lnc.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.lnc.hyp <- data.frame(tsi.lnc.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.lnc.hyp <- tsi.lnc.hyp[,-1]
  colnames(tsi.lnc.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.lnc.hyp <- tsi.lnc.hyp
  tsi.lnc.hyp <- data.frame(tsi.lnc.hyp=TSI(tsi.lnc.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.lnc.hyp
  p.pca.lnc.hyp
  
  # p.pca.lnc.nor|p.pca.lnc.hyp -> p.pca.lnc.nor.hyp.unite
  
  ##determine tissue specific lncRNA
  tsi.lnc.nor  <- data.frame(id=rownames(tsi.lnc.nor),tsi.lnc.nor,row.names = NULL)
  tsi.lnc.nor.spc <- tsi.lnc.nor[tsi.lnc.nor$tsi.lnc.nor>0.85,]
  tsi.lnc.hyp  <- data.frame(id=rownames(tsi.lnc.hyp),tsi.lnc.hyp,row.names = NULL)
  tsi.lnc.hyp.spc <- tsi.lnc.hyp[tsi.lnc.hyp$tsi.lnc.hyp>0.85,]
  
  ##determine expression of tissue specific lncRNA
  library(data.table)
  library(ggplot2)
  ##lnc对照组
  tsi.lnc.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.lnc.nor),avgexp.lnc.nor),tsi.lnc.nor.spc,by="id")
  
  #取前100
  tsi.lnc.nor.spc.long.1 <- head(tsi.lnc.nor.spc.long.1[order(tsi.lnc.nor.spc.long.1$tsi.lnc.nor, decreasing = T),],100)
  
  tsi.lnc.nor.spc.long <- melt(tsi.lnc.nor.spc.long.1,id=c("id","tsi.lnc.nor"))
  colnames(tsi.lnc.nor.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.lnc.nor.spc.long <- tsi.lnc.nor.spc.long[order(tsi.lnc.nor.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.lnc.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.lnc.nor.spc 
  p.lnc.nor.spc 
  
  ##lnc处理组
  tsi.lnc.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.lnc.hyp),avgexp.lnc.hyp),tsi.lnc.hyp.spc,by="id")
  
  #取前100
  tsi.lnc.hyp.spc.long.1 <- head(tsi.lnc.hyp.spc.long.1[order(tsi.lnc.hyp.spc.long.1$tsi.lnc.hyp, decreasing = T),],100)
  
  tsi.lnc.hyp.spc.long <- melt(tsi.lnc.hyp.spc.long.1,id=c("id","tsi.lnc.hyp"))
  colnames(tsi.lnc.hyp.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.lnc.hyp.spc.long <- tsi.lnc.hyp.spc.long[order(tsi.lnc.hyp.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.lnc.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.lnc.hyp.spc 
  p.lnc.hyp.spc 
  
  tsi.lnc.nor.spc.long <- data.frame(tsi.lnc.nor.spc.long,condition="Control")
  tsi.lnc.hyp.spc.long <- data.frame(tsi.lnc.hyp.spc.long,condition="Hypoxia")
  tsi.lnc.combine.spc.long <- rbind(tsi.lnc.nor.spc.long,tsi.lnc.hyp.spc.long)
  
  ggplot()+
    geom_point(data=tsi.lnc.combine.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    facet_wrap(~condition)+
    labs(x="",y="")+
    scale_size_continuous(range = c(-1,6))+
    viridis::scale_color_viridis()+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))->p.lnc.combine.spc 
  p.lnc.combine.spc 
  
 # ggsave(p.lnc.combine.spc,filename = "p.lnc.combine.spc.pdf",units = "cm",width = 16,height = 16,limitsize = F)
}##lncRNA

if (T) {
  setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  trf.in <- read.table("tRNA_all.quantified",header = T)
  
  trf.in.b <- data.frame(trf.in[,7:ncol(trf.in)],row.names = trf.in$Geneid)
  
  i=1
  for (i in 1:ncol(trf.in.b)) {
    colnames(trf.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(trf.in.b),"_")[[i]][1],strsplit(colnames(trf.in.b),"_")[[i]][2],sep = "_")]
  }
  
  trf.df <- trf.in.b
  ##trf处理成对照与处理组分开的表格
  trf.df.all <- data.frame(rep(0,nrow(trf.df)))
  trf.df.nor.all = data.frame(rep(0,nrow(trf.df)))
  trf.df.hyp.all = data.frame(rep(0,nrow(trf.df)))
  for (i in 1:8) {
    trf.df.nor.all <- cbind(trf.df.nor.all,trf.df[,c(i+8,i,i+16)])
  }
  
 
  for (i in 25:32) {
    trf.df.hyp.all <- cbind(trf.df.hyp.all,trf.df[,c(i,i+8,i+16)])
  }
  trf.df.nor.all <- trf.df.nor.all[,-1]
  trf.df.hyp.all <- trf.df.hyp.all[,-1]
 
  ##火山图
  { 
    trf.df.all <- cbind(trf.df.nor.all,trf.df.hyp.all)
    df.col.1 <- data.frame(id=colnames(trf.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                        "Liver","Lung","Intestines",
                                                                        "Heart","Spleen"),each=3)),2)),
                           condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
    i=7
    a=1
    p.trf.vol <- NULL  
    for (i in c(seq(1,22,3))) {
      df.col.tem <- df.col.1[c(i,i+1,i+2,i+24,i+25,i+26),]
      
      df.col.tem$condition <- relevel(as.factor(df.col.tem$condition), ref = "Normoxia")
      
      trf.df.all.tem <- trf.df.all[,c(i,i+1,i+2,i+24,i+25,i+26)]
      dds.1 <- DESeqDataSetFromMatrix(countData = data.frame(trf.df.all.tem),
                                      colData = df.col.tem%>%dplyr::select(-tissue), design = ~condition)
      res.tem <- na.omit(data.frame(results(DESeq(dds.1),pAdjustMethod = "BH")))
      
      res.vol.1 <- data.frame(res.tem)%>%filter(baseMean>1)
      
      res.vol.1 <- rownames_to_column(res.vol.1,var = "id")
      
      ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
             ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                    "Both")) -> res.vol.1$exin
      ifelse(res.vol.1$pvalue < .05,1.1,1) -> res.vol.1$size
      
      res.vol.1%>%filter(abs(log2FoldChange)>.5)->id.vol
      
      library(ggrepel)
      ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(pvalue)))+
        geom_point(aes(size=-log10(pvalue), color=exin))+
        scale_size_continuous(range = c(0,2.5),
                              # breaks = seq(0,10,2)
        )+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(pvalue)"))+
        geom_hline(yintercept = log10(10),linetype=5)+
        geom_vline(xintercept = .5,color="#f0a1a8")+
        geom_vline(xintercept = -.5,color="#4994c4")+
        theme_minimal()+
        geom_text_repel(data = id.vol%>%filter(-log10(pvalue)>1.5),aes(x=log2FoldChange,y=-log10(pvalue),label = id),max.overlaps = 1000,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
        labs(title = as.character(unique(df.col.tem$tissue)))+
        theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
        # theme(legend.position = legend.justification=c(0,0))+
        scale_y_continuous(limits = c(0,8))+
        scale_x_continuous(limits = c(-2,2))->vol.trf
      vol.trf
      
      ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
        geom_point(aes(size=-log10(pvalue), color=exin))+
        scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
        guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
        geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
        theme_light()+
        geom_hline(yintercept = .5,color="#f0a1a8")+
        geom_hline(yintercept = -.5,color="#4994c4")+
        scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
        theme(legend.position = "none",
              # axis.text = element_text(size = 10),
              # axis.title = element_text(size=10),
              # axis.line = element_line(color = "black"),
              # plot.background = element_rect(color="black",size=1),
              panel.background = element_rect(color = "white"),
              panel.grid = element_blank())->vol.2.trf
      vol.2.trf
      
      vol.trf+
        # theme(legend.position = "none")+
        annotation_custom(
          grob = ggplotGrob(vol.2.trf),
          xmin = -2.15,
          xmax = -.10,
          ymin = 4.6,
          ymax = 8.5
        )->uni.vol.trf
      uni.vol.trf
      
      p.trf.vol[[a]] <- uni.vol.trf
      a <- a+1
    }
    library(patchwork)
    
    p.trf.vol[[1]]+p.trf.vol[[2]]+p.trf.vol[[3]]+p.trf.vol[[4]]+
      p.trf.vol[[5]]+p.trf.vol[[6]]+p.trf.vol[[7]]+p.trf.vol[[8]]->a
    
    library(cowplot)
    plot_grid(p.trf.vol[[1]],p.trf.vol[[2]],p.trf.vol[[3]],p.trf.vol[[4]],
              p.trf.vol[[5]],p.trf.vol[[6]],p.trf.vol[[7]],p.trf.vol[[8]],
              ncol = 4,
              label_size = 20,
              labels = "AUTO")->a
    
    # ggsave(a,filename = "trf.vol.patch.pdf",path = "C:\\Users\\colet\\Desktop\\re snc atlas",units = "cm",width = 60,height = 32)
    
    ##火山图
    
    
    dds<-DESeqDataSetFromMatrix(countData = data.frame(trf.df.all), colData = df.col.1, design = ~tissue+condition)
    
    
    
    # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
    dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
    res = results(dds, pAdjustMethod = "BH")
    res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
    
    
    
    
    
    res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
    library(data.table)
    res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
    ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
           ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                  "Both")) -> res.vol.1$exin
    ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
    
    res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(padj)>1)->id.vol
    
    library(ggrepel)
    ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
      geom_point(aes(size=-log10(padj), color=exin))+
      scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
      scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
      geom_hline(yintercept = log10(10),linetype=5)+
      geom_vline(xintercept = .5,color="#f0a1a8")+
      geom_vline(xintercept = -.5,color="#4994c4")+
      theme_minimal()+
      geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "tRNA")+
      theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
      # theme(legend.position = legend.justification=c(0,0))+
      scale_y_continuous(limits = c(0,13.5))+
      scale_x_continuous(limits = c(-2,2))->vol.trf
    vol.trf
    
    ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
      geom_point(aes(size=-log10(padj), color=exin))+
      scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
      geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
      theme_light()+
      geom_hline(yintercept = .5,color="#f0a1a8")+
      geom_hline(yintercept = -.5,color="#4994c4")+
      scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
      theme(legend.position = "none",
            # axis.text = element_text(size = 10),
            # axis.title = element_text(size=10),
            # axis.line = element_line(color = "black"),
            # plot.background = element_rect(color="black",size=1),
            panel.background = element_rect(color = "white"),
            panel.grid = element_blank())->vol.2.trf
    vol.2.trf
    
    vol.trf+
      # theme(legend.position = "none")+
      annotation_custom(
        grob = ggplotGrob(vol.2.trf),
        xmin = -2.3,
        xmax = -.10,
        ymin = 7,
        ymax = 14.25
      )->uni.vol.trf
    
    uni.vol.trf
  }##火山图
  
  
  ##trf差异分析（对照对对照、处理对处理；组织别）
  
  if (T) {
    
    ##0.1all找常氧缺氧比对
    if (T) {
      
      ##0.1all找常氧缺氧比对
      # df.col.a <- data.frame(id=colnames(trf.df),tissue=c(rep(c("Thymus","Kidney","Brain",
      #                                                              "Liver","Lung","Intestines",
      #                                                              "Heart","Spleen"),6)),
      #                        condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),1)))
      df.col.1 <- data.frame(id=colnames(trf.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                          "Liver","Lung","Intestines",
                                                                          "Heart","Spleen"),each=3)),2)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
      dds<-DESeqDataSetFromMatrix(countData = data.frame(trf.df.all), colData = df.col.1, design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
      
      
      ##火山图
      { 
        df.col.1 <- data.frame(id=colnames(trf.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
                                                                            "Liver","Lung","Intestines",
                                                                            "Heart","Spleen"),each=3)),2)),
                               condition=c(rep(rep(c("Normoxia","Hypoxia"),each=24),1)))
        i=7
        a=1
        p.trf.vol <- NULL  
        for (i in c(seq(1,22,3))) {
          df.col.tem <- df.col.1[c(i,i+1,i+2,i+24,i+25,i+26),]
          
          df.col.tem$condition <- relevel(as.factor(df.col.tem$condition), ref = "Normoxia")
          
          trf.df.all.tem <- trf.df.all[,c(i,i+1,i+2,i+24,i+25,i+26)]
          dds.1 <- DESeqDataSetFromMatrix(countData = data.frame(trf.df.all.tem),
                                          colData = df.col.tem%>%dplyr::select(-tissue), design = ~condition)
          res.tem <- na.omit(data.frame(results(DESeq(dds.1),pAdjustMethod = "BH")))
          
          res.vol.1 <- data.frame(res.tem)%>%filter(baseMean>1)
          
          res.vol.1 <- rownames_to_column(res.vol.1,var = "id")
          
          ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
                 ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                        "Both")) -> res.vol.1$exin
          ifelse(res.vol.1$pvalue < .05,1.1,1) -> res.vol.1$size
          
          res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(pvalue)>1)->id.vol
          
          library(ggrepel)
          ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(pvalue)))+
            geom_point(aes(size=-log10(pvalue), color=exin))+
            scale_size_continuous(range = c(0,2.5),
                                  # breaks = seq(0,10,2)
            )+
            scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
            guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(pvalue)"))+
            geom_hline(yintercept = log10(10),linetype=5)+
            geom_vline(xintercept = .5,color="#f0a1a8")+
            geom_vline(xintercept = -.5,color="#4994c4")+
            theme_minimal()+
            geom_text_repel(data = id.vol%>%filter(-log10(pvalue)>1.5),aes(x=log2FoldChange,y=-log10(pvalue),label = id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
            labs(title = as.character(unique(df.col.tem$tissue)))+
            theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
            # theme(legend.position = legend.justification=c(0,0))+
            scale_y_continuous(limits = c(0,8))+
            scale_x_continuous(limits = c(-2,2))->vol.trf
          vol.trf
          
          ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
            geom_point(aes(size=-log10(pvalue), color=exin))+
            scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
            guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
            geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
            theme_light()+
            geom_hline(yintercept = .5,color="#f0a1a8")+
            geom_hline(yintercept = -.5,color="#4994c4")+
            scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
            theme(legend.position = "none",
                  # axis.text = element_text(size = 10),
                  # axis.title = element_text(size=10),
                  # axis.line = element_line(color = "black"),
                  # plot.background = element_rect(color="black",size=1),
                  panel.background = element_rect(color = "white"),
                  panel.grid = element_blank())->vol.2.trf
          vol.2.trf
          
          vol.trf+
            # theme(legend.position = "none")+
            annotation_custom(
              grob = ggplotGrob(vol.2.trf),
              xmin = -2.15,
              xmax = -.10,
              ymin = 4.6,
              ymax = 8.5
            )->uni.vol.trf
          uni.vol.trf
          
          p.trf.vol[[a]] <- uni.vol.trf
          a <- a+1
        }
        library(patchwork)
        
        p.trf.vol[[1]]+p.trf.vol[[2]]+p.trf.vol[[3]]+p.trf.vol[[4]]+
          p.trf.vol[[5]]+p.trf.vol[[6]]+p.trf.vol[[7]]+p.trf.vol[[8]]->a
        
        library(cowplot)
        plot_grid(p.trf.vol[[1]],p.trf.vol[[2]],p.trf.vol[[3]],p.trf.vol[[4]],
                  p.trf.vol[[5]],p.trf.vol[[6]],p.trf.vol[[7]],p.trf.vol[[8]],
                  ncol = 4,
                  label_size = 20,
                  labels = "AUTO")->a
        
        # ggsave(a,filename = "trf.vol.patch.pdf",path = "C:\\Users\\colet\\Desktop\\re snc atlas",units = "cm",width = 60,height = 32)
        
        ##火山图
        
        
        dds<-DESeqDataSetFromMatrix(countData = data.frame(trf.df.all), colData = df.col.1, design = ~tissue+condition)
        
        
        
        # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
        dds = DESeq(dds,test = "LRT", reduced = ~tissue, full = ~tissue+condition)
        res = results(dds, pAdjustMethod = "BH")
        res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
        
        
        
        
        
        res.vol = na.omit(data.frame(results(dds, pAdjustMethod = "BH")))
        library(data.table)
        res.vol.1 <- na.omit(rownames_to_column(res.vol,var = "id"))
        ifelse(res.vol.1$log2FoldChange > .5,"Hypoxia",
               ifelse(res.vol.1$log2FoldChange< -.5,"Normoxia",
                      "Both")) -> res.vol.1$exin
        ifelse(res.vol.1$padj < .01,1.1,1) -> res.vol.1$size
        
        res.vol.1%>%filter(abs(log2FoldChange)>.5&-log10(padj)>1)->id.vol
        
        library(ggrepel)
        ggplot(data = res.vol.1,aes(x=log2FoldChange,y=-log10(padj)))+
          geom_point(aes(size=-log10(padj), color=exin))+
          scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
          scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
          guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "log10(FDR)"))+
          geom_hline(yintercept = log10(10),linetype=5)+
          geom_vline(xintercept = .5,color="#f0a1a8")+
          geom_vline(xintercept = -.5,color="#4994c4")+
          theme_minimal()+
          geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
          labs(title = "tRNA")+
          theme(legend.justification=c(1,0), legend.position=c(1,.5),plot.title = element_text(size= 15, face = 2,hjust = .5),legend.title = element_text(face = 2))+
          # theme(legend.position = legend.justification=c(0,0))+
          scale_y_continuous(limits = c(0,13.5))+
          scale_x_continuous(limits = c(-2,2))->vol.trf
        vol.trf
        
        ggplot(data = res.vol.1%>%filter(baseMean>10),aes(x=log2(baseMean),y=log2FoldChange))+
          geom_point(aes(size=-log10(padj), color=exin))+
          scale_color_manual(values = c("grey70","#f0a1a8","#4994c4"))+
          guides(color=guide_legend(reverse = T,title = "Expressed in"),size=guide_legend(title = "FDR"))+
          geom_text_repel(data = id.vol%>%filter(baseMean>10),aes(size=8,x=log2(baseMean),y=log2FoldChange,label = id))+
          theme_light()+
          geom_hline(yintercept = .5,color="#f0a1a8")+
          geom_hline(yintercept = -.5,color="#4994c4")+
          scale_size_continuous(range = c(0,2.5),breaks = seq(0,7,2.5))+
          theme(legend.position = "none",
                # axis.text = element_text(size = 10),
                # axis.title = element_text(size=10),
                # axis.line = element_line(color = "black"),
                # plot.background = element_rect(color="black",size=1),
                panel.background = element_rect(color = "white"),
                panel.grid = element_blank())->vol.2.trf
        vol.2.trf
        
        vol.trf+
          # theme(legend.position = "none")+
          annotation_custom(
            grob = ggplotGrob(vol.2.trf),
            xmin = -2.3,
            xmax = -.10,
            ymin = 7,
            ymax = 14.25
          )->uni.vol.trf
        
        uni.vol.trf
      }##火山图
      
      
      
      trfdf.sig <- trf.df[rownames(trf.df)%in%rownames(res.df[abs(res.df$log2FoldChange)>0,]),]
      trfdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        trfdf.sig.mean.1  <-  rowMeans(matrix(trfdf.sig[,c(i:i+2)]))
        trfdf.sig.mean <- cbind(trfdf.sig.mean,trfdf.sig.mean.1)
      } 
      colnames(trfdf.sig.mean) <- colnames(trfdf.sig)[seq(1,48,3)]
      # TSI(trfdf.sig)
      trfdf.sig$id <- rownames(trfdf.sig)
      
      library(data.table)
      trfdf.sig.long <- melt(trfdf.sig,id="id")
      name.idx <- strsplit(as.character(trfdf.sig.long$variable),"_")
      for (i in 1:nrow(trfdf.sig.long)) {
        trfdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        trfdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      trfdf.sig.long$condition[trfdf.sig.long$condition=="NOR"] <- "Normoxia"
      trfdf.sig.long$condition[trfdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      trfdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))
      
      ##两个柱子挨在一起的图
      if (T) {
        library(ggplot2)
        # trfdf.sig.long%>%filter(value>0)%>%filter(!grepl('let', id))%>%
        trfdf.sig.long%>%filter(value>45&!grepl('let|3068|3084|23|26|30|342|21-',id))%>%
          ggplot(aes(x=id,y=log2(value)))+
          geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=condition))+
          stat_summary(fun = mean,
                       geom = "errorbar",
                       fun.max = function(x) mean(x) + sd(x),
                       fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
                       aes(group=condition,color=condition),
                       size=.5,
                       width=.25,
          )+
          # geom_vline(xintercept = .5)+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
          facet_wrap(~tissue,ncol=1,strip.position = "right")+
          labs(y="log2(CPM)",x=NULL,title = "tRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "2barEXP.trf.pdf",units = "cm",width = 16,height = 20,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子挨在一起的图
      
      ##两个柱子堆积成1
      if (T) {
        library(ggplot2)
        # trfdf.sig.long%>%filter(value>0)%>%
        trfdf.sig.long%>%filter(value>0&!grepl('let|3068|3084|23|26|30|342|21-',id))%>%
          ggplot(aes(y=id,x=log10(value)))+
          
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "tRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          geom_vline(xintercept = .5,linetype = 5)+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "fillbarEXP.trf.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (T) {
        library(ggplot2)
        # trfdf.sig.long%>%filter(value>0)%>%
        trfdf.sig.long%>%filter(value>0&!grepl('let|3068|3084|23|26|30|342|21-',id))%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "tRNA")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "stackbarEXP.trf.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
      
      # # Create a tibble for LRT results
      # res_LRT_tb <- res %>%
      #   data.frame() %>%
      #   rownames_to_column(var="gene") %>% 
      #   as_tibble()
      # # Subset to return genes with padj < 0.05
      # sigLRT_genes <- res_LRT_tb %>% 
      #   filter(padj < .05)
      # 
      # # Get number of significant genes
      # nrow(sigLRT_genes)
      # 
      # # Subset results for faster cluster finding (for classroom demo purposes)
      # clustering_sig_genes <- sigLRT_genes %>%
      #   arrange(padj) %>%
      #   head(n=1000)
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
      # library(DEGreport)
      # # install.packages('C:/Users/colet/Desktop/lasso2_1.2-22.tar.gz', repo=NULL)
      # clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
      # 
      # # Obtain rlog values for those significant genes
      # cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
      # 
      # BiocManager::install("DEGreport")
    }
    
    ##0.1all找组织比对(第一个是点图)
    if (T) {
      ##先算TSI
      df.col = data.frame(id = colnames(trf.df.nor.all),dex = c(rep("thymus",3),
                                                                rep("kidney",3),
                                                                rep("brain",3),
                                                                rep("liver",3),
                                                                rep("lung",3),
                                                                rep("intestines",3),
                                                                rep("heart",3),
                                                                rep("spleen",3)
      ))
      dds<-DESeqDataSetFromMatrix(countData = trf.df.nor.all, colData = df.col,design = ~dex)
      ## Prefiltering
      filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
      dds <- dds[!filt,]
      
      ## Perform DESeq2()
      dds = DESeq(dds)
      res = results(dds, pAdjustMethod = "BH")
      
      ##determine tissue specific(calculate TSI)
      # res.tsi <-data.frame(res)
      # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
      # res.sig <- na.omit(data.frame(res.sig))
      # res.sig.trf <- data.frame(merge(data.frame(id=rownames(trf.df.hyp.all),trf.df.nor.all),data.frame(id=rownames(res.sig))
      #                                ,by="id"),row.names = 1)
      
      tpm <- function(x){
        tpmres <- (x[,6:ncol(x)]/x[,5])/colSums(x[,6:ncol(x)]/x[,5])*1e6
        return(tpmres)
      }
      
      library(edgeR)
      # apply variance stabilizing transformation
      # vsted <- cpm(counts(dds))
      v = varianceStabilizingTransformation(dds, blind=FALSE)
      vsted = assay(v)
      
      ##calculate TSI
      tsi.trf.nor <- rep(0,nrow(vsted))
      for (i in c(1,4,7,10,13,16,19,21)) {
        tsi.trf.nor <- data.frame(tsi.trf.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
      }
      tsi.trf.nor <- tsi.trf.nor[,-1]
      colnames(tsi.trf.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
      avgexp.trf.nor <- tsi.trf.nor
      tsi.trf.nor <- data.frame(TSI=TSI(tsi.trf.nor,8))
      
      tsi.sig.id <- rownames_to_column(tsi.trf.nor,var="id")%>%filter(TSI>.85)
      
      ##0.1all找组织比对
      # df.col.a <- data.frame(id=colnames(trf.df.all),tissue=c(rep(c(rep(c("Thymus","Kidney","Brain",
      #                                                              "Liver","Lung","Intestines",
      #                                                              "Heart","Spleen"),each=3)),2)),
      #                        condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),2)))
      df.col.a <- data.frame(id=colnames(trf.df),tissue=c(rep(c("Thymus","Kidney","Brain",
                                                                "Liver","Lung","Intestines",
                                                                "Heart","Spleen"),6)),
                             condition=c(rep(rep(c("Normoxia","Hypoxia"),each=12),1)))
      dds<-DESeqDataSetFromMatrix(countData = trf.df, colData = df.col.a,design = ~tissue+condition)
      # dds = DESeq(dds,test = "LRT",reduced = ~condition,full = ~tissue+condition)
      dds = DESeq(dds,test = "LRT", reduced = ~condition, full = ~tissue+condition)
      res = results(dds, pAdjustMethod = "BH")
      res.df <- data.frame(res)%>%filter(baseMean>3&padj<.05)
      
      # res.df.1 <- res.df
      # res.df.1$id <- rownames(res.df.1)
      
      trfdf.sig <- trf.df[rownames(trf.df)%in%rownames(res.df[abs(res.df$log2FoldChange)>0,]),]
      trfdf.sig <- merge(rownames_to_column(trfdf.sig,"id"),tsi.sig.id,by="id")
      trfdf.sig <- column_to_rownames(trfdf.sig,"id")
      # rownames(trfdf.sig) <- trfdf.sig$id
      trfdf.sig.mean = NULL
      for (i in seq(1,48,3)) {
        trfdf.sig.mean.1  <-  rowMeans(matrix(trfdf.sig[,c(i:i+2)]))
        trfdf.sig.mean <- cbind(trfdf.sig.mean,trfdf.sig.mean.1)
      } 
      colnames(trfdf.sig.mean) <- colnames(trfdf.sig)[seq(1,48,3)]
      # TSI(trfdf.sig)
      trfdf.sig$id <- rownames(trfdf.sig)
      
      library(data.table)
      trfdf.sig.long <- melt(trfdf.sig,id=c("id","TSI"))
      name.idx <- strsplit(as.character(trfdf.sig.long$variable),"_")
      for (i in 1:nrow(trfdf.sig.long)) {
        trfdf.sig.long$tissue[i] <- str_to_title(name.idx[[i]][1])
        trfdf.sig.long$condition[i] <- name.idx[[i]][2]
      }
      trfdf.sig.long$condition[trfdf.sig.long$condition=="NOR"] <- "Normoxia"
      trfdf.sig.long$condition[trfdf.sig.long$condition=="HYP"] <- "Hypoxia"
      
      
      trfdf.sig.long%>%dplyr::select(-c(variable))%>%group_by(id,tissue,condition)%>%
        summarise(value_mean=mean(value),sd=sd(value),se=sd(value)/sqrt(n()))->p.dat
      trf.dot.dat <- data.frame(merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id"),RNA="tRNA")
      rbind(rny.dot.dat,trf.dot.dat)
      # dcast(data = p.dat,id~condition+tissue+value_mean+sd+se) 
      ##两个柱子挨在一起的图
      if (T) {
        library(ggplot2)
        # trfdf.sig.long%>%filter(value>0&condition=="Normoxia")%>%
        # p.dat%>%filter(value>0&condition=="Normoxia")%>%
        p.dat.log  <- merge(p.dat,data.frame(rownames_to_column(res.df,"id")),by="id")%>%
          filter(condition=="Normoxia")%>%
          # log10(p.dat.log$value_mean)
          # rbind(rny.dot.dat,trf.dot.dat)%>%
          #   filter(condition=="Normoxia")%>%
          ggplot(aes(y=id,x=tissue))+
          geom_point(aes(size=log(value_mean, 10),color=padj))+
          scale_size_continuous(range = c(-1,8),breaks = seq(.5,3,.5))+
          # facet_grid(~RNA)+
          # facet_wrap(~RNA,ncol = 1)+.
          # geom_point(fun = "mean", stat = "summary",aes(group=tissue,fill=condition,color=condition,size=log(value)),width = .65,position = "dodge")+
          # stat_summary(fun = "mean", geom = "point",aes(group=id,fill=condition,color=condition,size=log(value)))+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
        # )+
        # geom_curve()
        # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
        # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "dodge")+
        # facet_wrap(~tissue,nrow=1,strip.position = "top")+
        labs(y=NULL,x=NULL,title = "tRNA")+
          # scale_color_gradient(low = "blue",high = "red")+
          viridis::scale_color_viridis()+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          theme_gray()+
          guides(size=guide_legend(reverse = F,title = "log10(CPM)"),
                 # color=guide_legend(title = "p.adj")
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(),
                axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
          )->p.trf.tsi.dot
        p.trf.tsi.dot
        
        options(download.file.method = 'libcurl')
        options(url.method='libcurl')
        # BiocManager::install("lemon")
        library(lemon)
        #####
        #拼点图 TSI
        # patchwork::wrap_plots(guides = "collect",p.rny.tsi.dot, p.trf.tsi.dot)
        # grid_arrange_shared_legend(, p.trf.tsi.dot, p.pir.tsi.dot,nrow = 3,ncol = 1,position='right')
        
        p1=p.rny.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(rny.dot.dat$id))/1.8,"cm"))
        p2=p.mir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(mir.dot.dat$id))/1.8,"cm"))
        p3=p.pir.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(pir.dot.dat$id))/1.8,"cm"))
        p4=p.lnc.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(lnc.dot.dat$id))/1.8,"cm"))
        p5=p.trf.tsi.dot%>%egg::set_panel_size(width = unit(8,"cm"),height = unit(length(unique(trf.dot.dat$id))/1.8,"cm"))
        
        p.trf.tsi.dot
        
        
        # p <- grid_arrange_shared_legend(p1,p2,nrow = 2,ncol = 1,position='right')
        
        p <- patchwork::wrap_plots(guides = "collect",p1,p2,p3,p4,p5,ncol=5)
        
        ggsave(plot = p,filename = "test.pdf",units = "cm",width = 80,height = 450,path = "C:\\Users\\colet\\Desktop\\re snc atlas",limitsize = F)
        
        #拼点图 TSI
        #####      
        # ggsave(filename = "TIS.dot.EXP.trf.pdf",units = "cm",width = 16,height = 55,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子挨在一起的图
      
      ##两个柱子堆积成1
      if (F) {
        library(ggplot2)
        trfdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "fill")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="Relative Expression (scaled)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.fillbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积成1
      
      ##两个柱子堆积
      if (F) {
        library(ggplot2)
        trfdf.sig.long%>%filter(value>0)%>%
          ggplot(aes(y=id,x=log10(value)))+
          # geom_bar(fun = "mean", stat = "summary",aes(fill=condition,color=condition),width = .65,position = "stack")+
          stat_summary(fun = "mean", geom = "bar",aes(group=condition,fill=condition,color=condition),width=1,position = "stack")+
          # stat_summary(fun = mean,
          #              geom = "errorbar",
          #              fun.max = function(x) mean(x) + sd(x),
          #              fun.min = function(x) mean(x) - sd(x),position = position_dodge(.65),
          #              aes(group=condition,color=condition),
          #              size=.5,
          #              width=.25,
          # )+
          # geom_curve()
          # stat_summary(fun = "mean", geom = "line",aes(group=condition),color="white",alpha=50,linetype=3) +
          # geom_bar(stat = "identity",aes(fill=condition,color=condition),position = "fill")+
        facet_wrap(~tissue,nrow =1,strip.position = "top")+
          labs(y=NULL,x="log10(CPM)",title = "Rny")+
          # geom_hline(yintercept = .5,linetype=2,color="grey100")+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          scale_fill_manual(values = c("#f0a1a8","#4994c4"))+
          # scale_color_manual(values = c("#f0a1a8","#4994c4"))+
          scale_color_manual(values = c("white","white"))+
          theme_gray()+
          guides(fill=guide_legend(reverse = T,title = NULL),
                 color=F
          )+
          # scale_fill_manual(values = c("#f0c9cf","#a2d2e2"))+
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 8,angle = 90,vjust = .5)
          )
        # stat_summary(fun = mean,
        #              geom = "pointrange",
        #              fun.max = function(x) mean(x) + sd(x),
        #              fun.min = function(x) mean(x) - sd(x),position = "dodge")
        # ggsci::scale_fill_aaas()
        # theme(line = element_line(linetype = 1))
        # ggsave(filename = "TIS.stackbarEXP.rny.pdf",units = "cm",width = 24,height = 12,path = "C:\\Users\\colet\\Desktop\\re snc atlas")
      }##两个柱子堆积
    }
    
  }##新加入新分析
  
  
  
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(trf.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = trf.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  vsted.nor = data.frame(assay(v))
  
  ##calculate TSI
  tsi.trf.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.trf.nor <- data.frame(tsi.trf.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.trf.nor <- tsi.trf.nor[,-1]
  colnames(tsi.trf.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.trf.nor <- tsi.trf.nor
  tsi.trf.nor <- data.frame(tsi.trf.nor=TSI(tsi.trf.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.trf.nor
  p.pca.trf.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(trf.df.hyp.all),dex = c(rep("thymus",3),
                                                              rep("kidney",3),
                                                              rep("brain",3),
                                                              rep("liver",3),
                                                              rep("lung",3),
                                                              rep("intestines",3),
                                                              rep("heart",3),
                                                              rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = trf.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  vsted.hyp = data.frame(assay(v))
  
  ##calculate TSI
  tsi.trf.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.trf.hyp <- data.frame(tsi.trf.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.trf.hyp <- tsi.trf.hyp[,-1]
  colnames(tsi.trf.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.trf.hyp <- tsi.trf.hyp
  tsi.trf.hyp <- data.frame(tsi.trf.hyp=TSI(tsi.trf.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.trf.hyp
  p.pca.trf.hyp
  
  # p.pca.trf.nor|p.pca.trf.hyp -> p.pca.trf.nor.hyp.unite
  
  ##determine tissue specific trfRNA
  tsi.trf.nor  <- data.frame(id=rownames(tsi.trf.nor),tsi.trf.nor,row.names = NULL)
  tsi.trf.nor.spc <- tsi.trf.nor[tsi.trf.nor$tsi.trf.nor>0.85,]
  tsi.trf.hyp  <- data.frame(id=rownames(tsi.trf.hyp),tsi.trf.hyp,row.names = NULL)
  tsi.trf.hyp.spc <- tsi.trf.hyp[tsi.trf.hyp$tsi.trf.hyp>0.85,]
  
  ##determine expression of tissue specific trfRNA
  library(data.table)
  library(ggplot2)
  ##trf对照组
  tsi.trf.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.trf.nor),avgexp.trf.nor),tsi.trf.nor.spc,by="id")
  tsi.trf.nor.spc.long <- melt(tsi.trf.nor.spc.long.1,id=c("id","tsi.trf.nor"))
  colnames(tsi.trf.nor.spc.long) <- c("ID","TSI","Tissue","Exp")
  tsi.trf.nor.spc.long <- tsi.trf.nor.spc.long[order(tsi.trf.nor.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.trf.nor.spc.long,aes(x=Tissue,y=ID,size=Exp,color=TSI))+
    viridis::scale_color_viridis()->p.trf.nor.spc 
  p.trf.nor.spc 
  
  ##trf处理组
  tsi.trf.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.trf.hyp),avgexp.trf.hyp),tsi.trf.hyp.spc,by="id")
  tsi.trf.hyp.spc.long <- melt(tsi.trf.hyp.spc.long.1,id=c("id","tsi.trf.hyp"))
  colnames(tsi.trf.hyp.spc.long) <- c("ID","TSI","Tissue","Exp")
  tsi.trf.hyp.spc.long <- tsi.trf.hyp.spc.long[order(tsi.trf.hyp.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.trf.hyp.spc.long,aes(x=Tissue,y=ID,size=Exp,color=TSI))+
    viridis::scale_color_viridis()->p.trf.hyp.spc 
  p.trf.hyp.spc 
  
  #####
  ##trf热图
  
  #提取氨基酸名取并取唯一
  trf.aa.name.nor <- unique(substr(rownames(avgexp.trf.nor),str_length(rownames(avgexp.trf.nor))-5,str_length(rownames(avgexp.trf.nor))-3))
  trf.aa.name.hyp <- unique(substr(rownames(avgexp.trf.hyp),str_length(rownames(avgexp.trf.hyp))-5,str_length(rownames(avgexp.trf.hyp))-3))
  #给每行命名氨基酸名
  avgexp.trf.nor$aa <- substr(rownames(avgexp.trf.nor),str_length(rownames(avgexp.trf.nor))-5,str_length(rownames(avgexp.trf.nor))-3)
  avgexp.trf.hyp$aa <- substr(rownames(avgexp.trf.hyp),str_length(rownames(avgexp.trf.hyp))-5,str_length(rownames(avgexp.trf.hyp))-3)
  #trf.aa.name.order <- avgexp.trf.nor[order(avgexp.trf.nor$aa),]
  # trf.aa.merge <- left_join(aa,avgexp.trf.nor)
  #整理数据;对照组
  trf.aa.df = NULL
  for (i in trf.aa.name.nor) {
    trf.aa.df.1 <- t(data.frame(log10(colSums(avgexp.trf.nor[avgexp.trf.nor$aa==i,][,-ncol(avgexp.trf.nor)]))))
    trf.aa.df <- rbind(trf.aa.df,trf.aa.df.1)
  }
  rownames(trf.aa.df) <- trf.aa.name.nor
  trf.aa.df.nor <- trf.aa.df
  
  ggheatmap::ggheatmap(trf.aa.df.nor[rownames(trf.aa.df.nor)!="det",],scale = "row")+
    viridis::scale_fill_viridis()+
    ggprism::theme_prism()+
    theme(axis.title= element_text(size = 0),
          axis.text.x = element_text(size = 12,angle = 45,hjust = 1,vjust = 1),
          legend.position = "left",
          legend.key.size = unit(1,"lines")
    )+
    labs(fill = "log10(CPM)") -> p.trf.aa.hm.nor
  
  # trf.aa.df <- data.frame(aa=rownames(trf.aa.df),trf.aa.df)
  
  #整理数据;处理组
  trf.aa.df = NULL
  for (i in trf.aa.name.hyp) {
    trf.aa.df.1 <- t(data.frame(log10(colSums(avgexp.trf.hyp[avgexp.trf.hyp$aa==i,][,-ncol(avgexp.trf.hyp)]))))
    trf.aa.df <- rbind(trf.aa.df,trf.aa.df.1)
  }
  rownames(trf.aa.df) <- trf.aa.name.hyp
  trf.aa.df.hyp <- trf.aa.df
  
  ggheatmap::ggheatmap(trf.aa.df.hyp[rownames(trf.aa.df.hyp)!="det",],scale = "row")+
    viridis::scale_fill_viridis()+
    ggprism::theme_prism()+
    theme(axis.title= element_text(size = 0),
          axis.text.x = element_text(size = 12,angle = 45,hjust = 1,vjust = 1),
          legend.position = "left",
          legend.key.size = unit(1,"lines")
    )+
    labs(fill = "log10(CPM)") -> p.trf.aa.hm.hyp
  
  ##拼图
  p.trf.aa.hm.unite <- cowplot::plot_grid(p.trf.aa.hm.nor+labs(title = "Control"),p.trf.aa.hm.hyp+labs(title = "Hypoxia"),align = "h",labels = c("A","B"))
  # ggsave(filename = "p.hm.trf.aa.unite.pdf", p.trf.aa.hm.unite,units = "cm", width = 20,height = 15)
  
}##tRFs-1

if (T) {
  library(DESeq2)
  library(ggplot2)
  setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  trf.in <- read.table("tRNA_all.quantified",header = T)
  
  trf.in.b <- data.frame(trf.in[,7:ncol(trf.in)],row.names = trf.in$Geneid)
  
  i=1
  for (i in 1:ncol(trf.in.b)) {
    colnames(trf.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(trf.in.b),"_")[[i]][1],strsplit(colnames(trf.in.b),"_")[[i]][2],sep = "_")]
  }
  
  trf.df <- trf.in.b
  ##tRNA处理成对照与处理组分开的表格
  trf.df.all <- data.frame(rep(0,nrow(trf.df)))
  trf.df.nor.all = data.frame(rep(0,nrow(trf.df)))
  trf.df.hyp.all = data.frame(rep(0,nrow(trf.df)))
  for (i in 1:8) {
    trf.df.nor.all <- cbind(trf.df.nor.all,trf.df[,c(i+8,i,i+16)])
  }
  
  for (i in 25:32) {
    trf.df.hyp.all <- cbind(trf.df.hyp.all,trf.df[,c(i+8,i,i+16)])
  }
  trf.df.nor.all <- trf.df.nor.all[,-1]
  trf.df.hyp.all <- trf.df.hyp.all[,-1]
  
  ##trf差异分析（对照对对照、处理对处理；组织别）
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(trf.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = trf.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.trf <- data.frame(merge(data.frame(id=rownames(trf.df.hyp.all),trf.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.trf.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.trf.nor <- data.frame(tsi.trf.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.trf.nor <- tsi.trf.nor[,-1]
  colnames(tsi.trf.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.trf.nor <- tsi.trf.nor
  tsi.trf.nor <- data.frame(tsi.trf.nor=TSI(tsi.trf.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.trf.nor
  p.pca.trf.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(trf.df.hyp.all),dex = c(rep("thymus",3),
                                                              rep("kidney",3),
                                                              rep("brain",3),
                                                              rep("liver",3),
                                                              rep("lung",3),
                                                              rep("intestines",3),
                                                              rep("heart",3),
                                                              rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = trf.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.trf.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.trf.hyp <- data.frame(tsi.trf.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.trf.hyp <- tsi.trf.hyp[,-1]
  colnames(tsi.trf.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.trf.hyp <- tsi.trf.hyp
  tsi.trf.hyp <- data.frame(tsi.trf.hyp=TSI(tsi.trf.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.trf.hyp
  p.pca.trf.hyp
  
  # p.pca.trf.nor|p.pca.trf.hyp -> p.pca.trf.nor.hyp.unite
  
  ##determine tissue specific trfRNA
  tsi.trf.nor  <- data.frame(id=rownames(tsi.trf.nor),tsi.trf.nor,row.names = NULL)
  tsi.trf.nor.spc <- tsi.trf.nor[tsi.trf.nor$tsi.trf.nor>.9,]
  tsi.trf.hyp  <- data.frame(id=rownames(tsi.trf.hyp),tsi.trf.hyp,row.names = NULL)
  tsi.trf.hyp.spc <- tsi.trf.hyp[tsi.trf.hyp$tsi.trf.hyp>.9,]
  
  ##determine expression of tissue specific trfRNA
  library(data.table)
  library(ggplot2)
  ##trf对照组
  tsi.trf.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.trf.nor),avgexp.trf.nor),tsi.trf.nor.spc,by="id")
  tsi.trf.nor.spc.long <- melt(tsi.trf.nor.spc.long.1,id=c("id","tsi.trf.nor"))
  colnames(tsi.trf.nor.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.trf.nor.spc.long <- tsi.trf.nor.spc.long[order(tsi.trf.nor.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.trf.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.trf.nor.spc 
  p.trf.nor.spc 
  
  ##trf处理组
  tsi.trf.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.trf.hyp),avgexp.trf.hyp),tsi.trf.hyp.spc,by="id")
  tsi.trf.hyp.spc.long <- melt(tsi.trf.hyp.spc.long.1,id=c("id","tsi.trf.hyp"))
  colnames(tsi.trf.hyp.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.trf.hyp.spc.long <- tsi.trf.hyp.spc.long[order(tsi.trf.hyp.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.trf.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.trf.hyp.spc 
  p.trf.hyp.spc 
  
  tsi.trf.nor.spc.long <- head(data.frame(tsi.trf.nor.spc.long,condition="Control"),400)
  tsi.trf.hyp.spc.long <- head(data.frame(tsi.trf.hyp.spc.long,condition="Hypoxia"),400)
  tsi.trf.combine.spc.long <- rbind(tsi.trf.nor.spc.long,tsi.trf.hyp.spc.long)
  
  tsi.trf.combine.spc.long[-grep(pattern = "UndetNNN",tsi.trf.combine.spc.long$ID),]%>%
  ggplot()+
    geom_point(aes(x=Tissue,y=ID,size=VST,color=TSI))+
    labs(x="",y="")+
    facet_wrap(~condition)+
    scale_size_continuous(range = c(-1,6))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
    viridis::scale_color_viridis()->p.trf.combine.spc 
  p.trf.combine.spc 
  
  # ggsave(p.trf.combine.spc,filename = "p.trf.combine.spc.pdf",units = "cm",width = 16,height = 32,limitsize = F)
}##tRF-2(PCA与组织特异性)

if (T) {
  library(DESeq2)
  library(ggplot2)
  setwd("C:/Users/colet/Desktop/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  denovomir.in <- read.table("denovo_miR.quantified",header = T)
  
  denovomir.in.b <- data.frame(denovomir.in[,7:ncol(denovomir.in)],row.names = denovomir.in$Geneid)
  
  i=1
  for (i in 1:ncol(denovomir.in.b)) {
    colnames(denovomir.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(denovomir.in.b),"_")[[i]][1],strsplit(colnames(denovomir.in.b),"_")[[i]][2],sep = "_")]
  }
  
  denovomir.df <- denovomir.in.b[rowSums(denovomir.in.b>20),]
  ##denovomirNA处理成对照与处理组分开的表格
  denovomir.df.all <- data.frame(rep(0,nrow(denovomir.df)))
  denovomir.df.nor.all = data.frame(rep(0,nrow(denovomir.df)))
  denovomir.df.hyp.all = data.frame(rep(0,nrow(denovomir.df)))
  for (i in 1:8) {
    denovomir.df.nor.all <- cbind(denovomir.df.nor.all,denovomir.df[,c(i+8,i,i+16)])
  }
  
  for (i in 25:32) {
    denovomir.df.hyp.all <- cbind(denovomir.df.hyp.all,denovomir.df[,c(i+8,i,i+16)])
  }
  denovomir.df.nor.all <- denovomir.df.nor.all[,-1]
  denovomir.df.hyp.all <- denovomir.df.hyp.all[,-1]
  
  ##denovomir差异分析（对照对对照、处理对处理；组织别）
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(denovomir.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = denovomir.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ##determine tissue specific(calculate TSI)
  # res.tsi <-data.frame(res)
  # res.sig$group <- ifelse(res$log2FoldChange>1&res$pvalue<0.05,"Up",ifelse(res$log2FoldChange<1&res$pvalue<0.05,"Down","Not sig"))
  # res.sig <- na.omit(data.frame(res.sig))
  # res.sig.denovomir <- data.frame(merge(data.frame(id=rownames(denovomir.df.hyp.all),denovomir.df.nor.all),data.frame(id=rownames(res.sig))
  #                                ,by="id"),row.names = 1)
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.denovomir.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.denovomir.nor <- data.frame(tsi.denovomir.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.denovomir.nor <- tsi.denovomir.nor[,-1]
  colnames(tsi.denovomir.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.denovomir.nor <- tsi.denovomir.nor
  tsi.denovomir.nor <- data.frame(tsi.denovomir.nor=TSI(tsi.denovomir.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.denovomir.nor
  p.pca.denovomir.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(denovomir.df.hyp.all),dex = c(rep("thymus",3),
                                                              rep("kidney",3),
                                                              rep("brain",3),
                                                              rep("liver",3),
                                                              rep("lung",3),
                                                              rep("intestines",3),
                                                              rep("heart",3),
                                                              rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = denovomir.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind=FALSE)
  vsted = assay(v)
  
  ##calculate TSI
  tsi.denovomir.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.denovomir.hyp <- data.frame(tsi.denovomir.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.denovomir.hyp <- tsi.denovomir.hyp[,-1]
  colnames(tsi.denovomir.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.denovomir.hyp <- tsi.denovomir.hyp
  tsi.denovomir.hyp <- data.frame(tsi.denovomir.hyp=TSI(tsi.denovomir.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.denovomir.hyp
  p.pca.denovomir.hyp
  
  # p.pca.denovomir.nor|p.pca.denovomir.hyp -> p.pca.denovomir.nor.hyp.unite
  
  ##determine tissue specific denovomirRNA
  tsi.denovomir.nor  <- data.frame(id=rownames(tsi.denovomir.nor),tsi.denovomir.nor,row.names = NULL)
  tsi.denovomir.nor.spc <- tsi.denovomir.nor[tsi.denovomir.nor$tsi.denovomir.nor>0.85,]
  
  #加两个p值低的
  tsi.denovomir.nor.spc <- rbind(tsi.denovomir.nor[tsi.denovomir.nor$id%in%c("chr14_14682","chr17_20205"),],tsi.denovomir.nor.spc)
  
  
  # tsi.denovomir.nor[tsi.denovomir.nor$tsi.denovomir.nor<0.5,]
  tsi.denovomir.hyp  <- data.frame(id=rownames(tsi.denovomir.hyp),tsi.denovomir.hyp,row.names = NULL)
  tsi.denovomir.hyp.spc <- tsi.denovomir.hyp[tsi.denovomir.hyp$tsi.denovomir.hyp>0.85,]
  
  #加两个p值低的
  tsi.denovomir.hyp.spc <- rbind(tsi.denovomir.hyp[tsi.denovomir.hyp$id%in%c("chr14_14682","chr17_20205"),],tsi.denovomir.hyp.spc)
  
  ##determine expression of tissue specific denovomirRNA
  library(data.table)
  library(ggplot2)
  ##denovomir对照组
  tsi.denovomir.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.denovomir.nor),avgexp.denovomir.nor),tsi.denovomir.nor.spc,by="id")
  tsi.denovomir.nor.spc.long <- melt(tsi.denovomir.nor.spc.long.1,id=c("id","tsi.denovomir.nor"))
  colnames(tsi.denovomir.nor.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.denovomir.nor.spc.long <- tsi.denovomir.nor.spc.long[order(tsi.denovomir.nor.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.denovomir.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.denovomir.nor.spc 
  p.denovomir.nor.spc 
  
  ##denovomir处理组
  tsi.denovomir.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.denovomir.hyp),avgexp.denovomir.hyp),tsi.denovomir.hyp.spc,by="id")
  tsi.denovomir.hyp.unspc.long.1 <- merge(data.frame(id=rownames(avgexp.denovomir.hyp),avgexp.denovomir.hyp),tsi.denovomir.hyp[tsi.denovomir.hyp$id=="chr17_20207",],by="id")
  tsi.denovomir.hyp.spc.long.1 <- rbind(tsi.denovomir.hyp.spc.long.1,tsi.denovomir.hyp.unspc.long.1)
  ##转长表
  tsi.denovomir.hyp.spc.long <- melt(tsi.denovomir.hyp.spc.long.1,id=c("id","tsi.denovomir.hyp"))
  colnames(tsi.denovomir.hyp.spc.long) <- c("ID","TSI","Tissue","VST")

  tsi.denovomir.hyp.spc.long <- tsi.denovomir.hyp.spc.long[order(tsi.denovomir.hyp.spc.long$TSI,decreasing = T),]
  
  ggplot()+
    geom_point(data=tsi.denovomir.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.denovomir.hyp.spc 
  p.denovomir.hyp.spc 
  
  tsi.denovomir.nor.spc.long <- data.frame(tsi.denovomir.nor.spc.long,condition="Control")
  tsi.denovomir.hyp.spc.long <- data.frame(tsi.denovomir.hyp.spc.long,condition="Hypoxia")
  tsi.denovomir.combine.spc.long <- rbind(tsi.denovomir.nor.spc.long,tsi.denovomir.hyp.spc.long)
  tsi.denovomir.combine.spc.long <- tsi.denovomir.combine.spc.long[-grep("[.]",tsi.denovomir.combine.spc.long$ID),]
  ggplot()+
    geom_bar(data=tsi.denovomir.combine.spc.long,aes(y=abs(VST),x=Tissue,fill=TSI),stat = "identity",width = .5)+
    # facet_wrap(~condition+ID)+
    labs(y="DESeq2 Normalized Counts")+
    facet_grid(vars(ID),vars(condition))+
    scale_fill_viridis_b()+
    theme_light()+
    theme(axis.text.x = element_text(angle=45,hjust = 1))+
    # scale_size_continuous(range=c(-1,6))+
    viridis::scale_color_viridis()->p.denovomir.combine.spc 
  p.denovomir.combine.spc 
  # ggsave(filename = "p.denovomir.combine.spc.pdf",units = "cm",width = 16,height = 16)
  
  # ggsave(p.denovomir.combine.spc,filename = "p.2.denovomir.combine.spc.pdf",units = "cm",width = 16,height = 300,limitsize = F)
  
  ##画饼图
  denovomir.n <- nrow(denovomir.in.b)
  mir.n<- nrow(mir.in.b)
  df.pie <- data.frame(Type=c("Novel miRNA","Annotated miRNA"),percentage=c(paste0(round(denovomir.n/sum(denovomir.n,mir.n)*100,1),"%"),paste0(round(mir.n/sum(denovomir.n,mir.n)*100,1),"%")),n=c(denovomir.n/sum(denovomir.n,mir.n),mir.n/sum(denovomir.n,mir.n)),raw=c(mir.n,denovomir.n))
  library(ggplot2)
  library(ggrepel)
  ggplot(df.pie, aes(x = "", y = n, fill = Type)) +
    geom_col(color = "black") +
    geom_label(aes(label = percentage), color = c("black", "black"),
               position = position_stack(vjust = 0.5),
               show.legend = FALSE) +
    guides(fill = guide_legend(title = "miRNA")) +
    # scale_fill_viridis_d() +
    scale_fill_manual(values = c("#87CFEA", "#FFC0CB"))+
    # scale_fill_brewer(palette = "Pastel1")+ 
    # geom_label_repel(data = df.pie,
    #                  aes(y = n, label = percentage),
    #                  size = 4.5, nudge_x = 1, show.legend = FALSE, nudge_y = 2)+
    coord_polar(theta = "y") + 
    theme_void()
  
  library(plotrix)
  data <- c(19, 21, 54, 12, 36, 12)
  
  c(df.pie$raw,5)
  pie3D(c(df.pie$raw,5), mar = rep(1.75, 4),
        col = c("#87CFEA", "#FFC0CB","grey50"),
        labels = df.pie$percentage,
        explode =.1)
  }##denovo miR

if (T) {
  ##RNA的tsne3d
  
  #miRNA
  {
  df.all<- rbind(mirdf.sig)
  
  tsne.dat <- t(df.all%>%dplyr::select(-c("TSI","id")))
  
  # tsne <- data.frame(t(vsted[rownames(vsted)%in%rownames(data.frame(res)%>%filter(pvalue<0.05&log2FoldChange>.5)),]))
  library(Rtsne)
  res.tsne <- Rtsne(as.matrix(tsne.dat),perplexity = 15,dims = 2)
  colnames(res.tsne$Y) <- paste0("tsne",c("1","2"
                                          # "3"
                                          ))
  
  
  sta.tem <- rownames(tsne.dat)%>%strsplit2("_")
  sta.tem[,1] <- sta.tem[,1]%>%str_to_title()
  colnames(sta.tem) <- c("tissue","condition","sample")
  
  
  tsne.dat.res <- data.frame(res.tsne$Y,sta.tem)
  
  library(ggalt)
  ggplot(data = tsne.dat.res,aes(x=tsne1,y=tsne2,color=tissue,shape=condition))+
    geom_point()+
    guides(color=guide_legend(reverse = T,title = "Expressed in"),shape=guide_legend(title = "Condition",reverse = T))+
    theme_minimal()+
    # geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
    labs(title = "miRNA")+
    theme(
          # legend.justification=c(1,0), 
          # legend.position=c(1,.5),
          plot.title = element_text(size= 15, face = 2,hjust = .5),
          legend.title = element_text(face = 2))+
     # stat_ellipse(level = 0.95)
    geom_encircle(aes(), alpha = .2, show.legend = F)+
    xlim(-80,80)+
    ylim(-80,80)+
    scale_color_manual(values = c(viridis::viridis(8)))->p.mir.tsne
    # scale_fill_manual(values = c("#f0a1a8","#4994c4"))
    # scale_fill_manual(values = c("blue","red"))
  p.mir.tsne
  }
  
  #tRNA
  {
    df.all<- rbind(trfdf.sig)
    
    tsne.dat <- t(df.all%>%dplyr::select(-c("TSI","id")))
    
    # tsne <- data.frame(t(vsted[rownames(vsted)%in%rownames(data.frame(res)%>%filter(pvalue<0.05&log2FoldChange>.5)),]))
    library(Rtsne)
    res.tsne <- Rtsne(as.matrix(tsne.dat),perplexity = 15,dims = 2)
    colnames(res.tsne$Y) <- paste0("tsne",c("1","2"
                                            # "3"
    ))
    
    
    sta.tem <- rownames(tsne.dat)%>%strsplit2("_")
    sta.tem[,1] <- sta.tem[,1]%>%str_to_title()
    colnames(sta.tem) <- c("tissue","condition","sample")
    
    
    tsne.dat.res <- data.frame(res.tsne$Y,sta.tem)
    
    library(ggalt)
    ggplot(data = tsne.dat.res,aes(x=tsne1,y=tsne2,color=tissue,shape=condition))+
      geom_point()+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),shape=guide_legend(title = "Condition",reverse = T))+
      theme_minimal()+
      # geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "tRNA")+
      theme(
        # legend.justification=c(1,0), 
        # legend.position=c(1,.5),
        plot.title = element_text(size= 15, face = 2,hjust = .5),
        legend.title = element_text(face = 2))+
      # stat_ellipse(level = 0.95)
      geom_encircle(aes(), alpha = .2, show.legend = F)+
      xlim(-50,50)+
      ylim(-50,50)+
      scale_color_manual(values = c(viridis::viridis(8)))->p.trf.tsne
    p.trf.tsne
    # scale_fill_manual(values = c("#f0a1a8","#4994c4"))
    # scale_fill_manual(values = c("blue","red"))
  }
  
  #lncRNA
  {
    df.all<- rbind(lncdf.sig)
    
    tsne.dat <- t(df.all%>%dplyr::select(-c("TSI","id")))
    
    # tsne <- data.frame(t(vsted[rownames(vsted)%in%rownames(data.frame(res)%>%filter(pvalue<0.05&log2FoldChange>.5)),]))
    library(Rtsne)
    res.tsne <- Rtsne(as.matrix(tsne.dat),perplexity = 15,dims = 2)
    colnames(res.tsne$Y) <- paste0("tsne",c("1","2"
                                            # "3"
    ))
    
    
    sta.tem <- rownames(tsne.dat)%>%strsplit2("_")
    sta.tem[,1] <- sta.tem[,1]%>%str_to_title()
    colnames(sta.tem) <- c("tissue","condition","sample")
    
    
    tsne.dat.res <- data.frame(res.tsne$Y,sta.tem)
    
    library(ggalt)
    ggplot(data = tsne.dat.res,aes(x=tsne1,y=tsne2,color=tissue,shape=condition))+
      geom_point()+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),shape=guide_legend(title = "Condition",reverse = T))+
      theme_minimal()+
      # geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "lncRNA")+
      theme(
        # legend.justification=c(1,0), 
        # legend.position=c(1,.5),
        plot.title = element_text(size= 15, face = 2,hjust = .5),
        legend.title = element_text(face = 2))+
      # stat_ellipse(level = 0.95)
      geom_encircle(aes(), alpha = .2, show.legend = F)+
      xlim(-80,80)+
      ylim(-80,80)+
      scale_color_manual(values = c(viridis::viridis(8)))->p.lnc.tsne
    p.lnc.tsne
    # scale_fill_manual(values = c("#f0a1a8","#4994c4"))
    # scale_fill_manual(values = c("blue","red"))
  }
  
  #piRNA
  {
    df.all<- rbind(pirdf.sig)
    
    tsne.dat <- t(df.all%>%dplyr::select(-c("TSI","id")))
    
    # tsne <- data.frame(t(vsted[rownames(vsted)%in%rownames(data.frame(res)%>%filter(pvalue<0.05&log2FoldChange>.5)),]))
    library(Rtsne)
    res.tsne <- Rtsne(as.matrix(tsne.dat),perplexity = 15,dims = 2)
    colnames(res.tsne$Y) <- paste0("tsne",c("1","2"
                                            # "3"
    ))
    
    
    sta.tem <- rownames(tsne.dat)%>%strsplit2("_")
    sta.tem[,1] <- sta.tem[,1]%>%str_to_title()
    colnames(sta.tem) <- c("tissue","condition","sample")
    
    
    tsne.dat.res <- data.frame(res.tsne$Y,sta.tem)
    
    library(ggalt)
    ggplot(data = tsne.dat.res,aes(x=tsne1,y=tsne2,color=tissue,shape=condition))+
      geom_point()+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),shape=guide_legend(title = "Condition",reverse = T))+
      theme_minimal()+
      # geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "piRNA")+
      theme(
        # legend.justification=c(1,0), 
        # legend.position=c(1,.5),
        plot.title = element_text(size= 15, face = 2,hjust = .5),
        legend.title = element_text(face = 2))+
      # stat_ellipse(level = 0.95)
      geom_encircle(aes(), alpha = .2, show.legend = F)+
      scale_color_manual(values = c(viridis::viridis(8)))->p.pi.tsne
    p.pi.tsne
    # scale_fill_manual(values = c("#f0a1a8","#4994c4"))
    # scale_fill_manual(values = c("blue","red"))
  }
  
  library(patchwork)
  cowplot::plot_grid(p.mir.tsne,p.trf.tsne,p.lnc.tsne,p.pi.tsne,labels = "AUTO")->p.grid.tsne
  p.grid.tsne=patchwork::wrap_plots(guides = "collect",p.mir.tsne,p.trf.tsne,p.lnc.tsne,p.pi.tsne)
  
  ggsave(p.grid.tsne,filename = "tsne.pdf",path="C:\\Users\\colet\\Desktop\\re snc atlas",units = "cm",width = 16,height = 14)
  
  
  # grid_arrange_shared_legend(, p.lnc.tsi.dot, p.pir.tsi.dot,nrow = 3,ncol = 1,position='right')
  
  #yRNA(画不出)
  {
    df.all<- rbind(ysdf.sig)
    
    tsne.dat <- t(df.all%>%dplyr::select(-c("TSI","id")))
    
    # tsne <- data.frame(t(vsted[rownames(vsted)%in%rownames(data.frame(res)%>%filter(pvalue<0.05&log2FoldChange>.5)),]))
    library(Rtsne)
    res.tsne <- Rtsne(unique(tsne.dat),perplexity = 10,dims = 2)
    colnames(res.tsne$Y) <- paste0("tsne",c("1","2"
                                            # "3"
    ))
    
    
    sta.tem <- rownames(tsne.dat)%>%strsplit2("_")
    sta.tem[,1] <- sta.tem[,1]%>%str_to_title()
    colnames(sta.tem) <- c("tissue","condition","sample")
    
    
    tsne.dat.res <- data.frame(res.tsne$Y,sta.tem)
    
    library(ggalt)
    ggplot(data = tsne.dat.res,aes(x=tsne1,y=tsne2,color=tissue,shape=condition))+
      geom_point()+
      guides(color=guide_legend(reverse = T,title = "Expressed in"),shape=guide_legend(title = "Condition",reverse = T))+
      theme_minimal()+
      # geom_text_repel(data = id.vol,aes(x=log2FoldChange,y=-log10(padj),label = id.vol$id),max.overlaps = 100,arrow = arrow(length=unit(0.01, "npc")),force=1.5)+
      labs(title = "ysRNA")+
      theme(
        # legend.justification=c(1,0), 
        # legend.position=c(1,.5),
        plot.title = element_text(size= 15, face = 2,hjust = .5),
        legend.title = element_text(face = 2))+
      # stat_ellipse(level = 0.95)
      geom_encircle(aes(), alpha = .2, show.legend = F)+
      scale_color_manual(values = c(viridis::viridis(8)))->p.ys.tsne
    p.ys.tsne
    # scale_fill_manual(values = c("#f0a1a8","#4994c4"))
    # scale_fill_manual(values = c("blue","red"))
  }
  
    # geom_encircle(aes(group = tissue),expand=0.1,spread=0.5,s_shape= .5)
  
  
  
  df.all <- rbind(pir.df.all,mir.df.all,trf.df.all,lnc.df.all)
  
  
  
  id.tem=unique(colnames(df.all))
  col.dat <- data.frame(id=id.tem,
                        tissue=str_to_title(strsplit2(id.tem,split = "_")[,1]),
                        condition=str_to_title(strsplit2(id.tem,split = "_")[,2])
                        )
  
  for (i in 1:nrow(col.dat)) {
    col.dat$condition[i] <- ifelse(col.dat$condition[i]=="Nor","Normoxia","Hypoxia")  
  }
  
  col.dat$condition <- relevel(as.factor(col.dat$condition), ref = "Normoxia")
  
  dds <- DESeqDataSetFromMatrix(countData = df.all, colData = col.dat,design = ~tissue+condition)
  res.1 <- DESeq(dds,test = "Wald",fitType = "mean")
  res.res.1 <- data.frame(results(res.1,lfcThreshold = 1))%>%filter(pvalue<.05)
  
  res.2 <- DESeq(dds,test = "Wald",fitType = "local")
  res.res.2 <- data.frame(results(res.2,lfcThreshold = 1))%>%filter(pvalue<.05)
   
  nrow(res.res.1%>%filter(pvalue<.05))
  nrow(res.res.2%>%filter(pvalue<.05))
  
  # snc.df.nor.all <- rbind(pi.df.nor.all,mir.df.nor.all,trf.df.nor.all,denovomir.df.nor.all,lnc.df.nor.all)
  # snc.df.hyp.all <- rbind(pi.df.hyp.all,mir.df.hyp.all,trf.df.hyp.all,denovomir.df.hyp.all,lnc.df.hyp.all)
  # snc.df.nor.all <- rbind(trf.df.nor.all)
  # snc.df.hyp.all <- rbind(trf.df.hyp.all)
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(snc.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = snc.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v.nor = vst(dds, blind=FALSE)
  vsted = assay(v.nor)
  
  ##calculate TSI
  tsi.snc.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.snc.nor <- data.frame(tsi.snc.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.snc.nor <- tsi.snc.nor[,-1]
  colnames(tsi.snc.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.snc.nor <- tsi.snc.nor
  tsi.snc.nor <- data.frame(tsi.snc.nor=TSI(tsi.snc.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v.nor, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.snc.nor
  p.pca.snc.nor
  
  #2、实验组；组织别
  df.col = data.frame(id = colnames(snc.df.hyp.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = snc.df.hyp.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v.hyp = vst(dds, blind=FALSE)
  vsted = assay(v.hyp)
  
  ##calculate TSI
  tsi.snc.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.snc.hyp <- data.frame(tsi.snc.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.snc.hyp <- tsi.snc.hyp[,-1]
  colnames(tsi.snc.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.snc.hyp <- tsi.snc.hyp
  tsi.snc.hyp <- data.frame(tsi.snc.hyp=TSI(tsi.snc.hyp,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v.hyp, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.snc.hyp
  p.pca.snc.hyp
  
  
  
  
  cowplot::plot_grid(p.pca.snc.nor+ggtitle("Noncoding RNAs Control")+theme_prism(),p.pca.snc.hyp+ggtitle("Noncoding RNAs Hypoxia")+theme_prism(),
                     p.pca.mir.nor+ggtitle("miRNA Control")+theme_prism(),p.pca.mir.hyp+ggtitle("miRNA Hypoxia")+theme_prism(),
                     p.pca.pi.nor+ggtitle("piRNA Control")+theme_prism(),p.pca.pi.hyp+ggtitle("piRNA Hypoxia")+theme_prism(),
                     p.pca.rs.nor+ggtitle("rsRNA Control")+theme_prism(),p.pca.rs.hyp+ggtitle("rsRNA Hypoxia")+theme_prism(),
                     p.pca.rs.nor+ggtitle("ysRNA Control")+theme_prism(),p.pca.rs.hyp+ggtitle("ysRNA Hypoxia")+theme_prism(),
                     p.pca.snr.nor+ggtitle("snRNA Control")+theme_prism(),p.pca.snr.hyp+ggtitle("snRNA Hypoxia")+theme_prism(),
                     p.pca.sno.nor+ggtitle("snoRNA Control")+theme_prism(),p.pca.sno.hyp+ggtitle("snoRNA Hypoxia")+theme_prism(),
                     ncol = 2 , labels = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N"),align = "h"
  )->p.pca.all.rna
  # ggsave("p.pca.all.rna.pdf",units = "cm", width = 25,height = 65)
  snc.df.all.all <- cbind(snc.df.nor.all,snc.df.hyp.all)
  df.col = data.frame(id = colnames(snc.df.all.all),dex = c(rep("thymus-Control",3),
                                                            rep("kidney-Control",3),
                                                            rep("brain-Control",3),
                                                            rep("liver-Control",3),
                                                            rep("lung-Control",3),
                                                            rep("intestines-Control",3),
                                                            rep("heart-Control",3),
                                                            rep("spleen-Control",3),
                                                            rep("thymus-Hypoxia",3),
                                                            rep("kidney-Hypoxia",3),
                                                            rep("brain-Hypoxia",3),
                                                            rep("liver-Hypoxia",3),
                                                            rep("lung-Hypoxia",3),
                                                            rep("intestines-Hypoxia",3),
                                                            rep("heart-Hypoxia",3),
                                                            rep("spleen-Hypoxia",3)))
  df.col = data.frame(id = colnames(snc.df.all.all),dex = c(rep(c("Control","Hypoxia"),each=24)))
  dds<-DESeqDataSetFromMatrix(countData = snc.df.all.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds,fitType = "mean")
  res = results(dds)
  v.all = vst(dds, blind=FALSE)
  vsted = assay(v.all)
  
  tsne <- data.frame(t(trf.df.all))
  
  
  # tsne <- data.frame(t(vsted[rownames(vsted)%in%rownames(data.frame(res)%>%filter(pvalue<0.05&abs(log2FoldChange>1))),]))
  
  tsne <- data.frame(t(vsted[rownames(vsted)%in%rownames(data.frame(res)%>%filter(pvalue<0.05&log2FoldChange>.5)),]))
  library(Rtsne)
  res.tsne <- Rtsne(as.matrix(unique(tsne)),perplexity = 10,dims = 3)
  res.tsne.df <- data.frame(res.tsne$Y,group=c(rep(c("thymus-Control","kidney-Control","brain-Control","liver-Control","lung-Control","intestines-Control","heart-Control","spleen-Control",
                                                     "thymus-Hypoxia","kidney-Hypoxia","brain-Hypoxia","liver-Hypoxia","lung-Hypoxia","intestines-Hypoxia","heart-Hypoxia","spleen-Hypoxia"),each=3)),condition=rep(c("Control","Hypoxia"),each=24))
  
  
  {
    res.tsne.df.1 <- data.frame(res.tsne.df,
                                group.1=c(rep(rep(c(1:8),each=3),2)),
                                condition.1=c(rep(c(1,2),each=24)))
    
    res.tsne.df$condition
    shapes=c(16:18)
    shapes <- shapes[as.numeric( res.tsne.df.1$condition.1)]
    color <- viridis::viridis(8)[as.numeric(res.tsne.df.1$group.1)]
    
    library(scatterplot3d)
    scatterplot3d(res.tsne.df.1[,c(1:3)],pch = shapes,grid = T,box = F,angle = 45)
    addgrids3d(res.tsne.df, grid = c("xy", "xz", "yz"))
    
    # 1. 源函数
    source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
    # 2. 使用 pch="" 清空 3D 散点图
    s3d <- scatterplot3d(res.tsne.df.1[,c(1:3)], pch = "", grid=FALSE, box=FALSE,angle = 60)
    # 3. 添加网格
    addgrids3d(res.tsne.df.1[,c(1:3)], grid = c("xy", "xz", "yz"),angle = 60)
    # 4. 添加点
    s3d$points3d(res.tsne.df.1[,c(1:3)], pch = shapes,col = color)
    # 添加回归平面
    my.lm <- lm(pca.df$PC1 ~ pca.df$PC2 + pca.df$PC3)
    s3d$plane3d(my.lm)
    legend("right", legend = unique(pca.2$tissue),
           col = unique(color),
           pch = 16,
           # inset = -.4,
           xpd = T,
           # bty = "n",
           # bg = "transparent".
           # horiz = T
    )}
  
  
  
  
  
  
  
  
  
  
  library(tsne)
  X <- tsne(as.matrix(unique(tsne)), initial_config = NULL, k = 2, initial_dims = 30, perplexity = 30,
       max_iter = 1000, min_cost = 0, epoch_callback = NULL, whiten = TRUE,
       epoch=100)
  
  data.frame(X,group=c(rep(c("thymus-Control","kidney-Control","brain-Control","liver-Control","lung-Control","intestines-Control","heart-Control","spleen-Control",
                                     "thymus-Hypoxia","kidney-Hypoxia","brain-Hypoxia","liver-Hypoxia","lung-Hypoxia","intestines-Hypoxia","heart-Hypoxia","spleen-Hypoxia"),each=3)),condition=rep(c("Control","Hypoxia"),each=24))
  
  pca.2 <- data.frame(pca.df,group=c(rep(c("thymus-Control","kidney-Control","brain-Control","liver-Control","lung-Control","intestines-Control","heart-Control","spleen-Control",
                             "thymus-Hypoxia","kidney-Hypoxia","brain-Hypoxia","liver-Hypoxia","lung-Hypoxia","intestines-Hypoxia","heart-Hypoxia","spleen-Hypoxia"),each=3)),condition=rep(c("Control","Hypoxia"),each=24),
                      group.1=c(rep(rep(c(1:8),each=3),2)),condition.1=c(rep(c(1,2),each=24)))
  
  for (i in 1:nrow(pca.2)) {
    pca.2$tissue[i] <- strsplit(pca.2$group,"-")[[i]][1]
  }
  
  shapes=c(16:18)
  shapes <- shapes[as.numeric(pca.2$condition.1)]
  color <- viridis::viridis(8)[as.numeric(pca.2$group.1)]
  
  pca <- prcomp(tsne,center = T,scale. = T)
  pca.df <- as.data.frame(pca$x)[,c(1:3)]
  
  
  # install.packages("scatterplot3d")
  
  # library("scatterplot3d")
  
  
  scatterplot3d(pca.df,pch = shapes,grid = T,box = F,angle = 45)
  addgrids3d(pca.df, grid = c("xy", "xz", "yz"))
  
  # 1. 源函数
  source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
  # 2. 使用 pch="" 清空 3D 散点图
  s3d <- scatterplot3d(pca.df, pch = "", grid=FALSE, box=FALSE,angle = 60)
  # 3. 添加网格
  addgrids3d(pca.df, grid = c("xy", "xz", "yz"),angle = 60)
  # 4. 添加点
  s3d$points3d(pca.df, pch = shapes,col = color)
  # 添加回归平面
  my.lm <- lm(pca.df$PC1 ~ pca.df$PC2 + pca.df$PC3)
  s3d$plane3d(my.lm)
  legend("right", legend = unique(pca.2$tissue),
         col = unique(color),
         pch = 16,
         # inset = -.4,
         xpd = T,
         # bty = "n",
         # bg = "transparent".
         # horiz = T
         )

  
  
  
  # install.packages("umap")
  # BiocManager::install("umap")
  library(umap)
              
library(ggforce)
  ggplot(res.tsne.df,aes(X1,X2,group=group))+
    # geom_point(aes(X1,X2,fill = group),shape = 21,size=5,color="white")+
    # scale_fill_manual(values =  c(viridis::viridis(16)))+
    # geom_mark_hull(aes(fill = group, label = group),
    #                   con.cap = 0,label.fill='gray',
    #                   label.colour="black")+
    geom_mark_rect(aes(fill = group, label = group),
                   con.cap = 0,label.fill='gray',
                   label.colour="black",label.fontsize = 10)+
    scale_fill_manual(values =rep(c(viridis::viridis(8)),each=2))+
    geom_point(shape=21,aes(fill=group),colour="black",size=3) +
    theme_minimal()+
    theme(legend.position = "none")+
    labs(x="tSNE 1",y="tSNE 2")->p
  p
  # ggsave("tsne.all.pdf",units = "cm",width = 32,height = 16)
    # scale_fill_nejm() +
    # ggalt::geom_encircle(aes(fill=group), alpha = 0.1, show.legend = F)+
    # theme_prism()+
    # calculate_ellipse(aes(color = group,fill = group),
    #              geom = "polygon",
    #              alpha = 0.3,
    #              linetype = 2)+
    ggforce::geom_mark_ellipse(aes(fill = group,
                                  color = group))+
    theme_prism()+
    theme(legend.position = "topleft")
  
  
  # install.packages("basetheme")
  library(basetheme)
  pars <- basetheme("default")
  # pars$palette <- c("#2A363B", "#019875", "#99B898", "#FECEA8", "#FF847C", "#E84A5F", "#C0392B","#96281B")
  pars$palette <- c(viridis::viridis(5))
  pars$bg  <- "white"
  pars$fg  <- "gray20"
  pars$col <- "gray20"
  pars$col.main <- "black"
  pars$col.axis <- "gray20"
  pars$col.lab  <- "gray20"
  # pars$family   <-  "CascadiaCode-Regular" # 字体
  pars$lab      <-  c(10,10,7)
  pars$las      <-  1
  pars$rect.border <- "black"
  pars$rect.lwd    <- 4
  
  basetheme(pars)
  library(basetheme)
  basetheme("clean")
  cluster::clusplot(res.tsne.df[,1:2], res.tsne.df$group, color=T, shade=F, labels=2, lines=0)
  legend("topleft",
         legend = unique(res.tsne.df$group),
         col = lab2col(unique(res.tsne.df$group)),
         pch = par("pch"),
         cex = .5,
         # horiz = TRUE,
         bty = "n",
         # inset = c(0, 1),
         xpd = T)

    
  
  ## apply variance stabilizing transformation
  v.all = vst(dds, blind=FALSE)
  vsted = assay(v.all)
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v.all, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.snc
  p.pca.snc
}##所有非编码RNA.PCA或tsne

if (T) {
  setwd("C:/Users/colet/Desktop/re snc atlas/")
  all.snc.status <- read.csv("all.snc.status.csv")
  
  for (i in 1:nrow(all.snc.status)) {
    rownames(all.snc.status)[i] <- paste(all.snc.status$tissue[i],all.snc.status$condition[i],sep = "_")
  }
  all.snc.status <- data.frame(sample=rownames(all.snc.status),all.snc.status,row.names = NULL)
  #all.snc.df.t <- data.frame(t(all.snc.status[,c(2:11)]))
  
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(ggprism)
  all.snc.status.all <- melt(all.snc.status,id=c("sample","tissue","condition"))
  
  ##1.堆积柱状图，对照组
  ggplot(data = all.snc.status.all%>%filter(variable!="Other"),aes(x=tissue,y=value,fill=variable))+
    geom_bar(stat = "identity",position = "fill")+
    facet_wrap(~condition)+
    labs(x="",y="Ratio")+
    scale_x_discrete(guide = "prism_bracket")+
    scale_y_continuous(guide = "prism_offset_minor")+
    #  scale_fill_manual(values = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF","#FDE725FF","#B4DE2CFF"))+
    scale_fill_manual(values = c("#98d09d","#d7e698","#dadada","#fbf398","#9b8191","#8f888b","#f7a895","#e77381","#FDE725FF"))+
    theme_prism(axis_text_angle = 90) -> p.ratio.bycondition.palette.1
  #  coord_polar()
  #  ggsci::scale_fill_npg()
  #  viridis::scale_fill_viridis()+
  #  ylim(-1,NA)
  
  ##2.堆积柱状图，处理组
  ggplot(data = all.snc.status.all%>%filter(variable!="Other"),aes(x=tissue,y=value,fill=variable))+
    geom_bar(stat = "identity",position = "fill")+
    facet_wrap(~condition)+
    labs(x="",y="Ratio")+
    scale_x_discrete(guide = "prism_bracket")+
    scale_y_continuous(guide = "prism_offset_minor")+
    scale_fill_manual(values = c(viridis(9)))+
    #  scale_fill_manual(values = c("#98d09d","#d7e698","#dadada","#fbf398","#9b8191","#8f888b","#B4DE2CFF","#FDE725FF","#f7a895","#e77381"))+
    theme_prism(axis_text_angle = 90) -> p.ratio.bycondition.palette.2
  
  # ggsave(filename = "p.ratio.bycondition.palette.1.pdf", p.ratio.bycondition.palette.1, units = "cm", width = 16, height = 10)
  # ggsave(filename = "p.ratio.bycondition.palette.2.pdf", p.ratio.bycondition.palette.2, units = "cm", width = 16, height = 10)
  
  # library(DESeq2)
  # 
  # col <- data.frame(id = colnames(all.snc.df.t), dex = c(rep("Control",8),rep("Hypoxia",8)),stringsAsFactors = T)
  # dds <- DESeqDataSetFromMatrix(countData = all.snc.df.t,colData = col,design = ~dex)
  # dds <- DESeq(dds)
  # res = results(dds, pAdjustMethod = "bonferroni")
  # 
  # v = DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
  # vsted = data.frame(assay(v))
  # 
  # 
  # 
  #   a[[i]] <- cbind(a[[i]][,1],b)
  #   
  #   
  # for (i in 1:48) {
  #   spn <- sample$readable_name[sample$sample_name==gsub("-",".", paste(strsplit(listys[i],"_")[[1]][1],strsplit(listys[i],"_")[[1]][2],sep = "_"))]
  #   colnames(ys.in[[i]]) <- c("id",spn) 
  #   names(ys.in[i]) <- spn
  #   ys.in.b <- merge(ys.in.b,ys.in[[i]],by="id")
  # }
  all.snc.status.all.2 <- all.snc.status.all
  # rownames(all.snc.status.all.2) <- c("Sample","Tissue","")
  ##堆积点图
  ggplot(data = all.snc.status.all%>%filter(variable!="Other"), aes(y= tissue, x=value, color=variable))+
    geom_point(aes(size=value))+
    scale_color_manual(values = c(rev(viridis(n = 9, alpha = 0.5))))+
    labs(x="DESeq2 Normalized Counts",y="Tissue")+
    scale_x_continuous(guide = "prism_offset_minor")+
    scale_y_discrete(guide = "prism_bracket")+
    scale_size_area(max_size = 13)+
    facet_wrap(~condition,nrow = 2)+
    theme_minimal()+
    theme(legend.title = element_text(size = 0)) -> p.all.count.status.1
    p.all.count.status.1
    
    ggplot(data = all.snc.status.all%>%filter(variable!="Other"), aes(y= tissue, x=value, color=variable,alpha = 0.5))+
    geom_point(aes(size=value))+
    scale_color_manual(values = c(rev(viridis(n = 9))))+
    labs(x="DESeq2 Normalized Counts",y="Tissue")+
    scale_x_continuous(guide = "prism_offset_minor")+
    scale_y_discrete(guide = "prism_bracket")+
    scale_size_area(max_size = 13)+
    facet_wrap(~condition,ncol = 2)+
    theme_prism()+
    theme(axis.text.x = element_text(angle = 45,size = 10)) -> p.all.count.status.2
    p.all.count.status.2
    
  # ggsave(p.all.count.status.1,filename = "p.all.count.status.1.pdf",units = "cm", width = 12, height = 16)  
  # ggsave(p.all.count.status.2,filename = "p.all.count.status.2.pdf",units = "cm", width = 16, height = 14)
  
  
  
  
}##统计各种RNA占比

if (T) {
  library(dplyr)
  library(tidyverse)
  setwd("C:/Users/colet/Desktop/re snc atlas/length distribution/")
  ##read files as a large list
  listrs=list.files(pattern = "csv")
  
  ##rename every samples with readable names
  length.in = map(listrs, ~ read.csv(.,header = T))
  
  tissue.name <- c("length","thymus","kidney","brain","liver","lung","intestines","heart","spleen") 
  rna <- c("lncRNA","miRNA","piRNA","scaRNA","snoRNA","snRNA","tRFs")
  
  length.combine=NULL
  
  for (i in c(2:7)) {
    length.nor <- length.in[[i]][,c(1:9)]
    length.hyp <- length.in[[i]][,c(1,10:17)]
    colnames(length.nor) <- tissue.name
    colnames(length.hyp) <- tissue.name
    
    length.nor <- data.frame(length.nor,condition="Control",rna=rna[i])
    single.length.nor <- data.frame(length=unique(length.nor$length))
    length.nor <- left_join(single.length.nor,length.nor)
    
    length.hyp <- data.frame(length.hyp,condition="Hypoxia",rna=rna[i])
    single.length.hyp <- data.frame(length=unique(length.hyp$length))
    length.hyp <- left_join(single.length.hyp,length.hyp)
    
    length.combine <- rbind(length.combine,length.nor,length.hyp)
  }
  
  library(data.table)
  library(ggplot2)
  library(ggalt)
  library(ggbreak)
  library(ggprism)
  library(ggbump)
  library(viridis)
  length.combine.long <- melt(length.combine,id=c("length","condition","rna"))
  
  ggplot(data = length.combine.long%>%filter(length<100&value>20),aes(x=length,y=log10(value),color=rna))+
    #  geom_point()+
    geom_xspline(spline_shape = 0.5,size = 1)+
    coord_cartesian(xlim =c(0, 100),ylim = c(1,6))+
    labs(x="Length",y="log10(VST)")+
    scale_color_manual(values = c("#85BA8F", "#A3C8DC",
                                  "#349839","#EA5D2D",
                                  "#EABB77","#F09594","#2072A8"))+
    # scale_color_manual(values = c(viridis(5)))+
    facet_wrap(~variable+condition,nrow = 2)+
    guides(x="prism_offset_minor",y = "prism_offset_minor") +
    theme_prism()->p.length.distribution
  p.length.distribution
  # ggsave(filename = "p.length.distribution.pdf", p.length.distribution, units = "cm", width = 20,height = 16)
}##统计各种RNA长度分布（旧.匹配后基因的长度）

if (T) {
  setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##read files as a large list
  mir.in <- read.table("miR.quantified",header = T)
  
  mir.in.b <- data.frame(mir.in[,7:ncol(mir.in)],row.names = mir.in$Geneid)
  
  i=1
  for (i in 1:ncol(mir.in.b)) {
    colnames(mir.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(mir.in.b),"_")[[i]][1],strsplit(colnames(mir.in.b),"_")[[i]][2],sep = "_")]
  }
  
  mir.df <- mir.in.b
  ##miRNA处理成对照与处理组分开的表格
  mir.df.all <- data.frame(rep(0,nrow(mir.df)))
  mir.df.nor.all = data.frame(rep(0,nrow(mir.df)))
  mir.df.hyp.all = data.frame(rep(0,nrow(mir.df)))
  for (i in 1:8) {
    mir.df.nor.all <- cbind(mir.df.nor.all,mir.df[,c(i+8,i,i+16)])
  }
  i=9
  for (i in 25:32) {
    mir.df.hyp.all <- cbind(mir.df.hyp.all,mir.df[,c(i+8,i,i+16)])
  }
  mir.df.nor.all <- mir.df.nor.all[,-1]
  mir.df.hyp.all <- mir.df.hyp.all[,-1]
  
  ##mir差异分析（对照对对照、处理对处理；组织别）
  
  #1、对照组；组织别
  df.col = data.frame(id = colnames(mir.df.nor.all),dex = c(rep("thymus",3),
                                                            rep("kidney",3),
                                                            rep("brain",3),
                                                            rep("liver",3),
                                                            rep("lung",3),
                                                            rep("intestines",3),
                                                            rep("heart",3),
                                                            rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = mir.df.nor.all, colData = df.col,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  ## apply variance stabilizing transformation
  v = varianceStabilizingTransformation(dds, blind = FALSE)
  vsted = assay(v)
  vsted.nor = assay(v)
  
  ##calculate TSI
  tsi.mir.nor <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.mir.nor <- data.frame(tsi.mir.nor,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.mir.nor <- tsi.mir.nor[,-1]
  colnames(tsi.mir.nor) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.mir.nor <- tsi.mir.nor
  tsi.mir.nor <- data.frame(tsi.mir.nor=TSI(tsi.mir.nor,8))
  
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.mir.nor
  p.pca.mir.nor
  
  ##2、处理组；组织别
  df.col.2 = data.frame(id = colnames(mir.df.hyp.all),dex = c(rep("thymus",3),
                                                              rep("kidney",3),
                                                              rep("brain",3),
                                                              rep("liver",3),
                                                              rep("lung",3),
                                                              rep("intestines",3),
                                                              rep("heart",3),
                                                              rep("spleen",3)
  ))
  dds<-DESeqDataSetFromMatrix(countData = mir.df.hyp.all, colData = df.col.2,design = ~dex)
  ## Prefiltering
  filt <- rowSums(counts(dds) < 10) > dim(df.col.2)[1]*0.9
  dds <- dds[!filt,]
  
  ## Perform DESeq2()
  dds = DESeq(dds)
  res = results(dds, pAdjustMethod = "bonferroni")
  
  
  ##挑选显著的id
  ids <- rownames(res[data.frame(res)$pvalue<0.05,])
  mir.df.hyp.scale <- data.frame(scale(mir.df.hyp.all))
  mir.df.hyp.scale.1 <- mir.df.hyp.scale[rownames(mir.df.hyp.scale)%in%ids,]
  
  ## apply variance stabilizing transformation
  # v = varianceStabilizingTransformation(dds, blind=FALSE)
  # vsted = assay(v)
  # vsted.hyp = assay(v)
  
  
  vsted=mir.df.hyp.scale.1
  
  ##calculate TSI
  tsi.mir.hyp <- rep(0,nrow(vsted))
  for (i in c(1,4,7,10,13,16,19,21)) {
    tsi.mir.hyp <- data.frame(tsi.mir.hyp,rowMeans(vsted[,c(i,i+1,i+2)]))
  }
  tsi.mir.hyp <- tsi.mir.hyp[,-1]
  colnames(tsi.mir.hyp) <- c("thymus","kidney","brain","liver","lung","intestines","heart","spleen")
  avgexp.mir.hyp <- tsi.mir.hyp
  tsi.mir.hyp <- data.frame(tsi.mir.hyp=TSI(tsi.mir.hyp,8))
  
  ## Plot PCA of VST values
  DESeq2::plotPCA(v, intgroup=c("dex"))+
    ggalt::geom_encircle(aes(fill=dex), alpha = 0.1, show.legend = F)+
    geom_point(size = 3.5)+
    scale_fill_discrete(type = c(viridis::viridis(n=8)))+
    scale_colour_discrete(type = c(viridis::viridis(n=8)))+
    guides(fill = "none")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size = 15),
          axis.text = element_text(size = 11),axis.title = element_text(size = 13),
          legend.text = element_text(size = 11),legend.title = element_text(size = 0)
    )->p.pca.mir.hyp
  p.pca.mir.hyp
  
  #  p.pca.mir.nor|p.pca.mir.hyp -> p.pca.mir.nor.hyp.unite
  
  ##determine tissue specific mirRNA
  tsi.mir.nor  <- data.frame(id=rownames(tsi.mir.nor),tsi.mir.nor,row.names = NULL)
  tsi.mir.nor.spc <- tsi.mir.nor[tsi.mir.nor$tsi.mir.nor>0.85,]
  tsi.mir.hyp  <- data.frame(id=rownames(tsi.mir.hyp),tsi.mir.hyp,row.names = NULL)
  tsi.mir.hyp.spc <- tsi.mir.hyp[tsi.mir.hyp$tsi.mir.hyp>0.85,]
  
  ##在均值数据框中确定reads大于10的
  # avgexp.mir.nor <- avgexp.mir.nor[rowSums(avgexp.mir.nor)>10,]
  # avgexp.mir.hyp <- avgexp.mir.hyp[rowSums(avgexp.mir.hyp)>10,]
  ##determine expression of tissue specific mirRNA
  library(data.table)
  library(ggplot2)
  ##mir对照组
  tsi.mir.nor.spc.long.1 <- merge(data.frame(id=rownames(avgexp.mir.nor),avgexp.mir.nor),tsi.mir.nor.spc,by="id")
  tsi.mir.nor.spc.long <- melt(tsi.mir.nor.spc.long.1,id=c("id","tsi.mir.nor"))
  colnames(tsi.mir.nor.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.mir.nor.spc.long <- tsi.mir.nor.spc.long[order(tsi.mir.nor.spc.long$TSI,decreasing = T),]
  
  #####
  ##找出特异所在组织
  spc.tissue=data.frame(NULL)
  for (i in 1:nrow(tsi.mir.nor.spc.long.1)) {
    spc.tissue <- rbind(data.frame(id=tsi.mir.nor.spc.long.1$id[i],spcin=colnames(tsi.mir.nor.spc.long.1[i,-c(1,10)])[tsi.mir.nor.spc.long.1[i,-c(1,10)]==max(tsi.mir.nor.spc.long.1[i,-(c(1,10))])]),spc.tissue)
  }
  sig.tissue <- table(data.frame(spc.tissue,row.names = 1))
  sig.tissue
  sig.tissue <- data.frame(spc.tissue,row.names = 1)
  write.csv(sig.tissue, file = "sig.tissue.mir.csv", quote = F)
  # tsi.mir.nor.spc.long.2 <- merge(tsi.mir.nor.spc.long.1,spc.tissue,by="id")
  # tsi.mir.nor.spc.long.4 <- melt(tsi.mir.nor.spc.long.2,id=c("id","tsi.mir.nor","spcin"))
  # colnames(tsi.mir.nor.spc.long.4) <- c("ID","TSI","Speciticify","Tissue","VST")
  # tsi.mir.nor.spc.long.4 <- tsi.mir.nor.spc.long.4[order(tsi.mir.nor.spc.long.4$TSI,decreasing = T),]
  # ggplot()+
  #   geom_point(data=tsi.mir.nor.spc.long.4,aes(x=Tissue,y=ID,size=VST,color=TSI))+
  #   viridis::scale_color_viridis()+
  #   geom_point(data=tsi.mir.nor.spc.long.4,aes(x=9,y=ID,color=Speciticify))->p.mir.nor.spc 
  # p.mir.nor.spc 
  
  #####
  ggplot()+
    geom_point(data=tsi.mir.nor.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    viridis::scale_color_viridis()->p.mir.nor.spc 
  p.mir.nor.spc 
  
  ##mir处理组
  tsi.mir.hyp.spc.long.1 <- merge(data.frame(id=rownames(avgexp.mir.hyp),avgexp.mir.hyp),tsi.mir.hyp.spc,by="id")
  tsi.mir.hyp.spc.long <- melt(tsi.mir.hyp.spc.long.1,id=c("id","tsi.mir.hyp"))
  colnames(tsi.mir.hyp.spc.long) <- c("ID","TSI","Tissue","VST")
  tsi.mir.hyp.spc.long <- tsi.mir.hyp.spc.long[order(tsi.mir.hyp.spc.long$TSI,decreasing = T),]
  
  
  ggplot()+
    geom_point(data=tsi.mir.hyp.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    scale_size_continuous(range=c(0,6))+
    viridis::scale_color_viridis()->p.mir.hyp.spc 
  p.mir.hyp.spc 
  
  tsi.mir.nor.spc.long <- data.frame(tsi.mir.nor.spc.long,condition="Control")
  tsi.mir.hyp.spc.long <- data.frame(tsi.mir.hyp.spc.long,condition="Hypoxia")
  tsi.mir.combine.spc.long <- rbind(tsi.mir.nor.spc.long,tsi.mir.hyp.spc.long)
  
  ggplot()+
    geom_point(data=tsi.mir.combine.spc.long,aes(x=Tissue,y=ID,size=VST,color=TSI))+
    facet_wrap(~condition)+
    viridis::scale_color_viridis()->p.mir.combine.spc 
  p.mir.combine.spc
  
  # ggsave(filename = "p.mir.combine.spc.pdf",p.mir.combine.spc,units = "cm",width = 16,height = 150,limitsize = F)
  
  ####################3p/5p特异性
  ##5p/3p特异性
  i=1
  ##获取手臂名称和miR名称
  #1.对照组
  mir.df.vst.nor <- data.frame(vsted.nor)
  for (i in 1:nrow(mir.df.vst.nor)) {
    mir.df.vst.nor$arm[i] <- ifelse(
      strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.nor)[3],"-")[[1]])]=="3p","3p",
      ifelse(strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.nor)[3],"-")[[1]])]=="5p","5p","no"
      )
    )
    mir.df.vst.nor$mirname[i] <- paste(
      strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.nor)[3],"-")[[1]])-1)][1],
      strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.nor)[3],"-")[[1]])-1)][2],
      strsplit(rownames(mir.df.vst.nor)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.nor)[3],"-")[[1]])-1)][3],
      sep = "-")
  }
  
  mir.df.vst.nor <- na.omit(mir.df.vst.nor)
  
  
  mir.df.vst.nor.3p <- mir.df.vst.nor[mir.df.vst.nor$arm=="3p",]
  mir.df.vst.nor.5p <- mir.df.vst.nor[mir.df.vst.nor$arm=="5p",]
  
  
  mir.df.vst.nor[,c(1:24)] <- log2(mir.df.vst.nor[,c(1:24)])
  
  library(data.table)
  mir.df.vst.nor.long <- melt(mir.df.vst.nor,id=c("mirname","arm"))
  colnames(mir.df.vst.nor.long) <- c("miRNA","Arm","Sample","Reads") 
  
  ##
  for (i in 1:nrow(mir.df.vst.nor.long)) {
    mir.df.vst.nor.long$Tissue[i] <- strsplit(as.character(mir.df.vst.nor.long$Sample)[i],"_")[[1]][1]
    mir.df.vst.nor.long$Condition[i] <- strsplit(as.character(mir.df.vst.nor.long$Sample)[i],"_")[[1]][2]
  }
  
  # BiocManager::install("introdataviz")
  # BiocManager::install("CBNplot")
  library(tidyverse)
  # library(introdataviz)
  library(ggunchained)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(ggprism)
  
  mir.df.vst.nor.long.rmarmless <- mir.df.vst.nor.long[mir.df.vst.nor.long$Arm!="no",]
  ggplot(data = mir.df.vst.nor.long.rmarmless%>%filter(Condition=="NOR"),aes(x=Tissue,y=Reads,fill=Arm))+
    geom_split_violin()+
    facet_wrap(~miRNA,ncol = 4)+
    stat_summary(fun.data = "mean_sd", position = position_dodge(0.15), geom = "errorbar", width=.1)+
    stat_summary(fun = "mean", geom = "point", position = position_dodge(0.15), show.legend = F)+
    stat_compare_means(aes(group = mir.df.vst.nor.long.rmarmless$Arm),label = "p.signif", method = "t.test")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1 ))-> p.mir.arm.nor
  #  ggsave("violin.arm.nor.pdf",units = "cm", width = 100,height = 600,limitsize = F)
  
  ##2、处理组
  mir.df.vst.hyp <- data.frame(vsted.hyp)
  for (i in 1:nrow(mir.df.vst.hyp)) {
    mir.df.vst.hyp$arm[i] <- ifelse(
      strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.hyp)[3],"-")[[1]])]=="3p","3p",
      ifelse(strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.hyp)[3],"-")[[1]])]=="5p","5p","no"
      )
    )
    mir.df.vst.hyp$mirname[i] <- paste(
      strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.hyp)[3],"-")[[1]])-1)][1],
      strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.hyp)[3],"-")[[1]])-1)][2],
      strsplit(rownames(mir.df.vst.hyp)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.hyp)[3],"-")[[1]])-1)][3],
      sep = "-")
  }
  
  mir.df.vst.hyp <- na.omit(mir.df.vst.hyp)
  
  
  mir.df.vst.hyp.3p <- mir.df.vst.hyp[mir.df.vst.hyp$arm=="3p",]
  mir.df.vst.hyp.5p <- mir.df.vst.hyp[mir.df.vst.hyp$arm=="5p",]
  
  
  mir.df.vst.hyp[,c(1:24)] <- log2(mir.df.vst.hyp[,c(1:24)])
  
  library(data.table)
  mir.df.vst.hyp.long <- melt(mir.df.vst.hyp,id=c("mirname","arm"))
  colnames(mir.df.vst.hyp.long) <- c("miRNA","Arm","Sample","Reads") 
  
  ##
  for (i in 1:nrow(mir.df.vst.hyp.long)) {
    mir.df.vst.hyp.long$Tissue[i] <- strsplit(as.character(mir.df.vst.hyp.long$Sample)[i],"_")[[1]][1]
    mir.df.vst.hyp.long$Condition[i] <- strsplit(as.character(mir.df.vst.hyp.long$Sample)[i],"_")[[1]][2]
  }
  
  # BiocManager::install("introdataviz")
  # BiocManager::install("CBNplot")
  library(tidyverse)
  # library(introdataviz)
  library(ggunchained)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(ggprism)
  
  mir.df.vst.hyp.long.rmarmless <- mir.df.vst.hyp.long[mir.df.vst.hyp.long$Arm!="no",]
  ggplot(data = mir.df.vst.hyp.long.rmarmless,aes(x=Tissue,y=Reads,fill=Arm))+
    geom_split_violin()+
    facet_wrap(~miRNA,ncol = 4)+
    stat_summary(fun.data = "mean_sd", position = position_dodge(0.15), geom = "errorbar", width=.1)+
    stat_summary(fun = "mean", geom = "point", position = position_dodge(0.15), show.legend = F)+
    stat_compare_means(aes(group = mir.df.vst.hyp.long.rmarmless$Arm),label = "p.signif", method = "t.test")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1 ))-> p.mir.arm.hyp
  #ggsave("violin.arm.hyp.pdf",units = "cm", width = 100,height = 600,limitsize = F)
  
  ##3、缺氧后是否出现3p、5p转换
  vsted <- right_join(data.frame(id=rownames(vsted.nor),vsted.nor),data.frame(id=rownames(vsted.hyp),vsted.hyp))
  vsted <- data.frame(vsted,row.names = 1)
  for (i in 1:ncol(vsted)) {
    vsted[,i][is.na(vsted[,i])]=0
  }
  
  mir.df.vst.all <- data.frame(vsted)
  for (i in 1:nrow(mir.df.vst.all)) {
    mir.df.vst.all$arm[i] <- ifelse(
      strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.all)[3],"-")[[1]])]=="3p","3p",
      ifelse(strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][length(strsplit(rownames(mir.df.vst.all)[3],"-")[[1]])]=="5p","5p","no"
      )
    )
    mir.df.vst.all$mirname[i] <- paste(
      strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.all)[3],"-")[[1]])-1)][1],
      strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.all)[3],"-")[[1]])-1)][2],
      strsplit(rownames(mir.df.vst.all)[i],"-")[[1]][c(1:length(strsplit(rownames(mir.df.vst.all)[3],"-")[[1]])-1)][3],
      sep = "-")
  }
  
  mir.df.vst.all <- na.omit(mir.df.vst.all)
  
  
  mir.df.vst.all.3p <- mir.df.vst.all[mir.df.vst.all$arm=="3p",]
  mir.df.vst.all.5p <- mir.df.vst.all[mir.df.vst.all$arm=="5p",]
  
  
  mir.df.vst.all[,c(1:48)] <- log2(mir.df.vst.all[,c(1:48)])
  
  library(data.table)
  mir.df.vst.all.long <- melt(mir.df.vst.all,id=c("mirname","arm"))
  colnames(mir.df.vst.all.long) <- c("miRNA","Arm","Sample","Reads") 
  
  ##
  for (i in 1:nrow(mir.df.vst.all.long)) {
    mir.df.vst.all.long$Tissue[i] <- strsplit(as.character(mir.df.vst.all.long$Sample)[i],"_")[[1]][1]
    mir.df.vst.all.long$Condition[i] <- strsplit(as.character(mir.df.vst.all.long$Sample)[i],"_")[[1]][2]
  }
  
  # BiocManager::install("introdataviz")
  # BiocManager::install("CBNplot")
  library(tidyverse)
  # library(introdataviz)
  library(ggunchained)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(ggprism)
  
  mir.df.vst.all.long.rmarmless <- mir.df.vst.all.long[mir.df.vst.all.long$Arm!="no",]
  ggplot(data = mir.df.vst.all.long.rmarmless,aes(x=Condition,y=Reads,fill=Arm))+
    geom_split_violin()+
    facet_wrap(Tissue~miRNA,ncol = 4)+
    stat_summary(fun.data = "mean_sd", position = position_dodge(0.15), geom = "errorbar", width=.1)+
    stat_summary(fun = "mean", geom = "point", position = position_dodge(0.15), show.legend = F)+
    stat_compare_means(aes(group = mir.df.vst.all.long.rmarmless$Arm),label = "p.signif", method = "t.test")+
    #  theme_prism()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1 ))-> p.mir.arm.all
  
  ##拼图和导出
  # library(cowplot)
  # p.pca.mir.nor.prism <- p.pca.mir.nor+theme_prism()+guides(x = "prism_offset_minor", y = "prism_offset_minor")#+theme(legend.position = "none")
  # p.pca.mir.hyp.prism <- p.pca.mir.hyp+theme_prism()+guides(x = "prism_offset_minor", y = "prism_offset_minor")#+theme(legend.position = "none")
  # ggsave(plot = p.pca.mir.nor.prism,"p.pca.mir.nor.prism.pdf",units = "cm", width = 16,height = 16,limitsize = F)
  # ggsave(plot = p.mir.nor.spc+theme(axis.text.x = element_text(angle=45, hjust = 1)),"p.mir.nor.spc.pdf",units = "cm", width = 12,height = 125,limitsize = F)
  # ggsave(plot = p.mir.hyp.spc+theme(axis.text.x = element_text(angle=45, hjust = 1)),"p.mir.hyp.spc.pdf",units = "cm", width = 12,height = 125,limitsize = F)
  # ggsave(plot = p.pca.mir.hyp,"p.pca.mir.hyp.prism.pdf",units = "cm", width = 16,height = 16,limitsize = F)
  # plot_grid(p.pca.mir.nor.prism,p.pca.mir.hyp.prism)
  # ggsave(plot = p.mir.arm.nor,"violin.arm.nor.pdf",units = "cm", width = 75,height = 600,limitsize = F)
  # ggsave(plot = p.mir.arm.hyp,"violin.arm.hyp.pdf",units = "cm", width = 75,height = 600,limitsize = F)
}##miRNA.mmscale

if (T) {
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(readxl)
  library(ggprism)
  setwd("C:/Users/colet/Desktop/re snc atlas/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  setwd("C:/Users/colet/Desktop/re snc atlas/re length/")
  list.mir.len=list.files(pattern = "length")
  mir.len.in = map(list.mir.len, ~ read.table(.,header = F))
  
  for (i in 1:24) {
    mir.len.in[[i]] <- data.frame(length=mir.len.in[[i]],condition="Control",sample=sample$tissue[i])
  }
  for (i in 25:48) {
    mir.len.in[[i]] <- data.frame(length=mir.len.in[[i]],condition="Hypoxia",sample=sample$tissue[i])
  }
  mir.len.all=NULL
  for (i in 1:48) {
    mir.len.all <- rbind(mir.len.all,mir.len.in[[i]])
  }
  colnames(mir.len.all) <- c("Length","Condition","Sample")
  library(ggplot2)
  library(ggridges)
  head(mir.len.all)
  ggplot()+
    geom_density(data = mir.len.all,aes(x=Length,fill=Sample),lwd = 1, linetype = 1) + 
    scale_fill_manual(values = c(viridis::viridis(n=8,alpha = .5)))+
    facet_wrap(~Condition,nrow = 2)+
    xlim(0,100)+
    #  geom_boxplot(data = mir.len.all, aes(x=Length,y=0.15))+
    guides(x = "prism_offset_minor", y = "prism_offset_minor")+ 
    ggprism::theme_prism()->p
  # ggsave(p,filename = "t.len.trf.pdf",units = "cm",width = 16,height = 16)
  
}##tRF分布长度

if (T) {
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(ggprism)
  setwd("C:/Users/colet/Desktop/re snc atlas/quantified/")
  sample <- read_xlsx("sampleID.xlsx")
  ##读取tRFs
  setwd("C:/Users/colet/Desktop/re snc atlas/re length/trf/")
  list.trf.len=list.files(pattern = "length")
  trf.len.in = map(list.trf.len, ~ read.table(.,header = F))
  for (i in 1:24) {
    trf.len.in[[i]] <- data.frame(table(trf.len.in[[i]]),condition="Control",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.trf.len,"_")[[i]][2],strsplit(list.trf.len,"_")[[i]][3],sep="_"))])
  }
  for (i in 25:48) {
    trf.len.in[[i]] <- data.frame(table(trf.len.in[[i]]),condition="Hypoxia",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.trf.len,"_")[[i]][2],strsplit(list.trf.len,"_")[[i]][3],sep="_"))])
  }
  trf.len.all=NULL
  gc()
  for (i in 1:48) {
    trf.len.all <- rbind(trf.len.all,trf.len.in[[i]])
  }
  gc()
  trf.len.all <- data.frame(trf.len.all,rna="tRF")
  colnames(trf.len.all) <- c("Length","Counts","Condition","Sample","RNA")
  gc()
  len.all.1 <- trf.len.all
  
  ##读取miRNA长度
  setwd("C:/Users/colet/Desktop/re snc atlas/re length/mir/")
  list.mir.len=list.files(pattern = "length")
  mir.len.in = map(list.mir.len, ~ read.table(.,header = F))
  
  avg.len.mir=NULL
  for (i in 1:48) {
    avg.len.mir <- rbind(mir.len.in[[i]],avg.len.mir)  
  }
  avg.len.mir <- mean(na.omit(as.matrix(avg.len.mir)))
  avg.len.mir
  
  
  for (i in 1:24) {
    mir.len.in[[i]] <- data.frame(table(mir.len.in[[i]]),condition="Control",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.mir.len,"_")[[i]][3],strsplit(list.mir.len,"_")[[i]][4],sep="_"))])
  }
  for (i in 25:48) {
    mir.len.in[[i]] <- data.frame(table(mir.len.in[[i]]),condition="Hypoxia",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.mir.len,"_")[[i]][3],strsplit(list.mir.len,"_")[[i]][4],sep="_"))])
  }
  mir.len.all=NULL
  gc()
  for (i in 1:48) {
    mir.len.all <- rbind(mir.len.all,mir.len.in[[i]])
  }
  gc()
  mir.len.all <- data.frame(mir.len.all,rna="miRNA")
  colnames(mir.len.all) <- c("Length","Counts","Condition","Sample","RNA")
  gc()
  len.all.1 <- rbind(len.all.1,mir.len.all)
  
  
  
  ##读取piRNA
  setwd("C:/Users/colet/Desktop/re snc atlas/re length/pir/")
  list.pir.len=list.files(pattern = "length")
  pir.len.in = map(list.pir.len, ~ read.table(.,header = F))
  
  avg.len.pi=NULL
  for (i in 1:48) {
    avg.len.pi <- rbind(pir.len.in[[i]],avg.len.pi)  
  }
  avg.len.pi <- mean(na.omit(as.matrix(avg.len.pi)))
  avg.len.pi
  
  for (i in 1:24) {
    pir.len.in[[i]] <- data.frame(table(pir.len.in[[i]]),condition="Control",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.pir.len,"_")[[i]][3],strsplit(list.pir.len,"_")[[i]][4],sep="_"))])
  }
  for (i in 25:48) {
    pir.len.in[[i]] <- data.frame(table(pir.len.in[[i]]),condition="Hypoxia",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.pir.len,"_")[[i]][3],strsplit(list.pir.len,"_")[[i]][4],sep="_"))])
  }
  pir.len.all=NULL
  gc()
  for (i in 1:48) {
    pir.len.all <- rbind(pir.len.all,pir.len.in[[i]])
  }
  gc()
  pir.len.all <- data.frame(pir.len.all,rna="piRNA")
  colnames(pir.len.all) <- c("Length","Counts","Condition","Sample","RNA")
  gc()
  len.all.1 <- rbind(len.all.1,pir.len.all)
  
  ##读取snRNA长度
  setwd("C:/Users/colet/Desktop/re snc atlas/re length/sn/")
  list.snr.len=list.files(pattern = "length")
  snr.len.in = map(list.snr.len, ~ read.table(.,header = F))
  
  avg.len.snr=NULL
  for (i in 1:48) {
    avg.len.snr <- rbind(snr.len.in[[i]],avg.len.snr)  
  }
  avg.len.snr <- mean(na.omit(as.matrix(avg.len.snr)))
  avg.len.snr
  
  
  for (i in 1:24) {
    snr.len.in[[i]] <- data.frame(table(snr.len.in[[i]]),condition="Control",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.snr.len,"_")[[i]][3],strsplit(list.snr.len,"_")[[i]][4],sep="_"))])
  }
  for (i in 25:48) {
    snr.len.in[[i]] <- data.frame(table(snr.len.in[[i]]),condition="Hypoxia",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.snr.len,"_")[[i]][3],strsplit(list.snr.len,"_")[[i]][4],sep="_"))])
  }
  snr.len.all=NULL
  gc()
  for (i in 1:48) {
    snr.len.all <- rbind(snr.len.all,snr.len.in[[i]])
  }
  gc()
  snr.len.all <- data.frame(snr.len.all,rna="snRNA")
  colnames(snr.len.all) <- c("Length","Counts","Condition","Sample","RNA")
  gc()
  len.all.1 <- rbind(len.all.1,snr.len.all)
  
  # ##rsRNA
  # ##读取rsRNA长度
  # setwd("C:/Users/colet/Desktop/re snc atlas/re length/rsr/")
  # list.rsr.len=list.files(pattern = "length")
  # rsr.len.in = map(list.rsr.len, ~ read.table(.,header = F))
  # for (i in 1:24) {
  #   rsr.len.in[[i]] <- data.frame(table(rsr.len.in[[i]]),condition="Control",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.rsr.len,"_")[[i]][3],strsplit(list.rsr.len,"_")[[i]][4],sep="_"))])
  # }
  # for (i in 25:48) {
  #   rsr.len.in[[i]] <- data.frame(table(rsr.len.in[[i]]),condition="Hypoxia",sample=sample$tissue[sample$sample_name==gsub("-",".",paste(strsplit(list.rsr.len,"_")[[i]][3],strsplit(list.rsr.len,"_")[[i]][4],sep="_"))])
  # }
  # rsr.len.all=NULL
  # gc()
  # for (i in 1:48) {
  #   rsr.len.all <- rbind(rsr.len.all,rsr.len.in[[i]])
  # }
  # gc()
  # rsr.len.all <- data.frame(rsr.len.all,rna="rsRNA")
  # colnames(rsr.len.all) <- c("Length","Counts","Condition","Sample","RNA")
  # gc()
  # len.all.1 <- rbind(len.all.1,rsr.len.all)
  
  
  ##读取其他RNA
  setwd("C:/Users/colet/Desktop/re snc atlas/length distribution/")
  ##read files as a large list
  listrs=list.files(pattern = "csv")
  
  ##rename every samples with readable names
  length.in = map(listrs, ~ read.csv(.,header = T))
  
  tissue.name <- c("length","thymus","kidney","brain","liver","lung","intestines","heart","spleen") 
  rna <- c("lncRNA","miRNA","piRNA","scaRNA","snoRNA","snRNA","tRFs")
  
  length.combine=NULL
  
  for (i in c(2:7)) {
    length.nor <- length.in[[i]][,c(1:9)]
    length.hyp <- length.in[[i]][,c(1,10:17)]
    colnames(length.nor) <- tissue.name
    colnames(length.hyp) <- tissue.name
    
    length.nor <- data.frame(length.nor,condition="Control",rna=rna[i])
    single.length.nor <- data.frame(length=unique(length.nor$length))
    length.nor <- left_join(single.length.nor,length.nor)
    
    length.hyp <- data.frame(length.hyp,condition="Hypoxia",rna=rna[i])
    single.length.hyp <- data.frame(length=unique(length.hyp$length))
    length.hyp <- left_join(single.length.hyp,length.hyp)
    
    length.combine <- rbind(length.combine,length.nor,length.hyp)
  }
  
  library(data.table)
  library(ggplot2)
  library(ggalt)
  library(ggbreak)
  library(ggprism)
  library(ggbump)
  library(viridis)
  length.combine.long <- melt(length.combine,id=c("length","condition","rna"))
  other <- length.combine.long%>%filter(rna==c("snoRNA","scaRNA"))
  other.t <- other[,c(1,5,2,4,3)]
  colnames(other.t) <- c("Length","Counts","Condition","Sample","RNA")
  all.len <- rbind(len.all.1,other.t)
  all.len <- na.omit(all.len)
  all.len[,c(1,2)] <- data.frame(as.numeric(all.len$Length),as.numeric(all.len$Counts))
  write.csv(all.len,"all.len",row.names = F)
  all.len.2 <- read.csv("all.len")
  library(ggplot2)
  library(ggprism)
  library(viridis)
  ggplot()+
    #  geom_line(data = all.len.2%>%filter(Length>0&Counts>100),aes(x=Length,y=log10(Counts),color=RNA),lwd=1)+
    geom_bump(data = all.len.2%>%filter(Length<200&Counts>1),aes(x=Length,y=log10(Counts),color=RNA),lwd=1,smooth = 5)+
    # geom_bump(data = all.len.2%>%filter(RNA==c("miRNA","tRF","piRNA")&Length<100&Counts>10),aes(x=Length,y=log10(Counts),color=RNA),lwd=1,smooth = 5)+
    # geom_bump(data = all.len.2%>%filter(RNA!=c("miRNA","tRF","piRNA")&Length<100&Counts>10),aes(x=Length,y=log10(Counts),color=RNA),lwd=1,smooth = 5)+
    # geom_boxplot(data = all.len.2%>%filter(Length<200&Counts>100),aes(x=Length,y=1,color=RNA))+
    # geom_boxplot(data = all.len.2%>%filter(RNA!=c("miRNA","piRNA"),Length<100&Counts>100),aes(x=Length,y=0,color=RNA))+
    coord_cartesian(xlim =c(0, 150),ylim = c(0,6))+
    labs(x="Fragment Size (Bases)",y="log10(Reads)")+
    scale_color_manual(values = c("#85BA8F", "#A3C8DC",
                                  "#349839","#EA5D2D",
                                  "#EABB77","#F09594","#2072A8"))+
    # scale_color_manual(values = c(viridis(6)))+
    facet_wrap(~Condition+Sample,nrow = 2)+
    guides(x="prism_offset_minor",y = "prism_offset_minor") +
    theme_prism()->p.length.distribution
  p.length.distribution
  # ggsave(dpi=800,filename = "p.n.length.distribution.png", p.length.distribution, units = "cm", width = 38,height = 18)
}##新长度分布（base）

if (T) {
  ##注释占比横向堆积图
  #gtf文件的的RNA个数
  sno.gtf=1706
  snr.gtf=1515
  sca.gtf=74
  mir.gtf=1336
  trf.gtf=203395
  pi.gtf=1269304
  
  sno.a=nrow(sno.in)
  snr.a=nrow(snr.in)
  sca.a=nrow(sca.in)
  mir.a=1238
  trf.a=nrow(trf.in.b)
  pi.a=4152
  
  rate <- data.frame(id=c("snoRNA","snRNA","scaRNA","miRNA","tRF","piRNA"),
                     Dected=c(sno.a,snr.a,sca.a,mir.a,trf.a,pi.a),
                     Annotaded=c(sno.gtf,snr.gtf,sca.gtf,mir.gtf,trf.gtf,pi.gtf))
  
  rate.long <- melt(rate)
  
  library(ggbreak)
  ggplot()+
    geom_col(data = rate.long%>%filter(variable!="Dected"),aes(y=id, x = value,fill=variable))+
    geom_col(data = rate.long%>%filter(variable=="Dected"),aes(y=id, x = value,fill=variable))+
    scale_x_break(breaks = c(4500,70000),scales = "free")+
    labs(x="Number of RNA genes",y="RNA Type",fill="")+
    # scale_fill_discrete(breaks=c("gtf","dected"),
    #                     labels=c( "Annotaded","Dected in the\nprecent study"))+
    scale_fill_manual(values=c(alpha("#2fa1dd",.5),alpha("#f87669",.95)))+
    # scale_fill_npg()+
    theme_classic()+
    theme(text=element_text(family = "sans",colour ="gray30",size = 12),
          axis.text.x = element_text(angle = 45,hjust = 1),
          legend.text=element_text(colour ="gray30",size = 8),
          legend.title=element_text(colour ="gray30",size = 10),
          legend.key.size=unit(4,units = "mm"),
          legend.position=c(0.10,0.88),
          axis.line = element_line(size = 0.4,colour = "gray30"),
          axis.ticks = element_line(size = 0.4,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"),
          plot.margin=unit(x=c(.2,.2,.2,.2),
                           units="inches"))->p.dected.anno
  setwd("C:/Users/colet/Desktop/")
  # ggsave(p.dected.anno,filename = "p.dected.anno.pdf",units = "cm",width = 16,height = 8)
}##注释占比横向堆积图

library(readxl)
##miRNA靶点预测
setwd("C:/Users/colet/Desktop/quantified/quantified/20220605/quantified/")
sample <- read_xlsx("sampleID.xlsx")
##read files as a large list
mir.in <- read.table("miR.quantified",header = T)

#去除mir
mir.in <- mir.in[-grep("mir",mir.in$Geneid),]

mir.in.c=data.frame(mir.in,row.names = 1)
mir.in.b <- data.frame(mir.in[,7:ncol(mir.in)],row.names = mir.in$Geneid)

i=1
for (i in 1:ncol(mir.in.b)) {
  colnames(mir.in.b)[i] <- sample$readable_name[sample$sample_name==paste(strsplit(colnames(mir.in.b),"_")[[i]][1],strsplit(colnames(mir.in.b),"_")[[i]][2],sep = "_")]
}

mir.df=mir.in.b
all.mir.df=mir.in.b
mir.tpm.1 <- cbind(mir.in$Length,mir.in.b)
a=ncol(mir.tpm.1)
m.tpm <- (mir.tpm.1[,2:a]/mir.tpm.1[,1])/colSums(mir.tpm.1[,2:a]/mir.tpm.1[,1])*1e6
s.m.tpm <- na.omit(t(scale(t(m.tpm))))

# nor.m <- data.frame(s.m.tpm[rownames(s.m.tpm)%in%sig,1:24],condition="Normoxia")
# hyp.m <- data.frame(s.m.tpm[rownames(s.m.tpm)%in%sig,25:48],condition="Hypoxia")

nor.m <- data.frame(s.m.tpm[rownames(s.m.tpm),1:24],condition="Normoxia")
hyp.m <- data.frame(s.m.tpm[rownames(s.m.tpm),25:48],condition="Hypoxia")

colnames(nor.m) <- c(paste0(c("thymus","kidney","brain","liver", "lung","intestines","heart","spleen"),"-1"),
                     paste0(c("thymus","kidney","brain","liver", "lung","intestines","heart","spleen"),"-2"),
                     paste0(c("thymus","kidney","brain","liver", "lung","intestines","heart","spleen"),"-3"),"Condition"
                     )

colnames(hyp.m) <- c(paste0(c("thymus","kidney","brain","liver", "lung","intestines","heart","spleen"),"-1"),
                     paste0(c("thymus","kidney","brain","liver", "lung","intestines","heart","spleen"),"-2"),
                     paste0(c("thymus","kidney","brain","liver", "lung","intestines","heart","spleen"),"-3"),"Condition"
                     )



u.nor.hyp.m <- rbind(nor.m,hyp.m)
# u.nor.hyp.m <- u.nor.hyp.m[rownames(u.nor.hyp.m)%in%sig,]
library(ComplexHeatmap)
library(circlize)
##准备画热图
##1、准备行分裂（缺氧与否）
rowsplit <- list(Condition=unlist(strsplit(colnames(m.tpm),"_"))[seq(2,144,3)])
##2、准备列分裂（组织别）
colsplit <- unlist(strsplit(colnames(u.nor.hyp.m),"-"))[seq(1,48,2)]
# names(colsplit) <- colnames(u.nor.hyp.m)
# colsplit <- data.frame(unlist(strsplit(colnames(m.tpm),"_"))[seq(1,144,3)])
# colnames(colsplit)=NULL
##画热图
col = colorRamp2(breaks = c(-2, 0, 4), colors = c("#2fa1dd","white","#f87669"))
cl_col = cola:::brewer_pal_set2_col[1:8]
names(cl_col) <- c(colsplit[1:8])
rw_col <- cola:::brewer_pal_set1_col[1:2]
names(rw_col) <- c("Normoxia","Hypoxia")

rsp <- u.nor.hyp.m$Condition
names(rsp) <- rownames(u.nor.hyp.m)

ht = Heatmap(u.nor.hyp.m[,-ncol(u.nor.hyp.m)], name = "TPM",
             row_split = u.nor.hyp.m$Condition, 
             column_split = colsplit,
             # column_title = "D) Columns are pre-split by 'class'",
             show_row_names = FALSE, show_row_dend = FALSE, show_column_names = FALSE, 
             # bottom_annotation = HeatmapAnnotation(df = anno, col = anno_col),
             top_annotation = HeatmapAnnotation(Tissue = colsplit, name = "Tissue",
                                                col = list(Tissue = cl_col)
                                                ),
             col=colorRamp2(breaks = c(-1,0,3), colors = c("#2fa1dd","white","#f87669")),
             # col=colorRamp2(breaks = c(-2,-1,0,1,2,3,4), colors = c(viridis::viridis(7)))
             border = T
             )
print(ht)







# 
organlist <- list("thymus","kidney","brain","liver",
              "lung","intestines","heart","spleen")
i=1
for (i in 1:8) {
  tem.organ <- all.mir.df[,c(i+8,i,i+16,i+24,i+40,i+32)]
  coldat <- data.frame(id = colnames(tem.organ),dex = rep(c("Normoxia","Hypoxia"),each = 3))
  dds <- DESeqDataSetFromMatrix(countData = tem.organ, colData = coldat, design = ~dex)
  dds <- DESeq(dds)
  res <- data.frame(na.omit(results(dds)))
  organlist[[i]] <- res
}

sig = NULL
for (i in 1:8) {
  sig.1 <- rownames(data.frame(organlist[[i]])%>%filter(abs(log2FoldChange)>1))
  sig <- unlist(list(sig), list(sig.1))
}
sig <- unique(sig)



library(multiMiR)

tg.mir <- get_multimir(org = "rno", mirna = rownames(organlist[[1]]))
tg.mir <- data.frame(tg.mir@data)
sym.tg.mir <- unique(tg.mir$target_symbol)
library(clusterProfiler)
library(org.Rn.eg.db)
enrichGO(gene = sym.tg.mir,OrgDb = "org.Rn.eg.db",
         keyType = "SYMBOL",
         ont = "ALL",
         pvalueCutoff = 1,
         qvalueCutoff = 1)








##mir差异分析（对照对对照、处理对处理；组织别）

#1、对照组；组织别
df.col = data.frame(id = colnames(mir.df.nor.all),dex = c(rep("thymus",3),
                                                          rep("kidney",3),
                                                          rep("brain",3),
                                                          rep("liver",3),
                                                          rep("lung",3),
                                                          rep("intestines",3),
                                                          rep("heart",3),
                                                          rep("spleen",3)
))
dds<-DESeqDataSetFromMatrix(countData = mir.df.nor.all, colData = df.col,design = ~dex)
## Prefiltering
filt <- rowSums(counts(dds) < 10) > dim(df.col)[1]*0.9
dds <- dds[!filt,]

## Perform DESeq2()
dds = DESeq(dds)
res = results(dds, pAdjustMethod = "bonferroni")

########################再画图谱
mir.in.c=mir.in
for (i in 1:nrow(mir.in.c)) {
  mir.in.c$Chr[i] <- strsplit(mir.in.c$Chr,";")[[i]][1]
  mir.in.c$Start[i] <- strsplit(mir.in.c$Start,";")[[i]][1]
  mir.in.c$End[i] <- strsplit(mir.in.c$End,";")[[i]][1]
  mir.in.c$Strand[i] <- strsplit(mir.in.c$Strand,";")[[i]][1]
}
sub(mir.in.c$Chr,"NA","21")
mir.in.c$Chr[order(mir.in.c$Chr)]
# mir.in.c.nor$Chr <- factor(mir.in.c.nor$Chr,levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
#                                                        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chrX"))
mir.in.c.nor <- mir.in.c[,c(1:30)]
mir.in.c.nor$mean <- abs(rowMeans(mir.in.c.nor[,-c(1:6)]))
colors <-c("#FED439FF","#709AE1FF",
           "#D5E4A2FF","#197EC0FF","#F05C3BFF","#46732EFF",
           "#71D0F5FF","#370335FF","#075149FF","#C80813FF","#91331FFF",
           "#1A9993FF","#FD8CC1FF")
ggplot()+
  geom_bar(data = mir.in.c.nor,aes(x=Start,y=abs(log10(mean)),fill=Chr,color=Chr),stat="identity",alpha=0.8)+
  coord_polar(start = pi/3)+
  labs(x="",y="")+
  ylim(-7,4)+
  scale_fill_manual(values = c(viridis::viridis(21)))+
  scale_color_manual(values = c(viridis::viridis(21)))+
  theme_void()+
  theme(axis.text = element_text(size = 0))

bed.mir <- data.frame(mir.in.c.nor$Chr,mir.in.c.nor$Start,mir.in.c.nor$End,
                      value1=abs(log10(rowMeans(mir.in.c.nor[,-c(1:6)]))))
bed.mir$value1[bed.mir$value1=="Inf"]=0
colnames(bed.mir) <- c("chr","start","end","value1")
for (i in 1:nrow(bed.mir)) {
  # <- paste0("chr",bed.mir$chr[i])
  bed.mir$value1[i] <- (bed.mir$value1[i]-min(bed.mir$value1))/(max(bed.mir$value1)+min(bed.mir$value1))
}

library(circlize) 

bed = generateRandomBed(nr = 200)
head(bed)
t <- head(bed.mir,200)
circos.initializeWithIdeogram(species = "rn6")
circos.genomicTrack(t,
                    panel.fun = function(region, value, ...) {
                      
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = 1)})

circos.genomicTrack(bed.mir, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })




circos.clear()
circos.par("start.degree" = 90)
circos.par("gap.degree" = rep(c(2, 2), 12), ADD = TRUE)
circos.initializeWithIdeogram(species = "hg38", plotType = c("axis"))
text(0, 0, "Sample1", cex = 1)

####### 染色体使用不同颜色的方框表示：
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(24))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.6, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)

####### 添加SNV、INDEL的信息：
circos.genomicTrack(bed, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = 1)},
                    stack = F, track.height = 0.1)

####### 添加CNV：
circos.genomicTrack(bed, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "green"), ...)
                       circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    }, stack = F, track.height = 0.1)
# 大于阈值（例如0）即为红色柱子，否则为绿色柱子

####### 热图：
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
bed.1 = generateRandomBed(nr = 70, nc = 4)
# bed.1 = bed.1[1:10, ]
circos.genomicHeatmap(bed.1, col = col_fun, side = "inside", border = "white",
                      connection_height = mm_h(2),heatmap_height = 0.1,
                      line_col = as.numeric(factor(bed[[1]])))



#######################
  library(ggridges)
head(trf.len.all)
rm(trf.len.in)
gc()
ggplot()+
  geom_density(data = len.all.3%>%filter(Length<100),aes(x=Length,fill=RNA),lwd = 1, linetype = 1) + 
  scale_fill_manual(values = c(viridis::viridis(n=8,alpha = .5)))+
  facet_wrap(Condition~Sample,nrow = 2)+
  xlim(0,100)+
  #  geom_boxplot(data = mir.len.all, aes(x=Length,y=0.15))+
  guides(x = "prism_offset_minor", y = "prism_offset_minor")+ 
  ggprism::theme_prism()->p
p
gc()
# ggsave(p,filename = "t.len.trf.png",dpi=300,units = "cm",width = 16,height = 32)


length.combine.long%>%filter(rna=="snoRNA")

gc()
library(ggplot2)
library(ggridges)
head(trf.len.all)
rm(trf.len.in)
gc()
ggplot()+
  geom_density(data = trf.len.all%>%filter(Condition=="Control"&Length<100),aes(x=Length,fill=RNA),lwd = 1, linetype = 1) + 
  scale_fill_manual(values = c(viridis::viridis(n=8,alpha = .5)))+
  facet_wrap(~Sample,nrow = 2)+
  xlim(0,100)+
  #  geom_boxplot(data = mir.len.all, aes(x=Length,y=0.15))+
  guides(x = "prism_offset_minor", y = "prism_offset_minor")+ 
  ggprism::theme_prism()->p
p
gc()
# ggsave(p,filename = "t.len.trf.png",dpi=300,units = "cm",width = 16,height = 16)


#  scale_color_manual(values = cols)



data.frame(length.combine)

unique(length.combine$length)


list=list.files(pattern = "miR*")

##rename every samples with readable names
a = map(list, ~ read.table(.,header = T))





















## 建立一个list对象用于下面的按组织提取 20220530：首先完成denovo miRNAs部分(a[[1]],a[[2]])
kloop=7
e<-c(data.frame("thymus"),data.frame("kidney"),data.frame("brain"),
     data.frame("liver"),data.frame("lung"),
     data.frame("intestines"),data.frame("heart"),data.frame("spleen"))
deseq_done<-c(data.frame("thymus"),data.frame("kidney"),data.frame("brain"),
              data.frame("liver"),data.frame("lung"),
              data.frame("intestines"),data.frame("heart"),data.frame("spleen"))
for (i in 2:length(a)+1) {
  e[[i-1]] <- a[[kloop]][,c(1,i,i+8*1,i+8*2,i+8*3,i+8*4,i+8*5)] #1,thymus 2,kidney 3,brain 4,liver 5,lung 6,intestines 7,heart 8,spleen
  #  rownames(e[[i-1]]) <- e[[i-1]][,1]
  e[[i-1]] <- e[[i-1]][,-1]
}
i=2
e[[i-1]] <- a[[kloop]][,c(1,i,i+8*1,i+8*2,i+8*3,i+8*4,i+8*5)] #1,thymus 2,kidney 3,brain 4,liver 5,lung 6,intestines 7,heart 8,spleen
#rownames(e[[i-1]]) <- e[[i-1]][,1]
e[[i-1]] <- e[[i-1]][,-1]

##差异分析(20220530 denovo miRNA ：DONE lnc:done ,pi)
library(DESeq2)
i = 1
for (i in 1:length(e)) {
  COUNTDATA <- e[[i]]
  as.matrix(COUNTDATA)
  COLDATA <- data.frame(id=colnames(COUNTDATA),
                        dex=c(rep("control",3),rep("treat",3))
  )
  colnames(COUNTDATA) == COLDATA$id#判断是否一致，-/.
  colnames(COUNTDATA) <- COLDATA$id
  dds <- DESeqDataSetFromMatrix(countData = COUNTDATA, 
                                colData = COLDATA,
                                design = ~dex
  ) ## differential analysis
  dds <- DESeq(dds) 
  #  rbind(Nord,Nord1)
  res1 <- results(dds)
  head(res1)
  res2 <- data.frame(res1)
  deseq_done[[i]] <- na.omit(res2) #差异分析结果
  deseq_done[[i]]  <- data.frame(deseq_done[[i]],sample_tissue=strsplit(names(deseq_done[i]),split = "X.")[[1]][2])
  #  deseq_done[[i]] <- data.frame(deseq_done[[i]],sample_tissue = names(deseq_done[i]))
  ##重定向归一化数据至a
  Nord <- round(counts(dds, normalized=T),3) ## extract normalized counts , add them up
  e[[i]] <- Nord
}



##火山图##folded
if (FALSE) {
  bind2 <- rbind(deseq_done[[1]],deseq_done[[2]],deseq_done[[3]],
                 deseq_done[[4]],deseq_done[[5]],deseq_done[[6]],
                 deseq_done[[7]],deseq_done[[8]]
  )
  
  DF<-cbind(gene_name=row.names(bind2),bind2)
  row.names(DF)=NULL
  LOG2FC<-1
  PVALUE<-0.05
  XLINE<-c(-1,1)
  YLINE<-c(0.05)
  DF$group<-ifelse(DF$log2FoldChange>=LOG2FC& DF$pvalue<=PVALUE,"Up",
                   ifelse(DF$log2FoldChange<=-LOG2FC&DF$pvalue<=PVALUE,"Down","Not sig"))
  DF<-na.omit(DF)#去掉缺失值
  DF$dramatically_sig <- (-log10(DF$pvalue)) #挑选极为显著的基因
  DRAMATICALLY_sig<-DF[DF$dramatically_sig>=10,]
  DF_sig<-DF[-grep("Not sig",DF$group),]
  ##火山图
  #install.packages("remote")
  #remotes::install_github("tylermorganwall/rayshader")
  library(ggplot2)
  ggplot(data = DF,aes(x=log2FoldChange,y=-log10(pvalue)))+
    geom_point()+
    facet_wrap(.~sample_tissue,nrow=2)->p
  #火山图美化
  library(ggrepel)
  p+
    geom_point(aes(color=group))+
    scale_color_manual(values = c("#2fa1dd", "grey", "#f87669"))+
    #加极其显著标签
    geom_label_repel(data=DRAMATICALLY_sig,aes(x=log2FoldChange, y=-log10(pvalue),label=gene_name))+
    #坐标轴
    scale_x_continuous(limits = c(-6,6), breaks = c(-6,-4,-2,0,2,4,6))+
    #坐标轴标题
    # labs(x='logFC',y='-log10(pvalue)',title = '热图')+#坐标轴
    scale_y_continuous(limits = c(-0,5))+
    #阈值线
    geom_vline(xintercept = XLINE, linetype = 5,size=0.8)+
    geom_hline(yintercept=-log10(YLINE), linetype=5,size=0.8) +
    
    #主题
    theme_bw()+
    theme(#panel.background = element_blank(),#空白背景
      panel.border = element_rect(linetype = 1, color="black",size=1),
      #       axis.line = element_line(size = 1,colour = 'black'),#坐标轴线条
      axis.ticks = element_line(size = 1,colour = 'black'),#坐标轴格子
      axis.title = element_text(size = 12, color = 'black',face = "bold"),#坐标轴标题字体
      axis.text = element_text(size = 10,color = 'black',face = "bold"),#坐标轴字体
      legend.text = element_text(size = 10,color = 'black',face = "bold"),
      legend.title = element_text(size = 12,color = 'black',face = "bold")
      #       plot.title = element_text(size = 50 , hjust=0.5),#图表标题
      #       panel.border = element_rect(linetype = 1, color="black",size=1),
      #       legend.background = element_rect(fill = NULL, colour = NA)
    ) ->p1
  # ggsave(file= "miR_vocanol.png", plot = p1,dpi = 800,width=11, height=6)
}

#unlist
tissue <- c("thymus","kidney","brain",
            "liver","lung","intestines","heart","spleen")
library(dplyr)
library(data.table)
e[[1]] <- data.frame(e[[1]])
t <- data.frame(row.names = NULL,geneid=factor(rownames(e[[1]])),
                round(data.frame(NOR=rowMeans(e[[1]][,c(1:3)]),
                                 HYP=rowMeans(e[[1]][,c(4:6)])),3),tissue="thymus") %>% melt(geneid = geneid,tissue = tisste)
colnames(t) <- c("geneid","tissue","condition","counts")
for (i in 2:length(e)) {
  e[[i]] <- data.frame(e[[i]])
  u <- data.frame(row.names = NULL,geneid=factor(rownames(e[[i]])),
                  round(data.frame(NOR=rowMeans(e[[i]][,c(1:3)]),
                                   HYP=rowMeans(e[[i]][,c(4:6)])),3),tissue=tissue[[i]]) %>% melt(geneid = geneid,tissue = tisste)
  colnames(u) <- c("geneid","tissue","condition","counts")
  t <- rbind(t,u)
}
##配对箱线图
package.list=c("tidyverse","ggsignif","ggsci","ggprism")
for (package in package.list) {
  if (!require(package,character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

if (T) {
  t <- t[t$counts>1,]
  t[t$counts<10000,] %>%
    ggplot(aes(condition,counts)) +
    stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
    geom_boxplot(position=position_dodge(width =0.2),width=0.4)+
    geom_line(aes(group=geneid),position = position_dodge(0.2),color="grey80") +
    geom_point(aes(fill=condition,group=geneid,size=counts,alpha=counts),pch=21,
               position = position_dodge(0.2))+
    scale_size_continuous(range=c(1,3))+
    #  geom_signif(comparisons = list(c("NOR","HYP")),
    #              map_signif_level=T,vjust=0.5,color="black",
    #              textsize=5,test=wilcox.test,step_increase=0.1)+
    facet_wrap(.~tissue,nrow=1)+
    scale_fill_npg()+
    scale_x_discrete(guide = "prism_bracket")+
    scale_y_continuous(limits = c(0,25)#,minor_breaks = seq(0,50,5)
                       ,guide = "prism_offset_minor")+
    labs(x="conditon",y="counts")+
    theme_prism(base_line_size =0.5)+
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
          axis.line = element_line(color = "black",size = 0.4),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
          axis.text.y = element_text(color="black",size=10),
          axis.text.x = element_text(margin = margin(t = -5),color="black",size=10),
          legend.position = "none",
          panel.spacing = unit(0,"lines"))+
    coord_cartesian() -> p1
  # ggsave(file= "mir_overall_plot.png", plot = p1, dpi = 800,width=8, height=4,limitsize = FALSE)
}

##曼哈顿圈图
a = map(list, ~ read.table(.,header = T))
data.frame(list)
xxxx <- data.frame(a[[7]][,c(1:3)],rowMeans(a[[7]][,c(7:ncol(a[[7]]))]))
lnc=data.frame(rowMeans(a[[2]][,c(7:ncol(a[[2]]))]))
pi=data.frame(rowMeans(a[[3]][,c(7:ncol(a[[3]]))]))
trf=data.frame(rowMeans(a[[8]][,c(7:ncol(a[[8]]))]))
mir=data.frame(rowMeans(a[[7]][,c(7:ncol(a[[7]]))]))
sn=data.frame(rowMeans(a[[6]][,c(7:ncol(a[[6]]))]))
sno=data.frame(rowMeans(a[[5]][,c(7:ncol(a[[5]]))]))
sca=data.frame(rowMeans(a[[4]][,c(7:ncol(a[[4]]))]))
#整理格式（强行合并所有表格，以na补缺失）
yyyy<-dplyr::bind_rows(xxxx,lnc,pi,trf,mir,sno,sn,sca)
#na全部变为0
for (i in 1:ncol(yyyy)){
  yyyy[,i][is.na(yyyy[,i])] <- 0
}

xxxx <- merge(data.frame(a[[1:length(a)]][,c(1,2,3)]),data.frame(Geneid=row.names(e[[lloop]][e[[lloop]]>0.1,]),
                                                                 e[[kloop]][e[[kloop]]>0.1,]),
              by="Geneid")
#按比例缩放到（0，1）
normfun<-function(data,ymin=0.01,ymax=0.99){
  xmax=max(data)
  xmin=min(data)
  y = (ymax-ymin)*(data-xmin)/(xmax-xmin) + ymin
  return(y)
}
apply(yyyy[,c(-1,-2,-3)],2,function(x)normfun(x)) |>  as.data.frame() ->yyyy[,c(4:ncol(yyyy))]
##准备染色体标签
for (i in 1:length(strsplit(yyyy$Chr,';','',fixed = T))) {
  yyyy[i,2] <- strsplit(yyyy$Chr,';','',fixed = T)[[i]][1]
}

yyyy$Chr<-gsub(";*","",yyyy$Chr)
yyyy=yyyy[yyyy$Chr%in%paste("chr",c(1:18,"X","Y"),sep=""),]
colnames(yyyy)<-c("Geneid","Chr","Start","a","b","c","d","e","f","g")
xxxx$Chr=gsub("chr","",xxxx$Chr)
#xxxx[,c(4:9)]=xxxx[xxxx[,c(-1,-2,-3)]*1e-5!=0,]
#any(xxxx[,c(-1,-2,-3)]*1e-5<=0)
library(CMplot)
CMplot(yyyy,plot.type="c",chr.labels=paste("Chr",c(1:18,"X","Y"),sep=""),r=0.2,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col=c("skyblue", "grey", "red"),file="jpg",
       memo="",dpi=800,file.output=TRUE,verbose=TRUE,LOG10 = TRUE)
CMplot(pig60K,plot.type="d",bin.size=1e6,chr.den.col=c("#2fa1dd", "white", "#f87669"),file="jpg",memo="",dpi=800,
       file.output=TRUE)
CMplot(cattle50K,plot.type="c",LOG10=FALSE,outward=TRUE,col=matrix(c("#4DAF4A",NA,NA,"dodgerblue4",
                                                                     "deepskyblue",NA,"dodgerblue1", "olivedrab3", "darkgoldenrod1"), nrow=3, byrow=TRUE),
       chr.labels=paste("Chr",c(1:29),sep=""),threshold=NULL,r=1.2,cir.chr.h=1.5,cir.legend.cex=0.5,
       cir.band=1,file="jpg", memo="",dpi=300,chr.den.col="black",file.output=FALSE,verbose=TRUE)->p3
##draw rowMeans by tissues to final columns
i = 1
aa = 7
for (i in 1:length(a)) {
  b <- a[[i]]
  for (aa in as.numeric(c(7:14,31:38))) {
    b <- cbind(b,rowMeans(b[,c(aa,aa+8,aa+8+8)]))
    colnames(b)[ncol(b)] <- as.character(aa)
  }
  a[[i]] <- b
}
colnames(b)

##rename every samples
lis <- c("denovo_miRNA","lncRNA","piRNA","scaRNA","snoRNA","snRNA","tRFs","miRNA")
for (i in 1:length(a)) {
  b <- as.data.frame(a[[i]])
  c <- data.frame(b[,c(1,(x+1):(x+16))]) 
  colnames(c)<- c("Geneid","Thymus_NOR","kidney_NOR","brain_NOR","liver_NOR","lung_NOR","intestines_NOR","heart_NOR","spleen_NOR",
                  "Thymus_HYP","kidney_HYP","brain_HYP","liver_HYP","lung_HYP","intestines_HYP","heart_HYP","spleen_HYP"
  )
  a[[i]] <- c
}

##计算每种RNA在每一个样本中的总数
lnc_countsums <- colSums(a[[1]][,-1])
miR_countsums <- colSums(a[[2]][,-1])
tRFs_countsums <- colSums(a[[3]][,-1])
tRNA_countsums <- colSums(a[[4]][,-1])
denovo_miR_countsums <- colSums(a[[5]][,-1])


######################3
colnames(b)[as.character(paste(strsplit(colnames(b),'_') [[1]][1],
                               strsplit(colnames(b),'_') [[1]][2],sep = "_"))==as.character(sample$sample_name)] <- as.character(
                                 sample$readable_name[as.character(paste(strsplit(colnames(b),'_') [[1]][1],
                                                                         strsplit(colnames(b),'_') [[1]][2],sep = "_"))==as.character(sample$sample_name)])



###### requirement #####
library(tidyverse)
library(magrittr)
library(DESeq2)
library(limma)
library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(readxl)


##### 1106重做S4 #####
rs.df.all.1
mir.df.all.1 = mir.df.all %>% mutate(class = "mir")
pir.df.all.1 = pir.df.all %>% mutate(class = "pir")
rny.df.all.1 = rny.df.all %>% mutate(class = "rny")
trf.df.all.1

n_df_all = rbind(mir.df.all.1, pir.df.all.1, rny.df.all.1, trf.df.all.1, rs.df.all.1)

meta_1016 = colnames(n_df_all)[-49] %>% limma::strsplit2("_") %>% data.frame() %>% dplyr::select(tissue = 1, condition = 2) %>% mutate(id = colnames(n_df_all)[-49])
meta_1016$condition %<>% factor(levels = c("NOR", "HYP"))

res_all = data.frame()
for (rna in c("rs", "mir", "pir", "rny", "trf")) {
  rna_tmp = rna
  df_tmp = n_df_all %>% filter(class == rna_tmp) %>% dplyr::select(-class)
  for (tissue in unique(meta_1016$tissue)) {
    tissue_tmp = tissue
    meta_tmp = meta_1016 %>% filter(tissue == tissue_tmp)
    df_tmp.1 = df_tmp %>% dplyr::select(meta_tmp$id)
    dds = DESeqDataSetFromMatrix(df_tmp.1, meta_tmp, ~condition) 
    dds$condition %<>% relevel(ref = "NOR")
    res = dds %>% DESeq() %>% results() %>% data.frame()
    res %<>% mutate(change = case_when(
      log2FoldChange > 0.25 & pvalue < 0.1 ~ "UP",
      log2FoldChange < -0.25 & pvalue < 0.1 ~ "DOWN",
      T ~ "NOT_sig"
    ))
    res %<>% mutate(rna = rna_tmp, tissue = tissue_tmp)
    res %<>% rownames_to_column("id")
    res_all = res %>% rbind(res_all, .)
  }
}

res_all%>% view()

##### ceRNA #####
filter_mir = res_all %>% filter(tissue == "lung", rna == "mir", change != "NOT_sig") %>% na.omit()
miRchcange = res_all %>% filter(tissue == "lung", rna == "mir", change != "NOT_sig") %>% na.omit()
filter_mir = filter_mir$id

ceRNAdb = read_xlsx("/Users/jiahao_kuang/Library/Mobile Documents/com~apple~CloudDocs/sncRNA landscape in HPH/ceRNA/rno_MTI.xlsx")

# 创建一个空的数据框来存储结果
result <- data.frame()

# 循环遍历 a 中的每个元素
for (elem in filter_mir) {
  temp <- ceRNAdb[grep(elem, ceRNAdb$miRNA), ]
  result <- rbind(result, temp)
}
result <- unique(result)
result %<>% filter(miRNA %in% filter_mir)


##### DEmRNA 
luo.mRNA = read_xls("/Users/jiahao_kuang/Dropbox/邝嘉浩/sncRNA in tissue/iscience/Revision/ceRNA/6.罗灵杰-大鼠肺组织-2021-9-16.xls")
luo.meta = colnames(luo.mRNA)[1:17] %>% gsub("(\\D+)(\\d+).*", "\\1", .) %>% 
  data.frame(id = colnames(luo.mRNA)[1:17], group = .) %>% 
  filter(group %in% c("N", "H")) %>% 
  data.frame(row.names = 1)

luo.mRNA = luo.mRNA %>%  
  column_to_rownames("gene_id") %>% 
  dplyr::select(rownames(luo.meta))

dds <-  DESeqDataSetFromMatrix(countData = luo.mRNA, colData = luo.meta, design = ~group)
dds$group %<>% relevel(ref = "N")
res <-  dds %>% 
  DESeq() %>% 
  results() %>% 
  data.frame() %>% na.omit() %>% 
  mutate(change = case_when(log2FoldChange > 0.25 & pvalue < 0.1 ~ "UP", 
                            log2FoldChange < -0.25 & pvalue < 0.1 ~ "DOWN",
                            T ~ "NOT_sig"))

res <- rownames(res) %>% 
  clusterProfiler::bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Rn.eg.db, drop = F) %>%
  dplyr::distinct(ENSEMBL, .keep_all = T) %>% 
  data.frame(row.names = 1) %>%  merge(res, ., by = "row.names") %>% 
  data.frame(row.names = 1) %>% 
  na.omit()

res.1 <- res %>% filter(change !="NOT_sig", SYMBOL %in% result$`Target Gene`) 

network_data <- result %>% filter(`Target Gene` %in% res.1$SYMBOL) %>% dplyr::select("miRNA", mRNA = "Target Gene")


# 找到miRNA变化的方向
for (i in 1:nrow(network_data)) {
  network_data$mir_change[i] = miRchcange[miRchcange$id %in% network_data[i,1],]$change
}

# 找到mRNA变化的方向
for (i in 1:nrow(network_data)) {
  network_data$mrna_change[i] = res.1[res.1$SYMBOL %in% network_data[i,2],]$change
}

unique(network_data)

network_data %<>% filter(mrna_change != "NOT_sig")
network_data = unique(network_data)
network_data %<>% filter(network_data$mir_change != network_data$mrna_change)
network_data$miRNA %<>% gsub(pattern = "rno-", replacement = "")


# miRchcange$id %in% unique(network_data$miRNA)

setwd("/Users/jiahao_kuang/Library/Mobile Documents/com~apple~CloudDocs/sncRNA landscape in HPH/ceRNA")
write.table(network_data, file = "cenods.1.txt", row.names = F,sep = "\t",quote = F)

data.frame(rbind(data.frame(rna = network_data$miRNA, type = "mirna"), data.frame(rna = network_data$mRNA, type = "mrna"))) %>% write.table(., file = "attribute.txt", row.names = F,sep = "\t",quote = F)
data.frame(rna = network_data$mRNA, type = "mrna")

# PPI
gene = res %>% filter(change != "Not_sig") %>% rownames() %>% clusterProfiler::bitr(fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
network_data

gene = network_data$mRNA %>% unique() %>% clusterProfiler::bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)

library(STRINGdb)
?STRINGdb
string_db <- STRINGdb$new(version="11.0b", species=10116, 
                          score_threshold=500, input_directory="")
data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                      removeUnmappedRows = TRUE)
# string_db$plot_network( data_mapped$STRING_id )
string_db$plot_network(data_mapped$STRING_id)

# 使用get_interactions获取蛋白互作信息，以用于后续可视化。
hit <- data_mapped$STRING_id
info <- string_db$get_interactions(hit)

# 使用igraph和ggraph可视化蛋白互作网络图
# 转换stringID为Symbol，只取前两列和最后一列
links <- info %>%
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score)
# 节点数据
nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
# 创建网络图
# 根据links和nodes创建
net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
# 添加一些参数信息用于后续绘图
# V和E是igraph包的函数，分别用于修改网络图的节点（nodes）和连线(links)
igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
igraph::V(net)$size <- igraph::degree(net)/5 #
igraph::E(net)$width <- igraph::E(net)$weight/10


##### 1105 王老师要求的数据#####
### RNA的种类
rownames(df.all) = paste(df.all$class, rownames(df.all), sep = "_")

df.all.t = df.all %>% filter(class != "lnc") %>% dplyr::select(-"class") %>% t()

df.all.t = rownames(df.all.t) %>% limma::strsplit2("_") %>% data.frame(., id = rownames(df.all.t)) %>% dplyr::select(tissue = 1, 
                                                                                                                     condition = 2, 
                                                                                                                     id = id) %>% cbind(., df.all.t)
df.all.t.select = df.all.t[,-3] %>% group_by(tissue 
                                             # condition
                                             ) %>% summarise_all(mean)

# selector = function(x){ifelse(x > 20, x, NA)}
colname = paste(df.all.t.select$tissue, df.all.t.select$condition, sep = "_")
df.all.t.select = df.all.t.select[,-c(1:2)] %>% apply(., 1, function(x){ifelse(x > 20, x, NA)}) %>% data.frame()
colnames(df.all.t.select) = colname

# 一个统计每列NA数量的函数
# na.count = function(x){sum(is.na(x))}
na_counter = function(x){ifelse(is.na(x), 0, 1)}
na_counts = df.all.t.select %>% apply(., 2, function(x){ifelse(is.na(x), 0, 1)}) %>% data.frame() %>% mutate(class = limma::strsplit2(rownames(.), "_")[,1])

# 每种RNA的种类数量
na_counts %>% group_by(class) %>% summarise_all(sum)

##### 每种RNA变化的数量 #####
df_orign = df.all

df_all = readxl::read_xlsx("/Users/jiahao_kuang/Dropbox/邝嘉浩/sncRNA in tissue/iscience/Table S1 ExprSet gengerated in this study.xlsx")[,]
colnames(df_all) = df_all[1,]
df_all = df_all[-1,]
df.all = df_all %>% column_to_rownames("id")
df_all = df_all[!(df_all$id%in%rownames(lnc.df)),]

case_when(df_all$id %in% mir.df)




# meta_1105 = colnames(df.all)[-49] %>% limma::strsplit2("_") %>% data.frame() %>% dplyr::select(tissue = 1, condition = 2) %>% mutate(id = colnames(df.all)[-49])
meta_1105 = colnames(df.all) %>% limma::strsplit2("_") %>% data.frame() %>% dplyr::select(tissue = 1, condition = 2) %>% mutate(id = colnames(df.all)[-49])


library(DESeq2)
df.all = mir.df.all
res = data.frame()
for (tissue in unique(meta_1105$tissue)) {
  tissue_tmp = tissue
  meta_tmp = meta_1105 %>% filter(tissue == tissue_tmp)
  meta_tmp$condition %<>% factor(levels = c("NOR", "HYP"))
  df_tmp = df.all %>% dplyr::select(meta_tmp$id) 
  rowname_tmp = rownames(df_tmp)
  df_tmp %<>% apply(., 2, as.numeric) %>% data.frame()
  rownames(df_tmp) = rowname_tmp
  dds = DESeqDataSetFromMatrix(df_tmp, meta_tmp, ~condition)
  dds$condition %<>% relevel(ref = "NOR")
  res_tmp = DESeq(dds) %>% results() %>% data.frame() %>% rownames_to_column("id") %>% na.omit()
  res_tmp %<>% mutate(change = case_when(
    baseMean > 0  & log2FoldChange > 0.25  &pvalue  < 0.1 ~ "UP",
    baseMean > 0 &log2FoldChange < -0.25 & pvalue < 0.1 ~ "DOWN",
     # pvalue  < 0.05 ~ "UP",
     # pvalue < 0.05 ~ "DOWN", 
    T ~ "Not_sig"))
  res_tmp$tissue = tissue_tmp
  res = rbind(res, res_tmp)
}

res %>% view()


res %>% mutate(class = limma::strsplit2(rownames(res),"_")[,1]) %>% filter(change != "Not_sig") %>% group_by(tissue, class) %>% summarise(count = n()) %>% filter(class != "lnc") %>% reshape2::dcast(tissue ~ class, value.var = "count")


#### 显著变化的RNA数量 ####
tables5 = readxl::read_xlsx("/Users/jiahao_kuang/Library/Mobile Documents/com~apple~CloudDocs/sncRNA in tissue/iscience/Supplementary Table/Table S4 - new Differential Analysis Result.xlsx")[,]
colnames(tables5) = tables5[1,]
tables5 = tables5[-1,]

tables5 %<>% filter(rna != "lnc")
unique_class_type = tables5$rna %>% unique()
unique_tissue = tables5$tissue %>% unique()

mir_sig = (tables5 %>% filter(rna == "mir" & change != "Not_sig"))$tissue %>% table() %>% data.frame() %>% dplyr::rename(tissue = 1, mir = 2)
trf_sig = (tables5 %>% filter(rna == "trf" & change != "Not_sig"))$tissue %>% table() %>% data.frame() %>% dplyr::rename(tissue = 1, trf = 2)
pir_sig = (tables5 %>% filter(rna == "pir" & change != "Not_sig"))$tissue %>% table() %>% data.frame() %>% dplyr::rename(tissue = 1, pir = 2)
rny_sig = (tables5 %>% filter(rna == "rny" & change != "Not_sig"))$tissue %>% table() %>% data.frame() %>% dplyr::rename(tissue = 1, rny = 2)
rs_sig = (tables5 %>% filter(rna == "rs" & change != "Not_sig"))$tissue %>% table() %>% data.frame() %>% dplyr::rename(tissue = 1, rs = 2)

sig_count = Reduce(function(x, y) merge(x, y, by = "tissue", all = TRUE), list(mir_sig, trf_sig, pir_sig, rny_sig, rs_sig))



for (tissue in unique_tissue) {
  (tables5 %>% filter(rna == "mir" & change != "Not_sig"))$tissue %>% table()
}
  
}
tables5 %>% filter(rna == "mir" & change != "Not_sig" & tissue == "kidney") %>% view()


##### S4错误后的差异分析 #####
setwd("/Users/jiahao_kuang/Library/Mobile Documents/com~apple~CloudDocs/sncRNA landscape in HPH/")
load("20230224_final_SNC_image.RData")
##### miRNA Count Cor #####
library(magrittr)
library(tidyverse)
library(limma)
library(DESeq2)
library(plyr)
library(dplyr)
library(reshape2)

# miRNA
setwd("/Users/jiahao_kuang/Dropbox/邝嘉浩/Rat_plasma_miRNA_20230830/miRNA/")
miRlist <- list.files()

miRfile = NULL
x=1
for (i in miRlist) {
  miRfile[[x]] <- read.table(i) %>% data.frame()
  i %<>% gsub(pattern = "_miRNA_cal.txt", replacement = "")
  colnames(miRfile[[x]]) <- c(i, "id")
  x = x + 1
}

df.mir <- plyr::join_all(miRfile, by = "id", type = "full")
df.mir %<>% column_to_rownames("id")
df.mir[is.na(df.mir)] = 0

meta <- readxl::read_excel("../r_大鼠PAH模型血浆cfRNA测序信息表.xlsx")       

i=1
for (i in 1:ncol(df.mir)) {
  colnames(df.mir)[i] <- meta[meta$ID%in%colnames(df.mir)[i],]$造模编号
}

plasma_coldata <- data.frame(id = colnames(df.mir), group = toupper(limma::strsplit2(colnames(df.mir),split = "-")[,1]), tissue = "plasma")
tissue_coldata <- colnames(all.mir.df) %>% 
  strsplit2("_") %>% 
  data.frame() %>% 
  dplyr::rename(tissue = 1, group = 2) %>% 
  dplyr::select(2,1) %>% 
  data.frame(id = colnames(all.mir.df), .)
tissue_plasma.meta <- rbind(plasma_coldata, tissue_coldata)
tissue_plasma.meta$group %<>% gsub(pattern = "NC", replacement = "NOR")
tissue_plasma.meta.nor <- tissue_plasma.meta %>% filter(group %in% c("NOR"))

tissue_plasma.mir.df = merge(df.mir, all.mir.df, by = "row.names") %>% column_to_rownames("Row.names") 
tissue_plasma.mir.df %>% dim()
# ExprMat with counts
tissue_plasma.mir.df.1 <- tissue_plasma.mir.df[apply(tissue_plasma.mir.df, 1, sum) > 0 , ] %>% 
  # apply(. ,2, function(x) { log10(x/sum(x)*1000000) }) %>% 
  t() %>% cbind(tissue_plasma.meta) %>% dplyr::select(-"id")

# ExprMat with log10RPM
tissue_plasma.mir.df <- tissue_plasma.mir.df[apply(tissue_plasma.mir.df, 1, sum) > 0 , ] %>% 
  apply(. ,2, function(x) { log10(x/sum(x)*1000000) }) %>% 
  t() %>% cbind(tissue_plasma.meta) %>% dplyr::select(-"id")

library(GGally)
# ggpairs(tissue_plasma.mir.sum.nor)
# 
# nor.idx = tissue_plasma.meta %>% filter(group == "NOR")
# hyp.idx = tissue_plasma.meta %>% filter(group == "HYP")
# 
# data.frame(dplyr::select(tissue_plasma.mir.df, nor.idx$id))

tissue_plasma.meta
tissue_plasma.mir.sum = tissue_plasma.mir.df %>% group_by(group, tissue) %>% summarise_all(mean) %>% data.frame()
tissue_plasma.mir.sum[tissue_plasma.mir.sum == -Inf] = 0
rownames(tissue_plasma.mir.sum) = paste(tissue_plasma.mir.sum$group, tissue_plasma.mir.sum$tissue, sep = "_") 
tissue_plasma.mir.sum %<>% 
  dplyr::select(-c("tissue", "group")) %>% 
  t()

# tissue_plasma.mir.sum[tissue_plasma.mir.sum == -Inf] = 0
tissue_plasma.mir.sum %<>% data.frame()

meta.tmp <- colnames(tissue_plasma.mir.sum) %>% 
  limma::strsplit2("_") %>% data.frame() %>% dplyr::rename(group = 1, tissue = 2) %>% data.frame(id = colnames(tissue_plasma.mir.sum))

hyp.idx = meta.tmp %>% filter(group == "HYP")
nor.idx = meta.tmp %>% filter(group == "NOR")

tissue_plasma.mir.sum.nor = tissue_plasma.mir.sum[,nor.idx$id]
colnames(tissue_plasma.mir.sum.nor) = nor.idx$tissue
tissue_plasma.mir.sum.nor %<>% mutate(group = "NOR")

tissue_plasma.mir.sum.hyp = tissue_plasma.mir.sum[,hyp.idx$id]
colnames(tissue_plasma.mir.sum.hyp) = hyp.idx$tissue
tissue_plasma.mir.sum.hyp %<>% mutate(group = "HYP")

cor = rbind(tissue_plasma.mir.sum.nor, tissue_plasma.mir.sum.hyp)

library(GGally)
GGscatterPlot <- function(data, mapping, ..., 
                          method = "pearson") {
  
  #Get correlation coefficient
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor <- cor(x, y, method = method, use="pairwise.complete.obs")
  #Assemble data frame
  df <- data.frame(x = x, y = y)
  df <- na.omit(df)
  # PCA
  nonNull <- x!=0 & y!=0
  dfpc <- prcomp(~x+y, df[nonNull,])
  df$cols <- predict(dfpc, df)[,1]
  # Define the direction of color range based on PC1 orientation:
  dfsum <- x+y
  colDirection <- ifelse(dfsum[which.max(df$cols)] < 
                           dfsum[which.min(df$cols)],
                         1,
                         -1)
  #Get 2D density for alpha
  dens2D <- MASS::kde2d(df$x, df$y)
  df$density <- fields::interp.surface(dens2D ,df[,c("x", "y")])
  
  if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
  }
  #Prepare plot
  pp <- ggplot(df, aes(x=x, y=y, alpha = 1/density, color = cols)) +
    ggplot2::geom_point(shape=16, show.legend = FALSE) +
    ggplot2::scale_color_viridis_c(direction = colDirection) +
    ggplot2::scale_alpha(range = c(.05, .6)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = round(cor, digits = 3)),
      mapping = ggplot2::aes(x = xlabel,
                             y = ylabel,
                             label = lab),
      hjust = 0, vjust = 1,
      size = 3, fontface = "bold",
      inherit.aes = FALSE # do not inherit anything from the ...
    ) +
    theme_bw()
  return(pp)
}
colnames(cor) %<>% factor(levels = c("plasma", colnames(cor)[-7]))
ggpairs(cor, ggplot2::aes(color=group)) +
  # scale_color_manual(values = c("#f0a1a8","#4994c4")) +
  # scale_fill_manual(values = c("#f0a1a8","#4994c4")) +
  theme(axis.text = element_text(colour = "black", size = 11),
        strip.background = element_rect(fill = "#d63d2d"),
        strip.text = element_text(colour = "white", size = 12,
                                  face = "bold")) -> mir.plasmaCor
GGally::ggpairs(cor,
                # 1:4,
                mapping = aes(color = group),
                lower = list(continuous = wrap(GGscatterPlot, method="pearson")),
                upper = list(continuous = wrap(ggally_cor, align_percent = 0.8))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text =  element_text(color='black'),
        strip.background = element_rect(fill = "#4994c4"),
        strip.text = element_text(colour = "white", size = 10,
                                  face = "bold")) -> mir.plasmaCor
mir.plasmaCor

setwd("/Users/jiahao_kuang/Dropbox/邝嘉浩/Rat_plasma_miRNA_20230830/miRNA")
ggsave(filename = "miR_log10rpm_Cor.pdf", units = "cm", width = 26, height = 25)

##### tRNA Count Cor #####
plasma.trf = read.table("/Users/jiahao_kuang/Library/Mobile Documents/com~apple~CloudDocs/sncRNA landscape in HPH/cor with luo's dat/tRNA_plasma_ExprMat.txt", row.names = 1, header = T)
plasma.trf %<>% dplyr::select(-c(1:5))
i=1
for (i in 1:ncol(plasma.trf)) {
  colnames(plasma.trf)[i] %<>% gsub(pattern = "_cutada_trim_rmrsRNA_rmysRNA_R1.fastq.sam", replacement = "") %>% gsub(pattern = "\\.", replacement = "-")
  colnames(plasma.trf)[i] <- meta[meta$ID%in%colnames(plasma.trf)[i],]$造模编号
}


rownames(plasma.trf)


merge(plasma.trf, trf.df, by = "row.names")

plasma_coldata <- data.frame(id = colnames(plasma.trf), group = toupper(limma::strsplit2(colnames(df.mir),split = "-")[,1]), tissue = "plasma")
tissue_coldata <- colnames(trf.df) %>% 
  strsplit2("_") %>% 
  data.frame() %>% 
  dplyr::rename(tissue = 1, group = 2) %>% 
  dplyr::select(2,1) %>% 
  data.frame(id = colnames(trf.df), .)
tissue_plasma.meta <- rbind(plasma_coldata, tissue_coldata)
tissue_plasma.meta$group %<>% gsub(pattern = "NC", replacement = "NOR")
tissue_plasma.meta.nor <- tissue_plasma.meta %>% filter(group %in% c("NOR"))
# cpm <- apply(count ,2, function(x) { x/sum(x)*1000000 })
# countData <- count[apply(count, 1, sum) > 10 , ]
tissue_plasma.trf.df = merge(plasma.trf, trf.df, by = "row.names") %>% column_to_rownames("Row.names") 
tissue_plasma.trf.df.1 <- tissue_plasma.trf.df[apply(tissue_plasma.trf.df, 1, sum) > 15 , ] %>% 
  # apply(. ,2, function(x) { log10(x/sum(x)*1000000) }) %>% 
  t() %>% cbind(tissue_plasma.meta) %>% dplyr::select(-"id")
tissue_plasma.trf.df <- tissue_plasma.trf.df[apply(tissue_plasma.trf.df, 1, sum) > 15 , ] %>% 
  apply(. ,2, function(x) { log10(x/sum(x)*1000000) }) %>% 
  t() %>% cbind(tissue_plasma.meta) %>% dplyr::select(-"id")

library(GGally)
tissue_plasma.meta
tissue_plasma.trf.sum = tissue_plasma.trf.df %>% group_by(group, tissue) %>% summarise_all(mean) %>% data.frame()
tissue_plasma.trf.sum[tissue_plasma.trf.sum == -Inf] = 0
rownames(tissue_plasma.trf.sum) = paste(tissue_plasma.trf.sum$group, tissue_plasma.trf.sum$tissue, sep = "_") 
tissue_plasma.trf.sum %<>% 
  dplyr::select(-c("tissue", "group")) %>% 
  t()

# tissue_plasma.trf.sum[tissue_plasma.trf.sum == -Inf] = 0
tissue_plasma.trf.sum %<>% data.frame()

meta.tmp <- colnames(tissue_plasma.trf.sum) %>% 
  limma::strsplit2("_") %>% data.frame() %>% dplyr::rename(group = 1, tissue = 2) %>% data.frame(id = colnames(tissue_plasma.trf.sum))

hyp.idx = meta.tmp %>% filter(group == "HYP")
nor.idx = meta.tmp %>% filter(group == "NOR")

tissue_plasma.trf.sum.nor = tissue_plasma.trf.sum[,nor.idx$id]
colnames(tissue_plasma.trf.sum.nor) = nor.idx$tissue
tissue_plasma.trf.sum.nor %<>% mutate(group = "NOR")

tissue_plasma.trf.sum.hyp = tissue_plasma.trf.sum[,hyp.idx$id]
colnames(tissue_plasma.trf.sum.hyp) = hyp.idx$tissue
tissue_plasma.trf.sum.hyp %<>% mutate(group = "HYP")

cor = rbind(tissue_plasma.trf.sum.nor, tissue_plasma.trf.sum.hyp)

library(GGally)
colnames(cor) %<>% factor(levels = c("plasma", colnames(cor)[-7]))
ggpairs(cor, ggplot2::aes(color=group)) +
  # scale_color_manual(values = c("#f0a1a8","#4994c4")) +
  # scale_fill_manual(values = c("#f0a1a8","#4994c4")) +
  theme(axis.text = element_text(colour = "black", size = 11),
        strip.background = element_rect(fill = "#d63d2d"),
        strip.text = element_text(colour = "white", size = 12,
                                  face = "bold")) -> trf.plasmaCor

cor.1 = cor
tmp.1 = cor %>% dplyr::select(-"group")
# 
# cor.1[apply(tmp.1, 1, sum) > 1 , ] 
# 
# 
# tmp.1[rowSums(tmp.1) == 0]
# 
# tmp.1 = tmp.1+0.1
# cor.1 = cbind(tmp.1, group = cor.1$group)
cor.1[apply(tmp.1, 1, sum) > 10 , ] %>% 
  GGally::ggpairs(.,
                  # 1:4,
                  mapping = aes(color = group),
                  lower = list(continuous = wrap(GGscatterPlot, method="pearson")),
                  upper = list(continuous = wrap(ggally_cor, align_percent = 0.8))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text =  element_text(color='black'),
        strip.background = element_rect(fill = "#d63d2d"),
        strip.text = element_text(colour = "white", size = 10,
                                  face = "bold")) -> trf.plasmaCor
setwd("/Users/jiahao_kuang/Dropbox/邝嘉浩/Rat_plasma_miRNA_20230830/miRNA")
ggsave(trf.plasmaCor, filename = "tRNA_log10rpm_Cor.pdf", units = "cm", width = 26, height = 25)
ggsave(trf.plasmaCor, filename = "tRNA_log10rpm_Cor.png", dpi = 1000, units = "cm", width = 26, height = 25)

##### miRNA log2FC Cor #####
tissue_plasma.mir.ColDat = tissue_plasma.mir.df.1 %>% dplyr::select("group", "tissue")
tissue_plasma.mir.CountDat = tissue_plasma.mir.df.1 %>% dplyr::select(-c("group", "tissue")) %>% t()

colnames(tissue_plasma.mir.CountDat) == rownames(tissue_plasma.mir.ColDat)

# 血浆的miRNA做差异分析
coldata = tissue_plasma.mir.ColDat %>% rownames_to_column("id") %>% filter(tissue == "plasma") 
plasma_mir_res = NULL
x = 1
## DESEQ
# i = "HYP"
for (i in c("HYP")) {
  coldata.tmp <- coldata %>% filter(group %in% c("NOR", i))
  coldata.tmp$group <- factor(coldata.tmp$group, levels = c("NOR", i))
  countdata.tmp <- tissue_plasma.mir.CountDat[,coldata.tmp$id]
  # countdat.tem.1 = countdata.tmp
  # countdata.tmp.1 = countdata.tmp %>% data.frame() %>% dplyr::select("NC.4", "NC.5", "NC.9", "NC.2.7", "Hyp.1", "Hyp.4", "Hyp.2.1","Hyp.7")
  # colnames(countdata.tmp.1) %<>% gsub(pattern = "\\.", replacement = "-")
  # coldata.tmp.1 = coldata.tmp %>% filter(id %in% colnames(countdata.tmp.1))
  dds <- DESeqDataSetFromMatrix(countData = countdata.tmp, colData = coldata.tmp, design = ~group) 
  dds$group <- relevel(dds$group, ref = "NOR") 
  res <- DESeq(dds) %>% results() %>% data.frame(id = rownames(.), ., tissue = "plasma")
  res %<>% na.omit() %>% mutate(change =
                                  case_when(log2FoldChange > .25 & pvalue < 0.05 ~ "Up_regulated", 
                                            log2FoldChange < -.25 & pvalue < 0.05 ~ "Down_regulated",
                                            T ~ "Not_sig"
                                  )
  )
  plasma_mir_res[[x]] <- res
  names(plasma_mir_res)[x] <- paste(i, "Versus_NC", sep = "_")
  x = x + 1
}

# 火山图
if(T){
  
  # plasma_mir_res$HYP_Versus_NC[ce_target$SYMBOL,]
  
  plasma_mir_res$HYP_Versus_NC %>% 
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(pvalue), color = change)) + 
    geom_point(size = 1) + 
    scale_color_manual(values = c("#699ed4","grey90","#ef8183")) + 
    scale_fill_manual(values = c("#699ed4","#ef8183")) + 
    xlim(-2, 2) +
    # scale_x_continuous(limits = c(-8, 8),breaks = seq(-8, 8, by = 4)) +
    # scale_y_continuous(limits = c(0, 150),breaks = seq(0, 150, by = 50)) +
    geom_vline(xintercept = .25, color = "#ef8183") +
    geom_vline(xintercept = -.25, color = "#699ed4") +
    geom_hline(yintercept = -log10(0.05), color = "grey50") +
    # theme_classic()+
    theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = .75),
          panel.background = element_blank(), 
          axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold.italic")) -> vol.miRplasma
  vol.miRplasma
  
  tissue_plasma.mir.CountDat[rownames(tissue_plasma.mir.CountDat)%in% (plasma_mir_res$HYP_Versus_NC %>% filter(change!="Not_sig"))$id,coldata.tmp$id] %>% pheatmap::pheatmap(scale = "row")
  
  setwd("/Users/jiahao_kuang/Dropbox/邝嘉浩/sncRNA in tissue/iscience/Revision/cor with plasma/Plasma Append Figures")
  ggsave(plot = vol.miRplasma, filename = "DEmiR_plasma_Volcano.pdf", units = "cm", width = 10, height = 8)
  
  ce_target$SYMBOL
  
  
  plasma.sig = plasma_mir_res$HYP_Versus_NC %>% filter(change != "Not_sig")
  plasma.idx = coldata.tmp.1 %>% filter(tissue == "plasma" & group %in% c("NOR", "HYP"))
  plasma.idx %<>% column_to_rownames("id")
  plasma_hmDat = countdata.tmp.1[rownames(plasma.sig),rownames(plasma.idx)]
  plasma_hmDat %>% pheatmap::pheatmap(scale = "row", cellwidth = 13, cellheight = 13, border_color = NA, annotation = , color = colorRampPalette (c("#699ed4","white","#ef8183")) (100))
  
  
  (plasma_hmDat+1) %>% apply(. ,2, function(x) { log10(x/sum(x)*1000000) }) %>% 
    t() %>% scale() %>% data.frame() %>% merge(plasma.idx, by = "row.names") %>% 
    dplyr::select(-"tissue") %>% melt(id = c("group", "Row.names")) -> hm
  
  library(magrittr)
  library(reshape2)
  library(viridis)
  library(scales)
  
  getwd()
  setwd("/Users/jiahao_kuang/Dropbox/邝嘉浩/sncRNA in tissue/iscience/Revision/cor with plasma/Plasma Append Figures")
  # save.image("1108_SNC.RDA")
  # load("1108_SNC.RDA")
  
  library(ggplot2)
  library(tidyverse)
  library(magrittr)
  library(scales)
  hm$value %<>% as.numeric()
  hm$variable %<>% as.character() %>% gsub(pattern = "\\.", replacement = "-")
  hm$Row.names %<>% as.character()
  
  
  hm %>% data.frame() %>% 
    ggplot(data = ., aes(x = Row.names, y = variable, fill = value))+
    geom_tile(color = "black") +
    labs(x = NULL, y = NULL) +
    coord_fixed() + 
    # geom_text(aes(label = lab),color = "white") +
    scale_fill_gradientn(
      colours = c("#003366", "white", "#990033"),
      # values = rescale(c(-10, 0, 10)),
      # limits = c(-15, 15),
      oob = squish
    )+
    theme_minimal() +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    )+
    scale_y_discrete(position = "right") -> hm.feature
  
  library(ggtree)
  library(cowplot)
  library(aplot)
  
  gg <- hclust(dist(plasma_hmDat))    #对行聚类
  
  zz <- hclust(dist(t(plasma_hmDat)))    #对列聚类
  
  h <- ggtree(gg,layout = "rectangular",branch.length = "none") # 绘制列聚类树
  
  v <- ggtree(zz)+layout_dendrogram() # 绘制行聚类树
  
  library(aplot)
  hm.feature %>% insert_left(h,width = 0.3) 
 
  tsne.dat = (plasma_hmDat+1) %>% apply(. ,2, function(x) { log10(x/sum(x)*1000000) }) %>% 
    t() %>% scale() %>% data.frame() %>% merge(plasma.idx, by = "row.names") %>% column_to_rownames("Row.names") 
  
  tsne.dat %>% dim()
  library(Rtsne)
  tsne_out = Rtsne(
    tsne.dat[,c(1:8)],
    dims = 2,
    pca = T,
    max_iter = 1000,
    theta = 0.4,
    perplexity = 1,
    # verbose = F
  ) # 进行t-SNE降维分析 
  
  color = tsne.dat$group
  library(ggplot2)
  tsne_result = as.data.frame(tsne_out$Y)
  colnames(tsne_result) = c("tSNE1","tSNE2")
  ggplot(tsne_result,aes(tSNE1,tSNE2,color=color)) +
    geom_point() + 
    stat_ellipse() +
    coord_fixed(ratio = 2) +
    theme_bw()
  getwd()
  ggsave(filename = "DEmiR_tSNE.pdf", units = "cm", width = 8, height = 10)
  
  
  library(ggplot2)
  hm %<>% dplyr::rename(id=Row.names)
  hm$value %<>% as.numeric()
  hm$variable %<>% gsub(pattern = "\\.", replacement = "-")
  hm$id %<>% as.character()
  
  
  library(magrittr)
  library(reshape2)
  library(viridis)
  library(scales)
  hm %>% 
  ggplot(.,aes(x=id, y = variable, fill = value)) +
    geom_tile(color = "black") +
    labs(x = NULL, y = NULL) +
    # coord_fixed() + 
    # geom_text(aes(label = lab),color = "white") +
    scale_fill_gradientn(
      colours = c("#003366", "white", "#990033"),
      # values = rescale(c(-10, 0, 10)),
      limits = c(-.25, .25),
      oob = squish
    )+
    theme_minimal() +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    )+
    scale_y_discrete(position = "right") -> heatmap
  heatmap
library(ggtree)
gg <- hclust(dist(plasma_hmDat))    #对行聚类
zz <- hclust(dist(t(plasma_hmDat)))    #对列聚类
h <- ggtree(gg) # 绘制列聚类树
# v <- ggtree(zz)+layout_dendrogram() # 绘制行聚类树
library(aplot)
heatmap %>% 
  # insert_top(group,height = 0.02) %>% 
  # insert_top(v,height = 0.1) %>% 
  insert_left(h,width = 0.3)


group <- data.frame(
  id = hm$id %>% unique() ,
  x = ((hm$id %>% unique() %>% strsplit2("-") %>% data.frame() %>% dplyr::select(1))$X1 %>% factor(levels = c("NC", "Hyp"))),
  y = 1
  ) %>% 
# group <- data.frame(id = factor(c(rep("NOR", 14), rep("HYP", 11)), levels = c("NOR", "HYP")),
#            x = c(1:25),
#            y = 1
#            ) %>% 
# group <- rep(c("A","B"),each = 3) %>% data.frame(x = c(1:6),y = rep(1,6)) %>%
  ggplot()+
  geom_tile(aes(id,y,fill = x))+
  scale_fill_discrete(label = c("control","treat"))+
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank())



  
}

# 组织的miRNA做差异分析
coldata = tissue_plasma.mir.ColDat %>% rownames_to_column("id") %>% filter(tissue != "plasma") 
tissue_mir_res = NULL
# x = 1
## DESEQ
for (i in unique(coldata$tissue)) {
  coldata.tmp <- coldata %>% filter(tissue == i)
  # coldata.tmp <- coldata %>% filter(group %in% c("NOR", i))
  coldata.tmp$group <- factor(coldata.tmp$group, levels = c("NOR", "HYP"))
  countdata.tmp <- tissue_plasma.mir.CountDat[,coldata.tmp$id]
  dds <- DESeqDataSetFromMatrix(countData = countdata.tmp, colData = coldata.tmp, design = ~group) 
  dds$group <- relevel(dds$group, ref = "NOR") 
  res <- DESeq(dds) %>% results() %>% data.frame(id = rownames(.), ., tissue = i)
  res %<>% na.omit() %>% mutate(change =
                                  case_when(log2FoldChange > 0.25 & pvalue < 0.05 ~ "Up_regulated", 
                                            log2FoldChange < -0.25 & pvalue < 0.05 ~ "Down_regulated",
                                            T ~ "Not_sig"
                                  )
  )
  tissue_mir_res <- rbind(tissue_mir_res, res)
  # names(tissue_mir_res)[x] <- paste(i, "Versus_NC", sep = "_")
  # x = x + 1
}

# dcast(difres.trf, id~tissue, value.var = "log2FoldChange")
tissue_mir_res.wide = dcast(tissue_mir_res, id~tissue, value.var = "log2FoldChange")
plasma_mir_res.wide = plasma_mir_res$HYP_Versus_NC %>% dplyr::select("id", plasma = "log2FoldChange")

GGscatterPlot <- function(data, mapping, ..., 
                          method = "pearson") {
  
  #Get correlation coefficient
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor <- cor(x, y, method = method, use="pairwise.complete.obs")
  #Assemble data frame
  df <- data.frame(x = x, y = y)
  df <- na.omit(df)
  # PCA
  nonNull <- x!=0 & y!=0
  dfpc <- prcomp(~x+y, df[nonNull,])
  df$cols <- predict(dfpc, df)[,1]
  # Define the direction of color range based on PC1 orientation:
  dfsum <- x+y
  colDirection <- ifelse(dfsum[which.max(df$cols)] < 
                           dfsum[which.min(df$cols)],
                         1,
                         -1)
  #Get 2D density for alpha
  dens2D <- MASS::kde2d(df$x, df$y)
  df$density <- fields::interp.surface(dens2D ,df[,c("x", "y")])
  
  if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
  }
  #Prepare plot
  pp <- ggplot(df, aes(x=x, y=y, alpha = 1/density, color = cols)) +
    ggplot2::geom_point(shape=16, show.legend = FALSE) +
    ggplot2::scale_color_viridis_c(direction = colDirection) +
    ggplot2::scale_alpha(range = c(.05, .6)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = round(cor, digits = 3)),
      mapping = ggplot2::aes(x = xlabel,
                             y = ylabel,
                             label = lab),
      hjust = 0, vjust = 1,
      size = 3, fontface = "bold",
      inherit.aes = FALSE # do not inherit anything from the ...
    ) +
    theme_bw()
  return(pp)
}

cor.dat <- full_join(tissue_mir_res.wide, plasma_mir_res.wide, by = "id") %>%
  dplyr::select(-"id")
cor.dat[is.na(cor.dat)] = 0

GGally::ggpairs(cor.dat,
                # 1:4,
                lower = list(continuous = wrap(GGscatterPlot, method="pearson")),
                upper = list(continuous = wrap(ggally_cor, align_percent = 0.8))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text =  element_text(color='black'),
        strip.background = element_rect(fill = "#4994c4"),
        strip.text = element_text(colour = "white", size = 10,
                                  face = "bold"))

getwd()
# ggsave(filename = "miR_log2FC.pdf", units = "cm", height = 19, width = 20)

##### tRNA log2FC Cor #####
# 血浆的tRNA做差异分析
tissue_plasma.trf.df

tissue_plasma.trf.ColDat = tissue_plasma.trf.df.1 %>% dplyr::select("group", "tissue")
tissue_plasma.trf.CountDat = tissue_plasma.trf.df.1 %>% dplyr::select(-c("group", "tissue")) %>% t()

colnames(tissue_plasma.trf.CountDat) == rownames(tissue_plasma.trf.ColDat)

# 血浆的trfNA做差异分析
coldata = tissue_plasma.trf.ColDat %>% rownames_to_column("id") %>% filter(tissue == "plasma") 
plasma_trf_res = NULL
x = 1
## DESEQ
for (i in c("HYP", "MCT", "HYSU")) {
  coldata.tmp <- coldata %>% filter(group %in% c("NOR", i))
  coldata.tmp$group <- factor(coldata.tmp$group, levels = c("NOR", i))
  countdata.tmp <- tissue_plasma.trf.CountDat[,coldata.tmp$id]
  dds <- DESeqDataSetFromMatrix(countData = countdata.tmp, colData = coldata.tmp, design = ~group) 
  dds$group <- relevel(dds$group, ref = "NOR") 
  res <- DESeq(dds) %>% results() %>% data.frame(id = rownames(.), ., tissue = "plasma")
  res %<>% na.omit() %>% mutate(change =
                                  case_when(log2FoldChange > 0.25 & pvalue < 0.05 ~ "Up_regulated", 
                                            log2FoldChange < -0.25 & pvalue < 0.05 ~ "Down_regulated",
                                            T ~ "Not_sig"
                                  )
  )
  plasma_trf_res[[x]] <- res
  names(plasma_trf_res)[x] <- paste(i, "Versus_NC", sep = "_")
  x = x + 1
}

plasma_trf_res$HYP_Versus_NC %>% 
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(pvalue), color = change)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("#699ed4","grey90","#ef8183")) + 
  scale_fill_manual(values = c("#699ed4","#ef8183")) + 
  xlim(-2, 2) +
  # scale_x_continuous(limits = c(-8, 8),breaks = seq(-8, 8, by = 4)) +
  # scale_y_continuous(limits = c(0, 150),breaks = seq(0, 150, by = 50)) +
  geom_vline(xintercept = .25, color = "#ef8183") +
  geom_vline(xintercept = -.25, color = "#699ed4") +
  geom_hline(yintercept = -log10(0.05), color = "grey50") +
  # theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = .75),
        panel.background = element_blank(), 
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold.italic")) -> vol.tDFplasma
vol.tDFplasma
setwd("/Users/jiahao_kuang/Dropbox/邝嘉浩/sncRNA in tissue/iscience/Revision/cor with plasma/Plasma Append Figures")
ggsave(plot = vol.tDFplasma, filename = "DEdDF_plasma_Volcano.pdf", units = "cm", width = 10, height = 8)

plasma_mir_res %>% write.csv(x = ., file = "plasma_miR.csv")
plasma_trf_res$HYP_Versus_NC %>% write.csv(x = ., file = "plasma_tDR.csv")

# 组织的trfNA做差异分析
coldata = tissue_plasma.trf.ColDat %>% rownames_to_column("id") %>% filter(tissue != "plasma") 
tissue_trf_res = NULL
# x = 1
## DESEQ
for (i in unique(coldata$tissue)) {
  coldata.tmp <- coldata %>% filter(tissue == i)
  # coldata.tmp <- coldata %>% filter(group %in% c("NOR", i))
  coldata.tmp$group <- factor(coldata.tmp$group, levels = c("NOR", "HYP"))
  countdata.tmp <- tissue_plasma.trf.CountDat[,coldata.tmp$id]
  dds <- DESeqDataSetFromMatrix(countData = countdata.tmp, colData = coldata.tmp, design = ~group) 
  dds$group <- relevel(dds$group, ref = "NOR") 
  res <- DESeq(dds) %>% results() %>% data.frame(id = rownames(.), ., tissue = i)
  res %<>% na.omit() %>% mutate(change =
                                  case_when(log2FoldChange > 0.25 & pvalue < 0.05 ~ "Up_regulated", 
                                            log2FoldChange < -0.25 & pvalue < 0.05 ~ "Down_regulated",
                                            T ~ "Not_sig"
                                  )
  )
  tissue_trf_res <- rbind(tissue_trf_res, res)
  # names(tissue_trf_res)[x] <- paste(i, "Versus_NC", sep = "_")
  # x = x + 1
}

# dcast(difres.trf, id~tissue, value.var = "log2FoldChange")
tissue_trf_res.wide = dcast(tissue_trf_res, id~tissue, value.var = "log2FoldChange")
plasma_trf_res.wide = plasma_trf_res$HYP_Versus_NC %>% dplyr::select("id", plasma = "log2FoldChange")

GGscatterPlot <- function(data, mapping, ..., 
                          method = "pearson") {
  
  #Get correlation coefficient
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor <- cor(x, y, method = method, use="pairwise.complete.obs")
  #Assemble data frame
  df <- data.frame(x = x, y = y)
  df <- na.omit(df)
  # PCA
  nonNull <- x!=0 & y!=0
  dfpc <- prcomp(~x+y, df[nonNull,])
  df$cols <- predict(dfpc, df)[,1]
  # Define the direction of color range based on PC1 orientation:
  dfsum <- x+y
  colDirection <- ifelse(dfsum[which.max(df$cols)] < 
                           dfsum[which.min(df$cols)],
                         1,
                         -1)
  #Get 2D density for alpha
  dens2D <- MASS::kde2d(df$x, df$y)
  df$density <- fields::interp.surface(dens2D ,df[,c("x", "y")])
  
  if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
  }
  #Prepare plot
  pp <- ggplot(df, aes(x=x, y=y, alpha = 1/density, color = cols)) +
    ggplot2::geom_point(shape=16, show.legend = FALSE) +
    ggplot2::scale_color_viridis_c(direction = colDirection) +
    ggplot2::scale_alpha(range = c(.05, .6)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = round(cor, digits = 3)),
      mapping = ggplot2::aes(x = xlabel,
                             y = ylabel,
                             label = lab),
      hjust = 0, vjust = 1,
      size = 3, fontface = "bold",
      inherit.aes = FALSE # do not inherit anything from the ...
    ) +
    theme_bw()
  return(pp)
}

cor.dat <- full_join(tissue_trf_res.wide, plasma_trf_res.wide, by = "id") %>%
  dplyr::select(-"id")

# for (i in 1:ncol(cor.dat)) {
#   cor.dat[,i][is.na(cor.dat[,i])] <- runif(sum(is.na(cor.dat[,i])), min = -0.000005, max = 0.000005)
# }
# df$column_name[is.na(df$column_name)] <- runif(sum(is.na(df$column_name)), min = -0.5, max = 0.5)

# is.na(cor.dat) %>% length()
# rnorm(numeric(length(is.na(cor.dat))),mean=0,sd=.5)
# cor.dat[is.na(cor.dat)] = rnorm(numeric(length(is.na(cor.dat))),mean=0,sd=.05)

# ggpairs(cor.dat.1)
cor.dat[is.na(cor.dat)] = 0

cor.dat$intestines = cor.dat$intestines + rnorm(67563,mean = 0, sd=.001)
cor.dat %>% 
GGally::ggpairs(.,
                # 1:4,
                lower = list(continuous = wrap(GGscatterPlot, method="pearson")),
                upper = list(continuous = wrap(ggally_cor, align_percent = 0.8))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text =  element_text(color='black'),
        strip.background = element_rect(fill = "#d63d2d"),
        strip.text = element_text(colour = "white", size = 10,
                                  face = "bold")) -> tmp.2

getwd()
ggsave(tmp.2, filename = "trf_log2FC.pdf", units = "cm", height = 19, width = 20)
ggsave(tmp.2, filename = "trf_log2FC.png", dpi = 1000, units = "cm", height = 19, width = 20)

# ceRNA 靶基因的富集分析
setwd("/Users/jiahao_kuang/Library/Mobile Documents/com~apple~CloudDocs/sncRNA landscape in HPH/ceRNA/")

library(clusterProfiler)
library(org.Rn.eg.db)

ce_target = read.table("cenods.1.txt", header = T)$mRNA %>% 
  clusterProfiler::bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)


# ce_target = read.table("cenods.3_reverseChange.txt", header = T)$to %>% 
#   clusterProfiler::bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)

GO = enrichGO(gene = ce_target$ENTREZID, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = .05, qvalueCutoff = .1) %>% data.frame()

GO$GeneRatio %<>% gsub(pattern = "/13", replacement = "") %>% as.numeric()
GO %<>% filter(GeneRatio >= 3)
GO = GO[(GO$qvalue %>% order(decreasing = F)),] %>% head(20)

# goid = c("GO:0051020", "GO:0001216", "GO:0001221", "GO:0031968", "GO:0030315",
#          "GO:0016323", "GO:0036293", "GO:0002065", "GO:0014909", "GO:0060759",
#          "GO:0006936","GO:0048660", "GO:0048771", "GO:0043534", "GO:0071320", 
#          "GO:0016055", "GO:0001959","GO:1904019", "GO:0070305", "GO:0090257" )

library(ggplot2)
GO.draw = GO
# GO.draw = GO.draw[order(GO.draw$Count,decreasing = T),]  

GO.draw.bp = GO.draw %>% filter(ONTOLOGY == "BP")
GO.draw.mf = GO.draw %>% filter(ONTOLOGY == "MF")
GO.draw.cc = GO.draw %>% filter(ONTOLOGY == "CC")

GO.draw = rbind(GO.draw.bp, GO.draw.mf, GO.draw.cc)
GO.draw$Description %<>% factor(levels = rev(GO.draw$Description))

write.csv(GO.draw, "/Users/jiahao_kuang/Dropbox/邝嘉浩/sncRNA in tissue/iscience/Revision/Text Revision/Supplementary Table/Table S10 - ceRNA hub DEmRNA GO Enrichment.csv")


GO.draw %>% 
  ggplot(data = ., mapping = aes(y = Description, x = -log10(qvalue), fill = -log10(qvalue))) +
  # ggplot(data = ., mapping = aes(y = reorder(Description, Count, decreasing = F), x = Count, fill = Count)) +
  geom_col(width = .5, color = "grey") +
  # coord_polar() +
  # scale_fill_gradientn(colours = c("#003366", "grey90", "#990033")) +
  # scale_fill_gradientn(colours = c("#3b374c", "#44598e", "#64a0c0", "#7ec4b7", "#deebcd"))+
  scale_fill_gradientn(colours = rev(c("#073f82", "#1b71b4", "#58a4cf", "#a2cbe3", "#f2f9fe"))) +
  # scale_fill_gradientn(colours = c("#80ab1c", "#405335", "#99b69b", "#92e4ce", "#72c8b7")) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_text(color = "black", face = "italic"),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1))
getwd()
ggsave(filename = "GO_ceRNA.pdf", units = "cm", height = 14, width = 17.5)


KEGG = enrichKEGG(gene = ce_target$ENTREZID, organism = "rno", pvalueCutoff = .05, qvalueCutoff = .1)

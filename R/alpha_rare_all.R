# 绘制alpha多样性稀疏曲线 Alpha rarefracation curve
# This function named 'alpha_rare_all',
# which draw curve by sample group,and curve with standard error by group with otutab and metadata, and return a ggplot2 object and plot data.

#' @title Plotting rarefracation curve for each sample or group
#' @description Input otutab and metadata, and manual set metadata column names.
#' ggplot2 show lineplot, and/or standard error
#' @param otutab OTU/ASV table;
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param groupID column name for groupID, such as "Group".
#' @param start resampling OTU/ASV table with the start number sequence count;
#' @param step number of intervals for resampling
#' @param method method for culculate alpha diversity,including "observed", "chao1", "diversity_shannon", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_fisher",  "diversity_coverage", "evenness_camargo", "evenness_pielou", "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp",  "dominance_dmn", "dominance_absolute","dominance_relative", "dominance_simpson", "dominance_core_abundance" ,  "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_noncore_abundance",  "rarity_rare_abundance"
#' @details
#' By default, returns a list with the curve and plot data
#' \itemize{
#' \item{most used indices: "observed", "chao1", "diversity_shannon", "evenness_simpson"}
#' \item{other used indices: "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_fisher",  "diversity_coverage", "evenness_camargo", "evenness_pielou", "evenness_evar", "evenness_bulla", "dominance_dbp",  "dominance_dmn", "dominance_absolute","dominance_relative", "dominance_simpson", "dominance_core_abundance" ,  "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_noncore_abundance",  "rarity_rare_abundance"}
#' @return ggplot2 object.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso alpha_boxplot alpha_rare
#' @examples
#' # rare plot with otutab and map
#' result = alpha_rare_all(otu = otutab, map = metadata, group = "Group", method = "chao1", start = 200, step = 200)
#' result[[1]]# output sample curve plot
#' result[[2]]# data table
#' result[[3]]# output group curve plot
#' result[[4]]# output group curve with CI plot
#' # rare plot with phyloseq Object
#' library(phyloseq)
#' ps = phyloseq(otu_table(otutab, taxa_are_rows=TRUE), sample_data(metadata))
#' result = alpha_rare_all(ps = ps, group = "Group", method = "chao1", start = 1000, step = 1000)
#' result[[4]]# output group curve with CI plot
#' @export


#---------纯R语言绘制稀释曲线#----

# # ----# 测试代码
# otu = read.delim("./data/otutab.txt",row.names = 1)
# map = read.delim("./data/metadata.tsv",row.names = 1)
# result1 = alpha_rare_all(otu = otu,map = map,group = "Group", start = 100,step = 100,method = "chao1")
# # 提取图形
# p =  result1[[1]]
# p
# # 提取数据
# data = result1[[2]]
# head(data)
# # 提取分组稀释曲线
# p =  result1[[3]]
# p
# # 提取分组稀释曲线
# p =  result1[[4]]
# p

# ps = readRDS("./data/ps_liu.rds")
# result = alpha_rare_all(ps = ps, start = 100,step = 1000,method = "observed")
# # 提取图形
# p =  result[[1]]
# p
# # 提取数据
# data = result[[2]]
# head(data)
# # 提取分组稀释曲线
# p =  result[[3]]
# p
# # 提取分组稀释曲线
# p =  result[[4]]
# p

# 稀释曲线函数
alpha_rare_all =function(otu = otutab,map = metadata,ps = NULL,method = "Richness",group = "Group", start = 100,step = 100){
  # 所需R包
  library(vegan)
  library(microbiome)
  library(tidyverse)
  # 设置默认参数测试函数
  # otu = otutab
  # map = metadata
  # ps = NULL
  # method = "Richness"
  # group = "Group"
  # start = 100
  # step = 100
  # 构建所需子函数
  #----抽平函数，输入对象为Phyloseq对象和抽平数值#----
  phyRare = function(ps = ps,N = 3000){
    library(phyloseq)
    #---- 提取OTU表函数 #----
    vegan_otu = function(physeq){
      OTU =  otu_table(physeq)
      if(taxa_are_rows(OTU)){
        OTU =  t(OTU)
      }
      return(as(OTU,"matrix"))
    }
    #---- 提取OTU表函数
    otb = as.data.frame(t(vegan_otu(ps)))
    otb1 = vegan::rrarefy(t(otb), N)
    ps = phyloseq(otu_table(as.matrix(otb1),taxa_are_rows = F),
                  sample_data(ps)
    )
    ps
    return(ps)
  }

  #----数据导入phyloseq#----
  if (is.null(ps) ) {
    head(otu)
    otu = as.matrix(otu)
    str(otu)
    map = map[group]
    colnames(map) = "Group"
    map$Group = as.factor(map$Group)
    ps = phyloseq(otu_table(otu, taxa_are_rows=TRUE),sample_data(map))
  }
  #----分析phyloseq对象#----
  if (!is.null(ps) ) {
    ps = ps
    map = as.data.frame(sample_data(ps))
    map = map[, group]
    colnames(map) = "Group"
    map$Group = as.factor(map$Group)
    sample_data(ps) = map
    map = NULL
  }

  #---- 全部指标#----
  all = c("observed" , "chao1"  , "diversity_inverse_simpson" , "diversity_gini_simpson",
          "diversity_shannon"   ,   "diversity_fisher"   ,  "diversity_coverage"     ,    "evenness_camargo",
          "evenness_pielou"    ,   "evenness_simpson"       ,    "evenness_evar" ,   "evenness_bulla",
          "dominance_dbp"      ,  "dominance_dmn"        ,      "dominance_absolute"   ,      "dominance_relative",
          "dominance_simpson"      ,    "dominance_core_abundance" ,  "dominance_gini"  ,           "rarity_log_modulo_skewness",
          "rarity_low_abundance"   ,    "rarity_noncore_abundance",  "rarity_rare_abundance")

  #--- 运行计算#----
  for (i in seq(start,max(sample_sums(ps)), by = step) ) {
    psRe = phyRare(ps = ps, N = i)

    vegan_otu =  function(physeq){
      OTU =  otu_table(physeq)
      if(taxa_are_rows(OTU)){
        OTU =  t(OTU)
      }
      return(as(OTU,"matrix"))
    }

    if (method == "Richness") {
      count = as.data.frame(t(vegan_otu(psRe)))
      # head(count)
      x = t(count) ##转置，行为样本，列为OTU
      est = estimateR(x)
      index = est[1, ]
    }

    if (method %in% c("ACE")) {
      ap_phy = estimate_richness(psRe, measures =method)
      # head(ap_phy)
      index = ap_phy$ACE
    }

    if (method %in% all) {
      alp_mic = alpha(psRe,index=method)
      # head(alp_mic)
      index = alp_mic[,1]
    }

    tab = data.frame(ID = names(sample_sums(psRe)))
    #得到多样性的列
    tab = cbind(tab,i,index)
    # head(tab)
    if (i == start) {
      result = tab
    }
    if (i != start) {
      result = rbind(result,tab)
    }
  }

   #----稀释结果为整齐的表，为了对应map分组吗？#----
  for (ii in 1:length(sample_sums(ps))) {
    result$i[result$i > sample_sums(ps)[ii][[1]]]
    df_filter= filter(result, ID ==names(sample_sums(ps)[ii]) &i > sample_sums(ps)[ii][[1]])
    result$index
    result$index[result$i>sample_sums(ps)[ii][[1]]]
    a = result$i>sample_sums(ps)[ii][[1]]
    a[a == FALSE] = "a"
    b = result$ID == names(sample_sums(ps)[ii])
    b[b == FALSE] = "b"
    result$index[a== b] = NA
  }
  #----直接给列，添加分组#----
  map = as.data.frame(sample_data(ps))
  result$Group = map$Group

  ## 绘制稀释曲线
  library(ggplot2)
  main_theme =theme(panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.title = element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =element_text(size = 7,face = "bold",colour = "black"),
    axis.title.x =element_text(size = 7,face = "bold",colour = "black"),
    axis.text = element_text(size = 7,face = "bold"),
    axis.text.x = element_text(colour = "black",size = 7),
    axis.text.y = element_text(colour = "black",size = 7),
    legend.text = element_text(size = 7,face = "bold")
  )
  p = ggplot(data= result,aes(x = i,y = index,group = ID,colour = Group)) +
    geom_smooth(span = 0.7, se = FALSE, method = "loess") +
    labs(x= "",y=method,title="") +theme_bw()+main_theme

  #---分组求均值和标准误+se#---
  data = result
  groups= group_by(data, Group,i)
  data2 = summarise(groups , mean(index), sd(index))
  # head(data2)
  colnames(data2) = c(colnames(data2)[1:2],"mean","sd")
  # 按组均值绘图
  p2 = ggplot(data= data2,aes(x = i,y = mean,colour = Group)) +
    geom_smooth(span = 0.7,se = FALSE, method = "loess") +
    labs(x= "",y=method,title="") +theme_bw()+main_theme
  # p2

  # 组均值+标准差
  p4 = ggplot(data=data2,aes(x = i,y = mean,colour = Group)) +
    geom_errorbar(data = data2,aes(ymin=mean - sd, ymax=mean + sd,colour = Group),alpha = 0.4, width=.1)+labs(x= "",y=method,title="") +theme_bw()+main_theme
  p4
  return(list(p,table = result,p2,p4))
}

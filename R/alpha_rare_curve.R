# 基于USEARCH结果绘制alpha多样性稀疏曲线 Alpha rarefaction curve
# This function named 'alpha_rare_curve',
# which draw rarefaction curve and standard error with alpha_rare_curve and metadata, and return a ggplot2 object

#' @title Plotting rarefaction curve for each group
#' @description Input alpha rare and metadata, and manual set metadata column names.
#' ggplot2 show line plot, and standard error
#' @param alpha_rare alpha rarefaction matrix, typical output of usearch -alpha_div_rare;
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param groupID column name for groupID.
#' @details
#' By default, returns grouped curve, sample parameter to be continued
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai.
#' A practical guide to amplicon and metagenomic analysis of microbiome data.
#' Protein Cell, 2020(41), 1-16, DOI: \url{https://doi.org/10.1007/s13238-020-00724-8}
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso alpha_boxplot
#' @examples
#' # usearch -alpha_div_rare result(alpha_rare), metadata and groupID
#' alpha_rare_curve(alpha_rare, metadata, "Group")
#' # alpha_rare, metadata, and groupID using site
#' alpha_rare_curve(alpha_rare, metadata, groupID="Site")
#' @export



alpha_rare_curve <- function(alpha_rare, metadata, groupID = "Group") {

  # 依赖关系检测与安装
  p_list = c("ggplot2","reshape2")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  # # 测试默认参数
  # dir="http://210.75.224.110/github/MicrobiomeStatPlot/Data/Science2019/"
  # alpha_rare = read.table(paste0(dir, "alpha/alpha_rare.txt"), row.names= 1, header=T, sep="\t",  comment.char="", stringsAsFactors = F)
  # metadata = read.table(paste0(dir, "metadata.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
  # groupID = "Group"

  # 交叉筛选
  idx = rownames(metadata) %in% colnames(alpha_rare)
  metadata = metadata[idx,]
  alpha_rare = alpha_rare[,rownames(metadata)]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

  # 求各组均值
  # 默认步长为1，折线不平滑，改为4减少锯齿
  df =alpha_rare[(1:25)*4,]
  # 转置df表格与实验设计合并，并去除第一列样品名
  mat_t = merge(sampFile, t(df), by="row.names")[,-1]

  # 采用group_by按组合并
  library(dplyr)
  df_g = as.data.frame(mat_t %>% group_by(group) %>% summarise_all(mean))
  rownames(df_g)=df_g$group
  df_g = df_g[,-1]
  mat_mean_final = as.data.frame(t(df_g))

  # 按第一列合并求均值，aggregate合并+do.call合并
  # mat_mean = aggregate(mat_t[,-1], by=mat_t[1], FUN=mean)
  # # 修正行名
  # mat_mean_final = do.call(rbind, mat_mean)[-1,]
  # geno = mat_mean$group
  # colnames(mat_mean_final) = geno

  df=as.data.frame(round(mat_mean_final))
  df$x = rownames(df)
  df_melt = melt(df, id.vars=c("x"))

  # 求各组标准误
  # 转置df表格与实验设计合并，并去除第一列样品名
  se = function(x) sd(x)/sqrt(length(x)) # function for Standard Error
  # mat_se = aggregate(mat_t[,-1], by=mat_t[1], FUN=se)
  # mat_se_final = do.call(rbind, mat_se)[-1,]
  # colnames(mat_se_final) = geno

  df_g = as.data.frame(mat_t %>% group_by(group) %>% summarise_all(se))
  rownames(df_g)=df_g$group
  df_g = df_g[,-1]
  mat_se_final = as.data.frame(t(df_g))

  df_se=as.data.frame(round(mat_se_final))
  df_se$x = rownames(df_se)
  df_se_melt = melt(df_se, id.vars=c("x"))

  # 添加标准误到均值中se列
  df_melt$se=df_se_melt$value

  # 添加levels顺序，否则按字母顺序
  df_melt$x = factor(df_melt$x, levels=c(1:100))

  p = ggplot(df_melt, aes(x = x, y = value, group = variable, color = variable )) +
    geom_line()+
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.5) +
    labs(x="Percentage (%)", y=paste("Richness"), color=groupID)+theme_classic()+
    scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10)+
    theme(text=element_text(family="sans", size=7))
  p
}

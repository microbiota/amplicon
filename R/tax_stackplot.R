# 样本或组的物种组成堆叠柱状图 Stackplot of taxonomy for samples and groups
#
# This is the function named 'tax_stackplot'
# which draw stack plot, and return a ggplot2 object
#
#' @title Plotting stackplot of taxonomy for groups or samples
#' @description Input taxonomy composition, and metadata (SampleID and groupID). Then select top N high abundance taxonomy and group other low abundance. When Select samples can draw sample composition by facet groups. If used group can show mean of each group. Finally, return a ggplot2 object.
#' @param tax_sum composition matrix, like OTU table and rowname is taxonomy, typical output of usearch -sintax_summary;
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param topN Top N taxonomy to show, default 8, alternative 4, 6, 10 ...;
#' @param groupID column name for groupID;
#' @param style group or sample, default group
#' @param sorted Legend sorted type, default abundance, alternative alphabet
#' @details
#' By default, returns top 8 taxonomy and group mean stackplot
#' The available style include the following:
#' \itemize{
#' \item{group: group mean stackplot}
#' \item{sample: each sample stackplot and facet by group}
#' }
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai.
#' A practical guide to amplicon and metagenomic analysis of microbiome data.
#' Protein Cell, 2020, DOI: \url{https://doi.org/10.1007/s13238-020-00724-8}
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso tax_stackplot
#' @examples
#' # Taxonomy table in phylum level, rownames is Phylum, colnames is SampleID
#' data(tax_phylum)
#' # metadata, include SampleID, Group and Site
#' data(metadata)
#' # tax_sum and metadata as input, default include top 8 taxonomy, groupID is Group, show group mean and sorted by abundance
#' tax_stackplot(tax_sum = tax_phylum, metadata)
#' # Set top10 taxonomy, group by "Site", group mean, and sort by abundance
#' tax_stackplot(tax_sum = tax_phylum, metadata, topN = 10, groupID = "Site", style = "group", sorted = "abundance")
#' # Set top 10 taxonomy, group by "Group", and sample composition sorted by alphabet
#' tax_stackplot(tax_sum = tax_phylum, metadata, topN = 10, groupID = "Group", style = "sample", sorted = "alphabet")
#' @export
tax_stackplot <- function(tax_sum, metadata, topN = 8, groupID = "Group", style = "group", sorted = "abundance") {

  # 依赖关系检测与安装
  p_list = c("ggplot2", "reshape2")
  for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  # 测试默认参数
  # library(amplicon)
  # tax_sum = tax_phylum
  # topN = 8
  # groupID = "Group"
  # style = "group"
  # sorted = "abundance"

  # 交叉筛选
  idx = rownames(metadata) %in% colnames(tax_sum)
  metadata = metadata[idx,,drop=F]
  tax_sum = tax_sum[, rownames(metadata)]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

  # 按丰度降序排序
  mean_sort = as.data.frame(tax_sum[(order(-rowSums(tax_sum))), ])
  # 筛选前N类，其它归为Other，可设置不同组数
  other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
  mean_sort = mean_sort[1:(topN - 1), ]
  mean_sort = rbind(mean_sort,other)
  rownames(mean_sort)[topN] = c("Other")
  # 保存变量备份，并输出至文件
  merge_tax=mean_sort

  if (style == "sample"){
  # 添加分类学并设置排序方式，默认字母，abundancer按丰度
  mean_sort$Taxonomy = rownames(mean_sort)
  data_all = as.data.frame(melt(mean_sort, id.vars=c("Taxonomy")))
  if (sorted == "abundance"){
  data_all$Taxonomy  = factor(data_all$Taxonomy, levels=rownames(mean_sort))
  }
  # set group facet
  data_all = merge(data_all, sampFile, by.x="variable", by.y = "row.names")

  # 按group分面，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  # 关闭x轴刻度和标签

  p = ggplot(data_all, aes(x=variable, y = value, fill = Taxonomy )) +
    geom_bar(stat = "identity",position="fill", width=1)+
    scale_y_continuous(labels = scales::percent) +
    facet_grid( ~ group, scales = "free_x", switch = "x") +
    theme(strip.background = element_blank())+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
    xlab("Groups")+ylab("Percentage (%)")+
    theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))+
    theme(text=element_text(family="sans", size=7))
  p
  }else{
    # 按组合并求均值

    # 转置样品名添加组名，并去除多余的两个样品列
    mat_t = t(merge_tax)
    mat_t2 = merge(sampFile, mat_t, by="row.names")
    mat_t2 = mat_t2[,c(-1)]

    # 按组求均值，转置，再添加列名
    mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
    mat_mean_final = do.call(rbind, mat_mean)[-1,]
    geno = mat_mean$group
    colnames(mat_mean_final) = geno
    mean_sort=as.data.frame(mat_mean_final)

    # 添加分类学并设置排序方式，默认字母，abundancer按丰度
    mean_sort$Taxonomy = rownames(mean_sort)
    data_all = as.data.frame(melt(mean_sort, id.vars=c("Taxonomy")))
    if (sorted == "abundance"){
      data_all$Taxonomy  = factor(data_all$Taxonomy, levels=rownames(mean_sort))
    }
    data_all$value = as.numeric(data_all$value)
    p = ggplot(data_all, aes(x=variable, y = value, fill = Taxonomy )) +
      geom_bar(stat = "identity",position="fill", width=0.7)+
      scale_y_continuous(labels = scales::percent) +
      xlab("Groups")+ylab("Percentage (%)")+ theme_classic()+
      theme(text=element_text(family="sans", size=7))
    p
  }
}

# 绘制beta多样性组间统计
#
# This is the function named 'beta_pcoa_stat'
# Similar with beta_pcoa, and save a table with P-value
#
#' @title Calculate beta diversity p-value by Adonis
#' @description Input distance matrix and metadata, and manual set metadata column names. Save a table with P-value by Adonis.
#' @param dis_mat distance matrix, typical output of usearch -beta_div,
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param groupID column name for groupID.
#' @param pairwise group pairwise compare, only groups less than 5, T or F.
#' @param pairwise_list a file include pairwise list.
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Yong-Xin Liu, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, Liang Chen,
#' Xubo Qian, Jiao Xi, Hongye Lu, Huiluo Cao, Xiaoya Ma, Bian Bian, Pengfan Zhang,
#' Jiqiu Wu, Ren-You Gan, Baolei Jia, Linyang Sun, Zhicheng Ju, Yunyun Gao, Tao Wen,
#' Tong Chen. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research.
#' iMeta 2: e83. \url{https://doi.org/10.1002/imt2.83}
#'
#' @seealso beta_pcoa
#' @examples
#' # Set essential 3 parameters: distance matrix, metadata and groupID
#' beta_pcoa_stat(beta_bray_curtis, metadata, "Group", "beta_pcoa_stat.txt")
#' # Set 5 parameters: dis_mat, metadata, and groupID, pairwise FLASE, can set pairwise_list in file
#' beta_pcoa_stat(dis_mat, metadata, groupID="Group", pairwise=F, pairwise_list="vignettes/compare.txt")
#' @export

beta_pcoa_stat <- function(dis_mat, metadata, groupID = "Group", result = "beta_pcoa_stat.txt", pairwise = T, pairwise_list = "vignettes/compare.txt") {

  # # 测试默认参数
  # library(amplicon)
  # dis_mat = beta_unifrac
  # metadata = metadata
  # groupID = "Group"
  # result = "beta_pcoa_stat.txt"
  # pairwise = T
  # pairwise_list = "compare.txt"

  # 依赖关系检测与安装
  p_list = c("vegan")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  metadata$group = metadata[[groupID]]
  write.table(date(), file=result, append=T, sep="\t", quote=F, row.names=F, col.names=F)

  # 交叉筛选
  idx=rownames(metadata) %in% rownames(dis_mat)
  metadata=metadata[idx,,drop=F]
  dis_mat=dis_mat[rownames(metadata), rownames(metadata)]

  # Compare each group beta by vegan adonis in bray_curtis
  da_adonis = function(sampleV){
    # sampleV = tmp_compare
    sampleA = as.matrix(sampleV$sampA)
    sampleB = as.matrix(sampleV$sampB)
    design2 = subset(metadata, group %in% c(sampleA,sampleB))
    if (length(unique(design2$group))>1) {
      sub_dis_table = dis_table[rownames(design2),rownames(design2)]
      sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
      adonis_table = adonis2(sub_dis_table~group, data=design2, permutations = 10000)
      adonis_pvalue = adonis_table$`Pr(>F)`[1]
      # print(paste("In ",opts$type," pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
      adonis_pvalue = paste(sampleA, sampleB, adonis_pvalue, sep="\t")
      write.table(adonis_pvalue, file=result, append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    }
  }

  # loop for each group pair
  dis_table = as.matrix(dis_mat)
  if (pairwise) {
    compare_data = as.vector(unique(metadata$group))
    len_compare_data = length(compare_data)
    for(i in 1:(len_compare_data-1)) {
      for(j in (i+1):len_compare_data) {
        # i=1
        # j=2
        tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
        print(tmp_compare)
        da_adonis(tmp_compare)
      }
    }
  }else {
    compare_data = read.table(pairwise_list, sep="\t", check.names=F, quote='', comment.char="")
    colnames(compare_data) = c("sampA", "sampB")
    for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
  }
  write.table("", file=result, append=T, sep="\t", quote=F, row.names=F, col.names=F)

}

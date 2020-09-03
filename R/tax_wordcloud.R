# Wordcloud plot
#
# The function named 'tax_wordcloud'
# which wordcloud is another way to visualize compositional microbial community data
#
# You can learn more about package at:
#
#   https://github.com/microbiota/amplicon

#' @title Wordcloud plot relative abundance of taxonomy
#' @description Input otutab, metadata and taxonomy or phyloseq object; plot wordcloud to visualize compositional microbial community data
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param tax taxonomy file;
#' @param group group ID;
#' @param j taxonomy annotation classification level, typical including "Phylum","Class","Order","Family","Genus"
#' @param facet FALSE/F, if TURE/T, facet with Phylum in one of the "Class","Order","Family","Genus" level
#' @param abundance TRUE/T, show abundance in label, turn off in FLASE/F;
#' @param rand set reproducible layout, default 1, can set any integer;
#' @details
#' By default, input phyloseq object include metadata, otutab and metadata
#' The available classification level include the following:
#' \itemize{
#' \item{most used classification level was Phylum}
#' \item{other classification level include: Class, Order, Family, Genus }
#' }
#' @return list object including plot and data table
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai.
#' A practical guide to amplicon and metagenomic analysis of microbiome data.
#' Protein Cell, 2020(41), 1-16, DOI: \url{https://doi.org/10.1007/s13238-020-00724-8}
#'
#' @seealso tax_stackplot tax_circlize tax_maptree
#' @examples
#' # Input feature table, metadata, taxonomy, group column name, and taxonomic level
#' # Output each group Phylum wordcloud with relative abundance in bracket
#' tax_wordcloud(otu = otutab, map = metadata, tax = taxonomy, group = "Group", j = "Phylum")
#' # Input feature table, metadata, taxonomy, group column name, and taxonomic level
#' # in Class, facet by Phylum, not show abundance, layout in seed 2
#' tax_wordcloud(otu = otutab, map = metadata, tax = taxonomy, group = "Group", j = "Class", facet = T, abundance = F, rand = 2)
#' @export

tax_wordcloud = function(otu = otutab, map = metadata, tax = taxonomy, ps = NULL, j = "Phylum", group = "Group", facet = F, abundance = T, rand = 1){

  #----设置参数默认值测试函数#----
  # otu = otutab
  # map = metadata
  # tax = taxonomy
  # ps = NULL
  # j = "Phylum"
  # group = "Group"
  # facet = F
  # abundance = T
  # rand = 1

  #----加载R包#----
  suppressWarnings(suppressMessages(library(tidyverse)))
  suppressWarnings(suppressMessages(library(phyloseq)))
  suppressWarnings(suppressMessages(library(ggwordcloud)))

  #----自定义函数#----

  ## 从phyloseq对象中提取OTU表格
  vegan_otu =  function(physeq){
    OTU =  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU =  t(OTU)
    }
    return(as(OTU,"matrix"))
  }

  ### 从phyloseq对象中提取物种注释信息
  vegan_tax =  function(physeq){
    tax =  tax_table(physeq)
    return(as(tax,"matrix"))
  }

  #----数据导入PhyloSeq对象#----
  if (is.null(otu)&is.null(map)&is.null(tax)) {
    ps = ps
  }else{
    #导入otu表格
    otu = as.matrix(otu)
    # str(otu)
    #导入注释文件
    tax = as.matrix(tax)
    # taxa_names(tax)
    #导入分组文件
    map$Group =  as.factor(map[, group])
    # map$Group
    ps = phyloseq(otu_table(otu, taxa_are_rows=TRUE),
                   sample_data(map),
                   tax_table(tax)
    )
  }
  # ps


  #----按门/纲...水平分类汇总#----
  ps_sub = ps %>% tax_glom(taxrank = j) %>% transform_sample_counts(function(x) {x/sum(x)})
  # 为什么2631个ASV只剩下18个ASV?因为合并成门，但仍用ASV作行名
  # ps_sub
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  # head(otu_table)

  #----不同分组求均值#----
  # 可以简化为几行，用group_by
  count = t(otu_table)
  count2 = as.data.frame(count)
  # head(count)
  #提取分组文件
  sub_design = as.data.frame(sample_data(ps))
  #构造分组
  aa = sub_design[,"Group"]
  colnames(aa) = "Vengroup"
  iris.split = split(count2,as.factor(aa$Vengroup))
  #数据分组计算平均值
  iris.apply = lapply(iris.split,function(x)colMeans(x[]))
  # 组合结果
  iris.combine = do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  ven2 = as.data.frame(ven2)
  ven2$ID  =row.names(ven2)
  # 提取物种注释文件：
  tax_table = as.data.frame(vegan_tax(ps_sub))
  ven3 = cbind(tax_table,ven2)


  #----宽表格转换ggplot输入长表格#----
  library(reshape2)
  plotdata = melt(ven3,id.vars = c("ID",colnames(tax_table)),
                 variable.name='Group',
                 value.name='mean')
  #----物种名添加丰度#----
  if (abundance == T) {
    plotdata[,j] = paste(plotdata[,j],"(",round(plotdata[,"mean"],3),")",sep = "")
  }
  #----不同分组求均值#----
  # 分组展示词云
  set.seed(rand)
  p = ggplot(plotdata, aes_string(label = j, size = "mean",color = "Phylum")) +
    geom_text_wordcloud() +
    scale_size_area(max_size = 15) +
    theme_minimal() +
    facet_wrap(~Group)
  # p

  # 分组+门分面
  if(facet == TRUE){
    p = ggplot(plotdata, aes_string(label = j, size = "mean",color = "Phylum")) +
      geom_text_wordcloud() +
      scale_size_area(max_size = 15) +
      theme_minimal() +facet_grid(Phylum~Group)
    p
  }

  # 返回词云和出图数据
  return(list(p,plotdata))
}

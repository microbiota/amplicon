# 样本或组的物种组成弦图 Circlize of taxonomy for samples and gorups
#
# This is the function named 'data_to_maptree'
# which draw circle, and reture a circlize object
#
#' @title Plotting circlize of taxonomy for groups or samples
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
#' \item{group: group mean circlize}
#' \item{sample: each sample circlize}
#' }
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso data_to_maptree
#' @examples
#' # example data: feature table, rownames is OTU/taxonomy, colnames is SampleID
#' data(tax_phylum)
#' # example data: metadata or design, include SampleID, genotype and site
#' data(metadata)
#' # Set 4 parameters: set top 5 taxonomy, group by "genotype"
#' tax_circlize(tax_sum = tax_phylum, metadata, topN = 5, groupID = "genotype")
#' @export

##----构造高级节点注释文件----
tax_maptree = function(x){
  ps_sub = mapdata[[4]]
  vertices_t = mapdata[[3]]
  deg = mapdata[[2]]
  head(vertices_t)

  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)

    return(as(tax,"matrix"))
  }

  tax_table = as.data.frame(vegan_tax(ps_sub))
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  #计算全部均值
  otu_table$mean = rowMeans(otu_table)
  #组合结果
  inde = merge(otu_table,tax_table,by = "row.names", all = TRUE)
  head(inde)
  row.names(inde) = inde$Row.names
  inde$Row.names = NULL

  #添加高级注释信息进入函数
  row.names(inde) = paste("ASV",row.names(inde),sep = "--")
  ##将全部信息合并到节点属性中
  data_plot = merge(vertices_t,inde,by = "row.names",all = TRUE)
  row.names(data_plot) = data_plot$Row.names
  data_plot$Row.names = NULL

  #指定大小映射列将NA值做转换为0
  asa =data_plot$mean
  data_plot$mean[is.na(asa)] = 0
  #选择是否平方根标准化丰度
  # data_plot$mean[!is.na(asa)] = sqrt(data_plot$mean[!is.na(asa)])


  mygraph <- graph_from_data_frame( deg, vertices= data_plot )
  #-----------------------------------设置颜色映射参数-------------------------
  data = create_layout(mygraph, layout = 'circlepack',weight = mean, sort.by = NULL, direction = "out")
  #设置颜色
  data$Phylum = as.character(data$Phylum)
  data$Phylum[is.na(data$Phylum)] = "AA"
  data$Phylum = factor(data$Phylum,levels =unique(data$Phylum) )
  colbar <-length(unique(data$Phylum))
  fil_colors = colorRampPalette(c( "white","#CBD588", "#599861", "orange","#DA5724", "#508578", "#CD9BCD",
                                   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
                                   "#8569D5", "#5E738F","#D1A33D", "#8A7C64","black"))(colbar)
  names(fil_colors ) = unique(data$Phylum)
  fil_colors[1] = "white"
  p = ggraph(mygraph, layout = 'circlepack',weight = mean, sort.by = NULL, direction = "out") +
    geom_node_circle(aes(fill = as.factor(Phylum),color = as.factor(depth) ) ) +
    scale_fill_manual(values= fil_colors ) +
    scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5"="black", "6"="black", "7"="black") ) +
    geom_node_text( aes(label=new_label), size=6,repel = TRUE) +
    # geom_node_label( aes(label=new_label), size=3,repel = TRUE) +
    theme_void() +
    theme( legend.position="FALSE",plot.margin = unit(rep(0,4), "cm"))#
  # p
  return(list(p,mygraph,fil_colors))
}


# tax_maptree <- function(otutab, taxonomy, dot = 200) {
#   # otutab = as.matrix(otutab)
#   # taxonomy = as.matrix(taxonomy)
#   # dot = 200
#   mapdata = data_to_maptree (otutab,taxonomy,dot)
#   mapadd = maptree_add1_plot(mapdata)
#   p = mapadd[[1]]
# }

# 绘制beta多样性PCoA图+置信椭圆 Beta PCoA + stat ellipse
#
# This is the function named 'beta_pcoa'
# which draw PCoA scatter plot with stat ellipse, and return a ggplot2 object
#
#' @title Plotting beta diversity scatter plot
#' @description Input distance matrix and metadata, and manual set metadata column names.
#' Visualize PCoA with color and stat ellipse by ggplot2.
#' @param dis_mat distance matrix, typical output of usearch -beta_div,
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param groupID column name for groupID.
#' @param ellipse stat ellipse, T or F.
#' @param label sample name showing, T or F.
#' @param PCo principle coordinate used, default 12, alternative 13, or 23.
#' @details
#' By default, returns beta PCoA coordinate
#' The available diversity indices include the following:
#' \itemize{
#' \item{most used indices: bray_curtis, unifrac}
#' \item{other used indices: unifrac_binary, jaccard, euclidean, manhatten}
#' }
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
#' @seealso beta_cpcoa
#' @examples
#' # Set essential 3 parameters: distance matrix, metadata and groupID
#' beta_pcoa(beta_bray_curtis, metadata, "Group")
#' # Set full 6 parameters: distance matrix, metadata, and groupID as using "site",
#' beta_pcoa(dis_mat=beta_unifrac, metadata=metadata, groupID="Site", ellipse=F, label=T, PCo=13)
#' @export

beta_pcoa <- function(dis_mat, metadata, groupID="Group", ellipse=T, label=F, PCo=12) {
  # 依赖关系检测与安装
  p_list=c("ggplot2", "vegan", "ggrepel")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    suppressWarnings(suppressMessages(library(p)))}

  # 测试默认参数
  # dis_mat=beta_unifrac
  # metadata=metadata
  # groupID="Group"
  # ellipse=T
  # label=F
  # PCo=12

  # 交叉筛选
  idx=rownames(metadata) %in% rownames(dis_mat)
  metadata=metadata[idx,]
  dis_mat=dis_mat[rownames(metadata), rownames(metadata)]

  # 提取样品组信息,默认为group可指定
  sampFile=as.data.frame(metadata[, groupID],row.names=row.names(metadata))
  # colnames(sampFile)[1]="group"

  # PCoA
  pcoa=cmdscale(dis_mat, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points=as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  eig=pcoa$eig
  points=cbind(points, sampFile[rownames(points),])
  colnames(points)=c("x", "y", "z","group")

  # 按1、2轴绘图
  if (PCo == 12){
    p=ggplot(points, aes(x=x, y=y, color=group))  +
      labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按1、3轴绘图
  if (PCo == 13){
    p=ggplot(points, aes(x=x, y=z, color=group))  +
      labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按2、3轴绘图
  if (PCo == 23){
    p=ggplot(points, aes(x=y, y=z, color=group))  +
      labs(x=paste("PCo 2 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  p=p + geom_point(alpha=.7, size=2) + theme_classic() + theme(text=element_text(family="sans", size=7))
  # 是否添加置信椭圆
  if (ellipse == T){
    p=p + stat_ellipse(level=0.68)
  }
  # 是否显示样本标签
  if (label == T){
    p=p + geom_text_repel(label=paste(rownames(points)), colour="black", size=3)
  }
  p
}

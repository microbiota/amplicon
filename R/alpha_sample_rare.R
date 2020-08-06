# 基于特征表绘制样本量的alpha多样性稀释箱线图 Alpha sample rarefaction boxplot
# This function named 'alpha_sample_rare',
# which draw sample rarefaction boxplot, and return a ggplot2 object

#' @title Plotting sample rarefaction boxplot
#' @description Input feature (OTU/ASV) table
#' @param otutab OTU table
#' @param length Number of boxplot, recommend < 30 & sample number
#' @param count_cutoff cutoff of reads count , default is 1
#' @param rep replicate sample times, is the dot number in boxplot, default 30, recommend < 30 & sample number
#' @details
#' By default, returns ggplot2
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
#' @seealso alpha_boxplot alpha_rare_curve alpha_rare_all
#' @examples
#' # feature table, all parameters in default. Length is sample size times(boxplot number), rep is sampling times (dots in boxplot), count_cutoff is minimum reads count as detectable feature
#' alpha_sample_rare(otutab, length=18, rep=30, count_cutoff=1)
#' # Change the number of boxplot
#' alpha_sample_rare(otutab, length=9)
#' # Cutoff apparently affect trends
#' alpha_sample_rare(otutab, count_cutoff=9)
#' @export

alpha_sample_rare = function(otutab, length=18, rep=30, count_cutoff=1){

  # 设置初始值
  count=otutab
  # length=18
  # rep=30
  # count_cutoff=1

  # 转换特征表为二元值，计算丰富度
  # Set otu table to binary for easy calculate richness
  count[count < count_cutoff] = 0
  count[count >= count_cutoff] = 1
  # 样本量Sample number
  n = dim(count)[2]
  # 设置X轴各取样数量，即各步长样本量
  x = unique(as.integer(seq(1, n, length.out = length)))

  # 定义结果数据框
  result = data.frame(sample = c(0), richness = c(0))

  # 按步长统计多样性，每个步长也抽rep次
  for (i in x){
    m = choose(n, i)
    if (m > rep){m = rep}
    # loop list and calculate richness
    for(j in 1:m) {
      idx = sample(n, i)
      temp = count[, idx, drop = F]
      # subset non-zero row
      temp1=temp[rowSums(temp)>0, , drop=F]
      # row number is accumulative OTUs
      result = rbind(result, c(i, dim(temp1)[1]))
    }
  }
  # 移除空行 remove start 0,0
  result = result[-1,]
  # 设置X轴顺序 factor draw as each box
  result$sample = factor(result$sample, levels = unique(result$sample))
  # 设置绘图主题
  main_theme = theme(panel.background=element_blank(),
    panel.grid=element_blank(),
    axis.line.x=element_line(size=.5, colour="black"),
    axis.line.y=element_line(size=.5, colour="black"),
    axis.ticks=element_line(color="black"),
    axis.text=element_text(color="black", size=7),
    legend.position="right",
    legend.background=element_blank(),
    legend.key=element_blank(),
    legend.text= element_text(size=7),
    text=element_text(family="sans", size=7))

  # 绘图
  p = ggplot(result,aes(x=sample,y=richness, fill=sample))+
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.5, size = 0.1)+
    main_theme+ xlab("Sample size")+ylab("Feature richness") +
    theme(legend.position = "NA",axis.text.x = element_text(angle = 45,vjust=1, hjust=1))
  p
}

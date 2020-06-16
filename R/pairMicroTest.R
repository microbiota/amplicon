# Beta diversity statistics in pair
#
# The function named 'pairMicroTest'
# which call adonis, anosim or MRPP to test beta-diversity. It usually called by BetaDiv.
#
# You can learn more about package at:
#
#   https://github.com/microbiota/amplicon

#' @title Beta diversity statistics by adonis/anosim/MRPP in pair
#' @description Input phyloseq object, test method and distance type
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param ps alternative input;
#' @param tree tree/nwk file;
#' @param dist distance type, including "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski"  "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial"  "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co";
#' @param group group ID;
#' @param method DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA;
#' @param pvalue.cutoff Pvalue threshold, default in 0.05;
#' @param Micromet statistics default by adonis, alternative anosim or MRPP;

#' @details
#' By default, input phyloseq object include metadata and otutab
#' The available diversity indices include the following:
#' \itemize{
#' \item{most used indices: bray unifrac wunifrac}
#' \item{other used indices: dpcoa jsd manhattan euclidean canberra kulczynski jaccard gower altGower morisita horn mountford raup binomial chao cao w -1 c wb r I e t me j sor m -2 co}
#' }
#' @return stat table
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
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
#' @seealso BetaDiv beta_pcoa beta_cpcoa
#' @examples
#' # Input phyloseq format input, and options group, method and distance
#' # 生成phyloseq对象
#' ps=phyloseq(otu_table(otutab_rare, taxa_are_rows=TRUE), sample_data(metadata))
#' pairMicroTest (ps=ps, Micromet="anosim", dist="bray")
#' @export

pairMicroTest=function(ps=ps, Micromet="anosim", dist="bray"){

  if (!requireNamespace("vegan", quietly=TRUE))
    install.packages("vegan")
  library(vegan)
  # 安装Bioconductor的R包phyloseq
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  suppressWarnings(suppressMessages(library(BiocManager)))
  if (!requireNamespace("phyloseq", quietly=TRUE))
    BiocManager::install("phyloseq")
  library(phyloseq)

  # 生成phyloseq对象
  # ps=phyloseq(otu_table(otutab_rare, taxa_are_rows=TRUE),
  #               sample_data(metadata))
  ps1_rela =transform_sample_counts(ps, function(x) x / sum(x) );ps

  #-准备矩阵和分组文件
  map=as.data.frame(sample_data(ps1_rela))# ,stringsAsFactors=T
  # 转换为后，levels消失了，改用unique添加levels
  map$Group=factor(map$Group, levels=unique(map$Group))
  (aa=levels(map$Group))
  (aaa=combn(aa,2))
  dim(aaa)[2]

  # 构建三个空列
  ID=rep("a",dim(aaa)[2])
  R=rep("a",dim(aaa)[2])
  P=rep("a",dim(aaa)[2])
  # i=1
  for (i in 1:dim(aaa)[2]) {
    # print(i)
    Desep_group=aaa[,i]
    map=as.data.frame(sample_data(ps1_rela))
    # head(map)
    map$ID=row.names(map)
    # maps=dplyr::filter(map, Group %in% Desep_group)
    # 取子集 #----
    maps=subset(map, Group %in% Desep_group)
    row.names(maps)=maps$ID
    ps_sub=ps1_rela
    sample_data(ps_sub)=maps
    ps_sub=phyloseq::filter_taxa(ps_sub, function(x) sum(x ) > 0 , TRUE);ps_sub
    map=as.data.frame(sample_data(ps_sub))
    unif <- phyloseq::distance(ps_sub, method=dist)

    if (Micromet == "MRPP") {
      mrpp=vegan::mrpp(unif, map$Group)
      as1=round(mrpp$delta,3)
      R2 <- paste("MRPP.delta ",as1, sep="")
      # R[i]=R2
      R2
      p_v=paste("p: ",round(mrpp$Pvalue,3), sep="")
      # p_v
      # P[i]=p_v
    }

    if (Micromet == "anosim") {
      dat.ano=anosim(unif, map$Group)
      a=round(dat.ano$statistic,3)
      R2 <- paste("ANOSIM.r ",a, sep="")
      R[i]=R2
      p_v=paste("p: ",round(dat.ano$signif,3), sep="")
      # P[i]=p_v
    }

    gg =map$Group
    if (Micromet == "adonis") {
      ado= adonis(unif~gg,permutations=999)
      a=round(as.data.frame(ado$aov.tab[5])[1,1],3)
      R2 <- paste("adonis:R ",a, sep="")
      R[i]=R2
      b=as.data.frame(ado$aov.tab[6])[1,1]
      p_v=paste("p: ",b, sep="")
    }
    ID[i]=paste(Desep_group[1],Desep_group[2],sep="_VS_")
    P[i]=p_v
    R[i]=R2
  }
  # P
  # R
  result=data.frame(ID=ID,stat=R,p=P)
  result

  return(result)
}

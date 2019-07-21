# 生成LEfSe输入分类级汇总表 format to stamp
#
# This is the first function named 'format2lefse'

#' @title Format otutab and taxonomy into each level for STAMP
#' @description Input otutab, taxonomy and metadata, and manual set abundance threshold, metadata column names.
#' dplyr merge taxonomy.
#' @param otutab row or normalized OTU table
#' @param taxonomy Taxonomy include seven taxonomy level in tsv format
#' @param metadata matrix or dataframe, including sampleID and groupID
#' @param thre threshold of OTU abundance, especially for lefse create properly cladegram
#' @param groupID column name for groupID.
#' @details
#' By default, written 1 file
#' @return nothing, all table 0-8 saving in same directory.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#' Zhang, J., Zhang, N., Liu, YX. et al.
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage
#' Sci. China Life Sci. (2018) 61: 613.  DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#' @seealso alpha_rare
#' @examples
#' # Set six parameters: three input files, and three adjust parameters
#' format2lefse(otutab, taxonomy, metadata, thre = 0.01, groupID = "Group")
#' @export

format2lefse <- function(otutab, taxonomy, metadata, thre = 0.01, groupID = "Group", output = "LEfSe.txt") {

  # 依赖关系检测与安装
  p_list = c("dplyr")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  # 测试默认参数
  # thre = 0.01
  # groupID = "Group"

  # 标准化，并筛选高丰度菌均值最小百万分之一0.0001%
  norm = t(t(otutab)/colSums(otutab,na=T))*100
  dim(norm)
  colSums(norm)
  idx = rowMeans(norm) > thre
  HA = norm[idx,]
  dim(HA)
  colSums(HA)
  # 保在筛选后的OTU表
  # write.table(paste("OTU\t",  sep=""), file=paste(prefix, "_8OTU", thre, ".txt", sep = ""), append = F, sep="\t", quote=F, row.names=F, col.names=F, eol = "")
  # write.table(HA, file=paste(prefix, "_8OTU", thre, ".txt", sep = ""), append = T, sep="\t", quote=F, row.names=T, col.names=T)
  # 数据筛选并排序，要求每个OTU必须的注释，可以为空
  tax = tax[rownames(HA),]

  # 转换为等级|连接格式
  tax$Phylum=paste(tax$Kingdom,tax$Phylum,sep = "|")
  tax$Class=paste(tax$Phylum,tax$Class,sep = "|")
  tax$Order=paste(tax$Class,tax$Order,sep = "|")
  tax$Family=paste(tax$Order,tax$Family,sep = "|")
  tax$Genus=paste(tax$Family,tax$Genus,sep = "|")
  tax$Species=paste(tax$Genus,tax$Species,sep = "|")
  # head(tax)

  # 按Kingdom合并
  grp <- tax[rownames(tax), "Kingdom", drop=F]
  merge=cbind(HA, grp)
  HA_Kingdom = merge %>% group_by(Kingdom) %>% summarise_all(sum)
  # colnames(HA_Kingdom)[1]="Kingdom"
  # write.table(HA_Kingdom, file=paste(prefix, "_1Kingdom.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Kingdom)[1]="Class"

  # 按Phylum合并
  grp <- tax[rownames(tax), "Phylum", drop=F]
  merge=cbind(HA, grp)
  HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
  # colnames(HA_Phylum)[1]="Phylum"
  # write.table(HA_Phylum, file=paste(prefix, "_2Phylum.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Phylum)[1]="Class"

  # 按Class合并
  grp <- tax[rownames(tax), "Class", drop=F]
  merge=cbind(HA, grp)
  HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
  # colnames(HA_Class)[1]="Class"
  # write.table(HA_Class, file=paste(prefix, "_3Class.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Class)[1]="Class"

  # 按Order合并
  grp <- tax[rownames(tax), "Order", drop=F]
  merge=cbind(HA, grp)
  HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
  # colnames(HA_Order)[1]="Order"
  # write.table(HA_Order, file=paste(prefix, "_4Order.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Order)[1]="Class"

  # 按Family合并
  grp <- tax[rownames(tax), "Family", drop=F]
  merge=cbind(HA, grp)
  HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
  # colnames(HA_Family)[1]="Family"
  # write.table(HA_Family, file=paste(prefix, "_5Family.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Family)[1]="Class"

  # 按Genus合并
  grp <- tax[rownames(tax), "Genus", drop=F]
  merge=cbind(HA, grp)
  HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
  # colnames(HA_Genus)[1]="Genus"
  # write.table(HA_Genus, file=paste(prefix, "_6Genus.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Genus)[1]="Class"

  # 按Species合并
  # grp <- tax[rownames(tax), "Species", drop=F]
  # merge=cbind(HA, grp)
  # HA_Species = merge %>% group_by(Species) %>% summarise_all(sum)
  # colnames(HA_Species)[1]="Species"
  # write.table(HA_Species, file=paste(prefix, "_7Species.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  # colnames(HA_Species)[1]="Class"

  # 合并6个分类级
  all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus) # , HA_Species

  # 将选定的分组列统一命名为group
  metadata$group = metadata[, groupID]

  # 修改样品名为组名
  # 删除结尾的数字，兼容性不好
  # colnames(all) = gsub("\\d+$","",colnames(all),perl=TRUE)
  # 用实验设计中的group列替换样品名
  colnames(all)[2:dim(all)[2]] = as.character(metadata[colnames(all)[2:dim(all)[2]],]$group)

  # 保存结果为lefse
  write.table(all, file=paste(output, sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
}

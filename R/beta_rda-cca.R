
# RDA or CCA calculate
#
# The function named 'RDA_CCA'
# which do RDA or CCA calculate and plot, at the same time, this script would also culculate the  impacts of environmental indicators to Microbial community
# 
# You can learn more about package at:
#
#   https://github.com/microbiota/amplicon

#' @title RDA or CCA calculate
#' @description Input otutab, metadata and environmental indicators or phyloseq object; support Automatic select which method : RDA or CCA were used for analyse; output ggplot2 figure, data and statistical test result.
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param env environmental indicators table;
#' @param group group ID;
#' @details
#' By default, input phyloseq object include metadata, otutab and env tabe
#' @return list object including plot, stat table
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @examples
#'
#' RDA_CCA(otu = otu,map = map,ps = NULL,env = env,group = "Group")
#'
#' @export


# #导入otu表格
# otu = read.delim("../../data/otutab.txt",row.names = 1)
# #导入注释文件
# tax = read.delim("../../data/taxonomy.txt",row.names = 1)
# #导入分组文件
# map = read.delim("../../data/metadata.tsv",row.names = 1)
# head(map)
# 
# # 导入环境因子文件
# env = read.delim("../../data/env.txt",row.names = 1)
# # head(map)
# 
# result =RDA_CCA(otu = otu,map = map,ps = NULL,env = env,group = "Group")
# #=--提取图
# p = result[[1]]
# p
# #--提取作图数据
# plotdata = result[[2]]
# head(plotdata)
# 
# # 提取带有标记的图形
# p2 = result[[3]]
# p2
# # 提取环境因子同群落的差异检测
# table = result[[4]]
# head(table)
# # 如果两种模型都选择，则列表中5，6就被激活了#
# result[[5]]
# result[[6]]



RDA_CCA = function(otu = otutab,map = metadata,ps = NULL,env = env,group = "Group"){
  #----测试代码#----
  # otu = otu
  # map = map
  # ps = NULL
  # env = env
  # group = "Group"
  
  #----所需R包#----
  library("phyloseq")
  library("vegan")
  library("grid")
  library("ggplot2")
  #功能函数#----
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  
  #----数据导入--#-----

  if (is.null(ps) ) {
    
    head(otu)
    otu = as.matrix(otu)
    str(otu)
    
    head(tax)
    tax = as.matrix(tax)
    # taxa_names(tax)
    
    
    colnames(map) = gsub(group,"AA", colnames(map))
    
    map$Group = map$AA
    map$Group = as.factor(map$Group )
    map$Group 
    # #导入进化树
    # tree = read.tree("./otus.tree")
    # tree
    
    ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                   sample_data(map) 
                   # tax_table(tax)
                   # phy_tree(tree)
    )
  }
  
  if (!is.null(ps) ) {
    ps = ps
    
  }
  ps
  # 导入环境因子文件#----
  env.dat = env
  
  # 相对丰度标准化化#----
  ps_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
  #提取otu表格#----
  otu = as.data.frame(t(vegan_otu(ps_rela)))
  #-提取分组文件#----
  mapping = as.data.frame( sample_data(ps_rela))

  # match env and fg datasets
  samp.fg = colnames(otu)
  # 环境因子标准化#----
  env.st = decostand(env.dat, method="standardize", MARGIN=2,na.rm = TRUE)#
  samp.env= rownames(env.st)
  my.env = match(samp.fg, samp.env)
  env.st2 = na.omit(env.st[my.env, ])  # omit the NA rows if without fg data
  samp.env= rownames(env.st2)
  my.fg = match(samp.env, samp.fg)
  otu = otu[, my.fg]
  
  ##without Latitude and Longitude
  env.st3=env.st2
  
 # 专转置OTU表格#----
  otu = t(otu)
  #首先计算DCA排序#---
  DCA= decorana(otu)  
  #--提取DCA排序结果准备一下判断
  xxa = as.data.frame(DCA$rproj)
  max(xxa$DCA1)
  
  #--如果基于3-4之间，则CCA和RDA都可以使用#----
  if(max(xxa$DCA1)<=4 & max(xxa$DCA1)>=3){twochiose = "T"}else{twochiose = "F"}
  
 
  # 判断模型是否是线性模型还是单峰模型#----
  if(max(xxa$DCA1) > 4 | twochiose == "T") {
    ##choise CCA
    C.whole = cca(otu, env.st3)  ##rda(otu, env.st2)
    C.whole
    # for env selection by CCA inflation factors
    #Function vif.cca and alias.cca can be used to analyse linear dependencies among constraints and conditions.
    inf_factor = vif.cca(C.whole)
    
    # delete varable with max inflation factor
    na_env = which(is.na(inf_factor))
    if(isTRUE(length(na_env) > "0") ){
      inf_factor = inf_factor[-na_env]
    }
    
    max_env = which(inf_factor == max(inf_factor))
    env.st4 = env.st3
    while ( inf_factor[max_env] > 20){
      env.st4 = env.st4[,-max_env]
      C.reduced = cca(otu, env.st4)
      inf_factor = vif.cca(C.reduced)
      max_env = which(inf_factor == max(inf_factor))
    }
    output2 = inf_factor ;output2
    env.st4
    
    # for F and p values
    ind.p = array(0,dim=c(1,ncol(env.st4)))
    ind.F = array(0,dim=c(1,ncol(env.st4)))
    for(j in 1:ncol(env.st4)){
      ind.cca = cca(otu, env.st4[,j]) #ind.cca = cca(otu, env.st[,j], env.st[,-j])  #
      ind.sig = anova(ind.cca,step=1000)
      ind.p[1,j] = ind.sig$Pr[1]
      ind.F[1,j] = ind.sig$F[1]
    }
    
    colnames(ind.p) = colnames(env.st4)
    inf_Fp=rbind(output2,ind.F,ind.p)
    row.names(inf_Fp)=c("inf_factor","F","p")
    
    ##重新计算CCA
    C.whole = cca(otu, env.st4)  ##rda(otu, env.st3)
    x.sig = anova(C.whole)
    x.p = x.sig$Pr[1] ;x.p
    x.F = x.sig$F[1]  ;x.F
    F1 <- paste("anova F: ",round(x.F, 2), sep = "")
    pp1 = paste("p: ",round(x.p, 2), sep = "")
    title = paste(F1," ",pp1, sep = "")
    
    output1 = summary(C.whole)
    
    str(output1)
    a=output1$sites;a  ##样本坐标
    b=output1$cont$importance;b ##特征值，解释??? #eigenvals(C.whole)
    c=output1$biplot*5;c  ##环境因子坐标
    

    
    ca1=round(b[2,1],2);ca1
    ca2=round(b[2,2],2);ca2
    
    head(a)
    head(mapping)
    aa = merge(a,mapping, by="row.names",all=F);aa
    c = as.data.frame(c)
    
    library("ggplot2")
    library(ggrepel)
    # mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A","#E6AB02", "#B3DE69")
    p=ggplot()+
      geom_point(data=aa,aes(x=CCA1,y=CCA2, fill = aa$SampleType),pch = 21,colour = "black",size = 4)+
      geom_segment(data = c,aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
                   arrow = arrow(angle=22.5,length = unit(0.2,"cm")),linetype=1, size=0.6,colour = "black")+
      # stat_ellipse( data=aa,linetype = 2,level = 0.65,aes(x=CCA1,y=CCA2,group  =SampleType, colour =  SampleType))+
      geom_text_repel(data = c,aes(x=CCA1,y=CCA2, label = row.names(c)))+
      
      labs(x=paste("CCA 1 (", ca1*100, "%)", sep=""),
           y=paste("CCA 2 (", ca2*100, "%)", sep=""),
           title=title)
    
    p
    p = p+theme_bw()+
      
      #scale_y_continuous(expand = c(0,0))+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
      theme(
        
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        
        plot.title = element_text(vjust = -8.5,hjust = 0.1),
        axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
        axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
        axis.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 14),
        legend.text = element_text(size = 15,face = "bold")
        
        
      ) 
    p
    library(ggrepel) 
    p3 = p+geom_text_repel(data=aa, aes(x=CCA1,y=CCA2,label=aa$Row.names),size=4)#?stat_el
    p3
    
    #如果两者都做需要两个图片都输出#----
    if (twochiose == "T") {
      pRDA= p
      
    } 
    #--选择单峰模型#-----
  } else if (max(xxa$DCA1) < 3| twochiose == "T" ){
    ##choise RDA
    C.whole = rda(otu, env.st2)
    C.whole
    
    
    # for env selection by RDA inflation factors
    
    #Function vif.cca and alias.cca can be used to analyse linear dependencies among constraints and conditions.
    inf_factor = vif.cca(C.whole)
    inf_factor
    
    # delete varable with max inflation factor
    na_env = which(is.na(inf_factor))
    if(isTRUE(length(na_env) > "0") ){
      inf_factor = inf_factor[-na_env]
    }
    
    max_env = which(inf_factor == max(inf_factor))
    env.st4 = env.st3
    while ( inf_factor[max_env] > 20){
      env.st4 = env.st4[,-max_env]
      C.reduced = cca(otu, env.st4)
      inf_factor = vif.cca(C.reduced)
      max_env = which(inf_factor == max(inf_factor))
    }
    output2 = inf_factor ;output2
    
    # for F and p values
    ind.p = array(0,dim=c(1,ncol(env.st4)))
    ind.F = array(0,dim=c(1,ncol(env.st4)))
    for(j in 1:ncol(env.st4)){
      ind.cca = cca(otu, env.st4[,j]) #ind.cca = cca(otu, env.st[,j], env.st[,-j])  #
      ind.sig = anova(ind.cca,step=1000)
      ind.p[1,j] = ind.sig$Pr[1]
      ind.F[1,j] = ind.sig$F[1]
    }
    
    colnames(ind.p) = colnames(env.st4)
    inf_Fp=rbind(output2,ind.F,ind.p)
    row.names(inf_Fp)=c("inf_factor","F","p")
    
    ##重新计算rdA
    C.whole = rda(otu, env.st4)  ##rda(otu, env.st3)
    x.sig = anova(C.whole)
    x.p = x.sig$Pr[1] ;x.p
    x.F = x.sig$F[1]  ;x.F
    F1 <- paste("anova F: ",round(x.F, 2), sep = "")
    pp1 = paste("p: ",round(x.p, 2), sep = "")
    title = paste(F1," ",pp1, sep = "")
    
    output1 = summary(C.whole)
    
    str(output1)
    a=output1$sites;a  ##样本坐标
    b=output1$cont$importance;b ##特征值，解释??? #eigenvals(C.whole)
    c=output1$biplot;c  ##环境因子坐标
    
    
    # write.table(a,file="cca_site.txt",sep="\t",col.names=NA)
    # write.table(b,file="cca_evale.txt",sep="\t",col.names=NA)
    # write.table(c,file="cca_env.txt",sep="\t",col.names=NA)
    #保存环境因子的显著性检验结果
    

    
    
    ca1=round(b[2,1],2);ca1
    ca2=round(b[2,2],2);ca2
    
    # head(a)
    # head(mapping)
    aa = merge(a,mapping, by="row.names",all=F);aa
    c = as.data.frame(c)
    
    library("ggplot2")
    library(ggrepel)
    # mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A","#E6AB02", "#B3DE69")
    p=ggplot()+
      geom_point(data=aa,aes(x=RDA1,y=RDA2, fill = aa$Group),pch = 21,colour = "black",size = 4)+
      geom_segment(data = c,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                   arrow = arrow(angle=22.5,length = unit(0.2,"cm")),linetype=1, size=0.6,colour = "black")+
      # stat_ellipse( data=aa,linetype = 2,level = 0.65,aes(x=RDA1,y=RDA2,group  =SampleType, colour =  SampleType))+
      geom_text_repel(data = c,aes(x=RDA1,y=RDA2, label = row.names(c)))+
      
      labs(x=paste("RDA 1 (", ca1*100, "%)", sep=""),
           y=paste("RDA 2 (", ca2*100, "%)", sep=""),
           title=title)
    
    p
    head(aa)
    p = p+theme_bw()+
      
      #scale_y_continuous(expand = c(0,0))+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
      theme(
        
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        
        plot.title = element_text(vjust = -8.5,hjust = 0.1),
        axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
        axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
        axis.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 14),
        legend.text = element_text(size = 15,face = "bold")
        
        
      ) 
    p
    library(ggrepel) 
    p3 = p+geom_text_repel(data=aa, aes(x=RDA1,y=RDA2,label=aa$Row.names),size=4)#?stat_el
    p3
    
    #如果两者都做需要两个图片都输出CCA#----
    if (twochiose == "T") {
      pCCA= p
      
    }
    
    
    
    
  }
  #--如果只选择一个模型，那么需要设定另一种为无结果#-----
  if(max(xxa$DCA1) > 4|max(xxa$DCA1) < 3){
    pCCA= "no result"
    pRDA= "no result"
  }
  
  return(list(Plot = p,Plotdata = aa,Labelplot = p3,Difftest = inf_Fp,pRDA,pCCA))
  
}


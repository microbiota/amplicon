% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MicroTest.R
\name{MicroTest}
\alias{MicroTest}
\title{Beta diversity statistics by adonis/anosim/MRPP in all groups}
\usage{
MicroTest(
  otu = NULL,
  map = NULL,
  ps = ps,
  group = "Group",
  Micromet = "MRPP",
  dist = "bray"
)
}
\arguments{
\item{otu}{OTU/ASV table;}

\item{map}{Sample metadata;}

\item{ps}{alternative input;}

\item{group}{group ID;}

\item{Micromet}{statistics default by adonis, alternative anosim or MRPP;}

\item{dist}{distance type, including "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski"  "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial"  "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co";}

\item{tree}{tree/nwk file;}

\item{method}{DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA;}

\item{pvalue.cutoff}{Pvalue threshold, default in 0.05;}
}
\value{
stat table
}
\description{
Input phyloseq object, test method and distance type
}
\details{
By default, input phyloseq object include metadata and otutab
The available diversity indices include the following:
\itemize{
\item{most used indices: bray unifrac wunifrac}
\item{other used indices: dpcoa jsd manhattan euclidean canberra kulczynski jaccard gower altGower morisita horn mountford raup binomial chao cao w -1 c wb r I e t me j sor m -2 co}
}
}
\examples{
# Input otutab, metadata, and options group, method and distance
MicroTest(otu = otutab_rare, map = metadata, group = "Group", Micromet = "MRPP", dist = "bray")
# Input phyloseq format input, and options group, method and distance
MicroTest(ps = ps, group = "Group", Micromet = "MRPP", dist = "bray")
}
\references{
Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
}
\seealso{
BetaDiv beta_pcoa beta_cpcoa
}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Yong-Xin Liu \email{yxliu@genetics.ac.cn}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare.R
\name{compare}
\alias{compare}
\title{Difference comparison}
\usage{
compare(
  data = otutab,
  metadata = metadata,
  group = "Group",
  compare_pair = "KO-WT",
  method = "wilcox",
  RA = 0.01,
  pvalue = 0.05,
  fdr = 0.1,
  normalize = T
)
}
\arguments{
\item{data}{data matrix}

\item{metadata}{matrix or dataframe, including sampleID and groupID;}

\item{group}{column name for groupID.}

\item{compare_pair}{paired groups linked with dash.}

\item{method}{method of test.}

\item{RA}{threshold for relative abundance.}

\item{pvalue}{threshold of pvalue.}

\item{fdr}{threshold of fdr.}

\item{normalize}{data normalize to 100\%, default T, turn off by F.}
}
\value{
Difference comparison table.
}
\description{
Input feature table, group info and comparing pair
Features with relative abundance, p-value and FDR for pair.
}
\details{
By default, return table
The available test methods include the following:
\itemize{
\item{most used methods: wilcoxon test, t.test}
\item{other used indices: edgeR, DESeq2}
}
}
\examples{
# Table, metadata, groupID, compare pair and method (wilcoxon test)
compare(data = otutab, metadata = metadata, group = "Group", compare_pair = "KO-WT", method = "wilcox")
# Filter by 0.01\% relative abundance, Pvalue and FDR, not normalize
compare(data = otutab, metadata = metadata, group = "Group", compare_pair = "KO-WT", method = "wilcox", RA = 0.01, pvalue = 0.05, fdr = 0.1, normalize = F)
}
\references{
Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai.
A practical guide to amplicon and metagenomic analysis of microbiome data.
Protein Cell, 2020(41), 1-16, DOI: \url{https://doi.org/10.1007/s13238-020-00724-8}
}
\seealso{
compare_vocalno compare_heatmap compare_manhattan
}
\author{
Contact: Yong-Xin Liu \email{metagenome@126.com}
}

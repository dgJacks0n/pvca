\name{pvca-package}
\alias{pvca-package}
\alias{pvca}
\docType{package}
\title{A package that provides an approach to assess the source of batch effects in a microarray gene expression experiment
}
\description{
This package contains the function to assess the batch sources by fitting all "sources" as random effects including two-way interaction terms in the Mixed Model(depends on lme4 package) to selected principal components, which were obtained from the original data correlation matrix. This package accompanies the book "Batch Effects and Noise in Microarray Experiements, chapter 12.
}
\details{
\tabular{ll}{
Package: \tab pvca\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-09-11\cr
License: \tab LGPL (>= 2.0)\cr
}
library(golubEsets)

data(Golub_Merge)
pct_threshold <- 0.6
batch.factors <- c("ALL.AML", "BM.PB", "Source")

pvcaObj <- pvcaBatchAssess (Golub_Merge, batch.factors, pct_threshold) 
bp <- barplot(pvcaObj$dat,  xlab = "Effects", ylab = "Weighted average proportion variance", ylim= c(0,1.1),col = c("blue"), las=2, main="PVCA estimation bar chart")
axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8) 
print(sessionInfo())

}
\author{
Pierre Bushel  <bushel@niehs.nih.gov>

Maintainer: Jianying LI <li11@niehs.nih.gov>
}
\references{
Batch Effects and Noise in Microarray Experiments: Sources and Solutions. 2009 John Wiley & Sons, Ltd.
}
\keyword{ package }

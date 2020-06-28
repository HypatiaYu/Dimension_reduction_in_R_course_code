#install biobase 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biobase")
#--- load 
library(devtools)
library(Biobase)

par(pch = 19)
tropical = c('darkorange',
             'dodgerblue',
             'hotpink',
             'limegreen',
             'yellow')
palette(tropical)
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")

#extract phenotype data
load(file = con)
close(con)
mp = montpick.eset
pdata = pData(mp)
edata = as.data.frame(exprs(mp))
fdata = fData(mp)
dim(edata)
edata = log2(edata + 1)

#center the data
edata_centered = edata - rowMeans(edata)

svd1 = svd(edata_centered)

#plot singular values
plot(svd1$d, 
     ylab = "Singular values", 
     col = 2)

#plot the variance explained 
plot(svd1$d^2/sum(svd1$d^2),
ylab = "Percent Variance Explained",
col  = 2)

#what is the variable that explains most of the variance
par(mfrow = c(1, 2))

plot(svd1$v[,1],
     col = 2,
     ylab = "1st PC")
plot(svd1$v[,2],
     col = 2,
     ylab = "2nd PC")

#plot both against each other 
plot(svd1$v[,1],
     svd1$v[,2],
     ylab = "2nd PC",
     xlab = "1st PC",
     col = as.numeric(pdata$study))

boxplot(svd1$v[,1] ~ pdata$study,
        border = c(1, 2))

points(svd1$v[,1] ~ jitter(as.numeric(pdata$study)),
       col = as.numeric(pdata$study))

#plot the principle components
pc1 = prcomp(edata)
plot(pc1$rotation[,1],
     svd1$v[,1])

#scale data by subracting the column means
edata_centered2 = t(t(edata) - colMeans(edata))
svd2 = svd(edata_centered2)

plot(pc1$rotation[,1],
     svd2$v[,1])

#outliers can drive decomposition
edata_outlier = edata_centered
edata_outlier[6,] = edata_centered[6,] * 1000

svd3 = svd(edata_outlier)
plot(svd1$v[,1],
     svd3$v[,1],
     xlab = "Without outlier",
     ylab = "With outlier")

#-- set color pallette 

tropical = c('darkorange',
             'dodgerblue',
             'hotpink',
             'limegreen',
             'yellow')
palette(tropical)
par(pch = 19)

#-- load libraries 
library(devtools)
library(Biobase)
library(preprocessCore)

#--load dataset
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file = con)

#--extract phenotype and expression data

mp = montpick.eset
pdata = pData(mp)
edata = as.data.frame(exprs(mp))

edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 3, ]

dim(edata)

#--figure out distribution of each sample

colramp = colorRampPalette(c(3, 
                             "white"
                             , 2))(20)

#--plot density of each sample
plot(density(
  edata[,1]
),
col = colramp[1],
lwd = 3,
ylim = c(0, 
         .30))
for(i in 2:20){
  lines(
    density(
    edata[,i
          ]
  ),
  lwd = 3,
  col = colramp[i])
}

#quantile normalization to force 
#the distributions to be the same

norm_edata = normalize.quantiles(as.matrix(edata))
plot(
  density(
    norm_edata[, 1]),
  col = colramp[1],
  lwd = 3,
  ylim = c(0, .20))
for(i in 2:20){
  lines(
    density(
      norm_edata[,i]
    ),
    lwd = 3,
    col = colramp[i]
  )
}

dim(norm_edata)

#it has not removed gene by gene varibility

plot(norm_edata[1,],
     col = as.numeric(
       pdata$study))

#svd to see differences more clearly

svd1 = svd(norm_edata - rowMeans(norm_edata))
plot(svd1$v[,1], 
     svd1$v[,2],
     col = as.numeric(pdata$study))

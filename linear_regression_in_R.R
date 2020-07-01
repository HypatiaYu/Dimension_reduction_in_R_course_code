tropical = c(
  'darkorange',
  'dodgerblue',
  'hotpink',
  'limegreen',
  'yellow'
)
palette(tropical)
par(pch = 19)

#-- load packages
library(devtools)
library(Biobase)
library(broom)

#-- load dataset
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file = con)
close(con)

#-- define dataset variables
bm = bodymap.eset
pdata = pData(bm)
edata = as.data.frame(exprs(bm))
fdatq = fData(bm)
edata = as.matrix(edata)

lml = lm(edata[1, ] ~ pdata$age)
tidy(lml)

plot(pdata$age,
     edata[1,],
     col = 1)
abline(lml,
      col = 2,
      lwd = 3)

pdata$gender

boxplot(edata[1,] ~pdata$gender)
points(edata[1,] ~ jitter(
  as.numeric(
    pdata$gender
  )),
  col = as.numeric(
    pdata$gender
  )
)

dummy_m = pdata$gender == "M"

lm2 = lm(edata[1,] ~ pdata$gender)
tidy(lm2)

mod2 = model.matrix( ~ pdata$gender)

table(pdata$tissue.type)

pdata$tissue.type == 'adipose'

tidy(lm(edata[1,] ~pdata$tissue.type))

lm3 = lm(edata[1,] ~pdata$age + pdata$gender)
tidy(lm3)

lm4 = lm(edata[1,] ~pdata$age*pdata$gender)
tidy(lm4)

#-- overlay ontop of graph
lm4 = lm(edata[6,] ~pdata$age)
plot(pdata$age,
     edata[6,],
     col = 2)
abline(lm4,
       col = 1,
       lwd = 3)

#--- when outlier can plot 

lm5 = lm(edata[6,] ~ index)
index = 1 : 19
plot(index,
     edata[6,],
     col = 2)
abline(lm5,
       col = 1,
       lwd = 3
       )

#--- remove outlier and refit model

lm6 = lm(edata[6,
               -19] ~ index[-19])
abline(lm6,
       col = 3,
       lwd = 3
       )
legend(5,
       1000,
       c("With outlier",
         "Without outlier"),
       col = c(1, 3),
       lwd = 3)

par(mfrow = c(1, 2))
hist(lm6$residuals,
     col = 2)
hist(lm5$residuals,
     col = 3)

gene1 = log2(edata[1,]+1)
lm7 = lm(gene1 ~ index)
hist(lm7$residuals,
     col = 4)

#-- fit models with coefficeitns
lm8 = lm(gene1 ~ pdata$tissue.type + pdata$age)

tidy(lm8)

dim(model.matrix( ~ pdata$tissue.type + pdata$age))
colramp = colorRampPalette(1 : 4)(17)
lm9 = lm(edata[2,] ~ pdata$age)

plot(lm9$residuals,
     col = colramp[as.numeric(
       pdata$tissue.type
     )])

par(mfrow = c(1, 1))
colramp = colorRampPalette(1 : 4)(17)
lm9 = lm(edata[2, ] ~ pdata$age)

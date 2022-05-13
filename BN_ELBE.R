library(bnlearn)
library(caret)
library(readxl)
library(Rgraphviz)
library(writexl)
library(progress)

## Read the dataset
#setwd("C:/Users/Administer/Documents/GitHub/ElbeMetalAnalysis")
df = read_excel("BN_Elbe.xlsx")
df = as.data.frame(df)
bn_df = df[df$SubCat == 'D', 5:24]
#bn_df = df[5:24]

## Log-transform the data
bn_df = log(bn_df)
N = dim(bn_df)[1]

## Create Blacklist for structure constraint
bl <- data.frame(matrix(ncol = 2))
colnames(bl) <- c("from", "to")
rows = 1

for(i in 1:10)
{
  for(j in 11:20)
  {
    bl[rows, 1] = colnames(bn_df)[i]
    bl[rows, 2] = colnames(bn_df)[j]
    rows = rows + 1
  }
}

for(i in 1:10)
{
  for(j in 1:10)
  {
    bl[rows, 1] = colnames(bn_df)[i]
    bl[rows, 2] = colnames(bn_df)[j]
    rows = rows + 1
  }
}

for(i in 11:20)
{
  for(j in 11:20)
  {
    bl[rows, 1] = colnames(bn_df)[i]
    bl[rows, 2] = colnames(bn_df)[j]
    rows = rows + 1
  }
}

## Tabu-structure learning methods
net = tabu(bn_df, blacklist = bl)

## Strength analysis
strength = arc.strength(net, data=bn_df, criterion="bic-g")
#strength.plot(net, strength, shape='rectangle')
df_nodes = read.csv("d:/net_nodes.csv")
strength = merge(x=strength, y=df_nodes, by.x='from', by.y='x', all.x=TRUE)
names(strength) = c("from", "to", "strength", "edgeCat")
write.csv(strength, "d:/net_arcs.csv")

strengthmat = data.frame(matrix(ncol=11, nrow=10))
colnames(strengthmat) = c('factors', colnames(bn_df)[1:10])
strengthmat$factors = c(colnames(bn_df)[11:20])
for(i in 1:dim(strength)[1])
{
  row = which(strengthmat$factors == strength$from[i])
  col = which(colnames(strengthmat[-1]) == strength$to[i])
  strengthmat[row, col + 1] = abs(strength$strength[i])
}
strengthmat[is.na(strengthmat)] = '-'
#write.csv(strengthmat, "d:/strength.csv")
#write.csv(colnames(bn_df), "d:/net_nodes.csv")
#plot(net)

## Learn the model paramters
dfit = bn.fit(net, bn_df)

## Output the model parameters 
#df = data.frame(Factor = c('(Intercept)', 'As', 'Pb', 'Cd', 'Cr', 'Fe', 'Cu', 'Mn', 'Ni', 'Hg', 'Zn', 'TOC', 'EC', 'TP', 'TN', 'pH', 'RIR', 'TRF', 'ARG', 'NAT', 'POPD', 'RMSE', 'R2'))
df = data.frame(Factor = c('(Intercept)', 'TOC', 'EC', 'TP', 'TN', 'pH', 'RIR', 'TRF', 'ARG', 'NAT', 'POP', 'RMSE', 'R2'))
df$id = 1:nrow(df)
for(i in 1:10)
{
  target = dfit[[i]]$node
  df_tmp = data.frame(dfit[[i]]$coefficients)
  names(df_tmp) = target
  df_tmp$Factor = rownames(df_tmp)
  df = merge(x=df, y=df_tmp, by='Factor', sort = FALSE, all.x=TRUE)
}
df = df[order(df$id), ]
#df = df[-2]
write_xlsx(df, 'D:/bn_result.xlsx')

## 10-fold cross validation
r2list = c()
rmselist = c()
for(i in 1:10)
{ 
  target = dfit[[i]]$node
  mse = bn.cv(bn_df, net, method = "k-fold", loss = 'mse', loss.args = list(target = target))
  r2 =  bn.cv(bn_df, net, method = "k-fold", loss = 'cor', loss.args = list(target = target))
  rmselist = c(rmselist, sqrt(loss(mse)))
  r2list = c(r2list, loss(r2))
}



library(bnlearn)
library(caret)
library(readxl)
library(Rgraphviz)
library(writexl)
library(Metrics)
library(modelr)
#library(foreach)
#library(doParallel)

## Read the dataset
df_raw = read_excel("C:/Users/Administer/Documents/GitHub/ElbeMetalAnalysis/BN_Elbe.xlsx") # Read the data from the Excel table
df_raw = as.data.frame(df_raw)

setwd("C:/Users/Administer/Documents/GitHub/ElbeMetalAnalysis/BN_Unc/10-19") # set writing filepath

#bn_df.raw = df_raw[5:24] # Data in all subbasins
#bn_df.raw = df_raw[df_raw$Period == "10-19", 5:24] # Data in a period
bn_df.raw = df_raw[df_raw$Subcat == "U", 5:24] # Data in a subbasin

## Log-transform the data
bn_df = log(bn_df.raw)

## Create Blacklist for structure constraint
bl <- data.frame(matrix(ncol = 2))
colnames(bl) <- c("from", "to")
rows = 1

# No direction from metal to influential factors
for(i in 1:10)
{
  for(j in 11:20)
  {
    bl[rows, 1] = colnames(bn_df)[i]
    bl[rows, 2] = colnames(bn_df)[j]
    rows = rows + 1
  }
}

# No interactions between influential factors
for(i in 1:10)
{
  for(j in 1:10)
  {
    bl[rows, 1] = colnames(bn_df)[i]
    bl[rows, 2] = colnames(bn_df)[j]
    rows = rows + 1
  }
}

# No interactions between metals
for(i in 11:20)
{
  for(j in 11:20)
  {
    bl[rows, 1] = colnames(bn_df)[i]
    bl[rows, 2] = colnames(bn_df)[j]
    rows = rows + 1
  }
}

# Bootstrapped network structure learning
R.boot = 1000 # the number of bootstrap replicates
m = nrow(bn_df) # the size of bootstrap replicates

# Measure the strength of the probabilistic relationships by Bootstrapping tabu methods
tabu.boot.strength = boot.strength(bn_df, algorithm = "tabu", 
                                   R = R.boot, m = m, 
                                   algorithm.args = list(blacklist = bl))

# Network averaging
threshold = inclusion.threshold(tabu.boot.strength) # Estimate strength threshold
tabu.boot.net = averaged.network(tabu.boot.strength, threshold) # Average network structure
tbl.boot.confidence = tabu.boot.strength[(tabu.boot.strength$strength>threshold)
                                         & (tabu.boot.strength$direction>0.5), ] # Confidence table of the averaged network
tbl.boot.strength = arc.strength(tabu.boot.net, bn_df, criterion = "bic-g")
#bic.boot.net = score(tabu.boot.net, bn_df)

#graphviz.plot(bic.boot.net) # Plot the network

## k-fold cross validation
bn_df.cv = crossv_kfold(bn_df, k = 10)

# Cross-validation for RMSE and R2
rmse = c(NULL)
r2 = c(NULL)
p = c(NULL)
for(k in 1:10)
{
  train.idx = bn_df.cv$train[[k]]$idx
  test.idx = bn_df.cv$test[[k]]$idx
  dfit = bn.fit(tabu.boot.net, bn_df[train.idx, ])
  rmse_tmp = c(NULL)
  r2_tmp = c(NULL)
  p_tmp = c(NULL)
  for(i in 1:10)
  {
    target = dfit[[i]]$node
    pred = predict(dfit, data=bn_df[test.idx, ], node=target)
    rmse_tmp = c(rmse_tmp, Metrics::rmse(bn_df[test.idx, i], pred))
    res.cor = cor.test(bn_df.raw[test.idx, i], pred)
    r2_tmp = c(r2_tmp, res.cor$estimate[1])
    #r2_tmp = c(r2_tmp, stats::cor(bn_df.raw[test.idx, i], pred))
    p_tmp = c(p_tmp, res.cor$p.value)
  }
  rmse = c(rmse, list(rmse_tmp))
  r2 = c(r2, list(r2_tmp))
  p = c(p, list(p_tmp))
}
rmse_mean = rowMeans(sapply(rmse, function(x) x))
rmse_std = apply(sapply(rmse, function(x) x), 1, sd)
r2_mean = rowMeans(sapply(r2, function(x) x))
r2_std = apply(sapply(r2, function(x) x), 1, sd)
p_mean = rowMeans(sapply(p, function(x) x))
p_median = apply(sapply(p, function(x) x), 1, median)

# Conditional coefficient of the fitted model
dfit = bn.fit(tabu.boot.net, bn_df) # learn the model parameters

# Coefficient table
df = data.frame(Factor = c('(Intercept)', 'TOC', 'EC', 'TP', 'TN', 'pH', 'RIR', 'TRF', 'ARG', 'NAT', 'POP', 'RMSE', 'R2'))
df$id = 1:nrow(df)
for(i in 1:10)
{
  target = dfit[[i]]$node
  df_tmp = data.frame(round(dfit[[i]]$coefficients, 2))
  names(df_tmp) = target
  df_tmp$Factor = rownames(df_tmp)
  df = merge(x=df, y=df_tmp, by='Factor', sort = FALSE, all.x=TRUE)
}
df = df[order(df$id), ]
df[12, 3:12] = apply(data.frame(round(rmse_mean, 2), round(rmse_std, 2)), 1, paste, collapse="±")
df[13, 3:12] = apply(data.frame(round(r2_mean, 2), round(r2_std, 2)), 1, paste, collapse="±")

write_xlsx(df, 'bn_result.xlsx')

# Strength table for the CIRCOS diagram
strengthmat = data.frame(matrix(ncol=11, nrow=10))
colnames(strengthmat) = c('factors', colnames(bn_df)[1:10])
strengthmat$factors = c(colnames(bn_df)[11:20])
for(i in 1:dim(tbl.boot.strength)[1])
{
  row = which(strengthmat$factors == tbl.boot.strength$from[i])
  col = which(colnames(strengthmat[-1]) == tbl.boot.strength$to[i])
  strengthmat[row, col + 1] = abs(tbl.boot.strength$strength[i])
}
strengthmat[is.na(strengthmat)] = '-'
write.csv(strengthmat, "bn_strength.csv")

# Confidence table
confidencemat = data.frame(matrix(ncol=11, nrow=10))
colnames(confidencemat) = c('factors', colnames(bn_df)[1:10])
confidencemat$factors = c(colnames(bn_df)[11:20])
for(i in 1:dim(tbl.boot.confidence)[1])
{
  row = which(confidencemat$factors == tbl.boot.confidence$from[i])
  col = which(colnames(confidencemat[-1]) == tbl.boot.confidence$to[i])
  confidencemat[row, col + 1] = abs(tbl.boot.confidence$strength[i])
}
confidencemat[is.na(confidencemat)] = '-'
write.csv(confidencemat, "bn_confidence.csv")



# Uncertainty analysis caused by the data size
# R2 and RMSE
bootstrap.replicate = function(r, net, data, m, blacklist) {
  ## Bootstrapped sampling
  resampling = sample(nrow(data), m, replace = TRUE)
  replicate = data[resampling, , drop = FALSE]
  replicate.cv = modelr::crossv_kfold(replicate, k = 10)
  ## Uncertainty anaylsis
  bic = bnlearn::score(net, replicate, type="bic-g")
  
  rmse = c(NULL)
  r2 = c(NULL)
  for(k in 1:10)
  {
    train.idx = replicate.cv$train[[k]]$idx
    test.idx = replicate.cv$test[[k]]$idx
    rmse_tmp = c(NULL)
    r2_tmp = c(NULL)
    dfit = bnlearn::bn.fit(net, replicate[train.idx, ])
    for(i in 1:10)
    {
      target = dfit[[i]]$node
      res.predict = predict(dfit,data=replicate[test.idx, ],node=target)
      rmse_tmp = c(rmse_tmp, Metrics::rmse(replicate[test.idx, i], res.predict))
      res.cor = stats::cor.test(replicate[test.idx, i], res.predict)
      r2_tmp = c(r2_tmp, res.cor$estimate[1])
    }
    rmse = c(rmse, list(rmse_tmp))
    r2 = c(r2, list(r2_tmp))
  }
  
  rmse_mean = rowMeans(sapply(rmse, function(x) x))
  rmse_std = apply(sapply(rmse, function(x) x), 1, sd)
  r2_mean = rowMeans(sapply(r2, function(x) x))
  r2_std = apply(sapply(r2, function(x) x), 1, sd)

  res = list(m=m, bic=bic, rmse=rmse_mean, rmse_std=rmse_std, r2=r2_mean, r2_std = r2_std)
  return(res)
}

R.boot = 1000
interval = 200 # Interval for changing size of replicates

res.boot = data.frame(matrix(0, nrow = floor(nrow(bn_df)/interval) + 1, ncol = 22))
colnames(res.boot) = c("m", "bic", 
                       "As_R2", "Pb_R2", "Cd_R2", "Cr_R2", "Fe_R2", "Cu_R2", "Mn_R2", "Ni_R2", "Hg_R2", "Zn_R2", 
                       "As_RMSE", "Pb_RMSE", "Cd_RMSE", "Cr_RMSE", "Fe_RMSE", "Cu_RMSE", "Mn_RMSE", "Ni_RMSE", "Hg_RMSE", "Zn_RMSE")
row = 1
res.boot$m[row] = nrow(bn_df)
res.boot$bic[row] = score(tabu.boot.net, bn_df)
res.boot[row, 3:12] = apply(data.frame(round(r2_mean, 2), round(r2_std, 2)), 1, paste, collapse="±")
res.boot[row, 13:22] = apply(data.frame(round(rmse_mean, 2), round(rmse_std, 2)), 1, paste, collapse="±")


cluster = parallel::makeCluster(16)
for(m in rev(seq(interval,nrow(bn_df),by=interval)))
{
  row = row + 1
  res = as.list(seq(R))
  res = parallel::parLapplyLB(cluster, res, bootstrap.replicate, 
                              net = tabu.boot.net, data = bn_df, m = m, blacklist = bl)
  res.boot$m[row] = m
  res.boot$bic[row] = mean(sapply(res, function(x) x$bic))
  
  r2_mean = rowMeans(sapply(res, function(x) x$r2))
  r2_std = rowMeans(sapply(res, function(x) x$r2_std))
  rmse_mean = rowMeans(sapply(res, function(x) x$rmse))
  rmse_std = rowMeans(sapply(res, function(x) x$rmse_std))
  
  res.boot[row, 3:12] = apply(data.frame(round(r2_mean, 2), round(r2_std, 2)), 1, paste, collapse="±")
  res.boot[row, 13:22] = apply(data.frame(round(rmse_mean, 2), round(rmse_std, 2)), 1, paste, collapse="±")
}
parallel::stopCluster(cluster)

write_xlsx(res.boot, "bn_unc_evaluation_metrics.xlsx")


# confidence and strength
R.boot = 1000
interval = 200 # Interval for changing size of replicates

unc_confidence = tbl.boot.confidence[, 1:3]
unc_direction = tbl.boot.confidence[, c(1, 2, 4)]
unc_strength = tbl.boot.strength[, 1:3]

for(m in rev(seq(interval,nrow(bn_df),by=interval)))
{
  row = row + 1
  unc_tabu.boot.strength = boot.strength(bn_df, algorithm = "tabu", 
                                         R = R.boot, m = m, 
                                         algorithm.args = list(blacklist = bl))
  unc_threshold = inclusion.threshold(unc_tabu.boot.strength)
  unc_tabu.boot.net = averaged.network(unc_tabu.boot.strength, unc_threshold)
  unc_tbl.boot.strength = arc.strength(unc_tabu.boot.net, bn_df, criterion = "bic-g")

  unc_confidence = merge(unc_confidence, unc_tabu.boot.strength[, 1:3], by = c('from', 'to'), all.x = T)
  unc_direction = merge(unc_direction, unc_tabu.boot.strength[, c(1, 2, 4)], by = c('from', 'to'), all.x = T)
  unc_strength = merge(unc_strength, unc_tbl.boot.strength[, 1:3], by = c('from', 'to'), all.x = T)
}

write.csv(unc_confidence, "bn_unc_confidence.csv")
write.csv(unc_direction, "bn_unc_direction.csv")
write.csv(unc_strength, "bn_unc_strength.csv")

# Influence of different subbasins
subbasins = c("U", "M", "D")
tabu.boot.strength.list = c(NULL)
for(subbasin in subbasins)
{
  cat(crayon::green("Selection of Subbasin:", subbasin, "\n"))
  bn_df.raw = df_raw[df_raw$Period == period, 5:24] # Data in each subbasin
  bn_df = log(bn_df.raw)
  
  # Bootstrapped network structure learning
  R.boot = 1000 # the number of bootstrap replicates
  m = nrow(bn_df) # the size of bootstrap replicates
  
  tabu.boot.strength.sb = boot.strength(bn_df, algorithm = "tabu", 
                                        R = R.boot, m = m, 
                                        algorithm.args = list(blacklist = bl))
  tabu.boot.strength.list = c(tabu.boot.strength.list, list(tabu.boot.strength.sb))
}

compare.arc = data.frame(threshold = 0, U = 0, M = 0, D = 0)
compare.shd = data.frame(threshold = 0, U = 0, M = 0, D = 0)
row = 1

for(threshold.i in seq(0,1,0.05))
{
  compare.arc[row, 1] = threshold.i
  compare.shd[row, 1] = threshold.i
  for(i in 1:3)
  {
    tabu.boot.net.all = averaged.network(tabu.boot.strength, threshold.i)
    tabu.boot.net.sb = averaged.network(tabu.boot.strength.list[[i]], threshold.i)
    compare.res = compare(tabu.boot.net.all, tabu.boot.net.sb) # target, current; 
    TP = compare.res$tp
    FP = compare.res$fp
    FN = compare.res$fn
    TN = nrow(tabu.boot.strength) - nrow(bl) - TP - FP - FN
    a = nrow(tabu.boot.net.all)
    N = length(colnames(bn_df))
    i = N * (N - 1) / 2 - a
    
    
    compare.arc[row, i + 1] = 1 - shd(tabu.boot.net.all, tabu.boot.net.sb) / (nrow(tabu.boot.strength) - nrow(bl))
    compare.shd[row, i + 1] = shd(tabu.boot.net.all, tabu.boot.net.sb)
  }
  row = row + 1
}
compare.arc
compare.shd

# 90-99, 00-09, 10-19
periods = c("90-99", "00-09", "10-19")
tabu.boot.strength.list = c(NULL)
for(period in periods)
{
  cat(crayon::green("Selection of Subbasin:", period, "\n"))
  bn_df.raw = df_raw[df_raw$Period == period, 5:24] # Data in each subbasin
  bn_df = log(bn_df.raw)
  
  # Bootstrapped network structure learning
  R.boot = 1000 # the number of bootstrap replicates
  m = nrow(bn_df) # the size of bootstrap replicates
  
  tabu.boot.strength.p = boot.strength(bn_df, algorithm = "tabu", 
                                      R = R.boot, m = m, 
                                      algorithm.args = list(blacklist = bl))
  tabu.boot.strength.list = c(tabu.boot.strength.list, list(tabu.boot.strength.p))
}

compare.arc = data.frame(threshold = 0, p1 = 0, p2 = 0, p3 = 0)
compare.shd = data.frame(threshold = 0, p1 = 0, p2 = 0, p3 = 0)
row = 1

for(threshold.i in seq(0,1,0.05))
{
  compare.arc[row, 1] = threshold.i
  compare.shd[row, 1] = threshold.i
  for(i in 1:3)
  {
    tabu.boot.net.all = averaged.network(tabu.boot.strength, threshold.i)
    tabu.boot.net.p = averaged.network(tabu.boot.strength.list[[i]], threshold.i)
    compare.res = compare(tabu.boot.net.all, tabu.boot.net.p) # target, current; 

    compare.arc[row, i + 1] = 1 - shd(tabu.boot.net.all, tabu.boot.net.p) / (nrow(tabu.boot.strength) - nrow(bl))
    compare.shd[row, i + 1] = shd(tabu.boot.net.all, tabu.boot.net.p)
  }
  row = row + 1
}
compare.arc
compare.shd


# dn_df.raw = df_raw[df_raw$SubCat == "M", 5:24] # Data in each subbasin
# bn_df.raw = df_raw[5:24] # Data in all subbasins
# bn_df = log(bn_df.raw)
# 
# R.boot = 1000
# interval = 200 # Interval for changing size of replicates
# c(nrow(bn_df), rev(seq(interval,nrow(bn_df),by=interval)))

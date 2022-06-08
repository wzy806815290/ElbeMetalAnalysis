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
df = read_excel("C:/Users/Administer/Documents/GitHub/ElbeMetalAnalysis/BN_Elbe.xlsx") # Read the data from the Excel table
df = as.data.frame(df)

subbasins = c("all", "U", "M", "D")

for(subbasin in subbasins)
{
  print(subbasin)
  filepath = paste0("C:/Users/Administer/Documents/GitHub/ElbeMetalAnalysis/BN_Unc/", subbasin)
  setwd(filepath) # set writing filepath
  if(subbasin == "all") 
  {
    bn_df = df[5:24] # Data in all subbasins
  }
  else
  {
    bn_df = df[df$SubCat == subbasin, 5:24] # Data in each subbasin
  }

  ## Log-transform the data
  bn_df = log(bn_df)
  
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
  
  graphviz.plot(bic.boot.net) # Plot the network
  
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
      res.cor = cor.test(bn_df[test.idx, i], pred)
      r2_tmp = c(r2_tmp, res.cor$estimate[1])
      #r2_tmp = c(r2_tmp, stats::cor(bn_df[test.idx, i], pred))
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
  R.boot = 1000
  
  if(subbasin == "M")
  {
    interval = 100 # Interval for changing size of replicates
  }
  else
  {
    interval = 200 # Interval for changing size of replicates
  }

  res.boot = data.frame(matrix(0, nrow = floor(nrow(bn_df)/interval) + 1, ncol = 23))
  colnames(res.boot) = c("m", "bic", "threshold",
                         "As_R2", "Pb_R2", "Cd_R2", "Cr_R2", "Fe_R2", "Cu_R2", "Mn_R2", "Ni_R2", "Hg_R2", "Zn_R2", 
                         "As_RMSE", "Pb_RMSE", "Cd_RMSE", "Cr_RMSE", "Fe_RMSE", "Cu_RMSE", "Mn_RMSE", "Ni_RMSE", "Hg_RMSE", "Zn_RMSE")
  
  row = 1
  res.boot$m[row] = nrow(bn_df)
  res.boot$bic[row] = score(tabu.boot.net, bn_df)
  res.boot$threshold[row] = threshold
  res.boot[row, 4:13] = round(r2_mean, 2)
  res.boot[row, 14:23] = round(rmse_std, 2)
  
  unc_confidence = tbl.boot.confidence[, 1:3]
  unc_direction = tbl.boot.confidence[, c(1, 2, 4)]
  unc_strength = tbl.boot.strength[, 1:3]
  
  tbl.boot.confidence = tabu.boot.strength[(tabu.boot.strength$strength>threshold)
                                           & (tabu.boot.strength$direction>0.5), ] # Confidence table of the averaged network
  tbl.boot.strength = arc.strength(tabu.boot.net, bn_df, criterion = "bic-g")
  
  for(m in rev(seq(interval,nrow(bn_df),by=interval)))
  {
    row = row + 1
    unc_tabu.boot.strength = boot.strength(bn_df, algorithm = "tabu", 
                                           R = R.boot, m = m, 
                                           algorithm.args = list(blacklist = bl))
    unc_threshold = inclusion.threshold(unc_tabu.boot.strength)
    unc_tabu.boot.net = averaged.network(unc_tabu.boot.strength, unc_threshold)
    unc_tbl.boot.strength = arc.strength(unc_tabu.boot.net, bn_df, criterion = "bic-g")
    
    res.boot$m[row] = m
    res.boot$bic[row] = score(unc_tabu.boot.net, bn_df)
    
    bn_df.cv = crossv_kfold(bn_df, k = 10)
    
    for(k in 1:10)
    {
      train.idx = bn_df.cv$train[[k]]$idx
      test.idx = bn_df.cv$test[[k]]$idx
      dfit = bn.fit(unc_tabu.boot.net, bn_df[train.idx, ])
      rmse_tmp = c(NULL)
      r2_tmp = c(NULL)
      p_tmp = c(NULL)
      for(i in 1:10)
      {
        print(paste(m, k, i))
        target = dfit[[i]]$node
        pred = predict(dfit, data=bn_df[test.idx, ], node=target)
        rmse_tmp = c(rmse_tmp, Metrics::rmse(bn_df[test.idx, i], pred))
        res.cor = cor.test(bn_df[test.idx, i], pred)
        r2_tmp = c(r2_tmp, res.cor$estimate[1])
        #r2_tmp = c(r2_tmp, stats::cor(bn_df[test.idx, i], pred))
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
    res.boot[row, 4:13] = round(r2_mean, 2)
    res.boot[row, 14:23] = round(rmse_std, 2)
    
    unc_confidence = merge(unc_confidence, unc_tabu.boot.strength[, 1:3], by = c('from', 'to'), all.x = T)
    unc_direction = merge(unc_direction, unc_tabu.boot.strength[, c(1, 2, 4)], by = c('from', 'to'), all.x = T)
    unc_strength = merge(unc_strength, unc_tbl.boot.strength[, 1:3], by = c('from', 'to'), all.x = T)
  }
  
  write.csv(res.boot, "bn_unc_summary.csv")
  write.csv(unc_confidence, "bn_unc_confidence.csv")
  write.csv(unc_direction, "bn_unc_direction.csv")
  write.csv(unc_strength, "bn_unc_strength.csv")
}


  



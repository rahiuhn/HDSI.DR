#' The functions in this file are to be used for hyperparameter optimisation
#' Hyperparameters: q (k), n(number of components), cint, varmax

# m_createFolds = memoise::memoise(caret::createFolds)

folds_gen = function(df){
  # 5 fold OOB for 5 trials
  ## Create Five trials
  samples = nrow(df)
  trials=lapply(1:3, function(x) {set.seed(x); fold= caret::createFolds(1:samples, k=5)})
  trial_list = purrr::flatten(trials)
  trial_vector = unlist(trials)
  return(list(cv_list = trial_list, cv_vector = trial_vector))
}

#'@export
optim_hyperpara = function(optimtype = c("GA", "Grid"), cv= T, df_p, k = seq(5,varnum-1), ncomp= seq(2,10), summary_ci = seq(0.0, 0.99, 0.01), varmax =  seq(-2000, 2000), otherpara){
  if(cv){
    df=df_p
    # Generate CV list
    fold_list = folds_gen(df=df)
    trial_list = fold_list$cv_list
    trial_vector = fold_list$cv_vector
  }

  #optimisation
  # cat(k, ":", ncomp, ":", summary_ci, ":", varmax)
  paradf= expand.grid(k=k, ncomp=ncomp, summary_ci = summary_ci, varmax=varmax, stringsAsFactors = F)

  ## Constant Para
  {
   seeder = otherpara$seeder; i_numb=otherpara$i_numb; effectsize = otherpara$effectsize; model_perf = otherpara$model_perf; coarse_FS = otherpara$coarse_FS; Fine_FS = otherpara$Fine_FS;  strict_anova = otherpara$strict_anova; semi_strict = otherpara$semi_strict; kmean_type = otherpara$kmean_type; predict_method = otherpara$predict_method
  }

  if(optimtype == "GA"){
    nbits= length(GA::decimal2binary(nrow(paradf)))

    GA_OF=function(bin){
      dec_code = GA::binary2decimal(bin)
      # print(c(bin,dec_code))
      if(dec_code <= nrow(paradf)){
        para_k = paradf[dec_code, "k"]; para_n = paradf[dec_code, "ncomp"]; para_ci = paradf[dec_code, "summary_ci"]; para_vm = paradf[dec_code, "varmax"];
        if(cv){
          res_cv = future.apply::future_lapply(trial_list, function(x) {
            # cat(which(trial_list %in% list(x))," ")
            # Do training
            train_df <- df[-x,]
            # str(train_df[,1])
            test_df <- df[x,]
            cv_df = list(train_df, test_df)
            # print(kmean_type)
            fitvalue = Prop_mod(k=para_k, ncomp=para_n, summary_ci = para_ci, varmax=para_vm, df_p = cv_df, seeder=seeder,  i_numb=i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS,  strict_anova = strict_anova, semi_strict = semi_strict, prediction_type = "rmse", predict_method = predict_method, kmean_type = kmean_type)
            # print(fitvalue[[1]])
            res= - fitvalue[[1]]$RMSE
            #print(res)
            res
          }, future.seed = T)
          res=mean(unlist(res_cv), na.rm = T)
          # print(res)
        }
        else{fitvalue = Prop_mod(k=para_k, ncomp=para_n, summary_ci = para_ci, varmax=para_vm, df_p = df_p, seeder=seeder,  i_numb=i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS,  strict_anova = strict_anova, semi_strict = semi_strict, prediction_type = "rmse", predict_method = "aridge")
        res= - fitvalue[[1]]$RMSE
        }
        if(is.nan(res)){res=-10000}
      }else(res=-20000)

      # cat(round(res,2))
      return(Score=res)
    }
    # print(nbits)
    set.seed(3)
    GA <- GA::ga(type = "binary", fitness = GA_OF, nBits = nbits, popSize = 10, maxiter = 150, run=10, pcrossover = 0.8, pmutation = 0.5, parallel = F, monitor = T) #, suggestions = c(GA::decimal2binary(62703,nbits)) #150
    best_row = GA::binary2decimal(GA@solution[1,])
    best_para = paradf[best_row,]
  }
  else{
    Grid_fun = function(bin){
      para_k = paradf[bin, "k"]; para_n = paradf[bin, "ncomp"]; para_ci = paradf[bin, "summary_ci"]; para_vm = paradf[bin, "varmax"];
      if(cv){
        res_cv = future.apply::future_lapply(trial_list, function(x) {
          # Do training
          train_df <- df[-x,]
          # str(train_df)
          test_df <- df[x,]
          cv_df = list(train_df, test_df)
          fitvalue = Prop_mod(k=para_k, ncomp=para_n, summary_ci = para_ci, varmax=para_vm, df_p = cv_df, seeder=seeder,  i_numb=i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS,  strict_anova = strict_anova, semi_strict = semi_strict, prediction_type = "gmse", predict_method = "reg")
          res= - fitvalue[[1]]$RMSE
        }, future.seed = T)
        res=mean(unlist(res_cv), na.rm = T)
      }
      else{fitvalue = Prop_mod(k=para_k, ncomp=para_n, summary_ci = para_ci, varmax=para_vm, df_p = df_p, seeder=seeder,  i_numb=i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS,  strict_anova = strict_anova, semi_strict = semi_strict, prediction_type = "gmse", predict_method = "reg")
      res= - fitvalue[[1]]$RMSE
      }
      if(is.nan(res)){res=-10}
      #print(c(round(dec_code,0), round(res,2)))
      return(Score=res)
    }
    Grid_optim = future.apply::future_lapply(1:nrow(paradf), Grid_fun, future.seed=T)
    optim_row = which (Grid_optim == max(unlist(Grid_optim), na.rm = T))
    best_para = paradf[optim_row,]
  }
  return(best_para)
}

#'@export
cv_hyperpara = function(datatype=c("simulated"), param = list(miss = 33, corr = 0.9, cutter = 10, seed =1, data_code=1, test_train_ratio=0.2, varnum = 15, setting= "Correlation", main_var=10, var_effect=c(0.5, -0.5), correlation_var=15, correlation_val=5, high_dim=T, train_sample=500, var = "Mar"), seeder=1, i_numb=2, effectsize=13, model_perf="None", coarse_FS = T, Fine_FS =T, strict_anova=F, semi_strict=F, optimtype = "GA", cv = T, kmean_type = "Normal", predict_method = "aridge"){
  # Prepare the dataset
  df_p = data_fit(datatype = datatype, param = param)
  varnum = ncol(df_p[[1]])-1

  # Do para optimisation
  if(cv){df_p = df_p[[1]]}
  otherpara = list(seeder = seeder, i_numb = i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS, strict_anova=strict_anova, semi_strict=semi_strict, kmean_type=kmean_type, predict_method=predict_method)

  best_para = optim_hyperpara(optimtype = optimtype, cv= cv, df_p = df_p, k = seq(5,varnum-1), ncomp= seq(2,10), summary_ci = seq(0.0, 0.99, 0.01), varmax =  seq(1, 10), otherpara = otherpara)

  return(best_para)
}

#'@export
mloptim_hyperpara = function(optimtype = c("GA", "Grid"), cv= T, df_p, k = seq(5,varnum-1), ncomp= seq(2,10), summary_ci = seq(0.0, 0.99, 0.01), varmax =  seq(-2000, 2000), otherpara){
  if(cv){
    df=df_p
    # Generate CV list
    fold_list = folds_gen(df=df)
    trial_list = fold_list$cv_list
    trial_vector = fold_list$cv_vector
  }

  #optimisation
  # cat(k, ":", ncomp, ":", summary_ci, ":", varmax)
  paradf= expand.grid(k=k, ncomp=ncomp, summary_ci = summary_ci, varmax=varmax, stringsAsFactors = F)

  ## Constant Para
  {
    seeder = otherpara$seeder; i_numb=otherpara$i_numb; effectsize = otherpara$effectsize; model_perf = otherpara$model_perf; coarse_FS = otherpara$coarse_FS; Fine_FS = otherpara$Fine_FS;  strict_anova = otherpara$strict_anova; semi_strict = otherpara$semi_strict; kmean_type = otherpara$kmean_type; predict_method = otherpara$predict_method
  }

  if(optimtype == "GA"){
    nbits= length(GA::decimal2binary(nrow(paradf)))
    mldf <<- data.frame(q=integer(), lf=integer(), cl=integer(), cs=integer(), perf=integer())

    GA_OF=function(bin){
      dec_code = GA::binary2decimal(bin)
      # print(c(bin,dec_code))
      if(dec_code <= nrow(paradf)){
        para_k = paradf[dec_code, "k"]; para_n = paradf[dec_code, "ncomp"]; para_ci = paradf[dec_code, "summary_ci"]; para_vm = paradf[dec_code, "varmax"];

      qi = mlvalue = 0
      # ML based prediction
      if(nrow(mldf)>20){
        testdf = data.frame(q=para_k, lf = para_n, cl=para_ci, cs=para_vm)
        mlvalue = AIFS:::predict_function(df = testdf, ml_model=ml, ml_tech = "rf")
        qi = quantile(mldf$perf, .75, na.rm = T)
        res = mlvalue
        # cat(qi, " ")
        # print(mlvalue)
      }
      rand  = sample(c(T,F),1,prob = c(0.1,0.9))

      # Perform RHDSI
      if(mlvalue>=qi | rand){

        if(cv){
          # Create the Cross-Validation
          res_cv = future.apply::future_lapply(trial_list, function(x) {
            # cat(which(trial_list %in% list(x))," ")
            # Do training
            train_df <- df[-x,]
            # str(train_df[,1])
            test_df <- df[x,]
            cv_df = list(train_df, test_df)
            # print(kmean_type)
            fitvalue = Prop_mod(k=para_k, ncomp=para_n, summary_ci = para_ci, varmax=para_vm, df_p = cv_df, seeder=seeder,  i_numb=i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS,  strict_anova = strict_anova, semi_strict = semi_strict, prediction_type = "rmse", predict_method = predict_method, kmean_type = kmean_type)
            # print(fitvalue[[1]])
            res= - fitvalue[[1]]$RMSE
            #print(res)
            res
          }, future.seed = T)
          res=mean(unlist(res_cv), na.rm = T)
        }
        else{fitvalue = Prop_mod(k=para_k, ncomp=para_n, summary_ci = para_ci, varmax=para_vm, df_p = df_p, seeder=seeder,  i_numb=i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS,  strict_anova = strict_anova, semi_strict = semi_strict, prediction_type = "rmse", predict_method = "aridge")
        res= - fitvalue[[1]]$RMSE
        }
        mldf <<- rbind(mldf,c(para_k, para_n, para_ci, para_vm, res))
        names(mldf) <<- c("q", "lf", "cl","cs", "perf")
        # print(nrow(mldf))
      }

      # Perform ML Training
      if(is.nan(res)){res=-10000}
      else{
        # ML Training: If nrows is multiple of 20
        if(all.equal(nrow(mldf)/5, as.integer(nrow(mldf)/5)) == T & nrow(mldf)>19){
          ml <<- AIFS:::ML_funct(ml_tech = "rf", df=mldf, seed=2, init_wt=NULL)
          # str(ml)
        }
      }
      }
      else(res=-20000)

      # cat(round(res,2))
      return(Score=res)
    }
    # print(nbits)
    set.seed(3)
    GA <- GA::ga(type = "binary", fitness = GA_OF, nBits = nbits, popSize = 50, maxiter = 150, run=15, pcrossover = 0.8, pmutation = 0.5, parallel = F, monitor = T) #, suggestions = c(GA::decimal2binary(62703,nbits)) #150
    best_row = GA::binary2decimal(GA@solution[1,])
    best_para = paradf[best_row,]
  }
  else{
    Grid_fun = function(bin){
      para_k = paradf[bin, "k"]; para_n = paradf[bin, "ncomp"]; para_ci = paradf[bin, "summary_ci"]; para_vm = paradf[bin, "varmax"];
      if(cv){
        res_cv = future.apply::future_lapply(trial_list, function(x) {
          # Do training
          train_df <- df[-x,]
          # str(train_df)
          test_df <- df[x,]
          cv_df = list(train_df, test_df)
          fitvalue = Prop_mod(k=para_k, ncomp=para_n, summary_ci = para_ci, varmax=para_vm, df_p = cv_df, seeder=seeder,  i_numb=i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS,  strict_anova = strict_anova, semi_strict = semi_strict, prediction_type = "gmse", predict_method = "reg")
          res= - fitvalue[[1]]$RMSE
        }, future.seed = T)
        res=mean(unlist(res_cv), na.rm = T)
      }
      else{fitvalue = Prop_mod(k=para_k, ncomp=para_n, summary_ci = para_ci, varmax=para_vm, df_p = df_p, seeder=seeder,  i_numb=i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS,  strict_anova = strict_anova, semi_strict = semi_strict, prediction_type = "gmse", predict_method = "reg")
      res= - fitvalue[[1]]$RMSE
      }
      if(is.nan(res)){res=-10}
      #print(c(round(dec_code,0), round(res,2)))
      return(Score=res)
    }
    Grid_optim = future.apply::future_lapply(1:nrow(paradf), Grid_fun, future.seed=T)
    optim_row = which (Grid_optim == max(unlist(Grid_optim), na.rm = T))
    best_para = paradf[optim_row,]
  }
  return(best_para)
}

#'@export
mlcv_hyperpara = function(datatype=c("simulated"), param = list(miss = 33, corr = 0.9, cutter = 10, seed =1, data_code=1, test_train_ratio=0.2, varnum = 15, setting= "Correlation", main_var=10, var_effect=c(0.5, -0.5), correlation_var=15, correlation_val=5, high_dim=T, train_sample=500, var = "Mar"), seeder=1, i_numb=2, effectsize=13, model_perf="None", coarse_FS = T, Fine_FS =T, strict_anova=F, semi_strict=F, optimtype = "GA", cv = T, kmean_type = "Normal", predict_method = "aridge"){
  # Prepare the dataset
  df_p = data_fit(datatype = datatype, param = param)
  varnum = ncol(df_p[[1]])-1

  # Do para optimisation
  if(cv){df_p = df_p[[1]]}
  otherpara = list(seeder = seeder, i_numb = i_numb, effectsize = effectsize, model_perf = model_perf, coarse_FS = coarse_FS, Fine_FS = Fine_FS, strict_anova=strict_anova, semi_strict=semi_strict, kmean_type=kmean_type, predict_method=predict_method)
  varmax = seq(-varnum,varnum) #seq(1, 10)
  ncomp = seq(2,min(floor(nrow(df_p)*0.8/5), varnum)) #seq(2,10)

  best_para = mloptim_hyperpara(optimtype = optimtype, cv= cv, df_p = df_p, k = seq(5,varnum-1), ncomp= ncomp, summary_ci = seq(0.0, 0.99, 0.01), varmax = varmax , otherpara = otherpara)

  return(best_para)
}

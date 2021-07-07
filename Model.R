#' Prepare the model
#'
#'@ df_p: List containing training and test dataset. First element in list should be training dataset and second element in the list should be test dataset
#'@ k: This is "q" ,i.e., number of features to sample
#'@ ncomp: Number of latent features
#'@ summary_ci: confidence level (CL) used to filter the features in step 1. It can take value in between 0 and 1.
#'@ varmax: Number of features to be removed or added to the target cluster. It is "cs -  smallest mean CV value". 
#'@ seeder: It is used to control the randomness in the model 
#'@ i_numb: Level of interactions. Currently, it is fixed at 2.It is experimental and its values should not be changed. 
#'@ effectsize: What is the expected effectsize of features. "Large", "Small" and "Medium". It estimates the number of times a feature needs to be sampled. The code also allows to input user defined numeric value.
#'@ coarse_FS: if TRUE, stage two of feature selection, i.e., unsupervised learning is performed
#'@ Fine_FS: if TRUE, stage three of feature selection, i.e., supervised learning is performed
#'@ prediction_type: Prediction metric used for analysis. Currently, it is fixed to root mean squre error (RMSE)
#'@ predict_method: The technique used to build the final predictive model. Two options are defined: regression ("reg") and adaptive ridge regression ("aridge")
#'@ model_perf, strict_anova, semi_strict, kmean_type: Some experimental parameters explored. The values should not be changed.

#'@import plyr
Prop_mod = function(df_p, seeder=1, ncomp=3, i_numb=2, k=5, effectsize=13, model_perf="None", coarse_FS = T, Fine_FS =T, summary_ci = 0.95, strict_anova=F, semi_strict=F, varmax=10, prediction_type = "rmse", predict_method = "reg", kmean_type = "Normal"){
  # cat("ok")

  # Run proposed methods
  {
    # Generate methodlist
    df=df_p[[1]]
    outvar_col=ncol(df)
    methodlist=list(PLS=plsregfunct)
    {
      # Generate the bootstraps
      bootlist=mbootsample(type="normal", p=ncol(df)-1,k=k, rows=nrow(df), interaction_numb = i_numb,
                           plist = names(df)[-outvar_col], inputdf=df,seed_multiplier=seeder, effectsize = effectsize)

      ## Create technique result for each bootstrap
      res = lapply(bootlist, function(x) {
        inputdf=df[x[[2]],c(x[[1]],"y")]
        methodlist[[1]](df=inputdf, ncomp=ncomp, interactionterm=T, int_term=i_numb, model_perf=model_perf)
      })
      res_df=do.call(rbind,res)
      res_df$beta=as.numeric(res_df$beta)
      # str(res_df)
    }
  }

  para=c(summary_ci=summary_ci)
  # Coarse Variable Selection
  {
    if(coarse_FS){
      # Feature Selection Metric Beta Coefficient
      res_summary_beta = plyr::ddply(res_df[,c("Variable","imp")], .(Variable), function(x) {a=f_summary_beta(x$imp, cint=as.numeric(para['summary_ci']));unlist(a)})
      beta_sel_feature_list=variable_selection(df=res_summary_beta)
      sel_feature_list=beta_sel_feature_list
    }
    else{sel_feature_list = res_summary_beta = NA}
    # print(sel_feature_list)
  }

  # Fine Variable Selection
  {
    if(Fine_FS & coarse_FS){
      # Select Variable based on dispersion
      # print(T)
      if(length(sel_feature_list)>0){
        sel_feature_list=cv_filter(pre_filter_list = sel_feature_list, variable_summary_df=res_summary_beta, res_df=res_df,semi_strict=semi_strict, varmax=varmax, strict_anova=strict_anova, totalsteps = (ncol(df)-1), clusters = floor(2+(k/100)), inputdf = df, kmean_type = kmean_type)}
        # print(sel_feature_list)
    }
    else if (Fine_FS){
      sel_feature_list=cv_filter(pre_filter_list = NA, variable_summary_df=NA, res_df=res_df, cint=summary_ci,semi_strict=semi_strict, varmax=10, strict_anova=strict_anova, totalsteps = (ncol(df)-1), clusters = 10+(k/10), inputdf = df, kmean_type = kmean_type)
    }
    else{sel_feature_list}
  }

  ## Extract Results

  #flist <<- sel_feature_list
  # print(length(sel_feature_list))
  if(length(sel_feature_list)==0){ sel_feature_list=NA} #print(T);
  # print(sel_feature_list)
  if(is.na(sel_feature_list)){
    final_marvar=final_var=sel_feature_list
    MV=IV=0
    df_p_final = df_p
  }
  else{
    final_marvar=unique(unlist(strsplit(sel_feature_list, "_")))
    final_var=union(final_marvar, sel_feature_list)
    IV=length(grep("_",final_var))
    MV=length(final_var)-IV
    df_p_final = lapply(df_p, function(x) x[, c(final_marvar,"y")])
  }

  # print(final_var)

  Perf=prediction(varlist = final_var, df=df_p_final, outvar="y", intercept=T, i_number=i_numb, prediction_type = prediction_type, predict_method = predict_method)
  Perf[,c("MV", "IV", "method")] = list(MV, IV, "PLS_HDSI")
  Prop_Performance = Perf
  # print(Prop_Performance)
  return(list(Prop_Performance, final_var))
}

#'@export
fit_function = function(datatype=c("simulated", "real"), param = list(miss = 33, corr = 0.9, cutter = 10, seed =1, data_code=1, test_train_ratio=0.2, varnum = 15, setting= "Correlation", main_var=10, var_effect=c(0.5, -0.5), correlation_var=15, correlation_val=5, high_dim=T, train_sample=500, var = "Mar"), seeder=1, ncomp=3, i_numb=2, k=5, effectsize=13, model_perf="None", coarse_FS = T, Fine_FS =T, summary_ci = 0.95, strict_anova=F, semi_strict=F, varmax=10, kmean_type = "Normal", predict_method = "aridge"){
  # Prepare the dataset
  df_p = data_fit(datatype = datatype, param = param)
  # str(df_p[[1]][,"y"])
  # Perform the model analysis
  # print(c(kmean_type,predict_method))

  Performance = Prop_mod(df_p = df_p, seeder = seeder, ncomp=ncomp, i_numb=i_numb, k=k, effectsize=effectsize, model_perf=model_perf, coarse_FS = coarse_FS, Fine_FS =Fine_FS, summary_ci = summary_ci, strict_anova=strict_anova, semi_strict=semi_strict, varmax=varmax, kmean_type = kmean_type, predict_method = predict_method)
  return(Performance)
}

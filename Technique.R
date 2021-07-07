#' The functionS in this file define the function used for creating the model
#'
#'@export
# Define HDSI_PLS
plsregfunct=function(df, ncomp=3, interactionterm=T, int_term=3, model_perf=c("rsq", "rmse", "bic", NA)){
  if(interactionterm==T){
    if(int_term==3){
      f=as.formula(y ~ .*.*.)
    }else{f=as.formula(y ~ .*.)}
  }else{f=as.formula(y ~ .)}
  Matrix=model.matrix(f,df)[,-1]
  #str(Matrix); str(df)
  #cat(T)
  Varlist=stringr::str_replace_all(colnames(Matrix),"[:_*]", "_")
  #model= plsdepot::plsreg1(Matrix, df$y, comps = ncomp, crosval = F)
  model = tryCatch(plsdepot::plsreg1(Matrix, df$y, comps = ncomp, crosval = F), error = function(e){str(Matrix); str(df)})
  r_sq=sum(model$R2[1:ncomp])
  if(model_perf=="rmse"){MP=1/RMSE(model$y.pred, df$y)}
  else if (model_perf=="bic"){MP = -1*(-2*log(MSE(model$y.pred, df$y)) + ncol(Matrix) * log(nrow(Matrix)))}
  else if(model_perf=="rsq"){MP=r_sq}
  else(MP = model$reg.coefs[-1])

  beta_coef=model$reg.coefs[-1]
  df_out=data.frame(Variable=Varlist, beta=MP, imp = beta_coef, stringsAsFactors = F)
  return(df_out)
}
# mplsregfunct=memoise::memoise(plsregfunct)

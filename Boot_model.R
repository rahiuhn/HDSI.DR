bootridge = function(f, traindf, varlist, outvar){
  Matrix <- model.matrix(f, traindf[,c(varlist, outvar)])

  if(length(varlist) >1){Matrix= Matrix[,-1]}
  Y=as.matrix(traindf[,c(outvar)])

  # Find weight
  cv.ridge <- glmnet::cv.glmnet(Matrix, Y, alpha=0, standardize=F)
  dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[2:(ncol(Matrix)+1), 1]))^0.25 ## Using gamma = 1
  dif_wt[dif_wt[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999

  #Run Ridge
  fit = glmnet::cv.glmnet(Matrix, Y, alpha=0, penalty.factor=dif_wt, standardize=F) # get optimum lambda
  lambda.1se = fit$lambda.1se
  fit = glmnet::glmnet(Matrix, Y,alpha=0, penalty.factor=dif_wt,lambda = lambda.1se, standardize=F)
  Coef_beta = as.matrix(coef(fit,s=lambda.1se))
  Coef = data.frame(Coef_beta)
  Coef[,c("Variable", "Estimate")] = c(rownames(Coef), Coef[,1])
  Coef = Coef
  return(Coef)
}
bootreg = function(f, traindf, varlist, outvar){
  model=lm(f, data = traindf[,c(varlist, outvar)])
  Coef = data.frame(coef(summary(model)), stringsAsFactors = F)
  Coef$Variable = rownames(Coef)
  return(Coef)
}

#'@export
bootprediction = function(varlist=sel_feature_list, df=df_p, outvar="y", intercept=T, i_number=2, predict_method=c("aridge", "reg"), effectsize = 32, seed=1, k=16, boots=NA){
  # str(df)
  # Generate the bootstraps
  p = ncol(df[[1]])-1
  rows = nrow(df[[1]])
  inputdf = df[[1]]

  if(is.na(boots)){boots = p_est(p=p, k=k, rows=rows, interaction_numb=2, effectsize=effectsize)}
  if(boots != 0){
    rowlist = lapply(1:boots, function(x) {
    set.seed(x*1)
    random.samples=sample(rownames(inputdf), nrow(inputdf), replace=T)
    })
  }
  else{rowlist = list(rownames(inputdf))}

  # Generate the interaction based dataset
  if(!is.na(varlist)){varlist=var_organise(varlist)}

  # Create train and test df
  traindf=int_dfcreate(df[[1]],interaction_number=i_number)
  names(traindf)=var_organise(names(traindf))

  testdf=int_dfcreate(df[[2]], interaction_number=i_number)
  names(testdf)=var_organise(names(testdf))

  # Create general formula
  f=as.formula(y~.)

  if(length(varlist)>=1){
    if(!is.na(varlist)){
      # Run the Bootstrap model
      coeflist = lapply(rowlist, function(x){
        traindf = traindf[x,]
        if(predict_method=="aridge"){coefdf = bootridge(f=f, traindf=traindf, varlist=varlist, outvar=outvar)}
        else{coefdf = bootreg(f=f, traindf=traindf, varlist=varlist, outvar=outvar)
        coefdf}
      })

      coefdf = do.call(rbind, coeflist)

      Coef = ddply(coefdf, .(Variable), function(x) mean(x$Estimate))
      names(Coef) = c("Variable", "Estimate")

      cond = -grep("intercept|Intercept", Coef$Variable)
      intercept_val=ifelse(intercept==T, Coef$Estimate[Coef$Variable == "(Intercept)"],0)

      Coef=Coef[cond,c("Variable", "Estimate")]
      names(Coef)=c("Variable", "beta")

      Coef$beta=as.numeric(Coef$beta)
      intercept_val=as.numeric(intercept_val)
      var_effect=lapply(1:nrow(Coef), function(x){
        variable=strsplit(Coef$Variable[x], "_");
        colproduct(type = "column_constant", inputdf=testdf, var_list=variable, constant=Coef$beta[x])
      })
      y_est=apply(do.call(rbind,var_effect),2,function(x) sum(x, na.rm = T))
      y_est=y_est+intercept_val

      ## Estimate the performance
      Perf=tryCatch({predict_metric(actualy = testdf$y, predictedy = y_est)},
                    error=function(e){print("Error"); blank_metric(y=testdf$y)})
    }
    else{Perf = blank_metric(y=testdf$y)}}
  else{Perf = blank_metric(y=testdf$y)
  #data.frame(Corr=0, RMSE = sqrt(MLmetrics::MSE(testdf$y,mean(testdf$y, na.rm=T))), rsquare = 0, stringsAsFactors = F)
  }

  return(Perf)
}





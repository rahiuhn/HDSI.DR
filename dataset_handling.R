#' The function in this file prepare the dataset for the modeling. Currently test dataset size is fixed at 500 samples
#'
#'
#'@varnum: Number of features,p, in the model which is equal to sum of marginal target features and noise features
#'@setting: Should the features have any correlation or not. Two options: "No_Correlation" or "Correlation"
#'@var: Define the type of target features. "Mar" means both marginal and interaction terms are created as target features, "No_Mar" means only interaction terms are created as target features and "No_Var" mean null model with no target features is created.
#'@seed: Control the random seed generator 
#'@main_var: Number of target features in the model
#'@var_effect: Effect size of each target features. The input should be a vector.
#'@correlation_var: Number of marginal features with correlation.
#'@correlation_val: Covariance value for the correlated features. Perfect Positive Correlation is 10. The input should be a vector.
#'@high_dim: Different ways to just generate data. Use this to generate high dimensional data
#'@train_sample: Sample size of training sample

# Generate the artificial dataset
dataset=function(varnum, setting="No_Correlation", var=c("Mar", "No_Mar", "No_Var"), seed=2, main_var=10, var_effect=0.5, correlation_var=15, correlation_val=5, high_dim=T, train_sample=500){

  # Create the Covariance Matrix
  Sigma=matrix(rep(0,varnum), nrow=varnum, ncol=varnum, byrow=F)
  if(high_dim){
    for(i in 1:varnum){
      for(j in 1:varnum){
        if(i==j){Sigma[i,j]=10}
        else if(i<=correlation_var & j<=correlation_var & setting == "Correlation"){Sigma[i,j]=Sigma[j,i]=correlation_val}
        else{Sigma[i,j]=Sigma[j,i]=1}
      }
    }
  }
  else{
    for(i in 1:varnum){Sigma[i,i]=10}
    # Correlation Settings
    if(setting=="Correlation"){
      Sigma[1,2]=3;Sigma[1,3]=3;Sigma[1,4]=6;Sigma[1,5]=6
      Sigma[2,1]=3;Sigma[3,1]=3;Sigma[4,1]=6;Sigma[5,1]=6
      Sigma[2,3]=3;Sigma[2,4]=2;Sigma[2,5]=1
      Sigma[3,2]=3;Sigma[4,2]=2;Sigma[5,2]=1
      Sigma[3,4]=2;Sigma[3,5]=1
      Sigma[4,3]=2;Sigma[5,3]=1
      Sigma[4,5]=1
      Sigma[5,4]=1
    }
  }

  # Create the input dataset
  #print(Sigma)
  set.seed(seed)
  ta=data.frame(MASS::mvrnorm(n = train_sample+500, rep(0, varnum), Sigma/10))
  variablelist=list()
  for(i in 1:varnum){
    variablelist[[i]]=gsub(" ", "",paste("X",i))
    ta[,i]=mosaic::zscore(ta[,i])
  }
  # str(ta[,1:3])
  variablelist=unlist(variablelist)

  # Create the outcome Variable
  intercept=1
  if(high_dim){
    betas = c(rep(var_effect,(2*main_var)-1))
    beta_value = betas[1:((2*main_var)-1)]
    if(var=="Mar"){beta_value=beta_value*rep(1,(2*main_var)-1)}
    else if(var=="No_Mar"){beta_value=beta_value*c(rep(0,2), rep(1,main_var-2), rep(1,main_var-1))}
    else if(var=="No_Var"){beta_value=beta_value*c(rep(0,main_var), rep(0,main_var-1))}
    else if(var=="only_int"){beta_value=beta_value*c(rep(0,main_var), rep(1,main_var-1))}
    else {beta_value = beta_value*c(rep(1,main_var), rep(0,main_var-1))}

    # Generate theinteraction term
    mar_var= paste(c(names(ta)[1:main_var]), collapse = "+")
    int_var= paste(names(ta)[1:(main_var-1)],"*" , names(ta)[2:main_var], collapse = " + ")
    f = as.formula(paste("~",mar_var ," +" , int_var))
    # print(f)
    main_mat = model.matrix(f, ta)
    # print(main_mat[1:5,])
    # Get Outcome
    set.seed(2)
    random_value = rnorm(n=train_sample+500, mean=0, sd=0.25)
    # coef_value = apply(main_mat,1, function(x) {sum(x*c(intercept, beta_value))})

    ta$y = main_mat %*% c(intercept, beta_value) + random_value
    # ta$y  = coef_value + random_value
  }
  else{
    if(var=="Mar"){beta_a=1; beta_b=1}
    else if(var=="No_Mar"){beta_a=0; beta_b=1}
    else{beta_a=0; beta_b=0}
    b1=0.2*beta_a
    b2=0.3*beta_a
    b3=0.4*beta_b
    b4=0.3*beta_b
    set.seed(2)
    ta$y= intercept + (b1*ta$X1) + (b2*ta$X2) + (b3*ta$X3) + (b4*ta$X1*ta$X2) + rnorm(n=train_sample+500, mean=0, sd=0.25)
  }
  set.seed(2)
  index=sample(1:nrow(ta), train_sample, replace = F)
  traindf=ta[index,]
  validationdf=ta[-index,]
  return(list(train=traindf,test=validationdf))
}


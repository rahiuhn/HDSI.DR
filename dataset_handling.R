#' The functionS in this file prepare the dataset for the modeling
#'
#'
#'@export
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

# Process the real dataset
## Remove the variables with missing data
missing_data=function(inputdf, missingpercent){
    # Determine the missing data distribution across the dataset
    Missing_data=sapply(inputdf, function(x) 100*sum(is.na(x))/nrow(inputdf))
    # Remove variables with more than "missingpercent"
    var=attributes(Missing_data)$names[Missing_data>missingpercent]
    inputdf[,var]=NULL
    return(inputdf)
}
## Remove the variables with low variation
lowvar=function(inputdf, minfreq=100/10){
    al= unlist(apply(inputdf,2,function(x) caret::nearZeroVar(x, freqCut = minfreq, saveMetrics = T)[4]))
    drop_var_2=names(al)[which(al==TRUE)]
    drop_var_2 = stringr::str_replace(drop_var_2,".nzv","")
    newdf=inputdf[,setdiff(names(inputdf),drop_var_2)]
    return(newdf)
  }
## Remove the variables with high correlation To keep maximum variables
Correl=function(inputdf, cutoff){
    dfcor=cor(inputdf, use = "pairwise.complete.obs")
    highcor=apply(dfcor,1,function(x) length(which(x>cutoff | x< (-1*cutoff))))
    var=attributes(highcor)$names
    newdf=inputdf[,var]
    while (max(highcor>1)) {
      sort=order(-highcor)
      # drop the variable
      var=attributes(highcor)$names[sort[-1]]
      newdf=inputdf[,var]
      dfcor=cor(newdf, use = "pairwise.complete.obs")
      highcor=apply(dfcor,1,function(x) length(which(x>cutoff | x< (-1*cutoff))))
    }
    return(newdf)
}
## Remove the Outliers
out_rem=function(inputdf){
  Meanprofile=apply(inputdf, 2, function(x) mean(x, na.rm=T))
  sdprofile=apply(inputdf, 2, function(x) sd(x, na.rm=T))
  m_plus_s=Meanprofile+(4*sdprofile)
  m_minus_s=Meanprofile-(4*sdprofile)
  for(x in names(inputdf)){
    len=length(inputdf[which(inputdf[,x]>m_plus_s[x] | inputdf[,x]<m_minus_s[x] ),x])
    if(len>0){
      inputdf=inputdf[-which(inputdf[,x]>m_plus_s[x] | inputdf[,x]<m_minus_s[x]),]
    }
  }
  return(inputdf)
}

# Data Preparation Function
#' @importFrom magrittr %>%
df_process=function(inputdf,outcomedf,missingper=100,minvar=100/10,corr_cut=0.1, outcomevar=c("Unhealthy_Days"), cut_ratio=10,
                    rowsize=20, colsize=4, seeder=1, test_train_ratio=0.2){
  # Remove the independent variables with missing data and near zero variable frequency
  input_miss=missing_data(inputdf=inputdf, missingpercent=missingper)
  #str(input_miss)
  input_nozero=lowvar(inputdf=input_miss, minfreq=minvar)
  # str(input_nozero)
  if(typeof(input_nozero)=="integer"){return(NULL)}else{
    if(ncol(input_nozero)<2){return(NULL)}else{
      # Remove variables with high correlations
      input_nocorrel=Correl(inputdf=input_nozero, cutoff=corr_cut)
      #print(ncol(input_nocorrel))
      input_nocorrel=data.frame(input_nocorrel)
      #str(input_nocorrel)
      #rescale_f=function(x) scales::rescale(x, na.rm=T)
      rescale_f=function(x) {#hist(x);
        y=arm::rescale(x);
        #hist(y);
        return(y)} # Does Z scaling. It is desirable for models relying on Gaussian Distribution
      dataset=data.frame(apply(input_nocorrel,2, rescale_f))
      # str(dataset)
      dataset$y=outcomedf[,outcomevar]
      final_df=dataset[which(!is.na(dataset$y)),]
      # str(final_df)
      comp_set=final_df[complete.cases(final_df),]
      # str(comp_set)
      # Remove the empty columns, low variance columns and dataset with no "y" from final dataset
      comp_set_pre = comp_set %>% janitor::remove_empty("cols")
      comp_set_pre = lowvar(inputdf=comp_set_pre, minfreq=minvar)
      # str(comp_set_pre)
      if(ncol(comp_set_pre)>1){
        # print(1)
        cond_1= names(comp_set_pre)[ncol(comp_set_pre)]=="y"
        cond_2= nrow(comp_set_pre)/ncol(comp_set_pre)>cut_ratio
        cond_3= nrow(comp_set_pre)>=rowsize
        cond_4= ncol(comp_set_pre)>=colsize
        # cat(cond_1,cond_2,cond_3,cond_4)
        if(cond_1 & cond_2 & cond_3 & cond_4){
          #print(2)
          set.seed(seeder)
          index=sample(1:nrow(comp_set_pre),floor(test_train_ratio*nrow(comp_set_pre)))
          # print(index)
          comp_set_train=comp_set_pre[-index,]
          comp_set_test=comp_set_pre[index,]
          return(list(train=comp_set_train,test=comp_set_test))
        }else{return(NULL)}
      }else{return(NULL)}
    }}
}
mdf_process=memoise::memoise(df_process)

#'@export
data_prep=function(hyperparameter = data.frame(miss = 33, corr = 0.9, cutter = 10, seed =1, data_code=1), seed=1, test_train_ratio=0.2, file_source=c("D:/in_d", "D:/out_d")){
  data_prep=lapply(1:(nrow(hyperparameter)), function(x) {
    #print(c("Para",x))
    miss=hyperparameter[x,1]
    corr=hyperparameter[x,2]
    cutter=hyperparameter[x,3]
    if(is.na(seed)){seed=hyperparameter[x,4]}

    # Get the dataset
    {
      data_code=hyperparameter[x,5]
      in_data=  paste0(file_source[1],data_code,".csv")
      out_data= paste0(file_source[2],data_code,".csv")
      Outdata=read.csv(out_data)
      Inputdata=read.csv(in_data)
      Outdata=data.frame(out=Outdata[,2])
      Inputdata=Inputdata[,-1]
      # str(Inputdata)
      # print(ncol(Inputdata))
    }
    dp=mdf_process(missingper=miss, corr_cut=corr,seeder=seed, test_train_ratio=test_train_ratio,
                   inputdf=Inputdata,outcomedf=Outdata,minvar=100/10, outcomevar=c("out"), cut_ratio=cutter)
    dp
  })
  return(data_prep)
}

#'@export
data_fit =  function(datatype = c("simulated", "real"), param = list(miss = 33, corr = 0.9, cutter = 10, seed =1, data_code=1, test_train_ratio=0.2, varnum = 15, setting= "Correlation", main_var=10, var_effect=c(0.5, -0.5), correlation_var=15, correlation_val=5, high_dim=T, train_sample=500, var = "Mar")){
  if(datatype == "real"){
    # print("Enter")
    # print(param)
    hyperpara=data.frame(miss = param$miss, corr = param$corr, cutter = param$cutter, seed= param$seed, data_code = param$data_code)
    df_p=data_prep(hyperpara[1,], seed = param$seed, test_train_ratio = param$test_train_ratio)
    df_p=df_p[[1]]
  }
  else{
    df_p=dataset(varnum = param$varnum, setting= param$setting, var= param$var, main_var=param$main_var, var_effect=param$var_effect, correlation_var=param$correlation_var, correlation_val=param$correlation_val, high_dim=param$high_dim, train_sample=param$train_sample, seed = param$seed)
  }
  {
    # Generate Dataset
    realnames=names(df_p[[1]])
    # print(realnames)
    fakename=paste0("X",seq(1:(length(realnames)-1)))
    names(df_p[[1]])[1:(length(realnames)-1)]=fakename
    names(df_p[[2]])[1:(length(realnames)-1)]=fakename
  }
  return(df_p)
}

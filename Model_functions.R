#' The functionS in this file define the various functions used for running the model
#'


# Generate the bootstraps
p_est=function(p,k=NA,rows=NA, interaction_numb=2, effectsize=13){
  if(is.na(k)==T){
    optimal_k=function(x, r=rows){
      int=choose(x,interaction_numb)
      value=abs((int+x)-r)/r
      return(value)
    }
    k=optimize(optimal_k, interval = c(interaction_numb,p))$minimum
    k=floor(k)
  }
  denominator=choose(p,interaction_numb)
  numerator=choose(k,interaction_numb) # where, 2 is the order of interaction eg of 2 order interaction is X1_X2
  init_prob=numerator/denominator
  #print(init_prob)
  # Minimum bootstraps which will give minimum 13 occurrences (consider large effect) of an interaction variable with 99% confidence
  if(effectsize=="large"){prefer_occurence = 13}
  else if (effectsize == "medium"){prefer_occurence = 32}
  else if(effectsize == "small"){prefer_occurence = 200}
  else {prefer_occurence = effectsize}

  min_boot= ceiling(prefer_occurence/init_prob) # where 32 is minimum number of occurence desired for an interaction variable during bootstrapping.
  max_boot= 8*min_boot
  #print(c(min_boot, max_boot))
  f_opt=function(x,y=init_prob,max_x=max_boot,prefer=prefer_occurence){
    Val=qbinom(0.05, floor(x), y)
    value=(prefer-Val)^2+(x/max_x)
    return(value)
  }
  bootvalue=optimize(f_opt,interval = c(min_boot,max_boot))
  return(ceiling(bootvalue$minimum))
}
#'@export
bootsample=function(p,k,rows, plist, inputdf, intlist=FALSE, type="stratified", seed_multiplier=1,interaction_numb=2, effectsize=13){
  if(is.na(k)==T){
    optimal_k=function(x, r=rows){
      int=choose(x,2)
      value=abs((int+x)-r)/r
      return(value)
    }
    k=optimize(optimal_k,interval = c(2,p))$minimum
    k=floor(k)
  }
  #print(k)
  boots=p_est(p=p,k=k, rows=rows, interaction_numb=interaction_numb, effectsize=effectsize)

  ## Create Variable list for each bootstrap
  features=plist
  samples=1:rows
  res=lapply(1:boots, function(x) {
    set.seed(x*seed_multiplier)

    ### Create variable combinations
    random.features <- sample(features, k, replace = FALSE) # replaced samples with features/3
    feature=naturalsort::naturalsort(random.features)
    if(intlist==T){
      int_mat=lapply(1:(length(feature)-1), function(y){
        mat1=as.matrix(feature[1:(length(feature)-y)])
        mat2=as.matrix(feature[(1+y):(length(feature))])
        matname=paste(mat1,mat2,sep="_")
        matname
      })
      int_list=unlist(int_mat)
      varlist= union(random.features, int_list)
    }else{varlist=feature}

    ## Create Sample list for each bootstrap
    if(type=="stratified"){
      SURV_1=inputdf[which(inputdf$status==1),]
      S_1=sample(rownames(inputdf[which(inputdf$status==1),]), nrow(SURV_1), replace=T)
      SURV_0=inputdf[which(inputdf$status==0),]
      S_0=sample(rownames(inputdf[which(inputdf$status==0),]), nrow(SURV_0), replace=T)
      random.samples <- as.numeric(c(S_1,S_0))
    }else{random.samples=sample(rownames(inputdf), nrow(inputdf), replace=T)}

    list(varlist, random.samples)
  })
  return(res)
}
mbootsample=memoise::memoise(bootsample)

# Data processing functions
## Organise variables in two-way interactions
var_organise=function(inputlist, symbol="_"){
  org_var=lapply(inputlist, function(x) {
    a=strsplit(x, symbol);
    ifelse(length(a[[1]])>1, paste0(naturalsort::naturalsort(a[[1]]),collapse = symbol), a[[1]])})
  return(unlist(org_var))
}
var_combination=function(listofvar=list(), r=NULL){
  if(is.null(r)==TRUE){r=length(listofvar)-1}
  c=lapply(r, function(x) combn(listofvar, x, simplify = F))
  comb=c[[1]]
  if(length(r)>1){
    for(i in 2:length(r)){
      comb=append(comb,c[[i]])
    }
  }
  return(comb)
}
## Create the interaction list for simsurv df and model
interactionlist=function(listofvar=list(), interactionlength=2, comb=NULL, symbol="_"){
  if(is.null(comb)==TRUE){comb=var_combination(listofvar=listofvar, r=interactionlength)}
  i_list=lapply(comb, function(x) stringi::stri_paste(x, collapse = symbol))
  i_list=unlist(i_list)
  i_list=var_organise(i_list)
  return(i_list)
}
## Generate column products
colproduct=function(type=c("column", "column_constant", "column_list", "list"), ...){
  arglist=list(...)
  inputdf=arglist[[1]] #inputdf represents the df which contain the variable data,
  #print(names(inputdf))
  var_list=unlist(arglist[[2]]) # It provides list of all columns that need to be mulitplied together
  #print(var_list)
  if(type == "column"){
    inputdf$newprod=inputdf[,var_list[1]]
    if(length(var_list)>1){
      for(i in 2:length(var_list)){
        inputdf$newprod=inputdf$newprod*inputdf[,var_list[i]]
      }
    }
  }
  else if(type == "column_constant"){
    constant=unlist(arglist[[3]]) # It provides list of all constant that need to be mulitplied with the columns
    constantprod=prod(constant)
    inputdf$newprod=inputdf[,var_list[1]]*constantprod
    if(length(var_list)>1){
      for(i in 2:length(var_list)){
        inputdf$newprod=inputdf$newprod*inputdf[,var_list[i]]
      }
    }
  }
  return(inputdf$newprod)
}

# Summarize the bootstrap results
f_summary_rsquare=function(x, cint=0.95, min_max=c("ci", "quartile", "min")){
  #cat(x)
  quant=(1-cint)/2
  x=x[!is.na(x)]
  M=mean(x, na.rm=T); stdev=sd(x, na.rm=T); Med=median(x, na.rm = T)
  #coef_var=cv_versatile(x, na.rm = T, digits = 4)$statistics ; cqv=cqv_versatile(x, na.rm = T, digits = 4)$statistics
  quan=quantile(x, c(quant,1-quant)); lqi=quan[1]; uqi=quan[2]
  #gm=geometric.mean(sqrt(x^2), na.rm = T)
  #quan=MedianCI(x); lqi=quan[2]; uqi=quan[3]
  # Summarize the results based on the criteria: Confidene Interval, Quartile Interval or Minimum and Maximum
  if(min_max=="ci"){
    conf=unlist(list(attributes(suppressMessages(miscset:::confint.numeric(x, na.rm = T, level = cint)))));
    lci=conf[2]; uci=conf[3]; #cat(lci[[1]])
    out=c(mean=M, sd=stdev, median=Med, lowerci=lci, upperci=uci, lowerqi=lqi, upperqi=uqi, gm=0, coef_var=0, cqv=0)
  }else if(min_max=="quartile"){
    out=c(mean=M, sd=stdev, median=Med, lowerci=lqi, upperci=uqi, lowerqi=lqi, upperqi=uqi, gm=0, coef_var=0, cqv=0)
  }else{
    minimum=min(x, na.rm = T); maximum=max(x, na.rm = T)
    out=c(mean=M, sd=stdev, median=Med, lowerci=minimum, upperci=maximum, lowerqi=lqi, upperqi=uqi, gm=0, coef_var=0, cqv=0)
  }
  #out=c(mean=M, sd=stdev, median=Med, lowerci=lci, upperci=uci, lowerqi=lqi, upperqi=uqi, gm=gm)
  blank=c(mean=0, sd=0, median=0, lowerci=0, upperci=0, lowerqi=0, upperqi=0, gm=0, coef_var=0, cqv=0)
  if(any(is.na(out))){out<-blank}
  return(out)
}
f_summary_beta=function(x, cint=0.95){
  #cat(x)
  quant=(1-cint)/2
  x=x[!is.na(x)]
  M=mean(x, na.rm=T);  Med=median(x, na.rm = T)
  stdev=sd(x, na.rm=T);
  # conf=unlist(list(attributes(suppressMessages(miscset:::confint.numeric(x, na.rm = T, level = cint)))));
  # lci=conf[2]; uci=conf[3];
  quan=quantile(x, c(quant,1-quant)); lqi=quan[1]; uqi=quan[2]

  out=c(mean=M, sd=stdev, median=Med, lowerci=lqi, upperci=uqi, lowerqi=lqi, upperqi=uqi, gm=0, coef_var=0, cqv=0)
  #print(out)
  blank=c(mean=0, sd=0, median=0, lowerci=0, upperci=0, lowerqi=0, upperqi=0, gm=0, coef_var=0, cqv=0)
  if(any(is.na(out))){out<-blank}
  return(out)
}
f_summary_dispersion=function(x, cint=0.95){
    #cat(x)
    quant=(1-cint)/2
    x=x[!is.na(x)]
    M=mean(x, na.rm=T); stdev=sd(x, na.rm=T); Med=median(x, na.rm = T)
    coef_var=cv_versatile(x, na.rm = T, digits = 4)$statistics ; cqv=cqv_versatile(x, na.rm = T, digits = 4)$statistics

    conf=unlist(list(attributes(suppressMessages(miscset:::confint.numeric(x, na.rm = T, level = cint)))));
    lci=conf[2]; uci=conf[3];

    quan=quantile(x, c(quant,1-quant)); lqi=quan[1]; uqi=quan[2]
    gm=geometric.mean(sqrt(x^2), na.rm = T)

    out=c(mean=M, sd=stdev, median=Med, lowerci=lci, upperci=uci, lowerqi=lqi, upperqi=uqi, gm=gm, coef_var=coef_var, cqv=cqv)
    blank=c(mean=0, sd=0, median=0, lowerci=0, upperci=0, lowerqi=0, upperqi=0, gm=0, coef_var=0, cqv=0)
    if(any(is.na(out))){out<-blank}
    return(out)
}

# Perform Coarse feature Selection
variable_selection=function(df=res_summary){
  df$Selection=apply(df[,-1],1, function(x) {
      range= spatstat.utils::inside.range(0, r=c(min(x[4:7]), max(x[4:7])));
      ifelse(range==TRUE,0,x[1])})

  ## Extract selected features
  cond=df$Variable!="(Intercept)" & df$Selection!=0
  sel_feature=df[cond,]
  sel_list=unlist(stringr::str_split(sel_feature$Variable, "_"))
  sel_feature_list=union(sel_list, sel_feature$Variable)
  return(sel_feature_list)
}
mvariable_selection = memoise::memoise(variable_selection)

# Perform Fine Feature Selection (Clustering based)
{
  euclidean_dist=function(x,y){
    dist=sqrt(sum(c((x - y)^2)))
    return(dist)
  }
  constrain_kmean=function(dataset, steps, centers){
    #print(c("steps", steps))
    results=sapply(1:steps, function(x){
      #print(x)
      sort_list=dataset[order(dataset)]
      norm_sort_list=sort_list/min(sort_list)

      # Get the Main cluster points and distance
      main_cluster_center=norm_sort_list[1]
      main_cluster=norm_sort_list[1:x]
      main_cluster_dist=euclidean_dist(main_cluster_center, main_cluster)

      # Get the other cluster points and distance
      other_points=norm_sort_list[(x+1):length(norm_sort_list)]

      set.seed(1)
      other_cluster= kmeans(other_points, centers, nstart = 10)
      other_cluster_cluster = other_cluster$cluster #cutree(other_cluster, centers)
      other_cluster_center=as.numeric(other_cluster$centers)
      closest_center_index= which(other_cluster_center == min(other_cluster_center, na.rm = T))
      other_cluster_points_index=which(other_cluster_cluster == closest_center_index)

      # Find cluster closest to the main cluster
      close_cluster= other_points[other_cluster_points_index] #other_cluster_points[lapply(other_cluster_points,length)>0][[1]][[1]]
      #print(close_cluster)
      close_cluster_center= other_cluster_center[closest_center_index]#c(other_cluster_points[lapply(other_cluster_points,length)>0][[1]][[2]])
      #print(close_cluster_center)
      kmean_dist = euclidean_dist(close_cluster_center, close_cluster)

      total= main_cluster_dist + kmean_dist
      total})
    return(results)
  }
  get_mode=function(datalist){
    ux <- unique(datalist)
    tab <- tabulate(match(datalist, ux))
    mode_val=ux[tab == max(tab)]
    #print(mode_val)
    if(any(mode_val %in% 1)){
      ux_no1=setdiff(ux,1)
      tab <- tabulate(match(datalist, ux_no1))
      mode_val_no1 = ux_no1[tab == max(tab)]
      #print(mode_val_no1)
      mode_val=union(mode_val, mode_val_no1)
    }
    return(mode_val)
  }

#'@import plyr
  cv_filter=function(pre_filter_list=NA, variable_summary_df=res_summary, totalsteps=1/p, clusters=floor(3+(k/100)), varmax=4,
                     strict_anova=F, inputdf=df, res_df=res_df, cint= 0.95, semi_strict=T, kmean_type = c("Normal", "Constrained")){

    # print(T)
    if(is.na(variable_summary_df)){
      variable_summary_df = plyr::ddply(res_df[,c("Variable","imp")], .(Variable), function(x)
      {a=f_summary_beta(x$imp, cint = 0.95);unlist(a)})
    }
    if(!is.na(pre_filter_list)){
      res_summary_sel_fin <-variable_summary_df[variable_summary_df$Variable %in% pre_filter_list, ]
    }
    else{res_summary_sel_fin <-variable_summary_df}

    # Find Coefficient of Variation
    res_summary_sel_fin$coef_var = res_summary_sel_fin$sd / res_summary_sel_fin$mean
    res_summary_sel_fin$coef_var.est=abs(res_summary_sel_fin$coef_var) #abs(res_summary_sel_fin$coef_var.est)

    ## Remove Outliers
    outliers = boxplot(abs(res_summary_sel_fin$coef_var.est), plot = F)$out
    median_coef_var = median(res_summary_sel_fin$coef_var.est)
    outliers = outliers[which(outliers > median_coef_var*3)]
    # print(res_summary_sel_fin$Variable)
    #print(outliers)
    if(length(outliers) == 0){outlier_free=abs(res_summary_sel_fin$coef_var.est)}
    else{
      outlier_free = abs(res_summary_sel_fin$coef_var.est)[-which(abs(res_summary_sel_fin$coef_var.est) %in% outliers)]
      #print(outlier_free)
    }


    ## Do constrained clustering to select the optimal number of features
    #print(c("Clusters: ", clusters))
    if(kmean_type == "Constrained"){
      if(length(outlier_free)>10){
        steps = ifelse(length(outlier_free)>50, max(45,floor(length(outlier_free)/totalsteps)), length(outlier_free)-5)
        perf=lapply(1:min(clusters, length(outlier_free)-3), function(x) constrain_kmean(dataset = outlier_free, steps = steps, centers = x))
        #print(perf)
        min_list = sapply(perf, function(x) {
          get_min = min(x)
          min_index=which(x==get_min)
        })
        Potential_Variable_numb = get_mode(min_list)
        Potential_Variable_numb=Potential_Variable_numb[!is.na(Potential_Variable_numb)]
        #print(Potential_Variable_numb)
        var_potential=max(Potential_Variable_numb)*max(1,varmax)
      }
      else{var_potential=length(outlier_free)}
    }
    else{
      set.seed(1)
      outlier_free = outlier_free[order(outlier_free)]
      if(length(outlier_free)>2){
        km = kmeans(outlier_free, centers = 2, nstart = 10)
        min_clus_center = which(km$centers == min(km$centers))
        clus_points = which(km$cluster == min_clus_center)
        var_potential = length(clus_points) + max(varmax, 1-length(clus_points))
      }
      else(var_potential = length(outlier_free))
    }

    if(var_potential == 0){return(NA)}
    else{
      ## Perform selection
      res_summary_sel_fin = res_summary_sel_fin[order(res_summary_sel_fin$coef_var.est), c(1,10)]
      res_summary_sel_fin = res_summary_sel_fin[1:min(var_potential,nrow(res_summary_sel_fin), na.rm=T),]

      ### Remove the interaction terms
      marginals=strsplit(res_summary_sel_fin$Variable, "_")
      marginals = unique(unlist(marginals))
      # print(marginals)

      test_variable= stringr::str_replace_all(res_summary_sel_fin$Variable, "_", "*")

      ### Perform Stepwise selection
      {
        prep_data=inputdf[,c("y", marginals)]
        f=as.formula(paste("y~", paste0(test_variable, collapse ="+")))
        prep_mat=model.matrix(f,prep_data)[,-1, drop=FALSE]
        Y=as.matrix(prep_data[,"y"])

        colnames(prep_mat) = stringr::str_replace_all(colnames(prep_mat),"[:_*]", "_")
        prep_data[,colnames(prep_mat)] = list()
        prep_data[,colnames(prep_mat)] = prep_mat
        # str(prep_data)
        M = lm(y~1, data = prep_data)
        M1_f = as.formula(paste("y~", paste0(names(prep_data), collapse ="+")))
        model = step(M, direction = "both", k=log(nrow(prep_data)),trace=F, scope=M1_f)
        model_summary = summary(model)
        # str(model_summary)
      }

      ### Strict or Semi-strict criteria
      if(strict_anova | semi_strict){
        ### Do Final variable selection
        anova_df = model$anova
        # print(anova_df)
        anova_df = anova_df[,c(1,2,6)]
        remove_sign = do.call(rbind,strsplit(anova_df$Step, " "))
        # print(remove_sign)
        # cat(length(remove_sign))
        # print(is.empty.model(model))
        if(length(remove_sign)==0){return(NA)}
        anova_df[-1,1] = remove_sign[,2]
        anova_df[1,2] = 0
        #anova_df=anova_df[anova_df$Variable %in% names(model$coefficients),]

        clean_anova_df = anova_df[which(anova_df$Step %in% names(model$coefficients)[-1]),]
        var_freq=table(clean_anova_df$Step)
        get_duplicate_variable=names(var_freq)[var_freq>1]
        #print(get_duplicate_variable)
        get_duplicate_variable = get_duplicate_variable[!is.na(get_duplicate_variable)]
        removaldf=lapply(get_duplicate_variable, function(x)
        {min_value=min(clean_anova_df$AIC[clean_anova_df$Step %in% x], na.rm = T);
        #print(min_value)
        df <- clean_anova_df[clean_anova_df$AIC >min_value & clean_anova_df$Step %in% x,]})

        if(length(removaldf)>0){
          #print("Enter")
          removaldf=do.call(rbind, removaldf)
          #print(removaldf)
          clean_anova_df = clean_anova_df %>% dplyr::anti_join(removaldf)
        }
        clean_anova_df = clean_anova_df[clean_anova_df$Df <0, ]
        #print(clean_anova_df)
        clean_anova_df = rbind(anova_df[1,], clean_anova_df)
        #print(clean_anova_df)

        anova_value=clean_anova_df$AIC
        #print(anova_value)


        anova_sel=sapply(2: length(anova_value),function(x) {
          if(strict_anova){
            if(anova_value[x-1]<0){Old_var=1.01*anova_value[x-1]}
            else{Old_var=0.99*anova_value[x-1]}
          }else{Old_var=anova_value[x-1]}
          current_var=anova_value[x]
          #print(current_var)
          if(current_var<Old_var){TRUE}else{FALSE}})
        sel_feature_list=names(model$coefficients)[-1]
        #print(sel_feature_list)

        # Feature List
        semi_strict=semi_strict
        if (semi_strict){
          anova_diff=anova_value[-length(anova_value)]-anova_value[-1]
          norm_anova_diff = abs(anova_diff)/sum(abs(anova_diff))
          var_select_norm_anova_diff = norm_anova_diff >= 0.01

          cum_anova_diff = anova_value[1]-anova_value[-1]
          norm_cum_anova_diff = abs(cum_anova_diff)/sum(abs(anova_diff))
          var_select_norm_cum_anova_diff = norm_cum_anova_diff <= 0.99
          criteria_list=list(var_select_norm_anova_diff, var_select_norm_cum_anova_diff, anova_sel)
          #print(criteria_list)
        }
        else{criteria_list=list(anova_sel)}


        feature_list = lapply(criteria_list, function (x) sel_feature_list[x])
        feature_list = unlist(feature_list)
        sel_feature_list = unique(feature_list)
      }
      else{
        sel_feature_list = names(coef(model))[-1]
      }
    }
    return(sel_feature_list)
  }
}


# Performance Evaluation
rsquared=function(actualy, predictedy){
  mean_actual=mean(actualy, na.remove=T)
  rss=sum((actualy-predictedy)^2)
  ess=sum((mean_actual-predictedy)^2)
  tss=ess+rss
  rsquare=1-(rss/tss)
  output=c(rsquare,rss,ess)
  return(rsquare)
}
predict_metric=function(actualy=validationdf[,outname[x]], predictedy=y_estimate){
  Corr= cor(actualy,predictedy, method = "pearson")
  RMSE= sqrt(MLmetrics::MSE(actualy,predictedy))
  MSE = MLmetrics::MSE(actualy,predictedy)
  MAE = MLmetrics::MAE(actualy,predictedy)
  MAPE = MLmetrics::MAPE(y_true = actualy, y_pred = predictedy)
  rsquare=rsquared(actualy = actualy, predictedy = predictedy)
  res=data.frame(Corr=Corr, RMSE=RMSE, rsquare=rsquare, MSE=MSE, MAE=MAE, MAPE=MAPE, stringsAsFactors = F)
  return(res)
}
blank_metric = function(y){
  Corr= 0
  RMSE= sqrt(MLmetrics::MSE(y,mean(y, na.rm=T)))
  MSE = MLmetrics::MSE(y,mean(y, na.rm=T))
  MAE = MLmetrics::MAE(y,mean(y, na.rm=T))
  MAPE = MLmetrics::MAPE(y,mean(y, na.rm=T))
  rsquare=0
  res=data.frame(Corr=Corr, RMSE=RMSE, rsquare=rsquare, MSE=MSE, MAE=MAE, MAPE=MAPE, stringsAsFactors = F)
  return(res)
}

int_dfcreate=function(df, outvar="y", interaction_number=3){
  if(ncol(df)>2){
    ## Generate the interaction list
    if(interaction_number==3){int_f=as.formula(y~.*.*.)}else{int_f=as.formula(y~.*.)}
    Matrix=model.matrix(int_f, df)[, -1]
    colnames(Matrix)= stringr::str_replace_all(colnames(Matrix),"[:_*]", "_")
    var_list=setdiff(names(df),outvar)
    df[,var_list]=list()
    df[,colnames(Matrix)]=Matrix}
  return(df)
}
prediction = function(varlist=sel_feature_list, df=df_p, outvar="y", res=NA, intercept=T, i_number=2, predict_method=c("aridge", "reg"), prediction_type = c("rmse", "gmse", "bic")){

  if(!is.na(varlist)){varlist=var_organise(varlist)}
  # print(predict_method)
  # Create train and test df
  traindf=int_dfcreate(df[[1]],interaction_number=i_number)
  names(traindf)=var_organise(names(traindf))
  #str(traindf[,c(varlist, outvar)])
  testdf=int_dfcreate(df[[2]], interaction_number=i_number)
  names(testdf)=var_organise(names(testdf))

  f=as.formula(y~.)

  if(length(varlist)>=1){
    if(!is.na(varlist)){

      if(predict_method=="aridge"){
        Matrix <- model.matrix(f, traindf[,c(varlist, outvar)])
        # str(Matrix)
        if(length(varlist) >1){Matrix= Matrix[,-1]}
        Y=as.matrix(traindf[,c(outvar)])
        # Find weight
        cv.ridge <- glmnet::cv.glmnet(Matrix, Y, alpha=0, standardize=F)
        dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[2:(ncol(Matrix)+1), 1]))^0.25 ## Using gamma = 1
        dif_wt[dif_wt[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
        # print(dif_wt)
        #Run Ridge
        fit = glmnet::cv.glmnet(Matrix, Y, alpha=0, penalty.factor=dif_wt, standardize=F) # get optimum lambda
        lambda.1se = fit$lambda.1se
        fit = glmnet::glmnet(Matrix, Y,alpha=0, penalty.factor=dif_wt,lambda = lambda.1se, standardize=F)
        Coef_beta = as.matrix(coef(fit,s=lambda.1se))
        Coef = data.frame(Coef_beta)
        Coef[,c("Variable", "Estimate")] = c(rownames(Coef), Coef[,1])
        # print(Coef[Coef$Variable == "X15" | Coef$Variable == "X25"| Coef$Variable == "X15_25" | Coef$Variable == "X25_X15", ])
        Coef = Coef[,-1]
      }
      else{
        model=lm(f, data = traindf[,c(varlist, outvar)])
        # print(model)
        # newdata = data.frame(testdf[,varlist])
        # names(newdata) = varlist
        # predict_val = predict(model, newdata)
        Coef = data.frame(coef(summary(model)), stringsAsFactors = F)
        Coef$Variable = rownames(Coef)
      }

      cond = -grep("intercept|Intercept", Coef$Variable) #Coef$Variable != "(Intercept)" #grep("intercept", Coef$Variable)
      # print(cond)
      intercept_val=ifelse(intercept==T, Coef$Estimate[Coef$Variable == "(Intercept)"],0)
      # print(intercept_val)
      Coef=Coef[cond,c("Variable", "Estimate")]
      names(Coef)=c("Variable", "beta")
      # print(Coef)
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
                     error=function(e){print("Error");
                       blank_metric(y=testdf$y)
                       # data.frame(Corr=0, RMSE = sqrt(MLmetrics::MSE(testdf$y,mean(testdf$y, na.rm=T))), rsquare=0, stringsAsFactors = F)
                       })
    }
    else{
      Perf = blank_metric(y=testdf$y)
        #data.frame(Corr=0, RMSE = sqrt(MLmetrics::MSE(testdf$y,mean(testdf$y, na.rm=T))), rsquare = 0, stringsAsFactors = F)
      }
  }
  else{Perf = blank_metric(y=testdf$y)}
  # print(Perf)
  # Get the saturated model value
  if(prediction_type == "gmse"){
    actual_y_sat = testdf[,outvar]
    # Prepare Saturated Model
    if(predict_method=="aridge"){
      Matrix <- model.matrix(f, traindf)[, -1]
      Y=as.matrix(traindf[,c(outvar)])
      # Find weight
      cv.ridge <- glmnet::cv.glmnet(Matrix, Y, alpha=0, standardize=F)
      dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[2:(ncol(Matrix)+1), 1]))^0.25 ## Using gamma = 1
      dif_wt[dif_wt[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
      #Run Ridge
      fit = glmnet::cv.glmnet(Matrix, Y, alpha=0, penalty.factor=dif_wt, standardize=F) # get optimum lambda
      lambda.1se = fit$lambda.1se
      fit = glmnet::glmnet(Matrix, Y,alpha=0, penalty.factor=dif_wt,lambda = lambda.1se, standardize=F)
      newdata = model.matrix(f, testdf)[,-1]
      predict_val_sat = glmnet::predict.glmnet(fit, newx = newdata, s=lambda.1se, type = "response")
      #print(predict_val)
    }
    else{
      model=lm(f, data = traindf)
      newdata = testdf
      predict_val_sat = predict(model, newdata)
    }
    # Estimate test output
    sat_mse = MLmetrics::MSE(actual_y_sat,predict_val_sat)
    # Get performance
    var_num = ifelse(Perf$Corr != 0, length(varlist), 0)
    loss = Perf$RMSE
    SSE_m = (loss^2)*length(actual_y_sat)
    GMC_EST = ((SSE_m/sat_mse) + (var_num+1)*log(length(actual_y_sat)))
    Perf$RMSE = GMC_EST
  }
  if(prediction_type == "bic"){

    actual_y_sat = testdf[,outvar]
    # Prepare Saturated Model
    if(predict_method=="aridge"){
      Matrix <- model.matrix(f, traindf)[, -1]
      Y=as.matrix(traindf[,c(outvar)])

      # Find weight
      cv.ridge <- glmnet::cv.glmnet(Matrix, Y, alpha=0, standardize=F)
      dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[2:(ncol(Matrix)+1), 1]))^0.25 ## Using gamma = 1
      dif_wt[dif_wt[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999

      #Run Ridge
      fit = glmnet::cv.glmnet(Matrix, Y, alpha=0, penalty.factor=dif_wt, standardize=F) # get optimum lambda
      lambda.1se = fit$lambda.1se
      fit = glmnet::glmnet(Matrix, Y,alpha=0, penalty.factor=dif_wt,lambda = lambda.1se, standardize=F)
      newdata = model.matrix(f, testdf)[,-1]
      predict_val_sat = glmnet::predict.glmnet(fit, newx = newdata, s=lambda.1se, type = "response")
      #print(predict_val)
    }
    else{
      model=lm(f, data = traindf)
      newdata = testdf
      predict_val_sat = predict(model, newdata)
    }
    # Estimate test output
    sat_mse = MLmetrics::MSE(actual_y_sat,predict_val_sat)
    # Get performance
    var_num = ifelse(Perf$Corr != 0, length(varlist), 0)
    loss = Perf$RMSE
    SSE_m = (loss^2)*length(actual_y_sat)
    # print(T)
    # rmse_part = log(SSE_m/length(actual_y_sat))
    rmse_part = log(SSE_m/(length(actual_y_sat)-(var_num+1)))
    # feature_penalty = (var_num+1)*log(length(actual_y_sat))/length(actual_y_sat)
    # feature_penalty = (var_num+1)*2/length(actual_y_sat)
    feature_penalty=0
    # cat(rmse_part, " ")
    BIC_EST = (rmse_part + feature_penalty)*length(actual_y_sat)
    Perf$RMSE = BIC_EST
  }
  return(Perf)
}

# mprediction=memoise::memoise(prediction)

# Output processing functions
var_organise = function(inputlist, symbol="_"){
  org_var=lapply(inputlist, function(x) {
    a=strsplit(x, symbol);
    ifelse(length(a[[1]])>1, paste0(naturalsort::naturalsort(a[[1]]),collapse = symbol), a[[1]])})
  return(unlist(org_var))
}
#'@export
rmse_extractor = function(l=res, lindex=1){
  rmselist = lapply(l, function(x)x[[lindex]])
  rmsedf= data.frame(do.call(rbind, rmselist), stringsAsFactors = F)
  return(rmsedf)
}
#'@export
feature_extractor = function(l=res, lindex=2, iter=10){
  featurelist = lapply(l, function(x)x[[lindex]])
  featurevec= unlist(featurelist)
  # Organise the variable
  featurevec = var_organise(featurevec)
  feature_freq= table(featurevec)/iter
  return(feature_freq)
}

#'@export
feature_listcreator = function(l=res, lindex=2, iter=10){
  featurelist = lapply(l, function(x)x[[lindex]])
  return(featurelist)
}

#'@export
dftolist = function(df,dfpath=NA, variablename = "Variable", frequencyname = "Freq"){
  if(!is.na(dfpath)){df = read_excel(dfpath)}
  if(is.na(variablename) | is.na(frequencyname)){
    varname = df[,1]
    freq_value = df[,2]
  }
  else{
    varname = df[,variablename]
    freq_value = df[,frequencyname]
  }
  create_l = list()
  create_l[varname] <- lapply(1:nrow(df), function(x) freq_value[x])
  freq_list = unlist(create_l)
  return(freq_list)
}

freq_table = function(tablelist, varlist, target = c("X1", "X2", "X3", "X1_X2", "X5_X9"), freq_name = "RHDSI"){
  df= data.frame(var= names(tablelist), freq=as.numeric(tablelist), stringsAsFactors = F)
  names(df)[2] = freq_name
  other_var = setdiff(varlist, names(tablelist))
  len_othervar = length(other_var)
  df[(nrow(df)+1):(nrow(df)+len_othervar),] = c(other_var, rep(0,len_othervar))
  if(is.na(target)){
    target = unique(unlist(strsplit(varlist,"_")))
    df$vartype = ifelse(df$var %in% target, "MF", "IF")
  }
  else{
    df$vartype = ifelse(df$var %in% target, "T", "N")
  }

  return(df)
}
tn_table = function(varnum=15, vartype= c("int", "mar"), frequencylist, freq_name, target =  c("X1")){
  # Create list of variables
  {
    marvar = sapply(1:varnum, function(i) gsub(" ", "",paste("X",i)))
    intlist = var_combination(listofvar = marvar, r=2)
    intlist= sapply(intlist, function(x) paste(x, collapse = "_"))
    org_intlist = var_organise(inputlist = intlist)
    varlist = union(marvar, org_intlist)
  }
  # Create Frequency Table
  listsize = length(frequencylist)
  ## Check the list size and list name uniformity
  if(listsize != length(freq_name)) {print("Check the frequencylist and freq_name")}
  dflist= lapply(1:listsize, function(x) freq_table(tablelist = frequencylist[[x]], varlist = varlist, freq_name = freq_name[x], target = target))
  df = Reduce(function(x, y) merge(x, y, all=TRUE), dflist)
  return(df)
}
#'@export
tn_graph = function(varnum=15, vartype= c("int", "mar"), frequencylist, freq_name, target =  c("X1"), box = F){
  # Prepare df
  df = tn_table(varnum=varnum, vartype= vartype, frequencylist = frequencylist, freq_name = freq_name, target = target)
  longdf = reshape2::melt(df, id.vars = c("var","vartype"))
  # Prepare graph
  if(box){
    p <- ggplot2::ggplot(longdf, ggplot2::aes(variable, as.numeric(value))) +
      ggplot2::geom_point(ggplot2::aes(shape = vartype, colour = vartype, size = vartype)) +
      ggplot2::geom_boxplot(ggplot2::aes(colour=vartype))+
      ggplot2::labs(title="Target and Noise variable distribution",x= "Variable Type", y = "Frequency") +
      ggplot2::scale_shape_manual(values=c(2, 16))+
      ggplot2::scale_color_manual(values=c('#999999','#E69F00'))+
      ggplot2::scale_size_manual(values=c(4,2))+
      ggplot2::theme_bw()
  }
  else{
    p <-ggplot2::ggplot(longdf, ggplot2::aes(variable, as.numeric(value)*100)) +
      ggplot2::geom_point(ggplot2::aes(shape = vartype, colour = vartype, size = vartype)) +
      ggplot2::labs(x= "Method", y = "Frequency (percentage of Trials)", shape = "Feature",size = "Feature",colour = "Feature") +
      ggplot2::scale_shape_manual(values=c(2, 16))+
      ggplot2::scale_color_manual(values=c('#999999','#E69F00'))+
      ggplot2::scale_size_manual(values=c(4,2))+
      ggplot2::theme(panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(face = "bold", angle = 90))
  }
  return(p)
}



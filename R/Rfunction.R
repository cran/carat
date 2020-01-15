getData<-function(n,cov_num,level_num,pr,type,beta,mu1,mu2,sigma = 1,method = HuHuCAR,...){
  FUN = match.fun(paste(deparse(substitute(method)),"_getData",sep = ""))
  data = FUN(n,cov_num,level_num,pr,type,beta,mu1,mu2,sigma,...)
  names = rep(NA,cov_num+2)
  for(i in 1:cov_num){
    eval(parse(text = paste('names[i]=',"paste('covariate',i,sep='')")))
  }
  names[cov_num+1] = "assignment"
  names[cov_num+2] = "outcome"
  datafr = data.frame(data,row.names = names)
  return(datafr)
}

rand.test<-function(data,Pernum = 200,method = HuHuCAR,conf = 0.95,binwidth = 30,...){
  dname = deparse(substitute(data))
  FUN = match.fun(paste(deparse(substitute(method)),"_RT",sep = ""))
  result = FUN(data,Pernum,...)
  pval = result$pval
  testmethod<-"Randomization test"
  x<-result$Randata
  datanew = data.frame(x)
  estimate = result$estimate
  names(estimate)<-"difference for treatment effect"
  pic<-ggplot2::ggplot(data = datanew, ggplot2::aes(x = x))+ 
    ggplot2::geom_histogram(bins = binwidth) + 
    ggplot2::geom_vline(ggplot2::aes(xintercept = result$estimate),
                        colour = "#990000",linetype = "dashed")
  print(pic)
  rval<-list(p.value = pval,estimate = estimate,
             method = testmethod,data.name = dname)
  rval$plot = pic
  class(rval)<-"htest"
  return(rval)
}

boot.test<-function(data,B=200,method = HuHuCAR,conf = 0.95,...){
  dname = deparse(substitute(data))
  FUN = match.fun(paste(deparse(substitute(method)),"_BT",sep = ""))
  result = FUN(data,B,...)
  testmethod<-"Bootstrap t-test"
  estimate = result[1]
  stderr = result[2]
  tstat = result[3]
  pval = result[4]*2
  cint = c(estimate + stderr*qnorm((1-conf)/2),estimate - stderr*qnorm((1-conf)/2))
  attr(cint,"conf.level")<-conf
  names(tstat)<-"t"
  names(estimate)<-"difference for treatment effect"
  rval<-list(statistic = tstat,p.value = pval,conf.int = cint,
             estimate = estimate,stderr = stderr,
             method = testmethod,data.name = dname)
  class(rval)<-"htest"
  return(rval)
}

corr.test<-function(data,conf = 0.95){
  dname = deparse(substitute(data))
  result = CTT(data)
  testmethod<-"Corrected t-test"
  estimate = result[1]
  stderr = result[2]
  tstat = result[3]
  pval = result[4]*2
  cint = c(estimate + stderr*qnorm((1-conf)/2),estimate - stderr*qnorm((1-conf)/2))
  attr(cint,"conf.level")<-conf
  names(tstat)<-"t"
  names(estimate)<-"difference for treatment effect"
  rval<-list(statistic = tstat,p.value = pval,conf.int = cint,
             estimate = estimate,stderr = stderr,
             method = testmethod,data.name = dname)
  class(rval)<-"htest"
  return(rval)
}

evalPower<-function(n,cov_num,level_num,pr,type,beta,di = seq(0,0.5,0.1),sigma = 1,Iternum,sl = 0.05,method = HuHuCAR,test,plot = "TRUE",...){
  a = Sys.time()
  if(plot!= TRUE && plot!= FALSE){
    print("Please specify whether to plot or not! Enter ON or OFF")
    return(NULL)
  }
  else{
    mu2 = rep(0,length(di))
    if(deparse(substitute(test)) == "rand.test"){
      testname = "_RT_power"
    }
    else if(deparse(substitute(test)) == "boot.test"){
      testname = "_BT_power"
    }
    else if(deparse(substitute(test)) == "corr.test"){
      testname = "_CT_power"
    }
    else{stop("Please enter a valid test! rand.test, boot.test or corr.test")}
    FUN = match.fun(paste(deparse(substitute(method)),testname,sep = ""))
    result = FUN(n,cov_num,level_num,pr,type,beta,di,mu2,sigma,Iternum,sl,...)
    if(plot == TRUE){
      diff = di
      tp = result[1:length(di)]
      sd = round(result[(length(di)+1):(2*length(di))],3)
      tgg=data.frame(diff, tp, sd)
      pic = ggplot2::ggplot(tgg, ggplot2::aes(x=di, y=tp)) + 
        ggplot2::geom_line() + 
        ggplot2::geom_point(size=4, shape=20) + 
        ggplot2::labs(x = "mu1-mu2",y = "power")
      b = Sys.time()
      result = list(Powers = tgg,Plot = pic,Time = paste(paste("Execute time:",as.numeric(b-a)),units(b-a)))
      return(result)
    }
    else{
      diff = di
      tp = result[1:length(di)]
      sd = round(result[(length(di)+1):(2*length(di))],3)
      tgg=data.frame(diff, tp, sd)
      b = Sys.time()
      result = list(Powers = tgg,Time = paste(paste("Execute time:",round(as.numeric(b-a),2),units(b-a))))
      return(result)
    }
  }
}



compPower<-function(powers,diffs,testname){
  if(is.vector(testname) == FALSE || class(testname) != "character"){
    stop("Input of testname must be a vector of character!")
  }
  if(is.vector(diffs) == FALSE || is.numeric(diffs) == FALSE){
    stop("Input of testname must be a vector!")
  }
  if(class(powers) != "list"){
    stop("Input of powers must be a list!")
  }
  else if(length(powers) != length(testname)){
    stop("The length of powers must match that of testname!")
  }
  else{
    k = length(powers)
    l = length(powers[[1]]$Powers$tp)
    for(i in 2:k){
      if(length(powers[[k]]$Powers$tp)!=l){
        stop("The length of power vectors must match!")
      }
    }
    if(length(diffs)!=l){
      stop("The length of powers and diffs must match!")
    }
    supp = NULL
    popp = NULL
    popp_out = NULL
    power_temp = rep('',l)
    for(i in 1:k){
      supp = c(supp,rep(testname[i],l))
      for(j in 1:l){
        power_temp[j] = paste(powers[[i]]$Powers$tp[j],paste("(",paste(round(powers[[i]]$Powers$sd[j],digits = 3),")",sep = ''),sep = ''),sep = '')
      }
      popp_out = c(popp_out,power_temp)
      popp = c(popp,powers[[i]]$Powers$tp)
    }
    diffp = rep(diffs,k)
  }
  tgg = data.frame(supp, diffp, popp)
  pic = ggplot2::ggplot(tgg, ggplot2::aes(x = diffp,y = popp,
                                          color = supp,shape = supp)) + 
    ggplot2::geom_line() + 
    ggplot2::geom_point(size=4) + 
    ggplot2::labs(x = "mu1-mu2",y = "power")
  tpp = t(matrix(popp_out,nrow = l))
  rownames(tpp) = testname
  tpp = data.frame(tpp)
  colnames(tpp) = diffs
  result = list(powers = tpp,plot = pic)
  return(result)
}



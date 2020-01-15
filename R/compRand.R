#Rcpp::sourceCpp('C:/Users/Ye Xiaoqing/Desktop/CAR2.3/tes.cpp')
#compRand = function(...) UseMethod("compRand")

compRand.carcomp = function(...) UseMethod("carcomp")

compRand = function(...){
  Objects = list(...); 
  clch = as.character(lapply(Objects, class)); 
  if(length(which(clch != "careval")) >= 1){
    stop("Inputs must be of class 'SUM'!")
  }
  DataG = character(); 
  DataType = character(); 
  mechanism = character(); 
  leng = length(Objects); 
  bsize = vector(); 
  nvec = vector();
  Nvec = vector(); 
  cnvec = vector();
  for(j in 1 : leng){
    R = Objects[[j]]; 
    mechanism[j] = R$method; 
    bsize[j] = R$bsize; 
    nvec[j] = R$N; 
    Nvec[j] = R$iteration; 
    DataType[j] = R$`Data Type`;
    cnvec[j] = R$cov_num; 
    if(is.null(R$DataGeneration)){
      DataG[j] = NA;
    }else{
      DataG[j] = R$DataGeneration;
    }
  }
  if(length(unique(cnvec)) == 1){
    cov_num = cnvec[1]; 
    level_num = Objects[[1]]$level_num; 
  }else{
    warning("Results don't make sense: different dataframes are used in different methods."); 
    cov_num = cnvec; 
    level_num = NA; 
  }
  if(length(unique(DataType)) > 1){
    warning("Results don't make sense: comparison between simulated data and real data.")
  }
  if(length(unique(nvec)) > 1){
    n = NA; 
    warning("Results don't make sense: different sample sizes for different methods.")
  }else{
    n = nvec[1];
  }
  if(length(unique(Nvec)) > 1){
    N = NA; 
    warning("Results don't make sense: different iterations for different methods.")
  }else{
    N = Nvec[1]; 
  }
  bmax = max(bsize); 
  cname = vector()
  cname[1 : 5] = c("max", "95%-quan", "median", "mean", "loss");
  for(j in 1 : bmax){
    cname[5 + j] = paste("num", "=", j, seq = " ");
  }
  
  C_A = C_O = C_M = matrix(NA, nrow = leng, ncol = 5); 
  C_S = matrix(NA, nrow = leng, ncol = 5 + bmax);
  rownames(C_A) = rownames(C_O) = rownames(C_M) = rownames(C_S) = mechanism; 
  colnames(C_A) = colnames(C_O) = colnames(C_M) = c("max", "95%-quan", "median", "mean", "loss");
  colnames(C_S) = cname; 
  for(i in 1 : leng){
    R = Objects[[i]];
    C_A[i, ] = apply(R$Imb, 2, mean); 
    C_O[i, ] = R$Imb[1, ]; 
    C_M[i, ] = apply(R$Imb[(2 + R$strt_num) : (sum(R$level_num) + 1 + R$strt_num), ], 2, mean); 
    C_S[i, 1 : 5] = apply(R$Imb[2 : (1 + R$strt_num), ], 2, mean); 
    C_S[i, 6 : (5 + R$bsize)] = R$`Within-strt. by num of pats`; 
  }
  
  meth = rep(mechanism, times = 4);
  levels = rep(c("A", "O", "M", "S"), each = leng);
  crt = rep(1 : 4, each = leng);
  mean = c(C_A[, 4], C_O[, 4], C_M[, 4], C_S[, 4]);
  loss = c(C_A[, 5], C_O[, 5], C_M[, 5], C_S[, 5]);
  df1 = data.frame(meth, levels, crt, mean, loss); 
  
  p1 = ggplot2::ggplot(df1, ggplot2::aes(x = meth, y = mean, color = levels, group = crt,
                                         linetype = levels, shape = levels)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::xlab("method") +
    ggplot2::ylab("absolute mean") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                   legend.position = "bottom")
  
  val = vector()
  for(k in 1 : bmax){
    val = c(val, C_S[, 5 + k]);
  }
  num = rep(cname[6 : (5 + bmax)], each = leng);
  meths = rep(mechanism, times = bmax)
  df2 = data.frame("meth" = meths, "num" = num, "val" = val);
  
  p2 = ggplot2::ggplot(df2, ggplot2::aes(x = num, y = val, color = meth, group = meth,
                                         linetype = meth, shape = meth)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::xlab("numbers of patients for each strata") +
    ggplot2::ylab("absolute within-stratum mean") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                   legend.position = "bottom")
  
  p = gridExtra::grid.arrange(p1, p2, ncol = 1);
  
  RR = list("All Imbalances" = C_A, "Overall Imbalances" = C_O, 
            "Marginal Imbalances" = C_M, "Within-stratum Imbalances" = C_S);
  RR$plot = p;
  RR$mechanism = mechanism; 
  RR$N = n;
  RR$iteration = N;
  RR$cov_num = cov_num; 
  RR$level_num = level_num; 
  RR$DataType = DataType; 
  RR$DataGeneration = DataG; 
  
  class(RR) = "carcomp";
  return(RR)
}

fbagplot <- function(fd, inflation_factor = 2.58){
  #time is the domain of the observed data
  #fd are the functional data objects
  #inflation factor is for the bagplot inflation for outlier detection
  
  require(aplpack) #for the bagplot
  
  pca_fd <- pca.fd(fd)
  bagplot <- bagplot(pca_fd$scores[,1], pca_fd$scores[,2], inflation_factor,
                     create.plot = FALSE)
  if(is.null(bagplot$pxy.outlier)){ #if no outliers
    outliers <- NA
  }else{
    outliers <- which(pca_fd$scores[,1] %in% bagplot$pxy.outlier)
  }
  return(outliers)
}

#the tvdmss function is taken from the fdaoutlier package but with a bug fix
#for if the shape_boxstats$out < mean(mss) condition is false
tvdmss <- function (dts, emp_factor_mss = 1.5, emp_factor_tvd = 1.5, 
                    central_region_tvd = 0.5) {
  depths_mss <- total_variation_depth(dts = dts)
  tvd <- tvd_old <- depths_mss$tvd
  mss <- depths_mss$mss
  dta_dim <- dim(dts)
  n_curves <- dta_dim[1]
  n_points <- dta_dim[2]
  index <- (1:n_curves)
  n_central_obs <- ceiling(n_curves/2)
  shape_boxstats <- boxplot(mss, range = emp_factor_mss, plot = F)
  shape_outliers <- NULL
  if (length(shape_boxstats$out) != 0) {
    shape_outliers <- which(mss %in% shape_boxstats$out[shape_boxstats$out < 
                                                          mean(mss)])
    if(length(shape_outliers) == 0){  #no shape outliers less than mean(mss)
      shape_outliers <- NULL
    }else{
      dts <- dts[-shape_outliers, ]
      tvd <- tvd[-shape_outliers]
      index <- index[-shape_outliers]
    }
  }
  magnitude_outliers <- NULL
  outliers <- functional_boxplot(dts, depth_values = tvd, emp_factor = emp_factor_tvd, 
                                 central_region = central_region_tvd * n_curves/nrow(dts))$outliers
  if (length(outliers) != 0) 
    magnitude_outliers = index[outliers]
  return(list(outliers = sort(c(magnitude_outliers, shape_outliers)), 
              shape_outliers = shape_outliers, magnitude_outliers = magnitude_outliers, 
              tvd = tvd_old, mss = mss))
}

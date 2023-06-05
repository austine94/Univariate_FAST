#####
#Implements anomaly detection as in Huang (2019), using the fdaoutlier package

tvd_anomaly_detection <- function(fd_eval, mss_scale_factor = 3,
                                  tvd_scale_factor = 1.5){
  
  #fd_eval is the functional data evaluated at the time points
  #we input the data as in Ramsay fda package, so columns are obs rows are time
  #scale factor is the factor to scale the IQR by to set the detection limits,
  #Huang (2019) suggest 3
  
  #tvd and mss calculated as per the fdaoutlier package
  
  tvd <- total_variation_depth(t(fd_eval))
  
  mss_third <- quantile(tvd$mss, 0.75)
  mss_first <- quantile(tvd$mss, 0.25)
  mss_iqr <- mss_third - mss_first
  threshold_lower_mss <- mss_first - mss_scale_factor*mss_iqr  
  threshold_upper_mss <- mss_third + mss_scale_factor*mss_iqr
  
  tvd_third <- quantile(tvd$tvd, 0.75)
  tvd_first <- quantile(tvd$tvd, 0.25)
  tvd_iqr <- tvd_third - tvd_first
  threshold_lower_tvd <- tvd_first - tvd_scale_factor*tvd_iqr  
  threshold_upper_tvd <- tvd_third + tvd_scale_factor*tvd_iqr
  
  mss_anomalies <- which(tvd$mss < threshold_lower_mss) #Huang only uses lower threshold for shape
  tvd_anomalies <- which((tvd$tvd > threshold_upper_tvd) | (tvd$tvd < threshold_lower_tvd))
  
  return(list(mss_anomalies = mss_anomalies, tvd_anomalies = tvd_anomalies))
  
}
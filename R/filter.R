library(tseries)
library(lmtest)

wrapper.filter.splines <- function(raw_data, LMMSObject, stationnarity.test = FALSE, homoskedasticity = TRUE, MSE.filter = TRUE, 
                                   stationnarity.cutoff = 0.05, homoskedasticity.cutoff = 0.05){
  # raw_data :
  # ----------
  # sample_time in rows X molecules in columns
  # WARNING : sample_time
  # note to myself : should be in LMMS spline object
  
  # check for LMMS Spline object
  stopifnot(is(LMMSObject,"lmmspline"))
  
  result <- rownames(LMMSObject@predSpline) %>% as.data.frame() %>% purrr::set_names(c("molecule")) %>%
    mutate(modelsUsed = LMMSObject@modelsUsed)
  predSpline <- LMMSObject@predSpline %>% t %>% as.data.frame()
  
  # stationnarity.test
  # ------------------
  # ADF test : if Pvalue < signif thresold : non stationnory -> to keep
  # enough to remove straight lines and lots of splines
  # warning :  VERY stringent, pvalue correction applied
  if(stationnarity.test){
    ADF.res <- suppressWarnings(lapply(predSpline, tseries::adf.test)) %>% 
      lapply(function(x) x$p.value) # Extract pvalue
    result <- (ADF.res <= stationnarity.cutoff) %>% as.data.frame() %>% rownames_to_column() %>% purrr::set_names(c("molecule", "ADF.test")) %>% 
      right_join(result, by = c("molecule", "molecule")) # merge result
  }
  
  # homoskedasticity test
  # ---------------------
  #  Breusch-Pagan test. # WARNING : only for linear regression model (modelsUsed == 0)
  # * if pvalue < signif cutoff : heteroskedasticity
  # * if pvalue > signif cutoff : homooskedasticity  -> to keep, no pvalue correction
  if(homoskedasticity){
    models0 <- LMMSObject@models[LMMSObject@modelsUsed == 0]
    BP.res <- lapply(models0, lmtest::bptest) %>% 
      lapply(function(x) as.numeric(x$p.value)) %>% unlist() %>%
      as.data.frame() %>% set_names("BP.test") %>%
      mutate(molecule = rownames(LMMSObject@predSpline)[LMMSObject@modelsUsed == 0])
    result <- result %>% left_join(BP.res, by = c("molecule", "molecule")) %>% # merge result
      mutate(BP.test = ifelse(is.na(BP.test), 1, BP.test)) %>% # replace (model != 0) by  a p value of 1, must be homoskedastic if not, no spline
      mutate(BP.test = (BP.test >=  homoskedasticity.cutoff))  # TRUE/FALSE
  }
  
  # compute MSE in every case, always usefull in the event of a figure
  MSE.res <- get_MSE(raw_data, LMMSObject)
  # applied filter based on max MSE for model != 0
  if (MSE.filter) {
    MSE.cutoff <- MSE.res %>% filter(model_used != 0) %>% pull(MSE) %>% max
    result <- MSE.res %>% mutate(MSE.filter = (MSE <= MSE.cutoff)) %>% # TRUE/FALSE MSE.cutoff
      dplyr::select(molecule, MSE.filter) %>% right_join(result)
  }
  
  # conclusion
  OK <- c(stationnarity.test, homoskedasticity, MSE.filter) %>% as.numeric() %>% sum  # number of test that should be passed
  result %>% gather(test, res, -c(molecule, modelsUsed)) %>% 
    mutate(res = as.numeric(res)) %>%
    group_by(molecule) %>% 
    summarise(val = sum(res)) %>% # count number of positive test
    filter(val == OK) %>% # keep only molecule that passed all the test
    pull(molecule) -> to_keep
  
  return(list(to_keep =to_keep, res.filter = result))
}

get_MSE <- function(raw_data, LMMSObject){
  # sample_time in rows X molecules in columns
  # WARNING : sample_time
  # note to myself : should be in LMMS spline object
  
  # check for LMMS Spline object
  stopifnot(is(LMMSObject,"lmmspline"))
  
  raw_data %>% rownames_to_column("sample") %>%
    mutate(ID = sample %>% str_split("_") %>% map_chr(~.x[1])) %>%
    mutate(time = sample %>% str_split("_") %>% map_chr(~.x[2]) %>% as.numeric) %>%
    gather(molecule, Yi, -c(time, ID, sample)) %>%
    left_join(LMMSObject@predSpline %>% 
                rownames_to_column("molecule") %>% 
                mutate(model_used = factor(LMMSObject@modelsUsed)) %>%
                gather(time, Y_hat, -c(molecule, model_used)) %>% 
                mutate(time = as.numeric(time))) %>%
    mutate(error = (Yi-Y_hat)^2) %>% 
    group_by(molecule, model_used) %>%
    summarise(MSE = mean(error))
}

## Script for plotting output from OM & Diagnosing problems with discrepencies
## Brett Stacy
## 18_06_2019

#' Von Bertalanffy weight 2. This one doesn't require "object" like calc_VBweith()
#'
#' Growth may vary seasonally, therefore the whole (FLQuant-like)
#' object is passed on
#'
#' NOTE: params1 and params2 should reference names not position
#' @param age1 age
#' @param params1 von Bert parameters
#' @param params2 weight length parameters
#' @export
calc_VBweight2 <- function(age1,params1,params2) {
  Linf <- params1[1]; K   <- params1[2]; t0 <- params1[3]; CV <- params1[4]  # CV not used here
  WLa  <- params2[["a"]]; WLb <- params2[["b"]]
  res  <- WLa*(Linf*(1-exp(-K*(age1-t0))))^WLb
  return(res)
}

#' Calculate SSB
#'
#' Manually calculate SSB from age frequency data.
#'
#' This function calculates SSB from numbers at age, growth, weight at age, maturity at age.
#' The input data can be from the OM res object or from the casal output log.
#' @param type Type of object. One of c("res", "output_log").
#' @param age_data Object containting age frequency data.
#' @param VB_params vector containing parameters for length at age relationship, Von-Bertalanffy.
#' @param WL_params list containing parameters for weight at length relationship.
#' @param M_ogive type of ogive the maturity follows, probably "logistic".
#' @param M_params list containing parameters for the maturity at age relationship.
#' @return Returns a time series of SSB
#' @export
calc_SSB = function(type, age_data, VB_params, WL_params, M_ogive, M_params){
  switch(type,
         res = {
           wt  <- calc_VBweight2(as.numeric(rownames(age_data)), VB_params, WL_params); # calc weight from calc_VBweight
           m   <- ogive(type = M_ogive, ages = as.numeric(rownames(age_data)), params = M_params);
           wtm <- wt*m;
           SSB <- colSums(sweep(age_data, 1, wtm, "*"))
         },
         output_log = {
           wt  <- calc_VBweight2(seq(1, length(age_data[[1]])), VB_params, WL_params); # calc weight from calc_VBweight
           m   <- ogive(type = M_ogive, ages = seq(1, length(age_data[[1]])), params = M_params);
           wtm <- wt*m;
           SSB <- unlist(lapply(lapply(age_data, "*", wtm), sum))
         })

  return(SSB)
}
# VB_params = c(1565, 0.146, 0.015, 0.122)
# WL_params = list(a = 3.0088e-12, b = 3.2064)
# M_ogive   = "logistic"
# M_params  = list(x50 = 14.45, x95 = 6.5)
# age_data = res.TOA$pop$n[,,1,,]*2
# age_data = output.TOA$Numbers_at_age_R1



#' SSB Time Series
#'
#' Plot SSB from output from Planetfish2 operating model or assessment model.
#'
#' This function plots data from the \code{output} matrix. This can be either OM or AM output as SSB or REC etc.
#' The intention is to be able to plot an item by passing the unique column name identifier to \code{plot()}.
#' @param output The output matrix from a model run.
#' @param item = OM_ssb_R1 as default but can be any known column name type.
#' @param mean = Logical. Should the mean be plotted?
#' @return Returns a time series plot of item.
#' @export
plot_SSB = function(output, item = "OM_ssb_R1", mean = TRUE, ...){
  if(mean == T){
    get_item = output[, grepl(item, colnames(output))]
    if(NCOL(get_item) == 1){
      plot(get_item, xlab = "Year", ylab = item, type = "l", ...)
    }
    else {
      plot(colMeans(get_item), xlab = "Year", ylab = item, type = "l", ...)
    }
  }
  else {
    get_item = output[, grepl(item, colnames(output))]
    boxplot(get_item, xlab = "Year", ylab = item, type = "l", xaxt = "n", ...)
    axis(1, at = seq(10, round(dim(get_item)[2], -1), 10))
  }
}

# plot_OM_SSB(output)





#' Age structure of a population.
#'
#' Plot histogram of age frequencies for a particular region.
#'
#' This function plots age frequency data from the \code{output2} array. \code{output2} should be of dimention [n_iters, ages, years, regions]
#' @param output2 The output array from a model run.
#' @param region Which region should plotted?
#' @param skip Skip years if needed. = 1 means don't skip any years.
#' @return Returns a histogram of age structure for every year.
#' @export
plot_OM_Ages = function(output2, region, skip){
  par(mfrow = c(4,5), mar = c(.5, .5, .5, .5))
  for (i in seq(1, dim(output2)[3], skip)) {
    mean_iter = apply(output2, c(2,3,4), mean)
    barplot(mean_iter[,i,region], yaxt = "n", xaxt = "n")
    legend("topright", paste("year:", i), bty = "n")
  }
}

# plot_OM_Ages(output2, 1)


#' Relative Error SSB
#'
#' Calculate relative error for SSB
#'
#' This function needs a description.
#' @param output The output matrix from a model run.
#' @param type Switch. Type of SSB index. Options: "initial", "current", "status".
#' @return Returns a vector distribution of relative error.
#' @export
SSB_err = function(output, type){
  switch(type,
         initial = {om_output <- output[, grepl("OM_ssb_R1_1990", colnames(output))];
         am_output <- output[, grepl("AM_ssb_1990", colnames(output))]},

         current = {om_output <- output[, grepl("OM_ssb_R1", colnames(output))];
         om_output <- om_output[, ncol(om_output)];
         am_output <- output[, grepl("AM_ssb_", colnames(output))];
         am_output <- am_output[, ncol(am_output)]},

         status = {om_output <- output[, grepl("OM_ssb_R1", colnames(output))]/(output[, grepl("OM_ssb_R1_1990", colnames(output))]);
         om_output <- om_output[, ncol(om_output)];
         am_output <- output[, grepl("AM_ssb_", colnames(output))]/output[, grepl("AM_ssb_1990", colnames(output))];
         am_output <- am_output[, ncol(am_output)]})

  err = (am_output - om_output)/om_output

  return(err)

}

# plot_err(output, truth = "OM_ssb0", est = "AM_ssb0", half = "yes")


#' Relative Error Plot
#'
#' Plot relative error statistic between SSBs
#'
#' This function uses data from the \code{output} matrix to plot the relative error between operating and assessment model output quantities.
#' This will need adaptation as more results are acquired. The main one will be to have the option to vectorize the boxplots like Pauls.
#' @param output The output matrix from a model run.
#' @param type Switch. Type of SSB comparison. Options: "initial", "current", "status".
#' @return Returns a boxplot of relative error.
#' @export
plot_SSB_err = function(output, type, ...){

  switch(type,
         initial = {om_output <- output[, grepl("OM_ssb0", colnames(output))];
         am_output <- output[, grepl("AM_ssb0_", colnames(output))]},

         current = {om_output <- output[, grepl("OM_ssb_R1", colnames(output))];
         om_output <- om_output[, ncol(om_output)];
         am_output <- output[, grepl("AM_ssb_", colnames(output))];
         am_output <- am_output[, ncol(am_output)]},

         status = {om_output <- output[, grepl("OM_ssb_R1", colnames(output))]/(output[, grepl("OM_ssb0", colnames(output))]);
         om_output <- om_output[, ncol(om_output)];
         am_output <- output[, grepl("AM_ssb_", colnames(output))]/output[, grepl("AM_ssb0", colnames(output))];
         am_output <- am_output[, ncol(am_output)]})

  err = (am_output - om_output)/om_output

  boxplot(err, ...)


}




## Functions for checking parameter settings in para.

#' Check Match
#'
#' Check perfect knowledge between OM and AM
#'
#' check_match checks if certain parameters match between Planetfish2 and CASAL.
#' Perfect knowledge between the two models of parameters like life history values, fishing selectivity,
#' and movement are necessary for a baseline comparison.
#' @param para Parameter settings at the point just before running models.
#' @return Returns printed output listing matched/mismatched parameters
#' @export
#' @examples
#' # With default settings:
#' para = get_om_data()
#' para = get_casal_para(para)
#' check_match(para)
check_match = function(para){
  print("##### POPULATION ####")
  ifelse(identical(para$om$age, para$ass$age), print("age match"), print("age mismatch"))
  ifelse(identical(para$om$ages, para$ass$ages), print("ages match"), print("ages mismatch"))
  ifelse(identical(para$om$n_ages, length(para$ass$ages)), print("n_ages match"), print("n_ages mismatch"))
  ifelse(identical(para$om$move_rules, para$ass$move_rules), print("move_rules match"), print("move_rules mismatch"))

  ifelse(identical(para$om$growth[[1]], para$ass$estgrowth[[1]]), print("growth match"), print("growth mismatch"))

  ifelse(identical(unname(unlist(para$om$WL[[1]])), para$ass$estWL[[1]]), print("WL match"), print("WL mismatch"))

  ifelse(identical(para$om$pin_mat, para$ass$estpin.mat), print("maturity ogive match"), print("maturity ogive mismatch"))
  ifelse(identical(unlist(para$om$maturity[[1]], use.names = F), unlist(para$ass$estmaturity, use.names = F)), print("maturity match"), print("maturity mismatch"))
  ifelse(identical(c("allvalues ", round(ogive(para$om$pin_mat, para$om$ages, para$om$maturity[[1]]), 4)),
                   para$ass$maturity_props_all), print("maturity props match"), print("maturity props mismatch"))

  ifelse(identical(para$om$natM[1], para$ass$estnatM[[1]]), print("natM match"), print("natM mismatch"))

  ifelse(identical(para$om$rec_h, para$ass$rec_steepness), print("rec_h match"), print("rec_h mismatch"))

  print("##### TAGGING ####")
  ifelse(identical(para$sampling$tag_shedding, para$ass$tag_shedding_rate), print("tag shedding match"), print("tag shedding mismatch"))
  ifelse(identical(para$sampling$tag_mort, para$ass$tag_mortality), print("tag mortality match"), print("tag mortality mismatch"))



  print("##### SELECTIVITY ####")
  ifelse(all(sapply(list(unname(unlist(para$sampling$tag_select[1])), as.numeric(para$ass$selN_all[[1]][-1])), identical, unname(unlist(para$om$select[1])))),
         print("Selectivity Numbers match"), print("Selectivity Numbers mismatch"))
  ifelse(all(sapply(list(unname(unlist(para$sampling$pin_tag_sel[1])), para$ass$selN_all[[1]][1]), identical, unname(unlist(para$om$pin_sel[1])))),
         print("Selectivity Type match"), print("Selectivity Type mismatch"))



  ### Add these at some point if neccessary. Perhaps in a check_estim() function to check which parameters the user is estimating. (BS 3/5/19)
  # para$ass$estim_natural_mortality.all # what the heck is this?
  # para$ass$estim_size_at_age.cv
  # para$ass$estim_natural_mortality.ogive_all
  # para$ass$estimate_selectivity
}


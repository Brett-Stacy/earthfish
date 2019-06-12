## Function for running operating and assessment model.
## BS 6/5/19


#' Output SSB
#'
#' Run operating and assessment models to acquire SSB output.
#'
#' output_SSB runs the operating model and CASAL assessment model to generate data
#' for comparison. First and foremost, this is a trial function to streamline the process
#' of running the models and saving SSB output for Brett's PhD projects. Alterations will be
#' required if investigations outside the scope of Brett's PhD are desired. This function takes
#' as input the `para` parameter settings from a particular species. These settings should be the
#' "default" settings from which adjustment scenarios can be allocated. For example, a typical use could
#' use Patagonia Toothfish life history parameters, fishing, and sampling specifications reflecting
#' the stock on  the Kerguelen Plateau. Modifications to the baseline para object related to tagging
#' or otherwise can be made to facilitate inspection of the impact of hypothetical sampling scenarios.
#' @param baseline_para Load baseline parameterizaion para into environment.
#' @param n_tags Number of tags per year.
#' @param n_years_tags Number of years to apply tags.
#' @param iterations Number of iterations.
#' @param species Species of interest.
#' @param casal_path Point to folder with CASAL items.
#' @param save_output_path Point to the folder in which to save output.
#' @return Returns a matrix of SSB values generated from the OM and estimated from the AM.
#' @export
output_SSB = function(baseline_para, n_tags, n_years_tags, iterations, species, casal_path, save_output_path){  ## could change n_tags and n_years_tags to: "scenario_type", and "scenario_value"


  # library(planetfish2)
  # library(casal)


  # Establish file name
  file_name = paste0(species, "_tags_", n_tags, "_years_", n_years_tags, "_Niter_", iterations)

  # Set iterations
  n_iters = iterations

  # # Load up para into function environment
  # load("C:/Users/STA384/OneDrive - University of Tasmania/PhD/Year1/Work/Chapter_1_simulations/Patagonia_Toothfish/Baseline_Para/baseline_para_PT.Rda")
  para = baseline_para

  # Set Casal path
  para$control$casal_path = casal_path
  para$control$assess_type = "CASAL"

  # Set number of tags in OM
  para$sampling$tag_N <- c(n_tags, 0)


  # Set number of years to release tags
  tag_years = (para$om$year[2] - n_years_tags + 1):para$om$year[2] - 1
  n_years_aged = 10
  age_years = (para$om$year[2] - n_years_aged):para$om$year[2]
  para$ass$sample_years = am_sampling(years = para$om$years,
                                      ycurr = para$om$year[2],
                                      catchage_yrs = age_years,
                                      tagging_yrs = tag_years)$sample_years







  # Setup output object
  dim_names	<- c("OM_ssb0", paste0("OM_ssb_R1_", para$om$years), paste0("OM_ssb_R2_",para$om$years),
                 paste0("OM_rec_R1_", para$om$years), paste0("OM_rec_R2_", para$om$years),
                 "AM_ssb0_",paste0("AM_ssb_", para$om$years), paste0("AM_rec_", para$om$years))

  dim_length <- length(dim_names)

  ## construct the output array
  output <- array(data = 0, dim = c(n_iters, dim_length),
                  dimnames = list("Iter"=1:n_iters, dim_names))
  ## some conveniences for accessing arrays in the OM
  R1 <- 1
  R2 <- 2
  S1 <- para$om$season[2]
  ## paths for the CASAL outputs
  casal_path <- para[["control"]]$casal_path
  mpd_dat <- para[["control"]]$mpd_dat
  output_log <- para[["control"]]$output_log


  ## loop over the number of iterations
  for(i_iter in 1:n_iters){
    ## Set up om objects
    res	<- setup_om_objects(para=para)
    #### Populate om objects (mod) with biological data (M and fecundity), fishery effort & selectivity, and observation sample sizes
    res	<- populate_om_objects(para=para, res=res)
    #### Get initial population numbers and calculate SSB0
    res	<- get_initial_pop(para=para, res=res)
    #### Run Annual OM loops with/without assessment
    para[["control"]]$pin_casal_assess <- 1
    res <- run_annual_om(para=para, res=res) #, intern=TRUE) #


    ## OM calculated quantities
    ssb <- res$mod$ssb # ssb[quant, year, sex, season, area] # rec is the same
    rec <- res$mod$rec
    ## Spawning Biomass
    OM_SSB0 <- res$mod$ssb0
    OM_SSB_R1 <- apply(ssb[1,,,S1,R1],c(1),sum)
    OM_SSB_R2 <- apply(ssb[1,,,S1,R2],c(1),sum)
    ## Recruitment
    OM_Rec_R1 <- apply(rec[1,,,S1,R1],c(1),sum)
    OM_Rec_R2 <- apply(rec[1,,,S1,R2],c(1),sum)


    ## length of output
    om_ncols <- length(OM_SSB0) + length(OM_SSB_R1) + length(OM_SSB_R2) + length(OM_Rec_R1) +  length(OM_Rec_R2)

    output[i_iter, 1:om_ncols] <- c(OM_SSB0, OM_SSB_R1, OM_SSB_R2, OM_Rec_R1, OM_Rec_R2)



    ## add the AM output if it exists for this iteration
    if(file.exists(paste0(casal_path, para[["control"]]$mpd_dat)) &
       length(scan(file = paste0(casal_path, para[["control"]]$mpd_dat),
                   what = "character")) > 0){
      nname1 <- para[["control"]]$output_log
      casal_quants <- casal::extract.quantities(file=nname1, path=casal_path)
      casal_freeparams <- casal::extract.free.parameters(file=nname1, path=casal_path)
      ## quantities form the Assessment model
      AM_SSB0 <- casal_quants$B0
      AM_SSB_R1 <- casal_quants$SSBs$SSB
      AM_Rec_R1 <- casal_quants$recruitments$recruitment

      ## add the rest of the output
      output[i_iter, (om_ncols+1):ncol(output)] <- c(AM_SSB0, AM_SSB_R1, AM_Rec_R1)
    }
  }



  ## write to file
  write.csv(output, file=paste0(save_output_path, file_name,".csv"),
            quote=FALSE, na="NA", row.names=FALSE)



}


# baseline_para    = readRDS("C:/Users/STA384/OneDrive - University of Tasmania/PhD/Year1/Work/Chapter_1_simulations/Patagonia_Toothfish/Baseline_Para/baseline_para_PT_paul.Rds")
# n_tags           = 2500
# n_years_tags     = 5
# iterations       = 1
# species          = "PT"
# casal_path       = "C:/Users/STA384/Work/Chapter_1_output/output_using_function/"
# # casal_path       = "C:/Users/STA384/Work/simulations_test_PC/baseline_complementary_assessments_2/"
# save_output_path = "C:/Users/STA384/OneDrive - University of Tasmania/PhD/Year1/Work/Chapter_1_simulations/Patagonia_Toothfish/Scenarios/vary_tags/"
#
#
# output_SSB(baseline_para, n_tags, n_years_tags, iterations, species, casal_path, save_output_path)










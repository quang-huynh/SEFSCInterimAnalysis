######## R version 3.5.3

######## Install packages
# Note that this version of DLMtool is different than the version on CRAN.
# It has been modified so that the index follows vulnerable biomass instead of total biomass
# install.packages("DLMtool_5.3.1.tar.gz", repos = NULL)
# install.packages("MSEtool_1.4.1.tar.gz", repos = NULL)

library(MSEtool)
rm(list=ls())
t1 <- Sys.time()
set.seed(1358)

# Setup loops
scenario <- c("base", "hs", "hd", "lf", "epiM","dep")
OM_name <- c("OM_BlackSeaBass","OM_SnowyGrouper","OM_VermilionSnapper",
             "OM_RedPorgy","OM_GoldenTilefish","OM_Cobia",
             "OM_BluelineTilefish","OM_Gag","OM_GrayTriggerfish"
             )
seeds <- setNames(sample(1:10000,length(OM_name),replace = FALSE),OM_name)
ntrials <- 1500 # Sometimes this may need to be a really large number like 1500

DLMtool::setup(12) # Run in parallel over 12 cores


for(OM_name_k in "OM_Gag") { ######### Loop over operating model
  source('fn/iMP.R')

  for(scenario_i in scenario[2:6]) { ######### Loop over scenario
    myOM <- readRDS(paste0("OM/", OM_name_k, ".rds"))

    if(scenario_i == "hs") myOM@beta <- c(1/3, 2/3) # beta values below 1 lead to hyperstability
    if(scenario_i == "dep") myOM@cpars$D <- 0.3 * myOM@cpars$D
    if(scenario_i == "hd") myOM@beta <- c(1.5, 3)
    if(scenario_i == "lf") myOM@cpars$D <- 2 * myOM@cpars$D
    if(scenario_i == "epiM") {
      set.seed(seeds[OM_name_k])
      M_mult <- rbinom(myOM@proyears * myOM@nsim, 1, 0.1) * pmin(exp(rnorm(myOM@proyears * myOM@nsim, 0, 2)), 4)
      M_y <- myOM@M[1] * (1 + M_mult)
      M_array_hist <- array(myOM@M[1], dim = c(myOM@nsim, myOM@maxage, myOM@nyears))
      M_array_future <- aperm(array(M_y, dim = c(myOM@nsim, myOM@proyears, myOM@maxage)), perm = c(1, 3, 2))
      myOM@cpars$M_ageArray <- abind::abind(M_array_hist, M_array_future, along = 3)
    }

    ######## Fixed TAC MPs and Averaged Index MPs
    myOM@interval <- c(5, 10, 1, 1)
    # Save OM to object
    OM_name_ki <- paste0(OM_name_k, "_", scenario_i)
    assign(OM_name_ki,myOM)
    # Save OM to file
    saveRDS(get(OM_name_ki),
         file=paste0("OM_modified/",paste0(OM_name_ki, ".rds")))

    message(paste("Fixed TAC and Averaged MPs for:", OM_name_ki))
    t2 <- Sys.time()
    message(paste0("at: ",t2,".(",round(t2-t1,2)," since start)"))

    MSE_batch_1 <- runMSE(myOM, MPs = c("SCA_5", "SCA_10", "iMP_avg_5", "iMP_avg_10"),
                          parallel = TRUE, PPD = TRUE, ntrials = ntrials)
    t3 <- Sys.time()
    message(paste0("batch 1 finished at ",t3,".(",round(t3-t2,2)," duration"))

    ######## Annual Assessment MP
    myOM@interval <- 1
    MSE_batch_2 <- runMSE(myOM, MPs = "SCA_1", parallel = TRUE, PPD = TRUE, ntrials = ntrials)
    t4 <- Sys.time()
    message(paste0("batch 2 finished at ",t4,".(",round(t4-t3,2)," duration"))

    ####### Buffered interim MPs
    myOM@interval <- 1
    MSE_batch_3 <- runMSE(myOM, MPs = c("iMP_buffer_5", "iMP_buffer_10"), parallel = TRUE, PPD = TRUE, ntrials = ntrials)
    t5 <- Sys.time()
    message(paste0("batch 3 finished at ",t5,".(",round(t5-t4,2)," duration"))

    ######## Projection MPs
    myOM@interval <- 1
    MSE_batch_4 <- runMSE(myOM, MPs = c("pMP_5", "pMP_10"), parallel = TRUE, PPD = TRUE, ntrials = ntrials)
    t6 <- Sys.time()
    message(paste0("batch 4 finished at ",t6,".(",round(t6-t5,2)," duration"))

    ######## Merge MSE and save output
    res <- merge_MSE(MSE_batch_1,
                     MSE_batch_2,
                     MSE_batch_3,
                     MSE_batch_4
                     )
    saveRDS(res, file = paste0("MSE_obj/", OM_name_k, "_", scenario_i, ".rds"))

  }
}

sfStop()

# This is just a very slight modification to James Thorson's code of this function
# so that it can plot a single panel

Plot_taxa_NK <- 
function (Taxa,
          prob = 0.95,
          params = matrix(c("K", "M", 
                            "K", "Loo",
                            "tmax", "tm", 
                            #"Lm", "Temperature",
                            "ln_margsd", "rho", 
                            "logitbound_h", "ln_Fmsy",
                            "logitbound_h", "M"),
                          ncol = 2, byrow = TRUE),
          params.obs = NULL, # Optional observed values of params to plot
          legend.box.col=rgb(0,0,0,0.2),
          legend.bg=     rgb(1,1,1,0.9),
          Database = FishLife::FishBase_and_RAM,
          Cov_gjj = Database$Cov_gvv, 
          Mean_gj = Database$beta_gv,
          ParentChild_gz = Database$ParentChild_gz, 
          Y_ij = Database$Y_ij,
          g_i = Database$g_i,
          xlim = log(c(0.01,2)),
          ylim = xlim,
          ticks = c(0, 5),
          logticks = c(1, 2, 5),
          partial_match = FALSE,
          drop_pred = TRUE,
          mfrow = c(nrow(params), 1),
          legendnum = 1,
          verbose = FALSE,
          plot_lines = FALSE, 
          lcol = rainbow(length(Taxa)),
          lty = rep("solid", length(Taxa)), 
          xaxt = "s",
          yaxt = "s", 
          ...) 
{
  use_row = apply(params, MARGIN = 1, FUN = function(charvec) {
    all(charvec %in% colnames(Mean_gj))
  })
  params = params[which(use_row), ,drop=FALSE] # NK added drop=FALSE
  par(mfrow = mfrow, mar = c(3, 3, 0, 0), mgp = c(1.75, 0.25, 
                                                  0), tck = -0.02, oma = c(1, 1, 1, 1))
  Pred_taxa = NULL
  for (rowI in 1:nrow(params)) {
    for (uniqueI in 1:length(unique(Taxa))) {
      Pred_taxa[[uniqueI]] = Plot_trait(Taxon = Taxa[uniqueI], 
                                        params = params[rowI, ],
                                        Cov_gjj = Cov_gjj,
                                        Mean_gj = Mean_gj, 
                                        ParentChild_gz = ParentChild_gz,
                                        Y_ij = Y_ij, 
                                        g_i = g_i,
                                        add = ifelse(uniqueI == 1, FALSE, TRUE),
                                        xlim = quantile(Mean_gj[, params[rowI,1]], na.rm = TRUE, c(0, 1)),
                                        ylim = quantile(Mean_gj[,params[rowI, 2]], na.rm = TRUE, prob = c(0,1)),
                                        partial_match = partial_match,
                                        main = "", 
                                        lcol = lcol[uniqueI],
                                        ticks = ticks,
                                        logticks = logticks, 
                                        plot_lines = plot_lines,
                                        verbose = verbose,
                                        prob = prob, 
                                        lty = lty[uniqueI],
                                        xaxt = xaxt,
                                        yaxt = yaxt, 
                                        ...)
      # Add mean predictions for all taxa
      xy <- Pred_taxa[[uniqueI]]$Mean_pred[params[rowI, ]]
      points(x=xy[1],y=xy[2],pch=16,col=lcol[uniqueI])
      # Save mean predictions for finest taxon
      if(uniqueI==1){
        # Simple function to scale parameters to normal space
        scaleparam <- function(x,param){
          switch(param,
        "K"=exp(x),
        "M"=exp(x), 
        "Winfinity"=exp(x),
        "Loo"=exp(x),
        "tmax"=exp(x),
        "tm"=exp(x), 
        "Lm"=exp(x),
        "Temperature"=x,
        "ln_margsd"=exp(x),
        "rho"=x, 
        "logitbound_h"=exp(x)/(1+exp(x)),
        "h"=x,
        "ln_r"=exp(x),
        "ln_Fmsy"=exp(x))
        }
        
        parNameKey <- c(
        "K"="K", "M"="M", "Winfinity"="W_inf",
        "Loo"="L_inf", "tmax"="A_max", "tm"="A_mat", 
        "Lm"="L_mat","Temperature"="T", "ln_margsd"="Sigma_R",
        "rho"="rho", "logitbound_h"="h", "h"="h", "ln_r"="r","ln_Fmsy"="Fmsy")
        
        xy1 <- xy
        xy1sc <- sapply(1:length(xy),function(x){scaleparam(xy[x],names(xy)[x])})
        names(xy1sc) <- parNameKey[names(xy1sc)]
        
      }
      
      for (aI in 1:2) {
        Text = switch(params[rowI, aI], Loo = "Asymptotic length (L_inf)", 
                      K = "Relative growth rate (K)", Winfinity = "Asymptotic mass (W_inf)", 
                      tmax = "Maximum age (A_max)", tm = "Age at maturity (A_mat)", 
                      M = "Mortality rate (M)", Lm = "Length at maturity (L_mat)", 
                      Temperature = "Average temperature", 
                      ln_var = "Conditional recruitment variance", 
                      rho = "Recruitment autocorrelation (rho)", 
                      ln_MASPS = "Maximum annual spawners per spawner", 
                      ln_margsd = "SD of recruitment (Sigma_R)", 
                      h = "Steepness (h)", logitbound_h = "Steepness (h)", 
                      ln_Fmsy_over_M = "Ratio of F_msy and M", 
                      ln_Fmsy = "Fishing mortality rate at MSY", 
                      ln_r = "Intrinsic growth rate (r)", ln_G = "Generation time", 
                      r = "Intrinsic growth rate (r)", G = "Generation time", 
                      params[rowI, aI])
        mtext(side = aI, text = Text, line = 1.5)
      }
    }
    # Add mean prediction text for finest taxon
    legend("topright",
           legend=paste(names(xy1sc),signif(xy1sc,2),sep=" = "),
           #bty="n",
           text.col=lcol[1],
           title="Pre.",
           box.col=legend.box.col,
           bg=legend.bg,
           inset=0.02
    )
    
    # Add observed parameter values to plots
    if(!is.null(params.obs)){
      xy.obs <- params.obs_i[params[rowI, ]]
      points(x=xy.obs[1],y=xy.obs[2],pch=16)
      arrows(x0=xy1[1],y0=xy1[2],x1=xy.obs[1],y1=xy.obs[2],length=0.05)
      
      # Add mean prediction text for observed values
      xy.obssc <- sapply(1:length(xy.obs),function(x){scaleparam(xy.obs[x],names(xy.obs)[x])})
      names(xy.obssc) <- parNameKey[names(xy.obssc)]
      legend("right",
             legend=paste(names(xy.obssc),signif(xy.obssc,2),sep=" = "),
             #bty="n",
             text.col="black",
             title="Obs.",
             box.col=legend.box.col,
             bg=legend.bg,
             inset=0.02
      )
    }
    
    if (rowI %in% legendnum) {
      Legend = Taxa
      if (drop_pred == TRUE) 
        Legend = gsub(x = Taxa, pattern = "_predictive", 
                      replacement = "")
        Legend = lapply(strsplit(Legend,"_"),function(x){tail(x,1)}) # Just show finest taxon for each level
        #str_extract(Legend,"[A-Za-z]+$")  
        legend("topleft", legend = Legend, fill = lcol, 
             bty = "n")
    }
  }
  names(Pred_taxa) <- unique(Taxa)
  
  return(invisible(Pred_taxa))
}
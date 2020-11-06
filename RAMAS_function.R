
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###                                                                                       #
### This function runs RAMAS                                                              #
###                                                                                       #
### Inputs: asc habitat suitability map, parameter estimates, and working directory       #
### The working directory should contain the habitat suitability map                      #
### In this working directory, the function will create a new folder for the run          #
### Which, at the end of the run, will contain all the input and output RAMAS files.      #
### At the end of the run, working directory is set back to where it was before the run   #
###                                                                                       #
### Once the PVA is run, additional code at the end of this script                        #
###           can help extract information from RAMAS output files                        #
###           (Final population size, expected minimum abundance)                         #
###                                                                                       #
###                                                                                       #
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


### This function is valid if:
# Female only 
# No habitat dynamics
# overall initial abundance remain the same across different maps; it is the number of available cells that changes, and so do, therefore, local populations densities

############ Detailed inputs of the RAMAS() function 

#### Working environment 
## run = name of run/scenario/settings (name of the new folder, created within the function, containing all the RAMAS files)
## wd = working directory, where the habitat map is
## folder_path = folder where the RAMAS output files should be (i.e. "wd\\run") 
##               NOTE THAT THIS IS IN BATCH CODE, i.e. use \\ to separate directories,  avoid spaces, and start from the root "C:\\..."
## RAMAS_path <- path to the software 

#### Simulation length and replications
## n_repli <- how many replicate simulations? 1000 is a minimum, 10000 is advised
## time_step <- time step in years
## duration <- how many time steps do we want to run the simulation for?

#### Spatial program
# map <- map file name

# cell_length <- cell length in km
# HS <- HS threshold (Baseline should be representing 85% presence points on the habitat suitability map)
# ND <- Neighbourhood distance (in terms of cell length!)

####  coefficients for  the dispersal function between patches i and j: Mij = a*exp(-Dij^c/b) if Dij < Dmax, Mij = 0 if Dij>Dmax
# a <- y-intercept of the dispersal functino
# b <- average dipersal (km)
# c <- shape of curve
# Dmax <- maximum dispersal distance (km)
# CV_disp <- coefficient of variation for dispersal
# corr <- c(a, b, c) where a, b and c and coefficients for the environmental stochasticity correlation between patches  

#### Demographic parameters
# Rmax <- Rmax (maximum rate of population growth)
#         this value is used in some cases of density dependence (Contest, Scramble) as a deterministic growth rate of the population at certain levels of population abundances (depending on the carrying capacity)
# N0 <- How many FEMALE individuals are currently in the whole region?

#### carrying capacity function's terms
#   max_density <- maximum density of individuals per cell (i.e. carrying capacity)
#   min_patch_size <-  minimum size for a patch to be viable in terms of number of cells. default = 0
#   K_sd <- standard deviation of the carrying capacity

#### Stages
# breeding <- vector of length the number of stages, indicating the proportion of each stage breeding 
#             ordered by stage, from youngest to oldest
# rel_disp <- vector of the relative dispersal per stage 
# bool_DD <- TRUE/FALSE boolean vector if each stage is basis for density dependance
# weight <- average weight per stage 

# stageM =  stage matrix
#          needs to be a matrix object of characters and have 9 decimals
# sdM <- standard deviation matrix
#        same format as the stage matrix

#### Density dependence 
# DD_type = density dependance type : "BH" (contest = Beverton-Holt), "EXP" (exponential = no DD), "CE" (ceiling), "LO" (Scramble = Logistic)
# DD_vital_rates = what vital rates are affected by density dependence?
# stochasticity <- Environmental stochasticity - which variables (Fecundities, Survival, Karrying kapacity) should be correlated?
#                 "1 (F, S, K correlated)"             (default)
#                 "2 (F and S correlated)"
#                 "3 (F, S, K uncorrelated)"

RAMAS <- function(run, wd, folder_path, RAMAS_path,
                  time_step, duration, n_repli, 
                  map, cell_length, 
                  HS, 
                  a, b, c, Dmax, CV_disp, corr = c("0.50000", "10.00000", "1.00000"), 
                  ND, Rmax, N0, 
                  min_patch_size = 0, max_density, K_sd, 
                  stageM, sdM, 
                  rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
                  DD_type, DD_vital_rates = "all vital rates", 
                  stochasticity =  "1 (F, S, K correlated)"){
  
  # create new folder where all the RAMAS files will be stored for this run
  
  dir.create(paste0(wd,run))
  setwd(paste0(wd, run))
  
  
  # get habitat suitability map to calculate the number of habitable cells present in the landscape
  HS_map <- raster(paste0("../", map))
  # number of habitable cells
  ncell <- sum(getValues(HS_map)>= HS, na.rm = TRUE)
  
  
  #### Write base mp file ####
  sink(paste0(run, "_", "base.mp"))
  
  cat("Metapopulation input file (5.1) map=y\n")
  cat(paste0(run,"\n\n\n\n\n"))
  cat(paste0(n_repli,"\n"))     # number of replications
  cat(paste0(duration,"\n"))     # duration 
  cat("TRUE\n")
  cat(paste0(nrow(stageM), " FALSE\n\n\n")) 
  cat("Local\n\n")            # Catastrophes
  cat("not spread\n") 
  cat("0.000\n") 
  cat("0.000,0.000,0.000,0.000\n\n\n")
  cat("Local\n\n")
  cat("not spread\n") 
  cat("0.000\n")
  cat("0.000,0.000,0.000,0.000\n")
  cat("False,Zero\n")
  cat(paste(DD_vital_rates,"\n"))     # density dependence affects :
  cat("Lognormal,0\n")
  cat(CV_disp, "\n")           # CV for dispersal
  cat("count in total\n")
  cat(paste0(stochasticity,"\n"))   # stochasticity
  cat("No\n")
  cat("SelectedStages\n")      # Density dependence is based on the abundance of :
  cat("No\n")
  cat(paste0(DD_type,"\n\n"))
  cat(paste0(time_step,"\n"))
  cat("years\n")
  cat("OnlyFemale\n")          # sex structure, the model includes :
  cat("1\n")                   
  cat("Monogamous\n")
  cat("2.0\n")
  cat("2.0\n")
  cat("0.000\n")
  cat("0\n")
  cat(paste0("Pop 1,0.000,0.000,0,EX,", Rmax, ",0,", K_sd, ",0.0,,0.0,1,0,TRUE,1,1,1,0.0,1,0,1,0.0,0,0,1.0,1,1,\n"))
  cat("Migration\n")
  cat("FALSE\n")
  cat("0.000,0.00000,0.00000,0.00000\n")
  cat(" 0,\n")
  cat("Correlation\n")
  cat("FALSE\n")
  cat("0.000,0.00000,0.00000\n")
  cat(" 0,\n")
  cat("1 type(s) of stage matrix\n")
  cat("default\n")
  cat("1.000000\n")
  cat("1.000000\n")
  cat("0\n")
  for (i in 1:nrow(stageM)){
    for (j in 1:ncol(stageM)){
      cat(stageM[i, j],"") 
    }
    cat("\n")
  }
  cat("1 type(s) of st.dev. matrix\n")
  cat("default\n")
  for (i in 1:nrow(sdM)){
    for (j in 1:ncol(sdM)){
      cat(sdM[i, j],"") 
    }
    cat("\n")
  }
  cat("Constraints Matrix\n")
  for (i in 1:nrow(stageM)){
    for (j in 1:ncol(stageM)){
      if (i ==1){
        cat("0.000000","")
      }
      else{
        cat("1.000000", "")
      }
    }
    cat("\n")
  }
  for (i in 1:length(rel_disp)){
    cat(rel_disp[i],"")
  }
  cat("\n")
  for (i in 1:(2*nrow(stageM)+2)){
    for (j in 1:ncol(stageM)){
      cat("1.000000","") 
    }
    cat("\n")
  }
  for (i in 1:nrow(stageM)){
    cat("0 ")
  }
  cat("\n")                  # initial abundance per stage (in the base mp file, set zero, as this info is derived from the spatial data program)
  for (i in 1:nrow(stageM)){
    cat(paste("Stage", i, sep = "_"), "\n")
    cat(weight[i], "\n")          # average weight per stage - usually = 1
    cat("FALSE\n")                # exclude this stage from model
    cat(bool_DD[i], "\n")         # Basis for DD
    cat(breeding[i], "\n", sep = "")        # Proportion breeding
  }
  cat("0 (pop mgmnt)\n")
  cat(" 0.00000000000000E+0000\n")
  cat(" 0.00000000000000E+0000\n")
  cat("1\n")
  cat("-End of file-")
  
  sink()
  
  base_mp_path <- paste0(run, "_", "base.mp")
  
  #### Write ptc file ####
  
  ### determine number of habitable cells (to initialize local population sizes)
  
  
  
  ptc <- file(paste0(run, ".ptc"))      # where run = number of run 
  
  writeLines(c("Landscape input file (4.1) map=y",
               paste0(run, "\n\n\n\n"),    # title
               cell_length,                      # cell length in km 
               paste0("[",run,"]", "\n"),      # Habitat function: name of map between []
               HS,                             # habitat suitability threshold
               ND,                             # neighbourhood distance
               "Blue,False",                   
               "2",                            # decimals
               paste0("((ths)*", max_density, "/(noc))*thr((noc),", min_patch_size, ")"),    # Carrying capacity function                                      # carrying capacity function 
               Rmax,                                       
               paste0("(", N0, "/", ncell, ")*(noc)"),           # initial population per patch
               "1",                                        # relative fecundity 
               "1\n",                                     # relative survival 
               base_mp_path,                               # "other data from"
               "0",                                        # catastrophe 1 probability
               "1",                                        # catastrophe 1 multiplier
               "0",                                        # catastrophe 2 probability
               "1",                                        # catastrophe 2 probability
               "No\n",
               "Edge to edge",                 # distances between patches based on
               1,
               run,
               paste0("..\\", map),
               "ARC/INFO,ConstantMap",
               "Blue",
               "2845",                         # number of columns in map file
               paste0(",0.000,0.000,,", DD_type, ",,,",K_sd,",0.0,,0.0,1,0,TRUE,1,1,1,0.0,1,0,1,0,0,0,1.0,"),
               "Migration",
               "TRUE",                         # Dispersal distance function 
               paste(a, b, c, Dmax, sep = ","),
               "Correlation",
               "TRUE",                         # Correlation function 
               paste(corr[1], corr[2], corr[3], sep = ","),
               "-End of file-"), ptc)
  close(ptc)
  
  #### Spatial program batch run ####
  
  batch_ptc <- file(paste0(run, "Patch.bat"))
  writeLines(c(paste("CD ", folder_path),
               paste0("START /WAIT ", RAMAS_path, "\\Patch.exe ", paste0(run, ".ptc"), " /RUN=YES /TEX"),
               "CD .."), batch_ptc)
  close(batch_ptc)
  
  system(paste0(run, "Patch.bat")) #runs the batch file i.e. runs the spatial program
  
  
  #### Write pdy file #### 
  ### requires that the spatial data program has been run
  # there are no changes in habitat over time
  # we only run the Habitats Dynamics subprogram to automate the creation of mp file
  
  pdy <- file(paste0(run, ".pdy"))
  writeLines(c("Habitat Dynamics (version version.dll)", 
               paste0(run, "\n\n\n\n"),
               paste0(run, ".mp"), #edit a new .mp file for every different PVA run
               "POP", 
               "POP", 
               "POP",
               "2",
               paste0(run, ".ptc"), 
               "1", 
               "linear", 
               "linear", 
               "linear",
               paste0(run, ".ptc"), 
               "2", 
               "linear", 
               "linear", 
               "linear"), pdy)
  close(pdy)
  
  #### HabDyn program batch run ####
  ## write mp file  by running batch for habitat dynamics program 
  batch_pdy <- file(paste0(run, "HabDyn.bat"))
  writeLines(c(paste0("CD ", folder_path),
               paste0("START /WAIT ", RAMAS_path, "\\HabDyn.exe ", paste0(run, ".pdy"), " /RUN=YES /TEX"),
               "CD .."), batch_pdy)
  close(batch_pdy)
  
  system(paste0(run,"HabDyn.bat"))
  
  #### Metapop program batch run ####
  ## run metapopulation model -- write and execute batch file -- save outputs in new separate folder created in batch 
  batch_mp <- file(paste0(run, "Metapop.bat"))
  writeLines(c(paste0("CD ", folder_path),
               paste0("START /WAIT ", RAMAS_path, "\\Metapop.exe ", paste0(run, ".mp"), " /RUN=YES /TEX"),
               "CD .."), batch_mp)
  close(batch_mp)
  
  system(paste0(run,"Metapop.bat"))
  
  
  setwd("../../")
}


#### Extracting output data: EMA, final pop sizes ####


PVA_output <- function(wd, run, duration){
  # duration = number of time steps 
  setwd(paste0(wd, run))
  
  # create recipient dataframe
  df <- data.frame(Param_value = character(),
                   EMA = numeric(), 
                   PopF = numeric(), 
                   stringsAsFactors = F)
  
  # retrieve EMA
  IntExtRisk <- file("IntExtRisk.txt", 'r') #read ramas output file
  EMA <- readLines(IntExtRisk, n=14) 
  close(IntExtRisk)
  EMA <- EMA[14] # EMA appears at 14th line of that file
  EMA <- strsplit(EMA, split = " ")
  EMA <- as.numeric(EMA[[1]][5])
  
  # retrieve Final population abundances after the [duration] years simulation
  Abund <- file("Abund.txt", 'r') #read ramas output file
  PopF <- readLines(Abund)
  close(Abund)
  PopF <- PopF[duration + 16] # Final population appears on the (duration + 16)th line of that file
  PopF <- strsplit(PopF, split = "    ")
  PopF <- as.numeric(PopF[[1]][5])
  
  # put together
  df[1, ] <- c(run, EMA, PopF)
  
  return(df)
  
  
}




#### Sensitivity analysis ####
#### Or any analysis comparing different PVA results

PVA_SA <- function(spp,
                   run, wd, folder_path, RAMAS_path,
                   time_step, duration, n_repli, 
                   map, cell_length, 
                   HS, ncell,
                   a, b, c, Dmax, CV_disp, corr, 
                   ND, Rmax, N0, 
                   min_patch_size = 0, max_density, K_sd, 
                   stageM, sdM, 
                   rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
                   DD_type, DD_vital_rates = "all vital rates", 
                   stochasticity =  "1 (F, S, K correlated)"){
  
  #### EXECUTE RAMAS ####
  
  
  
  RAMAS(run, wd, folder_path, RAMAS_path,
        time_step, duration, n_repli, 
        map, cell_length, 
        HS, ncell,
        a, b, c, Dmax, CV_disp, corr, 
        ND, Rmax, N0, 
        min_patch_size = 0, max_density, K_sd, 
        stageM, sdM, 
        rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
        DD_type, DD_vital_rates = "all vital rates", 
        stochasticity =  "1 (F, S, K correlated)")
  
  
  #### Sensitivity analysis ####
  
  # Rmax +- 0.005 (otherwise too big a variation)
  for (p in c(Rmax-0.005, Rmax+0.005)){
    if (p == Rmax-0.005){
      run <- paste0("Rmax_", "min")
    }
    else if (p == Rmax+0.005) {
      run <- paste0("Rmax_", "max")
    }
    folder_path <- paste0("C:\\Users\\LOUISEMARIEL\\Dropbox\\QAEco_GreaterHunter\\", spp, "\\", run)
    RAMAS(run, wd, folder_path, RAMAS_path,
          time_step, duration, n_repli, 
          map, cell_length, 
          HS, ncell,
          a, b, c, Dmax, CV_disp, corr, 
          ND, p, N0, 
          min_patch_size = 0, max_density, K_sd, 
          stageM, sdM, 
          rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
          DD_type, DD_vital_rates = "all vital rates", 
          stochasticity =  "1 (F, S, K correlated)")
  }
  
  # habitat suitability threshold
  for (p in c(HS_80, HS_90)){
    if (p == HS_80){
      run <- paste0("HS_", "80")
    }
    else if (p == HS_90) {
      run <- paste0("HS_", "90")
    }
    folder_path <- paste0("C:\\Users\\LOUISEMARIEL\\Dropbox\\QAEco_GreaterHunter\\", spp, "\\", run)
    RAMAS(run, wd, folder_path, RAMAS_path,
          time_step, duration, n_repli, 
          map, cell_length, 
          p, ncell,
          a, b, c, Dmax, CV_disp, corr, 
          ND, Rmax, N0, 
          min_patch_size = 0, max_density, K_sd, 
          stageM, sdM, 
          rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
          DD_type, DD_vital_rates = "all vital rates", 
          stochasticity =  "1 (F, S, K correlated)")
  }
  
  
  # Neighbourhood distance 
  for (p in c(ND*(1-0.2), ND*(1+0.2))){
    if (p == ND*(1-0.2)){
      run <- paste0("ND_", "min")
    }
    else if (p == ND*(1+0.2)) {
      run <- paste0("ND_", "max")
    }
    folder_path <- paste0("C:\\Users\\LOUISEMARIEL\\Dropbox\\QAEco_GreaterHunter\\", spp, "\\", run)
    RAMAS(run, wd, folder_path, RAMAS_path,
          time_step, duration, n_repli, 
          map, cell_length, 
          HS, ncell,
          a, b, c, Dmax, CV_disp, corr, 
          p, Rmax, N0, 
          min_patch_size = 0, max_density, K_sd, 
          stageM, sdM, 
          rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
          DD_type, DD_vital_rates = "all vital rates", 
          stochasticity =  "1 (F, S, K correlated)")
  }
  
  # initial population size
  for (p in c(N0*(1-0.2), N0*(1+0.2))){
    if (p == N0*(1-0.2)){
      run <- paste0("N0_", "min")
    }
    else if (p == N0*(1+0.2)) {
      run <- paste0("N0_", "max")
    }
    folder_path <- paste0("C:\\Users\\LOUISEMARIEL\\Dropbox\\QAEco_GreaterHunter\\", spp, "\\", run)
    RAMAS(run, wd, folder_path, RAMAS_path,
          time_step, duration, n_repli, 
          map, cell_length, 
          HS, ncell,
          a, b, c, Dmax, CV_disp, corr, 
          ND, Rmax, p, 
          min_patch_size = 0, max_density, K_sd, 
          stageM, sdM, 
          rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
          DD_type, DD_vital_rates = "all vital rates", 
          stochasticity =  "1 (F, S, K correlated)")
  }
  
  ## dispersal average b 
  for (p in c(b*(1-0.2), b*(1+0.2))){
    if (p == b*(1-0.2)){
      run <- paste0("b_", "min")
    }
    else if (p == b*(1+0.2)) {
      run <- paste0("b_", "max")
    }
    folder_path <- paste0("C:\\Users\\LOUISEMARIEL\\Dropbox\\QAEco_GreaterHunter\\", spp, "\\", run)
    RAMAS(run, wd, folder_path, RAMAS_path,
          time_step, duration, n_repli, 
          map, cell_length, 
          HS, ncell,
          a, p, c, Dmax, CV_disp, corr, 
          ND, Rmax, N0, 
          min_patch_size = 0, max_density, K_sd, 
          stageM, sdM, 
          rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
          DD_type, DD_vital_rates = "all vital rates", 
          stochasticity =  "1 (F, S, K correlated)")
  }
  
  # dispersal max limit
  for (p in c(Dmax*(1-0.2), Dmax*(1+0.2))){
    if (p == Dmax*(1-0.2)){
      run <- paste0("Dmax_", "min")
    }
    else if (p == Dmax*(1+0.2)) {
      run <- paste0("Dmax_", "max")
    }
    folder_path <- paste0("D:\\Users\\LOUISEMARIEL\\Documents\\QAEco\\", spp, "\\", run)
    RAMAS(run, wd, folder_path, RAMAS_path,
          time_step, duration, n_repli, 
          map, cell_length, 
          HS, ncell,
          a, b, c, p, CV_disp, corr, 
          ND, Rmax, N0, 
          min_patch_size = 0, max_density, K_sd, 
          stageM, sdM, 
          rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
          DD_type, DD_vital_rates = "all vital rates", 
          stochasticity =  "1 (F, S, K correlated)")
    
  }
  
  # carrying capacity / maximum density estimate 
  for (p in c(max_density*(1-0.2), max_density*(1+0.2))){
    if (p == max_density*(1-0.2)){
      run <- paste0("max_density_", "min")
    }
    else if (p == max_density*(1+0.2)) {
      run <- paste0("max_density_", "max")
    }
    folder_path <- paste0("D:\\Users\\LOUISEMARIEL\\Documents\\QAEco\\", spp, "\\", run)
    RAMAS(run, wd, folder_path, RAMAS_path,
          time_step, duration, n_repli, 
          map, cell_length, 
          HS, ncell,
          a, b, c, p, CV_disp, corr, 
          ND, Rmax, N0, 
          min_patch_size = 0, max_density, K_sd, 
          stageM, sdM, 
          rel_disp, breeding, bool_DD = rep("TRUE", nrow(stageM)), weight = rep("1.0", nrow(stageM)), 
          DD_type, DD_vital_rates = "all vital rates", 
          stochasticity =  "1 (F, S, K correlated)")
    
  }
  
  
  
  ##### EXTRACTING THE  OUTPUTS DATA ######
  
  SA_outputs <- rbind(PVA_output(wd, run = paste0(spp, "_baseline")), 
                      PVA_output(wd, run = "Dmax_min"), 
                      PVA_output(wd, run = "Dmax_max"), 
                      PVA_output(wd, run = "Rmax_max"), 
                      PVA_output(wd, run = "Rmax_min"), 
                      PVA_output(wd, run = "HS_max"), 
                      PVA_output(wd, run = "HS_min"), 
                      PVA_output(wd, run = "b_max"), 
                      PVA_output(wd, run = "b_min"), 
                      PVA_output(wd, run = "ND_max"), 
                      PVA_output(wd, run = "ND_min"), 
                      PVA_output(wd, run = "N0_max"), 
                      PVA_output(wd, run = "N0_min"))
  
  
  library(tidyr)
  SA_output <- SA_outputs %>%
    separate(col = Param_value, into = c("parameter", "value"), sep = "_")
  
  
  #### % change from the baseline ####
  
  library(dplyr)
  
  SA_output$EMA <- as.numeric(SA_output$EMA)
  SA_output$PopF <- as.numeric(SA_output$PopF)
  
  EMA_base <- SA_output[1, 3]
  PopF_base <- SA_output[1, 4]
  
  SA_output <- SA_output %>%
    mutate(EMA_perc_change = (EMA - EMA_base) / EMA_base * 100) %>%
    mutate(PopF_perc_change = (PopF - PopF_base) / PopF_base * 100)
  
  
  #### bar plot ####
  library(ggplot2)
  
  #final pop sizes
  PopF_SA <- ggplot(data=SA_output[-1,], aes(x=parameter, y=PopF_perc_change)) +
    geom_bar(stat="identity",  fill = "steelblue", color = "black", width = 0.5) + 
    geom_hline(yintercept = 0, color = "red") +
    ggtitle(label = "Sensitivity of final population sizes to parameter variation",   subtitle = paste("for the", spp, "and original final abundances of", PopF_base, sep = " ")) + 
    xlab("Parameter") + 
    ylab("% change from original estimate")+ 
    ylim(-80, 80)
  
  #EMA
  EMA_SA <- ggplot(data=SA_output[-1,], aes(x=parameter, y=EMA_perc_change)) +
    geom_bar(stat="identity",  fill = "steelblue",  color = "black", width = 0.5) + 
    geom_hline(yintercept = 0, color = "red")+ 
    ggtitle(label = "Sensitivity of expected minimum abundance to parameter variation",   subtitle = paste("for the", spp, "and original EMA of", EMA_base, sep = " ")) +
    xlab("Parameter") + 
    ylab("% change from original estimate")+ 
    ylim(-80, 80)
  
  return(EMA_SA)
  return(PopF_SA)
  
}




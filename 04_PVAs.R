##### PVA_1 Squirrel Glider #####
##### BASELINE
library(raster)
wd <- paste0("~")
run <- paste0("BASELINE")  # this will be used, among other things, to name a separate folder for a particular PVA run
map <- print(paste0("SQG_BASELINE.asc"))
cell_length = 0.1
HS <- 0.405
#number of replications
n_repli <- 1000
# 1 time step in years
time_step <- 1
# how many time steps do you want to run the PVA for ?
duration <- 100
# path to folder for batch run 
folder_path <- paste0("~", run,"")
RAMAS_path <- "RAMAS\\RAMASGIS\\RAMASGIS"
# total number of cells in the area 
tc <- 454728
# stage matrix
stageM = matrix(c("0.000000000" , "0.970000000" , "1.050000000" ,      
                  "0.320000000" , "0.000000000" , "0.000000000" ,
                  "0.000000000" , "0.390000000" , "0.630000000" ),
                nrow=3,  
                ncol=3,  
                byrow = TRUE)           # fill matrix by rows 

# standard deviation matrix
sdM <- matrix(c("0.000000000", "0.120000000", "0.100000000",
                "0.120000000", "0.000000000", "0.000000000",
                "0.000000000", "0.060000000", "0.060000000"),
              nrow = 3, 
              ncol = 3,
              byrow = TRUE)

# proportion of each stage breeding 
breeding <- c("0.00", "0.95", "0.95")
max_density <- 0.228
min_patch_size <- 20
# K standard deviation
K_sd <- 0.1
# CV for dispersal 
CV_disp <- "0.250000"
# relative dispersal per stage 
rel_disp <- c("0.000000", "1.000000", "0.000000 ")
# Habitat suitability threshold 
# Neighbourhood distance
ND <- 1.50000
# Dispersal coefficients
a <- 1.000
b <- 0.75000
c <- 0.75000
Dmax <- 3.00000
# correlation coefficients 
corr <- c("1.00000", "3.00000", "1.00000")
# Rmax
Rmax <- 1.005
# Initial estimated abundance function
N0 <- 10000

# TRUE/FALSE boolean vector if each stage is basis for density dependance
# where the order of elements is ordered by stage, from youngest to oldest
bool_DD <- c("FALSE", "TRUE", "TRUE")

# average weight per stage 
weight <- c("1.0", "1.0", "1.0")

## density dependance type : BH = "contest", EXP = "exponential", CE = "ceiling"
DD_type <- "BH"

stochasticity <- "1 (F, S, K correlated)"
#ncell <- sum(getValues(map)>HS, na.rm = TRUE)
#### RUN ##

source("RAMAS_function.R")

RAMAS(run, wd, folder_path, RAMAS_path,
      time_step, duration, n_repli, 
      map, cell_length, 
      HS, 
      a, b, c, Dmax, CV_disp, corr, 
      ND, Rmax, N0, 
      min_patch_size, max_density, K_sd, 
      stageM, sdM, 
      rel_disp, breeding, bool_DD, weight, 
      DD_type, DD_vital_rates = "all vital rates", 
      stochasticity =  "1 (F, S, K correlated)")


##### PVA_2 Powerful Owl #####
#### BASELINE 

library(raster)
wd <- paste0("~")
run <- paste0("BASELINE")  # this will be used, among other things, to name a separate folder for a particular PVA run
map <- print(paste0("PO_BASELINE.asc"))

#### DEFINE PARAMETERS OF INTEREST 

cell_length = 0.1

HS <- 0.265

#number of replications
n_repli <- 1000

# 1 time step in years
time_step <- 1

# how many time steps do you want to run the PVA for ?
duration <- 100

# path to folder for batch run 

folder_path <- paste0("~", run,"")

RAMAS_path <- "C:\\RAMAS\\RAMASGIS\\RAMASGIS"

# total number of cells in the area 
tc <- 454728

# stage matrix
stageM = matrix(c("0.000000", "0.450000", "0.450000",      
                  "0.670000", "0.000000", "0.000000",
                  "0.000000", "0.670000", "0.890000"),
                nrow=3,  
                ncol=3,  
                byrow = TRUE)        # fill matrix by rows 

# standard deviation matrix
sdM <- matrix(c("0.000000000", "0.100000000", "0.100000000",
                "0.033000000", "0.000000000", "0.000000000",
                "0.000000000", "0.033000000", "0.400000000"),
              nrow = 3, 
              ncol = 3,
              byrow = TRUE)


# proportion of each stage breeding 
breeding <- c(0, 0.32, 0.63)

max_density <- 0.0035
min_patch_size <- 200

# K standard deviation
K_sd <- 0.1

# CV for dispersal 
CV_disp <- "0.100000"


# relative dispersal per stage 
rel_disp <- c("0.000000", "1.000000", "0.100000")

# Neighbourhood distance
ND <- 5.34000

# Dispersal coefficients
a <- 1.000
b <- 20.00000
c <- 0.70000
Dmax <- 50.00000

# correlation coefficients 
corr <- c("1.000", "50.00000", "1.50000")

# Rmax
Rmax <- 1.158

# Initial estimated abundance function
N0 <- 300


## TRUE/FALSE boolean vector if each stage is basis for density dependance
# where the order of elements is ordered by stage, from youngest to oldest
bool_DD <- c("FALSE", "TRUE", "TRUE")

# average weight per stage 
weight <- c("1.0", "1.0", "1.0")

## density dependance type : BH = "contest", EXP = "exponential", CE = "ceiling"
DD_type <- "BH"
DD_vital_rates <- "all vital rates"

stochasticity <- "1 (F, S, K correlated)"


source("RAMAS_function.R")

RAMAS(run, wd, folder_path, RAMAS_path,
      time_step, duration, n_repli, 
      map, cell_length, 
      HS, 
      a, b, c, Dmax, CV_disp, corr, 
      ND, Rmax, N0, 
      min_patch_size, max_density, K_sd, 
      stageM, sdM, 
      rel_disp, breeding, bool_DD, weight, 
      DD_type, DD_vital_rates = "all vital rates", 
      stochasticity =  "1 (F, S, K correlated)")


##### PVA_3 Northern Brown Bandicoot #####


wd <- paste0("~")
run <- paste0("BASELINE")  # this will be used, among other things, to name a separate folder for a particular PVA run
map <- print(paste0("NBB_BASELINE.asc"))

#### DEFINE PARAMETERS OF INTEREST 
cell_length = 0.1
HS <- 0.388
# number of replications
n_repli <- 1000
# how many time step
duration <- 200
# time step duration 
time_step <- 0.5
# path to folder for batch run 

folder_path <- paste0("~", run, "")

RAMAS_path <- "C:\\RAMAS\\RAMASGIS\\RAMASGIS"



# stage matrix using expert knowledge 
stageM = matrix(c("0.000000", "3.000000", "3.000000", "3.000000",      
                  "0.190000", "0.000000", "0.000000", "0.000000",
                  "0.000000", "0.500000", "0.000000", "0.000000",
                  "0.000000", "0.000000", "0.330000", "0.000000"),
                nrow=4,  
                ncol=4,  
                byrow = TRUE) 

# standard deviation matrix
sdM <- matrix(c("0.000000000", "0.500000000", "0.500000000", "0.500000000",
                "0.036100000", "0.000000000", "0.000000000", "0.000000000",
                "0.000000000", "0.095000000", "0.000000000", "0.000000000",
                "0.000000000", "0.000000000", "0.002700000", "0.000000000"),
              nrow = 4, 
              ncol = 4,
              byrow = TRUE)

# proportion of each stage breeding 
breeding <- c(0, 1, 1, 1)

# # carrying capacity function
# K <- "thr((ths*1.477),15)"
# original function, based on the southern brown bandicoot
max_density <- 0.985
min_patch_size <- 15
# K standard deviation
K_sd <- 0.2

# CV for dispersal 
CV_disp <- "0.100000"

# relative dispersal per stage 
rel_disp <- c("1.000000", "0.187500", "0.187500", "0.187500")
# rel_disp <- c("0.000000", "0.000000", "0.000000", "0.000000")

# Neighbourhood distance
ND <- 2.45

# Dispersal coefficients
a <- 0.3
b <- 0.15
c <- 2
Dmax <- 1

# correlation coefficients 
corr <- c("0.50000", "10.00000", "1.00000")

# Rmax
Rmax <- 1.5

# Initial estimated abundance function
N0 <- 10000

# TRUE/FALSE boolean vector if each stage is basis for density dependance
# where the order of elements is ordered by stage, from youngest to oldest
# density dependance affects all stages
bool_DD <- c("TRUE", "TRUE", "TRUE", "TRUE")

# average weight per stage 
weight <- c("1.0", "1.0", "1.0", "1.0")

## density dependance type : BH = "contest", EXP = "exponential", CE = "ceiling", LO = "SCRAMBLE"
DD_type <- "BH"

stochasticity <- "2 (F and S correlated)"

#### RUN ###

source("RAMAS_function.R")

RAMAS(run, wd, folder_path, RAMAS_path,
      time_step, duration, n_repli, 
      map, cell_length, 
      HS, 
      a, b, c, Dmax, CV_disp, corr, 
      ND, Rmax, N0, 
      min_patch_size, max_density, K_sd, 
      stageM, sdM, 
      rel_disp, breeding, bool_DD, weight, 
      DD_type, DD_vital_rates = "all vital rates", 
      stochasticity =  "1 (F, S, K correlated)")


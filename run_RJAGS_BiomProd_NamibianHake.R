# ---------------------------------------------------------------
# JAGS code for Dynamic Biomass Production Model
# (Model with Catches and AI and 5 years projection)
# V. November 2023
# E. Rivot
# ---------------------------------------------------------------

graphics.off()
rm(list = ls())
setwd("C:/Users/monau/Desktop/M2-S1/BAYES/_TD_BiomProd/_TD_BiomProd")
# ---------------------------------------------------------------

# Check the path to JAGS.exe

# Option
# If rjags can not connect directly to the JAGS.exe installed on your computer
# You will see that in the information message after loading the rjags library
# CAUTION : the path C:\\ is to be adapted to redirect where your JAGS.exe is installed
# the name JAGS-4.2.0 should also be changed if another version of JAGS is installed on your computer

# Sys.setenv(JAGS_HOME = "C:\\prog_td\\JAGS\\JAGS-4.2.0\\")

# ---------------------------------------------------------------

library(rjags)
library(coda)
library(MASS)


# Load the model
# ----------------------------------------

# Load the model, written in a .txt file
model <- "C:/Users/monau/Desktop/M2-S1/BAYES/_TD_BiomProd/_TD_BiomProd/model/modelJAGS_BiomProd_PriorParam.txt"

# Write the model specification into a virtual text file with name "model", to be found by JAGS
# The file "model" will be used by functions "jags.model"
source(model)
model <- textConnection(model_JAGS)


# Load and format data
# -----------------------------------------

# load data

data <- read.table(file = "C:/Users/monau/Desktop/M2-S1/BAYES/_TD_BiomProd/_TD_BiomProd/data/data_BiomProd.txt", header = TRUE, sep = "") 

n_obs <- length(data$Year)
n_proj <- 5
Year <- data$Year

# Vector of biomass - Needed to build the equilibrium curve
B_e <- seq(from = 0, to = 1500, by = 50)
n_equi <- length(B_e)

# Format data as a list to be read in JAGS

data <- list(
    "I_obs" = data$I,
    "C_obs" = data$C,
    "n_obs" = n_obs,
		"B_e" = B_e,
		"n_equi" = n_equi#,"Year"=Year",n_proj" = n_proj
  )

monitor <- c( "B",
              "h",
              "C",
              "I","r",
              "r_p",
              "K",
              "K_p",
              "q",
              "C_MSY",
              "C_MSY_p",
              "B_MSY",
              "B_MSY_p",
              "h_MSY", 
              "C_e"#"Over_C","I_pred","D","Over_h","risk","sigma2p","C_pred"
)
# MCMC options
# ----------------------------------------------------------------------

n.chains = 3

# Adapt the MCMC samplers with n.adapt iterations

n.adapt = 10000

# Iteration after adapting

n.burnin <- 1000
n.stored <- 1000
n.thin <- 10
n.iter <- n.stored*n.thin

# MANAGING NUMBER of ITERATIONS
# (works for both para = T and para = F)
# total number of REALIZED iterations = n.iter
# total number of STORED iterations = n.iter/n.thin



# Run the model to save MCMC results
# ---------------------------------------------------------------------

# Variable to store

monitor <- c( "B",
              "h",
              "C",
              "I","r",
		"r_p",
		"K",
		"K_p",
		"q",
		"C_MSY",
		"C_MSY_p",
		"B_MSY",
		"B_MSY_p",
		"h_MSY", 
	  "C_e"#"Over_C","I_pred","D","Over_h","risk","sigma2p","C_pred"
		)

# Compile the model, create a model object with the name "model.compiled.jags" and adapt the MCMC samplers
# with n.adapt iterations

print("adapting phase")
model <- "C:/Users/monau/Desktop/M2-S1/BAYES/_TD_BiomProd/_TD_BiomProd/model/modelJAGS_BiomProd_PriorParam.txt"
source(model)
model <- textConnection(model_JAGS)
model.compiled <- jags.model(file = model, data=data, n.chain=n.chains, n.adapt=n.adapt)

# Iteration after adapting

# Start to compute CPU time
ptm <- proc.time()
?proc.time
  # Burnin period (not stored)

print("burn-in")
update(model.compiled, n.iter=n.burnin)

# Store mcmc samples

print("mcmc stored for results")
mcmc <- coda.samples(model=model.compiled,variable.names=monitor,n.iter=n.iter,thin=n.thin)

time.to.run.mcmc <- proc.time() - ptm
print(time.to.run.mcmc)

# ----------------------------------------------------------------------------------------------


# Run the model to compute DIC
# ---------------------------------------------------------------------

# Start from a compiled model that has already been updated

n.iter.DIC <- 1000
dic.pD <- dic.samples(model.compiled, n.iter.DIC, "pD") 
dic.pD		# Deviance Information Criterion

# Alternative penalization of the Deviance
# dic.popt <- dic.samples(model.compiled.jags, n.iter, "popt")
# dic.popt
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# Calculate WAIC
# see also https://gist.github.com/oliviergimenez/68ad17910a62635ff6a062f8ec34292f
# Plummer - 2017 - I have added a WAIC monitor to the dic module in JAGS 4.3.0.
# Here the component "deviance" is the deviance statistic for each observed node and "WAIC" is the corresponding WAIC penalty.
# should add the two component to compute WAIC

# Note : need the module "dic"

mcmc.waic <- jags.samples(model = model.compiled, 
                          c("WAIC","deviance"), 
                          type = "mean", 
                          n.iter = 10000,
                          n.burnin = 1000,
                          n.thin = 10)
mcmc.waic$p_waic <- mcmc.waic$WAIC
mcmc.waic$waic <- mcmc.waic$deviance + mcmc.waic$p_waic
tmp <- sapply(mcmc.waic, sum)
waic <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic
# ----------------------------------------------------------------


# --------------------------------------------
# Work with mcmc.list
# --------------------------------------------

# "mcmc" is an object of the class "mcmc.list" (see package library(coda)
# to explore, plot ... mcmc objects

is(mcmc)

# Names of the variables stored in the mcmc list

varnames(mcmc)


# ---------------------------------------------------
# Convergence diagnostics and measure of performance
# ---------------------------------------------------

# Gelman-Rubin convergence diagnostics
# Point est. should be near 1

gelman.diag(mcmc[,c("B[1]")], confidence = 0.95, transform=TRUE, autoburnin=TRUE)
gelman.diag(mcmc[,c("q")], confidence = 0.95, transform=TRUE, autoburnin=TRUE)
gelman.diag(mcmc[,varnames(mcmc)[1:10]], confidence = 0.95, transform=TRUE, autoburnin=TRUE)

# Effective size
# An estimate of the number of "independent" samples

var = "B"
var.mcmc = mcmc[,which(substr(varnames(mcmc),1,nchar(var)+1)==paste(var,"[",sep=""))]
effectiveSize(var.mcmc) 



# --------------------------------------------
# plot results
# --------------------------------------------

# source("plot_results.R")

# Use this code to run the plot_results program with no errors
plot_results <- parse(file = "plot_results.R")
for (i in seq_along(plot_results)) {
  tryCatch(eval(plot_results[[i]]), 
           error = function(e) message("Oups!  ", as.character(e)))
}

#graphics.off()

#O


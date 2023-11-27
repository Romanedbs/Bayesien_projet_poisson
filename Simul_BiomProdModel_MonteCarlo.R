# Simulation from a dynamic Biomass Production Model
# Stochastic model with multiple trajectories

# E. Rivot - November 2023

# -------------------------------------------------
rm(list=ls())
# Number of years of projections
n <- 100

# Number of trajectories
n.traj <- 10

B <- array(data=NA, dim = c(n,1))

# Parameters for the Schaeffer model
K <- 100
C_MSY <- 10
r <- 4*C_MSY/K
B_MSY <- K/2
h_MSY <- C_MSY/B_MSY
h_MSY

# --------------------------------------------------
# --------------------------------------------------

# biomass first year
# B[1] = alpha*K
alpha <- 0.25

# --------------------------------------------------
# --------------------------------------------------

# Definition of a precautionary threshold
B_risk <- 0.20*K

# --------------------------------------------------
# --------------------------------------------------

# Harvest scenarios
# The values in h define a scenario (Ct = ht*Bt)
# In our case, h_MSY = 0.2
h <- 0.15

# Definition of harvest rate for t = 1 to n
h <- array(data = h, dim = c(n,n.traj))

# --------------------------------------------------
# --------------------------------------------------

# Noise
# 2 different levels of "noise"

# 1 - UNCERTAINTY ON PARAMETERS VALUES
# Values of parameters randomly drawn for each trajectory
# DOES not propagate through the dynamics

 sd.r <- 1
 sd.K <- 10

r.traj <- rnorm(mean=r,sd = sd.r, n=n.traj)
K.traj <- rnorm(mean=K,sd = sd.K, n=n.traj)

# 2 - ENVIRONMENTAL STOCHASTICITY
# Propagates through the dynamics

sd.noise = 0.15
# sd.noise = 0.01
# sd.noise = 0.05
# sd.noise = 0.10



# --------------------------------------------------
# --------------------------------------------------



# Definition of variable B[t] for t =1 to n
B <- array(NA,dim=c(n,n.traj))


# -----------------------------------------
# Loop on the trajectories

for (traj in 1:n.traj)
{

# Initial Biomass
# B[1,traj] <- alpha*K + rnorm(1,mean=0,sd=sd.noise)
B[1,traj] <- alpha*K*exp(rnorm(1,mean=0,sd=sd.noise))

# Dynamic evolution of the biomass (within one trajectory)
for (t in 1:(n-1))
{
# With a normal random noise
# Bnew <- B[t,traj] + r.traj[traj]*B[t,traj]*(1-B[t,traj]/K.traj[traj]) - h[t,traj]*B[t,traj] + rnorm(1,mean=0,sd=sd.noise)

# With a logNormal random noise
# Bnew <- (B[t,traj] + r.traj[traj]*B[t,traj]*(1-B[t,traj]/K.traj[traj]) - h[t,traj]*B[t,traj])*exp(rnorm(1,mean=0,sd=sd.noise))

# With LogNormal random noise + "trick" to avoid negative biomass
Bm_new <- max( (B[t,traj] + r.traj[traj]*B[t,traj]*(1-B[t,traj]/K.traj[traj]) - h[t,traj]*B[t,traj]) , 0.001*K )  
log_Bnew <- rnorm(1,mean=log(Bm_new),sd=sd.noise)
# Equivalent to log_Bnew <- log(Bm_new) + rnorm(1,mean=0,sd=sd.noise)
B[t+1,traj] <- exp(log_Bnew)
}

}
# -----------------------------------------



# Plot the biomass trajectories
# -----------------------------

# windows()
par(mfrow = c(3,1))

# Worm plot

plot(1:n, B[,1], type = "l", ylim = c(0,1.2*K), main = "Worm plot of Biomass trajectories")

for (traj in 2:n.traj) {
points(1:n, B[,traj], type = "l")
}
points(1:n,rep(K,times=n), col = "blue", type = "l", lwd = 2)
text(n/2,K,"Carrying capacity K", pos = 3, col = "blue")
points(1:n,rep(B_MSY,times=n), col = "green", type = "l", lwd = 2)
text(n/2,B_MSY,"B at MSY", pos = 1, col = "green")
points(1:n,rep(B_risk,times=n), col = "red", type = "l", lwd = 3)
text(n/2,B_risk,"Risk zone : B < 20% K", pos = 1, col="red")

# Boxplot

boxplot(as.data.frame(t(B)), main = "Marginal probability distribution of the biomass", ylim = c(0,1.2*K), outline=F)
points(1:n,rep(K,times=n), col = "blue", type = "l", lwd = 2)
text(n/2,K,"Carrying capacity K", pos = 3, col = "blue")
points(1:n,rep(B_MSY,times=n), col = "green", type = "l", lwd = 2)
text(n/2,B_MSY,"B at MSY", pos = 1, col = "green")
points(1:n,rep(B_risk,times=n), col = "red", type = "l", lwd = 3)
text(n/2,B_risk,"Risk zone : B < 20% K", pos = 1, col = "red")


# plot risk = P(B<B_risk)
# -----------------------
# windows()
# par(mfrow = c(1,1))
seuil <- B_risk
risk <- apply(  as.data.frame(B) < seuil  , MARGIN=1, FUN=mean)
plot(1:n,risk, ylim = c(0,1), main = "Risk = P(Biomass<20% K)", type = "l", col = "red")





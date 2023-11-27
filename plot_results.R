# ---------------------------------------------------------------
# Plot Results
# Dynamic Biomass Production Model
# E. Rivot - November 2023
# ---------------------------------------------------------------

library(rjags)
library(coda)
library(MASS)


# Source graphic functions
# ---------------------------------------------------------------

source("function/f.panel.cor.R")
source("function/f.panel.dens.R")
source("function/f.density.bivar.R")

# --------------------------------------------
# Work with mcmc.list
# --------------------------------------------

# "mcmc" is an object of the class "mcmc.list" (see package library(coda)
# to explore, plot ... mcmc objects

is(mcmc)

# Names of the variables stored in the mcmc list

varnames(mcmc)


# Extract the variable of interest from the mcmc list
# var.mcmc is still an mcmc list but contains only the variable of interest
# Works only for vectorial variable

var = "B"
var.mcmc = mcmc[,which(substr(varnames(mcmc),1,nchar(var)+1)==paste(var,"[",sep=""))]
varnames(var.mcmc)

# One dimensional variable can be extracted directly

r.mcmc <- mcmc[,c("r")]
B1.mcmc <- mcmc[,c("B[1]")]



# ---------------------------------------------------
# Traceplot (see also ?traceplot - coda package)
# ---------------------------------------------------

# Plot trace of ALL variables in var.mcmc-
# windows()
# plot(var.mcmc, trace = TRUE, density = FALSE)

windows()
par(mfrow=c(3,3))
traceplot(mcmc[,'K'],ylab="K")
traceplot(mcmc[,'r'],ylab="r")
traceplot(mcmc[,'C_MSY'],ylab="C_MSY")
traceplot(mcmc[,'q'],ylab="q")
traceplot(mcmc[,'sigma2p'],ylab="sigma2p")


windows()
par(mfrow=c(2,3))
traceplot(mcmc[,'B[1]'],ylab="B[1]")
traceplot(mcmc[,'B[5]'],ylab="B[5]")
traceplot(mcmc[,'B[10]'],ylab="B[10]")
traceplot(mcmc[,'B[15]'],ylab="B[15]")
traceplot(mcmc[,'B[20]'],ylab="B[20]")
traceplot(mcmc[,'B[25]'],ylab="B[25]")


# ---------------------------------------------------
# Density plot (see also ?densplot - coda package)
# ---------------------------------------------------

# Plot density of all variables in var.mcmc
# windows()
# plot(var.mcmc, trace = FALSE, density = TRUE)

windows()
par(mfrow=c(3,3))
densplot(mcmc[,'K'],xlab="K")
densplot(mcmc[,'r'],xlab="r")
densplot(mcmc[,'C_MSY'],xlab="C_MSY")
densplot(mcmc[,'q'],xlab="q")
densplot(mcmc[,'sigma2p'],xlab="sigma2p")


windows()
par(mfrow=c(2,3))
densplot(mcmc[,'B[1]'],xlab="B[1]")
densplot(mcmc[,'B[5]'],xlab="B[5]")
densplot(mcmc[,'B[10]'],xlab="B[10]")
densplot(mcmc[,'B[15]'],xlab="B[15]")
densplot(mcmc[,'B[20]'],xlab="B[20]")
densplot(mcmc[,'B[25]'],xlab="B[25]")


# --------------------------------------------
# Work with mcmc samples stored in TABLES
# --------------------------------------------


# Extract MCMC chains and store in a TABLE
# sometimes easier to manipulate than mcmc.list object
# each column = all mcmc samples for each variable

# as.matrix() does not work if mcmc contains a multidimensional variable
# but works after window is used

mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
head(mcmc.table)
dim(mcmc.table)


# Plot density from the mcmc.table

windows()
par(mfrow = c(2,3))
plot(density(mcmc.table$'K'), xlab = "K")
plot(density(mcmc.table$'r'), xlab = "r")
plot(density(mcmc.table$'C_MSY'), xlab = "C_MSY")
plot(density(mcmc.table$'q'), xlab = "q")
plot(density(mcmc.table$'sigma2p'), xlab = "sigma2p")



# ------------------------------------------------------------------
# Worm plot B
# ------------------------------------------------------------------

x = "B"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
x = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]
dim(x)

windows()
plot(1:dim(x)[2], x[1,], ylim=c(0,max(x)), type = "l", main = "Worm plot of Biomass trajectories (100 traj)", xlab = "time", ylab = "Biomass")

for (traj in 2:50) {
  points(1:dim(x)[2], x[traj,], type = "l")
}



# --------------------------------------------------------------------------
# Basic Boxplot of time series
# --------------------------------------------------------------------------
# Extract all variable "a" from the mcmc table --> put in a new table x

windows() ; par(mfrow = c(2,2))

x = "B"
title = "Biomass"

mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
x = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]

boxplot(x, outline = F, main = title, xlab = "Time", ylab = "Biomass")


x = "C"
title = "Catches"

mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
x = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]

boxplot(x, outline = F, main = title, xlab = "Time", ylab = "Catches")


x = "h"
title = "Harvest rates"

mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
x = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]

boxplot(x, outline = F, main = title, xlab = "Time", ylab = "Harvest rates")


x = "D"
title = "Depletion"

mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
x = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]

boxplot(x, outline = F, main = title, xlab = "Time", ylab = "Depletion")


# ------------------------------------------------------------------------------------
# Estimation of Risk
# ------------------------------------------------------------------------------------

var = "risk"
var.mcmc = mcmc[,which(substr(varnames(mcmc),1,nchar(var)+1)==paste(var,"[",sep=""))]
varnames(var.mcmc)
risk <- summary(var.mcmc)$statistics[,1]
windows()
plot(1:length(risk),risk, type = "l", col = "red", xlab = "Time", ylab = "Risk", main = "Risk measured as P(B<B_MSY)")



# ------------------------------------------------------------------------------------
# Customized nice graphs
# ------------------------------------------------------------------------------------

# ---------------------------------------------------
# Joint and marginal posterior distributions - 3 param
# ---------------------------------------------------

x1 <- mcmc.table[,'K']
x2 <- mcmc.table[,'r']
x3 <- mcmc.table[,'C_MSY']

def.par <- par(no.readonly = TRUE)
windows()

par(pch='.')

pairs( cbind(x1,x2,x3),labels=c("K","r",expression(C[MSY])), 
	 lower.panel=panel.smooth,
	 diag.panel=panel.dens, 
	 upper.panel=panel.cor,
	 cex.labels = 2.0, font.labels=1.2)


# ---------------------------------------------------
# Joint and marginal posterior distributions - 5 param
# ---------------------------------------------------

x1 <- mcmc.table[,'K']
x2 <- mcmc.table[,'r']
x3 <- mcmc.table[,'C_MSY']
x4 <- mcmc.table[,'q']
x5 <- mcmc.table[,'sigma2p']

def.par <- par(no.readonly = TRUE)
windows()

par(pch='.')

pairs( cbind(x1,x2,x3,x4,x5),labels=c("K","r",expression(C[MSY]),"q",expression(sigma[p]^2)), 
       lower.panel=panel.smooth,
       diag.panel=panel.dens, 
       upper.panel=panel.cor,
       cex.labels = 2.0, font.labels=1.2)


# ---------------------------------------------------
# 2D posterior distribution with contours
# ---------------------------------------------------

library(MASS) # necessaire pour la fonction kde2d qui permet de faire les contours

dim1 <- mcmc.table[,'K']
dim2 <- mcmc.table[,'r']

def.par <- par(no.readonly = TRUE)
windows()
f.density.bivar(dim1,dim2,nlevels=3,nb.points=2000)
par(def.par)




# --------------------------------------------------------------------------
# Nice graphs - Boxplot of time series
# --------------------------------------------------------------------------

n_obs <- data$n_obs
n_proj <- data$n_proj
n <- n_obs + n_proj
I_obs <- data$I_obs
Year = data$Year
Year = c(Year, seq(from = max(Year), to = max(Year)+n_proj-1, by=1))

size.text <- 1.3
size.labels <- 1
box.size <- 0.5
col.obs <- "white"
# col.proj <- "grey"
col.proj <- "grey75"
col.proj <- "pink"
# col2 <- "darkgrey"
col2 <- "grey30"
col <- c(rep(col.obs, times=n_obs),rep(col.proj, times=n_proj))
col2 <- rep(col2, times=n)


# ----------------------------------------------------------------------------
# Marginal posterior distribution - MSY
# ----------------------------------------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

# res <- 6
# name_figure <- "PostMSY.png"
# png(filename = name_figure, height = 500*res, width = 500*res, res=72*res)

par(mar=c(5,5,1,1), bty="n")

densMSY <- density(mcmc.table$C_MSY)
densMSYp <- density(mcmc.table$C_MSY_p)

plot(densMSY ,ylab = "", xlab = "", xaxt="n", yaxt="n", main="", xlim = c(0,1000), col = "red", type = "l", lty = 1, lwd = 2)
points(densMSYp ,ylab = "", xlab = "", xaxt="n", yaxt="n", main="", xlim = c(0,1000), col = "blue", type = "l", lty = 2, lwd = 2)

axis(side = 1, tick = T, cex.axis=1, las=1)
axis(side = 2, tick = T, lwd.ticks=0, labels=NA, cex.axis=1)
mtext(side=1, "MSY (x 1000 tons)", bty="n", line=3, cex = size.text)
mtext(side=2, "Density", bty="n", line=3, cex = size.text)

legend(	legend = c("Prior","Posterior"),
 		col = c("blue","red"), 	
 		lty = c(2,1), lwd = c(2,2), x = "topright", cex = size.text, bty ="n") 

par(def.par)

# dev.off()


# ----------------------------------------------
# Biomass
# ----------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

# res <- 6
# name_figure <- "Biomass.png"
# png(filename = name_figure, height = 500*res, width = 650*res, res=72*res)

par(mfrow = c(1,1), bty="n", mar=c(5,5,1,1))

x = "B"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]

x.min = 1 ; x.max = n+1
y.min = 0 ; y.max = 5000
line.y = y.max*2
y.label = "Biomass"
x.label = "Years"
title = "Biomass (X 1000 tons)"

boxplot(	X[,1], 
		xaxt = "n", yaxt = "n", 
		xlim = c(x.min,x.max), ylim = c(y.min,y.max), 
		at = 1, outpch = NA, boxwex = box.size, col = col[1]) 
for (i in 2:n) {
boxplot(X[,i], xaxt = "n", yaxt = "n", at = i, add = T, outpch = NA, boxwex = box.size, col = col[i]) }

points(x=x.min:(x.max-1),y=rep(line.y,n),col="red",type="l",lty=2,lwd="1")
axis(side =1, at=1:n, labels = Year, las=3, cex.axis=size.labels)
axis(side =2, cex.axis=size.labels)
mtext(x.label, line = 3, side = 1, las=1, cex = size.text)
mtext(y.label, line = 2, side = 2, cex = size.text)

q_B_MSY <- quantile(mcmc.table$B_MSY,probs=c(0.05,0.5,0.95))
abline(h=q_B_MSY[2], col = "green", lty=4, lwd = 3)
abline(h=q_B_MSY[1], col = "green", lty=4, lwd = 1)
abline(h=q_B_MSY[3], col = "green", lty=4, lwd = 1)

par(def.par)

# dev.off()



# ----------------------------------------------
# Harvest Rate
# ----------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

# res <- 6
# name_figure <- "Harvestrate.png"
# png(filename = name_figure, height = 500*res, width = 650*res, res=72*res)

par(mfrow = c(1,1), bty="n", mar=c(5,5,1,1))

x = "h"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]

x.min = 1 ; x.max = n+1
y.min = 0 ; y.max = 1.0
line.y = y.max*2
y.label = "Harvest Rate"
x.label = "Years"
title = "Harvest Rate (h=C/B)"

boxplot(	X[,1], 
		xaxt = "n", yaxt = "n", 
		xlim = c(x.min,x.max), ylim = c(y.min,y.max), 
		at = 1, outpch = NA, boxwex = box.size, col = col[1]) 
for (i in 2:n) {
boxplot(X[,i], xaxt = "n", yaxt = "n", at = i, add = T, outpch = NA, boxwex = box.size, col = col[i]) }

points(x=x.min:(x.max-1),y=rep(line.y,n),col="red",type="l",lty=2,lwd="1")
axis(side =1, at=1:n, labels = Year, las=3, cex.axis=size.labels)
axis(side =2, cex.axis=size.labels)
mtext(x.label, line = 3, side = 1, las=1, cex = size.text)
mtext(y.label, line = 2.5, side = 2, cex = size.text)

q_h_MSY <- quantile(mcmc.table$h_MSY,probs=c(0.05,0.5,0.95))
abline(h=q_h_MSY[2], col = "green", lty=4, lwd = 3)
abline(h=q_h_MSY[1], col = "green", lty=4, lwd = 1)
abline(h=q_h_MSY[3], col = "green", lty=4, lwd = 1)

par(def.par)

# dev.off()


# ----------------------------------------------------------------------------
# Equilibrium curve
# ----------------------------------------------------------------------------

windows()

x = "C_e"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]

col = c("lightgray","gray")

q.med <- apply(X,2,median)

q.inf <- NULL
q.sup <- NULL
for (k in 1:dim(X)[2])
{
q.inf[k] <- quantile(X[,k],probs=0.05)
q.sup[k] <- quantile(X[,k],probs=0.95)
}

q.inf1 <- NULL
q.sup1 <- NULL
for (k in 1:dim(X)[2])
{
q.inf1[k] <- quantile(X[,k],probs=0.25)
q.sup1[k] <- quantile(X[,k],probs=0.75)
}

plot(x=B_e,y=q.sup*2,ylab="Equilibrium catches",xlab="Biomass at equilibrium",type="n")
polygon(x=c(B_e,rev(B_e)),y=c(q.sup,rev(q.inf)),col=col[1],border=NA)
polygon(x=c(B_e,rev(B_e)),y=c(q.sup1,rev(q.inf1)),col = col[2],border=NA)
points(x=B_e,y=q.med,yaxt="n",col="red",lwd=2,type="l")

q_C_MSY <- quantile(mcmc.table$C_MSY,probs=c(0.05,0.5,0.95))
q_B_MSY <- quantile(mcmc.table$B_MSY,probs=c(0.05,0.5,0.95))

abline(v=q_B_MSY[2], col = "green", lty=4, lwd = 3)
abline(v=q_B_MSY[1], col = "green", lty=4, lwd = 1)
abline(v=q_B_MSY[3], col = "green", lty=4, lwd = 1)

abline(h=q_C_MSY[2], col = "green", lty=4, lwd = 3)
abline(h=q_C_MSY[1], col = "green", lty=4, lwd = 1)
abline(h=q_C_MSY[3], col = "green", lty=4, lwd = 1)

x = "B"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]
B.med <- apply(X,2,median)

points(B.med,data$C_obs[1:length(B.med)], type = "b", col="blue", pch=22)


# ----------------------------------------------
# Catches - predicted vs observed
# ----------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

par(mfrow = c(1,1), mar=c(5,5,1,1))

col = c("lightgray","gray")

# x = "C_pred"
x = "C"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]
n_plot <- dim(X)[2]

q.mean <- apply(X,2,mean)
q.inf <- NULL
q.sup <- NULL
for (k in 1:dim(X)[2])
{
  q.inf[k] <- quantile(X[,k],probs=0.05)
  q.sup[k] <- quantile(X[,k],probs=0.95)
}

q.inf1 <- NULL
q.sup1 <- NULL
for (k in 1:dim(X)[2])
{
  q.inf1[k] <- quantile(X[,k],probs=0.25)
  q.sup1[k] <- quantile(X[,k],probs=0.75)
}

plot(x=Year[1:n_plot],y=q.sup*1.2,ylim=c(0,max(q.sup*1.2)),ylab="Catches",type="n",las=1, main = "Fit Catches")
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup,rev(q.inf)),col = col[1],border=NA)
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup1,rev(q.inf1)),col = col[2],border=NA)
points(x=Year[1:n_plot], y=q.mean, col="red", lwd=2, type="l")
points(Year[1:n_plot],data$C_obs[1:n_plot],pch=20,cex=1,col="blue",type="p")

legend(	legend = c("Cred interv. (95%)","C pred.","C obs."),
        col = c("lightgray","red","blue"), 	
        lty = c(3,2,NA), pch=c(NA,NA,20), lwd = c(2,2,NA), x = "topright", cex = size.text, bty ="n") 

# dev.off()



# ----------------------------------------------
# Abundance indices - predicted vs observed
# ----------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

par(mfrow = c(1,1), mar=c(5,5,1,1))

col = c("lightgray","gray")

# x = "I_pred"
x = "I"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]
n_plot <- dim(X)[2]

q.mean <- apply(X,2,mean)
q.inf <- NULL
q.sup <- NULL
for (k in 1:dim(X)[2])
{
  q.inf[k] <- quantile(X[,k],probs=0.05)
  q.sup[k] <- quantile(X[,k],probs=0.95)
}

q.inf1 <- NULL
q.sup1 <- NULL
for (k in 1:dim(X)[2])
{
  q.inf1[k] <- quantile(X[,k],probs=0.25)
  q.sup1[k] <- quantile(X[,k],probs=0.75)
}

plot(x=Year[1:n_plot],y=q.sup*1.2,ylim=c(0,max(q.sup*1.2)),ylab="Abundance indices",type="n",las=1, main = "Fit Ab Indices")
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup,rev(q.inf)),col = col[1],border=NA)
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup1,rev(q.inf1)),col = col[2],border=NA)
points(x=Year[1:n_plot], y=q.mean, col="red", lwd=2, type="l")
points(Year[1:n_plot],data$I_obs[1:n_plot],pch=20,cex=1,col="blue",type="p")

legend(	legend = c("Cred interv. (95%)","IA pred.","IA obs."),
        col = c("lightgray","red","blue"), 	
        lty = c(3,2,NA), pch=c(NA,NA,20), lwd = c(2,2,NA), x = "topright", cex = size.text, bty ="n") 

# dev.off()


# ----------------------------------------------
# Biomass
# ----------------------------------------------

def.par <- par(no.readonly = TRUE)

windows()

par(mfrow = c(1,1), mar=c(5,5,1,1))

col = c("lightgray","gray")

# x = "C_pred"
x = "B"
mcmc <- window(mcmc)
mcmc.table <- as.data.frame(as.matrix(mcmc))
X = mcmc.table[,which(substr(colnames(mcmc.table),1,nchar(x)+1)==paste(x,"[",sep=""))]
n_plot <- dim(X)[2]

q.mean <- apply(X,2,mean)
q.inf <- NULL
q.sup <- NULL
for (k in 1:dim(X)[2])
{
  q.inf[k] <- quantile(X[,k],probs=0.05)
  q.sup[k] <- quantile(X[,k],probs=0.95)
}

q.inf1 <- NULL
q.sup1 <- NULL
for (k in 1:dim(X)[2])
{
  q.inf1[k] <- quantile(X[,k],probs=0.25)
  q.sup1[k] <- quantile(X[,k],probs=0.75)
}

plot(x=Year[1:n_plot],y=q.sup*1.2,ylim=c(0,max(q.sup*1.2)),ylab="Biomass",type="n",las=1, main = "Biomass")
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup,rev(q.inf)),col = col[1],border=NA)
polygon(x=c(Year[1:n_plot],rev(Year[1:n_plot])),y=c(q.sup1,rev(q.inf1)),col = col[2],border=NA)
points(x=Year[1:n_plot], y=q.mean, col="red", lwd=2, type="l", pch = "NA")

q_B_MSY <- quantile(mcmc.table$B_MSY,probs=c(0.05,0.5,0.95))
abline(h=q_B_MSY[2], col = "green", lty=4, lwd = 3)
abline(h=q_B_MSY[1], col = "green", lty=4, lwd = 1)
abline(h=q_B_MSY[3], col = "green", lty=4, lwd = 1)

legend(	legend = c("B-95%","B-50%", "B median"),
        col = c("lightgray","gray","red"), 	
        lty = c(3,2,1), pch=c(NA,NA,NA), lwd = c(2,2,2), x = "topright", cex = size.text, bty ="n") 

# dev.off()
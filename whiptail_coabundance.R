#load packages
library(R2jags)

#####load species data#####
Data <- read.csv("whiptail_data.csv", header=T)

#create matrices of lizard abundances by species
asma <- as.matrix(Data[,c("asma.01", "asma.02", "asma.03", "asma.04", "asma.05",
                          "asma.06", "asma.07", "asma.08", "asma.09", "asma.10")])

asun <- as.matrix(Data[,c("asun.01", "asun.02", "asun.03", "asun.04", "asun.05",
                          "asun.06", "asun.07", "asun.08", "asun.09", "asun.10")])

#####two-species n-mixture models:  model 1 = treatment effect#####

nsite <- 32
npair <- 16
nrep <- 10

#site covariates
pair <- Data$pair
trt <- Data$trt

#survey covariates
temp <- as.matrix(Data[,c("tmp.avg.c.01", "tmp.avg.c.02", "tmp.avg.c.03", "tmp.avg.c.04", "tmp.avg.c.05",
                          "tmp.avg.c.06", "tmp.avg.c.07", "tmp.avg.c.08", "tmp.avg.c.09", "tmp.avg.c.10")])
temp.s <- (temp-mean(temp))/sd(temp) #standardize covariates
temp.s2 <- temp.s^2 #quadratic term

humid <-as.matrix(Data[,c("rh.avg.p.01", "rh.avg.p.02", "rh.avg.p.03", "rh.avg.p.04", "rh.avg.p.05",
                          "rh.avg.p.06", "rh.avg.p.07", "rh.avg.p.08", "rh.avg.p.09", "rh.avg.p.10")])
humid.s <- (humid-mean(humid))/sd(humid)

wind <- as.matrix(Data[,c("win.avg.ms.01", "win.avg.ms.02", "win.avg.ms.03", "win.avg.ms.04", "win.avg.ms.05",
                          "win.avg.ms.06", "win.avg.ms.07", "win.avg.ms.08", "win.avg.ms.09", "win.avg.ms.10")])
wind.s <- (wind-mean(wind))/sd(wind)

time <- as.matrix(Data[,c("survtm.m.01", "survtm.m.02", "survtm.m.03", "survtm.m.04", "survtm.m.05",
                          "survtm.m.06", "survtm.m.07", "survtm.m.08", "survtm.m.09", "survtm.m.10")])
time.s <- (time-mean(time))/sd(time)

#bundle data
str(win.data <- list(y1 = asun, y2 = asma, nsite=nsite, nrep=dim(asun)[2], npair=npair,
                     pair=pair, trt=trt, 
                     temp=temp.s, temp2=temp.s2, humid=humid.s, wind=wind.s, time=time.s))

#specify model
sink("model_trt.txt")
cat("
    model{
    
    #priors
    for (i in 1:2) {
    a0[i] ~ dnorm(0,0.01)  #detection intercept
    a1[i] ~ dnorm(0,0.01)  #temp slope
    a2[i] ~ dnorm(0,0.01)  #temp^2 slope
    a3[i] ~ dnorm(0,0.01)  #humidity slope
    a4[i] ~ dnorm(0,0.01)  #wind slope
    a5[i] ~ dnorm(0,0.01)  #survey time slope
    b0[i] ~ dnorm(0,0.01)  #abundance intercept
    b1[i] ~ dnorm(0,0.01)  #treatment slope
    }
    
    b2 ~ dnorm(0,0.01)     #slope of y1 on y2
    b3 ~ dnorm(0,0.01)     #slope of y1 on y2*trt
    
    #random pair effect on abundance
    tau.b4 <- 1/(sd.b4*sd.b4)
    sd.b4 ~ dunif(0,10)
    tau.b5 <- 1/(sd.b5*sd.b5)
    sd.b5 ~ dunif(0,10)
    
    for (k in 1:npair) {
    b4[k] ~ dnorm(0,tau.b4)
    b5[k] ~ dnorm(0,tau.b5)
    }

    #site by survey heterogeneity in detection
    tau.a6 <- 1/(sd.a6*sd.a6)
    sd.a6 ~ dunif(0,2)
    tau.a7 <- 1/(sd.a7*sd.a7)
    sd.a7 ~ dunif(0,2)

    for(j in 1:nsite) {
      for(k in 1:nrep) {
        a6[j,k] ~ dnorm(0, tau.a6)
        a7[j,k] ~ dnorm(0, tau.a7)
      }
    }

    #ecological model
    for(j in 1:nsite){
    loglam1[j] <- b0[1] + b1[1]*trt[j] + b4[pair[j]]
    loglam1.lim[j] <- min(250, max(-250, loglam1[j])) #stabilize log
    lam1[j] <- exp(loglam1.lim[j])
    N1[j] ~ dpois(lam1[j])
    
    loglam2[j] <- b0[2] + b1[2]*trt[j] + b2*N1[j] + b3*N1[j]*trt[j] + b5[pair[j]] 
    lam2[j] <- exp(loglam2[j])
    N2[j] ~ dpois(lam2[j])
    
    #observation model for replicated counts
    for(k in 1:nrep){
    y1[j,k] ~ dbin(p1[j,k], N1[j])
    p1[j,k] <- exp(lp[j,k])/(1+exp(lp[j,k]))
    lp[j,k]<- a0[1] + a1[1]*temp[j,k] + a2[1]*temp2[j,k] + a3[1]*humid[j,k] + a4[1]*wind[j,k] + a5[1]*time[j,k] + a6[j,k]
    
    y2[j,k] ~ dbin(p2[j,k], N2[j])
    p2[j,k] <- exp(llp[j,k])/(1+exp(llp[j,k]))
    llp[j,k] <- a0[2] + a1[2]*temp[j,k] + a2[2]*temp2[j,k] + a3[2]*humid[j,k] + a4[2]*wind[j,k] + a5[2]*time[j,k] + a7[j,k]
    }
    }
    
    #derived quantities

    for(i in 1:2){
    mlambda[i] <- exp(b0[i]) #expected abundance on natural scale
    logit(mp[i]) <- a0[i]  #mean detection on natural scale
    }
    
    #mean abundance on treated and control sites
    MeanN1Trt <- (N1[2] + N1[4] + N1[6] + N1[8] + N1[10] + N1[12] + N1[14] + N1[16] + 
                  N1[18] + N1[20] + N1[22] + N1[24] + N1[26] + N1[28] + N1[30] + N1[32])/16
    
    MeanN1Con <- (N1[1] + N1[3] + N1[5] + N1[7] + N1[9] + N1[11] + N1[13] + N1[15] +
                  N1[17] + N1[19] + N1[21] + N1[23] + N1[25] + N1[27] + N1[29] + N1[31])/16

    MeanN2Trt <- (N2[2] + N2[4] + N2[6] + N2[8] + N2[10] + N2[12] + N2[14] + N2[16] + 
                  N2[18] + N2[20] + N2[22] + N2[24] + N2[26] + N2[28] + N2[30] + N2[32])/16
    
    MeanN2Con <- (N2[1] + N2[3] + N2[5] + N2[7] + N2[9] +  N2[11] + N2[13] + N2[15] +
                  N2[17] + N2[19] + N2[21] + N2[23] + N2[25] + N2[27] + N2[29] + N2[31])/16
    
    #difference in abundance between treated and control plots within pairs
    for(i in 1:npair) {
      N1dif[i] <- N1[i*2] - N1[(i*2)-1]
    }
  
    for(i in 1:npair) {
      N2dif[i] <- N2[i*2] - N2[(i*2)-1]
    }

    N1dif.mean <- mean(N1dif) 
    N2dif.mean <- mean(N2dif)
    
    }
    ",fill=TRUE)
sink()

#initial values
sp.inits = function() {
  list(a0=rnorm(2), a1=rnorm(2), a2=rnorm(2), a3=rnorm(2), a4=rnorm(2), a5=rnorm(2),
       b0=rnorm(2), b1=rnorm(2), b2=rnorm(1), b3=rnorm(1), b4=rnorm(16), b5=rnorm(16),
       N1 = as.vector(apply(asun,1,max, na.rm=T) + 5),   
       N2 = as.vector(apply(asma,1,max, na.rm=T) + 5)
  )     
}


#parameters monitored
params <- c("b0", "b1", "b2", "b3", "b4", "b5", "a0", "a1", "a2", "a3", "a4", "a5",
            "sd.b4", "sd.b5", "sd.a6", "sd.a7",
            "N1", "N2", "mlambda", "mp", 
            "MeanN1Trt", "MeanN1Con", "MeanN2Trt", "MeanN2Con", "N1dif", "N2dif", "N1dif.mean", "N2dif.mean")


#MCMC settings
ni <- 501000; nt <- 100; nb <- 1000; nc <- 3

m.trt <- jags(win.data, sp.inits, params, "model_trt.txt", n.chains = nc, n.thin = nt, 
              n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image("~/Dropbox/Documents/Research/Projects/RNM_AFRI/ExpPairs/sex_comp/ms/data_analysis/whiptail_coabundance.RData")



#####two-species n-mixture model; model2: vegetation PCs#####

nsite <- 32
npair <- 16
nrep <- 10

#site covariates
pair <- Data$pair
trt <- Data$trt
pc1.s <- (Data$PC1.six - mean(Data$PC1.six)) / sd(Data$PC1.six)
pc2.s <- (Data$PC2.six - mean(Data$PC2.six)) / sd(Data$PC2.six)

#survey covariates
temp <- as.matrix(Data[,c("tmp.avg.c.01", "tmp.avg.c.02", "tmp.avg.c.03", "tmp.avg.c.04", "tmp.avg.c.05",
                          "tmp.avg.c.06", "tmp.avg.c.07", "tmp.avg.c.08", "tmp.avg.c.09", "tmp.avg.c.10")])
temp.s <- (temp-mean(temp))/sd(temp) #standardize covariates
temp.s2 <- temp.s^2 #quadratic term

humid <-as.matrix(Data[,c("rh.avg.p.01", "rh.avg.p.02", "rh.avg.p.03", "rh.avg.p.04", "rh.avg.p.05",
                          "rh.avg.p.06", "rh.avg.p.07", "rh.avg.p.08", "rh.avg.p.09", "rh.avg.p.10")])
humid.s <- (humid-mean(humid))/sd(humid)

wind <- as.matrix(Data[,c("win.avg.ms.01", "win.avg.ms.02", "win.avg.ms.03", "win.avg.ms.04", "win.avg.ms.05",
                          "win.avg.ms.06", "win.avg.ms.07", "win.avg.ms.08", "win.avg.ms.09", "win.avg.ms.10")])
wind.s <- (wind-mean(wind))/sd(wind)

time <- as.matrix(Data[,c("survtm.m.01", "survtm.m.02", "survtm.m.03", "survtm.m.04", "survtm.m.05",
                          "survtm.m.06", "survtm.m.07", "survtm.m.08", "survtm.m.09", "survtm.m.10")])
time.s <- (time-mean(time))/sd(time)


#bundle data
str(win.data <- list(y1 = asun, y2 = asma, nsite=nsite, nrep=dim(asun)[2], npair=npair,
                     pair=pair, pc1=pc1.s, pc2=pc2.s,
                     temp=temp.s, temp2=temp.s2, humid=humid.s, wind=wind.s, time=time.s))

#specify model
sink("model_veg.txt")
cat("
    model{
    
    #priors

    for (i in 1:2) {
    a0[i] ~ dnorm(0,0.01)  #detection intercept
    a1[i] ~ dnorm(0,0.01) #temp slope
    a2[i] ~ dnorm(0,0.01)  #temp^2 slope
    a3[i] ~ dnorm(0,0.01)  #humidity slope
    a4[i] ~ dnorm(0,0.01)  #wind slope
    a5[i] ~ dnorm(0,0.01)  #survey time slope
    b0[i] ~ dnorm(0,0.01)  #abundance intercept
    b1[i] ~ dnorm(0,0.01)  #pc1
    b2[i] ~ dnorm(0,0.01)  #pc2
    }
    
    b3 ~ dnorm(0,0.01)  #slope of y1 on y2
    b4 ~ dnorm(0,0.01)  #pc2*N
    
    #random pair effect on abundance
    sd.b5 ~ dunif(0,10)
    tau.b5 <- 1/(sd.b5*sd.b5)
    
    sd.b6 ~ dunif(0,10)
    tau.b6 <- 1/(sd.b6*sd.b6)
    
    for (k in 1:npair) {
    b5[k] ~ dnorm(0,tau.b5)
    b6[k] ~ dnorm(0,tau.b6)
    }

    #site by survey heterogeneity in detection
    for(j in 1:nsite) {
    for(k in 1:nrep) {
    a6[j,k] ~ dnorm(0, tau.a6)
    a7[j,k] ~ dnorm(0, tau.a7)
    }
    }
    
    tau.a6 <- 1/(sd.a6*sd.a6)
    sd.a6 ~ dunif(0,2)
    
    tau.a7 <- 1/(sd.a7*sd.a7)
    sd.a7 ~ dunif(0,2)
    
    #ecological model
    for(j in 1:nsite){
    loglam1[j] <- b0[1] + b1[1]*pc1[j] + b2[1]*pc2[j] + b5[pair[j]]
    lam1[j] <- exp(loglam1[j])
    N1[j] ~ dpois(lam1[j])
    
    loglam2[j] <- b0[2] + b1[2]*pc1[j] + b2[2]*pc2[j] + b3*N1[j] + b4*N1[j]*pc2[j] + b6[pair[j]] 
    loglam2.lim[j] <- min(250, max(-250, loglam2[j])) #stabilize log
    lam2[j] <- exp(loglam2.lim[j])
    N2[j] ~ dpois(lam2[j])
    
    #observation model for replicated counts
    for(k in 1:nrep){
    y1[j,k] ~ dbin(p1[j,k], N1[j])
    p1[j,k] <- exp(lp[j,k])/(1+exp(lp[j,k]))
    lp[j,k]<- a0[1] + a1[1]*temp[j,k] + a2[1]*temp2[j,k] + a3[1]*humid[j,k] + a4[1]*wind[j,k] + a5[1]*time[j,k] + a6[j,k] 
    
    y2[j,k] ~ dbin(p2[j,k], N2[j])
    p2[j,k] <- exp(llp[j,k])/(1+exp(llp[j,k]))
    llp[j,k] <- a0[2] + a1[2]*temp[j,k] + a2[2]*temp2[j,k] + a3[2]*humid[j,k] + a4[2]*wind[j,k] + a5[2]*time[j,k] + a7[j,k]
    }
    }
    
    #derived quantities

    for(i in 1:2){
    mlambda[i] <- exp(b0[i]) #expected abundance on natural scale
    logit(mp[i]) <- a0[i]  #mean detection on natural scale
    }
    
    #mean abundance on treated and control sites
    MeanN1Trt <- (N1[2] + N1[4] + N1[6] + N1[8] + N1[10] + N1[12] + N1[14] + N1[16] + 
    N1[18] + N1[20] + N1[22] + N1[24] + N1[26] + N1[28] + N1[30] + N1[32])/16
    
    MeanN1Con <- (N1[1] + N1[3] + N1[5] + N1[7] + N1[9] + N1[11] + N1[13] + N1[15] +
    N1[17] + N1[19] + N1[21] + N1[23] + N1[25] + N1[27] + N1[29] + N1[31])/16
    
    MeanN2Trt <- (N2[2] + N2[4] + N2[6] + N2[8] + N2[10] + N2[12] + N2[14] + N2[16] + 
    N2[18] + N2[20] + N2[22] + N2[24] + N2[26] + N2[28] + N2[30] + N2[32])/16
    
    MeanN2Con <- (N2[1] + N2[3] + N2[5] + N2[7] + N2[9] +  N2[11] + N2[13] + N2[15] +
    N2[17] + N2[19] + N2[21] + N2[23] + N2[25] + N2[27] + N2[29] + N2[31])/16
    
    #difference in abundance between treated and control plots within pairs
    for(i in 1:npair) {
    N1dif[i] <- N1[i*2] - N1[(i*2)-1]
    }
    
    for(i in 1:npair) {
    N2dif[i] <- N2[i*2] - N2[(i*2)-1]
    }
    
    N1dif.mean <- mean(N1dif) 
    N2dif.mean <- mean(N2dif)
    
    
    }
    ",fill=TRUE)
sink()

#initial values
sp.inits = function() {
  list(a0=rnorm(2), a1=rnorm(2), a2=rnorm(2), a3=rnorm(2), a4=rnorm(2), a5=rnorm(2),
       b0=rnorm(2), b1=rnorm(2), b2=rnorm(2), b3=rnorm(1), b4=rnorm(1), b5=rnorm(16), b6=rnorm(16),
       N1 = as.vector(apply(asun,1,max, na.rm=T) + 5),   
       N2 = as.vector(apply(asma,1,max, na.rm=T) + 5)
  )     
}

#parameters monitored
params <- c("b0", "b1", "b2", "b3", "b4", "b5", "b6", "a0", "a1", "a2", "a3", "a4", "a5",
            "sd.b5", "sd.b6", "sd.a6", "sd.a7",
            "N1", "N2", "mlambda", "mp", 
            "MeanN1Trt", "MeanN1Con", "MeanN2Trt", "MeanN2Con", "N1dif", "N2dif", "N1dif.mean", "N2dif.mean")


#MCMC settings
ni <- 501000; nt <- 100; nb <- 1000; nc <- 3

m.veg <- jags(win.data, sp.inits, params, "model_veg.txt", n.chains = nc, n.thin = nt, 
              n.iter = ni, n.burnin = nb, working.directory = getwd())


save.image("~/Dropbox/Documents/Research/Projects/RNM_AFRI/ExpPairs/sex_comp/ms/data_analysis/whiptail_coabundance.RData")

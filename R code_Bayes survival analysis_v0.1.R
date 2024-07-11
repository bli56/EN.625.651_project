
library("BayesPPDSurv")
data("melanoma")


#########################
### Frequentist stats ###
#########################

hist <- melanoma[melanoma$study=="1684",]
current <- melanoma[melanoma$study=="1690",]
n.intervals <- c(4,3)
nMC <- 10000
nBI <- 200

historical <- list(list(time=hist$failtime,
                        event=hist$rfscens,
                        X=as.matrix(hist[, "trt"]),
                        S=hist$stratum))

set.seed(1)
result <- phm.fixed.a0(time=current$failtime,
                      event=current$rfscens,
                      X=as.matrix(current[,"trt"]),
                      S=current$stratum,
                      historical=historical,
                      a0=0.5,
                      n.intervals=n.intervals,
                      nMC=nMC,
                      nBI=nBI)

quantile(result$beta_samples)

colMeans(result$lambda_samples[[1]])
colMeans(result$lambda_samples[[2]])
colMeans(result$lambda0_samples[[1]])
colMeans(result$lambda0_samples[[2]])


######################
### Bayesian Power ###
######################

nMC <- 10000
nBI <- 200

set.seed(1)
samples <- phm.fixed.a0(historical=historical,
                        a0=1,
                        n.intervals=n.intervals,
                        current.data=FALSE,
                        nMC=nMC,
                        nBI=nBI)

beta_priors <- samples$beta_samples
DN_beta_samp_prior <- as.matrix(beta_priors[beta_priors[,1] > 0, ])
DA_beta_samp_prior <- as.matrix(beta_priors[beta_priors[,1] < 0, ])
lambda_samp_prior <- samples$lambda_samples



set.seed(1)

n.events <- 350
n.subjects <- n.events * 3

result <- power.phm.fixed.a0(historical = historical,
                              a0 = 0.6,
                              n.subjects = n.subjects,
                              n.events = n.events,
                              n.intervals = n.intervals,
                              samp.prior.beta = DA_beta_samp_prior,
                              samp.prior.lambda = lambda_samp_prior,
                              dist.enroll = "Uniform",
                              param.enroll = 4,
                              nMC = 10000,
                              nBI = 200,
                              delta = 0,
                              nullspace.ineq = ">",
                              N = 10)


result$'power/type I error'


############################################################################

rm(list=ls(all=TRUE))
cat("\014")

### Import & extract data

data_simu <- read.csv("D:\\Baini_files\\ZM\\Studies\\CB03-154-201\\Sample size\\Data_simulated\\data_simu.csv")
data_simu <- subset(data_simu, select = c(data_id, subject_id, trt, event_status, stop))
data_simu <- unique(data_simu)


data_simu_0 <- data_simu[data_simu$trt == 0, ]
data_simu_0 <- data_simu_0[order(data_simu_0$trt,
                                 data_simu_0$data_id,
                                 data_simu_0$subject_id),]
n_trt0 <- 49
data_simu_0 <- by(data_simu_0, data_simu_0["data_id"], head, n=n_trt0)
data_simu_0 <- Reduce(rbind, data_simu_0)



data_simu_1 <- data_simu[data_simu$trt == 1, ]
data_simu_1 <- data_simu_1[order(data_simu_1$trt,
                                 data_simu_1$data_id,
                                 data_simu_1$subject_id),]
n_trt1 <- 98
data_simu_1 <- by(data_simu_1, data_simu_1["data_id"], head, n=n_trt1)
data_simu_1 <- Reduce(rbind, data_simu_1)

hist <- rbind(data_simu_0, data_simu_1)
hist <- hist[order(hist$data_id,
                   hist$trt,
                   hist$subject_id),]
lngth <- length(hist$trt)

### Standardize the format
stratum <- replicate(lngth, 1)
historical <- list(list(X = as.matrix(hist[, "trt"]),
                        time = hist$stop,
                        event = hist$event_status,
                        S = stratum))
                   
### Bayesian power analysis

set.seed(1)
a0 <- replicate(100, 0.01)
nMC = 10000
nBI =  2000

samples <- phm.fixed.a0(historical = historical,
                        a0 = a0,
                        n.intervals = 1,
                        current.data = FALSE,
                        nMC = nMC,
                        nBI = nBI)

beta_priors <- samples$beta_samples
DN_beta_samp_prior <- as.matrix(beta_priors[beta_priors[,1] < 0, ])
DA_beta_samp_prior <- as.matrix(beta_priors[beta_priors[,1] > 0, ])
lambda_samp_prior <- samples$lambda_samples


n.subjects <- 147
n.events <- floor(n.subjects*1)
param.enroll <- 365*2

power <- power.phm.fixed.a0(historical = historical,
                            a0 = a0,
                            n.subjects = n.subjects,
                            n.events = n.events,
                            n.intervals = 1,
                            
                            samp.prior.beta = DA_beta_samp_prior,
                            samp.prior.lambda = lambda_samp_prior,
                            dist.enroll = "Uniform",
                            param.enroll = param.enroll,
                            nMC = nMC,
                            nBI = nBI,
                            delta = 0,
                            nullspace.ineq = "<",
                            gamma = 0.90,
                            dist.csr = "Constant",
                            param.csr = 334,
                            rand.prob = 0.6666667, 
                            N = 100)

power$'power/type I error'


















































































































































































































































































































































































































































































































































































































































































































































































































































































































































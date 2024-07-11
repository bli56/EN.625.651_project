
rm(list=ls(all=TRUE))
cat("\014")


set.seed(777)
library(survsim)


# Specify iterations
n_iter <- 100

# Specify the sample size per iteration
n_ttl <- 600

for(i_iter in 1: n_iter){

      ######################################################
      ### Generate survival time and indicator for death ###
      ######################################################
      
      # In AMX0035 study, death occurred in 5 subjects (6%) who received active drug, and in 2 subjects (4%) who received placebo. So, the overall death rate is (5+2)/(89+48) = 0.05
      # In filtered PRO-ACT data, considering the censoring time window [0, 330+5-1], death occurred in 1 subject (1/22=4.545455%) who received active drug, and in 15 subjects (15/149=10.06711%) who received placebo.
      # So, at here we set the pooled death rate around 11% (i.e., both median and mean >= 110)
  
      event_data <- simple.surv.sim(n=n_ttl,
                                    foltime=334,
                                    
                                    dist.ev=c('weibull'),                # time-to-EVENT dist (by day)
                                    anc.ev=c(20.3),                      # ancillary parameter for time-to-EVENT dist
                                    beta0.ev=c(5.9),                     # beta0 param for time-to-EVENT dist --> Mean(deathD) = 260.4375, so the initial setting of 'beta0.ev' of simple.surv.sim() = log(260.4375) = 5.562363. Start from here we will tune its value as well as 'anc.ev'
                                    
                                    dist.cens = 'unif',                  # time-to-CENSORING dist (by day)
                                    anc.cens=334,                        # ancillary parameter for time-to-CENSORING dist, or MAX in case of uniform dist 
                                    beta0.cens=324,                      # beta0 for time-to-CENSORING dist, or MIN in case of uniform dist
                                    
                                    z=NULL,                              # info on random effect(s)
                                    beta=list(c(0)),                     # covariate effect(s) --> treatment group ??no effect considered??
                                    x=list(c("bern", 0.6666667))         # covariate dist --> allocation ratio = 2:1
      )
      
      dataid <- rep(i_iter, n_ttl)
      event_data2 <- cbind(dataid, event_data)
      colnames(event_data2)[colnames(event_data2) == "nid"] ="subject_id"
      colnames(event_data2)[colnames(event_data2) == "status"] ="event_status"
      assign(paste("status_stat_", i_iter, sep=""), data.frame(i_iter, table(event_data2$event_status)))

      # summary(event_data2$stop)
      # sd(event_data2$stop)

      # data_names <- mget(ls(pattern = "^status_stat_.*"))
      # library(plyr)
      # combined <- plyr::rbind.fill(data_names)
      # combined2<-combined[!(combined$Var1==0),]
      # summary(combined2$Freq)
      
      
      #########################
      ### Generate ALSFRS_R ###
      #########################
      
      # Generate subject ID
      subject_id <- seq(from = 1, to = n_ttl)
      
      # Visit (by day)
      visit_day <- c(1, 8, 22, 50, 78, 106, 134, 162,   190, 218, 246, 274, 302, 330)
      
      # Visit (by week)
      visit_week <- c(0, 2, 4, 8, 12, 16, 20, 24,  28, 32, 36, 40, 44, 48)
      
      # Visit (by month)
      visit_month <- c(0, 0.5, 1, 2, 3, 4, 5, 6,  7, 8, 9, 10, 11, 12)
      
      # Assign visit
      visit <- visit_month
      visit_info <- as.data.frame(cbind(visit_day, visit_week, visit_month))
      
      # Total number of visits
      n_visit <- as.integer(length(visit))
      
      trt <- event_data2$x
      dispoistion_info <-as.data.frame(cbind(subject_id, trt))
      table(dispoistion_info$trt)
      
      # Specify fixed effect coefficients
      b_intercept <- 37.5393
      b_visit <- -1.2342 
      b_visit_trt <- 0.1992
      
      # Generate random effect coefficients
      b_intercept_rand <- rnorm(n_ttl, mean= 0, sd= 5.2072)
      b_visit_rand <- rnorm(n_ttl, mean= 0, sd= 0.6952)
      
      
      ##############
      ### Case b ###
      ##############
      
      ### Generate the post-baseline ALSFRS_R for each subject, as below:
      # ALSFRS_R <- b_intercept + b_visit*visit + b_visit_trt*trt*visit + e
      
      n_ALSFRS_R <- n_ttl*n_visit
      e <- array(-999, dim=c(n_ttl, n_visit))
      ALSFRS_R <- array(-999, dim=c(n_ttl, n_visit))
      
      for(i in 1: n_ttl){
        for(j in 1: n_visit){
          
          
          #??? Generate individual-level random error per each visit
          e[i,j] <- rnorm(1, mean=0, sd=1.6801010817)
          
          ALSFRS_R[i, j] <- (b_intercept+ b_intercept_rand[i]) + (b_visit+b_visit_rand[i])*visit[j] + b_visit_trt*trt[i]*visit[j] + e[i,j]
          
          # Re-set ALSFRS_R values
          if (ALSFRS_R[i, j] < 0) {
            ALSFRS_R[i, j] = 0
          } else if (ALSFRS_R[i, j] > 40) {
            ALSFRS_R[i, j] = 40
          }
        }
      }
      
      # Stacking
      ALSFRS_R2 <- as.data.frame(cbind(subject_id, ALSFRS_R))
      ALSFRS_R3 <- cbind(ALSFRS_R2[1], stack(ALSFRS_R2[2:ncol(ALSFRS_R2)]))
      
      ALSFRS_R3$visit_week <- 999
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V2', 0,   ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V3', 2,   ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V4', 4,   ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V5', 8,   ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V6', 12,  ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V7', 16,  ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V8', 20,  ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V9', 24,  ALSFRS_R3$visit_week)
      
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V10', 28,  ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V11', 32,  ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V12', 36,  ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V13', 40,  ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V14', 44,  ALSFRS_R3$visit_week)
      ALSFRS_R3$visit_week <- ifelse(ALSFRS_R3$ind == 'V15', 48,  ALSFRS_R3$visit_week)
     
      colnames(ALSFRS_R3)[colnames(ALSFRS_R3) == "values"] ="ALSFRS_R_Total"
      ALSFRS_R4 <- merge(ALSFRS_R3, visit_info, by=c("visit_week"))
      ALSFRS_R5 <- merge(ALSFRS_R4, dispoistion_info, by=c("subject_id"))
      ALSFRS_R6 <- ALSFRS_R5[order(ALSFRS_R5$subject_id, ALSFRS_R4$visit_week),]
      
      # Merge
      sim_data <- merge(ALSFRS_R6, event_data2, by=c("subject_id"))
      sim_data <- sim_data[,!names(sim_data) %in% c("ind", "start", "z", "x")]
      
      # Censor the record if it is later than day of death
      sim_data2 <- sim_data[!sim_data$visit_day > sim_data$stop,]
      sim_data3 <- sim_data2[order(sim_data2$subject_id, sim_data2$visit_week),]
      
      new_order <- c('dataid', 'subject_id', 'trt', 'visit_month', 'visit_week', 'visit_day', 'ALSFRS_R_Total', 'event_status', 'stop')
      sim_data4 <- sim_data3[, new_order]
      
      # Create a dataset per iteration
      assign(paste("data_", i_iter, sep=""), data.frame(sim_data4))
      
      rm(sim_data, sim_data2, sim_data3, sim_data4,
         ALSFRS_R, ALSFRS_R2, ALSFRS_R3, ALSFRS_R4, ALSFRS_R5, ALSFRS_R6,
         event_data, event_data2)
}

# Combine simu
data_names <- mget(ls(pattern = "^data_.*"))
library(plyr)
data_combined <- plyr::rbind.fill(data_names)
data_simu <- data_combined[order(data_combined$dataid, data_combined$subject_id, data_combined$visit_week),]
colnames(data_simu)[colnames(data_simu) == "dataid"] ="data_id"
rm(data_combined)

# Export data
library("writexl")
write.csv(data_simu, "D:\\Baini_files\\ZM\\Studies\\CB03-154-201\\Sample size\\Data_simulated\\data_simu.csv")














































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































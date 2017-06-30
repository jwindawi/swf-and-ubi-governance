library(survey)
library(dplyr)
library(convey)
library(parallel)
library(foreach)
library(doMC)

years <- cbind(2005:2015)


##   1. Generating initial files
totak_use_person <- read.csv("/Users/jasonwindawi/Dropbox/AK Transfer/IPUMS files/totak_use_person.csv") #need to adjust this for Nobel

totak_person_sim <- totak_use_person %>% select(ind, YEAR, SERIAL, NUMPREC, HHINCOME, HHWT, PERWT, DIVSIM, Poverty,  
                                                INCSIM, APFDperson, PUMA, eq_spm, eq_sq, PROPNSIM, starts_with("REPWT"))
totak_house_sim <- totak_person_sim %>% distinct(YEAR, SERIAL, .keep_all = TRUE)  

store <- list()
              
##   2. Generating random income data - once done, can save this object and call it later
              
set.seed = 123456
reps = 1000
incadj <- matrix(0, nrow(totak_house_sim), reps) # matrix for storing binomials

for (k in 1:reps) {
  for (i in 1:nrow(totak_house_sim)){
    incadj[i,k] <- (rbinom(1, totak_house_sim$DIVSIM[i], totak_house_sim$PROPNSIM[i])) * totak_house_sim$APFDperson[i]#Inner loop to assign Divs to remove
    }
} 
write.csv(incadj, file = "/Users/jasonwindawi/Dropbox/AK Transfer/IPUMS files/incadj.csv")

##   3. Simulating Dividend removal

incadj <- read.csv("/Users/jasonwindawi/Dropbox/AK Transfer/IPUMS files/incadj.csv")
incadj <- incadj[,-1] #to account for the added column created by making a csv

set.seed = 123456
registerDoMC(cores = 4)
store <- foreach (k = 1:1000, .packages = c('survey', 'convey', 'dplyr')) %dopar% {
  totak_house_sim <- mutate(totak_house_sim, 
                          Divadj = incadj[,k],
                          INCSIM1 = INCSIM - Divadj, # At Household level, deducting simulated Divs to create counterfactual 
                          inc_spm = INCSIM1/eq_spm,  # Equivalized income at SPM scale
                          inc_sq = INCSIM1/eq_sq,    # Equivalized income at sq rt scale, becomes "_sq" below
                          t = 1)                     # Vector of 1s to use in svyratio denominator
  
            totak_house_frag <- select(totak_house_sim, YEAR, SERIAL, INCSIM1, inc_spm, inc_sq, t)
            totak_person_temp <- merge(totak_person_sim, totak_house_frag, by = intersect(names(totak_person_sim), names(totak_house_frag)))
            totak_person_temp <- arrange(totak_person_temp, ind)
                            
                          totak_house_temp <- select(totak_house_sim, ind, YEAR, SERIAL, HHWT, PUMA, Poverty, Divadj,
                                                     INCSIM1, inc_spm, inc_sq, t, starts_with("REPWT"), -starts_with("REPWTP"))     #To use w/survey package
                                howts <- select(totak_house_temp, starts_with("REPWT"))
                          
                          totak_person_temp <- select(totak_person_temp, ind, YEAR, SERIAL, PERWT, PUMA, Poverty, 
                                                      INCSIM1, inc_spm, inc_sq, t, starts_with("REPWTP"))#To use w/survey package
                                pewts <- select(totak_person_temp, starts_with("REPWTP"))
                          
                          # Creating temporary survey objects to get medians and attaching to data frames
                                
                          hsurvtemp <- svrepdesign(data = totak_house_temp, repweights = howts, weights = totak_house_temp$HHWT,
                                   combined.weights=TRUE, type = "Fay", rho = 0.5, scale= 4/80, rscales = rep(1, 80), mse=TRUE) #Alaskan households
                          psurvtemp <- svrepdesign(data = totak_person_temp, repweights = pewts, weights = totak_person_temp$PERWT,
                                               combined.weights=TRUE, type = "Fay", rho = 0.5, scale= 4/80, rscales = rep(1, 80), mse=TRUE) #Persons within households

                          meds_h <- svyby(~inc_spm+~inc_sq, ~YEAR, hsurvtemp, svyquantile, interval.type = "quantile", na.rm=TRUE,
                                          quantile=0.5, keep.var=FALSE, multicore = TRUE)
                          meds_p <- svyby(~inc_spm + ~inc_sq, ~YEAR, psurvtemp, svyquantile, interval.type = "quantile", na.rm=TRUE,
                                          quantile=0.5, keep.var=FALSE, multicore = TRUE)
                          
                                    colnames(meds_h) <- c("YEAR", "hmed_spm", "hmed_sq") 
                                    colnames(meds_p) <- c("YEAR", "pmed_spm", "pmed_sq") 

                                    totak_person_temp <- merge(totak_person_temp, meds_p, by = "YEAR")
                                    totak_person_temp <- arrange(totak_person_temp, ind)
                          
                                    totak_house_temp <- merge(totak_house_temp, meds_h, by = "YEAR")
                                    totak_house_temp <- arrange(totak_house_temp, ind)
                          
                           # Generating survey objects to use in the remainder of processes
                                    
                           hsurv <- update(hsurvtemp, hmed_spm = totak_house_temp$hmed_spm, hmed_sq = totak_house_temp$hmed_sq)
                                    hsurv <- convey_prep(hsurv)
                                    hsurv4 <- subset(hsurv, YEAR > 2011 & PUMA == 400) #Subsistence Alaska, subsetting AFTER convey_prep only
                           
                           psurv <- update(psurvtemp, pmed_spm = totak_person_temp$pmed_spm, pmed_sq = totak_person_temp$pmed_sq)
                                    psurv <- convey_prep(psurv)
                                    psurv4 <- subset(psurv, YEAR > 2011 & PUMA == 400) #Subsistence Alaska, subsetting AFTER convey_prep only
                         
                           divs <- svyby(~Divadj, ~YEAR, design = hsurv, svytotal) #Measuring coverage of Dividends removed in sim
                          
                                    # Generating ginis using "convey" package
                                    # As with all subsequent stats, doing this for 8 sets:
                                    # Households vs. persons; Two equiv scales for each; All AK and Subsistence AK
                                    gini_spm <- svyby(~inc_spm, ~YEAR, hsurv, svygini, na.rm=TRUE, multicore = TRUE)
                                    gini_sq <- svyby(~inc_sq, ~YEAR, hsurv, svygini, na.rm=TRUE, multicore = TRUE)
                                    gini_spmp <- svyby(~inc_spm, ~YEAR, psurv, svygini, na.rm=TRUE, multicore = TRUE)
                                    gini_sqp <- svyby(~inc_sq, ~YEAR, psurv, svygini, na.rm=TRUE, multicore = TRUE)
                                    gini4_spm <- svyby(~inc_spm, ~YEAR, hsurv4, svygini, na.rm=TRUE, multicore = TRUE)
                                    gini4_sq <- svyby(~inc_sq, ~YEAR, hsurv4, svygini, na.rm=TRUE, multicore = TRUE)
                                    gini4_spmp <- svyby(~inc_spm, ~YEAR, psurv4, svygini, na.rm=TRUE, multicore = TRUE)
                                    gini4_sqp <- svyby(~inc_sq, ~YEAR, psurv4, svygini, na.rm=TRUE, multicore = TRUE)
          
                                    # Generating income quantiles
                                    survq_h <- svyby(~inc_spm + ~inc_sq, ~YEAR, hsurv, svyquantile, interval.type = "quantile", na.rm=TRUE,
                                                        quantiles=c(0.9, 0.75, 0.5, 0.25, 0.1), keep.var=TRUE, multicore = TRUE)
                                    survq_p <- svyby(~inc_spm + ~inc_sq, ~YEAR, psurv, svyquantile, interval.type = "quantile", na.rm=TRUE,
                                                        quantiles=c(0.9, 0.75, 0.5, 0.25, 0.1), keep.var=TRUE, multicore = TRUE)
                                    survq4_h <- svyby(~inc_spm + ~inc_sq, ~YEAR, hsurv4, svyquantile, interval.type = "quantile", na.rm=TRUE,
                                                         quantiles=c(0.9, 0.75, 0.5, 0.25, 0.1), keep.var=TRUE, multicore = TRUE)
                                    survq4_p <- svyby(~inc_spm + ~inc_sq, ~YEAR, psurv4, svyquantile, interval.type = "quantile", na.rm=TRUE,
                                                        quantiles=c(0.9, 0.75, 0.5, 0.25, 0.1), keep.var=TRUE, multicore = TRUE)

                          # Poverty
                          povH <- svyby(~(INCSIM1 < Poverty) + ~(INCSIM1 < 1.25*Poverty), ~YEAR, denominator=~t, svyratio, design = hsurv, multicore = TRUE)
                          povHtot <- svyby(~(INCSIM1 < Poverty) + ~(INCSIM1 < 1.25*Poverty), ~YEAR, denominator=~t, svytotal, design = hsurv, multicore = TRUE)
                          povH2 <- svyby(~(INCSIM1 < 1.5*Poverty) + ~(INCSIM1 < 2*Poverty), ~YEAR, denominator=~t, svyratio, design = hsurv, multicore = TRUE)
                          povH2tot <- svyby(~(INCSIM1 < 1.5*Poverty) + ~(INCSIM1 < 2*Poverty), ~YEAR, denominator=~t, svytotal, design = hsurv, multicore = TRUE)
                          povH_spm <- svyby(~(inc_spm < 0.25*hmed_spm) + (inc_spm < 0.5*hmed_spm), ~YEAR, denominator=~t, svyratio, design = hsurv, multicore = TRUE)
                          povH_sq <- svyby(~(inc_sq < 0.25*hmed_sq) + (inc_spm < 0.5*hmed_sq), ~YEAR, denominator=~t, svyratio, design = hsurv, multicore = TRUE)
                                    
                          povP <- svyby(~(INCSIM1 < Poverty) + ~(INCSIM1 < 1.25*Poverty), ~YEAR, denominator=~t, svyratio, design = psurv, multicore = TRUE)
                          povPtot <- svyby(~(INCSIM1 < Poverty) + ~(INCSIM1 < 1.25*Poverty), ~YEAR, denominator=~t, svytotal, design = psurv, multicore = TRUE)
                          povP2 <- svyby(~(INCSIM1 < 1.5*Poverty) + ~(INCSIM1 < 2*Poverty), ~YEAR, denominator=~t, svyratio, design = psurv, multicore = TRUE)
                          povP2tot <- svyby(~(INCSIM1 < 1.5*Poverty) + ~(INCSIM1 < 2*Poverty), ~YEAR, denominator=~t, svytotal, design = psurv, multicore = TRUE)
                          povP_spm <- svyby(~(inc_spm < 0.25*pmed_spm) + (inc_spm < 0.5*pmed_spm), ~YEAR, denominator=~t, svyratio, design = psurv, multicore = TRUE)
                          povP_sq <- svyby(~(inc_sq < 0.25*pmed_sq) + (inc_spm < 0.5*pmed_sq), ~YEAR, denominator=~t, svyratio, design = psurv, multicore = TRUE)
                                    
                          povH4 <- svyby(~(INCSIM1 < Poverty) + ~(INCSIM1 < 1.25*Poverty), ~YEAR, denominator=~t, svyratio, design = hsurv4, multicore = TRUE)
                          povH4tot <- svyby(~(INCSIM1 < Poverty) + ~(INCSIM1 < 1.25*Poverty), ~YEAR, denominator=~t, svytotal, design = hsurv4, multicore = TRUE)
                          povH42 <- svyby(~(INCSIM1 < 1.5*Poverty) + ~(INCSIM1 < 2*Poverty), ~YEAR, denominator=~t, svyratio, design = hsurv4, multicore = TRUE)
                          povH42tot <- svyby(~(INCSIM1 < 1.5*Poverty) + ~(INCSIM1 < 2*Poverty), ~YEAR, denominator=~t, svytotal, design = hsurv4, multicore = TRUE)
                          povH4_spm <- svyby(~(inc_spm < 0.25*hmed_spm) + (inc_spm < 0.5*hmed_spm), ~YEAR, denominator=~t, svyratio, design = hsurv4, multicore = TRUE)
                          povH4_sq <- svyby(~(inc_sq < 0.25*hmed_sq) + (inc_sq < 0.5*hmed_sq), ~YEAR, denominator=~t, svyratio, design = hsurv4, multicore = TRUE)
                                    
                          povP4 <- svyby(~(INCSIM1 < Poverty) + ~(INCSIM1 < 1.25*Poverty), ~YEAR, denominator=~t, svyratio, design = psurv4, multicore = TRUE)
                          povP4tot <- svyby(~(INCSIM1 < Poverty) + ~(INCSIM1 < 1.25*Poverty), ~YEAR, denominator=~t, svytotal, design = psurv4, multicore = TRUE)
                          povP42 <- svyby(~(INCSIM1 < 1.5*Poverty) + ~(INCSIM1 < 2*Poverty), ~YEAR, denominator=~t, svyratio, design = psurv4, multicore = TRUE)
                          povP42tot <- svyby(~(INCSIM1 < 1.5*Poverty) + ~(INCSIM1 < 2*Poverty), ~YEAR, denominator=~t, svytotal, design = psurv4, multicore = TRUE)
                          povP4_spm <- svyby(~(inc_spm < 0.25*pmed_spm) + (inc_spm < 0.5*pmed_spm), ~YEAR, denominator=~t, svyratio, design = psurv4, multicore = TRUE)
                          povP4_sq <- svyby(~(inc_sq < 0.25*pmed_sq) + (inc_sq < 0.5*pmed_sq), ~YEAR, denominator=~t, svyratio, design = psurv4, multicore = TRUE)
                  
                store[[k]] <- c(divs, gini_spm, gini_sq, gini_spmp, gini_sqp, gini4_spm, gini4_sq, gini4_spmp, gini4_sqp, survq_h, survq_p, survq4_h, survq4_p, povH, povHtot,
                                povH2, povH2tot, povH_spm, povH_sq, povP, povPtot, povP2, povP2tot, povP_spm, povP_sq, povH4, povH4tot, povH42, povH42tot, povH4_spm, povH4_sq, 
                                povP4,povP4tot, povP42, povP42tot, povP4_spm, povP4_sq)                            
}

save.image("/Users/jasonwindawi/Dropbox/AK Transfer/newsim.RData")








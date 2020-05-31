##############################################################################
#### EUGCOMM: A DATABASE OF EUGLOSSINE BEE ASSEMBLAGES ON FRAGRANCE BAITS ####
##############################################################################

# The following code reads the raw data and builds the R list object EUGCOMM

library(plyr)
library(reshape2)
library(stringr)

rm(list=ls())

#write.csv(eugcomm, "Database/R files/eugcomm.csv")

indat = read.csv("Database/R files/eugcomm.csv")

indat$Species = paste(indat$genus, indat$species)
studies = unique(indat$reference)

EUGCOMM = list()

for(i in 1:length(studies)){
  sub = indat[indat$reference==studies[i],]
  bSU = data.frame(baits=tapply(sub$baits_clean, factor(sub$SU_ID), function(x) as.character(unique(x))))
  df = data.frame(SU_ID=rep(rownames(bSU), times=str_count(bSU$baits, "_")+1), bait=unlist(strsplit(bSU$baits, split="_")))
  baitMat = dcast(df, SU_ID~bait, value.var="bait", fun.aggregate = length)
  
  EUGCOMM[[i]] = list(Reference = as.character(paste(unique(sub$reference), unique(sub$journal))),
                      Study_areas = as.character(unique(sub$study_area)),
                      Metadata = ddply(sub, .(SU_ID, study_area), summarize,
                                              sampling.unit = unique(sampling.unit),
                                              lat = mean(lat),
                                              long = mean(long),
                                              startyear = mean(startyear),
                                              startmonth = mean(startmonth),
                                              endyear = mean(endyear),
                                              endmonth = mean(endmonth),
                                              duration = mean(study_duration_m),
                                              sampling_times = mean(sampling_times),
                                              bait_number = mean(bait_number),
                                              baiting_method = unique(baiting_method),
                                              sample_duration = mean(sample_duration),
                                              station_number = mean(number_stations)
                                      ),
                      BaitData = baitMat,
                      SpeciesData = dcast(sub, SU_ID ~ Species, value.var = "number", fun.aggregate=sum)
                    )
}

names(EUGCOMM) = studies

EUGCOMM[[12]]

save(EUGCOMM, file="Database/R files/EUGCOMM.RData")



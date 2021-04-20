library(Rmisc)
library(labelled)
library(psych)
library(dplyr)

# NOTE: To load data, you must download both the extract's data and the DDI
# and also set the working directory to the folder with these files (or change the path below).

if (!require("ipumsr")) stop("Reading IPUMS data into R requires the ipumsr package. It can be installed using the following command: install.packages('ipumsr')")

ddi <- read_ipums_ddi("c:/yost2/pums/usa_00003.xml")
data <- read_ipums_micro(ddi)

data2 <- select(data,SERIAL,STATEFIP,PUMA,US2019A_GRNTP,US2019A_VALP,US2019A_HINCP,HHWT,
                PERNUM, PERWT, US2019A_SCHL)

data2 <- dplyr::rename(data2,rent=US2019A_GRNTP,inc=US2019A_HINCP,val=US2019A_VALP,
                       edu=US2019A_SCHL)

#data2 <- data2[!duplicated(data2$SERIAL), ]

data2$edu <- zap_labels(data2$edu)
data2$STATEFIP <- zap_labels(data2$STATEFIP)
data2$PUMA <- zap_labels(data2$PUMA)

data2$rent <- as.numeric(data2$rent)
data2$inc <- as.numeric(data2$inc)
data2$val <- as.numeric(data2$val)
data2$edu <- as.numeric(data2$edu)

#within each household, I need the maximum edu value

data2a <- group_by(as.data.frame(data2),SERIAL)
data2a <- dplyr::summarise(data2a, SERIAL=SERIAL, STATEFIP=STATEFIP, PUMA=PUMA, rent=rent,
  val=val, inc=inc, edu=max(edu), HHWT=HHWT, .groups="keep")

rm(ddi)
rm(data)
rm(data2)

data2a <- data2a[!duplicated(data2a$SERIAL), ]

memory.limit(size = 31000) 

sum(data2a$HHWT)

data3 <- data2a[rep(row.names(data2a), data2a$HHWT), 1:7]

rm(data2a)

set.seed(6)

###############
# State level #
###############

results <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(results) <- c('STATEFIP','mgr','mhi','mhv','edu','mgrrank','mhirank','mhvrank',
                       'edurank','ML1','row') 

data3 <- dplyr::select(data3,-PUMA)

#sample each column separately in order to exclude NAs

data3a <- as.data.frame(filter(data3,!is.na(rent)))
data3b <- as.data.frame(filter(data3,!is.na(inc)))
data3c <- as.data.frame(filter(data3,!is.na(val)))
data3d <- as.data.frame(filter(data3,!is.na(edu)))

#take a sample from the non-missing values in each variable
#equal to the number of values, rounded to the nearest million

data4a <- data3a %>% sample_n(size=trunc(nrow(data3a)/1000000)*1000000,replace=T)
rm(data3a)
data4b <- data3b %>% sample_n(size=trunc(nrow(data3b)/1000000)*1000000,replace=T)
rm(data3b)
data4c <- data3c %>% sample_n(size=trunc(nrow(data3c)/1000000)*1000000,replace=T)
rm(data3c)

#divide the million samples into 100 equal subsamples
  
y <- c(1:100)
samplen <- data.frame("samplen" = y)
  
data4a <- cbind(data4a,samplen)
data4b <- cbind(data4b,samplen)
data4c <- cbind(data4c,samplen)
  
#For each sample, calculate the median/mean value 
  
data4a <- group_by(as.data.frame(data4a),STATEFIP,samplen)
data4a <- dplyr::summarise(data4a, mgr=median(rent),.groups="keep")

data4b <- group_by(as.data.frame(data4b),STATEFIP,samplen)
data4b <- dplyr::summarise(data4b, mhi=median(inc),.groups="keep")
  
data4c <- group_by(as.data.frame(data4c),STATEFIP,samplen)
data4c <- dplyr::summarise(data4c, mhv=median(val),.groups="keep")
  
data5 <- inner_join(data4a,data4b,by=c("STATEFIP","samplen"))
data5 <- inner_join(data5,data4c,by=c("STATEFIP","samplen")) 

data5 <- as.data.frame(ungroup(data5))
  
#Calculate the resulting 'Yost Index' following the existing method
#using ranks as we have done previously

data6 <- data5 %>%
    group_by(samplen) %>%
    arrange(mgr) %>%
    mutate(mgrrank=rank(-mgr))
  
data6 <- data6 %>%
    group_by(samplen) %>%
    arrange(mhi) %>%
    mutate(mhirank=rank(-mhi))
  
data6 <- data6 %>%
    group_by(samplen) %>%
    arrange(mhi) %>%
    mutate(mhvrank=rank(-mhv))
  
data6 <- arrange(data6,STATEFIP)
  
data6 <- data.frame(ungroup(data6))

#Conduct a factor analysis on each subsample

for (i in 1:100){

print(i)
  
data7 <- filter(data6, samplen==i) 

data8 <- dplyr::select(data7, mgrrank, mhirank, mhvrank)
  
  data8f <- fa(data8,fm="ml")
  
  data9 <- mutate(arrange(cbind(data7,data8f$scores),ML1), row=row_number())
  
  results <- rbind(results,data9)
}

#Get statistics and present results

results2 <- results %>% 
  group_by(STATEFIP) %>%
  summarise(meanrank = mean(row), sdrank = sd(row), minrank = min(row), maxrank = max(row),
            lci=quantile(row,probs=c(.025)),uci=quantile(row,probs=c(.975)))

results2$state <- to_factor(results2$STATEFIP)
results2 <- arrange(results2, meanrank)
results2 <- dplyr::select(results2, state, meanrank, sdrank, lci, uci, minrank, maxrank)
results2$meanrank <- round(results2$meanrank,2)
results2$sdrank <- round(results2$sdrank,2)
results2$lci <- round(results2$lci,0)
results2$uci <- round(results2$uci,0)
knitr::kable(results2)







##############
# PUMA level #
##############

results <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(results) <- c('STATEFIP','PUMA','mgr','mhi','mhv','mgrrank','mhirank','mhvrank',
                       'ML1','row')  

#sample each column separately in order to exclude NAs

data3a <- as.data.frame(filter(data3,!is.na(rent)))
data3b <- as.data.frame(filter(data3,!is.na(inc)))
data3c <- as.data.frame(filter(data3,!is.na(val)))

#take a sample from the non-missing values in each variable
#equal to the number of values, rounded to the nearest million

data4a <- data3a %>% sample_n(size=trunc(nrow(data3a)/1000000)*1000000,replace=T)
rm(data3a)
data4b <- data3b %>% sample_n(size=trunc(nrow(data3b)/1000000)*1000000,replace=T)
rm(data3b)
data4c <- data3c %>% sample_n(size=trunc(nrow(data3c)/1000000)*1000000,replace=T)
rm(data3c)

#divide the million samples into 100 equal subsamples

y <- c(1:100)
samplen <- data.frame("samplen" = y)

data4a <- cbind(data4a,samplen)
data4b <- cbind(data4b,samplen)
data4c <- cbind(data4c,samplen)

#For each sample, calculate the median/mean value 

data4a <- group_by(as.data.frame(data4a),STATEFIP,PUMA,samplen)
data4a <- dplyr::summarise(data4a, mgr=median(rent),.groups="keep")

data4b <- group_by(as.data.frame(data4b),STATEFIP,PUMA,samplen)
data4b <- dplyr::summarise(data4b, mhi=median(inc),.groups="keep")

data4c <- group_by(as.data.frame(data4c),STATEFIP,PUMA,samplen)
data4c <- dplyr::summarise(data4c, mhv=median(val),.groups="keep")

data5 <- inner_join(data4a,data4b,by=c("STATEFIP","PUMA","samplen"))
data5 <- inner_join(data5,data4c,by=c("STATEFIP","PUMA","samplen")) 

data5 <- as.data.frame(ungroup(data5))

#Calculate the resulting 'Yost Index' following the existing method
#using ranks as we have done previously

data6 <- data5 %>%
  group_by(samplen) %>%
  arrange(mgr) %>%
  mutate(mgrrank=rank(-mgr))

data6 <- data6 %>%
  group_by(samplen) %>%
  arrange(mhi) %>%
  mutate(mhirank=rank(-mhi))

data6 <- data6 %>%
  group_by(samplen) %>%
  arrange(mhi) %>%
  mutate(mhvrank=rank(-mhv))

data6 <- arrange(data6,STATEFIP)

data6 <- data.frame(ungroup(data6))

#Conduct a factor analysis on each subsample

for (i in 1:100){
  
  print(i)
  
  data7 <- filter(data6, samplen==i) 
  
  data8 <- dplyr::select(data7, mgrrank, mhirank, mhvrank)
  
  data8f <- fa(data8,fm="ml")
  
  data9 <- dplyr::mutate(arrange(cbind(data7,data8f$scores),ML1), row=row_number())
  
  results <- rbind(results,data9)
}

#Get statistics and present results

results2 <- results %>% 
  group_by(STATEFIP,PUMA) %>%
  dplyr::summarise(meanrank = mean(row), sdrank = sd(row), minrank = min(row), maxrank = max(row),
            lci=quantile(row,probs=c(.025)),uci=quantile(row,probs=c(.975)))

results2$state <- to_factor(results2$STATEFIP)
results2 <- arrange(results2, meanrank)
results2 <- dplyr::select(results2, state, PUMA, meanrank, sdrank, lci, uci, minrank, maxrank)
results2$meanrank <- round(results2$meanrank,2)
results2$sdrank <- round(results2$sdrank,2)
results2$lci <- round(results2$lci,0)
results2$uci <- round(results2$uci,0)
knitr::kable(results2)

#####################################################
# is there a way to simulate the block group level? #
# frame it as a 'downscaling' problem?              #
#####################################################

test <- data3 %>% sample_n(size=100,replace=T)

#try to add education
if (!require("ipumsr")) stop("Reading IPUMS data into R requires the ipumsr package. It can be installed using the following command: install.packages('ipumsr')")

ddi <- read_ipums_ddi("c:/yost2/pums/usa_00003.xml")
data <- read_ipums_micro(ddi)







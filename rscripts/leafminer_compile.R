#Chronologies for leafminer analysis. Eight datasets: ANP, Howland, PEF, BAS sites w/ PIRU and THOC for each

#Load packages
install.packages("dfoliatR")

library(tidyverse)
library(dplR)
library(dfoliatR)
library(gridExtra)


# Load ringwidth data -----------------------------------------------------
getwd()
setwd("C:/01_NETN/Forest_Health/R_Dev/leafminer_data")

#Acadia
aPIRU <- read.rwl("ANP_PIRU_MDI_Mesic6.rwl") 
aTHOC <- read.tucson("ANP_THOC_No_Uplands.raw")

#Howland
hPIRU <- read.tucson("HOW_PIRU_BEST.txt") 
hTHOC <- read.tucson("HOW_THOC_BEST.txt")

#PEF
pPIRU <- read.tucson("PEF_PIRU.raw") 
pTHOC <- read.tucson("PEF_THOC.raw")

#BAS
bPIRU <- read.tucson("BAS_PIRU.raw") 
bTHOC <- read.tucson("BAS_THOC.raw")

getwd()
setwd("C:/01_NETN/Forest_Health/R_Dev/leafminer")

#Combine dataframes into a list
nnames <- c("aPIRU", "aTHOC", "hPIRU", "hTHOC", "pPIRU", "pTHOC", "bPIRU", "bTHOC")
rwl_data <- list(aPIRU, aTHOC, hPIRU, hTHOC, pPIRU, pTHOC, bPIRU, bTHOC)
rwl_data <- rwl_data %>% set_names(nm = nnames) #set names

# Summary Stats -----------------------------------------------------------
spag_plots <- rwl_data %>% map(spag.plot) #show up on the right, but can't list is NULL
rwl_stats <- rwl_data %>% map(rwl.stats)

# Detrend -----------------------------------------------------------------
dtrend_F <- rwl_data %>% map(detrend, method = c("Friedman")) #for climate growth relationship
dtrend_M <- rwl_data %>% map(detrend, method = c("Mean")) #for stand dynamics

# Chronologies ------------------------------------------------------------
nnames2 <- nnames <- c("aPIRU_F", "aTHOC_F", "hPIRU_F", "hTHOC_F", "pPIRU_F", "pTHOC_F", "bPIRU_F", "bTHOC_F")
nnames3 <- nnames <- c("aPIRU_M", "aTHOC_M", "hPIRU_M", "hTHOC_M", "pPIRU_M", "pTHOC_M", "bPIRU_M", "bTHOC_M")
crns_F <- dtrend_F %>% map(chron) %>% set_names(nm = nnames2)
crns_M <- dtrend_M %>% map(chron) %>% set_names(nm = nnames3)

# Plot chronologies -------------------------------------------------------
list2env(crns_F, envir = .GlobalEnv)
list2env(crns_M, envir = .GlobalEnv)

print_crn<-function(df_chr){
  site<-deparse(substitute(df_chr))
  fig_name<-paste0("./chrono_graphs/", site, '.jpg')
  ppi<-300
  jpeg(file = fig_name, units='px', width=10*ppi, height=7*ppi, res=300) 
  df.plot<-plot.crn(df_chr, add.spline=TRUE, nyrs=20, lab = c(7,7,7))
  dev.off()
}# function to print a chronology to specified folder; lab used to specify number of axis tick marks

#print to chrono graphs folder
print_crn(aPIRU_F)
print_crn(aTHOC_F)
print_crn(bPIRU_F)
print_crn(bTHOC_F)
print_crn(hPIRU_F)
print_crn(hTHOC_F)
print_crn(pPIRU_F)
print_crn(pTHOC_F)

print_crn(aPIRU_M)
print_crn(aTHOC_M)
print_crn(bPIRU_M)
print_crn(bTHOC_M)
print_crn(hPIRU_M)
print_crn(hTHOC_M)
print_crn(pPIRU_M)
print_crn(pTHOC_M)

# Messing with dfoliatR ---------------------------------------------------
#need detrended .rwl for host, and single chrono for nonhost
pPIRU_F2 <- pPIRU_F %>% select(1)
aPIRU_F2 <- aPIRU_F %>% select(1)
hPIRU_F2 <- hPIRU_F %>% select(1)
bPIRU_F2 <- bPIRU_F %>% select(1)

#PEF
pLM_defol <- defoliate_trees(host_tree = dtrend_F[["pTHOC"]], nonhost_chron = pPIRU_F2, bridge_events = TRUE)# did not work
PEF_defolEvents <- plot_defol(pLM_defol)
#Howland
hLM_defol <- defoliate_trees(host_tree = dtrend_F[["hTHOC"]], nonhost_chron = hPIRU_F2, bridge_events = TRUE)# did not work
H_defolEvents <- plot_defol(hLM_defol)
#ANP
aLM_defol <- defoliate_trees(host_tree = dtrend_F[["aTHOC"]], nonhost_chron = aPIRU_F2, bridge_events = TRUE)# did not work
A_defolEvents <- plot_defol(aLM_defol)
#BAS
bLM_defol <- defoliate_trees(host_tree = dtrend_F[["bTHOC"]], nonhost_chron = bPIRU_F2, bridge_events = TRUE)# did not work
B_defolEvents <- plot_defol(bLM_defol)

#outbreak analysis
#PEF
pLM_outbr <- outbreak(pLM_defol)
plot_outbreak(pLM_outbr)
#howland
hLM_outbr <- outbreak(hLM_defol)
plot_outbreak(hLM_outbr)
#ANP
aLM_outbr <- outbreak(aLM_defol)
plot_outbreak(aLM_outbr)
#BAS
bLM_outbr <- outbreak(bLM_defol)
plot_outbreak(bLM_outbr)

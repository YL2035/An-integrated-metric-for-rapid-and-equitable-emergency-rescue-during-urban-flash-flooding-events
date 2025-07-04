library(ggplot2)
library(dplyr)
library(tidyr)
library(igraph)

###################
#### Load Data ####
###################
# with powerline
file_dir <- ".../edge-withpowerline/"
# list of file names
file_list <- list.files(path = file_dir, pattern = "*.csv")
# read in all files
for (i in 1:length(file_list)){
  assign(file_list[i], read.csv(paste0(file_dir,file_list[i])))
}

# The SVI building location is different, need to rejoin
# Using spatial join, multiple buildings will have the same building location, need to remove them
node_res_svi <- read.csv(".../node_res_svi.csv")
colnames(node_res_svi)[61] <- "RID" # resident building ID
node_res_svi_R <- filter(node_res_svi, Type == "R")
dupRID <- node_res_svi_R$RID[duplicated(node_res_svi_R$RID)]
node_res_svi_dup <- filter(node_res_svi_R, !(RID %in% dupRID))
colnames(node_res_svi_dup)
# filter the selected RES ID
res_select <- node_res_svi_dup$RID
# block and RID
BRID <- node_res_svi_dup[,c(61,30,31)]
BRID$BID <- paste0(BRID$Block_y,BRID$Tract)

################################
#### Metric part 1: Pdamage ####
################################
PBdamage <- node_res_svi_dup[,c(61,28,30,31)]
PBdamage$BldgDmg <- PBdamage$BldgDmgPct/100
colnames(node_res_svi_dup)

############################
#### Metric part 3: SVI ####
############################
socioeconomic <- unique(node_res_svi_dup[,c(30,31,33:36)])
V1 <- socioeconomic[,c(1,2,3)]
V1$P65Over <- rank(V1$P65Over)/length(V1$P65Over)
V2 <- socioeconomic[,c(1,2,4)]
V2$P18Under <- rank(V2$P18Under)/length(V2$P18Under)
V3 <- socioeconomic[,c(1,2,5)]
V3$PNonWhite <- rank(V3$PNonWhite)/length(V3$PNonWhite)
V4 <- socioeconomic[,c(1,2,6)]
V4$PSingleParent <- rank(V4$PSingleParent)/length(V4$PSingleParent)

socioeconomic_PR <- left_join(V1, V2, by = c("Block_y", "Tract"))
socioeconomic_PR <- left_join(socioeconomic_PR, V3, by = c("Block_y", "Tract"))
socioeconomic_PR <- left_join(socioeconomic_PR, V4, by = c("Block_y", "Tract"))

# Sum the percentile
socioeconomic_PR$socioeconomic_PR_sum <- socioeconomic_PR$P65Over + socioeconomic_PR$P18Under + socioeconomic_PR$PNonWhite + socioeconomic_PR$PSingleParent
SVI_sum <- socioeconomic_PR[,c(1,2,7)]
SVI_sum$SVI_PR_total <- rank(SVI_sum$socioeconomic_PR_sum)/length(SVI_sum$socioeconomic_PR_sum)
SVI_final <- SVI_sum[,-3]

#############################
#### Metric (Same Weight)####
#############################
# load the efficiency data
file_dir <- ".../efficiency-withpowerline/"

# list of file names
file_list <- list.files(path = file_dir, pattern = "*.csv")
# read in all files
for (i in 1:length(file_list)){
  assign(file_list[i], read.csv(paste0(file_dir,file_list[i])))
}

colnames(efficiency1.csv)

# load the data for the network analysis
for (i in 1:26){
  # get the data frame
  name <- paste0("efficiency",as.character(i),".csv")
  efficiency_i <- get(name)
  efficiency_i <- efficiency_i[,c(5,8,11)]
  
  # calculate the metric
  PBdamage_efficiency <- left_join(PBdamage, efficiency_i, by = "RID")
  # average building damage and accessibility by block
  PBdamage_efficiency_block_avg <- PBdamage_efficiency %>%
    group_by(Block_y, Tract, BID) %>% 
    summarise(Pdamage_avg = mean(BldgDmg), inv_dist_avg = mean(min_inv_dist))
  # percentile ranking for building damage and efficiency
  PBdamage_efficiency_block_avg$Pdamage_PR <- rank(PBdamage_efficiency_block_avg$Pdamage_avg)/length(PBdamage_efficiency_block_avg$Pdamage_avg)
  PBdamage_efficiency_block_avg$Efficiency_PR <- rank(-PBdamage_efficiency_block_avg$inv_dist_avg)/length(PBdamage_efficiency_block_avg$inv_dist_avg)
  
  metric <- left_join(PBdamage_efficiency_block_avg, SVI_final, by = c("Block_y", "Tract"))
  
  # Sum the percentile
  metric$metric_sum <- metric$Pdamage_PR/3 + metric$Efficiency_PR/3 + metric$SVI_PR_total/3
  metric$metric_PR <- rank(metric$metric_sum)/length(metric$metric_sum)
  
  # Conventional ranking
  metric$metric_wosvi_sum <- metric$Pdamage_PR/2 + metric$Efficiency_PR/2
  metric$metric_wosvi_PR <- rank(metric$metric_wosvi_sum)/length(metric$metric_wosvi_sum)
  write.csv(metric, paste0(".../","_TID_",as.character(i),".csv"))
}

######################
#### No Powerline ####
######################
# without powerline
efficiency_i <- read.csv(".../efficiency_nopowerline.csv")
efficiency_i <- efficiency_i[,c(5,8,11)]
PBdamage_efficiency <- left_join(PBdamage, efficiency_i, by = "RID")
# average building damage and accessibility by block
PBdamage_efficiency_block_avg <- PBdamage_efficiency %>%
  group_by(Block_y, Tract, BID) %>% 
  summarise(Pdamage_avg = mean(BldgDmg), inv_dist_avg = mean(min_inv_dist))
# percentile ranking for building damage and efficiency
PBdamage_efficiency_block_avg$Pdamage_PR <- rank(PBdamage_efficiency_block_avg$Pdamage_avg)/length(PBdamage_efficiency_block_avg$Pdamage_avg)
PBdamage_efficiency_block_avg$Efficiency_PR <- rank(-PBdamage_efficiency_block_avg$inv_dist_avg)/length(PBdamage_efficiency_block_avg$inv_dist_avg)
metric <- left_join(PBdamage_efficiency_block_avg, SVI_final, by = c("Block_y", "Tract"))

# Sum the percentile
metric$metric_sum <- metric$Pdamage_PR/3 + metric$Efficiency_PR/3 + metric$SVI_PR_total/3
metric$metric_PR <- rank(metric$metric_sum)/length(metric$metric_sum)
# Conventional ranking
metric$metric_wosvi_sum <- metric$Pdamage_PR/2 + metric$Efficiency_PR/2
metric$metric_wosvi_PR <- rank(metric$metric_wosvi_sum)/length(metric$metric_wosvi_sum)

##################################
#### Metric (Different Weight)####
##################################
weight <- read.csv(".../weights.csv")

# load the efficiency data
file_dir <- ".../efficiency-withpowerline/"
# list of file names
file_list <- list.files(path = file_dir, pattern = "*.csv")
# read in all files
for (i in 1:length(file_list)){
  assign(file_list[i], read.csv(paste0(file_dir,file_list[i])))
}
colnames(efficiency1.csv)

# load the data for the network analysis
for (w in 1:nrow(weight)){
  w1 <- weight[w,2]
  w2 <- weight[w,3]
  w3 <- weight[w,4]
  
  for (i in 1:26){
    # get the data frame
    name <- paste0("efficiency",as.character(i),".csv")
    efficiency_i <- get(name)
    efficiency_i <- efficiency_i[,c(5,8,11)]
    
    # calculate the metric
    PBdamage_efficiency <- left_join(PBdamage, efficiency_i, by = "RID")
    # average building damage and accessibility by block
    PBdamage_efficiency_block_avg <- PBdamage_efficiency %>%
      group_by(Block_y, Tract, BID) %>% 
      summarise(Pdamage_avg = mean(BldgDmg), inv_dist_avg = mean(min_inv_dist))
    # percentile ranking for building damage and efficiency
    PBdamage_efficiency_block_avg$Pdamage_PR <- rank(PBdamage_efficiency_block_avg$Pdamage_avg)/length(PBdamage_efficiency_block_avg$Pdamage_avg)
    PBdamage_efficiency_block_avg$Efficiency_PR <- rank(-PBdamage_efficiency_block_avg$inv_dist_avg)/length(PBdamage_efficiency_block_avg$inv_dist_avg)
    
    metric <- left_join(PBdamage_efficiency_block_avg, SVI_final, by = c("Block_y", "Tract"))
    
    # Sum the percentile
    metric$metric_sum <- metric$Pdamage_PR*w1 + metric$Efficiency_PR*w2 + metric$SVI_PR_total*w3
    metric$metric_PR <- rank(metric$metric_sum)/length(metric$metric_sum)
    write.csv(metric, paste0(".../","TID",as.character(i),"WID",as.character(w),".csv"))
  }
}

######################
#### No Powerline ####
######################
for (w in 1:nrow(weight)){
  w1 <- weight[w,2]
  w2 <- weight[w,3]
  w3 <- weight[w,4]
  
  # without powerline
  efficiency_i <- read.csv(".../efficiency_nopowerline.csv")
  efficiency_i <- efficiency_i[,c(5,8,11)]
  
  # calculate the metric
  PBdamage_efficiency <- left_join(PBdamage, efficiency_i, by = "RID")
  # average building damage and accessibility by block
  PBdamage_efficiency_block_avg <- PBdamage_efficiency %>%
    group_by(Block_y, Tract, BID) %>% 
    summarise(Pdamage_avg = mean(BldgDmg), inv_dist_avg = mean(min_inv_dist))
  
  # percentile ranking for building damage and efficiency
  PBdamage_efficiency_block_avg$Pdamage_PR <- rank(PBdamage_efficiency_block_avg$Pdamage_avg)/length(PBdamage_efficiency_block_avg$Pdamage_avg)
  PBdamage_efficiency_block_avg$Efficiency_PR <- rank(-PBdamage_efficiency_block_avg$inv_dist_avg)/length(PBdamage_efficiency_block_avg$inv_dist_avg)
  metric <- left_join(PBdamage_efficiency_block_avg, SVI_final, by = c("Block_y", "Tract"))
  
  # Sum the percentile
  metric$metric_sum <- metric$Pdamage_PR*w1 + metric$Efficiency_PR*w2 + metric$SVI_PR_total*w3
  metric$metric_PR <- rank(metric$metric_sum)/length(metric$metric_sum)
  
  write.csv(metric, paste0(".../","WID",as.character(w),".csv"))
}
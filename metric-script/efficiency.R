library(ggplot2)
library(dplyr)
library(tidyr)
library(igraph)

###################
#### Load Data ####
###################
# with powerline
file_dir <- "/edge-withpowerline.csv"
# list of file names
file_list <- list.files(path = file_dir, pattern = "*.csv")
# read in all files
for (i in 1:length(file_list)){
  assign(file_list[i], read.csv(paste0(file_dir,file_list[i])))
}

node_res_svi <- read.csv("/node_res_svi.csv")
colnames(node_res_svi)[61] <- "RID"
node_res_svi_R <- filter(node_res_svi, Type == "R")
dupRID <- node_res_svi_R$RID[duplicated(node_res_svi_R$RID)]
node_res_svi_dup <- filter(node_res_svi_R, !(RID %in% dupRID))
colnames(node_res_svi_dup)
# filter the selected RES ID
res_select <- node_res_svi_dup$RID
# block and RID
BRID <- node_res_svi_dup[,c(61,30,31)]
BRID$BID <- paste0(BRID$Block_y,BRID$Tract)

#####################################
#### Efficiency (with Powerline) ####
#####################################
#### H1 network efficiency ####
#### node (the same for all networks)
node <- data.frame(RID = node.csv$OID_, type = node.csv$Type)
node$RID <- as.character(node$RID)
unique(node$type)
facility_ID <- filter(node, type %in% c("E","F"))$RID
facility_ID <- facility_ID[-4]

# a function to calculate the efficiency to reach essential facilities
efficiency_to_facility <- function(edge_h1.csv){
  # edge
  edge <- data.frame(from = edge_h1.csv$FID_From, to = edge_h1.csv$FID_To)
  edge_sort <- data.frame(t(apply(edge[1:2], 1, sort))) # Simplify edge 
  edge_unique <- unique(edge_sort)
  
  # create a graph
  g <- graph_from_data_frame(d = edge_unique, vertices = node, directed = F)
  # get the shortest distance between all pairs of nodes
  dist <- distances(g)
  # get the dist for the selected resID, and facility ID
  dist_res_facility <- dist[as.character(res_select),facility_ID]
  
  # get the closet dist and facility ID
  facility_min <- colnames(dist_res_facility)[apply(dist_res_facility,1,which.min)]
  dist_min <- apply(dist_res_facility,1,min)
  
  # final dist table
  dist_res_facility <- data.frame(dist_res_facility)
  dist_res_facility$RID <- rownames(dist_res_facility)
  dist_res_facility$RID <- as.numeric(dist_res_facility$RID)
  dist_res_facility$min_dist <- dist_min
  dist_res_facility$min_FID <- facility_min
  dist_res_facility$min_inv_dist <- 1/dist_res_facility$min_dist
  dist_res_facility <- left_join(dist_res_facility, BRID, by = "RID")
  colnames(dist_res_facility)[1:3] <- paste0("dist_F",facility_ID)
  return(dist_res_facility)
}

# edge with powerline
for (i in 1:26){
  # get the dataframe
  name <- paste0("edge_h",as.character(i),".csv")
  print(name)
  edge <- get(name)
  
  # efficiency
  efficiency_hr <- efficiency_to_facility(edge)
  # write the file to a folder, the output data were provided in the file
  write.csv(efficiency_hr, paste0("/efficiency-withpowerline/efficiency",as.character(i),".csv"))
}

######################
#### No Powerline ####
######################
# without powerline
edge_noboat <- read.csv("/edge_intact.csv")

# efficiency
efficiency_noboat <- efficiency_to_facility(edge_noboat)
# write the file to a folder, the output data were provided in the file
write.csv(efficiency_noboat, "/withoutpowerline/efficiency_nopowerline.csv")



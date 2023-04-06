
########## Get the file location right ###############################
#Creat a new Rstript file and save it to where you store the network data#####
#In Rstudio, click "Session" and "Set Working Directory" to "To Source File Location"


########## Use the network data to build Network1 #################
library(igraph)
library(visNetwork)
library(htmlwidgets)

netdata<-read.csv("day5_socialdynamics_Networkdata.csv", header = TRUE) 
Network1<-graph_from_data_frame(netdata, directed = FALSE)


########## Question 1 #######################################################
nodes <- data.frame(id=V(Network1)$name)
edges <- data.frame(from = netdata[,1],
                    to = netdata[,2])

set.seed(100)
VisNetwork1_layout<-visNetwork(nodes, edges, main = "Network1",
                               submain="Can zoom in/out to check the IDs and ties") %>%
  visIgraphLayout(layout = "layout_nicely",smooth =  FALSE) %>%
  visNodes(shape="circle",label = TRUE) %>% 
  visOptions(highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T)

VisNetwork1_layout
#Check out the network and find the weak ties and strong ties
#You can zoom in and zoom out to check the links; 
#By clicking a node, its id will be shown on the small window of left hand side; and its nearby nodes will be highlighted
#Recall the definitions of weak ties and strong ties in the lecture


########## Question 2 ####################################################

########## 2.1 Functions to build the IC model ########## 
#Note1: You can develop your own code according to the model description in Q2; 
#Note2: If you decide to use the code below, try to understand.

stopifnot(require(data.table))
stopifnot(require(Matrix))

calculate_value <- function(node, each_neighbors,Pprob){
  return(each_neighbors[[node]][ which(runif(length(each_neighbors[[node]]), 0, 1)<=Pprob)])
  #'runif' is a function to generate random number in R
}
#This function:
#1) searchs the neighbours of contangious node; 
#2) To those who are connected to a contagious node, generates a random number and compare to the 
#probability of p, if random number<p, this node will be infected and return the value of 1

IC<-function(node_seed,network,Pprob){
  
  #prepare input for the 'calculate_value' function#
  adj_matrix <- igraph::as_adjacency_matrix(network, type = 'both')
  each_neighbors <- which(adj_matrix > 0, arr.ind = TRUE)
  each_neighbors <- split(each_neighbors[, 2], each_neighbors[, 1]) #get the neigbhour list of each node
  
  nNode<-vcount(network)
  node_status <- rep.int(0, nNode) #start from a healthy population
  day_infected<-vector()#Total number of infected population
  new_infected <- list()  # Record the ID of person getting infected at each time step
  
  day<-1
  node_status[as.numeric(node_seed)] <- 1 # infected(value=1) health(value=0)
  day_infected[day] <- sum(node_status ) 
  new_infected[[day]]<-node_seed #The ID of the person infected in Day 1 (Patient Zero)

  #simulate the spread of virus with 6 weeks##
  for (day in c(2:28)){  
    ContagiousID<-which(node_status == 1) 
    infectedID<-unlist(lapply(ContagiousID,calculate_value,each_neighbors,Pprob))
    newinfectedID<- setdiff(infectedID, which(node_status == 1))
    
   #Update the node status and other variables
    node_status[newinfectedID] <- 1
    day_infected[day] <- length(newinfectedID)
    new_infected[[day]]<-newinfectedID
    
    day=day+1
  }
  return(day_infected)  #return the number of newly infected people by day 
  #return(list(day_infected,new_infected)) #if you want to see the ID of infected ppl in each day, use this command
}


########## 2.2 Apply the IC model to Network 1 #############
result1<- replicate(100, IC(5,Network1,0.35), simplify=FALSE)  #run 100 times since each IC model run has it own randomness 
result1<-do.call(rbind, result1)
result1_forplot<-colMeans(result1)
result1_forplot   # the number of newly infected people by day 

########## Question 3 ####################################################
### 3.1 Reconstruct the network by deleting 3 weak ties###

#Change the numbers to the ties that you want to delete,
#e.g., if you want to delete the tie between 22 and 46, change "1|2" to "22|46"
Network2<-Network1 %>%
  delete_edges("1|2") %>%   #!!!Change to your own numbers
  delete_edges("3|4") %>%   #!!!Change to your own numbers
  delete_edges("5|6")       #!!!Change to your own numbers


### 3.2 Apply the IC model to Network 2###
result2<- replicate(100, IC(5,Network2,0.35), simplify=FALSE)  #run 100 times since each IC model run has it own randomness 
result2<-do.call(rbind, result2)
result2_forplot<-colMeans(result2)  # the number of newly infected people by day 
result2_forplot

########## Question 4 ####################################################
### 4.1 Reconstruct the network by deleting 3 strong ties###

#Change the numbers to the ties that you want to delete,
#e.g., if you want to delete the tie between 22 and 46, change "1|2" to "22|46"
Network3<-Network1 %>%
  delete_edges("1|2") %>%   #!!!Change to your own numbers
  delete_edges("3|4") %>%   #!!!Change to your own numbers
  delete_edges("5|6")       #!!!Change to your own numbers

### 4.2 Apply the IC model to Network 3###
result3<- replicate(100, IC(5,Network3,0.35), simplify=FALSE)  #run 100 times since each IC model run has it own randomness 
result3<-do.call(rbind, result3)
result3_forplot<-colMeans(result3)  # the number of newly infected people by day 
result3_forplot


########## Question 5 #############################################
SumResult<-data.frame(matrix(nrow=28*3,ncol=3))
names(SumResult)=c("Day","Network","Infected")
SumResult[,1]=rep(c(1:28),3)
SumResult[,2]=c(rep("Network1",28),rep("Network2",28),rep("Network3",28))
SumResult[,3]=c(result1_forplot,result2_forplot,result3_forplot)

library(ggplot2)
Plot_p0<-ggplot(SumResult, aes(x=Day, y=Infected, fill=Network)) +
  geom_line(aes(color=Network))+geom_point(aes(color=Network))+
  labs(title = paste( "Daily infection curve when p =", 0.35 ))+
  ylab("Newly confirmed cases each day")
Plot_p0



########## Question 6 #############################################
transitivity(Network1,type="global") #global clustering coefficient
transitivity(Network2,type="global")
transitivity(Network3,type="global")

diameter(Network1) 
diameter(Network2) 
diameter(Network3) 

mean_distance(Network1)  #the average shortest path lengths of Network 2 change noticeably
mean_distance(Network2)
mean_distance(Network3)



########## Question 7 #############################################
################# 7.1 Higher p ##########
p1=0.7 #Or replace it to other high value 

result1_p1<- replicate(100, IC(5,Network1,p1), simplify=FALSE)  
result1_p1<-do.call(rbind, result1_p1)
result1_p1_forplot<-colMeans(result1_p1)

result2_p1<- replicate(100, IC(5,Network2,p1), simplify=FALSE)  
result2_p1<-do.call(rbind, result2_p1)
result2_p1_forplot<-colMeans(result2_p1)

result3_p1<- replicate(100, IC(5,Network3,p1), simplify=FALSE)  #run 100 times since each IC model run has it own randomness 
result3_p1<-do.call(rbind, result3_p1)
result3_p1_forplot<-colMeans(result3_p1)

SumResult_p1<-data.frame(matrix(nrow=28*3,ncol=3))
names(SumResult_p1)=c("Day","Network","Infected")
SumResult_p1[,1]=rep(c(1:28),3)
SumResult_p1[,2]=c(rep("Network1",28),rep("Network2",28),rep("Network3",28))
SumResult_p1[,3]=c(result1_p1_forplot,result2_p1_forplot,result3_p1_forplot)

Plot_p1<-ggplot(SumResult_p1, aes(x=Day, y=Infected, fill=Network)) +
  geom_line(aes(color=Network))+geom_point(aes(color=Network))+
  labs(title = paste( "Daily infection curve when p =", p1 ))+
  ylab("Newly confirmed cases each day")
Plot_p1



################# 7.2 Lower p ##########
p2=0.1 #Or replace it to other low value 

result1_p2<- replicate(100, IC(5,Network1,p2), simplify=FALSE)  
result1_p2<-do.call(rbind, result1_p2)
result1_p2_forplot<-colMeans(result1_p2)

result2_p2<- replicate(100, IC(5,Network2,p2), simplify=FALSE)  
result2_p2<-do.call(rbind, result2_p2)
result2_p2_forplot<-colMeans(result2_p2)

result3_p2<- replicate(100, IC(5,Network3,p2), simplify=FALSE)  #run 100 times since each IC model run has it own randomness 
result3_p2<-do.call(rbind, result3_p2)
result3_p2_forplot<-colMeans(result3_p2)

SumResult_p2<-data.frame(matrix(nrow=28*3,ncol=3))
names(SumResult_p2)=c("Day","Network","Infected")
SumResult_p2[,1]=rep(c(1:28),3)
SumResult_p2[,2]=c(rep("Network1",28),rep("Network2",28),rep("Network3",28))
SumResult_p2[,3]=c(result1_p2_forplot,result2_p2_forplot,result3_p2_forplot)

Plot_p2<-ggplot(SumResult_p2, aes(x=Day, y=Infected, fill=Network)) +
  geom_line(aes(color=Network))+geom_point(aes(color=Network))+
  labs(title = paste( "Daily infection curve when p =", p2 ))+
  ylab("Newly confirmed cases each day")
Plot_p2








################################# Afternoon Section ########################################################################

###################################### Question 8 ##################################################

##################### 8.1 Functions to build the threshold model ##########################
stopifnot(require(data.table))
stopifnot(require(Matrix))


calculate_adoptedNei <- function(node, node_status, each_neighbors){
  return(mean(node_status[each_neighbors[[node]]] == 1)) ### to calculate the percentage of adopted neigbhours
}


ThModel<-function(node_seed,network,threshold){ 
  #prepare input for the 'calculate_value' function#
  adj_matrix <- igraph::as_adjacency_matrix(network, type = 'both')
  each_neighbors <- which(adj_matrix > 0, arr.ind = TRUE)
  each_neighbors <- split(each_neighbors[, 2], each_neighbors[, 1]) #get the neigbhour list of each node
  
  nNode<-vcount(network)
  node_status <- rep.int(0, nNode) 
  neighbour_status<-rep.int(0, nNode)  ##percentage of adopted neighbours
  new_infected <- list()
  day_total_infected <- rep(0,28) ### Total number of active people by end of each day
  
  
  ### Day 1 ####
  day <- 1
  node_status[as.numeric((node_seed))] <- 1 
  new_infected[[day]] <-node_seed
  day_total_infected[day]=sum(node_status == 1)
  
  ####
  
  for (day in c(2:28)){
    NotAdopted <- which(node_status == 0)
    Adopted <- which(node_status == 1)
    
    neighbour_status[NotAdopted] <- unlist(lapply(NotAdopted, calculate_adoptedNei,
                                            node_status, each_neighbors))
    
    new_infected[[day]] <- setdiff(which(neighbour_status > threshold), Adopted)
    node_status[new_infected[[day]]] <- 1  #update the staus to 1 for those newly adopted
    day_total_infected[day] <- sum(node_status)
    
    day <- day + 1
  }
  #return(day_total_infected)
  return(list(day_total_infected,new_infected))
}


##################### 8.2 read threshold data ##########################
#since we don't have the 'real' threshold for the population, we will assign threshold by random number from 0 to 1
set.seed(10)  #don't change the seed number
threshold<-runif(64,min=0,max=1) 

##################### 8.3 Apply the threshold model to Network 1##########################
seeds<-c(1,2,5,11) #the persons who first join the campaign
ThM_Network1<-ThModel(seeds,Network1,threshold)
ThM_Network1[[1]]  #total numer of active people by day


#### Question 9 Visualize the seed status and discuss why network-wide congation is not happening #######
ThM_result1=rep("Inactive",64) 
ThM_result1[unlist(ThM_Network1[[2]])]<-"Active"
ThM_result1[seeds]<-"Seeds"

nodes <- data.frame(id=V(Network1)$name,group=ThM_result1)
set.seed(100)
ThM_status1<-visNetwork(nodes, edges, main = "Outcome by seeding Node 1, 2, 5",
                        submain="Seeds in red, Active nodes in orange, Inactive nodes in grey") %>%
  visIgraphLayout(layout = "layout_nicely",smooth =  FALSE) %>%
  visGroups(groupname = "Seeds", color ="#D62728FF", shape = "circle",label = TRUE) %>%
  visGroups(groupname = "Active", color = "#FF7F0EFF", shape = "circle",label = TRUE) %>%
  visGroups(groupname = "Inactive", color = "#8A9197FF", shape = "circle",label = TRUE)%>% 
  visOptions(highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T)


ThM_status1 #Inactive nodes in grey, seed nodes in red, active ndoes in orange



#################### Question 10 high-degree seeding ###################################
degree(Network1) #degree by nodes
nSeed=4 #can hire 4 ambassadors


#Find out 4 high-Degree seeds and use them as seeds: "seedHD"


ThM_NetworkHD<-ThModel(seedHD,Network1,threshold)  #seedHD is the high-degree nodes that you want to seed
ThM_NetworkHD[[1]]  #total numer of active people by day

#### Question 9 Visualize the seed status and discuss why network-wide congation is not happening #######
ThM_resultHD=rep("Inactive",64) 
ThM_resultHD[unlist(ThM_NetworkHD[[2]])]<-"Active"
ThM_resultHD[seedHD]<-"Seeds"

nodes <- data.frame(id=V(Network1)$name,group=ThM_resultHD)
set.seed(100)
ThM_statusHD<-visNetwork(nodes, edges, main = "Outcome of high-degree seeding",
                         submain="Seeds in red, Active nodes in orange, Inactive nodes in grey") %>%
  visIgraphLayout(layout = "layout_nicely",smooth =  FALSE) %>%
  visGroups(groupname = "Seeds", color ="#D62728FF", shape = "circle",label = TRUE) %>%
  visGroups(groupname = "Active", color = "#FF7F0EFF", shape = "circle",label = TRUE) %>%
  visGroups(groupname = "Inactive", color = "#8A9197FF", shape = "circle",label = TRUE)%>% 
  visOptions(highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T)

ThM_statusHD  #Outcome of high-degree seeding



#################### Question 11 Greedy algorithm ###################################

greedy_ThM<-function(network,threshold,k){#if you want to select 5 seeds, k=5
  nNode<-vcount(network)
  SeedCan<-c(1:nNode)  # nodes that can still be used as seeds
  Seedset<-c()  # nodes that will be used as seeds
  
  for (s in 1:k){
    best_seed<--1  # initial setting of each round
    best_spread <- -Inf
    
    N<-length(SeedCan)
    for (i in 1:N){
      current_seed<-SeedCan[i]
      current_seedset<-c(Seedset,current_seed)
      current_spread<-ThModel(current_seedset,network,threshold)[[1]]
      if (sum(current_spread)>sum(best_spread)){
        best_seed<-current_seed
        best_spread<-current_spread
      }
    }
    Seedset<-c(Seedset,best_seed)
    SeedCan<-setdiff(SeedCan,Seedset) #exclude the seed nodes in next round seed selection
  }
  return(Seedset)
}

###### 11.1 the four ambassadors suggested by the greedy algorithm#########
greedy_ThM4<-greedy_ThM(Network1,threshold,4) #return the seeds suggested by greedy algorithm


###### 11.2 seeding outcome of the greedy algorithm############
ThModel(greedy_ThM4,Network1,threshold)[[1]] #Diffusion outcome of the seeds suggested by greedy algorithm


######visualize the seeding outcome##########
ThM_resultGreedy=rep("Inactive",64) 
ThM_resultGreedy[unlist(ThModel(greedy_ThM4,Network1,threshold)[[2]])]<-"Active"
ThM_resultGreedy[greedy_ThM4]<-"Seeds"

nodes <- data.frame(id=V(Network1)$name,group=ThM_resultGreedy)
set.seed(100)
ThM_statusGreedy<-visNetwork(nodes, edges, main = "Seeding outcome of greedy algorithm",
                             submain="Seeds in red, Active nodes in orange, Inactive nodes in grey") %>%
  visIgraphLayout(layout = "layout_nicely",smooth =  FALSE) %>%
  visGroups(groupname = "Seeds", color ="#D62728FF", shape = "circle",label = TRUE) %>%
  visGroups(groupname = "Active", color = "#FF7F0EFF", shape = "circle",label = TRUE) %>%
  visGroups(groupname = "Inactive", color = "#8A9197FF", shape = "circle",label = TRUE)%>% 
  visOptions(highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T)

ThM_statusGreedy


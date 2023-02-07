library(network)
library(igraph)
library(iGraphMatch)
library(dplyr)
library(mclust) 
library(irlba)
library(multiRDPG)
library(ggplot2)
library(foreach)

source("functions/functions.R")
source("functions/mase.R")
source("functions/parametric-bootstrap.R")
source("functions/sample_graphs.R")
source("functions/jointembedding.R")

array_split <- function(data, m) {
  rowIdx <- seq_len(nrow(data))    
  lapply(split(rowIdx,cut(rowIdx, m)), function(x) data[x, ])
}

# example
#brain.dat.1 <- read_graph("https://www.cis.jhu.edu/~parky/Microsoft/JHU-MSR/ZMx2/BNU1/DS01216/subj1-scan1.graphml",format="graphml")
#brain.dat.2<- read_graph("https://www.cis.jhu.edu/~parky/Microsoft/JHU-MSR/ZMx2/BNU1/DS01216/subj1-scan2.graphml",format="graphml")

# read data---------
subjs <- c(1:57)
scans <- c(1:2)
f <- function(i, n, e) paste0("https://www.cis.jhu.edu/~parky/Microsoft/JHU-MSR/ZMx2/BNU1/DS01216/subj", i,"-scan",n, e)
urls <- sapply(".graphml", function(.e) outer(subjs, scans, f, e = .e))

listOfFiles <- lapply(urls, function(x) read_graph(x, format = "graphml")) 
namelist <- lapply(urls, function(x) substr(x, 68, 79)) 

# list of graphs----
listOfGraphs <- lapply(listOfFiles, function(x) 
  {delete_vertices(x, V(x)[hemisphere == "none"]) })

# remove none tissues
listOfGraphs2 <- lapply(listOfGraphs, function(x) 
{delete_vertices(x, V(x)[tissue == "none"]) })


table(vertex_attr(test,"region"), vertex_attr(test,"hemisphere"))

# list of adjs----
brain_adjs <- lapply(listOfGraphs2, function (x) as_adjacency_matrix(x, sparse = F)) 
names(brain_adjs) <- namelist

# check number of vertex----
number_of_v <- sapply(brain_adjs, function(x) dim(x)[1])
number_of_v <- sapply(listOfGraphs, function(x) length(V(x)))

number_of_v[1:57]-number_of_v[58:114]
# check region 
number_of_region <- sapply(listOfGraphs, function(x) length(unique(vertex_attr(x, "region"))))
# check tissue
number_of_tissue <-  sapply(listOfGraphs, function(x) length(unique(vertex_attr(x, "tissue"))))
# check hemisphere
num_of_hemi <- sapply(listOfGraphs, function(x) length(unique(vertex_attr(x, "hemisphere"))))

# get vertex list
# vertexnames_list <- list()
# vertexnames_list <- lapply(1:length(brain_adjs), function(i){
#   colnames(brain_adjs[[i]])
# })

# vertex names for all graphs
vertexnames_list <- lapply(1:length(listOfGraphs2), function(i) (
  V(listOfGraphs2[[i]])$name
))

# choose graphs---- to embed together and get common vertex----
commonnumber <- c(1:114)
length(Reduce(intersect, vertexnames_list[commonnumber]))
commonvertex <- Reduce(intersect, vertexnames_list[commonnumber])

# common subgraphs----
brain_commongraphs <- list()
for (i in 1:length(commonnumber)) {
  brain_commongraphs[[i]] <- induced_subgraph(listOfGraphs[[commonnumber[i]]], vids=commonvertex)
} # not in order yet

# real region
realregion <- vertex_attr(brain_commongraphs[[1]],"region")
# real hemisphere
realhemisphere <- vertex_attr(brain_commongraphs[[1]],"hemisphere")
# real tissue
realtissue <-vertex_attr(brain_commongraphs[[1]],"tissue")

# get common adjs----
brain_commonadjs <- list()
for (k in 1:length(commonnumber)) {
  brain_commonadjs[[k]] <- as_adjacency_matrix(brain_commongraphs[[k]],sparse = F)
  brain_commonadjs[[k]] <- brain_commonadjs[[k]][commonvertex,commonvertex]
}
#check elbow----
d<-getElbows(irlba(brain_commonadjs[[1]],120)$d)[2]
d

m = length(commonnumber) # number of graphs
n = dim(brain_commonadjs[[3]])[1] # number of vertex
i = 26 # elbow dim to embed

# mclust num----
Mnum = 4  #2 # 70 # 4
reallabels = realtissue
# calculate preds 

#calculate_preds <- function(m=5,i=2) {
  # Omni----
  OmniM<-buildOmni(brain_commonadjs)
  Xomni <- ase2(OmniM,d=i)$X
  # get average and prediction

  OmniXlist <-array_split(data=Xomni, m = m)
  
  AvgOmni <-Reduce("+", OmniXlist) / m
  preds <- Mclust(AvgOmni,Mnum) 
  Omni <- preds$classification
  OmniARI<- adjustedRandIndex(Omnipreds,reallabels) 
  
  #------------MASE----
  scaled_MASE <- mase(brain_commonadjs, d=i)
  V <- scaled_MASE$V
  # R list
  R_list <- scaled_MASE$R
  # sig list 'Eigen-decomp of R_1, R_2'
  sigs_list <- list()
  for (k in 1:m) {
    sigs_list[[k]] <-  matrix(0, nrow = i, ncol =i)
    diag(sigs_list[[k]]) <- eigen(R_list[[k]])$values
  }
    
  'Procrustes'
  new_eig_vec_lsit <- list()
  for (k in 1:(m-1)) {
    W <- procrustes(eigen(R_list[[k]])$vectors, eigen(R_list[[m]])$vectors)
    new_eig_vec_lsit[[k]] <- eigen(R_list[[k]])$vectors%*%W$W
  }
  new_eig_vec_lsit[[m]] <-eigen(R_list[[m]])$vectors
  'Estimated latent positions'
  Xs <- list()
  for (k in 1:m) {
    Xs[[k]]  <- V%*%new_eig_vec_lsit[[k]]%*%(abs(sigs_list[[k]])^{1/2})
  }

  # get average and prediction
  AvgMase <-Reduce("+", Xs) / m
  preds <- Mclust(AvgMase,Mnum)
  Mase <- preds$classification
  MaseARI<- adjustedRandIndex(Masepreds,reallabels) 
  
  #-----------MRDPG----
  scaled_MRDPG <- multiRDPG(brain_commonadjs, d=i)
  U <- scaled_MRDPG$U
  'Estimated latent positions'
  Xs <- list()
  for (k in 1:m) {
    Xs[[k]] <- scaled_MRDPG$Lambda[[k]]
    Xs[[k]] <- U%*%Xs[[k]]^{1/2}
  }
  # get average
  AvgMRDPG <- Reduce("+", Xs) / m
  preds <- Mclust(AvgMRDPG,Mnum)
  MRDPG <- preds$classification
  MRDPGARI<- adjustedRandIndex(MRDPGpreds,reallabels) 
  
  #-----------joint embedding----
  jointemb <- multidembed(brain_commonadjs,d=i)
  jointemb <- jointemb$h
  preds <- Mclust(jointemb,Mnum)
  jointemb <- preds$classification
  jointembARI<- adjustedRandIndex(jointembpreds,reallabels) 
  
  #-----------separate ASE, procrustes, average ----
  Xs <- list()
  for (k in 1:m) {
    Xs[[k]] <- ase2(brain_commonadjs[[k]],d=i)$X
  }
  "Procrustes"
  Xs_new <- list()
  for (k in 1:(m-1)) {
    W <- procrustes(Xs[[k]],Xs[[m]])
    Xs_new[[k]] <- Xs[[k]]%*%W$W
  }
  Xs_new[[m]] <- Xs[[m]]

  # get average
  AvgProc <- Reduce("+", Xs_new) / m
  preds <- Mclust(AvgProc,Mnum)
  Proc <- preds$classification
  ProcARI<- adjustedRandIndex(Procpreds,reallabels) 
  
  #-----------average first, then embed----
  avgadj <- Reduce("+", brain_commonadjs) / m
  ASEavg <- ase2(avgadj,d=i)$X
  preds <- Mclust(ASEavg,Mnum)
  Avg<- preds$classification
  AvgARI <- adjustedRandIndex(Avgpreds,reallabels) 
  
  #print(c(OmniARI, MaseARI, MRDPGARI,jointembARI,ProcARI,AvgARI))
  #return(c(OmniARI, MaseARI, MRDPGARI,jointembARI,ProcARI,AvgARI))
#}
  #results[nrow(results) + 1,] <- c(OmniARI, MaseARI, MRDPGARI,jointembARI,ProcARI,AvgARI)
  
  #row.names(results)[5] <- "tissue_1_25"
#result <- data.frame(OmniARI, MaseARI, MRDPGARI,jointembARI,ProcARI,AvgARI)

# add real on heatmap data ----
  heatmapdata <- matrix(rep(0,147),ncol = 3)
  colnames(heatmapdata) <- c("colname","rowname","value")
  heatmapdata[,1]<-rep(c("Omni","Mase","MRDPG","jointemb","Proc","Avg","realhemisphere"),7)
  heatmapdata[,2]<- c(rep("Omni",7),rep("Mase",7),rep("MRDPG",7),rep("jointemb",7),
                      rep("Proc",7),rep("Avg",7),rep("realhemisphere",7))
  
  for (i in 1:nrow(heatmapdata)) {
    heatmapdata[i,3]<- adjustedRandIndex(get(heatmapdata[i,1]),get(heatmapdata[i,2]))
  }
  heatmapdata <- as.data.frame(heatmapdata)
  heatmapdata$value <- as.numeric(heatmapdata$value)

  #return(heatmapdata)

# plot ----
ggplot(heatmapdata, aes(x = rowname, y = colname, fill = value)) +
  geom_tile(color="black")+
  geom_text(aes(label = round(value, 2)),colour = "darkgoldenrod1")+
  scale_fill_gradient(low="white", high="darkgreen",limits= c(0,1)) +
  xlab("methods")+
  ylab("methods")+
  theme_bw()+
  ggtitle("Brain data 1-20 with real hemisphere, d=26")


# SAME algo different dim ----
  d = seq(2,26, 1)
  namesvec <- c()
  for(i in 2:26) { 
    namesvec[i-1] <- paste("dim", i, sep = "_")
  }
  namesvec[length(d)+1] <- "realtissue" #8888888
  
  Mnum = NULL #4 # 70 # 2
  reallabels = realtissue #8888888
  
  ## Omni----
  preds_Omni = list()
  #for (i in d) {
  foreach(i = 2:26) %do%{
    OmniM<-buildOmni(brain_commonadjs)
    Xomni <- ase2(OmniM,d=i)$X
    # get average and prediction
    OmniXlist <-array_split(data=Xomni, m = m)
    
    AvgOmni <-Reduce("+", OmniXlist) / m
    preds <- Mclust(AvgOmni,Mnum) 
    Omni <- preds$classification
    
    preds_Omni[[i-1]]<- Omni
  }
  preds_Omni[[length(d)+1]] <- realtissue #8888888
  names(preds_Omni) <- namesvec
  ## heatmap Omni ----
  heatmap_Omni <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_Omni) <- c("colname","rowname","value")
  heatmap_Omni[,1] <- rep(namesvec,26)
  heatmap_Omni[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_Omni)) {
    heatmap_Omni[i,3]<- adjustedRandIndex(preds_Omni[[heatmap_Omni[i,1]]],preds_Omni[[heatmap_Omni[i,2]]])
  }
  
  ## MASE----
  preds_MASE = list()
  #for (i in d) {
  foreach(i = 2:26) %do%{
    scaled_MASE <- mase(brain_commonadjs, d=i)
    V <- scaled_MASE$V
    # R list
    R_list <- scaled_MASE$R
    # sig list 'Eigen-decomp of R_1, R_2'
    sigs_list <- list()
    for (k in 1:m) {
      sigs_list[[k]] <-  matrix(0, nrow = i, ncol =i)
      diag(sigs_list[[k]]) <- eigen(R_list[[k]])$values
    }
    'Procrustes'
    new_eig_vec_lsit <- list()
    for (k in 1:(m-1)) {
      W <- procrustes(eigen(R_list[[k]])$vectors, eigen(R_list[[m]])$vectors)
      new_eig_vec_lsit[[k]] <- eigen(R_list[[k]])$vectors%*%W$W
    }
    new_eig_vec_lsit[[m]] <-eigen(R_list[[m]])$vectors
    'Estimated latent positions'
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]]  <- V%*%new_eig_vec_lsit[[k]]%*%(abs(sigs_list[[k]])^{1/2})
    }
    # get average and prediction
    AvgMase <-Reduce("+", Xs) / m
    preds <- Mclust(AvgMase,Mnum)
    Mase <- preds$classification
    
    preds_MASE[[i-1]]<- Mase
  }
  preds_MASE[[length(d)+1]] <- realtissue #8888888
  names(preds_MASE) <- namesvec
  ## heatmap mase
  heatmap_MASE <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_MASE) <- c("colname","rowname","value")
  heatmap_MASE[,1] <- rep(namesvec,26)
  heatmap_MASE[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_MASE)) {
    heatmap_MASE[i,3]<- adjustedRandIndex(preds_MASE[[heatmap_MASE[i,1]]],preds_MASE[[heatmap_MASE[i,2]]])
  }
  
  
  ## MRDPG----
  preds_MRDPG = list()
  #for (i in d) {
  foreach(i = 2:26) %do%{
    scaled_MRDPG <- multiRDPG(brain_commonadjs, d=i)
    U <- scaled_MRDPG$U
    'Estimated latent positions'
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]] <- scaled_MRDPG$Lambda[[k]]
      Xs[[k]] <- U%*%Xs[[k]]^{1/2}
    }
    # get average
    AvgMRDPG <- Reduce("+", Xs) / m
    preds <- Mclust(AvgMRDPG,Mnum)
    MRDPG <- preds$classification
    
    preds_MRDPG[[i-1]]<- MRDPG
  }
  preds_MRDPG[[length(d)+1]] <- realtissue #8888888
  names(preds_MRDPG) <- namesvec
  ## heatmap MRDPG
  heatmap_MRDPG <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_MRDPG) <- c("colname","rowname","value")
  heatmap_MRDPG[,1] <- rep(namesvec,26)
  heatmap_MRDPG[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_MRDPG)) {
    heatmap_MRDPG[i,3]<- adjustedRandIndex(preds_MRDPG[[heatmap_MRDPG[i,1]]],preds_MRDPG[[heatmap_MRDPG[i,2]]])
  }
  
  ## Joint embedding ----
  preds_Joint = list()
  foreach(i = 2:26) %do%{
    jointemb <- multidembed(brain_commonadjs,d=i)
    jointemb <- jointemb$h
    preds <- Mclust(jointemb,Mnum)
    jointemb <- preds$classification
    
    preds_Joint[[i-1]] <- jointemb
  }
  preds_Joint[[length(d)+1]] <- realtissue #8888888
  names(preds_Joint) <- namesvec
  ## heatmap Joint
  heatmap_Joint <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_Joint) <- c("colname","rowname","value")
  heatmap_Joint[,1] <- rep(namesvec,26)
  heatmap_Joint[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_Joint)) {
    heatmap_Joint[i,3]<- adjustedRandIndex(preds_Joint[[heatmap_Joint[i,1]]],preds_Joint[[heatmap_Joint[i,2]]])
  }
  
  
  ## separate Procruste----
  preds_Proc = list()
  foreach(i = 2:26) %do%{
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]] <- ase2(brain_commonadjs[[k]],d=i)$X
    }
    "Procrustes"
    Xs_new <- list()
    for (k in 1:(m-1)) {
      W <- procrustes(Xs[[k]],Xs[[m]])
      Xs_new[[k]] <- Xs[[k]]%*%W$W
    }
    Xs_new[[m]] <- Xs[[m]]
    
    # get average
    AvgProc <- Reduce("+", Xs_new) / m
    preds <- Mclust(AvgProc,Mnum)
    Proc <- preds$classification
    
    preds_Proc[[i-1]] <- Proc
  }
  preds_Proc[[length(d)+1]] <- realtissue #8888888
  names(preds_Proc) <- namesvec
  ## heatmap Proc
  heatmap_Proc <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_Proc) <- c("colname","rowname","value")
  heatmap_Proc[,1] <- rep(namesvec,26)
  heatmap_Proc[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_Proc)) {
    heatmap_Proc[i,3]<- adjustedRandIndex(preds_Proc[[heatmap_Proc[i,1]]],preds_Proc[[heatmap_Proc[i,2]]])
  } 
  
  ## AVG ----
  preds_Avg = list()
  foreach(i = 2:26) %do%{
    avgadj <- Reduce("+", brain_commonadjs) / m
    ASEavg <- ase2(avgadj,d=i)$X
    preds <- Mclust(ASEavg,Mnum)
    Avg<- preds$classification
    
    preds_Avg[[i-1]] = Avg
  }
  preds_Avg[[length(d)+1]] <- realtissue #8888888
  names(preds_Avg) <- namesvec
  ## heatmap Avg
  heatmap_Avg <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_Avg) <- c("colname","rowname","value")
  heatmap_Avg[,1] <- rep(namesvec,26)
  heatmap_Avg[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_Avg)) {
    heatmap_Avg[i,3]<- adjustedRandIndex(preds_Avg[[heatmap_Avg[i,1]]],preds_Avg[[heatmap_Avg[i,2]]])
  }

  heatmaps <- c("heatmap_Omni","heatmap_MASE","heatmap_MRDPG",
               "heatmap_Joint","heatmap_Avg","heatmap_Proc")
  results <- list()
  plots <- list()
  for (j in heatmaps) {
    df = get(j)
    df = data.frame(df)
    df$value = as.numeric(df$value)
    order = 1:nrow(df)
    p = ggplot(df, aes(x = reorder(rowname,order), y = reorder(colname,order), fill = value)) +
      geom_tile(color="black")+
      #geom_text(aes(label = round(value, 2)),colour = "darkgoldenrod1")+
      scale_fill_gradient(low="white", high="darkgreen",limits= c(0,1)) +
      xlab("dimensions")+
      ylab("dimensions")+
      theme_bw()+
      ggtitle(paste(j,"with different embedding dimensions"))
    best <- df %>% filter(colname == "realtissue") %>% arrange(value) %>%  #8888888
      summarise(maxvalue = max(value[value != max(value)])) 
    best = df[which(df$value==best$maxvalue),]
    plots[[which(heatmaps == j)]] <- p
    results[[which(heatmaps == j)]] <- best
  }
  names(results) <- heatmaps
  
  Omni=preds_Omni[['dim_2']]
  Mase=preds_MASE[['dim_2']]
  MRDPG=preds_MRDPG[['dim_2']]
  jointemb=preds_Joint[['dim_2']]
  Proc=preds_Proc[['dim_2']]
  Avg=preds_Avg[['dim_2']]
  
  
  #--- get heatmap data list ----
  heatmapdata <- matrix(rep(0,147),ncol = 3)
  colnames(heatmapdata) <- c("colname","rowname","value")
  heatmapdata[,1]<-rep(c("Omni","Mase","MRDPG","jointemb","Proc","Avg","realtissue"),7) #8888888
  heatmapdata[,2]<- c(rep("Omni",7),rep("Mase",7),rep("MRDPG",7),rep("jointemb",7),
                      rep("Proc",7),rep("Avg",7),rep("realtissue",7)) #8888888
  
  for (i in 1:nrow(heatmapdata)) {
    heatmapdata[i,3]<- adjustedRandIndex(get(heatmapdata[i,1]),get(heatmapdata[i,2]))
  }
  #heatmapdata <- as.data.frame(heatmapdata)
  heatmapdata = data.frame(heatmapdata) %>% mutate(dim=i)
  heatmapdata$value <- as.numeric(heatmapdata$value)
  ggplot(heatmapdata, aes(x = rowname, y = colname, fill = value)) +
    geom_tile(color="black")+
    scale_fill_gradient(low="white", high="darkgreen",limits= c(0,1)) +
    xlab("dimensions")+
    ylab("dimensions")+
    ggtitle("Brain data tissue")

  
  # KMEAN algo different dim ----
  d = seq(2,26, 1)
  namesvec <- c()
  for(i in 2:26) { 
    namesvec[i-1] <- paste("dim", i, sep = "_")
  }
  namesvec[length(d)+1] <- "realtissue" #8888888
  
  Mnum =  4 # 70 # 2
  reallabels = realtissue #8888888
  
  ## Omni----
  preds_Omni = list()
  #for (i in d) {
  foreach(i = 2:26) %do%{
    OmniM<-buildOmni(brain_commonadjs)
    Xomni <- ase2(OmniM,d=i)$X
    # get average and prediction
    OmniXlist <-array_split(data=Xomni, m = m)
    
    AvgOmni <-Reduce("+", OmniXlist) / m
    preds <- kmeans(AvgOmni,Mnum) 
    Omni <- preds$cluster
    
    preds_Omni[[i-1]]<- Omni
  }
  preds_Omni[[length(d)+1]] <- realtissue #8888888
  names(preds_Omni) <- namesvec
  ## heatmap Omni
  heatmap_Omni <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_Omni) <- c("colname","rowname","value")
  heatmap_Omni[,1] <- rep(namesvec,26)
  heatmap_Omni[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_Omni)) {
    heatmap_Omni[i,3]<- adjustedRandIndex(preds_Omni[[heatmap_Omni[i,1]]],preds_Omni[[heatmap_Omni[i,2]]])
  }
  
  ## MASE----
  preds_MASE = list()
  #for (i in d) {
  foreach(i = 2:26) %do%{
    scaled_MASE <- mase(brain_commonadjs, d=i)
    V <- scaled_MASE$V
    # R list
    R_list <- scaled_MASE$R
    # sig list 'Eigen-decomp of R_1, R_2'
    sigs_list <- list()
    for (k in 1:m) {
      sigs_list[[k]] <-  matrix(0, nrow = i, ncol =i)
      diag(sigs_list[[k]]) <- eigen(R_list[[k]])$values
    }
    'Procrustes'
    new_eig_vec_lsit <- list()
    for (k in 1:(m-1)) {
      W <- procrustes(eigen(R_list[[k]])$vectors, eigen(R_list[[m]])$vectors)
      new_eig_vec_lsit[[k]] <- eigen(R_list[[k]])$vectors%*%W$W
    }
    new_eig_vec_lsit[[m]] <-eigen(R_list[[m]])$vectors
    'Estimated latent positions'
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]]  <- V%*%new_eig_vec_lsit[[k]]%*%(abs(sigs_list[[k]])^{1/2})
    }
    # get average and prediction
    AvgMase <-Reduce("+", Xs) / m
    preds <- kmeans(AvgMase,Mnum)
    Mase <- preds$cluster
    
    preds_MASE[[i-1]]<- Mase
  }
  preds_MASE[[length(d)+1]] <- realtissue #8888888
  names(preds_MASE) <- namesvec
  ## heatmap mase
  heatmap_MASE <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_MASE) <- c("colname","rowname","value")
  heatmap_MASE[,1] <- rep(namesvec,26)
  heatmap_MASE[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_MASE)) {
    heatmap_MASE[i,3]<- adjustedRandIndex(preds_MASE[[heatmap_MASE[i,1]]],preds_MASE[[heatmap_MASE[i,2]]])
  }
  
  
  ## MRDPG----
  preds_MRDPG = list()
  #for (i in d) {
  foreach(i = 2:26) %do%{
    scaled_MRDPG <- multiRDPG(brain_commonadjs, d=i)
    U <- scaled_MRDPG$U
    'Estimated latent positions'
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]] <- scaled_MRDPG$Lambda[[k]]
      Xs[[k]] <- U%*%Xs[[k]]^{1/2}
    }
    # get average
    AvgMRDPG <- Reduce("+", Xs) / m
    preds <- kmeans(AvgMRDPG,Mnum)
    MRDPG <- preds$cluster
    
    preds_MRDPG[[i-1]]<- MRDPG
  }
  preds_MRDPG[[length(d)+1]] <- realtissue #8888888
  names(preds_MRDPG) <- namesvec
  ## heatmap MRDPG
  heatmap_MRDPG <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_MRDPG) <- c("colname","rowname","value")
  heatmap_MRDPG[,1] <- rep(namesvec,26)
  heatmap_MRDPG[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_MRDPG)) {
    heatmap_MRDPG[i,3]<- adjustedRandIndex(preds_MRDPG[[heatmap_MRDPG[i,1]]],preds_MRDPG[[heatmap_MRDPG[i,2]]])
  }
  
  ## Joint embedding ----
  preds_Joint = list()
  foreach(i = 2:26) %do%{
    jointemb <- multidembed(brain_commonadjs,d=i)
    jointemb <- jointemb$h
    preds <- kmeans(jointemb,Mnum)
    jointemb <- preds$cluster
    
    preds_Joint[[i-1]] <- jointemb
  }
  preds_Joint[[length(d)+1]] <- realtissue #8888888
  names(preds_Joint) <- namesvec
  ## heatmap Joint
  heatmap_Joint <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_Joint) <- c("colname","rowname","value")
  heatmap_Joint[,1] <- rep(namesvec,26)
  heatmap_Joint[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_Joint)) {
    heatmap_Joint[i,3]<- adjustedRandIndex(preds_Joint[[heatmap_Joint[i,1]]],preds_Joint[[heatmap_Joint[i,2]]])
  }
  
  
  ## separate Procruste----
  preds_Proc = list()
  foreach(i = 2:26) %do%{
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]] <- ase2(brain_commonadjs[[k]],d=i)$X
    }
    "Procrustes"
    Xs_new <- list()
    for (k in 1:(m-1)) {
      W <- procrustes(Xs[[k]],Xs[[m]])
      Xs_new[[k]] <- Xs[[k]]%*%W$W
    }
    Xs_new[[m]] <- Xs[[m]]
    
    # get average
    AvgProc <- Reduce("+", Xs_new) / m
    preds <- kmeans(AvgProc,Mnum)
    Proc <- preds$cluster
    
    preds_Proc[[i-1]] <- Proc
  }
  preds_Proc[[length(d)+1]] <- realtissue #8888888
  names(preds_Proc) <- namesvec
  ## heatmap Proc
  heatmap_Proc <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_Proc) <- c("colname","rowname","value")
  heatmap_Proc[,1] <- rep(namesvec,26)
  heatmap_Proc[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_Proc)) {
    heatmap_Proc[i,3]<- adjustedRandIndex(preds_Proc[[heatmap_Proc[i,1]]],preds_Proc[[heatmap_Proc[i,2]]])
  } 
  
  ## AVG ----
  preds_Avg = list()
  foreach(i = 2:26) %do%{
    avgadj <- Reduce("+", brain_commonadjs) / m
    ASEavg <- ase2(avgadj,d=i)$X
    preds <- kmeans(ASEavg,Mnum)
    Avg<- preds$cluster
    
    preds_Avg[[i-1]] = Avg
  }
  preds_Avg[[length(d)+1]] <- realtissue #8888888
  names(preds_Avg) <- namesvec
  ## heatmap Avg
  heatmap_Avg <- matrix(rep(0,2028),ncol = 3)
  colnames(heatmap_Avg) <- c("colname","rowname","value")
  heatmap_Avg[,1] <- rep(namesvec,26)
  heatmap_Avg[,2] <- rep(namesvec,each = 26)
  for (i in 1:nrow(heatmap_Avg)) {
    heatmap_Avg[i,3]<- adjustedRandIndex(preds_Avg[[heatmap_Avg[i,1]]],preds_Avg[[heatmap_Avg[i,2]]])
  }
  
  heatmaps <- c("heatmap_Omni","heatmap_MASE","heatmap_MRDPG",
                "heatmap_Joint","heatmap_Avg","heatmap_Proc")
  results <- list()
  plots <- list()
  for (j in heatmaps) {
    df = get(j)
    df = data.frame(df)
    df$value = as.numeric(df$value)
    order = 1:nrow(df)
    p = ggplot(df, aes(x = reorder(rowname,order), y = reorder(colname,order), fill = value)) +
      geom_tile(color="black")+
      #geom_text(aes(label = round(value, 2)),colour = "darkgoldenrod1")+
      scale_fill_gradient(low="white", high="darkgreen",limits= c(0,1)) +
      xlab("dimensions")+
      ylab("dimensions")+
      theme_bw()+
      ggtitle(paste(j,"with different embedding dimensions"))
    best <- df %>% filter(colname == "realtissue") %>% arrange(value) %>%  #8888888
      summarise(maxvalue = max(value[value != max(value)])) 
    best = df[which(df$value==best$maxvalue),]
    plots[[which(heatmaps == j)]] <- p
    results[[which(heatmaps == j)]] <- best
  }
  names(results) <- heatmaps
  
  Omni=preds_Omni[['dim_2']]
  Mase=preds_MASE[['dim_2']]
  MRDPG=preds_MRDPG[['dim_2']]
  jointemb=preds_Joint[['dim_2']]
  Proc=preds_Proc[['dim_2']]
  Avg=preds_Avg[['dim_2']]
  
  
  #--- get heatmap data list ----
  heatmapdata <- matrix(rep(0,147),ncol = 3)
  colnames(heatmapdata) <- c("colname","rowname","value")
  heatmapdata[,1]<-rep(c("Omni","Mase","MRDPG","jointemb","Proc","Avg","realtissue"),7) #8888888
  heatmapdata[,2]<- c(rep("Omni",7),rep("Mase",7),rep("MRDPG",7),rep("jointemb",7),
                      rep("Proc",7),rep("Avg",7),rep("realtissue",7)) #8888888
  
  for (i in 1:nrow(heatmapdata)) {
    heatmapdata[i,3]<- adjustedRandIndex(get(heatmapdata[i,1]),get(heatmapdata[i,2]))
  }
  #heatmapdata <- as.data.frame(heatmapdata)
  heatmapdata = data.frame(heatmapdata) %>% mutate(dim=i)
  heatmapdata$value <- as.numeric(heatmapdata$value)
  ggplot(heatmapdata, aes(x = rowname, y = colname, fill = value)) +
    geom_tile(color="black")+
    scale_fill_gradient(low="white", high="darkgreen",limits= c(0,1)) +
    xlab("dimensions")+
    ylab("dimensions")+
    ggtitle("Brain data tissue")
  
  
 
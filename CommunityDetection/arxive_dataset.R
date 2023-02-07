# arxiv dataset from https://manliodedomenico.com/data.php
library(network)
library(igraph)
library(iGraphMatch)
library(dplyr)
library(mclust) 
library(irlba)
library(multiRDPG)
library(ggplot2)

setwd("/Users/tongqi/Desktop/papersandbooks")
source("codes/functions.R")
source("repos/mase/R/mase.R")
source("repos/mase/R/parametric-bootstrap.R")
source("repos/mase/R/sample_graphs.R")
source("codes/Algcompare/jointembedding.R")

arxiv<- read.table("./codes/Algcompare/arXiv-Netscience_Multiplex_Coauthorship/Dataset/arxiv_netscience_multiplex.edges", header=TRUE,
                   sep=" ")
colnames(arxiv) = c("layerID","nodeID","nodeID2","weight")
layers <- unique(arxiv$layerID)

arxivs <- list()
middles <- list()
arxiv_adjs <- list()
dims<-data.frame()
for (i in layers) {
  arxivs[[i]] <- subset(arxiv,layerID==i)
  arxivs[[i]] <- arxivs[[i]][,c(2,3)]
  middles[[i]] <- dplyr::inner_join(arxivs[[i]], arxivs[[i]], by = "nodeID")[,-1]
  middles[[i]] <- apply(middles[[i]], 2, as.character)
  arxiv_adjs[[i]] <- network(middles[[i]], directed = FALSE)
  arxiv_adjs[[i]] <- as.sociomatrix(arxiv_adjs[[i]])
  # arxiv_adjs[[i]] <- data.matrix(arxiv_adjs[[i]])
  print(dim(arxiv_adjs[[i]]))
  dims <- rbind(dims, dim(arxiv_adjs[[i]]))
}

colnames_list <- list()
colnames_list <- lapply(1:length(layers), function(i){
  colnames(arxiv_adjs[[i]])
})
# graphs ----
commonnumber <- c(2,3,6,12)
length(Reduce(intersect, colnames_list[commonnumber]))
commonvertex <- Reduce(intersect, colnames_list[commonnumber])

arxiv_commonadjs <- list()
for (i in 1:length(commonnumber)) {
  arxiv_commonadjs[[i]] <- arxiv_adjs[[commonnumber[i]]][commonvertex,commonvertex]
}

d<-getElbows(irlba(arxiv_commonadjs[[1]],30)$d)[2]
d

n = length(commonvertex) # number of vertex
m = length(commonnumber) # number of graphs
i=2 # number of dimension to embed

# only two graphs
getdatas2 <- function(i=2) {
  OmniM<-buildOmni(arxiv_commonadjs[c(1,2)])
  Xomni <- ase(OmniM,d=i)
  # get average and prediction
  AvgOmni <- (Xomni[1:n,]+Xomni[(n+1):(2*n),])/m
  preds <- Mclust(AvgOmni)
  Omnipreds <- preds$classification
  #OmniARI<- adjustedRandIndex(Omnipreds,real_class)
  
  #------------MASE
  scaled_MASE <- mase(arxiv_commonadjs[c(1,2)], d=i)
  V <- scaled_MASE$V
  R1 <- scaled_MASE$R[[1]]
  R2 <- scaled_MASE$R[[2]]
  'Eigen-decomp of R_1, R_2'
  sig1 <- matrix(0, nrow = i, ncol =i)
  diag(sig1) <- eigen(R1)$values
  sig2 <- matrix(0, nrow = i, ncol =i)
  diag(sig2) <- eigen(R2)$values
  'Procrustes'
  eig_vec1 <-  eigen(R1)$vectors
  eig_vec2 <-  eigen(R2)$vectors
  W<-procrustes(eig_vec1,eig_vec2)
  eig_vec1_new<-eig_vec1%*%W$W
  'Estimated latent positions'
  X1 <- V%*%eig_vec1_new%*%(abs(sig1)^{1/2})
  X2 <- V%*%eig_vec2%*%(abs(sig2)^{1/2})
  # get average and prediction
  AvgMase <- (X1+X2)/2
  preds <- Mclust(AvgMase)
  Masepreds <- preds$classification
  #MaseARI<- adjustedRandIndex(Masepreds,real_class)
  
  #-----------MRDPG
  scaled_MRDPG <- multiRDPG(arxiv_commonadjs[c(1,2)], d=i)
  U <- scaled_MRDPG$U
  L1 <- scaled_MRDPG$Lambda[[1]]
  L2 <- scaled_MRDPG$Lambda[[2]]
  'Estimated latent positions'
  X1 <- U%*%L1^{1/2}
  X2 <- U%*%L2^{1/2}
  # get average
  AvgMRDPG <- (X1+X2)/2
  preds <- Mclust(AvgMRDPG)
  MRDPGpreds <- preds$classification
  #MRDPGARI<- adjustedRandIndex(MRDPGpreds,real_class)
  
  #-----------joint embedding
  jointemb <- multidembed(arxiv_commonadjs[c(1,2)],d=i)
  jointemb <- jointemb$h
  preds <- Mclust(jointemb)
  jointembpreds <- preds$classification
  #jointembARI<- adjustedRandIndex(jointembpreds,real_class)
  
  #-----------separate ASE, procrustes, average 
  X1 <- ase(arxiv_commonadjs[[1]],d=i)
  X2 <- ase(arxiv_commonadjs[[2]],d=i)
  "Procrustes"
  W<-procrustes(X1,X2)
  X1_new<-X1%*%W$W
  # get average
  AvgProc <- (X1_new+X2)/2
  preds <- Mclust(AvgProc)
  Procpreds <- preds$classification
  #ProcARI<- adjustedRandIndex(Procpreds,real_class)
  
  #-----------average first, then embed
  avgadj <- (arxiv_commonadjs[[1]]+arxiv_commonadjs[[2]])/2
  ASEavg <- ase(avgadj,d=i)
  preds <- Mclust(ASEavg)
  Avgpreds<- preds$classification
  #AvgARI <- adjustedRandIndex(Avgpreds,real_class)
  
  heatmatrix <- matrix(rep(0,36), ncol = 6)
  colnames(heatmatrix) <- c("Omnipreds","Masepreds","MRDPGpreds","jointembpreds","Procpreds","Avgpreds")
  rownames(heatmatrix) <- c("Omnipreds","Masepreds","MRDPGpreds","jointembpreds","Procpreds","Avgpreds")
  
  heatmapdata <- matrix(rep(0,108),ncol = 3)
  colnames(heatmapdata) <- c("colname","rowname","value")
  heatmapdata[,1]<-rep(c("Omnipreds","Masepreds","MRDPGpreds","jointembpreds","Procpreds","Avgpreds"),6)
  heatmapdata[,2]<- c(rep("Omnipreds",6),rep("Masepreds",6),rep("MRDPGpreds",6),rep("jointembpreds",6),
                      rep("Procpreds",6),rep("Avgpreds",6))
  
  for (i in 1:nrow(heatmapdata)) {
    heatmapdata[i,3]<- adjustedRandIndex(get(heatmapdata[i,1]),get(heatmapdata[i,2]))
  }
  heatmapdata <- as.data.frame(heatmapdata)
  heatmapdata$value <- as.numeric(heatmapdata$value)
  return(heatmapdata)
}
# plot ----
ggplot(heatmapdata, aes(x = rowname, y = colname, fill = value)) +
  geom_tile(color="black")+
  geom_text(aes(label = round(value, 2)),colour = "darkgoldenrod1")+
  scale_fill_gradient(low="white", high="darkgreen",limits= c(0.2,1)) +
  xlab("methods")+
  ylab("methods")+
  ggtitle("arxiv_23_dim9")

# m}
#--------------------------------------------

# all the graphs---- 
m=length(commonnumber)

getdatas4<- function(i=2) {
    OmniM<-buildOmni(arxiv_commonadjs)
    Xomni <- ase(OmniM,d=i)
    # get average and prediction
    AvgOmni <- (Xomni[1:n,]+Xomni[(n+1):(2*n),])/m
    preds <- Mclust(AvgOmni)
    Omni <- preds$classification
    #OmniARI<- adjustedRandIndex(Omnipreds,real_class)
    
    #------------MASE
    scaled_MASE <- mase(arxiv_commonadjs, d=i)
    V <- scaled_MASE$V
    R1 <- scaled_MASE$R[[1]]
    R2 <- scaled_MASE$R[[2]]
    R3 <- scaled_MASE$R[[3]]
    R4 <- scaled_MASE$R[[4]]
    'Eigen-decomp of R_1, R_2'####
    sig1 <- matrix(0, nrow = i, ncol =i)
    diag(sig1) <- eigen(R1)$values
    sig2 <- matrix(0, nrow = i, ncol = i)
    diag(sig2) <- eigen(R2)$values
    sig3 <- matrix(0, nrow = i, ncol = i)
    diag(sig3) <- eigen(R3)$values
    sig4 <- matrix(0, nrow = i, ncol = i)
    diag(sig4) <- eigen(R4)$values
    'Procrustes'
    eig_vec1 <-  eigen(R1)$vectors
    eig_vec2 <-  eigen(R2)$vectors
    eig_vec3 <-  eigen(R3)$vectors
    eig_vec4 <-  eigen(R4)$vectors
    W1<-procrustes(eig_vec1,eig_vec4)
    W2<-procrustes(eig_vec2,eig_vec4)
    W3<-procrustes(eig_vec3,eig_vec4)
    eig_vec1_new<-eig_vec1%*%W1$W
    eig_vec2_new<-eig_vec2%*%W2$W
    eig_vec3_new<-eig_vec3%*%W3$W
    'Estimated latent positions'
    X1 <- V%*%eig_vec1_new%*%(abs(sig1)^{1/2})
    X2 <- V%*%eig_vec2_new%*%(abs(sig2)^{1/2})
    X3 <- V%*%eig_vec3_new%*%(abs(sig3)^{1/2})
    X4 <- V%*%eig_vec4%*%(abs(sig4)^{1/2})
    # get average and prediction
    AvgMase <- (X1+X2+X3+X4)/m
    preds <- Mclust(AvgMase)
    Mase <- preds$classification
    #MaseARI<- adjustedRandIndex(Masepreds,real_class)
    
    #-----------MRDPG
    scaled_MRDPG <- multiRDPG(arxiv_commonadjs, d=i)
    U <- scaled_MRDPG$U
    L1 <- scaled_MRDPG$Lambda[[1]]
    L2 <- scaled_MRDPG$Lambda[[2]]
    L3 <- scaled_MRDPG$Lambda[[3]]
    L4 <- scaled_MRDPG$Lambda[[4]]
    'Estimated latent positions'
    X1 <- U%*%L1^{1/2}
    X2 <- U%*%L2^{1/2}
    X3 <- U%*%L3^{1/2}
    X4 <- U%*%L4^{1/2}
    # get average
    AvgMRDPG <- (X1+X2+X3+X4)/m
    preds <- Mclust(AvgMRDPG)
    MRDPG <- preds$classification
    #MRDPGARI<- adjustedRandIndex(MRDPGpreds,real_class)
    
    #-----------joint embedding
    jointemb <- multidembed(arxiv_commonadjs,d=i)
    jointemb <- jointemb$h
    preds <- Mclust(jointemb)
    jointemb <- preds$classification
    #jointembARI<- adjustedRandIndex(jointembpreds,real_class)
    
    #-----------separate ASE, procrustes, average 
    X1 <- ase(arxiv_commonadjs[[1]],d=i)
    X2 <- ase(arxiv_commonadjs[[2]],d=i)
    X3 <- ase(arxiv_commonadjs[[3]],d=i)
    X4 <- ase(arxiv_commonadjs[[4]],d=i)
    "Procrustes"
    W1<-procrustes(X1,X2)
    W3<-procrustes(X3,X2)
    W4<-procrustes(X4,X2)
    X1_new<-X1%*%W1$W
    X3_new<-X3%*%W3$W
    X4_new<-X4%*%W4$W
    # get average
    AvgProc <- (X1_new+X2+X3_new+X4_new)/m
    preds <- Mclust(AvgProc)
    Proc <- preds$classification
    #ProcARI<- adjustedRandIndex(Procpreds,real_class)
    
    #-----------average first, then embed
    avgadj <- (arxiv_commonadjs[[1]]+arxiv_commonadjs[[2]]+
                 arxiv_commonadjs[[3]]+arxiv_commonadjs[[4]])/m
    ASEavg <- ase(avgadj,d=i)
    preds <- Mclust(ASEavg)
    Avg<- preds$classification
    #AvgARI <- adjustedRandIndex(Avgpreds,real_class)
    
    heatmapdata2 <- matrix(rep(0,108),ncol = 3)
    colnames(heatmapdata2) <- c("colname","rowname","value")
    heatmapdata2[,1]<-rep(c("Omni","Mase","MRDPG","jointemb","Proc","Avg"),6)
    heatmapdata2[,2]<- c(rep("Omni",6),rep("Mase",6),rep("MRDPG",6),rep("jointemb",6),
                         rep("Proc",6),rep("Avg",6))
    
    for (i in 1:nrow(heatmapdata2)) {
      heatmapdata2[i,3]<- adjustedRandIndex(get(heatmapdata2[i,1]),get(heatmapdata2[i,2]))
    }
    
    heatmapdata2 <- data.frame(heatmapdata2)
    heatmapdata2$value <- as.numeric(heatmapdata2$value)
    p = ggplot(heatmapdata2, aes(x = rowname, y = colname, fill = value)) +
      geom_tile(color="black")+
      geom_text(aes(label = round(value, 2)),colour = "darkgoldenrod1")+
      scale_fill_gradient(low="white", high="darkgreen", limits= c(0.25,1)) +
      xlab("methods")+
      ylab("methods")+
      ggtitle("arxiv_23612_dim2")
    
    print(p)
    
    return(heatmapdata2)
}

library(viridis)
library(RColorBrewer)
heatmapdata2_2 <- getdatas4(i=2)

# plot ----
p1= ggplot(heatmapdata2_2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile(color="black")+
  geom_text(aes(label = round(value, 2)),colour = "darkgoldenrod1")+
  scale_fill_gradient(low="white", high="darkgreen", limits= c(0.25,1)) +
  xlab(" ")+
  ylab("Methods")+
  ggtitle("arxiv_4graphs_dim2")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.5)) 
  #scale_fill_distiller(palette = "Blues",direction = 1)+
  

heatmapdata2_4 <- getdatas4(i=4)

# plot ----
p2= ggplot(heatmapdata2_4, aes(x = rowname, y = colname, fill = value)) +
  geom_tile(color="black")+
  geom_text(aes(label = round(value, 2)),colour = "darkgoldenrod1")+
  scale_fill_gradient(low="white", high="darkgreen", limits= c(0.25,1)) +
  xlab("Methods")+
  ylab(" ")+
  ggtitle("arxiv_4graphs_dim4")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))
  

heatmapdata2_9 <- getdatas4(i=9)

# plot ----
p3= ggplot(heatmapdata2_9, aes(x = rowname, y = colname, fill = value)) +
  geom_tile(color="black")+
  geom_text(aes(label = round(value, 2)),colour = "darkgoldenrod1")+
  scale_fill_gradient(low="white", high="darkgreen", limits= c(0.25,1)) +
  xlab(" ")+
  ylab(" ")+
  ggtitle("arxiv_4graphs_dim9")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))
  




library(ggpubr)
arxiv = ggarrange(
  p1,p2,p3,
  ncol = 3, nrow = 1, #labels = c("A","B"),
  common.legend = TRUE, legend = "right"
) 
arxiv
  
  





  
  
  colname <-c("Omnipreds","MRDPGpreds","jointembpreds","Procpreds","Avgpreds")
  lapply(colname, function(x){
    data.frame(value = get(x)) %>% mutate(variable = x)-> df
    return(df)
  }) %>% bind_rows() -> data


  lapply(colname, function(cvar) {
    lapply(colname, function(nvar){

      data %>% filter(variable %in% cvar) %>%
        rename(M_cal = variable, value1 = value) %>%
        bind_cols(
          data %>% filter(variable %in% nvar) %>%
            rename(M_row = variable, value2 = value))

    }) %>% bind_rows()
  }) %>% bind_rows() %>%
    group_by(M_cal, M_row) %>%
    summarise(value = adjustedRandIndex(value1, value2), .groups = "drop") ->
    heatmapdata

  ggplot(heatmapdata, aes(x = M_row, y = M_cal, fill = value)) +
    geom_tile(color="black")+
    scale_fill_gradient(low="white", high="darkgreen") +
    xlab("methods")+
    ylab("methods")
  

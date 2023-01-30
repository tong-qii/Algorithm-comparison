setwd("/Users/tongqi/Desktop/papersandbooks")
source("codes/functions.R")
source("repos/mase/R/mase.R")
source("repos/mase/R/parametric-bootstrap.R")
source("repos/mase/R/sample_graphs.R")
source("codes/Algcompare/jointembedding.R")
library(igraph)
library(iGraphMatch)
library(irlba)
library(multiRDPG)
library(mclust) 
library(dplyr)
library(ggplot2)

# Laplancian function ----
makeLaplancian <- function(A){
  r<- rowSums(A)
  smartroot<-function(y){if(y>0){y=y^(-1/2); return(y)}else {return(y)}}
  rr<-sapply(r,smartroot)
  D <- diag(rr)
  L <- D%*%A%*%D
}

# # generate SBM
# P = matrix(rep(0.2,4), ncol = 2)
# epsilons = seq(0.001,0.35, by=0.03)
# qr = matrix(c(0.7,0.1,0.1,0.3),ncol = 2)
# 
# B = matrix(c(0.9,0.2,0.2,0.5), ncol = 2)
# n = 200 # for each block
# m = 2 # number of graphs
# sbm <- sample_sbm(n,pref.matrix=B, block.sizes=c(n/2,n/2))
# sudo_classes <- c(rep(1,50),rep(2,50))
# # Z <- matrix(0, n, 2)
# # Z[1:(n/2), 1] = Z[(n/2 + 1):n,2] = 1
# # P1 = Z %*% B %*% t(Z)
# # P2 = Z %*% B %*% t(Z)


# adj matrix----
ARIs <- list()
df = list()
epsilons = seq(0.05,0.8, by=0.08)
n = 120 # for total number of vertex
m = 2 # number of graphs
blocks = 2 # number of blocks
sudo_classes <- c(rep(1,60),rep(2,60))
#change the value of epsilon----
for (j in 1:length(epsilons)) {
  print(epsilons[j])
  #P = matrix(c(0.8,0.5,0.5,0.8), ncol = 2)
  P = matrix(rep(0.1,4), ncol = 2)
  #P = matrix(rep(0.1,9), ncol = 3)
  #P = matrix(rep(0.1,16), ncol = 4)
  qr = matrix(c(0.9,0.01,0.01,0.3),ncol = 2)
  #qr = matrix(c(0.9,0.01,0.01,0.01,0.5,0.01,0.01,0.01,0.3),ncol = 3)
  #qr = matrix(rep(0.01,16),ncol = 4)
  #diag(qr)<- c(0.9,0.6,0.4,0.2)
  B = P
  diag(B) <- diag(B)-epsilons[i]
  
  S = 20
  seeds = 1000:(1000-1+S)
  ARIs[[j]] <- sapply(1:S, function(i) {
    set.seed(seeds[i])
    graphs_adjs = list()
    graphs_laps = list()
    for (k in 1:m) {
      graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=c(n/2,n/2))
      #graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=c(n/3,n/3,n/3))
      #graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=c(n/4,n/4,n/4,n/4))
      graphs_adjs[[k]] = as_adjacency_matrix(graphs_adjs[[k]], sparse = F)
      graphs_laps[[k]] = makeLaplancian(graphs_adjs[[k]])
    }

    #-----------Omni----
    OmniM<-buildOmni(graphs_adjs)
    Xomni <- ase2(OmniM,m)$X
    # get average and prediction
    AvgOmni <- (Xomni[1:n,]+Xomni[(n+1):(2*n),])/m
    preds <- Mclust(AvgOmni)
    Omnipreds <- preds$classification
    OmniARI<- adjustedRandIndex(Omnipreds,sudo_classes)
    
    #------------MASE----
    scaled_MASE <- mase(graphs_adjs, d=m)
    V <- scaled_MASE$V
    R1 <- scaled_MASE$R[[1]]
    R2 <- scaled_MASE$R[[2]]
    'Eigen-decomp of R_1, R_2'
    sig1 <- cbind( c( eigen(R1)$values[1],0), c(0, eigen(R1)$values[2]) )
    sig2 <- cbind( c( eigen(R2)$values[1],0), c(0, eigen(R2)$values[2]) )
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
    MaseARI<- adjustedRandIndex(Masepreds,sudo_classes)
    
    #-----------MRDPG----
    scaled_MRDPG <- multiRDPG(graphs_adjs, d=m)
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
    MRDPGARI<- adjustedRandIndex(MRDPGpreds,sudo_classes)
    
    #-----------joint embedding----
    jointemb <- multidembed(graphs_adjs,d=m)
    jointemb <- jointemb$h
    preds <- Mclust(jointemb)
    jointembpreds <- preds$classification
    jointembARI<- adjustedRandIndex(jointembpreds,sudo_classes)
    
    #-----------separate ASE, procrustes, average ----
    X1 <- ase2(graphs_adjs[[1]],d=m)$X
    X2 <- ase2(graphs_adjs[[2]],d=m)$X
    "Procrustes"
    W<-procrustes(X1,X2)
    X1_new<-X1%*%W$W
    # get average
    AvgProc <- (X1_new+X2)/2
    preds <- Mclust(AvgProc)
    Procpreds <- preds$classification
    ProcARI<- adjustedRandIndex(Procpreds,sudo_classes)
    
    #-----------average first, then embed----
    avgadj <- (graphs_adjs[[1]]+graphs_adjs[[2]])/2
    ASEavg <- ase2(avgadj,d=m)$X
    preds <- Mclust(ASEavg)
    Avgpreds<- preds$classification
    AvgARI <- adjustedRandIndex(Avgpreds,sudo_classes)
    

    #-----------Omni_lap----
    OmniM<-buildOmni(graphs_laps)
    Xomni <- ase2(OmniM,m)$X
    # get average and prediction
    AvgOmni <- (Xomni[1:n,]+Xomni[(n+1):(2*n),])/m
    preds <- Mclust(AvgOmni)
    Omnipreds <- preds$classification
    Omni_lapARI<- adjustedRandIndex(Omnipreds,sudo_classes)

    #MASE_lap----
    scaled_MASE <- mase(graphs_laps, d=m)
    V <- scaled_MASE$V
    R1 <- scaled_MASE$R[[1]]
    R2 <- scaled_MASE$R[[2]]
    'Eigen-decomp of R_1, R_2'
    sig1 <- cbind( c( eigen(R1)$values[1],0), c(0, eigen(R1)$values[2]) )
    sig2 <- cbind( c( eigen(R2)$values[1],0), c(0, eigen(R2)$values[2]) )
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
    Mase_lapARI<- adjustedRandIndex(Masepreds,sudo_classes)


    # #-----------jointembedding_lap----
    # jointemb <- multidembed(graphs_laps,d=m)
    # jointemb <- jointemb$h
    # preds <- Mclust(jointemb)
    # jointembpreds <- preds$classification
    # jointemb_lapARI<- adjustedRandIndex(jointembpreds,sudo_classes)

    #-----------lap_separate ASE, procrustes, average ----
    X1 <- ase2(graphs_laps[[1]],d=m)$X
    X2 <- ase2(graphs_laps[[2]],d=m)$X
    "Procrustes"
    W<-procrustes(X1,X2)
    X1_new<-X1%*%W$W
    # get average
    AvgProc <- (X1_new+X2)/2
    preds <- Mclust(AvgProc)
    Procpreds <- preds$classification
    Proc_lapARI<- adjustedRandIndex(Procpreds,sudo_classes)

    #-----------lap_average first, then embed----
    avgadj <- (graphs_laps[[1]]+graphs_laps[[2]])/2
    ASEavg <- ase2(avgadj,d=m)$X
    preds <- Mclust(ASEavg)
    Avgpreds<- preds$classification
    Avg_lapARI <- adjustedRandIndex(Avgpreds,sudo_classes)
    #

    return(c(OmniARI, MaseARI, MRDPGARI,jointembARI,ProcARI,AvgARI,
             Omni_lapARI,Mase_lapARI,Proc_lapARI,Avg_lapARI))
    
  })
  df[[j]] = data.frame(t(ARIs[[j]])) %>% 
    rename(OmniARI=X1,MaseARI=X2,MRDPGARI=X3,jointembARI=X4,ProcARI=X5,AvgARI=X6,
           Omni_lapARI=X7,Mase_lapARI=X8,Proc_lapARI=X9,Avg_lapARI=X10) %>% 
    tidyr::gather(variable,value) %>% 
    mutate(epsilons=epsilons[j])
}
  
# plot the result----
finaldf <- df %>% bind_rows() %>% 
  group_by(epsilons,variable) %>% 
  summarise_at(vars(value), list(Min = min, Mean = mean, Max = max, Sd = sd)) %>% 
  mutate(se = Sd/sqrt(S)) %>% 
  mutate(AdjLap = ifelse(variable %in% c("Avg_lapARI","Mase_lapARI",
                                         "Omni_lapARI","Proc_lapARI"),  "Lap", "Adj"))



finaldf %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = variable)) +
  #geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-Sd, 
                    #ymax = Mean+Sd,color=variable,linetype=variable),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() + 
  scale_x_continuous(name = "Epsilons", breaks = epsilons) +
  ylab( "Mean of simulation of ARIs") +
  ggtitle("Four blocks SBM with same P") +
  scale_color_discrete(name = "Variable") +
  scale_fill_discrete(name = "Variable") +
  scale_linetype_discrete(name = "Variable")

# Laplacian matrix----
ARIs2 <- list()
df2 = list()
epsilons = seq(0.01,0.1, by=0.008)
n = 100 # for each block
m = 2 # number of graphs
#----change the value of epsilon
for (j in 1:length(epsilons)) {
  print(epsilons[j])
  P = matrix(rep(0.1,4), ncol = 2)
  qr = matrix(c(0.9,0.01,0.01,0.3),ncol = 2)
  B = P+epsilons[j]*qr
  
  S = 100
  seeds = 1000:(1000-1+S)
  ARIs2[[j]] <- sapply(1:S, function(i) {
    set.seed(seeds[i])
    print(i)
    graphs_adjs2 = list()
    for (k in 1:m) {
      graphs_adjs2[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=c(n/2,n/2))
      graphs_adjs2[[k]] = laplacian_matrix(graphs_adjs2[[k]], sparse = F)
    }
    
    #-----------Omni_lap----
    OmniM<-buildOmni(graphs_laps)
    Xomni <- ase2(OmniM,m)$X
    # get average and prediction
    AvgOmni <- (Xomni[1:n,]+Xomni[(n+1):(2*n),])/m
    preds <- Mclust(AvgOmni)
    Omnipreds <- preds$classification
    OmniARI<- adjustedRandIndex(Omnipreds,sudo_classes)
    
    #------------MASE_lap----
    scaled_MASE <- mase(graphs_laps, d=dims)
    V <- scaled_MASE$V
    R1 <- scaled_MASE$R[[1]]
    R2 <- scaled_MASE$R[[2]]
    'Eigen-decomp of R_1, R_2'
    sig1 <- cbind( c( eigen(R1)$values[1],0), c(0, eigen(R1)$values[2]) )
    sig2 <- cbind( c( eigen(R2)$values[1],0), c(0, eigen(R2)$values[2]) )
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
    MaseARI<- adjustedRandIndex(Masepreds,sudo_classes)
    
    # #-----------MRDPG----
    # scaled_MRDPG <- multiRDPG(graphs_adjs2, d=m)
    # U <- scaled_MRDPG$U
    # L1 <- scaled_MRDPG$Lambda[[1]]
    # L2 <- scaled_MRDPG$Lambda[[2]]
    # 'Estimated latent positions'
    # X1 <- U%*%L1^{1/2}
    # X2 <- U%*%L2^{1/2}
    # # get average
    # AvgMRDPG <- (X1+X2)/2
    # preds <- Mclust(AvgMRDPG)
    # MRDPGpreds <- preds$classification
    # MRDPGARI<- adjustedRandIndex(MRDPGpreds,sudo_classes)
    # 
    #-----------joint embedding----
    jointemb <- multidembed(graphs_laps,d=m)
    jointemb <- jointemb$h
    preds <- Mclust(jointemb)
    jointembpreds <- preds$classification
    jointembARI<- adjustedRandIndex(jointembpreds,sudo_classes)
    
    
    #-----------separate ASE, procrustes, average---- 
    X1 <- ase(graphs_laps[[1]],d=m)
    X2 <- ase(graphs_laps[[2]],d=m)
    "Procrustes"
    W<-procrustes(X1,X2)
    X1_new<-X1%*%W$W
    # get average
    AvgProc <- (X1_new+X2)/2
    preds <- Mclust(AvgProc)
    Procpreds <- preds$classification
    ProcARI<- adjustedRandIndex(Procpreds,sudo_classes)
    
    #-----------average first, then embed----
    avgadj <- (graphs_adjs2[[1]]+graphs_adjs2[[2]])/2
    ASEavg <- ase(avgadj,d=m)
    preds <- Mclust(ASEavg)
    Avgpreds<- preds$classification
    AvgARI <- adjustedRandIndex(Avgpreds,sudo_classes)
    
    return(c(OmniARI, MaseARI,jointembARI,ProcARI,AvgARI))
    
  })
  df2[[j]] = data.frame(t(ARIs[[j]])) %>% 
    rename(OmniARI=X1,MaseARI=X2,MRDPGARI=X3,jointembARI=X4,ProcARI=X5,AvgARI=X6) %>% 
    tidyr::gather(variable,value) %>% 
    mutate(epsilons=epsilons[j])
}

# plot the result
finaldf2 <- df2 %>% bind_rows() %>% 
  group_by(epsilons,variable) %>% 
  summarise_at(vars(value), list(Min = min, Mean = mean, Max = max, Sd = sd))

finaldf2 %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = variable)) +
  geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-Sd, 
                    ymax = Mean+Sd,color=variable,linetype=variable),width = 0.03)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() + 
  scale_x_continuous(name = "Epsilons", breaks = epsilons) +
  ylab( "Mean of simulation of ARIs") +
  ggtitle("Two block SBM with same P (Laplacian matrix)") +
  scale_color_discrete(name = "Variable") +
  scale_fill_discrete(name = "Variable") +
  scale_linetype_discrete(name = "Variable")








# change m ----

ARIs <- list()
df2 = list()
#epsilons = seq(0.005,0.2, by=0.02) #old
#epsilons = seq(0.05,0.8, by=0.08)
epsilons = seq(0.005,0.9, by=0.14) #new
n = 120 # for total number of vertex
m = 25 # number of graphs
blocks = 4 # number of blocks
dims = blocks # dimension
#sudo_classes <- c(rep(1,n/2),rep(2,n/2))
sudo_classes <- c(rep(1,n/4),rep(2,n/4),rep(3,n/4),rep(4,n/4))
#change the value of epsilon----
for (j in 1:length(epsilons)) {
  print(epsilons[j])
  #P = matrix(rep(0.1,4), ncol = 2)
  #P = matrix(rep(0.1,9), ncol = 3)
  #P = matrix(rep(0.1,16), ncol = 4)
  #qr = matrix(c(0.9,0.01,0.01,0.3),ncol = 2)
  #qr = matrix(c(0.9,0.01,0.01,0.01,0.5,0.01,0.01,0.01,0.3),ncol = 3)
  # qr = matrix(rep(0.01,16),ncol = 4)
  # diag(qr)<- c(0.9,0.6,0.4,0.2)
  # B = P+epsilons[j]*qr
  
  ##---- different graphs probs ----
  P = matrix(rep(0.1,16), ncol = 4)
  qr1 = matrix(rep(0.01,16),ncol = 4)
  diag(qr1)<- c(0.9,0.6,0.4,0.2)
  qr2 = matrix(rep(0.04,16),ncol = 4)
  diag(qr2)<- c(0.4,0.3,0.8,0.9)
  B1 = P+epsilons[j]*qr1
  B2 = P+epsilons[j]*qr2

  S = 100
  seeds = 1000:(1000-1+S)
  ARIs[[j]] <- sapply(1:S, function(i) {
    set.seed(seeds[i])
    graphs_adjs = list()
    graphs_laps = list()
    for (k in 1:12) {
      #graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=c(n/2,n/2))
      #graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=c(n/3,n/3,n/3))
      graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B1, block.sizes=c(n/4,n/4,n/4,n/4))
      graphs_adjs[[k]] = as_adjacency_matrix(graphs_adjs[[k]], sparse = F)
      graphs_laps[[k]] = makeLaplancian(graphs_adjs[[k]])
    }
    for (k in 13:m) {
      #graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=c(n/2,n/2))
      #graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=c(n/3,n/3,n/3))
      graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B2, block.sizes=c(n/4,n/4,n/4,n/4))
      graphs_adjs[[k]] = as_adjacency_matrix(graphs_adjs[[k]], sparse = F)
      graphs_laps[[k]] = makeLaplancian(graphs_adjs[[k]])
    }
    
    #-----------Omni----
    OmniM<-buildOmni(graphs_adjs)
    Xomni <- ase2(OmniM,dims)$X
    # get average and prediction
    array_split <- function(data, m) {
      rowIdx <- seq_len(nrow(data))    
      lapply(split(rowIdx,cut(rowIdx, m)), function(x) data[x, ])
    }
    OmniXlist <-array_split(data=Xomni, m = m)
    AvgOmni <-Reduce("+", OmniXlist) / m
    preds <- Mclust(AvgOmni) 
    Omnipreds <- preds$classification
    OmniARI<- adjustedRandIndex(Omnipreds,sudo_classes)
    
    #------------MASE----
    scaled_MASE <- mase(graphs_adjs, d=dims)
    V <- scaled_MASE$V
    # R list
    R_list <- scaled_MASE$R
    # sig list 'Eigen-decomp of R_1, R_2'
    sigs_list <- list()
    for (k in 1:m) {
      sigs_list[[k]] <-  matrix(0, nrow = dims, ncol =dims)
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
    preds <- Mclust(AvgMase)
    Masepreds <- preds$classification
    MaseARI<- adjustedRandIndex(Masepreds,sudo_classes)
    
    #-----------MRDPG----
    scaled_MRDPG <- multiRDPG(graphs_adjs, d=dims)
    U <- scaled_MRDPG$U
    'Estimated latent positions'
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]] <- scaled_MRDPG$Lambda[[k]]
      Xs[[k]] <- U%*%Xs[[k]]^{1/2}
    }
    # get average
    AvgMRDPG <- Reduce("+", Xs) / m
    preds <- Mclust(AvgMRDPG)
    MRDPGpreds <- preds$classification
    MRDPGARI<- adjustedRandIndex(MRDPGpreds,sudo_classes)
    
    #-----------joint embedding----
    # jointemb <- multidembed(graphs_adjs,d=dims)
    # jointemb <- jointemb$h
    # preds <- Mclust(jointemb)
    # jointembpreds <- preds$classification
    # jointembARI<- adjustedRandIndex(jointembpreds,sudo_classes)
    jointemb <- multidembed(graphs_adjs,d=dims)
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]]<- matrix(0,dims,dims)
      diag(Xs[[k]])<- abs(jointemb$lambda[k,])^{1/2}
      Xs[[k]] <- jointemb$h%*%Xs[[k]]^{1/2}
    }
    
    # get average
    Avgjoint <- Reduce("+", Xs) / m
    preds <- Mclust(Avgjoint)
    jointpreds <- preds$classification
    jointembARI<- adjustedRandIndex(jointpreds,sudo_classes)
    
    #-----------separate ASE, procrustes, average ----
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]] <- ase2(graphs_adjs[[k]],d=dims)$X
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
    preds <- Mclust(AvgProc)
    Procpreds <- preds$classification
    ProcARI<- adjustedRandIndex(Procpreds,sudo_classes)
    
    # #-----------average first, then embed----
    # avgadj <- Reduce("+", graphs_adjs) / m
    # ASEavg <- ase2(avgadj,d=dims)$X
    # preds <- Mclust(ASEavg)
    # Avgpreds<- preds$classification
    # AvgARI <- adjustedRandIndex(Avgpreds,sudo_classes)
    
    
    #-----------Omni_lap----
    OmniM<-buildOmni(graphs_laps)
    Xomni <- ase2(OmniM,m)$X
    OmniXlist <-array_split(data=Xomni, m = m)
    AvgOmni <-Reduce("+", OmniXlist) / m
    preds <- Mclust(AvgOmni)
    Omnipreds <- preds$classification
    Omni_lapARI<- adjustedRandIndex(Omnipreds,sudo_classes)
    
    #------------MASE_lap----
    scaled_MASE <- mase(graphs_laps, d=dims)
    V <- scaled_MASE$V
    # R list
    R_list <- scaled_MASE$R
    # sig list 'Eigen-decomp of R_1, R_2'
    sigs_list <- list()
    for (k in 1:m) {
      sigs_list[[k]] <-  matrix(0, nrow = dims, ncol =dims)
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
    preds <- Mclust(AvgMase)
    Masepreds <- preds$classification
    Mase_lapARI<- adjustedRandIndex(Masepreds,sudo_classes)
    
    
    # #-----------jointembedding_lap----
    # jointemb <- multidembed(graphs_laps,d=dims)
    # jointemb <- jointemb$h
    # preds <- Mclust(jointemb)
    # jointembpreds <- preds$classification
    # jointemb_lapARI<- adjustedRandIndex(jointembpreds,sudo_classes)
    
    #-----------lap_separate ASE, procrustes, average ----
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]] <- ase2(graphs_laps[[k]],d=dims)$X
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
    preds <- Mclust(AvgProc)
    Procpreds <- preds$classification
    Proc_lapARI<- adjustedRandIndex(Procpreds,sudo_classes)
    
    # #-----------lap_average first, then embed----
    # avgadj <- Reduce("+", graphs_laps) / m
    # ASEavg <- ase2(avgadj,d=dims)$X
    # preds <- Mclust(ASEavg)
    # Avgpreds<- preds$classification
    # Avg_lapARI <- adjustedRandIndex(Avgpreds,sudo_classes)
    # #
    
    return(c(OmniARI, MaseARI, MRDPGARI,jointembARI,ProcARI,
             Omni_lapARI,Mase_lapARI,Proc_lapARI))
    
  })
  df2[[j]] = data.frame(t(ARIs[[j]])) %>% 
    rename(OmniARI=X1,MaseARI=X2,MRDPGARI=X3,jointembARI=X4,ProcARI=X5,
           Omni_lapARI=X6,Mase_lapARI=X7,Proc_lapARI=X8) %>% 
    tidyr::gather(variable,value) %>% 
    mutate(epsilons=epsilons[j])
}

# plot the result----
finaldf2 <- df2 %>% bind_rows() %>% 
  group_by(epsilons,variable) %>% 
  summarise_at(vars(value), list(Min = min, Mean = mean, Max = max, Sd = sd)) %>%
  mutate(se = Sd/sqrt(S)) %>% 
  mutate(AdjLap = ifelse(variable %in% c("Avg_lapARI","Mase_lapARI",
                                         "Omni_lapARI","Proc_lapARI"),  "Lap", "Adj"))

finaldf2 %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = AdjLap)) +
  geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-se, 
                    ymax = Mean+se,color=variable,linetype=AdjLap),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() + 
  scale_x_continuous(name = "Epsilons", breaks = epsilons) +
  ylab( "Mean of simulation of ARIs") +
  ggtitle("Four blocks SBM with different P (12/13), n=120, m=25") +
  scale_color_discrete(name = "Variable") +
  scale_fill_discrete(name = "Variable") +
  scale_linetype_discrete(name = "Variable")







#-----------------SBM_corr figure
  
  
m2 = readRDS("./codes/Algcompare/SBM_corr/vertex120.rds")
m4 = readRDS("./codes/Algcompare/SBM_corr/m4.rds")
m6 = readRDS("./codes/Algcompare/SBM_corr/m6.rds")
ver220 = readRDS("./codes/Algcompare/SBM_corr/vertex220.rds")
ver320 = readRDS("./codes/Algcompare/SBM_corr/vertex320.rds")
ver420 = readRDS("./codes/Algcompare/SBM_corr/vertex420.rds")

p1 = m2 %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = AdjLap)) +
  geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-se, 
              ymax = Mean+se,color=variable,linetype=AdjLap),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "", breaks = seq(min(m2$epsilons), max(m2$epsilons), by = 0.04)) +
  ylab( "Mean of simulation of ARIs") +
  #ylab("")+
  xlab("")+
  #ggtitle("Two blocks SBM with same P, n=120, m=2") +
  ggtitle("n = 120, m = 2")+ 
  scale_color_discrete(name = "Methods") +
  scale_fill_discrete(name = "Methods") +
  scale_linetype_discrete(name = "Adj/Lap")+
  theme(plot.margin=unit(c(-0.05,0.05,-0.05,0.05), "cm"))
p1
p2 = m4 %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = AdjLap)) +
  geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-se, 
                    ymax = Mean+se,color=variable,linetype=AdjLap),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "", breaks = seq(min(m2$epsilons), max(m2$epsilons), by = 0.04)) +
  #ylab( "Mean of simulation of ARIs") +
  ylab("")+
  xlab("")+
  #ggtitle("Two blocks SBM with same P, n=120, m=4") +
  ggtitle("n = 120, m = 4")+ 
  scale_color_discrete(name = "Variable") +
  scale_fill_discrete(name = "Variable") +
  scale_linetype_discrete(name = "Variable")+
  theme(plot.margin=unit(c(-0.05,0.05,-0.05,0.05), "cm"))
p2
p3 = m6 %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = AdjLap)) +
  geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-se, 
                    ymax = Mean+se,color=variable,linetype=AdjLap),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "", breaks = seq(min(m2$epsilons), max(m2$epsilons), by = 0.04)) +
  #ylab( "Mean of simulation of ARIs") +
  ylab("")+
  xlab("")+
  #ggtitle("Two blocks SBM with same P, n=120, m=6") +
  ggtitle("n = 120, m = 6")+ 
  scale_color_discrete(name = "Variable") +
  scale_fill_discrete(name = "Variable") +
  scale_linetype_discrete(name = "Variable")+
  theme(plot.margin=unit(c(-0.05,0.05,-0.05,0.05), "cm"))
p3
p4 = ver220 %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = AdjLap)) +
  geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-se, 
                    ymax = Mean+se,color=variable,linetype=AdjLap),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "", breaks = seq(min(m2$epsilons), max(m2$epsilons), by = 0.04)) +
  ylab( "Mean of simulation of ARIs") +
  #ylab("")+
  xlab("")+
  #ggtitle("Two blocks SBM with same P, n=220, m=2") +
  ggtitle("n = 220, m = 2")+ 
  scale_color_discrete(name = "Variable") +
  scale_fill_discrete(name = "Variable") +
  scale_linetype_discrete(name = "Variable")+
  theme(plot.margin=unit(c(-0.05,0.05,-0.05,0.05), "cm"))
p4
p5 = ver320 %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = AdjLap)) +
  geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-se, 
                    ymax = Mean+se,color=variable,linetype=AdjLap),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Epsilons", breaks = seq(min(m2$epsilons), max(m2$epsilons), by = 0.04)) +
  #ylab( "Mean of simulation of ARIs") +
  ylab("")+
  xlab("")+
  #ggtitle("Two blocks SBM with same P, n=320, m=2") +
  ggtitle("n = 320, m = 2")+ 
  scale_color_discrete(name = "Variable") +
  scale_fill_discrete(name = "Variable") +
  scale_linetype_discrete(name = "Variable")+
  theme(plot.margin=unit(c(-0.05,0.05,-0.05,0.05), "cm"))
p5
p6 = ver420 %>% ggplot() +
  geom_line(aes(x = epsilons, y = Mean, color = variable, linetype = AdjLap)) +
  geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-se, 
                    ymax = Mean+se,color=variable,linetype=AdjLap),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "", breaks = seq(min(m2$epsilons), max(m2$epsilons), by = 0.04)) +
  #ylab( "Mean of simulation of ARIs") +
  ylab("")+
  xlab("")+
  #ggtitle("Two blocks SBM with same P, n=420, m=2") +
  ggtitle("n = 420, m = 2")+ 
  scale_color_discrete(name = "Variable") +
  scale_fill_discrete(name = "Variable") +
  scale_linetype_discrete(name = "Variable")+ 
  theme(plot.margin=unit(c(-0.05,0.05,-0.05,0.05), "cm")) 
p6

library(ggpubr)
cluster_samep = ggarrange(
  p1, p2, p3,p4,p5,p6,
  ncol = 3, nrow = 2, #labels = c("A","B"),
  common.legend = TRUE, legend = "right"
) 

annotate_figure(cluster_samep, 
                top = text_grob(" ", 
                               color = "black",face = "bold", size = 3),
                bottom = text_grob(" ", 
                                color = "black",face = "bold", size = 3)
                #left = text_grob("Mean of simulation of ARIs", color = "black", rot = 90),
                #bottom = text_grob("Epsilons")
                )













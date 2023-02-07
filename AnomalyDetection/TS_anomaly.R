source("./functions/functions.R")
source("./functions/mase.R")
source("./functions/parametric-bootstrap.R")
source("./functions/sample_graphs.R")
source("./functions/jointembedding.R")
library(igraph)
library(iGraphMatch)
library(irlba)
library(multiRDPG)
library(mclust) 
library(dplyr)
library(ggplot2)

# Laplancian function
makeLaplancian <- function(A){
  r<- rowSums(A)
  smartroot<-function(y){if(y>0){y=y^(-1/2); return(y)}else {return(y)}}
  rr<-sapply(r,smartroot)
  D <- diag(rr)
  L <- D%*%A%*%D
}
# split matrix
array_split <- function(data, m) {
  rowIdx <- seq_len(nrow(data))    
  lapply(split(rowIdx,cut(rowIdx, m)), function(x) data[x, ])
}

# BBs----
n = 200
m = 20 # number of graphs
j = 10 # the jth one is anomalous
d = 10 # dim to embed
epsilons = seq(0.05,0.89, by=0.1)

dflist <- list()
for (e in 1:length(epsilons)) {
  B = matrix(rep(0.1,100), ncol = 10)
  B_j <- B
  # blocks ----
  diag(B_j)[1:1]<- 0.1+epsilons[e]
  
  print(epsilons[e])
  
  Ds <- list()
  S = 100
  seeds = 1000:(1000-1+S)
  Ds <- lapply(1:S, function(i) {  
    set.seed(seeds[i])
    graphs_adjs = list()
    graphs_laps = list()
    for (k in 1:m) {
      graphs_adjs[[k]] = sample_sbm(n,pref.matrix=B, block.sizes=rep(n/d,d))
      graphs_adjs[[k]] = as_adjacency_matrix(graphs_adjs[[k]], sparse = F)
      graphs_laps[[k]] = makeLaplancian(graphs_adjs[[k]])
    }
    
    # generate anomalous graph ----
    graphs_adjs[[j]] = sample_sbm(n,pref.matrix=B_j, block.sizes=rep(n/d,d))
    #graphs_adjs[[j]] = sample_gnp(n, p=0.5)
    graphs_adjs[[j]] = as_adjacency_matrix(graphs_adjs[[j]], sparse = F)
    graphs_laps[[j]] = makeLaplancian(graphs_adjs[[j]])
    # compute ds----
    #-----------Omni----
    OmniM<-buildOmni(graphs_adjs)
    Xomni <- ase2(OmniM,d)$X
    OmniXlist <-array_split(data=Xomni, m = m)
    
    Dists_Omni = c()
    for (aa in 2:m) {
      Dists_Omni[aa-1] <- norm(OmniXlist[[aa]]-OmniXlist[[aa-1]], "F")
    }
    
    # Omni_Lap ----
    OmniM<-buildOmni(graphs_laps)
    Xomni <- ase2(OmniM,d)$X
    OmniXlist <-array_split(data=Xomni, m = m)
    
    Dists_OmniLap = c()
    for (aa in 2:m) {
      Dists_OmniLap[aa-1] <- norm(OmniXlist[[aa]]-OmniXlist[[aa-1]], "F")
    }
    
    #-------MASE-------
    scaled_MASE <- mase(graphs_adjs, d=d)
    V <- scaled_MASE$V
    # R list
    R_list <- scaled_MASE$R
    # sig list 'Eigen-decomp of R_1, R_2'
    sigs_list <- list()
    for (k in 1:m) {
      sigs_list[[k]] <-  matrix(0, nrow = d, ncol =d)
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
    
    Dists_MASE = c()
    for (aa in 2:m) {
      Dists_MASE[aa-1] <- norm(Xs[[aa]]-Xs[[aa-1]], "F")
    }
    
    #-----------MRDPG----
    scaled_MRDPG <- multiRDPG(graphs_adjs, d=d)
    U <- scaled_MRDPG$U
    'Estimated latent positions'
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]] <- scaled_MRDPG$Lambda[[k]]
      Xs[[k]] <- U%*%Xs[[k]]^{1/2}
    }
    Dists_MRDPG = c()
    for (aa in 2:m) {
      Dists_MRDPG[aa-1] <- norm(Xs[[aa]]-Xs[[aa-1]], "F")
    }
    
    #-----------joint embedding----
    jointemb <- multidembed(graphs_adjs,d=d)
    Xs <- list()
    for (k in 1:m) {
      Xs[[k]]<- matrix(0,d,d)
      diag(Xs[[k]])<- abs(jointemb$lambda[k,])^{1/2}
      Xs[[k]] <- jointemb$h%*%Xs[[k]]^{1/2}
    }
    
    Dists_jointemb = c()
    for (aa in 2:m) {
      Dists_jointemb[aa-1] <- norm(Xs[[aa]]-Xs[[aa-1]], "F")
    }
    
    # sequential procrustes ----
    sep_ases<- list()
    for (k in 1:m) {
      sep_ases[[k]] <- ase2(graphs_adjs[[k]],d=d)$X
    }
    new_Xs <- list()
    for (k in 2:m) {
      W <- procrustes(sep_ases[[k]], sep_ases[[k-1]])
      new_Xs[[k]] <- sep_ases[[k]]%*%W$W
    }
    new_Xs[[1]]<- sep_ases[[1]]
    Dists_sqProc =c()
    for (aa in 2:m) {
      Dists_sqProc[aa-1] <- norm(new_Xs[[aa]]-sep_ases[[aa-1]], "F")
    }
    
    # seqPro Lap ----
    sep_ases<- list()
    for (k in 1:m) {
      sep_ases[[k]] <- ase2(graphs_laps[[k]],d=d)$X
    }
    new_Xs <- list()
    for (k in 2:m) {
      W <- procrustes(sep_ases[[k]], sep_ases[[k-1]])
      new_Xs[[k]] <- sep_ases[[k]]%*%W$W
    }
    new_Xs[[1]]<- sep_ases[[1]]
    Dists_sqProcLap =c()
    for (aa in 2:m) {
      Dists_sqProcLap[aa-1] <- norm(new_Xs[[aa]]-sep_ases[[aa-1]], "F")
    }
    
    
    return(data.frame(Dists_Omni,Dists_MASE,Dists_MRDPG,Dists_jointemb,Dists_sqProc,
                      Dists_OmniLap,Dists_sqProcLap))
    
  })
  
  # average the Ds ----
  AvgDs <- Reduce("+", Ds) / S
  AvgDs$m <- as.numeric(c(1:19))
  dflist[[e]] = AvgDs %>% 
    rename(jointemb=Dists_jointemb, MASE=Dists_MASE,MRDPG=Dists_MRDPG, 
           Omni=Dists_Omni,SeqProc=Dists_sqProc,
           SeqProcLap=Dists_sqProcLap, OmniLap =Dists_OmniLap) %>%
    tidyr::gather("methods", "value", 1:7) %>% 
    mutate(epsilons=epsilons[e])
}


# plot ----
bigdf <- dflist %>%  bind_rows()
bigdf %>% ggplot(aes(x = m, y=value, color = methods,linetype = methods)) +
  geom_line()+
  facet_wrap(~ epsilons,  scales = "free_y", ncol = 3) +
  #geom_errorbar(aes(x = epsilons, y = Mean, ymin = Mean-Sd, 
  #                  ymax = Mean+Sd,color=variable,linetype=variable),width=0.005)+
  #geom_ribbon(aes(x = epsilons, y = Mean,ymax=Max, ymin=Min,fill = variable),alpha=0.2) +
  theme_bw() +
  #scale_x_continuous(name = "Epsilons", breaks = epsilons) +
  ylab( "Distances between consecutive latent positions of the graphs") +
  ggtitle("SBMs with different P on 1 block, n=200, d=10") 











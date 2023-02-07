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
library(network)
library(HelpersMG)

library(AnomalyDetection)

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


# ENRON data ----

enron <- read.csv("https://www.cis.jhu.edu/~parky/Enron/execs.email.linesnum",
                  sep = "",header = F )
colnames(enron) <- c("time","from","to")

library("lubridate")

enron$time <- as_datetime(enron$time)
enron$time <- as.Date(enron$time)
enron$time <- format(as.Date(enron$time), "%Y-%m")
table(enron$time)

#ggplot(enron, aes(x=time)) + geom_histogram(color="black", fill="white")
enron$ind <- 1:nrow(enron)
split_enron <- split(enron, format(as.Date(enron$time), "%Y-%m")) # split by month/week

# get subset for specific time
sub_2002 <- split_enron[25:44]
sub_2002 <- data.frame(bind_rows(sub_2002))
sub_2002$time <- as.Date(sub_2002$time)

sub_2002_split <- split(sub_2002,format(as.Date(sub_2002$time), "%Y-%m"))


# get the adjs
middles <- list()
enron_adjs <- list()
dims<-data.frame()
for (i in 1:length(sub_2002_split)) {
  sub_2002_split[[i]] <- sub_2002_split[[i]][,c(2,3)]
  middles[[i]] <- as.matrix(get.adjacency(graph.data.frame(sub_2002_split[[i]])))
  enron_adjs[[i]]<- ifelse(middles[[i]]>1,1,middles[[i]])
  enron_adjs[[i]] <- symmetricize(enron_adjs[[i]],"max")
  print(dim(enron_adjs[[i]]))
  dims <- rbind(dims, dim(enron_adjs[[i]]))
}

colnames_list <- list()
colnames_list <- lapply(1:length(enron_adjs), function(i){
  colnames(enron_adjs[[i]])
})

# select range of adjs
commonnumber <- c(1:18)
length(commonnumber)
#length(Reduce(intersect, colnames_list[commonnumber]))
commonvertex <- Reduce(intersect, colnames_list[commonnumber])
length(commonvertex)

enron_commonadjs <- list()
for (i in 1:length(commonnumber)) {
  enron_commonadjs[[i]] <- enron_adjs[[commonnumber[i]]][commonvertex,commonvertex]
}


m = length(enron_commonadjs) # number of graphs
d = 2 # number of dim to embed

enronlaps <- list()
for (k in 1:m) {
  enronlaps[[k]] = makeLaplancian(enron_commonadjs[[k]])
}


#-----------Omni----
OmniM<-buildOmni(enron_commonadjs)
Xomni <- ase2(OmniM,d)$X
OmniXlist <-array_split(data=Xomni, m = m)

Dists_Omni = c()
for (aa in 2:m) {
  Dists_Omni[aa-1] <- norm(OmniXlist[[aa]]-OmniXlist[[aa-1]], "F")
}

# Omni_Lap ----
OmniM<-buildOmni(enronlaps)
Xomni <- ase2(OmniM,d)$X
OmniXlist <-array_split(data=Xomni, m = m)

Dists_OmniLap = c()
for (aa in 2:m) {
  Dists_OmniLap[aa-1] <- norm(OmniXlist[[aa]]-OmniXlist[[aa-1]], "F")
}

#-------MASE-------
scaled_MASE <- mase(enron_commonadjs, d=d)
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
scaled_MRDPG <- multiRDPG(enron_commonadjs, d=d)
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
jointemb <- multidembed(enron_commonadjs,d=d)
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
  sep_ases[[k]] <- ase2(enron_commonadjs[[k]],d=d)$X
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
  sep_ases[[k]] <- ase2(enronlaps[[k]],d=d)$X
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


Dists <- data.frame(Dists_Omni,Dists_MASE,Dists_MRDPG,Dists_jointemb,Dists_sqProc,
                    Dists_OmniLap,Dists_sqProcLap)
Dists$m <- as.numeric(c(1:(length(enron_commonadjs)-1)))
Distsdf <- Dists %>% 
  rename(jointemb=Dists_jointemb, MASE=Dists_MASE,MRDPG=Dists_MRDPG, 
         Omni=Dists_Omni,SeqProc=Dists_sqProc,
         SeqProcLap=Dists_sqProcLap, OmniLap =Dists_OmniLap) %>% 
  tidyr::gather("methods", "value", 1:7)

Distsdf %>% ggplot(aes(x = m, y=value, color = methods,linetype = methods)) +
  geom_line()+
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 9),
        text = element_text(size = 9))+
  scale_x_continuous(label = scales::label_comma(accuracy = 1),breaks = Dists$m)+
  ylab( "Distances between consecutive latent positions") +
  ggtitle("Enron Email data") 








# Zebrafish Brain data ----
load("./codes/Algcompare/TS/Data_and_Code/tsg-11-30-469.Rbin")

set.seed(12345)

# Zebrafish Brain Data 
Zeb<-list()
for(i in 1:20){
  A<-as.matrix(tsg21[[i]])
  g<-graph.adjacency(A,mode=c("undirected"))
  #D[[i]]<-1-similarity(g,method=c("jaccard"))
  Zeb[[i]] <- as.matrix(g[])
}

m = length(Zeb) # number of graphs
d = 2 # number of dim to embed

Zeblaps <- list()
for (k in 1:m) {
  Zeblaps[[k]] = makeLaplancian(Zeb[[k]])
}

    #-----------Omni----
    OmniM<-buildOmni(Zeb)
    Xomni <- ase2(OmniM,d)$X
    OmniXlist <-array_split(data=Xomni, m = m)
    
    Dists_Omni = c()
    for (aa in 2:m) {
      Dists_Omni[aa-1] <- norm(OmniXlist[[aa]]-OmniXlist[[aa-1]], "F")
    }
    
    # Omni_Lap ----
    OmniM<-buildOmni(Zeblaps)
    Xomni <- ase2(OmniM,d)$X
    OmniXlist <-array_split(data=Xomni, m = m)
    
    Dists_OmniLap = c()
    for (aa in 2:m) {
      Dists_OmniLap[aa-1] <- norm(OmniXlist[[aa]]-OmniXlist[[aa-1]], "F")
    }
    
    #-------MASE-------
    scaled_MASE <- mase(Zeb, d=d)
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
    scaled_MRDPG <- multiRDPG(Zeb, d=d)
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
    jointemb <- multidembed(Zeb,d=d)
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
      sep_ases[[k]] <- ase2(Zeb[[k]],d=d)$X
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
      sep_ases[[k]] <- ase2(Zeblaps[[k]],d=d)$X
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
    
    
    Dists <- data.frame(Dists_Omni,Dists_MASE,Dists_MRDPG,Dists_jointemb,Dists_sqProc,
                        Dists_OmniLap,Dists_sqProcLap)
    Dists$m <- as.numeric(c(1:(length(Zeb)-1)))
    Distsdf <- Dists %>% 
      rename(jointemb=Dists_jointemb, MASE=Dists_MASE,MRDPG=Dists_MRDPG, 
             Omni=Dists_Omni,SeqProc=Dists_sqProc,
             SeqProcLap=Dists_sqProcLap, OmniLap =Dists_OmniLap) %>% 
      tidyr::gather("methods", "value", 1:7)
    
    Distsdf %>% ggplot(aes(x = m, y=value, color = methods,linetype = methods)) +
      geom_line()+
      theme_bw() +
      scale_x_continuous(label = scales::label_comma(accuracy = 1),breaks = Dists$m)+
      ylab( "Distances between \nconsecutive latent positions") +
      ggtitle("Zebrafish brain data") 
    



# MSG data ----
msg <- read.delim("./codes/Algcompare/TS/CollegeMsg.txt", header=F, sep=" ")
colnames(msg) <- c("SRC","TGT","timestamp")

library("lubridate")

msg$timestamp <- as_datetime(msg$timestamp)
msg$timestamp <- as.Date(msg$timestamp)

#table(msg$timestamp)
#ggplot(msg, aes(x=timestamp)) + geom_histogram(color="black", fill="white")

split_data <- split(msg, cut(msg$timestamp, 14))

middles <- list()
msg_adjs <- list()
dims<-data.frame()
for (i in 1:length(split_data)) {
  split_data[[i]] <- split_data[[i]][,c(1,2)]
  #middles[[i]] <- dplyr::inner_join(split_data[[i]], split_data[[i]], by = "SRC")#[,-1]
  #middles[[i]] <- apply(middles[[i]], 2, as.character)
  middles[[i]] <- as.matrix(get.adjacency(graph.data.frame(split_data[[i]])))
  msg_adjs[[i]]<- ifelse(middles[[i]]>1,1,middles[[i]])
  msg_adjs[[i]] <- symmetricize(msg_adjs[[i]],"max")

  print(dim(msg_adjs[[i]]))
  dims <- rbind(dims, dim(msg_adjs[[i]]))
}

colnames_list <- list()
colnames_list <- lapply(1:length(split_data), function(i){
  colnames(msg_adjs[[i]])
})

commonnumber <- c(1:11)
#length(commonnumber)
length(Reduce(intersect, colnames_list[commonnumber]))
commonvertex <- Reduce(intersect, colnames_list[commonnumber])
length(commonvertex)

msg_commonadjs <- list()
for (i in 1:length(commonnumber)) {
  msg_commonadjs[[i]] <- msg_adjs[[commonnumber[i]]][commonvertex,commonvertex]
}

Distsdflist = list()

# injection

n=dim(msg_commonadjs[[1]])[1]
p=seq(0.1,0.9,by=0.1)
set.seed(118)
for (e in 1:length(p)) {
  
  msg_commonadjs <- list()
  for (i in 1:length(commonnumber)) {
    msg_commonadjs[[i]] <- msg_adjs[[commonnumber[i]]][commonvertex,commonvertex]
  }
  
  xs = matrix(0, ncol = n, nrow = n)
  xs[lower.tri(xs)] = rbinom(length(xs[lower.tri(xs)]),1,p[e])
  xs=symmetricize(xs,"max")
  # change the 9th graph
  zz=8
  msg_commonadjs[[zz]] = msg_commonadjs[[zz]]*(1-xs)+xs*(1-msg_commonadjs[[zz]])

# # random matrix
# test= matrix(sample(c(0,1), replace=TRUE, prob=c(0.9,0.1), size=n*n), nrow = n, ncol = n)
# test2 = symmetricize(test,"max")
# diag(test2) <- 0
# colnames(test2) = colnames(msg_commonadjs[[1]])
# rownames(test2) = rownames(msg_commonadjs[[1]])
# msg_commonadjs <- append(msg_commonadjs, list(test2),7 )


  m = length(msg_commonadjs) # number of graphs
  d = 2 # number of dim to embed
  
  msglaps <- list()
  for (k in 1:m) {
    msglaps[[k]] = makeLaplancian(msg_commonadjs[[k]])
  }
          
        #-----------Omni----
        OmniM<-buildOmni(msg_commonadjs)
        Xomni <- ase2(OmniM,d)$X
        OmniXlist <-array_split(data=Xomni, m = m)
        
        Dists_Omni = c()
        for (aa in 2:m) {
          Dists_Omni[aa-1] <- norm(OmniXlist[[aa]]-OmniXlist[[aa-1]], "F")
        }
        
        # Omni_Lap ----
        OmniM<-buildOmni(msglaps)
        Xomni <- ase2(OmniM,d)$X
        OmniXlist <-array_split(data=Xomni, m = m)
        
        Dists_OmniLap = c()
        for (aa in 2:m) {
          Dists_OmniLap[aa-1] <- norm(OmniXlist[[aa]]-OmniXlist[[aa-1]], "F")
        }
        
        #-------MASE-------
        scaled_MASE <- mase(msg_commonadjs, d=d)
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
        scaled_MRDPG <- multiRDPG(msg_commonadjs, d=d)
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
        jointemb <- multidembed(msg_commonadjs,d=d)
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
          sep_ases[[k]] <- ase2(msg_commonadjs[[k]],d=d)$X
        }
        new_Xs <- list()
        for (k in 2:m) {
          W <- procrustes(sep_ases[[k]], sep_ases[[k-1]])
          new_Xs[[k]] <- sep_ases[[k]]%*%W$W
        }
        new_Xs[[1]]<- sep_ases[[1]]
        Dists_sqProc =c()
        for (aa in 2:m) {
          #Dists_sqProc[aa-1] <- norm(new_Xs[[aa]]-new_Xs[[aa-1]], "F")
          Dists_sqProc[aa-1] <- norm(new_Xs[[aa]]-sep_ases[[aa-1]], "F")
        }
        
        # seqPro Lap ----
        sep_ases<- list()
        for (k in 1:m) {
          sep_ases[[k]] <- ase2(msglaps[[k]],d=d)$X
        }
        new_Xs <- list()
        for (k in 2:m) {
          W <- procrustes(sep_ases[[k]], sep_ases[[k-1]])
          new_Xs[[k]] <- sep_ases[[k]]%*%W$W
        }
        new_Xs[[1]]<- sep_ases[[1]]
        Dists_sqProcLap =c()
        for (aa in 2:m) {
          #Dists_sqProcLap[aa-1] <- norm(new_Xs[[aa]]-new_Xs[[aa-1]], "F")
          Dists_sqProcLap[aa-1] <- norm(new_Xs[[aa]]-sep_ases[[aa-1]], "F")
        }
        
        
        Dists <- data.frame(Dists_Omni,Dists_MASE,Dists_MRDPG,Dists_jointemb,Dists_sqProc,
                            Dists_OmniLap,Dists_sqProcLap)
        Dists$m <- as.numeric(c(1:(length(msg_commonadjs)-1)))
        Distsdflist[[e]] <- Dists %>% 
          rename(MASE=Dists_MASE,MRDPG=Dists_MRDPG,Jointemb=Dists_jointemb,
                 Omni=Dists_Omni,SeqProc=Dists_sqProc,
                 SeqProcLap=Dists_sqProcLap, OmniLap =Dists_OmniLap) %>% 
          #tidyr::gather("methods", "value", 1:7) %>% 
          mutate(probs=p[e])
        
}
        
bigdf <- Distsdflist %>%  bind_rows() %>% 
  tidyr::gather("methods", "value", 1:7)

bigdf %>% ggplot(aes(x = m, y=value, color = methods,linetype = methods)) +
          geom_line()+
          facet_wrap(~ probs,  scales = "free_y", ncol = 3) +
          theme_bw() +
          scale_x_continuous(label = scales::label_comma(accuracy = 1),breaks = Dists$m)+
          labs( y="Distances between consecutive latent positions",color = "Methods",linetype = "Methods") +
          ggtitle("College Message data") 


# anomaly detection ----
test = split(bigdf, f = bigdf$methods)
res = AnomalyDetectionVec(test[[2]][, 3], max_anoms=0.1, direction='pos', period=10,plot=TRUE)
res = AnomalyDetectionVec(bigdf[, 3], max_anoms=0.1, direction='both', period=10,plot=TRUE)
res$plot
        
        
# BOSTON marathon ----

boston <- read.delim("./codes/Algcompare/TS/BostonBomb2013_activity.txt", header=F, sep=" ")
colnames(boston) <- c("from","to","timestamp","layer")
        
library("lubridate")
        
boston$timestamp <- as_datetime(boston$timestamp)
boston$timestamp <- as.Date(boston$timestamp)
        
table(boston$timestamp)
table(boston$layer)

# subset one layer
boston1 <- subset(boston, boston$layer =="RE")
split_boston <- split(boston1, format(as.Date(boston1$time), "%Y-%m-%d"))

split_boston <- split(boston1, cut(boston1$timestamp,10))

# get the adjs
middles <- list()
boston_adjs <- list()
dims<-data.frame()
for (i in 1:length(split_boston)) {
  split_boston[[i]] <- split_boston[[i]][,c(1,2)]
  middles[[i]] <- as.matrix(get.adjacency(graph.data.frame(split_boston[[i]])))
  boston_adjs[[i]]<- ifelse(middles[[i]]>1,1,middles[[i]])
  boston_adjs[[i]] <- symmetricize(boston_adjs[[i]],"max")
  print(dim(boston_adjs[[i]]))
  dims <- rbind(dims, dim(boston_adjs[[i]]))
}

colnames_list <- list()
colnames_list <- lapply(1:length(enron_adjs), function(i){
  colnames(enron_adjs[[i]])
})

# select range of adjs
commonnumber <- c(1:18)
length(commonnumber)
#length(Reduce(intersect, colnames_list[commonnumber]))
commonvertex <- Reduce(intersect, colnames_list[commonnumber])
length(commonvertex)

enron_commonadjs <- list()
for (i in 1:length(commonnumber)) {
  enron_commonadjs[[i]] <- enron_adjs[[commonnumber[i]]][commonvertex,commonvertex]
}


m = length(enron_commonadjs) # number of graphs
d = 2 # number of dim to embed

enronlaps <- list()
for (k in 1:m) {
  enronlaps[[k]] = makeLaplancian(enron_commonadjs[[k]])
}

        








#' 
#' Function to sample an independent-edge graph given a probability matrix P
#' 
#' @param P matrix encoding the probabilities of an edge in the graph
#' 
#' @return The adjacency matrix of a random graph
#' 
sample_from_P <- function(P) {
  n = ncol(P)
  A = Matrix(0, n, n)
  A[upper.tri(A)] <- 1*(runif(sum(upper.tri(P))) < P[upper.tri(P)])
  A = A + t(A)
  return(A)
}


plot_adjmatrix <- function(edgevalues, type=c("undirected", "directed"),
                           edgetype = c("real", "prob", "binary"),
                           communities = NULL, 
                           community_labels = NULL, 
                           main= "", axislabel = "Nodes",
                           colorlims = NULL) {
  type <- match.arg(type)
  require(lattice)
  require(Matrix)
  edgetype <- match.arg(edgetype)
  if(is.null(dim(edgevalues))) {
    if(type=="undirected") {
      NODES <- (1+sqrt(1+8*length(edgevalues)))/2
    }else{if(type=="directed") {
      NODES <- (1+sqrt(1+8*length(edgevalues)/2))/2
    }else{
      stop("The value of type should be \"undirected\" or \"directed\"")
    }}
    Adj_matr <- as.matrix(get_matrix(edgevalues, type))
  }else{
    Adj_matr <- as.matrix(edgevalues)
    NODES <- ncol(Adj_matr)
  }
  
  tckspace <- round(NODES/5, -floor(log10(NODES/5)))
  
  cuts <- 100
  colorkey <- TRUE
  #edgetype = real----------------------------------------
  atneg <- 0
  atpos <- 0
  col.regions.neg <- rgb(red = 1,green = 1, blue = 1)
  col.regions.pos <- rgb(red = 1,green = 1, blue = 1)
  if(is.null(colorlims)) {
    min_Adj <- min(Adj_matr[!is.na(Adj_matr)])
    max_Adj <- max(Adj_matr[!is.na(Adj_matr)])
  }else{
    min_Adj <- colorlims[1]
    max_Adj <- colorlims[2]
  }
  if(min_Adj <0 ) {
    atneg <- seq(min_Adj ,max(-1e-16,min(Adj_matr)), length.out =cuts)
    col.regions.neg <- rgb((cuts:0)/cuts, green = 0, blue = 1,red = 0)
  }
  if(max_Adj > 0) {
    atpos <- seq(1e-16,max_Adj,length.out =cuts)
    col.regions.pos <- rgb((0:cuts)/cuts, green = 0, blue = 0,red = 1)
  }
  atval = unique(c(atneg, 0, atpos))
  col.vals = unique(c(col.regions.neg, rgb(red = 1,green = 1, blue = 1), col.regions.pos))
  #edgetype = prob----------------------------------------
  if(edgetype == "prob") {
    atval = seq(0,1,length.out =cuts)
    col.vals = rgb((0:cuts)/cuts, green = 0, blue = 1,red = 0)
  }
  #edgetype = binary -------------------------------------
  if(edgetype == "binary") {
    colorkey <- FALSE
    Adj_matr <- 1*(Adj_matr != 0)
    atval <- seq(0,1,length.out =cuts)
    col.vals <- rgb((0:cuts)/cuts, green = 0, blue = 1,red = 0)
  }
  
  # Plot community groupings --------------------------------
  if(!is.null(communities)) {
    lengths_coms = sapply(communities, length)
    scales_list = list(tck = c(0,0),
                       x=list(at=Reduce('+',c(0,lengths_coms[1:(length(lengths_coms)-1)]),accumulate = T ) + 
                                lengths_coms/2, 
                              labels=community_labels),
                       y=list(at=Reduce('+',c(0,rev(lengths_coms)[1:(length(lengths_coms)-1)]),accumulate = T ) + 
                                rev(lengths_coms)/2, 
                              labels=rev(community_labels)))
    panel_func = function(...){ panel.levelplot(...)
      if(type=="prob_cells") {
        select_list <- which(sel_cells,arr.ind = T)
        for(cell in 1:nrow(select_list)) {
          fill_block(select_list[cell,1],select_list[cell,2], communities)
        }
      }
      
      for(u in Reduce('+',lengths_coms,accumulate = T)) {
        panel.abline(v = u+0.5)
      }
      for(u in Reduce('+',rev(lengths_coms),accumulate = T)) {
        panel.abline(h = u+0.5)
      }
    }
    nodeorder = unlist(communities)
    Adj_matr <- Adj_matr[nodeorder, nodeorder]
    Adj_matr <- Adj_matr[,seq(from=ncol(Adj_matr),to=1,by=-1)] #reversing the columns
    levelplot(Adj_matr, at = atval,
              xlab = axislabel, ylab = axislabel,
              main = main,
              colorkey = colorkey,
              col.regions = col.vals,
              panel = panel_func,
              scales=scales_list)
  }else{
    Adj_matr <- Adj_matr[,seq(from=ncol(Adj_matr),to=1,by=-1)] #reversing the columns
    levelplot(Adj_matr, at = atval,
              xlab = axislabel, ylab = axislabel,
              col.regions = col.vals,
              main = main,
              colorkey = colorkey,
              scales = list(tck = c(1,0), 
                            x = list(at=seq(0,ncol(Adj_matr), by = tckspace)),
                            y = list(at = NODES-tckspace- seq(0,ncol(Adj_matr), by = tckspace),
                                     labels = (seq(tckspace,ncol(Adj_matr),tckspace)))))
  }
  
}
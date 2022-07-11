## Mean generality of species in the network
Generality <- function(M){
  return(sum(colSums(M))/sum((colSums(M)!=0)));
}

## Mean vulnerability of the species in the network
Vulnerability <- function(M){
  return(sum(rowSums(M))/sum((rowSums(M)!=0)));
}

## Standard deviation of generality:
SDGenerality <- function(M){
  return(sd(colSums(M) / (sum(M)/dim(M)[1]) ));
}

## Standard deviation of vulnerability:
SDVulnerability <- function(M){
  return(sd(rowSums(M) / (sum(M)/dim(M)[1]) ));
}


## In-degree or number of prey of all species in the network
InDegree <- function(M){
  return(colSums(M));
}

## Out-degree or number of predators of all species in the network
OutDegree <- function(M){
  return(rowSums(M));
}

FractionOfBasal <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));
  
  return(b_sps / dim(M)[1]);
}

## Fraction of top predator species
FractionOfTop <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));
  
  return(t_sps / dim(M)[1]);
}

## Fraction of intermediate consumer species
FractionOfIntermediate <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;
  
  i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));
  
  return(i_sps / dim(M)[1]);
}



## ---------------------------
##
## Script name: space.r
##
## Purpose of script: This script implements functions to generate random spatial layouts of networks of meta-communities
##
## Author: Dr Miguel Lurgi
## Lecturer in Biosciences (Computational Ecology)
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
## and
##
## Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## Date Created: November-2017
##
## Copyright (c) Miguel Lurgi, 2017-2021
## Email: miguel.lurgi@swansea.ac.uk
##
## ---------------------------
##
## Notes:
##
## This script is provided as supplementary material for the paper:
## Galiana, N., Lurgi, M., Claramunt-López, B. et al. The spatial scaling of species interaction networks. 
## Nature Ecology & Evolution 2, 782–790 (2018). https://doi.org/10.1038/s41559-018-0517-3
##
## ---------------------------


## required libraries
require(igraph)

space <- function(N, C)  {
  X = runif(N,0,1)
  Y = runif(N,0,1)
  XY = cbind(X,Y)
  distMat = as.matrix(dist(XY,method = "euclidean", upper = T, diag = T))
  
  distMat[distMat > C] <- 0
  adj_mat <- distMat
  adj_mat[adj_mat!=0] <- 1
  
  
  graf <- graph.adjacency(adj_mat);
  y <- igraph::no.clusters(graf, mode='weak');
  while(y > 1){
    X = runif(N,0,1)
    Y = runif(N,0,1)
    XY = cbind(X,Y)
    distMat = as.matrix(dist(XY,method = "euclidean", upper = T, diag = T))
    
    distMat[distMat > C] <- 0
    adj_mat <- distMat
    adj_mat[adj_mat!=0] <- 1
    
    graf <- graph.adjacency(adj_mat);
    y <- igraph::no.clusters(graf, mode='weak');
    
    print('finding connected landscape')
  }
  
  draw_landscape(XY, adj_mat)
  
  return(list(XY,distMat)) 
}

dispersal <- function(d, BS, beta){
  
  return( exp( (-d/(BS^beta) ) ) )
  
}

draw_landscape <- function(XY, adj_mat){
  quartz(height = 5, width = 6)
  par(mar=c(5,6,2,1))
  
  plot(XY[,1],XY[,2],xlab = "East - West", ylab = "South - North",cex.lab = 1.5, cex.axis = 1.25, main="Meta-community structure") # xy coordinates of each patch
  ConVec = stack(as.data.frame(adj_mat))[,1]
  XX = expand.grid(XY[,1],XY[,1]) #this compares each x coordinate of each patch with all the x coordinates of the others patches. 25*25-1 to avoid putting together the same coordinate with itself
  YY = expand.grid(XY[,2],XY[,2]) #same with y coordinates. es per posar els patches en l'espai en relació al primer
  XX = subset(XX,ConVec==1) #subset of the patches that are connected to each other to draw the arrow
  YY = subset(YY,ConVec==1)
  arrows(x0 = XX[,1],x1=XX[,2],y0 = YY[,1], y1 = YY[,2], length = 0,lwd = 0.1)
  #points(XY[,1],XY[,2],pch=21,bg=col.vec)  
}


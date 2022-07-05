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
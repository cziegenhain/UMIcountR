# backend utils -----------------------------------------------------------

.ham_mat <- function(umistrings) {
  X<- matrix(unlist(strsplit(umistrings, "")),ncol = length(umistrings))
  #function below thanks to Johann de Jong
  #https://goo.gl/u8RBBZ
  uniqs <- unique(as.vector(X))
  U <- X == uniqs[1]
  H <- t(U) %*% U
  for ( uniq in uniqs[-1] ) {
    U <- X == uniq
    H <- H + t(U) %*% U
  }
  nrow(X) - H
}

.hammingFilter<-function(umiseq, edit=1, collapse_mode = c("adjacency","adjacency_directional","adjacency_singleton","cluster")){
  # umiseq a vector of umis, one per read
  uc <- data.table(us = umiseq)[, .N, by = "us"] # normal UMI counts
  setorder(uc, us) #order by sequence
  
  if(length(uc$us)>1 && length(uc$us)<45000){ #prevent use of > 100Gb RAM
    #Sys.time()
    umi <-  .ham_mat(uc$us) #construct pairwise UMI distances
    umi[upper.tri(umi,diag=T)] <- NA #remove upper triangle of the output matrix
    umi <-  data.table(
      row = rep(seq(nrow(umi)), ncol(umi)),
      col = rep(seq(ncol(umi)), each = nrow(umi)),
      value = as.vector(umi)
    )[value <= edit ] #make a long data frame and filter according to cutoff
    umi[, "n.1" := uc[row]$N ][
      , "n.2" := uc[col]$N ] #add in observed freq
    
    
    #adjacency singleton: keep only minor UMI with readcount 1
    if(collapse_mode == "adjacency_singleton"){
      umi <- umi[n.1 <= 1 | n.2 <= 1]
    }
    
    #directional adjacency: keep only minor UMI with readcount 1
    if(collapse_mode == "adjacency_directional"){
      umi <- umi[n.1 <= n.2/2 | n.2 <= n.1/2]
    }
    
    
    umi_out <- copy(umi)
    if(nrow(umi_out) == 0){
      return(umi_out)
    }
    
    if(collapse_mode == "cluster"){
      umi_out[, falseUMI := ifelse( n.1>=n.2, col, row ) ][
        , n.false := ifelse( n.1>=n.2, n.2, n.1 )]
    }else{
      umi_out[, falseUMI := ifelse( n.1>n.2, col, row ) ][
        , n.false := ifelse( n.1>n.2, n.2, n.1 )]
    }
    umi_out[, trueUMI := ifelse( n.1<n.2, col, row ) ][
      , n.true := ifelse( n.1<n.2, n.2, n.1 )][
        , falseUMI := uc[falseUMI]$us ][
          , trueUMI  := uc[trueUMI ]$us][
            , c("row", "col", "value", "n.1", "n.2") := NULL]
    
    if(collapse_mode != "cluster"){
      umi_out <- umi_out[!falseUMI == trueUMI]
    }
    
    dup_daughters <- unique(umi_out[which(duplicated(falseUMI))]$falseUMI)
    if(length(dup_daughters>0)){
      umi_out[,rem := FALSE]
      setorder(umi_out, falseUMI, -n.true)
      setkey(umi_out, falseUMI)
      for(i in dup_daughters){
        umi_out[ i, rem := TRUE ] #remove duplicates
        umi_out[ i, mult = "first" , rem := FALSE] # keep the most frequent parent UMI
      }
      umi_out <- umi_out[rem == FALSE]
      umi_out[, rem := NULL]
    }
    
    
    non_true_UMIs <- unique(umi_out[trueUMI %in% umi_out$falseUMI]$trueUMI)
    real_true_UMIs <- unique(umi_out[!trueUMI %in% umi_out$falseUMI]$trueUMI)
    if(length(non_true_UMIs>0)){
      setkey(umi_out, falseUMI)
      for(i in non_true_UMIs){
        true_parent_UMI <- umi_out[i][!trueUMI %in% non_true_UMIs]$trueUMI
        if(length(true_parent_UMI)==0){#find closest match in case there is no clear parent UMI!
          true_parent_UMI <- real_true_UMIs[stringdist::amatch(umi_out[i][1]$trueUMI, real_true_UMIs, method = "hamming", maxDist=edit, nthread = 1)[1]]
        }
        if(length(true_parent_UMI)>1){ #take a random good parent UMI if more possibilities exist
          true_parent_UMI <- true_parent_UMI[1]
        }
        umi_out[trueUMI == i, trueUMI := true_parent_UMI]
      }
    }
    umi_out[,dist := stringdist::stringdist(falseUMI,trueUMI,method = "hamming", nthread = 1)]
    
    if(collapse_mode != "cluster"){
      umi_out <- umi_out[dist <= edit]
    }
    
    umi_out[, c("n.false","n.true","dist") := NULL]
    return(umi_out)
    
  }else{
    return(NULL)
  }
}
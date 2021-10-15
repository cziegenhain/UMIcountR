# backend utils -----------------------------------------------------------

# .ham_mat <- function(umistrings) {
#   X<- matrix(unlist(strsplit(umistrings, "")),ncol = length(umistrings))
#   #function below thanks to Johann de Jong
#   #https://goo.gl/u8RBBZ
#   uniqs <- unique(as.vector(X))
#   U <- X == uniqs[1]
#   H <- t(U) %*% U
#   for ( uniq in uniqs[-1] ) {
#     U <- X == uniq
#     H <- H + t(U) %*% U
#   }
#   nrow(X) - H
# }

# .hammingFilter<-function(umiseq, edit=1, collapse_mode = c("adjacency","adjacency_directional","adjacency_singleton","cluster")){
#   # umiseq a vector of umis, one per read
#   uc <- data.table(us = umiseq)[, .N, by = "us"] # normal UMI counts
#   setorder(uc, us) #order by sequence
#   
#   if(length(uc$us)>1 && length(uc$us)<45000){ #prevent use of > 100Gb RAM
#     #Sys.time()
#     umi <-  .ham_mat(uc$us) #construct pairwise UMI distances
#     umi[upper.tri(umi,diag=T)] <- NA #remove upper triangle of the output matrix
#     umi <-  data.table(
#       row = rep(seq(nrow(umi)), ncol(umi)),
#       col = rep(seq(ncol(umi)), each = nrow(umi)),
#       value = as.vector(umi)
#     )[value <= edit ] #make a long data frame and filter according to cutoff
#     umi[, "n.1" := uc[row]$N ][
#       , "n.2" := uc[col]$N ] #add in observed freq
#     
#     
#     #adjacency singleton: keep only minor UMI with readcount 1
#     if(collapse_mode == "adjacency_singleton"){
#       umi <- umi[n.1 <= 1 | n.2 <= 1]
#     }
#     
#     #directional adjacency: keep only minor UMI with readcount 1
#     if(collapse_mode == "adjacency_directional"){
#       umi <- umi[n.1 <= n.2/2 | n.2 <= n.1/2]
#     }
#     
#     
#     umi_out <- copy(umi)
#     if(nrow(umi_out) == 0){
#       return(umi_out)
#     }
#     
#     if(collapse_mode == "cluster"){
#       umi_out[, falseUMI := ifelse( n.1>=n.2, col, row ) ][
#         , n.false := ifelse( n.1>=n.2, n.2, n.1 )]
#     }else{
#       umi_out[, falseUMI := ifelse( n.1>n.2, col, row ) ][
#         , n.false := ifelse( n.1>n.2, n.2, n.1 )]
#     }
#     umi_out[, trueUMI := ifelse( n.1<n.2, col, row ) ][
#       , n.true := ifelse( n.1<n.2, n.2, n.1 )][
#         , falseUMI := uc[falseUMI]$us ][
#           , trueUMI  := uc[trueUMI ]$us][
#             , c("row", "col", "value", "n.1", "n.2") := NULL]
#     
#     if(collapse_mode != "cluster"){
#       umi_out <- umi_out[!falseUMI == trueUMI]
#     }
#     
#     dup_daughters <- unique(umi_out[which(duplicated(falseUMI))]$falseUMI)
#     if(length(dup_daughters>0)){
#       umi_out[,rem := FALSE]
#       setorder(umi_out, falseUMI, -n.true)
#       setkey(umi_out, falseUMI)
#       for(i in dup_daughters){
#         umi_out[ i, rem := TRUE ] #remove duplicates
#         umi_out[ i, mult = "first" , rem := FALSE] # keep the most frequent parent UMI
#       }
#       umi_out <- umi_out[rem == FALSE]
#       umi_out[, rem := NULL]
#     }
#     
#     
#     non_true_UMIs <- unique(umi_out[trueUMI %in% umi_out$falseUMI]$trueUMI)
#     real_true_UMIs <- unique(umi_out[!trueUMI %in% umi_out$falseUMI]$trueUMI)
#     if(length(non_true_UMIs>0)){
#       setkey(umi_out, falseUMI)
#       for(i in non_true_UMIs){
#         true_parent_UMI <- umi_out[i][!trueUMI %in% non_true_UMIs]$trueUMI
#         if(length(true_parent_UMI)==0){#find closest match in case there is no clear parent UMI!
#           true_parent_UMI <- real_true_UMIs[stringdist::amatch(umi_out[i][1]$trueUMI, real_true_UMIs, method = "hamming", maxDist=edit, nthread = 1)[1]]
#         }
#         if(length(true_parent_UMI)>1){ #take a random good parent UMI if more possibilities exist
#           true_parent_UMI <- true_parent_UMI[1]
#         }
#         umi_out[trueUMI == i, trueUMI := true_parent_UMI]
#       }
#     }
#     umi_out[,dist := stringdist::stringdist(falseUMI,trueUMI,method = "hamming", nthread = 1)]
#     
#     if(collapse_mode != "cluster"){
#       umi_out <- umi_out[dist <= edit]
#     }
#     
#     umi_out[, c("n.false","n.true","dist") := NULL]
#     return(umi_out)
#     
#   }else{
#     return(NULL)
#   }
# }

.hammingFilter_bktree <-function(umiseq, edit=1, collapse_mode = c("adjacency","adjacency_directional","adjacency_singleton","cluster")){
  # umiseq a vector of umis, one per read
  uc <- data.table(us = umiseq)[, .N, by = "us"] # normal UMI counts
  setorder(uc, us) #order by sequence
  setkey(uc,us) #index by sequence
  
  umi <- return_dist_frame(uc$us, edit)
  if(length(umi) == 0){ # leave function as soon as possible if nothing to do.
    return(NULL)
  }
  umi <- as.data.table(matrix(unlist(lapply(umi, unlist)), byrow = T, ncol = 3)) #reformat list into correct shape
  #setDT(umi)
  setnames(umi, c("umi.1","value","umi.2"))
  umi[, value := as.integer(value)]
  
  umi[, "n.1" := uc[umi.1]$N ][
    , "n.2" := uc[umi.2]$N ] #add in observed freq
  
  #next, keep only unique pairings
  umi[, idx := seq(.N)][, pairings := paste(sort(c(umi.1,umi.2)),collapse="_"), by = idx][, idx := NULL]
  umi <- unique(umi, by="pairings")
  umi[,pairings := NULL]
  
  #adjacency singleton: keep only minor UMI with readcount 1
  if(collapse_mode == "adjacency_singleton"){
    umi <- umi[n.1 <= 1 | n.2 <= 1]
  }
  
  #directional adjacency: keep only minor UMI with readcount 1
  if(collapse_mode == "adjacency_directional"){
    umi <- umi[n.1 <= n.2/2 | n.2 <= n.1/2]
  }
  
  if(nrow(umi) == 0){
    return(umi)
  }
  
  if(collapse_mode == "cluster"){
    umi[, falseUMI := ifelse( n.1>=n.2, umi.2, umi.1 ) ][
      , n.false := ifelse( n.1>=n.2, n.2, n.1 )]
  }else{
    umi[, falseUMI := ifelse( n.1>n.2, umi.2, umi.1 ) ][
      , n.false := ifelse( n.1>n.2, n.2, n.1 )]
  }
  umi[, trueUMI := ifelse( n.1<n.2, umi.2, umi.1 ) ][
    , n.true := ifelse( n.1<n.2, n.2, n.1 )][
      , c("umi.1", "umi.2", "value", "n.1", "n.2") := NULL]
  
  umi_out <- copy(umi)
  
  if(collapse_mode != "cluster"){
    umi_out <- umi_out[!falseUMI == trueUMI]
  }
  
  if(sum(duplicated(umi_out$falseUMI))>0){ #is the collapsing ambigous?
    setorder(umi_out, falseUMI, -n.true) # order by most abundant true observation
    umi_out <- umi_out[, head(.SD,1), by = falseUMI] #select the first row
  }
  
  setorder(umi_out, -n.true)
  non_true_UMIs <- unique(umi_out[trueUMI %in% umi_out$falseUMI]$trueUMI)
  real_true_UMIs <- unique(umi_out[!trueUMI %in% umi_out$falseUMI]$trueUMI)
  if(length(non_true_UMIs>0)){
    #first, see if there is an alternative target UMI to collapse into, sorted by biggest difference between n.true & n.false
    #alternate_pairings <- umi[falseUMI %in% umi_out[trueUMI %in% non_true_UMIs]$falseUMI][trueUMI %in% real_true_UMIs][order(n.true-n.false)][,head(.SD,1), by = "falseUMI"]
    alternate_pairings <- umi[falseUMI %in% umi_out[trueUMI %in% non_true_UMIs]$falseUMI][trueUMI %in% real_true_UMIs][!(trueUMI %in% umi_out$falseUMI)][!(trueUMI %in% falseUMI)][order(n.true-n.false)][,head(.SD,1), by = "falseUMI"]
    if(nrow(alternate_pairings)>0){
      #if there are some, swap them into the collapse table to keep
      umi_out <- rbind(umi_out[!(falseUMI %in% alternate_pairings$falseUMI)], alternate_pairings)
    }
    #update if there is still work to do
    non_true_UMIs <- unique(umi_out[trueUMI %in% umi_out$falseUMI]$trueUMI)
    #real_true_UMIs <- unique(umi_out[!trueUMI %in% umi_out$falseUMI]$trueUMI)
    umi_out[, diff := n.true-n.false]
    
    if(collapse_mode == "cluster"){
      #if cluster, pull the chained true UMI identity
      chain_info <- umi_out[falseUMI %in% non_true_UMIs,c("falseUMI","trueUMI"), with = F]
      setnames(chain_info, c("trueUMI","real_trueUMI"))
      umi_out <- merge(umi_out, chain_info, by = "trueUMI", all.x = T)
      #replace the wrong true UMI with the resolved chain target
      umi_out[,trueUMI := ifelse(is.na(real_trueUMI), trueUMI, real_trueUMI)][,real_trueUMI:= NULL]
    }else{
      #if not cluster, deal with chaining differently.
      # 1. don't collapse if that UMI serves as true UMI for several other collapses.
      #umi_out <- umi_out[!(falseUMI %in% umi_out[trueUMI %in% non_true_UMIs,.N, by = "trueUMI"][N>1]$trueUMI )]
      #non_true_UMIs <- unique(umi_out[trueUMI %in% umi_out$falseUMI]$trueUMI)
      #real_true_UMIs <- unique(umi_out[!trueUMI %in% umi_out$falseUMI]$trueUMI)
      #check if there is still issues
      #if(length(non_true_UMIs)>0){
      #Next, select those collapses, that remove the nodes with least read support
      #}
      #easier solution: if not cluster, collapse the node with the least read support, which means remove all collapses that would 
      umi_out <- umi_out[!(falseUMI %in% trueUMI)]
    }
    #previously, I double checked if all edit distances are as asked by the user, but it is not necessary any more
    #umi_out[,dist := stringdist::stringdist(falseUMI,trueUMI,method = "hamming", nthread = 1)]
    #if(collapse_mode != "cluster"){
    #  umi_out <- umi_out[dist <= edit]
    #}
  }
  umi_out[, c("n.false","n.true") := NULL]
  return(umi_out)
  
}

.onLoad <- function(libname, pkgname) {
  reticulate::source_python(system.file("python/bktreeumis_ngram_nodepend.py",package="UMIcountR"),envir=globalenv())
}
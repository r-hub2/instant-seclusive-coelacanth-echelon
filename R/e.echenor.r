#-------------------------------------------------------------------------------
# Echelon scan based on Normal model
#-------------------------------------------------------------------------------
echenor <- 
  function(echelon.obj, val, weight = NULL, K = length(val)/2, Kmin = 2, n.sim = 99,
           cluster.type = "high", cluster.legend.pos = "bottomleft",
           dendrogram = TRUE, cluster.info = FALSE, coo = NULL, ...
  ){


#-------------------------------------------------------------------------------
# Check of echelon class
#-------------------------------------------------------------------------------
  if(!inherits(echelon.obj, what = "echelon")) stop(paste("The class 'echelon' is incorrect\n\n"))

  x <- echelon.obj$x
  rin <- echelon.obj$rin
  locs <- echelon.obj$locs
  peaks <- echelon.obj$peaks
  separates <- echelon.obj$separates
  c_separates <- echelon.obj$c_separates
  parents <- echelon.obj$parents


#-------------------------------------------------------------------------------
# Echelon scan
#-------------------------------------------------------------------------------
  if(Kmin < 2) stop("For the Normal model, 'Kmin' must be set to 2 or greater.\n")

  if(!is.null(weight)){
    if(length(val) != length(weight)) stop("length(val) and length(weight) must be of the same size\n")
    if(any(weight <= 0) || any(is.na(weight))) stop("Negative values or NA must not be assigned to the argument 'weight'\n")
  }
  else weight <- rep(1, times = length(val))

  if(K <= 0) K <- floor(length(x)/2)
  else if(K > length(x)) stop("Please check the argument 'K'. It must satisfy 'K' <= length(x)\n")
  else{
    K <- floor(K)
    reg_data <- e.scan(x, locs, peaks, c_separates, parents, K, par = NULL)
  }

  if(ncol(reg_data) < Kmin) stop("No clusters found! Please Try again by modifying the argument 'x', 'K', 'Kmin' or 'cluster.type'\n")
  reg_data <- reg_data[which(!is.na(reg_data[,floor(Kmin)])),]
  if(is.null(dim(reg_data))) reg_data <- t(reg_data)


#-------------------------------------------------------------------------------
# Statistic calculation
#-------------------------------------------------------------------------------

  TotalWeight <- sum(weight)
  WeightedMean <- sum(val*weight)/TotalWeight
  WeightedVar <- sum(weight*(val - WeightedMean)^2)/(TotalWeight - mean(weight))

  log.lambda <- numeric(nrow(reg_data))
  for(i in 1:nrow(reg_data)){
    Z <- reg_data[i,][!is.na(reg_data[i,])]
    WeightedMeanIn <- sum(val[Z] * weight[Z])/sum(weight[Z])
    WeightedMeanOut <- sum(val[-Z] * weight[-Z])/(TotalWeight - sum(weight[Z]))

    if(cluster.type == "high" && WeightedMeanIn < WeightedMeanOut){
      log.lambda[i] <- 0
      next
    }
    if(cluster.type == "low" && WeightedMeanIn > WeightedMeanOut){
      log.lambda[i] <- 0
      next
    }

    WeightedVarInOut <- sum(c(weight[Z] * (val[Z] - WeightedMeanIn)^2, weight[-Z]*(val[-Z] - WeightedMeanOut)^2))/(TotalWeight - mean(weight))
    log.lambda[i] <- length(val)/2 * (log(WeightedVar) - log(WeightedVarInOut))
  }


#-------------------------------------------------------------------------------
# Cluster decision
#-------------------------------------------------------------------------------
  temp <- e.cluster.decision(reg_data, log.lambda)
  cluster_reg <- temp$cluster_reg
  cluster_log.lambda <- temp$cluster_log.lambda


#-------------------------------------------------------------------------------
# Monte Carlo Estimation
#-------------------------------------------------------------------------------
  p_rank <- NULL
  if(n.sim > 0){
    n.sim <- floor(n.sim)
    cat(paste("Starting", n.sim, "Monte Carlo replications...\n"), sep="")

    nulldata <- replicate(n.sim, sample(length(x)))
    type <- ifelse(cluster.type == "high", 31, 32)

    sim_lambda <- e_monteCPP(x = nulldata, rin = rin, K = K, Kmin = floor(Kmin),
                             par1 = val, par2 = weight, type = type)
    monte_lost <- which(sim_lambda < 0)
    if(length(monte_lost) != 0) sim_lambda <- sim_lambda[-monte_lost]
    p_rank <- (n.sim + 1 - length(monte_lost)) - findInterval(cluster_log.lambda, sort(sim_lambda))
    cat("\n\n")
  }
  else{
    n.sim <- 0
    monte_lost <- NULL
  }


#-------------------------------------------------------------------------------
# Echelon dendrogram
#-------------------------------------------------------------------------------
  if(dendrogram){
    temp <- e.cluster.dendrogram(echelon.obj, n.sim, cluster.legend.pos, cluster_reg, p_rank, para = list(...))
    coord <- temp$coord
  }
  else coord <- NULL


#-------------------------------------------------------------------------------
# Cluster map
#-------------------------------------------------------------------------------
  if(!is.null(coo)){
    if(nrow(coo) != length(x)) stop("'length(x)' and 'nrow(coo)' must be same size\n\n")
    e.cluster.map(x, c_separates, locs, coo, rin, p_rank, cluster_reg, n.sim, cluster.type)
  }


#-------------------------------------------------------------------------------
# Out put
#-------------------------------------------------------------------------------
  if(cluster.info){
    cat("------------- CLUSTERS DETECTED -------------\n")
    cat(paste("Number of locations ......: ",length(x), " region\n", sep=""))

    cat(paste("Limit length of cluster ..: ", K, " regions\n", sep=""))
    if(Kmin != 1) cat(paste("Minimum length of cluster : ", Kmin, " regions\n", sep=""))

    cat(paste("Total values ..............: ", round(sum(val), digits=4), "\n", sep=""))
    if(!is.null(weight)) cat(paste("Total weights ............: ", round(TotalWeight, digits=4), "\n", sep=""))
    cat(paste("Mean .....................: ", round(mean(val), digits=4), "\n", sep=""))
    cat(paste("Variance .................: ", round(var(val), digits=4), "\n", sep=""))
    cat(paste("Standard deviation .......: ", round(sqrt(var(val)), digits=4), "\n", sep=""))
    if(!is.null(weight)) cat(paste("Weighted Mean ............: ", round(WeightedMean, digits=4),"\n", sep=""))
    if(!is.null(weight)) cat(paste("Weighted Variance ........: ", round(WeightedVar, digits=4),"\n", sep=""))
    if(!is.null(weight)) cat(paste("Weighted Std deviation....: ", round(sqrt(WeightedVar), digits=4), "\n",sep=""))


    if(cluster.type == "high") cat(paste("Scan for Area with .......: High Values\n"))
    else if(cluster.type == "low") cat(paste("Scan for Area with .......: Low Values\n"))
    cat(paste("Number of Replications ...: ", n.sim, sep=""))
    if(length(monte_lost) != 0) cat(paste("(No solution found ", length(monte_lost), " times)", sep=""))
    cat("\n")
    cat("Model ....................: Normal\n\n")
    cat("---------------------------------------------\n")
  }

  Z1 <- cluster_reg[1,][!is.na(cluster_reg[1,])]
  cat(paste("MOST LIKELY CLUSTER -- ", length(Z1), " regions\n Cluster regions included : ", sep=""))
  cat(echelon.obj$reg_name[Z1], sep=", ")
  MLC <- list(regionsID = Z1)

  cat(paste("\n Mean inside .............: ", round(mean(val[Z1]), digits=4), sep=""))
  cat(paste("\n Mean outside ............: ", round(mean(val[-Z1]), digits=4), sep=""))
  var_Z <- sum((val[Z1] - mean(val[Z1]))^2, (val[-Z1] - mean(val[-Z1]))^2)/(length(x) - 1)
  cat(paste("\n Variance ................: ", round(var_Z, digits=4), sep=""))
  cat(paste("\n Standard deviation ......: ", round(sqrt(var_Z), digits=4), sep=""))

  WeightedMeanIn <- sum(val[Z1] * weight[Z1])/sum(weight[Z1])
  WeightedMeanOut <- sum(val[-Z1] * weight[-Z1])/(TotalWeight - sum(weight[Z1]))
  WeightedVarInOut <- sum(c(weight[Z1] * (val[Z1] - WeightedMeanIn)^2, weight[-Z1]*(val[-Z1] - WeightedMeanOut)^2))/(TotalWeight - mean(weight))
  if(!is.null(weight)) cat(paste("\n Weighted mean inside ......: ", round(WeightedMeanIn, digits=4), sep=""))
  if(!is.null(weight)) cat(paste("\n Weighted mean outside .....: ", round(WeightedMeanOut, digits=4), sep=""))
  if(!is.null(weight)) cat(paste("\n Weighted variance .........: ", round(WeightedVarInOut, digits=4), sep=""))
  if(!is.null(weight)) cat(paste("\n Weighted std deviation ....: ", round(sqrt(WeightedVarInOut), digits=4), sep=""))

  cat(paste("\n Log likelihood ratio ....: ", round(cluster_log.lambda[1], digits=4), "", sep=""))
  MLC <- c(MLC, list(mean_inZ = mean(val[Z1]), mean_outZ = mean(val[Z1]), var_Z = var_Z,
    wmean_inZ = WeightedMeanIn, wmean_outZ = WeightedMeanOut, wvar_Z = WeightedVarInOut,
    LLR = cluster_log.lambda[1]))

  if(n.sim != 0){
    cat(paste("\n Monte Carlo rank ........: ", p_rank[1], "/", (n.sim + 1 - length(monte_lost)), "", sep=""))
    cat(paste("\n P-value .................: ", round(p_rank[1]/(n.sim + 1 - length(monte_lost)), digits=nchar(as.character(n.sim))+1), "", sep=""))
    MLC <- c(MLC,list(p=p_rank[1]/(n.sim+1)))
  }
  cat("\n\n")
  clusters <- MLC

  if(cluster.info) cat("----------------------------------------------\n")

  if(nrow(cluster_reg) != 1){
    if(cluster.info) cat("SECONDARY CLUSTERS\n")

    if(nrow(cluster_reg) > 5) len2C <- 5
    else len2C <- nrow(cluster_reg)

    for(i in 2:len2C){
      Z <- cluster_reg[i,][!is.na(cluster_reg[i,])]
      secondC <- NULL
      if(cluster.info){
        cat(paste(i, " -- ", length(Z), " regions\n Cluster regions included : ", sep=""))
        cat(echelon.obj$reg_name[Z], sep=", ")
      }
      secondC <- list(regionsID = Z)

      if(cluster.info){
        cat(paste("\n Mean inside .............: ", round(mean(val[Z]), digits=4), sep=""))
        cat(paste("\n Mean outside ............: ", round(mean(val[-Z]), digits=4), sep=""))
        var_Z <- sum((val[Z] - mean(val[Z]))^2, (val[-Z] - mean(val[-Z]))^2)/(length(x) - 1)
        cat(paste("\n Variance ................: ", round(var_Z, digits=4), sep=""))
        cat(paste("\n Standard deviation ......: ", round(sqrt(var_Z), digits=4), sep=""))

        WeightedMeanIn <- sum(val[Z] * weight[Z])/sum(weight[Z])
        WeightedMeanOut <- sum(val[-Z] * weight[-Z])/(TotalWeight - sum(weight[Z]))
        WeightedVarInOut <- sum(c(weight[Z] * (val[Z] - WeightedMeanIn)^2, weight[-Z]*(val[-Z] - WeightedMeanOut)^2))/(TotalWeight - mean(weight))
        if(!is.null(weight)) cat(paste("\n Weighted mean inside ....: ", round(WeightedMeanIn, digits=4), sep=""))
        if(!is.null(weight)) cat(paste("\n Weighted mean outside ...: ", round(WeightedMeanOut, digits=4), sep=""))
        if(!is.null(weight)) cat(paste("\n Weighted variance .......: ", round(WeightedVarInOut, digits=4), sep=""))
        if(!is.null(weight)) cat(paste("\n Weighted std deviation ..: ", round(sqrt(WeightedVarInOut), digits=4), sep=""))

        cat(paste("\n Log likelihood ratio ....: ", round(cluster_log.lambda[i], digits=4), "", sep=""))
      }
      secondC <- c(secondC, list(mean_inZ = mean(val[Z]), mean_outZ = mean(val[Z]), var_Z = var_Z,
        wmean_inZ = WeightedMeanIn, wmean_outZ = WeightedMeanOut, wvar_Z = WeightedVarInOut,
        LLR = cluster_log.lambda[i]))

      if(n.sim > 0){
        if(cluster.info){
          cat(paste("\n Monte Carlo rank ........: ", p_rank[i], "/", (n.sim + 1 - length(monte_lost)), "", sep=""))
          cat(paste("\n P-value .................: ", round(p_rank[i]/(n.sim + 1 - length(monte_lost)), digits=nchar(as.character(n.sim))+1), "", sep=""))
        }
        secondC <- c(secondC, list(p = p_rank[i]/(n.sim+1)))
      }
      if(i == 2) clusters <- list(clusters, secondC)
      else clusters[[i]] <- secondC
      if(cluster.info) cat("\n\n")
    }
    if(cluster.info) cat("----------------------------------------------\n")

    if(nrow(cluster_reg) > 5){
      if(cluster.info) cat("Display only the top 5 clusters. See object 'clusters' for more details\n\n")
      for(i in 6:nrow(cluster_reg)){
        Z <- cluster_reg[i,][!is.na(cluster_reg[i,])]
        secondC <- NULL
        secondC <- list(regionsID = Z)
        WeightedMeanIn <- sum(val[Z] * weight[Z])/sum(weight[Z])
        WeightedMeanOut <- sum(val[-Z] * weight[-Z])/(TotalWeight - sum(weight[Z]))
        WeightedVarInOut <- sum(c(weight[Z] * (val[Z] - WeightedMeanIn)^2, weight[-Z]*(val[-Z] - WeightedMeanOut)^2))/(TotalWeight - mean(weight))
        secondC <- c(secondC, list(mean_inZ = mean(val[Z]), mean_outZ = mean(val[Z]), var_Z = var_Z,
          wmean_inZ = WeightedMeanIn, wmean_outZ = WeightedMeanOut, wvar_Z = WeightedVarInOut,
          LLR = cluster_log.lambda[i]))
        if(n.sim != 0) secondC <- c(secondC,list(p = p_rank[i]/(n.sim + 1 - length(monte_lost))))
        clusters[[i]] <- secondC
      }
    }
  }
  else clusters <- list(clusters, "not detected")

  result <- c(clusters = list(clusters), list(scanned.regions = reg_data))
  if(n.sim > 0) result <- c(result, list(simulated.LLR = sim_lambda))
  if(dendrogram) result <- c(result, list(coord = coord, regions.value = echelon.obj$regions.value,
                             regions.name = echelon.obj$regions.name))
  invisible(result)
}

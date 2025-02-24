###########################################
### Echelon scan based on Binomial model ###
###########################################

echebin<-function(echelon.obj,cas,ctl,K=length(cas)/2,Kmin=1,n.sim=99,
									cluster.type="high",cluster.legend.pos="bottomleft",
									dendrogram=TRUE,cluster.info=FALSE,coo=NULL,...){


##############################
### Check of echelon class ###
##############################

	if(!inherits(echelon.obj, what="echelon"))
		stop(paste("The class 'echelon' is incorrect\n\n"))

	x <- echelon.obj$x
	rin <- echelon.obj$rin
	locs <- echelon.obj$locs
	peaks <- echelon.obj$peaks
	c_separates <- echelon.obj$c_separates
	parents <- echelon.obj$parents
	pop <- cas + ctl

####################
### Echelon scan ###
####################

	if(is.null(cas)||is.null(ctl)) stop("The Binomial model requires two arguments 'cas' and 'ctl'\n")
	if(K <= 0) K <- floor(length(x)/2)
	else if(K < 1) reg_data <- e.scan(x,locs,peaks,c_separates,parents,K,par=pop)
	else if(K > length(x)) stop("Please check the argument 'K'. It must satisfy 'K' <= length(x)\n")
	else{
		K <- floor(K)
		reg_data <- e.scan(x,locs,peaks,c_separates,parents,K, par=NULL)
	}

	if(Kmin < 1) Kmin <- 1
	if(Kmin == 1){
		if(is.null(reg_data)) stop("No clusters found! Please Try again by modifying the argument 'x', 'K' or 'cluster.type'\n")
	}
	else{
		if(ncol(reg_data) < Kmin) stop("No clusters found! Please Try again by modifying the argument 'x', 'K', 'Kmin' or 'cluster.type'\n")
		reg_data <- reg_data[which(!is.na(reg_data[,floor(Kmin)])),]
		if(is.null(dim(reg_data))) reg_data <- t(reg_data)
	}

#############################
### Statistic calculation ###
#############################

	casg <- sum(cas)
	ctlg <- sum(ctl)
	popg <- casg + ctlg

	casz <- apply(array(cas[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
	ctlz <- apply(array(ctl[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
	popz <- casz + ctlz
	log.lambda <- (casz*log(casz/popz) + ctlz*log(ctlz/popz) + (casg-casz)*log((casg-casz)/(popg-popz))
		+ (ctlg-ctlz)*log((ctlg-ctlz)/(popg-popz)) -casg*log(casg/popg) -ctlg*log(ctlg/popg))

	if(cluster.type == "high") log.lambda[which(casz/popz < casg/popg)] <- 0
	if(cluster.type == "low") log.lambda[which(casz/popz > casg/popg)] <- 0


########################
### Cluster decision ###
########################

	temp <- e.cluster.decision(reg_data,log.lambda)
	cluster_reg <- temp$cluster_reg
	cluster_log.lambda <- temp$cluster_log.lambda


##############################
### Monte Carlo Estimation ###
##############################

	p_rank <- NULL
	if(n.sim > 0){
		n.sim <- floor(n.sim)
		cat(paste("Starting",n.sim,"Monte Carlo replications...\n"),sep="")

		nulldata <- rmultinom(n.sim,round(sum(cas)),prob=pop)
		type <- ifelse(cluster.type == "high", 21, 22)

		sim_lambda <- e_monteCPP(x=nulldata,rin=rin, K=K, Kmin = floor(Kmin), par1=numeric(length(x)), par2=pop, type=type)
		monte_lost <- which(sim_lambda < 0)
		if(length(monte_lost) != 0) sim_lambda <- sim_lambda[-monte_lost]
		p_rank <- (n.sim + 1 - length(monte_lost)) - findInterval(cluster_log.lambda, sort(sim_lambda))
	}
	else{
		n.sim <- 0
		monte_lost <- NULL
	}


##########################
### Echelon dendrogram ###
##########################

	if(dendrogram){
		temp <- e.cluster.dendrogram(echelon.obj,n.sim,cluster.legend.pos,cluster_reg,p_rank,para=list(...))
		coord <- temp$coord
	}
	else coord <- NULL


###################
### Cluster map ###
###################

	if(!is.null(coo)){
		if(nrow(coo) != length(x)) stop("length(x) and nrow(coo) must have same size\n\n")
		e.cluster.map(x,c_separates,locs,coo,rin,p_rank,cluster_reg,n.sim,cluster.type)
	}


###############
### Out put ###
###############

	if(cluster.info){
		cat("------------- CLUSTERS DETECTED -------------\n")
		cat(paste("Number of locations ......: ",length(x)," region\n",sep=""))

		if(K >= 1) cat(paste("Limit length of cluster ..: ",K," regions\n",sep=""))
		else cat(paste("Limit length of cluster ..: ",K*100," percent of population\n",sep=""))
		if(Kmin != 1) cat(paste("Minimum length of cluster : ",Kmin," regions\n",sep=""))

		cat(paste("Total cases ..............: ",sum(cas),"\n",sep=""))
		cat(paste("Total population .........: ",sum(pop),"\n",sep=""))
		if(cluster.type == "high") cat(paste("Scan for Area with .......: High Rates\n"))
		if(cluster.type == "low") cat(paste("Scan for Area with .......: Low Rates\n"))
		cat(paste("Number of Replications ...: ",n.sim,"\n",sep=""))
		if(length(monte_lost) != 0) cat(paste("(No solution found ",length(monte_lost)," times)",sep=""))
		cat("\n")
		cat("Model ....................: Binomial\n\n")
		cat("---------------------------------------------\n")
	}

	cat(paste("MOST LIKELY CLUSTER -- ",length(cluster_reg[1,][!is.na(cluster_reg[1,])])," regions\n Cluster regions included : ",sep=""))
	cat(echelon.obj$reg_name[cluster_reg[1,][!is.na(cluster_reg[1,])]],sep=", ")
	MLC <- list(regionsID=cluster_reg[1,][!is.na(cluster_reg[1,])])

	pop_inZ <- sum(pop[cluster_reg[1,]],na.rm=TRUE)
	cat(paste("\n Population ..............: ",pop_inZ,sep=""))
	MLC <- c(MLC,list(pop_inZ=pop_inZ))

	cas_inZ <- sum(cas[cluster_reg[1,]],na.rm=TRUE)
	cas_outZ <- casg-cas_inZ
	ex_inZ <- sum(pop[cluster_reg[1,]]*casg/popg,na.rm=TRUE)
	ex_outZ <- casg-ex_inZ

	cat(paste("\n Number of cases .........: ",cas_inZ,sep=""))
	cat(paste("\n Expected cases ..........: ",round(ex_inZ,digits=4),sep=""))
	cat(paste("\n Observed / expected .....: ",round(cas_inZ/ex_inZ,digits=4),sep=""))
	cat(paste("\n Relative risk ...........: ",round((cas_inZ/ex_inZ)/(cas_outZ/ex_outZ),digits=4),sep=""))
	cat(paste("\n Log likelihood ratio ....: ",round(cluster_log.lambda[1],digits=4),"",sep=""))
	MLC <- c(MLC,list(cas_inZ=cas_inZ,ex_inZ=ex_inZ,LLR=cluster_log.lambda[1]))

	if(n.sim != 0){
		cat(paste("\n Monte Carlo rank ........: ",p_rank[1],"/",(n.sim + 1 - length(monte_lost)),"",sep=""))
		cat(paste("\n P-value .................: ",round(p_rank[1]/(n.sim + 1 - length(monte_lost)),digits=nchar(as.character(n.sim))+1),"",sep=""))
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
			secondC <- NULL
			if(cluster.info){
			  cat(paste(i," -- ",length(cluster_reg[i,][!is.na(cluster_reg[i,])])," regions\n Cluster regions included : ",sep=""))
				cat(echelon.obj$reg_name[cluster_reg[i,][!is.na(cluster_reg[i,])]],sep=", ")
			}
			secondC <- list(regionsID=cluster_reg[i,][!is.na(cluster_reg[i,])])

			if(!is.null(pop)){
				pop_inZ <- sum(pop[cluster_reg[i,]],na.rm=TRUE)
				if(cluster.info) cat(paste("\n Population ..............: ",pop_inZ,sep=""))
				secondC <- c(secondC,list(pop_inZ=pop_inZ))
			}

			cas_inZ <- sum(cas[cluster_reg[i,]],na.rm=TRUE)
			cas_outZ <- casg-cas_inZ
			ex_inZ <- sum(pop[cluster_reg[i,]]*casg/popg,na.rm=TRUE)
			ex_outZ <- casg-ex_inZ

			if(cluster.info){
			  cat(paste("\n Number of cases .........: ",cas_inZ,sep=""))
				cat(paste("\n Expected cases ..........: ",round(ex_inZ,digits=4),sep=""))
				cat(paste("\n Observed / expected .....: ",round(cas_inZ/ex_inZ,digits=4),sep=""))
				cat(paste("\n Relative risk ...........: ",round((cas_inZ/ex_inZ)/(cas_outZ/ex_outZ),digits=4),sep=""))
				cat(paste("\n Log likelihood ratio ....: ",round(cluster_log.lambda[i],digits=4),"",sep=""))
			}
			secondC <- c(secondC,list(cas_inZ=cas_inZ,ex_inZ=ex_inZ,LLR=cluster_log.lambda[i]))

			if(n.sim > 0){
			  if(cluster.info){
			    cat(paste("\n Monte Carlo rank ........: ",p_rank[i],"/",(n.sim + 1 - length(monte_lost)),"",sep=""))
					cat(paste("\n P-value .................: ",round(p_rank[i]/(n.sim + 1 - length(monte_lost)),digits=nchar(as.character(n.sim))+1),"",sep=""))
				}
			  secondC <- c(secondC,list(p=p_rank[i]/(n.sim+1)))
			}
			if(i == 2) clusters <- list(clusters,secondC)
			else clusters[[i]] <- secondC
			if(cluster.info) cat("\n\n")
		}
		if(cluster.info) cat("----------------------------------------------\n")

		if(nrow(cluster_reg) >5){
		  if(cluster.info) cat("Display only the top 5 clusters. See object 'clusters' for more details\n\n")
			for(i in 6:nrow(cluster_reg)){
				secondC <- NULL
				secondC <- list(regionsID=cluster_reg[i,][!is.na(cluster_reg[i,])])
				if(!is.null(pop)){
					pop_inZ <- sum(pop[cluster_reg[i,]],na.rm=TRUE)
					secondC <- c(secondC,list(pop_inZ=pop_inZ))
				}
				cas_inZ <- sum(cas[cluster_reg[i,]],na.rm=TRUE)
				ex_inZ <- sum(pop[cluster_reg[i,]]*casg/popg,na.rm=TRUE)
				secondC <- c(secondC,list(cas_inZ=cas_inZ,ex_inZ=ex_inZ,LLR=cluster_log.lambda[i]))
				if(n.sim != 0){
					secondC <- c(secondC,list(p=p_rank[i]/(n.sim + 1 - length(monte_lost))))
				}
				clusters[[i]] <- secondC
			}
		}
	}
	else clusters <- list(clusters, "not detected")

	result <- c(clusters=list(clusters),list(scanned.regions=reg_data))
	if(n.sim > 0) result <- c(result,list(simulated.LLR=sim_lambda))
	if(dendrogram) result <- c(result,list(coord=coord,regions.value=echelon.obj$regions.value,
															regions.name=echelon.obj$reg_name[locs]))
	invisible(result)
}

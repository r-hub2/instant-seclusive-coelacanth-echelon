#############################################
### Main and scan program for Monte Carlo ###
#############################################

e.main.monte <- function(x, rin, K, Kmin, par1, par2, type){
	
	if(type == 11){
		cas <- x
		ex <- par1
		x <- x/ex
	}
	else if(type == 12){
		cas <- x
		ex <- par1
		x <- -x/ex
	}
	else if(type == 21){
		cas <- x
		ctl <- par2 - cas
		x <- x/ctl
	}
	else if(type == 22){
		cas <- x
		ctl <- par2 - cas
		x <- -x/ctl
	}
	dat <- x

### Echelon analysis ###
	total_eche_locs <- numeric(length(dat))
	include_eche_locs <- NULL
	eche_separates <- NULL
	eche_separates_include <- NULL
	eche_peaks <- NULL
	pare_temp <- NULL

	while(any(!is.na(dat))){
		now_val <- max(dat,na.rm=TRUE)
		now_loc <- which(dat==now_val)[1]
		
		t_dummy <- now_loc

		now_eche_locs <- now_loc
		p_dummy <- 1

		while(!is.na(t_dummy)){
			now_nei <- as.vector(rin[now_loc,])
			now_nei <- now_nei[!is.na(now_nei)]
			if(any(is.element(now_nei, now_loc))) now_nei <- now_nei[-which(is.element(now_nei, now_loc))]
			now_nei <- unique(now_nei)

			if(any(is.element(total_eche_locs, now_nei))){
				temp <- rep(seq_along(eche_separates_include), times = eche_separates_include)
				del_eche_separates_include <- unique(temp[is.element(include_eche_locs, now_nei)])
				now_loc <- c(include_eche_locs[is.element(temp, del_eche_separates_include)], now_loc)

				include_eche_locs <- include_eche_locs[!is.element(include_eche_locs, now_loc)]
				include_eche_locs <- c(include_eche_locs, now_loc)

				if(p_dummy == 1) eche_separates_include <- eche_separates_include[-del_eche_separates_include]
				else eche_separates_include <- eche_separates_include[-c(del_eche_separates_include, length(eche_separates_include))]
				eche_separates_include <- c(eche_separates_include, length(now_loc))

				now_nei <- as.vector(rin[now_loc,])
				now_nei <- now_nei[!is.na(now_nei)]
				if(any(is.element(now_nei, now_loc))) now_nei <- now_nei[-which(is.element(now_nei, now_loc))]
				now_nei <- unique(now_nei)

				p_dummy <- p_dummy + length(del_eche_separates_include)
		  }

			if(any(x[now_nei] == x[t_dummy])){
				now_loc <- c(now_loc, now_nei[which(x[now_nei] == x[t_dummy])])
				now_eche_locs <- c(now_eche_locs, now_nei[which(x[now_nei] == x[t_dummy])])
				
				if(p_dummy != 1){
					include_eche_locs <- c(include_eche_locs, now_nei[which(x[now_nei] == x[t_dummy])])
					eche_separates_include[length(eche_separates_include)] <- eche_separates_include[length(eche_separates_include)] + length(now_nei[which(x[now_nei] == x[t_dummy])])
				}
			}
			else t_dummy <- NA
		}

		i_dummy <- now_eche_locs
		s_dummy <- 0

		while(s_dummy == 0){
			now_nei <- as.vector(rin[now_loc,])
			now_nei <- now_nei[!is.na(now_nei)]
			if(any(is.element(now_nei, now_loc))) now_nei <- now_nei[-which(is.element(now_nei, now_loc))]
			now_nei <- unique(now_nei)

			if(length(now_nei) == 0){
				dat[now_eche_locs] <- NA
				total_eche_locs[which(total_eche_locs==0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs
				if(p_dummy == 1) eche_peaks <- c(eche_peaks, 1)
				else eche_peaks <- c(eche_peaks, 0)
				eche_separates <- c(eche_separates, length(now_eche_locs))
				now_eche_locs <- NULL
				pare_temp <- c(pare_temp, 0)
				break
			}

			now_nei_max_val <- max(dat[now_nei])
			now_nei_max_locs <- now_nei[which(now_nei_max_val == dat[now_nei])]
			t_dummy <- now_nei_max_val

			while(!is.na(t_dummy)){
				sub_nei <- c(now_nei,as.vector(rin[now_nei_max_locs,]))
				sub_nei <- sub_nei[!is.na(sub_nei)]
				sub_nei <- sub_nei[-which(is.element(sub_nei, c(now_loc, now_nei_max_locs)))]
				sub_nei <- unique(sub_nei)

				if(any(x[sub_nei] == t_dummy)){
					now_nei <- c(now_nei, sub_nei[x[sub_nei] == t_dummy])
					now_nei_max_locs <- c(now_nei_max_locs, sub_nei[x[sub_nei] == t_dummy])
				}
				
				else t_dummy <- NA
			}

			if(length(sub_nei) == 0){
				now_eche_locs <- c(now_eche_locs, now_nei_max_locs)
				dat[now_eche_locs] <- NA
				total_eche_locs[which(total_eche_locs==0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs
				if(p_dummy == 1){
					eche_peaks <- c(eche_peaks, 1)
					include_eche_locs <- c(include_eche_locs, now_eche_locs)
					eche_separates_include <- c(eche_separates_include, length(now_eche_locs))
				}
				else{
					eche_peaks <- c(eche_peaks, 0)
					include_eche_locs <- c(include_eche_locs, now_eche_locs[-c(seq_along(i_dummy))])
					eche_separates_include[length(eche_separates_include)] <- eche_separates_include[length(eche_separates_include)] + length(now_eche_locs[-c(seq_along(i_dummy))])
				}
				eche_separates <- c(eche_separates, length(now_eche_locs))
				now_eche_locs <- NULL
				pare_temp <- c(pare_temp, 0)
				break
			}

			else if(now_nei_max_val >= max(x[sub_nei])){
				now_eche_locs <- c(now_eche_locs, now_nei_max_locs)
				now_loc <- unique(c(now_loc, now_eche_locs))
			}

			else{
				total_eche_locs[which(total_eche_locs==0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs
				dat[now_eche_locs] <- NA
				pare_temp <- c(pare_temp, now_nei_max_locs[1])

				if(p_dummy == 1){
					eche_peaks <- c(eche_peaks, 1)
					include_eche_locs <- c(include_eche_locs, now_eche_locs)
					eche_separates_include <- c(eche_separates_include, length(now_eche_locs))
				}
				else{
					eche_peaks <- c(eche_peaks, 0)
					include_eche_locs <- c(include_eche_locs, now_eche_locs[-c(seq_along(i_dummy))])
					eche_separates_include[length(eche_separates_include)] <- eche_separates_include[length(eche_separates_include)] + length(now_eche_locs[-c(seq_along(i_dummy))])
				}
				eche_separates <- c(eche_separates, length(now_eche_locs))
				now_eche_locs <- NULL
				s_dummy <- 1
			}
		}
	}

	eche_parent <- rep(seq_along(eche_separates), times = eche_separates)[match(pare_temp, total_eche_locs)]
	eche_parent[is.na(eche_parent)] <- 0
	

### Echelon scan ###
	e_num <- 1
	reg_data <- NULL
	Kout <- 0
	c_separates <- c(0,cumsum(eche_separates))

### K>=1 ###
	if(K >= 1){
		while(e_num <= length(eche_peaks)){
			if(eche_peaks[e_num] == 1){
				areas <- total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]

				if(length(areas) > K){
					i <- 1
					temp <- NULL

					while(length(temp) < K){
						t_check <- x[areas[i]]
						t_check2 <- which(x[areas] == t_check)
						if(length(temp)+length(t_check2) <= K){
							temp <- c(temp,areas[t_check2])
							i <- length(temp)+1
						}
						else{
						  if(is.null(temp)) Kout <- 1
							break
						}
					}
					if(!is.null(temp)) areas <- temp
				}

				if(Kout == 0){
					temp <- t(array(areas,c(length(areas),length(areas))))
					temp[upper.tri(temp)] <- NA

					if(any(duplicated(x[areas]))){
						del_tie <- which(duplicated(x[areas]))-1
						temp <- temp[-del_tie,,drop=FALSE]
					}

					if(e_num == 1){
						reg_data <- temp
						reg_data_col <- ncol(reg_data)
					}
					else if(ncol(temp) == reg_data_col) reg_data <- rbind(reg_data,temp)
					else if(ncol(temp) > reg_data_col){
						reg_data <- cbind(reg_data,array(NA,dim=c(nrow(reg_data),ncol(temp)-reg_data_col)))
						reg_data <- rbind(reg_data,temp)
						reg_data_col <- ncol(reg_data)
					}
					else{
						temp <- cbind(temp,array(NA,dim=c(nrow(temp),reg_data_col-ncol(temp))))
						reg_data <- rbind(reg_data,temp)
					}
				}
				else Kout <- 0
			}

			else{
				areas <- total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]
				now_eches <- e_num
				temp <- NULL

				while(any(eche_peaks[now_eches] == 0)){
					child_eches <- which(is.element(eche_parent,now_eches))
					for(i in 1:length(child_eches)) temp <- c(temp,total_eche_locs[(c_separates[child_eches[i]]+1):c_separates[child_eches[i]+1]])

					e_top <- x[total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]]
					if((length(temp)+length(which(e_top[1] == e_top))) > K){
						Kout <- 1
						break
					}
					now_eches <- child_eches
				}

				if(Kout == 0){
					if(length(temp)+length(areas) > K){
						i <- 1
						temp2 <- NULL
						while(length(temp)+length(temp2) <= K){
							t_check <- x[areas[i]]
							t_check2 <- which(x[areas] == t_check)
							if(length(temp)+length(temp2)+length(t_check2) <= K){
								temp2 <- c(temp2,areas[t_check2])
								i <- length(temp2)+1
							}
							else break
						}
						areas <- temp2
					}

					temp2 <- t(array(temp,c(length(temp),length(areas))))
					temp <- t(array(areas,c(length(areas),length(areas))))
					temp[upper.tri(temp)] <- NA
					temp <- cbind(temp2,temp)

					if(any(duplicated(x[areas]))){
						del_tie <- which(duplicated(x[areas]))-1
						temp <- temp[-del_tie,,drop=FALSE]
					}
					if(ncol(temp) == reg_data_col) reg_data <- rbind(reg_data,temp)
					else if(ncol(temp) > reg_data_col){
						reg_data <- cbind(reg_data,array(NA,dim=c(nrow(reg_data),ncol(temp)-reg_data_col)))
						reg_data <- rbind(reg_data,temp)
						reg_data_col <- ncol(reg_data)
					}
					else{
						temp <- cbind(temp,array(NA,dim=c(nrow(temp),reg_data_col-ncol(temp))))
						reg_data <- rbind(reg_data,temp)
					}
				}
				else Kout<-0
			}
			e_num <- e_num+1
		}
	}
### K<1 ###
	else if(K < 1){
		ng <- sum(par2)

		while(e_num <= length(eche_peaks)){
			if(eche_peaks[e_num] == 1){
				areas <- total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]

				if(sum(par2[areas]) > ng*K){
					i<-1
					temp <- NULL
					while(sum(par2[temp]) <= ng*K){
						t_check <- x[areas[i]]
						t_check2 <- which(x[areas] == t_check)
						if(sum(par2[temp])+sum(par2[areas[t_check2]]) <= ng*K){
							temp <- c(temp,areas[t_check2])
							i <- length(temp)+1
						}
						else{
					  	if(is.null(temp)) Kout <-1
							break
						}
					}
					if(!is.null(temp)) areas <- temp
				}

				if(Kout == 0){
					temp <- t(array(areas,c(length(areas),length(areas))))
					temp[upper.tri(temp)] <- NA

					if(any(duplicated(x[areas]))){
						del_tie <- which(duplicated(x[areas]))-1
						temp <- temp[-del_tie,,drop=FALSE]
					}

					if(is.null(reg_data)){
						reg_data <- temp
						reg_data_col <- ncol(reg_data)
					}
					else if(ncol(temp) == reg_data_col) reg_data <- rbind(reg_data,temp)
					else if(ncol(temp) > reg_data_col){
						reg_data <- cbind(reg_data,array(NA,dim=c(nrow(reg_data),ncol(temp)-reg_data_col)))
						reg_data <- rbind(reg_data,temp)
						reg_data_col <- ncol(reg_data)
					}
					else{
						temp <- cbind(temp,array(NA,dim=c(nrow(temp),reg_data_col-ncol(temp))))
						reg_data <- rbind(reg_data,temp)
					}
				}
				else Kout <- 0
			}

			else{
				areas <- total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]
				now_eches <- e_num
				temp <- NULL

				while(any(eche_peaks[now_eches] == 0)){
					child_eches <- which(is.element(eche_parent,now_eches))
					for(i in 1:length(child_eches)) temp <- c(temp,total_eche_locs[(c_separates[child_eches[i]]+1):c_separates[child_eches[i]+1]])

					e_top <- x[total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]]
					if(sum(par2[temp])+sum(par2[areas[which(e_top[1]==e_top)]]) > ng*K){
						Kout <- 1
						break
					}

					now_eches <- child_eches
				}

				if(Kout == 0){
					if(sum(par2[temp])+sum(par2[areas]) > ng*K){
						i <- 1
						temp2 <- NULL
						while(sum(par2[temp])+sum(par2[temp2]) <= ng*K){
							t_check <- x[areas[i]]
							t_check2 <- which(x[areas] == t_check)
							if(sum(par2[temp])+sum(par2[temp2])+sum(par2[areas[t_check2]]) <= ng*K){
								temp2 <- c(temp2,areas[t_check2])
								i <- length(temp2)+1
							}
							else break
						}
						areas <- temp2
					}

					temp2 <- t(array(temp,c(length(temp),length(areas))))
					temp <- t(array(areas,c(length(areas),length(areas))))
					temp[upper.tri(temp)] <- NA
					temp <- cbind(temp2,temp)

					if(any(duplicated(x[areas]))){
						del_tie <- which(duplicated(x[areas]))-1
						temp <- temp[-del_tie,,drop=FALSE]
					}

					if(ncol(temp) == reg_data_col) reg_data <- rbind(reg_data,temp)
					else if(ncol(temp) > reg_data_col){
						reg_data <- cbind(reg_data,array(NA,dim=c(nrow(reg_data),ncol(temp)-reg_data_col)))
						reg_data <- rbind(reg_data,temp)
						reg_data_col <- ncol(reg_data)
					}
					else{
						temp <- cbind(temp,array(NA,dim=c(nrow(temp),reg_data_col-ncol(temp))))
						reg_data <- rbind(reg_data,temp)
					}
				}
				else Kout <- 0
			}
			e_num <- e_num+1
		}
	}

#	list(reg_data=reg_data)
	if(!is.null(reg_data)){
		if(is.vector(reg_data)) reg_data <- t(reg_data)
		if(ncol(reg_data) < Kmin) reg_data <- NULL
		else if(Kmin != 1) reg_data <- reg_data[which(!is.na(reg_data[,Kmin])),]
		if(is.vector(reg_data)) reg_data <- t(reg_data)
	
### LLR ###
		if(!is.null(reg_data)){

### Poisson ###
			if(type == 11 || type == 12){
				cg <- sum(cas)
				eg <- sum(ex)
				cz <- apply(array(cas[reg_data], dim(reg_data)),1,sum,na.rm=TRUE)
				ez <- apply(array(ex[reg_data], dim(reg_data)),1,sum,na.rm=TRUE)

				if(type == 11) temp <- which(cz>ez)
				else if(type == 12) temp <- which(cz<ez)

				cz <- cz[temp]
				ez <- ez[temp]
				if(length(cz) == 0) maxLLR <- 0
				else maxLLR <- max(cz*log(cz/ez)+(cg-cz)*log((cg-cz)/(eg-ez)), na.rm=TRUE)
			}

### Binomial ###
			if(type == 21 || type == 22){
				casg <- sum(cas)
				ctlg <- sum(ctl)
				casz <- apply(array(cas[reg_data], dim(reg_data)),1,sum,na.rm=TRUE)
				ctlz <- apply(array(ctl[reg_data], dim(reg_data)),1,sum,na.rm=TRUE)

				if(type == 21) temp <- which(casz/(casz + ctlz) > casg/(casg + ctlg))
				else if(type == 22) temp <- which(casz/(casz + ctlz) < casg/(casg + ctlg))

				casz <- casz[temp]
				if(length(casz) == 0) maxLLR <- 0
				else{
					ctlz <- ctlz[temp]
					popz <- casz + ctlz
					popg <- casg + ctlg

					maxLLR <- max(casz*log(casz/popz) + ctlz*log(ctlz/popz) + (casg-casz)*log((casg-casz)/(popg-popz))
						+ (ctlg-ctlz)*log((ctlg-ctlz)/(popg-popz)) -casg*log(casg/popg) -ctlg*log(ctlg/popg), na.rm=TRUE)
				}
			}
### Normal ###

		}
		else maxLLR <- -1
	}
	else maxLLR <- -1

	return(maxLLR)
}


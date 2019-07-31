

##################################################################################



	
	
	################################################################################
	#### Plot observed and fitted observations

	#res_fits		<- extract.fits(file=paste(rootdir,outputprefix,"output.log",sep=""), path="") 
	#names(res_fits)
	#plot_fits_in_trellis(dat=res_fits, dat_name="f1_LLs1r1_catchA", xlab="Age", ylab="Proportion")
	#plot_fits_in_trellis(dat=res_fits, dat_name="f1_LLs1r1_catchS", xlab="Size", ylab="Proportion")
	
	# unpaste for R (from NZ CASAL code ('MCMC plot functions.s)
	unpaste <- function (str, sep = "/") {
		w <- strsplit(str, sep)
		w <- matrix(unlist(w), ncol = length(str))
		nr <- nrow(w)
		ans <- vector("list", nr)
		for (j in 1:nr) ans[[j]] <- w[j, ]
		ans
	}
	
	#####################################
	#### Bubble plots of age observations		
	plot_YCS_bubbles <- function(est_file, casal_path, yline, dat_name) {
		dat	<- extract.csl.file(est_file,casal_path)
		dat_df <- NULL
		
		for (i in 1:length(dat_name)) {
			if(dat_name[i]=="Surv1A")  {
				a1	<- dat[[paste("relative_numbers_at[",dat_name[i],"]",sep="")]]
				if(is.null(a1)) a1	<- dat[[paste("proportions_at[",dat_name[i],"]",sep="")]]
			} else { 
				a1	<- dat[[paste("catch_at[",dat_name[i],"]",sep="")]]
			}
			if(!is.null(a1)) {
				as.numeric(a1$years)
				a	<- lapply(a1[a1$years],as.numeric)	# For this to work, data must be space separated (not tab!)
				a	<- lapply(a,function(x) x/sum(x))
				a	<- unlist(a)
				age	<- as.numeric(a1$"min_class"):as.numeric(a1$"max_class")
				a1	<- data.frame(year=as.numeric(substring(names(a),1,4)),age=rep(age,length(a1$years)),n=a,row.names=NULL)
				a1$year.class	<- a1$year-a1$age-1
				a1$dat			<- dat_name[i]		# Name of data
				if(dat_name[i] %in% c("Catch_LL1Aearly")) a1$dat <- "Catch_LL1A"
				if(dat_name[i] %in% c("Catch_LL2Aearly")) a1$dat <- "Catch_LL2A"
				a1$small		<- unlist(tapply(a1$n,a1$age,function(x) ifelse(x > median(x), F, T)))	# Obs smaller/bigger than median
				dat_df <- rbind(dat_df,a1)
			}
		}
		
		dat_name1 <- dat_name[! dat_name %in% c("Catch_LL1Aearly","Catch_LL2Aearly")]
		d1 <- dat_df[dat_df$dat %in% dat_name1[1],]		

		par(xaxs="i",yaxs="i")
		xlim <- c(min(dat_df$year)-0.5,max(dat_df$year)+1)
		# ylim <- c(min(dat_df$year.class)-0.5,max(dat_df$year.class)+1)
		ylim <- c(1978.5,max(dat_df$year.class)+1)
		plot(d1$year,d1$year.class, xlab="Year", ylab="Year class",xaxt="n",type="n",xlim=xlim,ylim=ylim,font.lab=2)	
		axis(1,at=dat_df$year,label=dat_df$year)
		abline(h=c(1980,2010,yline),lty=8,col="grey")
		for (i in 1:length(dat_name1)) {
			d1 <- dat_df[dat_df$dat %in% dat_name1[i],]
			if(nrow(d1)>0) {
				xshift <- 0.17 * (i-1)				
				if(i==1) symbols(d1$year+xshift,d1$year.class,circles=sqrt(d1$n),inches=0.1,add=T,fg="#80808080", bg="#FF000080")	# Survey = red
				if(i==2) symbols(d1$year+xshift,d1$year.class,circles=sqrt(d1$n),inches=0.1,add=T,fg="#80808080", bg="#0000FF80")	# Trawl  = blue
				if(i==3) symbols(d1$year+xshift,d1$year.class,circles=sqrt(d1$n),inches=0.1,add=T,fg="#80808080", bg="#0000FF80")	# Trawl  = blue
				if(i==4) symbols(d1$year+xshift,d1$year.class,circles=sqrt(d1$n),inches=0.1,add=T,fg="#80808080", bg="#80808080")	# LL1    = grey	
				if(i==5) symbols(d1$year+xshift,d1$year.class,circles=sqrt(d1$n),inches=0.1,add=T,fg="#80808080", bg="#80808080")	# LL2    = grey	
				if(i==6) symbols(d1$year+xshift,d1$year.class,circles=sqrt(d1$n),inches=0.1,add=T,fg="#80008080", bg="#80008080")	# Trap   = purple	
				
				# symbols(d1$year+xshift,d1$year.class,circles=sqrt(d1$n),inches=0.1,add=T,fg="#80808080", bg="#80808080")
				# symbols(d1$year+xshift,d1$year.class,circles=sqrt(d1$n),inches=0.1,add=T,fg="#80808080", bg=ifelse(a1$small,"#FF000080","#0000FF80"))  
			}
			if( length(dat_name1)==4) legend(y=ylim[2]-4,x=min(dat_df$year),pt.bg=c("#FF000080","#0000FF80","#80808080"),col=c("red","blue","grey"),# bg="white",
											pch=21,bty="l",box.col="white", legend=c("Survey","Trawl","Longline"),cex=1.2)
			if( length(dat_name1)==5) legend(y=ylim[2]-4,x=min(dat_df$year),pt.bg=c("#FF000080","#0000FF80","#80808080","#80008080"),col=c("red","blue","grey","purple"),# bg="white",
											pch=21,bty="l",box.col="white", legend=c("Survey","Trawl","Longline","Trap"),cex=1.2)
		}	
	}
	

	###########################################################################
	plot_fits_survey <- function(dat, dat_name, xlab, ylab, mcmcdat=NULL) {
		obs 	<- dat[[dat_name]]$obs
		pred 	<- dat[[dat_name]]$fits
		yr		<- dat[[dat_name]]$year
		cv		<- dat[[dat_name]]$error.value
		sigma2 	<- log(cv^2 + 1)				# sigma2 = log(CV^2+1) for lognormal
		mu		<- log(obs)-sigma2/2			# Mean of transformed distribution
		obs_upp	<- exp(mu + 1.96 * sqrt(sigma2))	# Quantiles (& percentiles) are preserved under transformations like exp(x) & ln(x) 
		obs_low	<- exp(mu - 1.96 * sqrt(sigma2))	# so 5% & 95% of X will be the 55 & 95% when transformed through exp.
		# obs_upp1 	<- qlnorm(0.975, meanlog=log(obs), sdlog=sqrt(sigma2))	# Lognormal
        # obs_low1	<- qlnorm(0.025, meanlog=log(obs), sdlog=sqrt(sigma2))
		data_df <- data.frame(pred=pred, obs=obs, year=yr, obs_upp, obs_low) #, obs_upp1, obs_low1)
		# If there is any MCMC data: calculate CI
		if(!is.null(mcmcdat) & length(mcmcdat>0)) {	
			ci		<- apply(mcmcdat, 2, quantile, c(0.025,0.5,0.975))	
			nn		<- unpaste(dimnames(ci)[[2]],sep="\\[")
			nn[[2]]	<- as.numeric(sub(pattern="]",replacement="",x=nn[[2]]))	# Remove ']'
			nn[[3]]	<- as.numeric(sub(pattern="]",replacement="",x=nn[[3]]))	# Remove ']'
			ci1		<- as.data.frame(cbind(t(ci),Year=nn[[2]],Age=nn[[3]]))
			ci1		<- ci1[order(ci1$Age),]
			data_df <- cbind(data_df,ci1)
			#data_df$pred1	<- -100000			# Make MPD fit negative so it will not be visible in the plots
		} else {
			data_df <- cbind(data_df,'2.5%'=-100000,'50%'=-100000,'97.5%'=-100000)	# Make it negative so it will not be visible in the plots
		}
		## Plot
		par(mfrow=c(1,1)); par(mar=c(4,4,2,1))
		maxY 	<- 1.1*max(c(data_df$pred,data_df$obs_upp,data_df$'97.5%'))
		plot(y=data_df$obs, x=data_df$year, ylim=c(0,max(maxY)), ylab=ylab, xlab=xlab, lwd=1,
				type="b",col="black", bg="black", pch=19, font.lab=2,xaxt="n")
		for (ii in 1:length(data_df$year)) 
			segments(data_df$year[ii],data_df$obs_low[ii],data_df$year[ii],data_df$obs_upp[ii],col="black",lwd=1)
		# for (ii in 1:length(data_df$year)) 
		# 	segments(data_df$year[ii]+0.1,data_df$obs_low1[ii],data_df$year[ii]+0.1,data_df$obs_upp1[ii],col="grey",lwd=1)
		axis(side=1,at=min(data_df$year):max(data_df$year),labels=min(data_df$year):max(data_df$year))
		# MCMC fits
		polygon(x=c(data_df$year,rev(data_df$year)),y=c(data_df$'2.5%',rev(data_df$'97.5%')),
				  border="red", col="#FF000080", lty="solid",lwd=1)
		lines(y=data_df$pred, x=data_df$year,col="red",bg="red",pch=19,type="b",lwd=1)	
		legend(x=max(data_df$year)-3,y=0.95*maxY,legend=c("Observed","MPD fits"),pch=c(19,19),#fill=c("black","red"),
			col=c("black","red"),lty=c(1,1),lwd=c(1,1),bty="l",box.col="white")
	}

	
	plot_fits_survey_biomass <- function(dat,dat_name,xlab,ylab, mcmcdat=NULL) {
		data_df <- NULL
		for (i in 1: length(dat_name)) {
			type	<- substring(dat_name[i],first=nchar(dat_name[i]),last=nchar(dat_name[i]))
			obs 	<- dat[[dat_name[i]]]$obs
			pred 	<- dat[[dat_name[i]]]$fits
			cv		<- dat[[dat_name[i]]]$error.value 	# CV: Use only for abundance data (otherwise absolute N for proportions)					
			if(!is.null(obs)) data_df <- rbind(data_df,cbind(create_dataframe(pred=pred, obs=obs, cv=cv),type))
		}
		# se1 is CV, sigma2 = log(CV^2+1) for lognormal
		data_df$sigma2 <- log(data_df$cv^2 + 1)		# Variance
		data_df$sigma <- sqrt(data_df$sigma2)		# Standard deviation
		# Calculate mu
		data_df$mu	<- log(data_df$obs)-data_df$sigma2/2
		# data_df$Median	<- exp(data_df$mu)		
		## Approximate 95% CI 
		data_df$obs_upp	<- round(exp(data_df$mu + 1.96*data_df$sigma),0)
		data_df$obs_low	<- round(exp(data_df$mu - 1.96*data_df$sigma),0)
		data_df$obs_low[data_df$obs_low<0]	<- 0
	
		## Calculate Weights
		casalpop	<- extract.csl.file(paste(casal_path,runN,"_",pop_csl,sep=""))
		# lenatage	<- as.numeric(casalpop$size_at_age$all_1982)
		ages		<- casalpop$min_age$value:casalpop$max_age$value
		vb 			<- casalpop$size_at_age
		lw			<- as.numeric(c(casalpop$size_weight$a,casalpop$size_weight$b))
		lenatage	<- as.numeric(vb$Linf)*(1-exp(-as.numeric(vb$k)*((ages+0.5)-as.numeric(vb$t0))))		
			
		data_df$len[data_df$type=="A"]	<- lenatage[data_df$bin[data_df$type=="A"]+0.5]
		data_df$len[data_df$type=="L"]	<- data_df$bin[data_df$type=="L"]
		data_df$wt	<- (lw[1] * data_df$len^lw[2])/1000 # divide by 1000 to convert to tonnes

		## Convert to Weights
		data_df$w_obs		<- data_df$wt * data_df$obs		
		data_df$w_pred    	<- data_df$wt * data_df$pred # / survq
		data_df$w_obs_upp 	<- data_df$wt * data_df$obs_upp
		data_df$w_obs_low 	<- data_df$wt * data_df$obs_low

		Bobs	<- aggregate(data_df$w_obs,by=list(Year=data_df$year),sum)
		Bpred	<- aggregate(data_df$w_pred,by=list(Year=data_df$year),sum)
		Bobsupp	<- aggregate(data_df$w_obs_upp,by=list(Year=data_df$year),sum)
		Bobslow	<- aggregate(data_df$w_obs_low,by=list(Year=data_df$year),sum)

		B	<- cbind(Year=as.numeric(as.character(Bobs$Year)),Bobs=Bobs$x,Bpred=Bpred$x,Bobsupp=Bobsupp$x,Bobslow=Bobslow$x)
		B 	<- B[order(as.numeric(as.character(B[,"Year"]))),]

		## If there is any MCMC data: calculate CI
		if(!is.null(mcmcdat) & length(mcmcdat)>0) {	
			nn			<- unpaste(dimnames(mcmcdat)[[2]],sep="\\[")
			nn[[2]]		<- as.numeric(sub(pattern="]",replacement="",x=nn[[2]]))	# Remove ']', Year
			nn[[3]]		<- as.numeric(sub(pattern="]",replacement="",x=nn[[3]]))	# Remove ']', Ages (or lengths)
			# Multiply age numbers by weight
			mcmc_len 	<- as.numeric(vb$Linf)*(1-exp(-as.numeric(vb$k)*(nn[[3]]-as.numeric(vb$t0))))	# Treat all as ages
			len 		<- unique(data_df$bin[data_df$type == "L"])					# Correct length for survey length data
			mcmc_len[data_df$type == "L"]	<- len[nn[[3]]][data_df$type == "L"]
			mcmc_wt		<- (lw[1] * mcmc_len^lw[2])/1000			
			mcmc_wtN	<- sweep(mcmcdat,2,mcmc_wt,"*")
			bb			<- apply(mcmc_wtN,1,function(x) aggregate(x,by=list(Year=nn[[2]]),sum))	# Sum up by MCMC sample and year 
			bb1			<- NULL
			for (i in 1:length(bb)) bb1	<- cbind(bb1,bb[[i]]$x)	# All results together
			rownames(bb1)	<- bb[[1]]$Year
			ci			<- apply(bb1, 1, quantile, c(0.025,0.5,0.975))	# Calculate quantiles
			B 			<- cbind(B,t(ci))
			#B[,"Bpred"]		<- -100000			# Make MPD fit negative so it will not be visible in the plots
		} else {
			B <- cbind(B,'2.5%'=-100000,'50%'=-100000,'97.5%'=-100000)	# Make it negative so it will not be visible in the plots
		}
	
		## Plot
		maxY 	<- max(c(B[,"Bpred"],B[,"Bobsupp"],B[,"97.5%"]))	
		par(mfrow=c(1,1)); par(mar=c(4,4,2,1))
		plot(y=B[,"Bobs"], x=B[,"Year"], ylim=c(0,max(maxY)), ylab="Survey biomass (t)", xlab="Year", lwd=1,
				type="b",col="black", bg="black", pch=19, font.lab=2,xaxt="n")
		for (ii in 1:length(B[,"Year"])) segments(B[ii,"Year"],B[ii,"Bobslow"],B[,"Year"][ii],B[ii,"Bobsupp"],col="black",lwd=1)
		axis(side=1,at=min(B[,"Year"]):max(B[,"Year"]),labels=min(B[,"Year"]):max(B[,"Year"]))
		# MCMC fits
		polygon(x=c(B[,"Year"],rev(B[,"Year"])),y=c(B[,"2.5%"],rev(B[,"97.5%"])),border="red",col="#FF000080",lty="solid",lwd=1)
		lines(y=B[,"Bpred"], x=B[,"Year"],col="red",bg="red",pch=19,type="b",lwd=1)	
		legend(x=max(B[,"Year"])-3,y=0.95*maxY,legend=c("Observed","MPD fits"),pch=c(19,19),
			col=c("black","red"),lty=c(1,1),lwd=c(1,1),bty="l",box.col="white")
	}


	plot_fits_in_trellis <- function(dat, dat_name, xlab, ylab, mcmcdat=NULL) {
		## Master file to plot MPD data fits in trellis 
		data_df <- NULL
		for (i in 1: length(dat_name)) {
			obs 	<- dat[[dat_name[i]]]$obs
			pred 	<- dat[[dat_name[i]]]$fits
			cv		<- dat[[dat_name[i]]]$error.value 	# CV: Use only for abundance data (otherwise absolute N for proportions)					
			if(!is.null(obs)) data_df <- rbind(data_df,create_dataframe(pred=pred, obs=obs, cv=cv))
		}
		# se1 is CV, sigma2 = log(CV^2+1) for lognormal
		data_df$sigma2 <- log(data_df$cv ^ 2 +1)		# Variance
		data_df$sigma <- sqrt(data_df$sigma2)			# Standard deviation
		# Calculate mu
		data_df$mu	<- log(data_df$obs)-data_df$sigma2/2
		# data_df$Median	<- exp(data_df$mu)		
		## Approximate 95% CI 
		data_df$obs_upp	<- round(exp(data_df$mu + 1.96*data_df$sigma),0)
		data_df$obs_low	<- round(exp(data_df$mu - 1.96*data_df$sigma),0)
		data_df$obs_low[data_df$obs_low<0]	<- 0

		# If there is any MCMC data: calculate CI
		if(!is.null(mcmcdat)) {	
			ci		<- apply(mcmcdat, 2, quantile, c(0.025,0.5,0.975))	
			nn		<- unpaste(dimnames(ci)[[2]],sep="\\[")
			nn[[2]]	<- as.numeric(sub(pattern="]",replacement="",x=nn[[2]]))	# Remove ']'
			nn[[3]]	<- as.numeric(sub(pattern="]",replacement="",x=nn[[3]]))	# Remove ']'
			ci1		<- as.data.frame(cbind(t(ci),Year=nn[[2]],Age=nn[[3]]))
			ci1		<- ci1[order(ci1$Age),]
			data_df <- cbind(data_df,ci1)
			#data_df$pred	<- -100000			# Make MPD fit negative so it will not be visible in the plots
		} else {
			data_df <- cbind(data_df,'2.5%'=-100000,'50%'=-100000,'97.5%'=-100000)	# Make it negative so it will not be visible in the plots
		}
		
		bin <- sort(unique(data_df$bin))
		years <- sort(unique(data_df$year))
		if(ylab=="Proportion") maxY 	<- max(c(data_df$obs,data_df$pred,data_df$'97.5%'))
		if(ylab=="Numbers")    maxY 	<- max(c(data_df$obs_upp,data_df$pred,data_df$'97.5%'))
		print(xyplot(data_df$obs + data_df$pred ~ data_df$bin |  data_df$year, 
					xlab = xlab, ylab = ylab, ylim=c(0,maxY), main=dat_name[1],
					type=c("l","l"),lty=c(1,1), col=c("black","red"),lwd=c(1,2),
					key=list(lines = list(col=c("black","red"),lty=c(1,1)),lwd=c(1,2), background="white",
						text=list(lab=c("Observed","MPD fits")), columns=1, cex=0.7, x=0.73, y=0.95),
					panel=function(...) { 
						panel.xyplot(...)		
						if(ylab=="Numbers") for(i in 1:length(bin))  
							llines(x=c(bin[i],bin[i]),y=c(data_df$obs_low[data_df$year==years[panel.number()]][i],
														  data_df$obs_upp[data_df$year==years[panel.number()]][i]), col="black", lty="solid")
						if(ylab=="Proportion")
							ltext(x=0.9*max(data_df$bin),y=0.2*maxY,labels=data_df$'cv'[data_df$year==years[panel.number()]][1],cex=0.8,font=2)	# ESS
						# MCMC fits
						lpolygon(x=c(bin,rev(bin)),y=c(data_df$'2.5%'[data_df$year==years[panel.number()]],rev(data_df$'97.5%'[data_df$year==years[panel.number()]])),
								  border="red", col="#FF000080", lty="solid",lwd=1)
						llines(x=bin,y=data_df$pred[data_df$year==years[panel.number()]], col="red", lty="solid",lwd=1)
						llines(x=c(5,5), y=c(0,10000000), col="grey", lty="dotted")
						llines(x=c(10,10), y=c(0,10000000), col="grey", lty="dotted")
						llines(x=c(500,500), y=c(0,10000000), col="grey", lty="dotted")
						llines(x=c(1000,1000), y=c(0,10000000), col="grey", lty="dotted")
						llines(x=c(1000,1000), y=c(0,10000000), col="grey", lty="dotted")
					}))  					
	}
		
		
	create_dataframe <- function(pred,obs,cv) {
		## Creates a data frame (suitable format for plotting in trellis)
		dims 	<- dim(obs)
		row.names <- rownames(obs)
		bin.names <- colnames(obs)
		obs <- as.vector(unlist(obs))
		substring(bin.names,first=1) <- " "
		bin <- as.integer(bin.names)
		bin <- bin + (bin[2]-bin[1])/2
		bin <- rep(bin,each=dims[1])
		year <- factor(rep(row.names,times=length(bin.names)))
		pred <- as.vector(unlist(pred))
		cv <- as.vector(unlist(cv))
		data_df <- data.frame(year,bin,obs,pred,cv)
	}

	

	#### Plots mean ages to all items in dat_name
	boxplot_fits_age <- function(dat,dat_name,ylab,mcmcdat=NULL) {
		allobs <- NULL; allpred <- NULL
		for (i in 1: length(dat_name)) {			
			obs 	<- dat[[dat_name[i]]]$obs
			pred 	<- dat[[dat_name[i]]]$fits
			yr		<- dat[[dat_name[i]]]$year
			if(!is.null(obs)) {
				bin.names <-  names(obs)
				substring(bin.names,first=1) <- " "
				bin <- as.integer(bin.names)		
				
				obs 	<- obs/rowSums(obs)*100
				pred 	<- pred/rowSums(pred)*100
				
				for (y in 1: length(yr)) {
					yrobs <- NULL; yrpred <- NULL
					for (ii in 1:length(bin)) {
						yrobs  <- c(yrobs, rep(bin[ii],round(obs[y,ii])))
						yrpred <- c(yrpred,rep(bin[ii],round(pred[y,ii])))
					}
					allobs  <- rbind(allobs,cbind(Fishery=rep(dat_name[i],length(yrobs)),Year=rep(yr[y],length(yrobs)), Obs=yrobs))
					allpred <- rbind(allpred,cbind(Fishery=rep(dat_name[i],length(yrobs)),Year=rep(yr[y],length(yrpred)), Pred=yrpred))
				}
			}
		}
		allobs  <- as.data.frame(allobs)
		allpred <- as.data.frame(allpred)
		
		
		# If there is any MCMC data: calculate CI
		if(!is.null(mcmcdat) & length(mcmcdat>0)) {	
			nn		<- unpaste(dimnames(mcmcdat)[[2]],sep="\\[")
			nn[[2]]	<- as.numeric(sub(pattern="]",replacement="",x=nn[[2]]))	# Remove ']'
			nn[[3]]	<- as.numeric(sub(pattern="]",replacement="",x=nn[[3]]))	# Remove ']'
			# Calculate cumulative sum and median position (code is a bit messy...)
			bb		<- apply(mcmcdat,1,function(x) aggregate(x,by=list(Fishery=nn[[1]],Year=nn[[2]]),cumsum))	# Sum up by MCMC sample and year 
			bb1		<- NULL
			for (i in 1:length(bb)) {
				z	<- bb[[i]][,"x"]
				zrel<- lapply(z, function(x) {x/max(x)})  
				bb1	<- cbind(bb1,sapply(zrel, function(y) {which(y>=0.50)[1]})) # 'Locations of median', not actual ages
			}
			ci		<- apply(bb1, 1, quantile, c(0.025,0.5,0.975))	# Calculate quantiles
			ci		<- cbind(Year=bb[[1]]$Year,t(ci))
			rownames(ci)	<- bb[[1]]$Fishery
			# Convert 'location' to actual ages
			ffish	<- unique(nn[[1]])
			for (f in ffish) {
				fages	<- nn[[3]][nn[[1]] %in% f]
				# print(paste(f, min(fages),max(fages)))
				ci[, "2.5%"][rownames(ci) %in% f] <- seq(min(fages),max(fages))[ci[, "2.5%"][rownames(ci) %in% f]]
				ci[,  "50%"][rownames(ci) %in% f] <- seq(min(fages),max(fages))[ci[,  "50%"][rownames(ci) %in% f]]
				ci[,"97.5%"][rownames(ci) %in% f] <- seq(min(fages),max(fages))[ci[,"97.5%"][rownames(ci) %in% f]]
			}
		}

		## Plot
		par(windows(width=8,height=8))
		par(mfrow=c(2,ceiling(length(dat_name)/2)))
		
		for (i in 1: length(dat_name)) {
			allobsF		<- allobs[allobs$Fishery==dat_name[i],]
			allpredF	<- allpred[allpred$Fishery==dat_name[i],]
			if(!is.null(mcmcdat)) ci1	<- ci[rownames(ci) %in% paste(dat_name[i],".fits",sep=""),]
			boxplot(as.numeric(as.character(allobsF$Obs)) ~ as.numeric(as.character(allobsF$Year)),main=dat_name[i],xlab="Year",
				ylab=ylab, font.lab=2,ylim=c(1,25))			
			median_age <- aggregate(as.numeric(as.character(allpredF$Pred)),by=list(Year=as.numeric(as.character(allpredF$Year))),median)
			if(is.null(mcmcdat))  lines(median_age$x,type="l",col="red",lwd=2)
			if(!is.null(mcmcdat)) polygon(x=c(1:nrow(ci1),nrow(ci1):1),y=c(ci1[,"2.5%"],rev(ci1[,"97.5%"])), border="red", col="#FF000080", lty="solid",lwd=1)

		}
	}


	#### Plot residuals from model age/length data fits
	plot_residuals <- function(dat, dat_name, ylab, plotdims=c(8,8)) {
		# Plot Pearson's residuals for length/age composition data
		allobs 	<- NULL 
		for (i in 1: length(dat_name)) {			
			obs <- dat[[dat_name[i]]]$pearson.resids
			yr	<- dat[[dat_name[i]]]$year
			if(!is.null(obs)) {
				bin.names <-  names(obs)
				substring(bin.names,first=1) <- " "
				bin <- as.integer(bin.names)					
				for (y in 1: length(yr)) {
					yrobs  <- as.numeric(obs[y,])
					allobs  <- rbind(allobs,cbind(Fishery=rep(dat_name[i],length(yrobs)),Year=rep(yr[y],length(yrobs)),Bin=bin,Obs=yrobs))
				}
			}
		}
		allobs  		<- as.data.frame(allobs)
		allobs$Year		<- as.numeric(as.character(allobs$Year))
		allobs$Bin		<- as.numeric(as.character(allobs$Bin))
		allobs$Obs		<- as.numeric(as.character(allobs$Obs))		
		allobs$pos		<- ifelse(allobs$Obs >= 0, TRUE, FALSE)
		## Pool fisheries
		# allobs$Fishery[allobs$Fishery %in% "Catch_LL1Aearly"] <- "Catch_LL1A"
		# allobs$Fishery[allobs$Fishery %in% "Catch_LL2Aearly"] <- "Catch_LL2A"
		dat_name 	<- as.vector(unique(allobs$Fishery))
		
		# Plot
		par(windows(width=plotdims[1],height=plotdims[2]))
		par(mfrow=c(1,1))
		if(length(dat_name)>1) par(mfrow=c(1,2))
		par(mar=c(4,4,2,1))		# inner figure margins
		for (i in 1:length(dat_name)) {
			d1	<- allobs[allobs$Fishery==dat_name[i],]
			xlim 	<- c(min(d1$Year),max(d1$Year))
			ylim	<- c(min(d1$Bin),max(d1$Bin))
			plot(d1$Bin ~ d1$Year, cex=d1$Obs,xlab="Year",xaxt="n",ylab=ylab,type="n",xlim=xlim,ylim=ylim,font.lab=2)
			symbols(d1$Year,d1$Bin,circles=abs(d1$Obs),inches=0.1,add=T,fg=ifelse(d1$pos,"blue","red"), bg=ifelse(d1$pos,"blue","red"))
			mtext(dat_name[i], side=3, line=0.5, font=2)
			axis(1,at=min(d1$Year):max(d1$Year),label=min(d1$Year):max(d1$Year),tick = TRUE)
			#savePlot(filename = paste(casal_path,runN,"_Plot Residuals ",dat_name[i],".png",sep=""),type = "png")
		}
	}


	
	#### Residual plots (code by Timothy Earle 2015-08)
	## Each age by year
	plot_diagnostics_age_by_year <- function(dat, dat_name, trans = rep(identity, length(dat_name)), new_device=FALSE) 	{
		names(trans) <- dat_name
		for(i in dat_name) { 		# All fits not to tags	  
			var_data <- lapply(dat[[i]],as.matrix)
			ages	<- colnames(var_data$obs)
			substring(ages,first=1) <- " "
			ages <- as.integer(ages)	
			for (q in 1:ncol(var_data$obs)) {	    
				lims <- range(c(var_data$obs[,q],var_data$fits[,q]),na.rm=TRUE)
				par(mfrow=c(3,2))
				# Plot 1: Observations by year
				plot(var_data$year, var_data$obs[,q], pch=19, cex=0.8, ylim=lims, xlab="Year", ylab="Observed", main="Observations by year")
				lines(var_data$year, var_data$fits[,q], col='red')
				# Plot 2: Observed versus Fitted
				plot(var_data$obs[,q],var_data$fits[,q], cex=0.8, xlim=lims, ylim=lims, xlab="Observed", ylab="Fitted", main="Observed versus fitted")
				# Plot 3: Residuals
				plot(var_data$year, var_data$pearson.resids[,q], pch=19, cex=0.8, xlab="Year", ylab="Residuals", main="Pearson residuals by year")
				abline(h=0,col="blue",lty=2)
				# Plot 4: QQNorm
				if(!is.na(trans[[i]])) {
					## browser()
					f <- get(trans[i])
					y <- f(var_data$resids[,q])
				    qqnorm(y=y, main="Normal Q-Q Plot")
				    qqline(y=y, col='blue')
				} else {
					plot.new()
				}
				# Plot 4: ACF
				if(all(diff(var_data$obs[,q])==0)) { 	# Avoid error when perfect residuals
					plot.new()
				} else {
					acf(var_data$obs[,q], main="ACF")
				}
				var_title <- paste0("\n",i)
				if (ncol(var_data$obs)>1) var_title <- paste(var_title, ": Age",ages[q])
				title(var_title,outer=TRUE)
				if(new_device) dev.new()
			}
		}
	}
	
	## Each year by age
	plot_diagnostics_year_by_age <- function(dat, dat_name, trans = rep(identity, length(dat_name)), new_device=FALSE) 	{
		names(trans) <- dat_name
		for(i in dat_name) { 		# All fits not to tags	  
			var_data <- lapply(dat[[i]],as.matrix)
			for (q in 1:nrow(var_data$obs)) {	    
				lims <- range(c(var_data$obs[q,],var_data$fits[q,]),na.rm=TRUE)
				ages	<- colnames(var_data$obs)
				substring(ages,first=1) <- " "
				ages <- as.integer(ages)		
				par(mfrow=c(3,2))
				# Plot 1: Observations by age
				plot(ages, var_data$obs[q,], pch=19, cex=0.8, ylim=lims, xlab="Ages", ylab="Observed", main="Observations by age")
				lines(ages, var_data$fits[q,], col='red')
				# Plot 2: Observed versus Fitted
				plot(var_data$obs[q,],var_data$fits[q,], cex=0.8, xlim=lims, ylim=lims, xlab="Observed", ylab="Fitted", main="Observed versus fitted")
				lines(x=c(0,100000000), y=c(0,100000000),col="blue",lty=2)
				# Plot 3: Residuals
				plot(x=ages, y=var_data$pearson.resids[q,], pch=19, cex=0.8, xlab="Ages", ylab="Residuals", main="Pearson residuals by age")
				abline(h=0,col="blue",lty=2)
				# Plot 4: QQNorm
				if(!is.na(trans[[i]])) {
					## browser()
					f <- get(trans[i])
					y <- f(var_data$resids[q,])
				    qqnorm(y=y, main="Normal Q-Q Plot")
				    qqline(y=y, col='blue')
				} else {
					plot.new()
				}
				# Plot 5: ACF
				if(all(diff(var_data$obs[q,])==0)) { 	# Avoid error when perfect residuals
					plot.new()
				} else {
					acf(var_data$obs[q,], main="ACF")
				}
				var_title <- paste0("\n",i)
				if (ncol(var_data$obs)>1) var_title <- paste(var_title, ": Year", var_data$year[q])
				title(var_title,outer=TRUE)
				if(new_device) dev.new()
			}
		}
	}
	
		
	
	
	###### TAG DATA FITS ##################  
	#### Plot tag-recapture fits (by release year, recapture year and length bin)
	plot_fits_tags <- function (dat, plotdims=c(12,12)) {
		tagfits <- NULL
		tagnames <- names(dat)[which(substr(names(dat),1,4)=="Tags")]
		for (i in 1:length(tagnames)) {
			tagy	<- as.integer(substr(tagnames[i],5,8))
			for (ii in 1:length(dat[[tagnames[i]]])) {
				recapy		<- as.integer(names(dat[[tagnames[i]]][ii]))
				dd			<- dat[[tagnames[i]]][[ii]]
				scanned 	<- dd["scanned",]
				recap 		<- dd["recaptured",]
				expected	<- dd["expected_prop",]
				LogLik		<- dd["neg_log_likelihood",]
				bin.names 	<- as.integer(names(scanned))
				tagfits		<- rbind(tagfits,data.frame(tagy,recapy,lengthbin=bin.names,scanned,recap,expected,LogLik))
			}		
		}
		Ntagy	<- min(tagfits$tagy):max(tagfits$tagy)
		Nrecapy	<- min(tagfits$recapy):max(tagfits$recapy)
		
		par(windows(width=plotdims[1],height=plotdims[2]))
		par(mfrow=c(length(Ntagy),length(Nrecapy)))
		par(oma=c(2,6,4,1))		# outer margins
		par(mar=c(0,0,0,0))		# inner figure margins

		ymax	<- max(tagfits$recap,tagfits$expected*tagfits$scanned)
		for (y in 1:length(Ntagy)) {
			count <- 1
			for (x in 1:length(Nrecapy)) {
				d1 <- tagfits[tagfits$tagy==Ntagy[y] & tagfits$recapy==Nrecapy[x],]
				if(nrow(d1)>0) { 
					plot(d1$recap ~ d1$lengthbin,type=c("b"),pch=c(21),lty=c(1),lwd=c(1),col=c("black"),bg=c("black"),
							xlab='Length (mm)',ylab='Recaptures',font.lab=2,ylim=c(0,ymax),xaxt="n",yaxt="n") 
					lines((d1$expected*d1$scanned) ~ d1$lengthbin,type=c("l"),lty=c(1),lwd=c(2),col=c("red"),bg=c("red"))
					abline(v=c(500,1000),lty="dotted",col="grey")
					if(count==1) axis(1, at=d1$lengthbin, labels=d1$lengthbin, tick = TRUE)
					if(count >1) axis(1, at=d1$lengthbin, labels=FALSE, tick = TRUE)
					if(count==1) axis(2, at=NULL, tick = TRUE)					# Tick marks are recycles with at=NULL
					if(count >1) axis(2, at=NULL, labels=FALSE, tick = TRUE)
					count <- count + 1
				} else {		# Blank plot
					plot(1,1, type="n", axes=F, ylim=c(0,1), xlim=c(0,1),xlab="", ylab="")		# or frame()
				}
				if(x==1) mtext(Ntagy[y], side=2, line=3, font=2)
				if(y==1) mtext(Nrecapy[x], side=3, line=0.5, font=2)
			}
		}
		mtext("Releases",   side=2, line=1, font=2, outer=TRUE)
		mtext("Recaptures", side=3, line=2, font=2, outer=TRUE)
	}
	# Multi-area: 
	plot_fits_tags_area <- function (dat, plotdims=c(12,12),casal_path,runN) {
		tagfits 	<- NULL
		tagnames 	<- names(dat)[which(substr(names(dat),1,4)=="Tags")]
		tagnames1 	<- strsplit(tagnames,"_")
		for (i in 1:length(tagnames)) {
			tagy	<- as.integer(substr(tagnames[i],5,8))
			for (ii in 1:length(dat[[tagnames[i]]])) {
				recapy		<- as.integer(names(dat[[tagnames[i]]][ii]))
				dd			<- dat[[tagnames[i]]][[ii]]
				scanned 	<- dd["scanned",]
				recap 		<- dd["recaptured",]
				expected	<- dd["expected_prop",]
				LogLik		<- dd["neg_log_likelihood",]
				bin.names 	<- as.integer(names(scanned))
				relarea		<- tagnames1[[i]][2]
				recarea		<- tagnames1[[i]][3]
				tagfits		<- rbind(tagfits,data.frame(tagy,recapy,recap,expected,scanned,lengthbin=bin.names,LogLik,relarea,recarea))
			}		
		}
		rela	<- unique(tagfits[,"relarea"])
		reca	<- unique(tagfits[,"recarea"])
		Ntagy	<- min(tagfits$tagy):max(tagfits$tagy)
		Nrecapy	<- min(tagfits$recapy):max(tagfits$recapy)
		ymax	<- max(tagfits$recap,tagfits$expected*tagfits$scanned)
		
		for (ll in 1:length(rela)){		# One plot for each combination of release and recapture area 
			for (cc in 1:length(reca)){
				tagfits1	<- tagfits[tagfits[,"relarea"] %in% rela[ll] & tagfits[,"recarea"] %in% reca[cc],]
				
				par(windows(width=plotdims[1],height=plotdims[2]))
				par(mfrow=c(length(Ntagy),length(Nrecapy)))
				par(oma=c(4,6,6,1))		# outer margins
				par(mar=c(0,0,0,0))		# inner figure margins
				for (y in 1:length(Ntagy)) {
					count <- 1
					for (x in 1:length(Nrecapy)) {
						d1 <- tagfits1[tagfits1$tagy==Ntagy[y] & tagfits1$recapy==Nrecapy[x],]
						if(nrow(d1)>0) { 
							plot(d1$recap ~ d1$lengthbin,type=c("b"),pch=c(21),lty=c(1),lwd=c(1),col=c("black"),bg=c("black"),
									xlab='Length (mm)',ylab='Recaptures',font.lab=2,ylim=c(0,ymax),xaxt="n",yaxt="n") 
							lines((d1$expected*d1$scanned) ~ d1$lengthbin,type=c("l"),lty=c(1),lwd=c(2),col=c("red"),bg=c("red"))
							abline(v=c(500,1000),lty="dotted",col="grey")
							if(count==1) axis(1, at=d1$lengthbin, labels=d1$lengthbin, tick = TRUE)
							if(count >1) axis(1, at=d1$lengthbin, labels=FALSE, tick = TRUE)
							if(count==1) axis(2, at=NULL, tick = TRUE)					# Tick marks are recycles with at=NULL
							if(count >1) axis(2, at=NULL, labels=FALSE, tick = TRUE)
							count <- count + 1
						} else {		# Blank plot
							plot(1,1, type="n", axes=F, ylim=c(0,1), xlim=c(0,1),xlab="", ylab="")		# or frame()
						}
						if(x==1) mtext(Ntagy[y], side=2, line=3, font=2)
						if(y==1) mtext(Nrecapy[x], side=3, line=0.5, font=2)
					}
				}
				mtext("Releases",   side=2, line=5, font=2, outer=TRUE)
				mtext("Recaptures", side=3, line=2, font=2, outer=TRUE)
				mtext("Length(mm)", side=1, line=3, font=2, outer=TRUE)
				mtext(paste(rela[ll]," to ",reca[cc],sep=""), side=3, outer=TRUE, line=4, ,cex=1.5,font=2)
				savePlot(filename = paste(casal_path,runN,"_Plot Fits Tag ByLength ",rela[ll],"to",reca[cc],".png",sep=""),type = "png")
			}
		}
	}


	#### Plot tag-recapture Numbers (by release and recapture year)
	plot_fits_tagN <- function (dat,plotdims=c(8,8)) {
		tagfits <- NULL
		tagnames <- names(dat)[which(substr(names(dat),1,4)=="Tags")]
		for (i in 1:length(tagnames)) {
			tagy	<- as.integer(substr(tagnames[i],5,8))
			for (ii in 1:length(dat[[tagnames[i]]])) {
				recapy		<- as.integer(names(dat[[tagnames[i]]][ii]))
				recap 		<- sum(dat[[tagnames[i]]][[ii]]["recaptured",])
				expected	<- sum(dat[[tagnames[i]]][[ii]]["expected_prop",]*dat[[tagnames[i]]][[ii]]["scanned",])
				tagfits		<- rbind(tagfits,data.frame(tagy,recapy,recap,expected))
			}		
		}
		Ntagy	<- min(tagfits$tagy):max(tagfits$tagy)
		Nrecapy	<- min(tagfits$recapy):max(tagfits$recapy)
		
		par(windows(width=plotdims[1],height=plotdims[2]))
		par(mfcol=c(ceiling(length(Ntagy)/2),2))
		par(oma=c(2,1,1,1))		# outer margins
		par(mar=c(2,2.5,3,0))		# inner figure margins
		ymax	<- max(tagfits$recap,tagfits$expected)
		for (y in 1:length(Ntagy)) {
			d1 <- tagfits[tagfits$tagy==Ntagy[y],]
			if(nrow(d1)>0) { 
				plot(d1$recap ~ d1$recapy,type=c("b"),pch=c("o"),lty=c(1),lwd=c(1),col=c("black"),bg=c("black"),
						xlab='',ylab='',ylim=c(0,ymax),xlim=c(min(Nrecapy),max(Nrecapy)),xaxt="n",yaxt="n") 
				lines((d1$expected) ~ d1$recapy,type=c("b"),pch=c("e"),lty=c(1),lwd=c(2),col=c("red"),bg=c("red"))
				#abline(v=c(500,1000),lty="dotted",col="grey")
				axis(1, at=Nrecapy, labels=Nrecapy, tick = TRUE)
				axis(2, at=NULL, tick = TRUE)					# Tick marks are recycles with at=NULL
				if(y<=ceiling(length(Ntagy)/2)) mtext("Numbers", side=2, line=2, font=2)
			} else {		# Blank plot
				plot(1,1, type="n", axes=F, ylim=c(0,1), xlim=c(0,1),xlab="", ylab="")		# or frame()
			}
			mtext(paste("Release year: ",Ntagy[y],sep=""), side=3, line=0.5, font=2)
		}
	}
	# Multi-area: 
	plot_fits_tagN_area <- function (dat,plotdims=c(8,8),casal_path,runN) {
		tagfits 	<- NULL
		tagnames 	<- names(dat)[which(substr(names(dat),1,4)=="Tags")]
		tagnames1 	<- strsplit(tagnames,"_")
		for (i in 1:length(tagnames)) {
			tagy	<- as.integer(substr(tagnames[i],5,8))
			for (ii in 1:length(dat[[tagnames[i]]])) {
				recapy		<- as.integer(names(dat[[tagnames[i]]][ii]))
				recap 		<- sum(dat[[tagnames[i]]][[ii]]["recaptured",])
				expected	<- sum(dat[[tagnames[i]]][[ii]]["expected_prop",]*dat[[tagnames[i]]][[ii]]["scanned",])
				relarea		<- tagnames1[[i]][2]
				recarea		<- tagnames1[[i]][3]
				tagfits		<- rbind(tagfits,data.frame(tagy,recapy,recap,expected,relarea,recarea))
			}		
		}
		rela	<- unique(tagfits[,"relarea"])
		reca	<- unique(tagfits[,"recarea"])
		Ntagy	<- min(tagfits$tagy):max(tagfits$tagy)
		Nrecapy	<- min(tagfits$recapy):max(tagfits$recapy)
		ymax	<- max(tagfits$recap,tagfits$expected)
		for (ll in 1:length(rela)){		# One plot for each combination of release and recapture area 
			for (cc in 1:length(reca)){
				tagfits1	<- tagfits[tagfits[,"relarea"] %in% rela[ll] & tagfits[,"recarea"] %in% reca[cc],]
				
				par(windows(width=plotdims[1],height=plotdims[2]))
				par(mfcol=c(ceiling(length(Ntagy)/2),2))
				par(oma=c(2,1,1,1))		# outer margins
				par(mar=c(2,2.5,3,0))		# inner figure margins
				for (y in 1:length(Ntagy)) {
					d1 <- tagfits1[tagfits1$tagy==Ntagy[y],]
					if(nrow(d1)>0) { 
						plot(d1$recap ~ d1$recapy,type=c("b"),pch=c("o"),lty=c(1),lwd=c(1),col=c("black"),bg=c("black"),
								xlab='',ylab='',ylim=c(0,ymax),xlim=c(min(Nrecapy),max(Nrecapy)),xaxt="n",yaxt="n") 
						lines((d1$expected) ~ d1$recapy,type=c("b"),pch=c("e"),lty=c(1),lwd=c(2),col=c("red"),bg=c("red"))
						#abline(v=c(500,1000),lty="dotted",col="grey")
						axis(1, at=Nrecapy, labels=Nrecapy, tick = TRUE)
						axis(2, at=NULL, tick = TRUE)					# Tick marks are recycles with at=NULL
						if(y<=ceiling(length(Ntagy)/2)) mtext("Numbers", side=2, line=2, font=2)
					} else {		# Blank plot
						plot(1,1, type="n", axes=F, ylim=c(0,1), xlim=c(0,1),xlab="", ylab="")		# or frame()
					}
					mtext(paste("Release year: ",Ntagy[y],sep=""), side=3, line=0.5, font=2)
				}
				mtext(paste(rela[ll]," to ",reca[cc],sep=""), side=3, outer=T, line=-1, ,cex=1.5,font=2)
				savePlot(filename = paste(casal_path,runN,"_Plot Fits Tag Numbers ",rela[ll],"to",reca[cc],".png",sep=""),type = "png")
			}
		}
	}
	



	#### Plot expected - Observed tag-recapture Numbers in one graph by recapture year
	plot_fits_tags_Exp_Obs <- function (dat,plotdims=c(9,6)) {
		tagfits <- NULL
		tagnames <- names(dat)[which(substr(names(dat),1,4)=="Tags")]
		for (i in 1:length(tagnames)) {
			tagy	<- as.integer(substr(tagnames[i],5,8))
			for (ii in 1:length(dat[[tagnames[i]]])) {
				recapy		<- as.integer(names(dat[[tagnames[i]]][ii]))
				recap 		<- sum(dat[[tagnames[i]]][[ii]]["recaptured",])
				expected	<- sum(dat[[tagnames[i]]][[ii]]["expected_prop",]*dat[[tagnames[i]]][[ii]]["scanned",])
				tagfits		<- rbind(tagfits,data.frame(tagy,recapy,recap,expected))
			}		
		}
		tagfits$diff		<- tagfits$expected-tagfits$recap 
		Ntagy	<- min(tagfits$tagy):max(tagfits$tagy)
		Nrecapy	<- min(tagfits$recapy):max(tagfits$recapy)
		ally	<- seq(min(Ntagy,Nrecapy),max(Ntagy,Nrecapy))
		# Colours to use from cyan to blue
		ramp <- colorRamp(c("cyan", "blue"))
		cols <- rgb(ramp(seq(0, 1, length = length(Ntagy))), max = 255)
		par(windows(width=plotdims[1],height=plotdims[2]))
		par(mfrow=c(1,1))
		par(mar=c(5,5,3,2))		# inner figure margins
		ylims	<- c(min(tagfits$diff),max(tagfits$diff))
		plot(y=rep(0,length(ally)),x=ally,type="l",lty=1,lwd=1,col=c("grey"),ylim=ylims,font.lab=2,xaxt="n",
				xlab='Year of recapture',ylab='Expected minus observed recaptures') 
		axis(1,at=ally,labels=ally)
		for (y in Ntagy) {
			d1 <- tagfits[tagfits$tagy==y,]
			if(nrow(d1)>0) { 
				lines(d1$diff ~ d1$recapy,type="b",lty=1,lwd=1,pch=as.character(y-floor(y/10)*10),cex=1,col=cols[y+1-Ntagy[1]])
			}
		}
	}
	# Multi-area: 
	plot_fits_tags_Exp_Obs_area <- function (dat,plotdims=c(9,6),casal_path,runN) {
		tagfits 	<- NULL
		tagnames 	<- names(dat)[which(substr(names(dat),1,4)=="Tags")]
		tagnames1 	<- strsplit(tagnames,"_")
		for (i in 1:length(tagnames)) {
			tagy	<- as.integer(substr(tagnames[i],5,8))
			for (ii in 1:length(dat[[tagnames[i]]])) {
				recapy		<- as.integer(names(dat[[tagnames[i]]][ii]))
				recap 		<- sum(dat[[tagnames[i]]][[ii]]["recaptured",])
				expected	<- sum(dat[[tagnames[i]]][[ii]]["expected_prop",]*dat[[tagnames[i]]][[ii]]["scanned",])
				relarea		<- tagnames1[[i]][2]
				recarea		<- tagnames1[[i]][3]
				tagfits		<- rbind(tagfits,data.frame(tagy,recapy,recap,expected,relarea,recarea))
			}		
		}
		tagfits$diff		<- tagfits$expected-tagfits$recap 
		rela	<- unique(tagfits[,"relarea"])
		reca	<- unique(tagfits[,"recarea"])
		Ntagy	<- min(tagfits$tagy):max(tagfits$tagy)
		Nrecapy	<- min(tagfits$recapy):max(tagfits$recapy)
		ally	<- seq(min(Ntagy,Nrecapy),max(Ntagy,Nrecapy))
		# Colours to use from cyan to blue
		ramp <- colorRamp(c("cyan", "blue"))
		cols <- rgb(ramp(seq(0, 1, length = length(Ntagy))), max = 255)
		ylims	<- c(min(tagfits$diff),max(tagfits$diff))
		
		par(windows(width=plotdims[1],height=plotdims[2]))
		par(mfcol=c(ceiling(length(Ntagy)/2),2))
		par(oma=c(0,0,1,1))		# outer margins
		par(mar=c(5,5,3,2))		# inner figure margins

		for (ll in 1:length(rela)){		# One plot for each combination of release and recapture area 
			for (cc in 1:length(reca)){
				tagfits1	<- tagfits[tagfits[,"relarea"] %in% rela[ll] & tagfits[,"recarea"] %in% reca[cc],]
				plot(y=rep(0,length(ally)),x=ally,type="l",lty=1,lwd=1,col=c("grey"),ylim=ylims,font.lab=2,xaxt="n",
						xlab='Year of recapture',ylab='Expected minus observed recaptures') 
				axis(1,at=ally,labels=ally)
				mtext(paste(rela[ll]," to ",reca[cc],sep=""), side=3, outer=F, line=1, ,cex=1,font=2)
				for (y in Ntagy) {
					d1 <- tagfits1[tagfits1$tagy==y,]
					if(nrow(d1)>0) { 
						lines(d1$diff ~ d1$recapy,type="b",lty=1,lwd=1,pch=as.character(y-floor(y/10)*10),cex=1,col=cols[y+1-Ntagy[1]])
					}
				}
			}
		}
		savePlot(filename = paste(casal_path,runN,"_Plot Fits Tag Exp-Obs",".png",sep=""),type = "png")
	}
	
	#### Plot median length of tag-recapture fits (by release year, recapture year and length bin)
	plot_fits_tags_MedianL <- function (dat, plotdims=c(12,12)) {
		tagfits <- NULL
		tagnames <- names(dat)[which(substr(names(dat),1,4)=="Tags")]
		for (i in 1:length(tagnames)) {
			tagy	<- as.integer(substr(tagnames[i],5,8))
			for (ii in 1:length(dat[[tagnames[i]]])) {
				recapy		<- as.integer(names(dat[[tagnames[i]]][ii]))
				dd			<- dat[[tagnames[i]]][[ii]]
				bin.names 	<- as.integer(names(dd["scanned",]))
				recap 		<- mean(rep(bin.names,times=dd["recaptured",]))						# Mean observed
				expected	<- mean(rep(bin.names,dd["expected_prop",] * dd["scanned",] * 100))	# Mean expected
				tagfits		<- rbind(tagfits,data.frame(tagy,recapy,recap,expected))
			}		
		}		
		Ntagy	<- min(tagfits$tagy):max(tagfits$tagy)
		Nrecapy	<- min(tagfits$recapy):max(tagfits$recapy)
		par(windows(width=plotdims[1],height=plotdims[2]))
		par(mfcol=c(ceiling(length(Ntagy)/2),2))
		par(oma=c(2,1,1,1))			# outer margins
		par(mar=c(2,2.5,3,0))		# inner figure margins
		ymax	<- max(tagfits$recap,tagfits$expected)*1.2
		for (y in 1:length(Ntagy)) {
			d1 <- tagfits[tagfits$tagy==Ntagy[y],]
			if(nrow(d1)>0) { 
				plot(d1$recap ~ d1$recapy,type=c("b"),pch=c("o"),lty=c(1),lwd=c(1),col=c("black"),bg=c("black"),
						xlab='',ylab='',ylim=c(0,ymax),xlim=c(min(Nrecapy),max(Nrecapy)),xaxt="n",yaxt="n") 
				lines((d1$expected) ~ d1$recapy,type=c("b"),pch=c("e"),lty=c(1),lwd=c(2),col=c("red"),bg=c("red"))
				#abline(v=c(500,1000),lty="dotted",col="grey")
				axis(1, at=Nrecapy, labels=Nrecapy, tick = TRUE)
				axis(2, at=NULL, tick = TRUE)					# Tick marks are recycles with at=NULL
				if(y<=ceiling(length(Ntagy)/2)) mtext("Mean length (mm)", side=2, line=2, font=2)
			} else {		# Blank plot
				plot(1,1, type="n", axes=F, ylim=c(0,1), xlim=c(0,1),xlab="", ylab="")		# or frame()
			}
			mtext(paste("Release year: ",Ntagy[y],sep=""), side=3, line=0.5, font=2)
		}
	}
	# Multi-area: 
	plot_fits_tags_MedianL_area <- function (dat, plotdims=c(12,12),casal_path,runN) {
		tagfits 	<- NULL
		tagnames 	<- names(dat)[which(substr(names(dat),1,4)=="Tags")]
		tagnames1 	<- strsplit(tagnames,"_")
		for (i in 1:length(tagnames)) {
			tagy	<- as.integer(substr(tagnames[i],5,8))
			for (ii in 1:length(dat[[tagnames[i]]])) {
				recapy		<- as.integer(names(dat[[tagnames[i]]][ii]))
				dd			<- dat[[tagnames[i]]][[ii]]
				bin.names 	<- as.integer(names(dd["scanned",]))
				recap 		<- mean(rep(bin.names,times=dd["recaptured",]))						# Mean observed
				expected	<- mean(rep(bin.names,dd["expected_prop",] * dd["scanned",] * 100))	# Mean expected
				relarea		<- tagnames1[[i]][2]
				recarea		<- tagnames1[[i]][3]
				tagfits		<- rbind(tagfits,data.frame(tagy,recapy,recap,expected,relarea,recarea))
			}		
		}		
		rela	<- unique(tagfits[,"relarea"])
		reca	<- unique(tagfits[,"recarea"])
		Ntagy	<- min(tagfits$tagy):max(tagfits$tagy)
		Nrecapy	<- min(tagfits$recapy):max(tagfits$recapy)
		ymax	<- max(tagfits$recap,tagfits$expected,na.rm=T)*1.2
		for (ll in 1:length(rela)){		# One plot for each combination of release and recapture area 
			for (cc in 1:length(reca)){
				tagfits1	<- tagfits[tagfits[,"relarea"] %in% rela[ll] & tagfits[,"recarea"] %in% reca[cc],]
				
				par(windows(width=plotdims[1],height=plotdims[2]))
				par(mfcol=c(ceiling(length(Ntagy)/2),2))
				par(oma=c(2,1,1,1))		# outer margins
				par(mar=c(2,2.5,3,0))		# inner figure margins

				for (y in 1:length(Ntagy)) {
					d1 <- tagfits1[tagfits1$tagy==Ntagy[y],]
					if(nrow(d1)>0) { 
					plot(d1$recap ~ d1$recapy,type=c("b"),pch=c("o"),lty=c(1),lwd=c(1),col=c("black"),bg=c("black"),
							xlab='',ylab='',ylim=c(0,ymax),xlim=c(min(Nrecapy),max(Nrecapy)),xaxt="n",yaxt="n") 
					lines((d1$expected) ~ d1$recapy,type=c("b"),pch=c("e"),lty=c(1),lwd=c(2),col=c("red"),bg=c("red"))
					#abline(v=c(500,1000),lty="dotted",col="grey")
					axis(1, at=Nrecapy, labels=Nrecapy, tick = TRUE)
					axis(2, at=NULL, tick = TRUE)					# Tick marks are recycles with at=NULL
					if(y<=ceiling(length(Ntagy)/2)) mtext("Mean length (mm)", side=2, line=2, font=2)
				} else {		# Blank plot
					plot(1,1, type="n", axes=F, ylim=c(0,1), xlim=c(0,1),xlab="", ylab="")		# or frame()
				}
				mtext(paste("Release year: ",Ntagy[y],sep=""), side=3, line=0.5, font=2)
				}
				mtext(paste(rela[ll]," to ",reca[cc],sep=""), side=3, outer=T, line=-1, ,cex=1.5,font=2)
				savePlot(filename = paste(casal_path,runN,"_Plot Fits Tag MedianLength ",rela[ll],"to",reca[cc],".png",sep=""),type = "png")
			}
		}
	}
	


	############################################################
	#### Extra plots
	## Population numbers-at-age
	PopN_at_age <- function(dat) {
		natage		<- dat$Numbers_at_age
		age_year 	<- as.data.frame(cbind(Year=as.numeric(rep(names(natage),each=length(natage[[1]]))),
									  Age=(rep(1:length(natage[[1]]),length(names(natage)))),Numbers=unlist(natage)))
		# Median age
		cumsums <- aggregate(age_year$Numbers,by=list(Year=age_year$Year), cumsum)
		cumsums  <- as.matrix(cumsums)
		cumsums[,2:ncol(cumsums)]	<- cumsums[,2:ncol(cumsums)]/cumsums[,ncol(cumsums)]
		median_age <- unique(age_year$Age)[apply(cumsums[,2:ncol(cumsums)]<=0.5,1,sum)+1] # First Age with cumsum > 0.5
		
		par(mfrow=c(1,1))
		symbols(age_year[,"Year"], age_year[,"Age"], circles=age_year[,"Numbers"],inches=0.1,fg="black",bg="black",xaxt="n",
				xlab="Year",ylab="Age (year)",xlim=c(1980,max(age_year$Year)),ylim=c(0,35),font.lab=2,main="Population numbers at age")		
		axis(1,seq(1980,max(age_year$Year),1))
		lines(unique(age_year[,"Year"]),median_age, col="red",lwd=2)
	}
	
	#### Stock-recruitment relationship
	stock_rec <- function(dat) {
		SSB 		<- dat$SSBs
		YCS 		<- dat$YCS
		true_YCS 	<- dat$true_YCS
		# Select years: All years where YCS was estimated
		YCS_years	<- YCS$year[YCS$YCS!=1]
		# SSB_years	<- SSB$year[SSB$SSB != SSB$SSB[2]]
		# est_years	<- SSB_years[SSB_years %in% YCS_years]	 # Years where SSB is not SSB0
		# Select data
		est_SSB		<- SSB$SSB[SSB$year %in% YCS_years]
		est_YCS		<- true_YCS$true_YCS[true_YCS$year %in% YCS_years] * dat$R0/1000000
		# Plot
		plot(est_YCS ~ est_SSB, type="p", pch=21,col="red",bg="red",tcl=-0.3,
				lwd=1,lty="solid",ylab="Recruitment (Mio)",xlab="SSB (t)",ylim=c(0,max(est_YCS)),xlim=c(0,max(est_SSB)),font.lab=2)
		text(est_SSB, est_YCS, labels=YCS_years, cex= 0.7,pos=3,col="grey")
	}
	
	
	cryptic_Biomass <- function(dat=res_quant, dat_name, ages, WatA, mat, firsty, plotdims=c(10,6)) {
		# dat 		<- res_quant
		# dat_name	<- c("Sel_Trawl1","Sel_LL1","Sel_LL2","Sel_Trap")
		# firsty 	<- 2000
		# ages 		<- cpara$age[1]:cpara$age[2]
		# WatA: Weight at age

		## Numbers-at-age to be used
		N <- dat$'Numbers_at_age'
		yrs	<- which(as.numeric(names(N)) >= firsty)
	
		## Select selectivities from dat_names
		ogives	<- dat$'Ogive parameter values'
		og_names <- which(names(ogives) %in% c(paste("selectivity[",dat_name,"].all",sep="")))
		ogives 	<- ogives[og_names]
		# Combined LL selectivity
		if(sum(c("Sel_LL1","Sel_LL2") %in% dat_name)>0) {
			ll <- rbind(ogives[["selectivity[Sel_LL1].all"]],ogives[["selectivity[Sel_LL2].all"]])
			sel_LL	<- apply(ll,2,max)
		}
		# Combined selectivity of all gear types
		all 	<- matrix(unlist(ogives),nrow=length(ogives),ncol=length(ages),byrow=TRUE)
		sel_all	<- apply(all,2,max)
		
		par(windows(width=9,height=6))
		par(mar=c(4,4,1,1))	
		plot(y=sel_all,x=ages, col="black",type="l",lwd=2,xlab="Ages (y)",ylab="Selectivity",ylim=c(0,1),font.lab=2)
		lines(y=ogives[["selectivity[Sel_LL1].all"]],x=ages, col="grey",lwd=1,lty=1)
		lines(y=ogives[["selectivity[Sel_LL2].all"]],x=ages, col="grey",lwd=1,lty=1)
		lines(y=ogives[["selectivity[Sel_Trap].all"]],x=ages, col="grey",lwd=1,lty=2)
		lines(y=sel_LL,x=ages, col="blue",lwd=2)
		lines(y=mat,x=ages, col="red",lwd=1,lty=2)
		legend(y=0.95, x=max(ages)-6,col=c("black","blue","grey","grey","red"),lty=c(1,1,1,2,2),lwd=c(2,2,1,1,1),# bg="white",
					bty="l",box.col="white", legend=c("All","All longline","Longline","Trap","Maturity"),cex=0.9)
		# savePlot(filename = paste(casal_path,runN,"_Plot Selectivities Combined",".png",sep=""),type = "png")
		
		## Plot Numbers-at-age by year
		par(windows(width=plotdims[1],height=plotdims[2]))
		par(mfrow=c(ceiling(length(yrs)/4),4)); par(mar=c(4,4,2,1))
		par(oma=c(4.5,4.5,1,1))		# outer margins
		par(mar=c(0,0,0,0))		# inner figure margins
		
		ymax <- max(unlist(N[yrs])* mat/1000)
		for (y in 1:length(yrs)) {
			Ny	<- N[[yrs[y]]] * mat/1000
			plot(y=Ny, x=ages, type="l", col="black", tcl=0.2, lwd=1, lty="solid",ylim=c(0,ymax),
					xaxt="n",yaxt="n",ylab="",xlab="")
			text(x=1,y=ymax*0.9,labels=names(N)[yrs[y]], pos=4,font=2,cex=1.4)
			if(ceiling((y+3)/4) == floor((y+3)/4)) { axis(2, at=NULL, tick = TRUE) # Tick marks are recycles with at=NULL
				} else {axis(2, at=NULL, labels=FALSE, tick=TRUE, tcl=0.2)}
			if(y >= length(yrs)-3) { axis(1, at=NULL, tick=TRUE) # Tick marks are recycles with at=NULL
				} else {axis(1, at=NULL, labels=FALSE, tick=TRUE, tcl=0.2)}
			
			lines(y=sel_LL* Ny,x=ages, col="blue",lwd=2)
			lines(y=sel_all*Ny,x=ages, col="red",lwd=2)
			lines(y=Ny,x=ages, col="black",lwd=2)
			mtext("Age (y)", side=1, line=3, font=2, outer=TRUE)
			mtext("Numbers (in '000s)", side=2, line=3, font=2, outer=TRUE)
			if(y==1) legend(y=ymax, x=max(ages)-12,col=c("black","red","blue"),lty=c(1,1,1),lwd=c(2,2,2),# bg="white",
					bty="l",box.col="white", legend=c("Est","Obs All","Obs LL"),cex=1)
		}
	
		## Calculations: 
		# Proportions of total SSB observed 
		props	<- matrix(0,nrow=length(yrs),ncol=5,dimnames=list(names(N)[yrs],c("Total","VulnAll","VulnLL","PropAll","PropLL")))
		for (y in 1:length(yrs)) {
			Ny	<- N[[yrs[y]]] * mat * WatA
			Ny_all	<- Ny* sel_all
			Ny_LL	<- Ny* sel_LL
			props[y,]	<- c(round(c(sum(Ny), sum(Ny_all), sum(Ny_LL)),0), round(c(sum(Ny_all)/sum(Ny), sum(Ny_LL)/sum(Ny)),3))
		}		
		return(props)
	}

	
	calc_SSB_by_region	<- function(casal_path,pop_csl,res_quant) {
		## Calculate SSB by regions from numbers at age, maturity-at-age, length-at-age and length-weight
		## and match with CASAL calculated total SSB
		pop.file	<- extract.csl.file(paste(casal_path,pop_csl,sep=""))
		ages 		<- seq(as.numeric(pop.file[["min_age"]]$value),as.numeric(pop.file[["max_age"]]$value))
		mat 		<- pop.file[["maturity_props"]]$all
		wt			<- as.numeric(pop.file[["size_weight"]]$a)
		wt[2]		<- as.numeric(pop.file[["size_weight"]]$b)
		vb 			<- as.numeric(pop.file[["size_at_age"]]$Linf)
		vb[2] 		<- as.numeric(pop.file[["size_at_age"]]$k)
		vb[3] 		<- as.numeric(pop.file[["size_at_age"]]$t0)
		vb[4] 		<- as.numeric(pop.file[["size_at_age"]]$cv)
		
		if(mat[1] == "allvalues_bounded") {	# c("allvalues_bounded", 11,17, 0.0,0.1667,0.3333,0.5000,0.6667,0.8333,1.0000) 
			maturity	<- c(rep(0,as.numeric(mat[2])-1),mat[4:length(mat)])
			maturity	<- c(as.numeric(maturity),rep(1,length(ages)-length(mat)))
		}
		if(mat[1] == "allvalues") maturity <- as.numeric(mat[2:length(mat)])
		if(mat[1] == "logistic")  maturity <- 1/(1+19^((as.numeric(mat[2])-ages)/as.numeric(mat[3])))
		growth	<- vb[1] * (1-exp(-vb[2]*(ages-vb[3])))
		# Weight in tonnes by length-at-age (without cv in length-at-age)
		weight	<- wt[1] * growth^wt[2]	
		## But: Need to acount for cv in VB
		dens_val	<- seq(0.5,1.5,0.1)
		vb_density	<- dnorm(dens_val,mean=1,sd=vb[4])
		vb_density	<- vb_density/sum(vb_density)
		
		a1 <- res_quant[substring(names(res_quant),1,14)== "Numbers_at_age"]
		a_reg	<- NULL
		a_res	<- NULL
		for (rr in 1:length(a1)) {
			a2 <- a1[[rr]]
			a_region <- substring(names(a1[rr]),nchar(names(a1[rr]))-1,nchar(names(a1[rr])))
			a3 <- a2[[1]]	# Select ages of first year
			# Numbers of mature fish at age
			a_mat	<- a3 * maturity
			# Weight of mature fish at age (without cv in VB)
			# a_weight <- a_mat * weight
			a_weight <- 0
			for (aa in 1:length(ages)) {
				aN 		<- a_mat[aa] * vb_density 	# Numbers of fish in each bin
				lb 		<- dens_val*growth[aa]		# Length of fish in each bin
				wgt		<- aN * (wt[1] * lb^wt[2])	# Weight of all fish in each bin
				atot	<- sum(wgt)					# Total weight of all fish of age aa
				a_weight <- a_weight + atot
			}
			a_reg <- rbind(a_reg,a_region)
			a_res <- rbind(a_res,round(sum(a_weight),0))
		}
		a_reg <- rbind(a_reg,"Total"); 
		sum_a_res <- sum(a_res); a_res <- rbind(a_res,sum_a_res)
		print(paste("Total SSB: ", round(res_quant$B0,0)," tonnes",sep=""))
		print(cbind(a_reg,a_res))		
		# Approximately the same as SSB, thus probably similar calculations as within CASAL

		## Correct a_res
		print("Regional SSB corrected to exactly match CASAL B0:")
		a_res_corr <- round(a_res/sum_a_res*res_quant$B0,0)
		print(cbind(a_reg,a_res_corr))
	}
	
	
	
	
	
	############################################################################
	####   2. CASAL.proj_boxplot
	
	CASAL_box	<- function(Plot_mat, YearRange, Ylimits, ablines_h, XYlabels, current_y="") {
		## Create box plots given CASAL output matrix with rows = trials & cols = years	
		## Input parameters:
		#    Plot_mat      Input data matrix with rows = trials & cols = years
		#    YearRange     First and last years of time seriesc - (yr_first,yr_last) 
		#    XYlabels      Axes labels - c(Xlabel,Ylabel) 
		#    Xlimits       X-axis limits - c(Xmin,Xmax)
		#    Ylimits       Y-axis limits - c(Ymin,Ymax) 
		#    ablines_h     Vector for horizontal lines - c(line1,line2...)
		#    ablines_v     Vector for vertical lines - c(line1,line2...)

		# par(windows(width=9,height=6))
		colnames(Plot_mat) <- seq(YearRange[1],YearRange[2])
		nyr	<- seq(YearRange[1],YearRange[2])- YearRange[1]+1
		# Set the quantiles
		kk 	<- apply(Plot_mat, 2, quantile, c(0.0,0.1,0.5,0.9,1.0))	
		print(kk[3,])
		yrs	<- seq(YearRange[1],YearRange[2])
		plot(kk[3,] ~ yrs, type="l", pch=21,col="black",bg="black",tcl=-0.3,
			lwd=1,lty="solid",xlab=XYlabels[1],ylab=XYlabels[2],ylim=Ylimits,font.lab=2)
		polygon(x=c(yrs,rev(yrs)),y=c(kk[1,],rev(kk[5,])),col="lightgrey",border="lightgrey")
		polygon(x=c(yrs,rev(yrs)),y=c(kk[2,],rev(kk[4,])),col="darkgrey",border="darkgrey")
		points(kk[3,] ~ yrs, type="l", pch=21,col="black",bg="black",lwd=2,lty="solid")
		abline(h=ablines_h, col="gray60", lty=2)
		abline(v=(current_y), col="blue", lty=2)
	}

					CASAL_box_old	<- function(Plot_mat, YearRange, Ylimits, ablines_h, XYlabels, current_y="") {
						## Boxplot
						## Create box plots given CASAL output matrix with rows = trials & cols = years	
						## Input parameters:
						#    Plot_mat      Input data matrix with rows = trials & cols = years
						#    YearRange     First and last years of time seriesc - (yr_first,yr_last) 
						#    XYlabels      Axes labels - c(Xlabel,Ylabel) 
						#    Xlimits       X-axis limits - c(Xmin,Xmax)
						#    Ylimits       Y-axis limits - c(Ymin,Ymax) 
						#    ablines_h     Vector for horizontal lines - c(line1,line2...)
						#    ablines_v     Vector for vertical lines - c(line1,line2...)

						par(windows(width=9,height=6))
						colnames(Plot_mat) <- seq(YearRange[1],YearRange[2])
						nyr			<- seq(YearRange[1],YearRange[2])- YearRange[1]+1
						kk 			<- boxplot(Plot_mat, plot=FALSE) 
						# Set the quantiles
						kk$stats 	<- apply(Plot_mat, 2, quantile, c(0.1,0.25,0.5,0.75,0.9))	
						# Plot all outliers outside the 10-90 Percentile
						outlierU	<- sweep(Plot_mat,2,apply(Plot_mat, 2, quantile, c(0.1)),"<")	
						outlierL	<- sweep(Plot_mat,2,apply(Plot_mat, 2, quantile, c(0.9)),">")
						kk$out		<- c(as.vector(unlist(Plot_mat[outlierL])),as.vector(unlist(Plot_mat[outlierU])))
						yrs			<- matrix(rep(nyr,each=nrow(Plot_mat)),nrow=nrow(Plot_mat), ncol=ncol(Plot_mat))
						kk$group	<- c(as.vector(unlist(yrs[outlierL])),as.vector(unlist(yrs[outlierU])))
						bxp(kk, ylim=Ylimits, xlab=XYlabels[1], ylab=XYlabels[2],outcol="darkgrey",outcex=0.5)
						bxp(kk, add=TRUE, outline=FALSE)		# Make sure that the ends of the whiskers are on top
						abline(h=ablines_h, col="gray60", lty=2)
						abline(v=(current_y-YearRange[1]+1), col="blue", lty=2)
						if(FALSE) {
							MatBox		<- NULL
							YearBox		<- c(YearRange[1]:YearRange[2])
							for (i in 1:ncol(Plot_mat)) MatBox <- rbind(MatBox,cbind(YearBox[i],Plot_mat[,i]))
							par(windows(width=9,height=6))
							# boxplot(MatBox[,2]~MatBox[,1], ylim=Ylimits, xlab=XYlabels[1], ylab=XYlabels[2])
							abline(h=ablines_h, col="gray60", lty=2)
						}
					}


	
	
	
	
	
	
					################################################################################################  
					################################################################################################  
					################################################################################################  
					################################################################################################  
					
					#### Code from Rob Scott (Cefas) ####
					plot_fits_tags_RobScott <- function (dat, dat_name, finalyr, xlab, ylab) {
						
						#### From Rob Scott, 2014 Punta Arenas
						firstl		<- gregexpr("Tags",dat_name,fixed=TRUE)[[1]][1]			# First letter for 'Tags'
						firstly		<- gregexpr("x",dat_name,fixed=TRUE)[[1]][1]			# First letter for year
						tagnames	<- names(dat)[substr(names(dat),firstl,firstl+3) == "Tags"]
						tagevents	<- as.numeric(substr(tagnames,firstly,firstly+3))
						allyears	<- seq(tagevents[1],finalyr)
						tag.years 	<- seq(tagevents[1],finalyr-1)
						
						length.bins  <- colnames(dat[[tagnames[1]]][[1]])
						tag.events   <- c(rep(tag.years,seq(length(tag.years),1)))
						recap.events <- NULL
						for (i in 1:length(allyears)) recap.events <-  c(recap.events,allyears[-1:-i])
						
						# tag.events   <- c(rep(2005,7),rep(2006,6),rep(2007,5),rep(2008,4),rep(2009,3), rep(2010,2), 2011)
						# tag.years    <- seq(2005,2012)
						# recap.events <- c(tag.years[-1], tag.years[-1:-2], tag.years[-1:-3], tag.years[-1:-4], tag.years[-1:-5], tag.years[-1:-6], tag.years[-1:-7])
					
						# length.bins  <- seq(400,1700,by=100)
						# tag.events   <- c(2012,2012,2013)
						# recap.events <- c(2013,2014,2014)

						scanned <- recaptured <- expected <- LogLik <- NULL
						for(aa in 1:length(tag.events)){
							scanned    <- c(scanned,    el(el(dat, paste(as.character(tag.events[aa]),"Tags",sep="")),as.character(recap.events[aa]))[1,])
							recaptured <- c(recaptured, el(el(dat, paste(as.character(tag.events[aa]),"Tags",sep="")),as.character(recap.events[aa]))[2,])
							expected   <- c(expected,   el(el(dat, paste(as.character(tag.events[aa]),"Tags",sep="")),as.character(recap.events[aa]))[3,])
							LogLik     <- c(LogLik,     el(el(dat, paste(as.character(tag.events[aa]),"Tags",sep="")),as.character(recap.events[aa]))[4,])
						}

						tag.fits.df<- data.frame(tag.events  =rep(tag.events, each=length(length.bins)), 
												 recap.events=rep(recap.events,each=length(length.bins)), 
												 length.bins =length.bins,
												 scanned     =scanned,
												 recaptured  =recaptured,
												 expected    =expected,
												 LogLik      =LogLik)
						tag.fits.df <- tag.fits.df[tag.fits.df$recap.events>2006,]
						tag.fits.df$recaptured[tag.fits.df$recaptured==0] <- NA
						xyplot(recaptured+(expected*scanned)~length.bins|as.character(recap.events)*as.character(tag.events), data=tag.fits.df, 
							   type=c("p","l"), pch=c(3,NA), lty=c(1,1), xlab='Length (cm)', ylab='Recaptures')	
					}
						

					#### Tag fits - not used 
					plot_fits_in_trellis_tags <- function(dat, dat_name, xlab, ylab) {
						data_df <- NULL
						for (tt in 1:length(names(dat[[dat_name]]))) {		# Data are by year
							yyear <- names(dat[[dat_name]])[tt]
							dat_y 	<- dat[[dat_name]][[yyear]]
							data_df <- rbind(data_df, create_dataframe_tags(dat_y=dat_y, yyear=yyear))
							data_df[data_df==0]	<- NA
						}
						
						xlim	<- c(0.9*min(data_df$bin),1.1*max(data_df$bin))
						ylim	<- c(0,max(c(data_df$obs,data_df$pred))*1.05)
						background <- trellis.par.get("background")
						background$col <- "white"
						trellis.par.set("background",background)
						print(xyplot(data_df$obs + data_df$pred ~ data_df$bin |  data_df$year.f, 
									xlab = xlab, ylab = ylab, ylim=ylim, xlim=xlim, 
									type=c("b","b"),lty=c(1,2),	pch=c(21,21), col=c("black","grey"),fill=c("black","grey"),lwd=c(1,1),main=dat_name,
									key=list(lines = list(col=c("black","grey"),lty=c(1,2)),type=c("b","b"),lwd=c(1,1), background="white",
										pch=c(21,21), col=c("black","grey"),fill=c("black","grey"),
										text=list(lab=c("Observed","Fitted")), columns=1, cex=0.7, x=0.73, y=0.95)))  					
					}
					
					create_dataframe_tags <- function(dat_y,yyear) {
						dims 	<- dim(dat_y)
						obs 	<- as.vector(unlist(dat_y["recaptured",]))
						#if(sum(obs) > 0) obs <- obs /sum(obs)
						pred1 	<- as.vector(unlist(dat_y["expected_prop",]*dat_y["scanned",]))
						bin.names <-  names(as.data.frame(dat_y))
						bin		<- as.integer(bin.names)
						bin 	<- bin + (bin[2]-bin[1])/2
						year 	<- factor(rep(yyear,times=length(bin.names)))
						data_df <- data.frame(year,bin,obs,pred)
					}



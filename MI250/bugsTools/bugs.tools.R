#$ 5/20/2008 WRG: Fixed bug in save.model. Added n.chains argument.
## 7/2/2010 WRG: Added na.rm = TRUE to summary(x1) statement

myshell = function(cmd){
# runs DOS command "cmd" without opening a command window
	system(paste(Sys.getenv("COMSPEC"),"/c",cmd), intern = FALSE, wait = TRUE, input = NULL,
       show.output.on.console = FALSE,
       minimized = FALSE, invisible = TRUE)
}

mcmc.history = function(sims.array){

	n1 = dim(sims.array)[1]
	n2 = dim(sims.array)[2]
	n3 = dim(sims.array)[3]
	x = data.frame(value=as.vector(sims.array),
		chain=rep(rep(1:n2,ea=n1),n3),
		parameter=rep(dimnames(sims.array)[[3]],ea=n1*n2))
	print(xyplot(value~rep(1:n1,n2*n3)|parameter,x,groups=chain,
		panel=panel.superpose,type="l",col=c(1,2,3,4),
		layout=c(1,6),scales=list(cex=1,y=list(relation="free")),
		xlab=list(label="sample",cex=1.2),ylab=list(label="value",cex=1.2),
		par.strip.text=list(cex=1),strip = function(...) strip.default(..., style = 1)))
	NULL
}

mcmc.density = function(sims.array,n=50){

# n = number of points at which the density is calculated

	n1 = dim(sims.array)[1]
	n2 = dim(sims.array)[2]
	n3 = dim(sims.array)[3]
	x = data.frame(value=as.vector(sims.array),
		chain=rep(rep(1:n2,ea=n1),n3),
		parameter=rep(dimnames(sims.array)[[3]],ea=n1*n2))
	x = data.frame(parameter=rep(unique(x$parameter),ea=n),
		value=as.vector(sapply(unique(x$parameter),function(x,n,sim.list){
		density(sim.list$value[sim.list$parameter==x],
		n=n,na.rm=T)$x},n=n,sim.list=x)),
		frequency=as.vector(sapply(unique(x$parameter),function(x,n,sim.list){density(
		sim.list$value[sim.list$parameter==x],
		n=n,na.rm=T)$y},n=n,sim.list=x)))

	print(xyplot(frequency~value|parameter,x,scales="free",type="l",col=1,
		layout=c(0,min(16,length(unique(x$parameter)))),
		par.strip.text=list(cex=1),strip = function(...) strip.default(..., style = 1)))
	NULL
}

save.model = function(bugs.output,model.name,parent.dir=".",n.chains=1){

# assumes model script and R script are both named "model.name" with the extensions .txt and .R, respectively

	file.move = function(from, to, overwrite=F){
		file.copy(from, to, overwrite=overwrite)
		file.remove(from)
	}

	dir.create(paste(parent.dir,"/",model.name,sep=""))
	files = eval(parse(text=paste("c(\"codaIndex.txt\",",paste("\"coda",1:n.chains,".txt\",",
		sep="",collapse=""),paste("\"inits",1:n.chains,".txt\",",sep="",
		collapse=""),"\"data.txt\",\"log.odc\",\"log.txt\",\"script.txt\")",sep="")))
	sapply(files,function(x) file.move(x,paste(parent.dir,"/",model.name,sep=""),overwrite=T))
	file.copy(paste(parent.dir,"/",model.name,".txt",sep=""),paste(parent.dir,"/",model.name,sep=""),overwrite=T)
	file.copy(paste(parent.dir,"/",model.name,".R",sep=""),paste(parent.dir,"/",model.name,sep=""),overwrite=T)
	bugs.fit = bugs.output # looks silly but it works and appears to be necessary
	save(bugs.fit,file=paste(parent.dir,"/",model.name,"/",model.name,".fit.Rsave",sep=""))

}

summary.mcmc.list <- function (object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
    ...) 
{
    x <- mcmc.list(object)
    statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
    varstats <- matrix(nrow = nvar(x), ncol = length(statnames), 
        dimnames = list(varnames(x), statnames))
    xtsvar <- matrix(nrow = nchain(x), ncol = nvar(x))
    if (is.matrix(x[[1]])) {
        for (i in 1:nchain(x)) for (j in 1:nvar(x)) xtsvar[i, 
            j] <- coda:::safespec0(x[[i]][, j])
        xlong <- do.call("rbind", x)
    }
    else {
        for (i in 1:nchain(x)) xtsvar[i, ] <- coda:::safespec0(x[[i]])
        xlong <- as.matrix(x)
    }
    xmean <- apply(xlong, 2, mean, na.rm = TRUE)
    xvar <- apply(xlong, 2, var, na.rm = TRUE)
    xtsvar <- apply(xtsvar, 2, mean, na.rm = TRUE)
    varquant <- t(apply(xlong, 2, quantile, quantiles, na.rm = TRUE))
    varstats[, 1] <- xmean
    varstats[, 2] <- sqrt(xvar)
    varstats[, 3] <- sqrt(xvar/(niter(x) * nchain(x)))
    varstats[, 4] <- sqrt(xtsvar/(niter(x) * nchain(x)))
    varquant <- drop(varquant)
    varstats <- drop(varstats)
    out <- list(statistics = varstats, quantiles = varquant, 
        start = start(x), end = end(x), thin = thin(x), nchain = nchain(x))
    class(out) <- "summary.mcmc"
    return(out)
}

parameter.plot.table = function(parameter.array){
# create history, density and Gelman-Rubin-Brooks plots
# return value is a table of summary stats for the parameters

	# create history, density and Gelman-Rubin-Brooks plots
	mcmc.history(parameter.array)
	mcmc.density(parameter.array,n=50)
	if(length(dim(parameter.array))<3){
		n.chains = 1
	}else{
		n.chains = dim(parameter.array)[2]
	}
	x1 = mcmc.list(lapply(1:n.chains,function(i) mcmc(parameter.array[,i,]))) # format required by CODA
	try(gelman.plot(x1, ask = FALSE),TRUE)

	# summary stats on parameters
	psummary = summary(x1, na.rm = TRUE)
	if(is.vector(psummary$statistics)){ # when there is only one monitored parameter
		psummary$statistics = matrix(psummary$statistics,nrow=1,dimnames=list(NULL,names(psummary$statistics)))
		psummary$quantiles = matrix(psummary$quantiles,nrow=1,dimnames=list(NULL,names(psummary$quantiles)))
	}
	ptable = cbind(psummary$statistics,psummary$quantiles)
	neff = try(effectiveSize(x1),TRUE)
	if(exists("neff")){
		ptable = cbind(psummary$statistics,psummary$quantiles,neff)
		dimnames(ptable)[[2]][ncol(ptable)] = "Effective N"
	}
	ptable
}

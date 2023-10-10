# hierarchical dose-response model

modelName <- "dose-response1"

if(.Platform$OS.type == "windows") Sys.setenv(HOME = substr(R.home(), 1, 2))

courseDir <- file.path(Sys.getenv("HOME"), "MI250")
exampleDir <- file.path(courseDir, modelName)
toolsDir <- file.path(courseDir, "bugsTools") # change to match your setup
bugsDir <- "c:/Program Files/WinBUGS14"
if(.Platform$OS.type == "windows") 
   bugsDir <- file.path(Sys.getenv("HOME"), "Program Files/WinBUGS14")
wineBin <- "/opt/local/bin" # directory containing wine and winepath programs.
                                        # only relevant on unix or Mac OS X
setwd(exampleDir)
library(R2WinBUGS)
library(coda)
library(lattice)
library(nlme)
source(file.path(toolsDir, "bugs.tools.R"))
source(file.path(toolsDir, "bgillespie.utilities.R"))

set.seed(10271998) # not required but assures repeatable results

# get data file
xdata = read.table("dose-response1.csv",sep=",",header=T)

# create WinBUGS data set
bugsdata = list(
	narm = nrow(xdata),
	nstudy = length(unique(xdata$study)),
	study = xdata$study,
	dose = xdata$dose,
	n = xdata$n,
	resp = xdata$resp)

# create initial estimates

# method 1: explicitly specify initial estimates as a list of lists (1 per chain)
bugsinit = list(
	list(
		mu.emax = 10,
		mu.e0 = 2,
		log.ed50 = log(6),
		omega.e0 = 0.1,
		omega.emax = 0.1,
		0.15),
	list(
		mu.emax = 10,
		mu.e0 = 2,
		log.ed50 = log(6),
		omega.e0 = 0.1,
		omega.emax = 0.1,
		0.15),
	list(
		mu.emax = 10,
		mu.e0 = 2,
		log.ed50 = log(6),
		omega.e0 = 0.1,
		omega.emax = 0.1,
		0.15))

# method 2: specify initial estimates with a random generation function
fit = gnls(resp~e0+emax*dose/(ed50+dose),xdata,start=list(e0=2,emax=10,ed50=6))
se = sqrt(diag(fit$varBeta))

bugsinit = function() list(
	mu.emax = rnorm(1,log(fit$coefficients["emax"]),4*se["emax"]/fit$coefficients["emax"]),
	mu.e0 = rnorm(1,log(fit$coefficients["e0"]),4*se["emax"]/fit$coefficients["e0"]),
	log.ed50 = rnorm(1,log(fit$coefficients["ed50"]),4*se["ed50"]/fit$coefficients["ed50"]),
	omega.e0 = exp(rnorm(1,log(0.2),0.5)),
	omega.emax = exp(rnorm(1,log(0.2),0.5)),
	sigma = exp(rnorm(1,log(fit$sigma),0.5)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("med.emax","med.e0","ed50","omega.emax","omega.e0","sigma")

## Additional variables to monitor
otherRVs <- c("emax","e0","resp.cond","resp.pred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("deviance", parametersToPlot)

################################################################################################
# run WinBUGS

n.chains = 3
n.iter = 10000
n.burnin = 4000
n.thin = 5
bugs.fit <- bugs(data = bugsdata, inits = bugsinit, parameters.to.save = parameters,
                 model.file = file.path(exampleDir, paste(modelName, ".txt", sep="")),
                 n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                 n.thin = n.thin, clearWD = T,
                 bugs.directory = bugsDir, working.directory = getwd(),
                 useWINE = (.Platform$OS.type == "unix"),
                 WINE = file.path(wineBin,"wine"),
                 newWINE = (.Platform$OS.type == "unix"),
                 WINEPATH = file.path(wineBin,"winepath"))

# save scripts, data and results to a directory

save.model(bugs.fit, modelName)
#load(paste(modelName,"/",modelName,".fit.Rsave",sep=""))

# rename and reformat MCMC results to facilitate later calculations and plots

sims.array = bugs.fit$sims.array
posterior = array(as.vector(sims.array),dim=c(prod(dim(sims.array)[1:2]),dim(sims.array)[3]),
	dimnames=list(NULL,dimnames(sims.array)[[3]]))

################################################################################################
## posterior distributions of parameters

## open graphics device
pdf(file = file.path(exampleDir, paste(modelName,"/",modelName,"Plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

## subset of sims.array containing selected variables
x1 <- sims.array[,,unlist(sapply(c(paste("^",parametersToPlot,"$",sep=""),
	paste("^",parametersToPlot,"\\[",sep="")),grep,x=dimnames(sims.array)[[3]]))]


## create history, density and Gelman-Rubin-Brooks plots, and a table of summary stats
ptable <- parameter.plot.table(x1)
write.csv(signif(ptable,3),paste(modelName,"/",modelName,".summary.csv",sep=""))

################################################################################################
# posterior predictive distributions

# prediction of future observations in the same studies, i.e., posterior predictions
# conditioned on observed data from the same study

pred = posterior[,grep("resp.cond\\[",dimnames(posterior)[[2]])]

x1 = xdata
x1$type =rep("observed",nrow(x1))
x2 = rbind(x1,x1,x1)
x2$resp = as.vector(t(apply(pred,2,quantile,probs=c(0.05,0.5,0.95))))
x2$type = rep(c("5%ile","median","95%ile"),ea=nrow(x1))
x1 = rbind(x1,x2)

xyplot(resp~dose|as.factor(study),x1,groups=type,
	panel=panel.superpose,type=c("l","l","l","p"),lty=c(3,3,1,0),col=c(2,2,4,1),
	pch=c(NA,NA,NA,1),lwd=2,cex=0.5,scales=list(cex=1),
	xlab=list(label="dose",cex=1.2),ylab=list(label="response",
	cex=1.2),par.strip.text=list(cex=1),
	strip = function(...) strip.default(..., style = 1))

# prediction of future observations in a new study, i.e., posterior predictive distributions

pred = posterior[,grep("resp.pred\\[",dimnames(posterior)[[2]])]

x1 = xdata
x1$type =rep("observed",nrow(x1))
x2 = rbind(x1,x1,x1)
x2$resp = as.vector(t(apply(pred,2,quantile,probs=c(0.05,0.5,0.95))))
x2$type = rep(c("5%ile","median","95%ile"),ea=nrow(x1))
x1 = rbind(x1,x2)

xyplot(resp~dose|as.factor(study),x1,groups=type,
	panel=panel.superpose,type=c("l","l","l","p"),lty=c(3,3,1,0),col=c(2,2,4,1),
	pch=c(NA,NA,NA,1),lwd=2,cex=0.5,scales=list(cex=1),
	xlab=list(label="dose",cex=1.2),ylab=list(label="response",
	cex=1.2),par.strip.text=list(cex=1),
	strip = function(...) strip.default(..., style = 1))

dev.off()



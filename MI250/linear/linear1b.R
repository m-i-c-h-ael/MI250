# Simple linear model example

modelName <- "linear1b"
exampleDir <- "linear"

if(.Platform$OS.type == "windows") Sys.setenv(HOME = substr(R.home(), 1, 2))

courseDir <- file.path(Sys.getenv("HOME"), "MI250")
exampleDir <- file.path(courseDir, exampleDir)
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
source(file.path(toolsDir, "bugs.tools.R"))
source(file.path(toolsDir, "bgillespie.utilities.R"))

set.seed(10271998) # not required but assures repeatable results

# create WinBUGS data set
bugsdata = list(
	x = c(1,2,3,4,5,6,7,8,9,10),
	y = c(5.19,6.56,9.19,8.09,7.6,7.08,6.74,9.3,8.98,11.5)
)

# create initial estimates
fit = lm(y~x,bugsdata)
se = summary(fit)$coefficients[,2]
sigma = summary(fit)$sigma
bugsinit = function() list(
	a = rnorm(1,coef(fit)[1],5*se[1]),
	b = rnorm(1,coef(fit)[2],5*se[2]),
	sigma = exp(rnorm(1,log(sigma),log(sigma)))
)

## Specify the variables for which you want history and density plots
parametersToPlot <- c("a","b","sigma")

## Additional variables to monitor
otherRVs <- c("ypred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("deviance", parametersToPlot)

################################################################################################
# run WinBUGS

n.chains = 3
n.iter = 1000
n.burnin = 100
n.thin = 1
bugs.fit <- bugs(data = bugsdata, inits = bugsinit, parameters.to.save = parameters,
                 model.file = file.path(exampleDir, paste(modelName, ".txt", sep="")),
                 n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                 n.thin = n.thin, clearWD = FALSE, debug = FALSE,
                 bugs.directory = bugsDir, working.directory = getwd(),
                 useWINE = (.Platform$OS.type == "unix"),
                 WINE = file.path(wineBin,"wine"),
                 newWINE = (.Platform$OS.type == "unix"),
                 WINEPATH = file.path(wineBin,"winepath"))

## save scripts, data and results to a directory

save.model(bugs.fit, modelName, n.chains = n.chains)
#load(paste(modelName,"/",modelName,".fit.Rsave",sep=""))

## rename and reformat MCMC results to facilitate later calculations and plots

sims.array <- bugs.fit$sims.array
posterior <- array(as.vector(sims.array),dim=c(prod(dim(sims.array)[1:2]),dim(sims.array)[3]),
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

pred = posterior[,grep("ypred\\[",dimnames(posterior)[[2]])]

x1 = data.frame(x=bugsdata$x,y=bugsdata$y)
x1$type =rep("observed",nrow(x1))
x2 = rbind(x1,x1,x1)
x2$y = as.vector(t(apply(pred,2,quantile,probs=c(0.05,0.5,0.95))))
x2$type = rep(c("5%ile","median","95%ile"),ea=nrow(x1))
x1 = rbind(x1,x2)

xyplot(y~x,x1,groups=type,panel=panel.superpose,
	type=c("l","l","l","p"),lty=c(3,3,1,0),col=c(2,2,4,1),
	pch=c(NA,NA,NA,19),cex=1.5,lwd=3,
	scales=list(cex=1),xlab=list(xlab="x",cex=1.2),
	ylab=list(ylab="y",cex=1.2))

dev.off()


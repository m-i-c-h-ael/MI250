modelName <- "poppk" # root names of model file

if(.Platform$OS.type == "windows") Sys.setenv(HOME = substr(R.home(), 1, 2))

courseDir <- file.path(Sys.getenv("HOME"), "MI250")
exampleDir <- file.path(courseDir, modelName)
toolsDir <- file.path(courseDir, "bugsTools") # change to match your setup
bugsDir <- "c:/Program Files/BlackBoxWinBUGS"
if(.Platform$OS.type == "windows") 
   bugsDir <- file.path(Sys.getenv("HOME"), "Program Files/BlackBoxWinBUGS")
wineBin <- "/opt/local/bin" # directory containing wine and winepath programs.
                                     # only relevant on unix or Mac OS X
setwd(exampleDir)
library(R2WinBUGS)
library(coda)
library(lattice)
source(file.path(toolsDir, "bugs.tools.R"))
source(file.path(toolsDir, "bgillespie.utilities.R"))

## Get data files
pkData <- read.csv(file.path(exampleDir, "censoredData.csv"), as.is = TRUE)

## Logical indicating BLQ data
pkData$blq <- (pkData$DV == "BLQ")

lq <- 0.1 # Lower limit of quantitation
big <- 1.0E20 # a big number

## Convert AMT and DV to numeric. This will also convert "." and "BLQ" to NA
pkData$AMT <- as.numeric(pkData$AMT)
pkData$DV <- as.numeric(pkData$DV)

## Make sequential ID
pkData$seqID <- as.numeric(as.factor(pkData$ID))

## Row index for start of each individual's data
nObs <- nrow(pkData)
start <- (1:nObs)[!duplicated(pkData$seqID)]

## Format data for WinBUGS
bugsdata <- list(
                 nObs = nObs,
                 nSubjects = length(unique(pkData$seqID)),
                 start = start,
                 end = c(start[-1]-1,nObs),
                 subject = pkData$seqID,
                 time = pkData$TIME,
                 amt = ifelse(is.na(pkData$AMT), 0, pkData$AMT),
                 rate = rep(0, nObs),
                 ii = rep(0, nObs),
                 evid = pkData$EVID,
                 cmt = ifelse(pkData$EVID > 0, 1, 2),
                 addl = rep(0, nObs),
                 ss = rep(0, nObs),
                 logCobs = log(pkData$DV),
                 loglq = ifelse(pkData$blq, log(lq), big),
                 omegaInvPrior = 5 * diag(rep(0.05,5))
                 )

## Create initial estimates
bugsinit = function() list(
	logCLHat = rnorm(1,log(25),0.5),
	logQHat = rnorm(1,log(10),0.5),
	logV1Hat = rnorm(1,log(50),0.5),
	logV2Hat = rnorm(1,log(75),0.5),
	logDkaHat = rnorm(1,log(0.25),0.5),
	sigma = runif(1,0.05,0.5),
	omegaInv = solve(diag(exp(2 * rnorm(5, log(0.2), 0.5)))))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "DkaHat",
                      "sigma", "omega")

## Additional variables to monitor
otherRVs <- c("logCobsCond", "logCobsPred", "theta")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("deviance", parametersToPlot)

################################################################################################
# run WinBUGS

n.chains <- 3
n.iter <- 10000
n.burnin <- 4000
n.thin <- 5
bugs.fit <- bugs(data = bugsdata, inits = bugsinit, parameters.to.save = parameters,
                 model.file = file.path(exampleDir, paste(modelName, ".txt", sep="")),
                 n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                 n.thin = n.thin, clearWD = F,
                 bugs.directory = bugsDir, working.directory = getwd(),
                 useWINE = (.Platform$OS.type == "unix"),
                 WINE = file.path(wineBin,"wine"),
                 newWINE = (.Platform$OS.type == "unix"),
                 WINEPATH = file.path(wineBin,"winepath"))

## save scripts, data and results to a directory
save.model(bugs.fit, modelName, n.chains = n.chains)
##load(paste(modelName,"/",modelName,".fit.Rsave",sep=""))

## rename and reformat MCMC results to facilitate later calculations and plots
sims.array = bugs.fit$sims.array
posterior = array(as.vector(sims.array),dim=c(prod(dim(sims.array)[1:2]),dim(sims.array)[3]),
  dimnames=list(NULL,dimnames(sims.array)[[3]]))
                                 
################################################################################################
## posterior distributions of parameters

## open graphics device
pdf(file = file.path(exampleDir, paste(modelName,"/",modelName,"Plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

## subset of sims.array containing selected variables
x1 = sims.array[,,unlist(sapply(c(paste("^",parametersToPlot,"$",sep=""),
	paste("^",parametersToPlot,"\\[",sep="")),grep,x=dimnames(sims.array)[[3]]))]

## create history, density and Gelman-Rubin-Brooks plots, and a table of summary stats
ptable = parameter.plot.table(x1)
write.csv(ptable, paste(modelName,"/",modelName,".summary.csv",sep=""))

################################################################################################
# posterior predictive distributions

## prediction of future observations in the same subjects, i.e., posterior predictions
## conditioned on observed data from the same subject

pred <- exp(posterior[,grep("logCobsCond\\[",dimnames(posterior)[[2]])])

## Calculate posterior median, 5th and 95th percentiles for each observation

x1 <- pkData
x2 <- rbind(x1,x1,x1)
x2$DV <- as.vector(t(apply(pred,2,quantile,probs=c(0.05,0.5,0.95))))
x2$type <- rep(c("5%ile","median","95%ile"),ea=nrow(x1))
x1$type <- rep("observed",nrow(x1))
x1 <- rbind(x1,x2)

xyplot(DV ~ TIME | as.factor(ID), x1,
       groups = type,
       panel = panel.superpose, type = c("l","l","l","p"),
       lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
       pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5,
       scales = list(cex = 1),
       xlab = list(label = "time (h)", cex=1.2),
       ylab = list(label = "plasma concentration (mg/L)", cex=1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1),
       main =  list(label = "individual predictions", cex=1.2))

xyplot(log(DV) ~ TIME | as.factor(ID), x1,
       groups = type,
       panel = panel.superpose, type = c("l","l","l","p"),
       lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
       pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5,
       scales = list(cex = 1),
       xlab = list(label = "time (h)", cex=1.2),
       ylab = list(label = "log(plasma concentration (mg/L))", cex=1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1),
       main =  list(label = "individual predictions", cex=1.2))

################################################################################################
# posterior predictive distributions

## prediction of future observations in new subjects

pred <- exp(posterior[,grep("logCobsPred\\[",dimnames(posterior)[[2]])])

## Calculate posterior median, 5th and 95th percentiles for each observation

x1 <- pkData
x2 <- rbind(x1,x1,x1)
x2$DV <- as.vector(t(apply(pred,2,quantile,probs=c(0.05,0.5,0.95))))
x2$type <- rep(c("5%ile","median","95%ile"),ea=nrow(x1))
x1$type <- rep("observed",nrow(x1))
x1 <- rbind(x1,x2)

xyplot(DV ~ TIME | as.factor(ID), x1,
       groups = type,
       panel = panel.superpose, type = c("l","l","l","p"),
       lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
       pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5,
       scales = list(cex = 1),
       xlab = list(label = "time (h)", cex=1.2),
       ylab = list(label = "plasma concentration (mg/L)", cex=1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1),
       main =  list(label = "population predictions", cex=1.2))

xyplot(log(DV) ~ TIME | as.factor(ID), x1,
       groups = type,
       panel = panel.superpose, type = c("l","l","l","p"),
       lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
       pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5,
       scales = list(cex = 1),
       xlab = list(label = "time (h)", cex=1.2),
       ylab = list(label = "log(plasma concentration (mg/L))", cex=1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1),
       main =  list(label = "population predictions", cex=1.2))

################################################################################################
## random effect plots

effects <- "theta"
effectsPost <- posterior[,unlist(sapply(c(paste("^",effects,"$",sep=""),
	paste("^",effects,"\\[",sep="")),grep,x=dimnames(posterior)[[2]]))]

x1 <- pkData[!duplicated(pkData$seqID),]

## rename to improve readability
effects <- c("logCL", "logQ", "logV1", "logV2", "logDka")
effectNames <- colnames(effectsPost)
for(id in x1$seqID){
  effectNames <- gsub(paste("theta\\[", id, ",1\\]", sep = ""),
                 paste("logCL\\[", id, "\\]", sep = ""), effectNames)
  effectNames <- gsub(paste("theta\\[", id, ",2\\]", sep = ""),
                 paste("logQ\\[", id, "\\]", sep = ""), effectNames)
  effectNames <- gsub(paste("theta\\[", id, ",3\\]", sep = ""),
                 paste("logV1\\[", id, "\\]", sep = ""), effectNames)
  effectNames <- gsub(paste("theta\\[", id, ",4\\]", sep = ""),
                 paste("logV2\\[", id, "\\]", sep = ""), effectNames)
  effectNames <- gsub(paste("theta\\[", id, ",5\\]", sep = ""),
                 paste("logDka\\[", id, "\\]", sep = ""), effectNames)
}
colnames(effectsPost) <- effectNames

medianEffects = sapply(effects,
  function(effect, posterior, n){
    eta <- t(posterior[,grep(paste(effect,"\\[",sep=""),dimnames(posterior)[[2]])])
    apply(eta,
          1, median)
  }, posterior = effectsPost)
colnames(medianEffects) <- effects
medianEffects <- cbind(x1, medianEffects)

for(effect in effects){

  print(histogram(~ medianEffects[effect], medianEffects,
                  ylab = list(cex=1.2), xlab=list(label = effect, cex=1.2),
                  scales=list(cex=1.2)))
  print(qqmath(~ medianEffects[effect], medianEffects,
               prepanel = prepanel.qqmathline, distribution = qnorm,
               panel=function(x, distribution){
                 panel.qqmath(x,cex=1.2,pch=16, distribution = distribution)
                 panel.qqmathline(x,distribution=distribution,col=3,lwd=3)
               },
               ylab=list(cex=1.2),
               xlab=list(label = effect, cex=1.2),
               scales=list(cex=1.2)))

}

#########################################################################################
## residuals

## prediction of future observations in the same subjects, i.e., posterior predictions
## conditioned on observed data from the same subject

pred = posterior[,grep("logCobsCond\\[",dimnames(posterior)[[2]])]

cPred <- apply(pred, 2, median, na.rm = T)

resid <- log(pkData$DV) - cPred

xyplot(resid ~ TIME, pkData, panel =
       function(x, y, ...){
         panel.xyplot(jitter(x, 2), y, ...)
         panel.abline(h = 0, col = "red", lty = 2, lwd = 2)
       },
       cex=0.5,scales=list(cex=1),
       xlab=list(label="time (h)",cex=1.2),
       ylab=list(label="log(plasma concentration) residual", cex=1.2),
       par.strip.text=list(cex=1),
       strip = function(...) strip.default(..., style = 1))

xyplot(resid[!is.na(pkData$DV)] ~ cPred[!is.na(pkData$DV)], panel =
       function(x, y, ...){
         panel.xyplot(x, y, ...)
         panel.abline(h = 0, col = "red", lty = 2, lwd = 2)
       },
       cex=0.5,scales=list(cex=1),
       xlab=list(label="log(plasma concentration) predicted",cex=1.2),
       ylab=list(label="log(plasma concentration) residual", cex=1.2),
       par.strip.text=list(cex=1),
       strip = function(...) strip.default(..., style = 1))

histogram(~ resid,
          xlab=list(label="log(plasma concentration) residual",cex=1.2),
          ylab = list(cex=1.2),
          scales=list(cex=1.2))
qqmath(~ resid,
       prepanel = prepanel.qqmathline, distribution = qnorm,
       panel=function(x, distribution){
         panel.qqmath(x,cex=1.2,pch=16, distribution = distribution)
         panel.qqmathline(x,distribution=distribution,col=3,lwd=3)
       },
       ylab=list(label="log(plasma concentration) residual",cex=1.2),
       xlab=list(cex=1.2),
       scales=list(cex=1.2))

dev.off()


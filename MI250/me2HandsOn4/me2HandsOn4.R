modelName = "me2HandsOn4" # root names of model file

if(.Platform$OS.type == "windows") Sys.setenv(HOME = substr(R.home(), 1, 2))

courseDir <- file.path(Sys.getenv("HOME"), "MI250")
exampleDir <- file.path(courseDir, modelName)
toolsDir <- file.path(courseDir, "bugsTools") # change to match your setup
bugsDir <- "c:/Program Files/BlackBoxWinBUGS"
if(.Platform$OS.type == "windows") 
   bugsDir <- file.path(Sys.getenv("HOME"), "Program Files/BlackBoxWinBUGS")
wineBin <- "/opt/local/bin" # directory containing wine and winepath programs.
                                     # only relevant on unix or Mac OS X
memory.limit(2048)

setwd(exampleDir)
library(R2WinBUGS)
##source(file.path(Sys.getenv("HOME"),"bugsTools", "myR2WinBUGS.R"))
library(coda)
library(lattice)
source(file.path(toolsDir, "bugs.tools.R"))
source(file.path(toolsDir, "bgillespie.utilities.R"))
source(paste(modelName, "Sim.R", sep = ""))

set.seed(10271998) # not required but assures repeatable results

lq <- 10 # Lower limit of quantitation
big <- 1.0E20 # a big number

## get the NONMEM-formatted data file
xdata <- read.csv("me2HandsOn4NONMEMData.csv", as.is = TRUE)

## Let's just do a subset for now to illustrate the method
## Drop placebo arm
##xdata <- xdata[xdata$DOSE > 0,]
## Select first half of each treatment arm
subjects <- with(xdata, tapply(ID, list(DOSE), unique))
subjects <- unlist(lapply(subjects, function(x) x[1:(length(x)%/%2)]))
xdata <- xdata[xdata$ID %in% subjects,]

## convert those pesky "."'s to NA's
xdata$COBS <- as.numeric(xdata$COBS)
xdata$NEUT <- as.numeric(xdata$NEUT)

nobs <- nrow(xdata)
## Create unique sequential ID's for all subjects. This allows the NONMEM practice of
## using the same ID's for different subjects
xdata$origID <- xdata$ID
subject <- 1
id <- xdata$ID
xdata$ID[1] <- subject
for(i in 2:nobs){
  if(id[i] != id[i-1]) subject <- subject + 1
  xdata$ID[i] <- subject
}

## Row index for start of each individual's data
nobs <- nrow(xdata)
start <- (1:nobs)[!duplicated(xdata$ID)]

## Parameters for informative priors

circ0Hat <- 5.4
circ0CV <- 0.20
mttHat <- 110
mttCV <- 0.16
gammaHat <- 0.16
gammaCV <- 0.16

# The following refer to priors for IIV parameters expressed as precisions (1/variance)
circ0PrecMean <- 11
circ0PrecSD <- 21
mttPrecMean <- 37
mttPrecSD <- 61

## create WinBUGS data set
bugsdata <- list(
                 nobs = nobs,
                 nsub = length(unique(xdata$ID)),
                 start = start,
                 end = c(start[-1] - 1, nobs),
                 subject = xdata$ID,
                 weight = xdata$WEIGHT,
                 time = xdata$TIME,
                 amt = xdata$AMT,
                 rate = xdata$RATE,
                 ii = xdata$II,
                 evid = xdata$EVID,
                 cmt = xdata$CMT,
                 addl = xdata$ADDL,
                 ss = xdata$SS,
                 logCobs = ifelse(xdata$COBS <= 0, NA, log(xdata$COBS)),
                 loglq = ifelse(xdata$BLQ == 1, log(lq), big),
                 omegaPKInvPrior = 5 * diag(rep(0.05, 5)),
                 logNeut = log(xdata$NEUT),
                 logCirc0HatPriorMean = log(circ0Hat),
                 logCirc0HatPriorTau = 1 / (circ0CV^2),
                 logMttHatPriorMean = log(mttHat),
                 logMttHatPriorTau = 1 / (mttCV^2),
                 logGammaPriorMean = log(gammaHat),
                 logGammaPriorTau = 1 / (gammaCV^2),
                 tauCirc0PriorAlpha = (circ0PrecMean / circ0PrecSD)^2,
                 tauCirc0PriorBeta = circ0PrecMean / circ0PrecSD^2,
                 tauMttPriorAlpha = (mttPrecMean / mttPrecSD)^2,
                 tauMttPriorBeta = mttPrecMean / mttPrecSD^2
                 )

## create initial estimates
bugsinit <- function(){
list(
     logCLHat = rnorm(1, log(10), 0.2),
     logQHat = rnorm(1, log(15), 0.2),
     logV1Hat = rnorm(1, log(35), 0.2),
     logV2Hat = rnorm(1, log(106), 0.2),
     logDkaHat = rnorm(1, log(1), 0.2),
     omegaPKInv = solve(diag(exp(2 * rnorm(5, log(c(0.07, 0.18, 0.057, 0.16, 0.063)) / 2, 0.2)))),
     tauPK = 1/(runif(1, 0.1, 2)^2),
     alphaHat = exp(rnorm.trunc(1, log(5.0E-4), 0.5, upper = log(1))),
     omegaAlpha = exp(rnorm.trunc(1, log(0.5), 0.5, upper = log(5))),
     logMttHat = rnorm(1, log(mttHat), mttCV),
     logCirc0Hat = rnorm(1, log(circ0Hat), circ0CV),
     logGamma = rnorm(1, log(gammaHat), gammaCV),
     tauMtt = rgamma(1, bugsdata$tauMttPriorAlpha, bugsdata$tauMttPriorBeta),
     tauCirc0 = rgamma(1, bugsdata$tauCirc0PriorAlpha, bugsdata$tauCirc0PriorBeta),
     tauNeut = 1 / (runif(1, 0.1, 1)^2),
##     dfNeut = runif(1, 2, 20),
     logtheta = matrix(rep(log(c(10, 15, 35, 106, 1)),
       ea = bugsdata$nsub), nrow = bugsdata$nsub),
     logMtt = rep(log(mttHat), bugsdata$nsub),
     logCirc0 = rep(log(circ0Hat), bugsdata$nsub),
     logAlpha = rep(log(5.0E-4), bugsdata$nsub)
     )
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "kaHat",
                      "sigmaPK", "omegaPK",
                      "alphaHat", "circ0Hat", "mttHat", "gamma",
                      "omegaAlpha", "omegaCirc0", "omegaMtt",
                      "sigmaNeut")##, "dfNeut")

## Additional variables to monitor
##otherRVs <- c("logCobsCond", "logCobsPred", "thetaPD",
##              "logNeutCond", "logNeutPred")
otherRVs <- c("logCobsCond", "thetaPD",
              "logNeutCond")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("deviance", parametersToPlot)

################################################################################################
## run WinBUGS

n.chains <- 3
n.iter <- 10000
n.burnin <- 5000
n.thin <- 10
system.time(
bugs.fit <- bugs(data = bugsdata, inits = bugsinit, parameters.to.save = parameters,
                 model.file = file.path(exampleDir, paste(modelName, ".txt", sep="")),
                 n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                 n.thin = n.thin, clearWD = TRUE, debug = FALSE,
                 bugs.directory = bugsDir, working.directory = getwd(),
                 useWINE = (.Platform$OS.type == "unix"),
                 WINE = file.path(wineBin,"wine"),
                 newWINE = (.Platform$OS.type == "unix"),
                 WINEPATH = file.path(wineBin,"winepath"))
)

## Uncomment the following to relaod from coda files
##bugs.fit <- R2WinBUGS:::bugs.sims(parameters.to.save = parameters, n.chains = n.chains,
## n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = TRUE)

## save scripts, data and results to a directory
save.model(bugs.fit, modelName)

## Uncomment the following to relaad a previously saved run
## load(paste(modelName,"/",modelName,".fit.Rsave",sep=""))
## load(paste("../solutionsComplete/",modelName,"/",modelName,".fit.Rsave",sep=""))

## rename and reformat MCMC results to facilitate later calculations and plots
sims.array <- bugs.fit$sims.array
posterior <- array(as.vector(sims.array), dim = c(prod(dim(sims.array)[1:2]), dim(sims.array)[3]),
                   dimnames = list(NULL,dimnames(sims.array)[[3]]))

################################################################################################
## posterior distributions of parameters

## open graphics device
pdf(file = file.path(exampleDir, paste(modelName, "/", modelName, "Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

## subset of sims.array containing selected variables
x1 <- sims.array[,,unlist(sapply(c(paste("^", parametersToPlot, "$", sep=""),
                                   paste("^", parametersToPlot, "\\[",sep="")),
                                 grep, x=dimnames(sims.array)[[3]]))]

## create history, density and Gelman-Rubin-Brooks plots, and a table of summary stats
ptable <- parameter.plot.table(x1)
write.csv(signif(ptable, 3), paste(modelName, "/", modelName, ".summary.csv", sep=""))

################################################################################################
## posterior predictive distributions of plasma concentrations

## prediction of future observations in the same subjects, i.e., posterior predictions
## conditioned on observed data from the same subject

pred <- exp(posterior[,grep("logCobsCond\\[",dimnames(posterior)[[2]])])

x1 <- xdata
x1$COBS <- as.numeric(x1$COBS)
x1$type <-rep("observed", nrow(x1))
x2 <- rbind(x1, x1, x1)
x2$COBS <- as.vector(t(apply(pred, 2, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)))
x2$type <- rep(c("5%ile", "median", "95%ile"), ea = nrow(x1))
x1 <- rbind(x1, x2)

doses <- sort(unique(x1$DOSE))
for(dose in doses){
  print(xyplot(COBS ~ TIME | as.factor(ID), x1[x1$DOSE == dose,],
               groups = type,
               panel = panel.superpose, type = c("l", "l", "l", "p"),
               lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
               pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "concentration", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste(dose, " mg", sep = ""), cex = 1.5),
               sub = list(label = "individual predictions", cex = 1.5)))
}

pred <- exp(posterior[,grep("logNeutCond\\[",dimnames(posterior)[[2]])])

x1 <- xdata
x1$NEUT <- as.numeric(x1$NEUT)
x1$type <-rep("observed", nrow(x1))
x2 <- rbind(x1, x1, x1)
x2$NEUT <- as.vector(t(apply(pred, 2, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)))
x2$type <- rep(c("5%ile", "median", "95%ile"), ea = nrow(x1))
x1 <- rbind(x1, x2)

doses <- sort(unique(x1$DOSE))
for(dose in doses){
  print(xyplot(NEUT ~ TIME | as.factor(ID), x1[x1$DOSE == dose,],
               groups = type,
               panel = panel.superpose, type = c("l", "l", "l", "p"),
               lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
               pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "ANC", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste(dose, " mg", sep = ""), cex = 1.5),
               sub = list(label = "individual predictions", cex = 1.5)))
}

## prediction of future observations in a new subject, i.e., posterior predictive distributions

## This takes a while, so let's do a subset for now.
predAll <- apply(posterior[sample(nrow(posterior), size = 100),], 1,
                 get(paste(modelName, "Sim", sep = "")), bugsdata=bugsdata)
save(predAll, file = file.path(exampleDir, modelName, "simsForPPC.Rsave"))
nobs <- sum(xdata$EVID == 0)
predConc <- predAll[1:nobs,]
predNeut <- predAll[nobs + 1:nobs,]

x1 <- xdata[xdata$EVID == 0,]
x1$COBS <- as.numeric(x1$COBS)
x1$type <-rep("observed", nrow(x1))
x2 <- rbind(x1, x1, x1)
x2$COBS <- as.vector(t(apply(predConc, 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)))
x2$type <- rep(c("5%ile", "median", "95%ile"), ea = nrow(x1))
x1 <- rbind(x1, x2)

doses <- sort(unique(x1$DOSE))
for(dose in doses){
  print(xyplot(COBS ~ TIME | as.factor(ID), x1[x1$DOSE == dose,],
               groups = type,
               panel = panel.superpose, type = c("l", "l", "l", "p"),
               lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
               pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "concentration", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste(dose, " mg", sep = ""), cex = 1.5),
               sub = list(label = "population predictions", cex = 1.5)))
}


x1 <- xdata[xdata$EVID == 0,]
x1$NEUT <- as.numeric(x1$NEUT)
x1$type <-rep("observed", nrow(x1))
x2 <- rbind(x1, x1, x1)
x2$NEUT <- as.vector(t(apply(predNeut, 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)))
x2$type <- rep(c("5%ile", "median", "95%ile"), ea = nrow(x1))
x1 <- rbind(x1, x2)

doses <- sort(unique(x1$DOSE))
for(dose in doses){
  print(xyplot(NEUT ~ TIME | as.factor(ID), x1[x1$DOSE == dose,],
               groups = type,
               panel = panel.superpose, type = c("l", "l", "l", "p"),
               lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
               pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "ANC", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste(dose, " mg", sep = ""), cex = 1.5),
               sub = list(label = "population predictions", cex = 1.5)))
}

## Posterior predictive checks for 5th, 50th and 95th population percentiles
## Plasma concentration data

xobs <- xdata[xdata$EVID == 0,]
tSim <- sort(unique(xobs$TIME[!is.na(xobs$COBS)])) # times for simulated data

## Calculate 5th, 50th and 95th percentiles of observed data at specified times

predConc <- predConc[with(xobs, order(DOSE, TIME)),]
xobs <- xobs[with(xobs, order(DOSE, TIME)),]

pctObs <- NULL
  doses <- sort(unique(xobs$DOSE))
  for(dose in doses){
  
  cobs <- apply(sapply(unique(xobs$ID), function(id, times, data){
    data <- data[data$ID == id & !is.na(data$COBS),]
    if(nrow(data) < 2)
      return(rep(NA, length(times)))
    else{
      return(approx(data$TIME, data$COBS, xout = times)$y)
    }
  }, times = tSim, data = xobs[xobs$DOSE == dose,
                     c("ID", "TIME", "COBS")]),
               1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)

  pctObs <- rbind(pctObs,
                  data.frame(DOSE = rep(dose, length(tSim) * 3),
                             TIME = rep(tSim, 3),
                             COBS = as.vector(t(cobs)),
                             type = rep(c("obs5", "obs50", "obs95"), ea = length(tSim))))
}

## Calculate 95% credible intervals for predicted values

concSimByTreatment <- NULL
doses <- sort(unique(xobs$DOSE))
for(dose in doses){

  cpred <- apply(apply(predConc[xobs$DOSE == dose,],
                       2, function(cobs, times, tobs, subject){
                         t(apply(sapply(unique(subject), function(i, cobs, times, tobs, subject){
                           tobs <- tobs[subject == i & !is.na(cobs)]
                           cobs <- cobs[subject == i & !is.na(cobs)]
                           if(length(tobs) < 2)
                             return(rep(NA, length(times)))
                           else{
                             return(approx(tobs, cobs, xout = times)$y)
                           }
                         }, times = times, cobs = cobs, tobs = tobs, subject = subject),
                                 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))
                       }, times = tSim, subject = xobs$ID[xobs$DOSE == dose],
                       tobs = xobs$TIME[xobs$DOSE == dose]),
                 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  nTimes <- length(tSim)
  concSimByTreatment <- rbind(concSimByTreatment,
                              data.frame(DOSE = rep(dose, nTimes * 6),
                                         TIME = rep(tSim, 6),
                                         COBS = as.vector(t(cpred)),
                                         type = rep(c("lci5", "lci50", "lci95", "uci5", "uci50", "uci95"),
                                           ea = length(tSim))))

}

obs <- xobs[c("DOSE", "TIME", "COBS")]
obs$type <- rep("observed", nrow(obs))
ppc <- rbind(obs, pctObs[names(obs)], concSimByTreatment)
ppc <- ppc[!is.na(ppc$COBS),]

  
xyplot(COBS ~ TIME | as.factor(DOSE),
       ppc, groups = type,
       panel = function(x, y, subscripts, groups, ...){
         panel.polygon(c(x[groups[subscripts] == "lci5"], rev(x[groups[subscripts] == "uci5"])),
                       c(y[groups[subscripts] == "lci5"], rev(y[groups[subscripts] == "uci5"])),
                       col = "light blue", alpha = 0.33)
         panel.polygon(c(x[groups[subscripts] == "lci50"], rev(x[groups[subscripts] == "uci50"])),
                       c(y[groups[subscripts] == "lci50"], rev(y[groups[subscripts] == "uci50"])),
                       col = "light blue", alpha = 0.33)
         panel.polygon(c(x[groups[subscripts] == "lci95"], rev(x[groups[subscripts] == "uci95"])),
                       c(y[groups[subscripts] == "lci95"], rev(y[groups[subscripts] == "uci95"])),
                       col = "light blue", alpha = 0.33)
         panel.superpose(ifelse(groups[subscripts] == "observed", jitter(x), x),
                         y, subscripts, groups, ...)
       },
       type = c(rep("l", 6), "p", rep("l", 3)),
       lty = c(rep(1, 6), 0, rep(1, 3)), pch = c(rep(NA, 6), 1, rep(NA, 3)),
       col = c(rep("transparent", 3), rep("red", 3), "black", rep("transparent", 3)),
       lwd = 2,
       scales = list(cex = 1, y = list(relation = "free")),
       xlab = list(label = "time (d)", cex=1.2),
       ylab = list(label = "ME-2 plasma concentration (ng/mL)", cex=1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1))

## Posterior predictive checks for 5th, 50th and 95th population percentiles
## ANC data

xobs <- xdata[xdata$EVID == 0,]
tSim <- sort(unique(xobs$TIME[!is.na(xobs$NEUT)])) # times for simulated data

## Calculate 5th, 50th and 95th percentiles of observed data at specified times

predNeut <- predNeut[with(xobs, order(DOSE, TIME)),]
xobs <- xobs[with(xobs, order(DOSE, TIME)),]

pctObs <- NULL
  doses <- sort(unique(xobs$DOSE))
  for(dose in doses){
  
  yobs <- apply(sapply(unique(xobs$ID), function(id, times, data){
    data <- data[data$ID == id & !is.na(data$NEUT),]
    if(nrow(data) < 2)
      return(rep(NA, length(times)))
    else{
      return(approx(data$TIME, data$NEUT, xout = times)$y)
    }
  }, times = tSim, data = xobs[xobs$DOSE == dose,
                     c("ID", "TIME", "NEUT")]),
               1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)

  pctObs <- rbind(pctObs,
                  data.frame(DOSE = rep(dose, length(tSim) * 3),
                             TIME = rep(tSim, 3),
                             NEUT= as.vector(t(yobs)),
                             type = rep(c("obs5", "obs50", "obs95"), ea = length(tSim))))
}

## Calculate 95% credible intervals for predicted values

simByTreatment <- NULL
doses <- sort(unique(xobs$DOSE))
for(dose in doses){

  ypred <- apply(apply(predNeut[xobs$DOSE == dose,],
                       2, function(yobs, times, tobs, subject){
                         t(apply(sapply(unique(subject), function(i, yobs, times, tobs, subject){
                           tobs <- tobs[subject == i & !is.na(yobs)]
                           yobs <- yobs[subject == i & !is.na(yobs)]
                           if(length(tobs) < 2)
                             return(rep(NA, length(times)))
                           else{
                             return(approx(tobs, yobs, xout = times)$y)
                           }
                         }, times = times, yobs = yobs, tobs = tobs, subject = subject),
                                 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))
                       }, times = tSim, subject = xobs$ID[xobs$DOSE == dose],
                       tobs = xobs$TIME[xobs$DOSE == dose]),
                 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  nTimes <- length(tSim)
  simByTreatment <- rbind(simByTreatment,
                              data.frame(DOSE = rep(dose, nTimes * 6),
                                         TIME = rep(tSim, 6),
                                         NEUT = as.vector(t(ypred)),
                                         type = rep(c("lci5", "lci50", "lci95", "uci5", "uci50", "uci95"),
                                           ea = length(tSim))))

}

obs <- xobs[c("DOSE", "TIME", "NEUT")]
obs$type <- rep("observed", nrow(obs))
ppc <- rbind(obs, pctObs[names(obs)], simByTreatment)
ppc <- ppc[!is.na(ppc$NEUT),]

  
xyplot(NEUT ~ TIME | as.factor(DOSE),
       ppc, groups = type,
       panel = function(x, y, subscripts, groups, ...){
         panel.polygon(c(x[groups[subscripts] == "lci5"], rev(x[groups[subscripts] == "uci5"])),
                       c(y[groups[subscripts] == "lci5"], rev(y[groups[subscripts] == "uci5"])),
                       col = "light blue", alpha = 0.33)
         panel.polygon(c(x[groups[subscripts] == "lci50"], rev(x[groups[subscripts] == "uci50"])),
                       c(y[groups[subscripts] == "lci50"], rev(y[groups[subscripts] == "uci50"])),
                       col = "light blue", alpha = 0.33)
         panel.polygon(c(x[groups[subscripts] == "lci95"], rev(x[groups[subscripts] == "uci95"])),
                       c(y[groups[subscripts] == "lci95"], rev(y[groups[subscripts] == "uci95"])),
                       col = "light blue", alpha = 0.33)
         panel.superpose(ifelse(groups[subscripts] == "observed", jitter(x), x),
                         y, subscripts, groups, ...)
       },
       type = c(rep("l", 6), "p", rep("l", 3)),
       lty = c(rep(1, 6), 0, rep(1, 3)), pch = c(rep(NA, 6), 1, rep(NA, 3)),
       col = c(rep("transparent", 3), rep("red", 3), "black", rep("transparent", 3)),
       lwd = 2,
       scales = list(cex = 1, y = list(relation = "free")),
       xlab = list(label = "time (d)", cex=1.2),
       ylab = list(label = "ANC", cex=1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1))

dev.off()


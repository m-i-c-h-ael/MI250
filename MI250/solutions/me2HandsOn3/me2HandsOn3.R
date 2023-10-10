modelName = "me2HandsOn3" # root names of model file

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

set.seed(10271998) # not required but assures repeatable results

lq <- 10 # Lower limit of quantitation
big <- 1.0E20 # a big number

## get the NONMEM-formatted data file
xdata <- read.csv("fxaNONMEMData.csv", as.is = TRUE)

## Let's just do a subset for now to illustrate the method
##subjectsPhase1 <- with(xdata[xdata$STUDY %in% 1:2,],
##                       tapply(ID, list(STUDY, DOSE), unique))
##subjectsPhase2 <- with(xdata[xdata$STUDY == 3,],
##                       tapply(ID, list(STUDY, DOSE), unique))
## Select first half of each Phase 1 treatment arm and a quarter of the
## Phase 2 patients
##subjects <- c(unlist(lapply(subjectsPhase1, function(x) x[1:(length(x)%/%4)])),
##              unlist(lapply(subjectsPhase2, function(x) x[1:(length(x)%/%4)])))
##xdata <- xdata[xdata$ID %in% subjects,]

## convert those pesky "."'s to NA's
xdata$DV <- as.numeric(xdata$DV)

nobs <- nrow(xdata)
## Create unique sequential ID's for all subjects. This allows the NONMEM practice of
## using the same ID's for different subjects
xdata$origID <- xdata$ID
subject <- 1
id <- xdata$ID
for(i in 2:nobs){
  if(id[i] != id[i-1]) subject <- subject + 1
  xdata$ID[i] <- subject
}

## Row index for start of each individual's data
nobs <- nrow(xdata)
start <- (1:nobs)[!duplicated(xdata$ID)]

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
                 logCobs = ifelse(xdata$DV <= 0, NA, log(xdata$DV)),
                 loglq = ifelse(xdata$BLQ == 1, log(lq), big),
                 omegaInvPrior = 5 * diag(rep(0.05, 5))
                 )

## create initial estimates
bugsinit <- function() 
    list(logCLHat = rnorm(1, log(10), 0.2),
         logQHat = rnorm(1, log(20), 0.2),
         logV1Hat = rnorm(1, log(70), 0.2),
         logV2Hat = rnorm(1, log(70), 0.2),
         logDkaHat = rnorm(1, log(1), 0.2),
         omegaInv = solve(diag(exp(2 * rnorm(5, log(0.25), 0.5)))),
         tau = 1/(runif(1, 0.1, 2)^2)
##         sigma = runif(1, 0.1, 2)
         )

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "kaHat",
                      "sigma", "omega")

## Additional variables to monitor
otherRVs <- c("logCobsCond", "logCobsPred", "theta")

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
x1$DV <- as.numeric(x1$DV)
x1$type <-rep("observed", nrow(x1))
x2 <- rbind(x1, x1, x1)
x2$DV <- as.vector(t(apply(pred, 2, quantile, probs = c(0.05, 0.5, 0.95))))
x2$type <- rep(c("5%ile", "median", "95%ile"), ea = nrow(x1))
x1 <- rbind(x1, x2)

doses <- sort(unique(x1$DOSE[x1$STUDY == 1]))
for(dose in doses){
  print(xyplot(DV ~ TIME | as.factor(ID), x1[x1$STUDY == 1 & x1$DOSE == dose,],
               groups = type,
               panel = panel.superpose, type = c("l", "l", "l", "p"),
               lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
               pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "concentration", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste("study 1  ", dose, " mg", sep = ""), cex = 1.5),
               sub = list(label = "individual predictions", cex = 1.5)))
}

doses <- sort(unique(x1$DOSE[x1$STUDY == 2]))
for(dose in doses){
  print(xyplot(DV ~ TIME | as.factor(ID), x1[x1$STUDY == 2 & x1$DOSE == dose,],
               groups = type,
               panel = panel.superpose, type = c("l", "l", "l", "p"),
               lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
               pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "concentration", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste("study 2  ", dose, " mg", sep = ""), cex = 1.5),
               sub = list(label = "individual predictions", cex = 1.5)))
}

print(xyplot(DV ~ TIME | as.factor(ID), x1[x1$STUDY == 3,],
             groups = type, layout = c(5, 5),
             panel = panel.superpose, type = c("l", "l", "l", "p"),
             lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
             pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
             xlab = list(label = "time (h)", cex = 1.2),
             ylab = list(label = "concentration", cex = 1.2),
             par.strip.text = list(cex = 1),
             strip = function(...) strip.default(..., style = 1),
             main = list(label = paste("study 2 20 mg", sep = ""), cex = 1.5),
             sub = list(label = "individual predictions", cex = 1.5)))

## prediction of future observations in a new subject, i.e., posterior predictive distributions

pred <- exp(posterior[,grep("logCobsPred\\[", dimnames(posterior)[[2]])])

x1 <- xdata
x1$DV <- as.numeric(x1$DV)
x1$type <-rep("observed", nrow(x1))
x2 <- rbind(x1, x1, x1)
x2$DV <- as.vector(t(apply(pred, 2, quantile, probs=c(0.05, 0.5, 0.95))))
x2$type <- rep(c("5%ile", "median", "95%ile"), ea = nrow(x1))
x1 <- rbind(x1, x2)

doses <- sort(unique(x1$DOSE[x1$STUDY == 1]))
for(dose in doses){
  print(xyplot(DV ~ TIME | as.factor(ID), x1[x1$STUDY == 1 & x1$DOSE == dose,],
               groups = type,
               panel = panel.superpose, type = c("l", "l", "l", "p"),
               lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
               pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "concentration", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste("study 1  ", dose, " mg", sep = ""), cex = 1.5),
               sub = list(label = "population predictions", cex = 1.5)))
}

doses <- sort(unique(x1$DOSE[x1$STUDY == 2]))
for(dose in doses){
  print(xyplot(DV ~ TIME | as.factor(ID), x1[x1$STUDY == 2 & x1$DOSE == dose,],
               groups = type,
               panel = panel.superpose, type = c("l", "l", "l", "p"),
               lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
               pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "concentration", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste("study 2  ", dose, " mg", sep = ""), cex = 1.5),
               sub = list(label = "population predictions", cex = 1.5)))
}

print(xyplot(DV ~ TIME | as.factor(ID), x1[x1$STUDY == 3,],
             groups = type, layout = c(5, 5),
             panel = panel.superpose, type = c("l", "l", "l", "p"),
             lty = c(3, 3, 1, 0), col = c("red", "red", "blue", "black"),
             pch = c(NA, NA, NA, 19), lwd = 2, cex = 0.5, scales = list(cex = 1),
             xlab = list(label = "time (h)", cex = 1.2),
             ylab = list(label = "concentration", cex = 1.2),
             par.strip.text = list(cex = 1),
             strip = function(...) strip.default(..., style = 1),
             main = list(label = paste("study 2 20 mg", sep = ""), cex = 1.5),
             sub = list(label = "population predictions", cex = 1.5)))

################################################################################################
## random effect plots

effects <- "theta"
effectsPost <- posterior[,unlist(sapply(c(paste("^", effects, "$", sep = ""),
	paste("^", effects, "\\[", sep = "")), grep, x = dimnames(posterior)[[2]]))]

x1 <- xdata[!duplicated(xdata$ID),]

## rename to improve readability
effects <- c("logCL", "logQ", "logV1", "logV2", "logKa")
effectNames <- colnames(effectsPost)
for(id in x1$ID){
  effectNames <- gsub(paste("theta\\[", id, ",1\\]", sep = ""),
                 paste("logCL\\[", id, "\\]", sep = ""), effectNames)
  effectNames <- gsub(paste("theta\\[", id, ",2\\]", sep = ""),
                 paste("logQ\\[", id, "\\]", sep = ""), effectNames)
  effectNames <- gsub(paste("theta\\[", id, ",3\\]", sep = ""),
                 paste("logV1\\[", id, "\\]", sep = ""), effectNames)
  effectNames <- gsub(paste("theta\\[", id, ",4\\]", sep = ""),
                 paste("logV2\\[", id, "\\]", sep = ""), effectNames)
  effectNames <- gsub(paste("theta\\[", id, ",5\\]", sep = ""),
                 paste("logKa\\[", id, "\\]", sep = ""), effectNames)
}
colnames(effectsPost) <- effectNames

medianEffects <- sapply(effects,
                        function(effect, posterior, n){
                          eta <- t(posterior[,grep(paste(effect, "\\[", sep = ""),
                                                   dimnames(posterior)[[2]])])
                          apply(eta,
                                1, median)
                        }, posterior = effectsPost)
colnames(medianEffects) <- effects
medianEffects <- cbind(x1, medianEffects)

for(effect in effects){

  print(histogram(~ medianEffects[effect], medianEffects,
                  ylab = list(cex = 1.2), xlab = list(label = effect, cex = 1.2),
                  scales = list(cex = 1.2)))
  print(qqmath(~ medianEffects[effect], medianEffects,
               prepanel = prepanel.qqmathline, distribution = qnorm,
               panel=function(x, distribution){
                 panel.qqmath(x, cex = 1.2, pch = 16, distribution = distribution)
                 panel.qqmathline(x, distribution = distribution, col = 3, lwd = 3)
               },
               ylab = list(cex = 1.2),
               xlab = list(label = effect, cex = 1.2),
               scales = list(cex = 1.2)))

}

#########################################################################################
## residuals

## prediction of future observations in the same subjects, i.e., posterior predictions
## conditioned on observed data from the same subject

pred <- posterior[, grep("logCobsCond\\[", dimnames(posterior)[[2]])]

cPred <- apply(pred, 2, median, na.rm = T)

resid <- log(xdata$DV) - cPred

xyplot(resid ~ TIME, xdata, panel =
       function(x, y, ...){
         panel.xyplot(jitter(x, 2), y, ...)
         panel.abline(h = 0, col = "red", lty = 2, lwd = 2)
       },
       cex = 0.5, scales = list(cex = 1),
       xlab = list(label = "time (h)", cex = 1.2),
       ylab = list(label = "log(plasma concentration) residual", cex = 1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1))

xyplot(resid[!is.na(xdata$DV)] ~ cPred[!is.na(xdata$DV)], panel =
       function(x, y, ...){
         panel.xyplot(x, y, ...)
         panel.abline(h = 0, col = "red", lty = 2, lwd = 2)
       },
       cex = 0.5, scales = list(cex = 1),
       xlab = list(label = "log(plasma concentration) predicted", cex = 1.2),
       ylab = list(label = "log(plasma concentration) residual", cex = 1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1))

histogram(~ resid,
          xlab=list(label="log(plasma concentration) residual",cex=1.2),
          ylab = list(cex=1.2),
          scales=list(cex=1.2))
qqmath(~ resid,
       prepanel = prepanel.qqmathline, distribution = qnorm,
       panel = function(x, distribution){
         panel.qqmath(x, cex = 1.2, pch = 16, distribution = distribution)
         panel.qqmathline(x, distribution = distribution, col = 3, lwd = 3)
       },
       ylab = list(label = "log(plasma concentration) residual", cex = 1.2),
       xlab = list(cex = 1.2),
       scales = list(cex = 1.2))

dev.off()


## population PK-PD model for relationship between factor Xa inhibition
## and ME-2 plasma concentration

modelName <- "me2HandsOn1" # root names of model file

if(.Platform$OS.type == "windows") Sys.setenv(HOME = substr(R.home(), 1, 2))

courseDir <- file.path(Sys.getenv("HOME"), "project/MI250")
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
source(file.path(toolsDir, "bugs.tools.R"))
source(file.path(toolsDir, "bgillespie.utilities.R"))

set.seed(10271998) # not required but assures repeatable results

## create WinBUGS data set
bugsdata <- list(
                 nobs=8.00E+01,
                 cobs=c(4.10E+00, 3.77E+00, 4.61E+00, 3.45E+00, 4.56E+00, 3.74E+00,
                   3.58E+00, 4.63E+00, 1.90E+01, 1.78E+01, 1.51E+01, 1.14E+01,
                   1.25E+01, 1.87E+01, 1.97E+01, 2.10E+01, 2.01E+01, 2.97E+01,
                   3.87E+01, 3.48E+01, 4.20E+01, 2.46E+01, 2.84E+01, 4.74E+01,
                   6.18E+01, 6.23E+01, 2.73E+01, 3.84E+01, 3.44E+01, 5.05E+01,
                   5.99E+01, 4.23E+01, 1.07E+02, 9.44E+01, 5.17E+01, 7.43E+01,
                   1.15E+02, 7.21E+01, 6.08E+01, 4.98E+01, 1.92E+02, 7.41E+01,
                   9.48E+01, 1.00E+02, 6.65E+01, 1.40E+02, 1.68E+02, 1.02E+02,
                   1.11E+02, 1.08E+02, 1.65E+02, 1.52E+02, 1.65E+02, 1.12E+02,
                   1.33E+02, 1.05E+02, 1.41E+02, 1.75E+02, 2.66E+02, 2.05E+02,
                   2.84E+02, 2.64E+02, 1.47E+02, 2.19E+02, 2.30E+02, 1.97E+02,
                   3.79E+02, 2.57E+02, 3.71E+02, 2.83E+02, 2.02E+02, 2.69E+02,
                   0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,
                   0.00E+00, 0.00E+00),
                 fxa=c(2.45E+00, 6.05E+00, 5.02E+00, 1.20E+00, 4.98E+00, 3.99E+00, 2.85E+00,
                   3.93E+00, 1.72E+01, 1.48E+01, 8.60E+00, 4.18E+00, 8.25E+00, 1.99E+01,
                   1.29E+01, 1.76E+01, 6.56E+00, 7.38E+00, 2.31E+01, 1.42E+01, 3.26E+01,
                   1.66E+01, 2.34E+01, 3.33E+01, 2.89E+01, 3.08E+01, 2.31E+01, 2.70E+01,
                   1.49E+01, 3.06E+01, 2.84E+01, 1.71E+01, 4.75E+01, 3.81E+01, 2.67E+01,
                   3.37E+01, 5.68E+01, 3.46E+01, 2.85E+01, 2.88E+01, 4.45E+01, 2.69E+01,
                   4.94E+01, 4.37E+01, 2.10E+01, 4.76E+01, 4.59E+01, 3.61E+01, 4.91E+01,
                   2.74E+01, 6.44E+01, 6.00E+01, 5.20E+01, 3.89E+01, 5.95E+01, 3.54E+01,
                   5.00E+01, 5.78E+01, 7.72E+01, 5.89E+01, 5.75E+01, 6.20E+01, 4.83E+01,
                   6.84E+01, 6.02E+01, 4.65E+01, 6.84E+01, 6.13E+01, 6.15E+01, 6.60E+01,
                   5.49E+01, 6.09E+01, 1.20E+01, -4.43E+00, -2.43E+00, -3.28E+00, 1.35E+00,
                   4.20E-01, -1.70E+00, -1.00E+00))

## create initial estimates
bugsinit <- list(
                 list(emax=75, log.ec50=4, gamma=1, sigma=10),
                 list(emax=60, log.ec50=6, gamma=2, sigma=5),
                 list(emax=95, log.ec50=5, gamma=0.5, sigma=20)
                 )

## Specify the variables for which you want history and density plots
parametersToPlot <- c("emax", "ec50", "gamma", "sigma")

## Additional variables to monitor
otherRVs <- c("fxa.pred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("deviance", parametersToPlot)

################################################################################################
## run WinBUGS

n.chains <- 3
n.iter <- 10000
n.burnin <- 4000
n.thin <- 5
bugs.fit <- bugs(data = bugsdata, inits = bugsinit, parameters.to.save = parameters,
                 model.file = file.path(exampleDir, paste(modelName, ".txt", sep="")),
                 n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                 n.thin = n.thin, clearWD = TRUE, debug = FALSE,
                 bugs.directory = bugsDir, working.directory = getwd(),
                 useWINE = (.Platform$OS.type == "unix"),
                 WINE = file.path(wineBin,"wine"),
                 newWINE = (.Platform$OS.type == "unix"),
                 WINEPATH = file.path(wineBin,"winepath"))

## save scripts, data and results to a directory

save.model(bugs.fit, modelName)
## load(paste(modelName,"/",modelName,".fit.Rsave",sep=""))

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
                                   paste("^",parametersToPlot,"\\[",sep="")),
                                 grep,x=dimnames(sims.array)[[3]]))]

## create history, density and Gelman-Rubin-Brooks plots, and a table of summary stats
ptable <- parameter.plot.table(x1)
write.csv(signif(ptable,3),paste(modelName,"/",modelName,".summary.csv",sep=""))

################################################################################################
## posterior predictive distributions

pred <- posterior[, grep("fxa.pred\\[", dimnames(posterior)[[2]])]

x1 <- as.data.frame(bugsdata[c("cobs", "fxa")])
x1$type <- rep("observed",nrow(x1))
x2 <- rbind(x1,x1,x1)
x2$fxa <- as.vector(t(apply(pred,2,quantile,probs=c(0.05,0.5,0.95))))
x2$type <- rep(c("5%ile","median","95%ile"),ea=nrow(x1))
x1 <- rbind(x1,x2)

x1 <- x1[order(x1$type, x1$cobs),]
xyplot(fxa ~ cobs, x1, groups = type,
               panel = panel.superpose, type = c("l","l","l","p"), lty = c(3,3,1,0), col = c(2,2,4,1),
               pch = c(NA,NA,NA,1), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "time-averaged ME-2 plasma concentration (ng/mL)", cex = 1.2),
               ylab = list(label = "time-averaged inhibition of factor Xa activity (%)", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1))

dev.off()


## population PK-PD model for relationship between factor Xa inhibition
## and ME-2 plasma concentration

modelName <- "me2HandsOn2" # root names of model file

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

## get data file
xdata <- read.csv(paste(exampleDir,"/fxa.data.csv",sep=""))

## create WinBUGS data set
bugsdata <- list(
                 nobs = nrow(xdata),
                 nsub = length(unique(xdata$subject)),
                 subject = xdata$subject,
                 cobs = xdata$cobs,
                 fxa = xdata$fxa.inh.obs,
                 cmin = 0, cmax = 1600, nsim = 201)

## create initial estimates
bugsinit <- function() list(
                            emax = min(100, max(0, exp(rnorm(1, log(90), 0.1)))),
                            log.ec50.hat = rnorm(1, log(200), 1),
                            gamma = runif(1, 0.25, 3),
                            sigma = exp(rnorm(1, log(10), 1)),
                            omega.ec50 = exp(rnorm(1, log(0.2), 1)))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("emax", "ec50.hat", "gamma", "sigma", "omega.ec50")

## Additional variables to monitor
otherRVs <- c("ec50", "fxa.cond", "fxa.pred", "fxa.pred2")

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

## prediction of future observations in the same subjects, i.e., posterior predictions
## conditioned on observed data from the same subject

pred <- posterior[, grep("fxa.cond\\[", dimnames(posterior)[[2]])]

x1 <- xdata
x1$type <- rep("observed",nrow(x1))
x2 <- rbind(x1,x1,x1)
x2$fxa.inh.obs <- as.vector(t(apply(pred,2,quantile,probs=c(0.05,0.5,0.95))))
x2$type <- rep(c("5%ile","median","95%ile"),ea=nrow(x1))
x1 <- rbind(x1,x2)

doses <- sort(unique(x1$dose))
for(dose in doses){
  xplot <- x1[x1$dose == dose,]
  print(xyplot(fxa.inh.obs ~ time | as.factor(subject), xplot, groups = type,
               panel = panel.superpose, type = c("l","l","l","p"), lty = c(3,3,1,0), col = c(2,2,4,1),
               pch = c(NA, NA, NA, 1), lwd = 2, cex = 0.5, scales = list(cex=1),
               xlab = list(label = "time (h)", cex = 1.2),
               ylab = list(label = "inhibition of factor Xa activity (%)", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste(dose, "mg"), cex = 1.2),
               sub = list(label = "individual predictions", cex = 1.2)))
}

x1 <- x1[order(x1$subject, x1$type, x1$cobs),]
doses <- sort(unique(x1$dose))
for(dose in doses){
  xplot <- x1[x1$dose == dose,]
  print(xyplot(fxa.inh.obs ~ cobs | as.factor(subject), xplot, groups = type,
               panel = panel.superpose, type = c("l","l","l","p"), lty = c(3,3,1,0), col = c(2,2,4,1),
               pch = c(NA,NA,NA,1), lwd = 2, cex = 0.5, scales = list(cex = 1),
               xlab = list(label = "ME-2 plasma concentration (ng/mL)", cex = 1.2),
               ylab = list(label = "inhibition of factor Xa activity (%)", cex = 1.2),
               par.strip.text = list(cex = 1),
               strip = function(...) strip.default(..., style = 1),
               main = list(label = paste(dose, "mg"), cex = 1.2),
               sub = list(label = "individual predictions", cex = 1.2)))
}

## prediction of future observations in a new subject, i.e., posterior predictive distributions

pred <- posterior[, grep("fxa.pred2\\[", dimnames(posterior)[[2]])]

x1 <- rbind(data.frame(conc = rep(bugsdata$cmin + (0:(bugsdata$nsim-1)) *
                         (bugsdata$cmax - bugsdata$cmin) / (bugsdata$nsim - 1), 3),
                       fxa = as.vector(t(apply(pred, 2, quantile, probs = c(0.05, 0.5, 0.95)))),
                       type = rep(c("5%ile", "median", "95%ile"), ea = bugsdata$nsim)),
            data.frame(conc = xdata$cobs, fxa = xdata$fxa.inh.obs,
                       type = rep("observed", nrow(xdata))))

xyplot(fxa ~ conc, x1, groups = type,
       panel = panel.superpose, type = c("l","l","l","p"), lty = c(3,3,1,0), col = c(2,2,4,1),
       pch = c(NA, NA, NA, 1), lwd = 2, cex = 0.5, scales = list(cex=1),
       xlab = list(label = "ME-2 plasma concentration (ng/mL)", cex = 1.2),
       ylab = list(label = "inhibition of factor Xa activity (%)", cex = 1.2),
       par.strip.text = list(cex = 1),
       strip = function(...) strip.default(..., style = 1),
       main = list(label = "population predictions", cex = 1.2))

################################################################################################
## examples of posterior predictive checking

pred <- posterior[, grep("fxa.pred\\[", dimnames(posterior)[[2]])]

## Prediction of mean peak inhibition of factor Xa activity
epeak <- aggregate(xdata["fxa.inh.obs"], list(subject = xdata$subject, dose = xdata$dose), max)
epeak <- aggregate(epeak["fxa.inh.obs"], list(dose = epeak$dose), mean)
x1 <- xdata[!duplicated(xdata$subject), c("subject", "dose")]
x1 <- x1[order(x1$subject),]
x2 <- sapply(sort(unique(xdata$subject)), function(subject, xdata, pred)
	apply(pred[, xdata$subject == subject], 1, max), xdata = xdata, pred = pred)
doses <- sort(unique(x1$dose))
epeak.pred <- data.frame(dose = rep(doses, ea = nrow(pred)),
	epeak = as.vector(sapply(doses, function(dose, y)
          apply(y[, x1$dose == dose], 1, mean), y = x2)))

## density plot 
densityplot(~ epeak | ordered(paste(dose, "mg"),
                              levels = paste(sort(unique(epeak.pred$dose)), "mg")),
            epeak.pred, subscripts = T,
            panel = function(x, y, subscripts, ...){
              panel.densityplot(x, ...)
              dose = epeak.pred$dose[subscripts[1]]
              epeak.obs = epeak$fxa.inh.obs[epeak$dose == dose]
              panel.abline(v = epeak.obs, col = "red", lwd = 2, lty = 2)
              axis.limits = current.panel.limits()
              ltext(x = axis.limits$xlim[2], y = axis.limits$ylim[2], adj = c(1.05, 1.05),
                    labels=paste("Pr(pred >= obs)\n=",
                      signif(sum(x >= epeak.obs) / length(x), 3), sep = ""))
            },
            xlab = list(label = "mean peak inhibition of factor Xa activity (%)", cex = 1.25),
            ylab = list(cex = 1.25),
            scales = list(x = list(relation = "free"), cex = 1.25),
            par.strip.text = list(cex = 1.25))

## histogram
histogram(~ epeak | ordered(paste(dose, "mg"),
                              levels = paste(sort(unique(epeak.pred$dose)), "mg")),
            epeak.pred, subscripts = T,
            panel = function(x, y, subscripts, ...){
              panel.histogram(x, ...)
              dose = epeak.pred$dose[subscripts[1]]
              epeak.obs = epeak$fxa.inh.obs[epeak$dose == dose]
              panel.abline(v = epeak.obs, col = "red", lwd = 2, lty = 2)
              axis.limits = current.panel.limits()
              ltext(x = axis.limits$xlim[2], y = axis.limits$ylim[2], adj = c(1.05, 1.05),
                    labels=paste("Pr(pred >= obs)\n=",
                      signif(sum(x >= epeak.obs) / length(x), 3), sep = ""))
            },
            xlab = list(label = "mean peak inhibition of factor Xa activity (%)", cex = 1.25),
            ylab = list(cex = 1.25),
            scales = list(x = list(relation = "free"), cex = 1.25), breaks = NULL,
            par.strip.text = list(cex = 1.25))

## Prediction of mean inhibition of factor Xa activity AUC
eauc <- aggregate(1:nrow(xdata), list(subject = xdata$subject, dose = xdata$dose),
                    function(i, xdata)
                    trapez(xdata$time[i], xdata$fxa.inh.obs[i]), xdata = xdata)
names(eauc)[ncol(eauc)] <- "fxa.inh.obs"
eauc <- aggregate(eauc["fxa.inh.obs"], list(dose = eauc$dose), mean)
x1 <- xdata[!duplicated(xdata$subject), c("subject", "dose")]
x1 <- x1[order(x1$subject),]
x2 <- sapply(sort(unique(xdata$subject)),
             function(subject, xdata,pred)
             apply(pred[, xdata$subject == subject], 1,
                   function(x,xx)
                   trapez(xx,x), xx = xdata$time[xdata$subject == subject]),
             xdata = xdata, pred = pred)

doses <- sort(unique(x1$dose))
eauc.pred <- data.frame(dose = rep(doses, ea = nrow(pred)),
	eauc = as.vector(sapply(doses,
          function(dose,y)
          apply(y[, x1$dose==dose], 1, mean), y = x2)))

## density plot
densityplot(~ eauc | ordered(paste(dose, "mg"),
                             levels = paste(sort(unique(eauc.pred$dose)), "mg")),
            eauc.pred, subscripts = T,
            panel = function(x, y, subscripts, ...){
              panel.densityplot(x, ...)
              dose = eauc.pred$dose[subscripts[1]]
              eauc.obs = eauc$fxa.inh.obs[eauc$dose == dose]
              panel.abline(v = eauc.obs, col = "red", lwd = 2, lty = 2)
              axis.limits = current.panel.limits()
              ltext(x = axis.limits$xlim[2], y = axis.limits$ylim[2], adj = c(1.05, 1.05),
                    labels = paste("Pr(pred >= obs)\n=",
                      signif(sum(x >= eauc.obs) / length(x), 3), sep = ""))
            },
            xlab = list(label = "mean inhibition of factor Xa activity AUC", cex = 1.25),
            ylab = list(cex = 1.25),
            scales = list(x = list(relation = "free"), cex = 1.25),
            par.strip.text = list(cex = 1.25))

## histogram
histogram(~ eauc | ordered(paste(dose, "mg"),
                             levels = paste(sort(unique(eauc.pred$dose)), "mg")),
            eauc.pred, subscripts = T,
            panel = function(x, y, subscripts, ...){
              panel.histogram(x, ...)
              dose = eauc.pred$dose[subscripts[1]]
              eauc.obs = eauc$fxa.inh.obs[eauc$dose == dose]
              panel.abline(v = eauc.obs, col = "red", lwd = 2, lty = 2)
              axis.limits = current.panel.limits()
              ltext(x = axis.limits$xlim[2], y = axis.limits$ylim[2], adj = c(1.05, 1.05),
                    labels = paste("Pr(pred >= obs)\n=",
                      signif(sum(x >= eauc.obs) / length(x), 3), sep = ""))
            },
            xlab = list(label = "mean inhibition of factor Xa activity AUC", cex = 1.25),
            ylab = list(cex = 1.25),
            scales = list(x = list(relation = "free"), cex = 1.25), breaks = NULL,
            par.strip.text = list(cex = 1.25))

################################################################################################
## Residual plots

pred <- posterior[, grep("fxa.cond\\[", dimnames(posterior)[[2]])]

x2 <- rbind(xdata, xdata, xdata)
x2$fxa.inh.pred <- as.vector(t(apply(pred, 2, quantile, probs = c(0.05, 0.5, 0.95))))
x2$type <- rep(c("5%ile", "median", "95%ile"), ea = nrow(xdata))
x2$res <- x2$fxa.inh.obs - x2$fxa.inh.pred

xyplot(res ~ fxa.inh.pred, x2, groups = type,
       panel = function(x, y, subscripts, groups,...){
         panel.segments(x[groups[subscripts] == "5%ile"], y[groups[subscripts] == "5%ile"],
                        x[groups[subscripts] == "95%ile"], y[groups[subscripts] == "95%ile"],
                        lty = 1, col = "black")
         panel.xyplot(x[groups[subscripts] == "median"], y[groups[subscripts] == "median"],
                      pch = 19, cex = 0.5)
         panel.abline(h = 0, lty = 2, col = "red", lwd = 2)
       },
       xlab = list(label = "predicted inhibition of factor Xa activity (%)", cex = 1.25),
       ylab = list(label = "residuals", cex = 1.25),
       scales = list(x = list(relation = "free"), cex = 1.25),
       par.strip.text = list(cex = 1.25),
       main = list(label = "individual predictions", cex = 1.5))

xyplot(res ~ fxa.inh.pred | ordered(paste(dose, "mg"), levels = paste(sort(unique(eauc.pred$dose)), "mg")),
       x2, groups = type,
       panel = function(x, y, subscripts,groups, ...){
         panel.segments(x[groups[subscripts] == "5%ile"], y[groups[subscripts] == "5%ile"],
                        x[groups[subscripts] == "95%ile"], y[groups[subscripts] == "95%ile"],
                        lty = 1, col = "black")
         panel.xyplot(x[groups[subscripts] == "median"], y[groups[subscripts] == "median"],
                      pch = 19, cex = 0.5)
         panel.abline(h = 0, lty = 2, col = "red", lwd = 2)
       },
       xlab = list(label = "predicted inhibition of factor Xa activity (%)", cex = 1.25),
       ylab = list(label = "residuals", cex = 1.25),
       scales = list(x = list(relation = "free"), cex = 1.25),
       par.strip.text = list(cex = 1.25),
       main = list(label = "individual predictions", cex = 1.5))

xyplot(res ~ time | ordered(paste(dose, "mg"), levels = paste(sort(unique(eauc.pred$dose)), "mg")),
       x2, groups=type,
       panel=function(x, y, subscripts, groups, ...){
         panel.segments(x[groups[subscripts] == "5%ile"], y[groups[subscripts] == "5%ile"],
                        x[groups[subscripts] == "95%ile"], y[groups[subscripts] == "95%ile"],
                        lty = 1, col = "black")
         panel.xyplot(x[groups[subscripts] == "median"], y[groups[subscripts] == "median"],
                      pch = 19, cex = 0.5)
         panel.abline(h = 0, lty = 2, col = "red", lwd = 2)
       },
       xlab = list(label = "time (h)", cex = 1.25),
       ylab = list(label = "residuals", cex = 1.25),
       scales = list(x = list(relation = "free"), cex = 1.25),
       par.strip.text = list(cex = 1.25),
       main = list(label = "individual predictions", cex = 1.5))

pred <- posterior[,grep("fxa.pred\\[",dimnames(posterior)[[2]])]

x2 <- rbind(xdata, xdata, xdata)
x2$fxa.inh.pred <- as.vector(t(apply(pred, 2, quantile, probs = c(0.05, 0.5, 0.95))))
x2$type <- rep(c("5%ile", "median", "95%ile"), ea = nrow(xdata))
x2$res <- x2$fxa.inh.obs - x2$fxa.inh.pred

xyplot(res ~ fxa.inh.pred, x2, groups = type,
       panel = function(x, y, subscripts, groups, ...){
         panel.segments(x[groups[subscripts] == "5%ile"], y[groups[subscripts] == "5%ile"],
                        x[groups[subscripts] == "95%ile"], y[groups[subscripts] == "95%ile"],
                        lty = 1, col = "black")
         panel.xyplot(x[groups[subscripts] == "median"], y[groups[subscripts] == "median"],
                      pch = 19, cex = 0.5)
         panel.abline(h = 0, lty = 2, col = "red", lwd = 2)
       },
       xlab = list(label = "predicted inhibition of factor Xa activity (%)", cex = 1.25),
       ylab = list(label = "residuals", cex = 1.25),
       scales = list(x = list(relation = "free"), cex = 1.25),
       par.strip.text = list(cex = 1.25),
       main=list(label = "population predictions", cex = 1.5))

xyplot(res ~ fxa.inh.pred | ordered(paste(dose, "mg"), levels = paste(sort(unique(eauc.pred$dose)), "mg")),
       x2, groups = type,
       panel = function(x, y, subscripts, groups,...){
         panel.segments(x[groups[subscripts] == "5%ile"], y[groups[subscripts] == "5%ile"],
                        x[groups[subscripts] == "95%ile"], y[groups[subscripts] == "95%ile"],
                        lty = 1, col = "black")
         panel.xyplot(x[groups[subscripts] == "median"], y[groups[subscripts] == "median"],
                      pch = 19, cex = 0.5)
         panel.abline(h = 0, lty = 2, col = "red", lwd = 2)
       },
       xlab = list(label = "predicted inhibition of factor Xa activity (%)", cex = 1.25),
       ylab = list(label = "residuals", cex = 1.25),
       scales = list(x = list(relation = "free"), cex = 1.25),
       par.strip.text = list(cex = 1.25),
       main = list(label = "population predictions", cex = 1.5))

xyplot(res ~ time | ordered(paste(dose, "mg"), levels = paste(sort(unique(eauc.pred$dose)), "mg")),
       x2, groups = type,
       panel = function(x, y, subscripts, groups, ...){
         panel.segments(x[groups[subscripts] == "5%ile"], y[groups[subscripts] == "5%ile"],
                        x[groups[subscripts] == "95%ile"], y[groups[subscripts] == "95%ile"],
                        lty = 1, col = "black")
         panel.xyplot(x[groups[subscripts] == "median"], y[groups[subscripts] == "median"],
                      pch = 19, cex = 0.5)
         panel.abline(h = 0, lty = 2, col = "red", lwd = 2)
       },
       xlab = list(label = "time (h)", cex = 1.25),
       ylab = list(label = "residuals", cex = 1.25),
       scales = list(x = list(relation = "free"), cex = 1.25),
       par.strip.text = list(cex = 1.25),
       main = list(label = "population predictions", cex = 1.5))

dev.off()


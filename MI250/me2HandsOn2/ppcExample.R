## Posterior predictive checking example for MI250

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
library(lattice)
source(file.path(toolsDir, "bgillespie.utilities.R"))

################################################################################################
## Use data set and MCMC samples from handOn2 example

## get data file
xdata <- read.csv(paste(exampleDir,"/fxa.data.csv",sep=""))

 load(paste(modelName,"/",modelName,".fit.Rsave",sep=""))

## rename and reformat MCMC results to facilitate later calculations and plots

sims.array <- bugs.fit$sims.array
posterior <- array(as.vector(sims.array),dim=c(prod(dim(sims.array)[1:2]),dim(sims.array)[3]),
                   dimnames=list(NULL,dimnames(sims.array)[[3]]))

################################################################################################

## open graphics device
pdf(file = file.path(exampleDir, "ppcPlots%03d.pdf"),
    width = 6, height = 6, onefile = F)

calcFxaInh <- function(p, conc){
  ec50 <- exp(rnorm(1, log(p["ec50.hat"]), p["omega.ec50"]))
  fxa.hat <- p["emax"] * conc^p["gamma"] / (ec50^p["gamma"] + conc^p["gamma"])
  rnorm(length(conc), fxa.hat, p["sigma"])
}

cmin <- 0
cmax <- 1600
nsim <- 201
cpred <- cmin + ((1:nsim) - 1) * (cmax - cmin) / (nsim - 1)

fxa.pred <- apply(posterior, 1, calcFxaInh, conc = cpred)

x1 <- rbind(data.frame(conc = rep(cpred, 3),
                       fxa = as.vector(t(apply(fxa.pred, 1, quantile, probs = c(0.05, 0.5, 0.95)))),
                       type = rep(c("5%ile", "median", "95%ile"), ea = nsim)),
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

## Get predicted data
pred <- posterior[, grep("fxa.pred\\[", dimnames(posterior)[[2]])]

## Calculate mean AUC factor Xa inhibition by treatment arm
stat <- aggregate(1:nrow(xdata), list(subject = xdata$subject, dose = xdata$dose),
                  function(i, xdata)
                  trapez(xdata$time[i], xdata$fxa.inh.obs[i]), xdata = xdata)
names(stat)[ncol(stat)] <- "stat"
stat <- aggregate(stat["stat"], list(dose = stat$dose), mean, na.rm = TRUE)

## Now do it for each MCMC sample
x1 <- xdata[!duplicated(xdata$subject), c("subject", "dose")]
x1 <- x1[order(x1$subject),]
x2 <- sapply(sort(unique(xdata$subject)),
             function(subject, xdata, pred)
             apply(pred[, xdata$subject == subject], 1,
                   function(x, time)
                   trapez(time, x), time = xdata$time[xdata$subject == subject]),
             xdata = xdata, pred = pred)

doses <- sort(unique(x1$dose))
stat.pred <- data.frame(dose = rep(doses, ea = nrow(pred)),
                        stat = as.vector(sapply(doses, function(dose,y)
                          apply(y[, x1$dose==dose], 1, mean), y = x2)))

## density plot 
densityplot(~ stat | as.factor(dose),
            stat.pred, subscripts = T,
            panel = function(x, y, subscripts, ...){
              panel.densityplot(x, ...)
              dose = stat.pred$dose[subscripts[1]]
              stat.obs = stat$stat[stat$dose == dose]
              panel.abline(v = stat.obs, col = "red", lwd = 2, lty = 2)
              axis.limits = current.panel.limits()
              ltext(x = axis.limits$xlim[2], y = axis.limits$ylim[2], adj = c(1.05, 1.05),
                    labels=paste("Pr(pred >= obs)\n=",
                      signif(sum(x >= stat.obs) / length(x), 3), sep = ""), cex = 0.75)
            },
            xlab = list(label = "mean AUC of factor Xa % inhibition", cex = 1.25),
            ylab = list(cex = 1.25),
            scales = list(x = list(relation = "free"), cex = 1.25),
            par.strip.text = list(cex = 1.25))

## histogram
histogram(~ stat | as.factor(dose),
            stat.pred, subscripts = T,
            panel = function(x, y, subscripts, ...){
              panel.histogram(x, ...)
              dose = stat.pred$dose[subscripts[1]]
              stat.obs = stat$stat[stat$dose == dose]
              panel.abline(v = stat.obs, col = "red", lwd = 2, lty = 2)
              axis.limits = current.panel.limits()
              ltext(x = axis.limits$xlim[2], y = axis.limits$ylim[2], adj = c(1.05, 1.05),
                    labels=paste("Pr(pred >= obs)\n=",
                      signif(sum(x >= stat.obs) / length(x), 3), sep = ""), cex = 0.75)
            },
            xlab = list(label = "mean AUC of factor Xa % inhibition", cex = 1.25),
            ylab = list(cex = 1.25),
            scales = list(x = list(relation = "free"), cex = 1.25), breaks = NULL,
            par.strip.text = list(cex = 1.25))

#####################################################################################################

## Calculate mean peak factor Xa inhibition by treatment arm
stat <- aggregate(1:nrow(xdata), list(subject = xdata$subject, dose = xdata$dose),
                  function(i, xdata)
                  max(xdata$fxa.inh.obs[i]), xdata = xdata)
names(stat)[ncol(stat)] <- "stat"
stat <- aggregate(stat["stat"], list(dose = stat$dose), mean, na.rm = TRUE)

## Now do it for each MCMC sample
x1 <- xdata[!duplicated(xdata$subject), c("subject", "dose")]
x1 <- x1[order(x1$subject),]
x2 <- sapply(sort(unique(xdata$subject)),
             function(subject, xdata, pred)
             apply(pred[, xdata$subject == subject], 1, max),
             xdata = xdata, pred = pred)

doses <- sort(unique(x1$dose))
stat.pred <- data.frame(dose = rep(doses, ea = nrow(pred)),
                        stat = as.vector(sapply(doses, function(dose,y)
                          apply(y[, x1$dose==dose], 1, mean), y = x2)))

## density plot 
densityplot(~ stat | as.factor(dose),
            stat.pred, subscripts = T,
            panel = function(x, y, subscripts, ...){
              panel.densityplot(x, ...)
              dose = stat.pred$dose[subscripts[1]]
              stat.obs = stat$stat[stat$dose == dose]
              panel.abline(v = stat.obs, col = "red", lwd = 2, lty = 2)
              axis.limits = current.panel.limits()
              ltext(x = axis.limits$xlim[2], y = axis.limits$ylim[2], adj = c(1.05, 1.05),
                    labels=paste("Pr(pred >= obs)\n=",
                      signif(sum(x >= stat.obs) / length(x), 3), sep = ""), cex = 0.75)
            },
            xlab = list(label = "mean peak factor Xa % inhibition", cex = 1.25),
            ylab = list(cex = 1.25),
            scales = list(x = list(relation = "free"), cex = 1.25),
            par.strip.text = list(cex = 1.25))

## histogram
histogram(~ stat | as.factor(dose),
            stat.pred, subscripts = T,
            panel = function(x, y, subscripts, ...){
              panel.histogram(x, ...)
              dose = stat.pred$dose[subscripts[1]]
              stat.obs = stat$stat[stat$dose == dose]
              panel.abline(v = stat.obs, col = "red", lwd = 2, lty = 2)
              axis.limits = current.panel.limits()
              ltext(x = axis.limits$xlim[2], y = axis.limits$ylim[2], adj = c(1.05, 1.05),
                    labels=paste("Pr(pred >= obs)\n=",
                      signif(sum(x >= stat.obs) / length(x), 3), sep = ""), cex = 0.75)
            },
            xlab = list(label = "mean peak factor Xa % inhibition", cex = 1.25),
            ylab = list(cex = 1.25),
            scales = list(x = list(relation = "free"), cex = 1.25), breaks = NULL,
            par.strip.text = list(cex = 1.25))

dev.off()

dydx2cptNeut <- function(t, y, p){
  ## PK model: 2 compartment with 1st order absorption
  ## PD model: neutropenia model of Friberg & Karlsson
  ## Assumes p is a named vector with elements
  ## named cl, cld, v1, v2, ka, alpha, mtt, circ0 and gamma
  
  k10 <- p["CL"] / p["V1"]
  k12 <- p["Q"] / p["V1"]
  k21 <- p["Q"] / p["V2"]
  ktr <- 4 / p["mtt"]
  conc <- 1000 * y[2] / p["V1"]
  EDrug <- p["alpha"] * conc
  dydx <- double(8)
  K1 <- matrix(c(-p["ka"], p["ka"], 0,
                 0, -(k10 + k12), k12,
                 0, k21, -k21),
               ncol = 3)
  K2 <- matrix(c(ktr, 0, 0, 0,
                 -ktr, ktr, 0, 0,
                 0, -ktr, ktr, 0,
                 0, 0, -ktr, ktr,
                 0, 0, 0, -ktr),
               ncol = 5)
  dydx[1:3] <- K1 %*% y[1:3]
  dydx[4] <- ktr * y[4] * ((1 - EDrug) * (p["circ0"] / y[8])^p["gamma"] - 1)
  dydx[5:8] <- K2 %*% y[4:8]

  list(dydx)
}

getDoseRecords <- function(nmdata){
  ## Expands NONMEM dose records to account for multiple doses
  ## specified by addl and ii.
  ## Only returns the dose records.
  doseRecords <- nmdata[nmdata$evid == 1,]
  for(i in 1:nrow(doseRecords)){
    if(doseRecords$addl[i] > 0 & doseRecords$ii[i] > 0){
      newRecords <- do.call(rbind, lapply(1:doseRecords$addl[i], function(x) doseRecords[i,]))
      newRecords$addl <- rep(0, doseRecords$addl[i])
      newRecords$ii <- rep(0, doseRecords$addl[i])
      newRecords$time <- doseRecords$time[i] +
        cumsum(rep(doseRecords$ii[i], doseRecords$addl[i]))
      doseRecords <- rbind(doseRecords, newRecords)
    }
    doseRecords$addl[i] <- 0
    doseRecords$ii[i] <- 0
  }
  doseRecords[order(doseRecords$subject, doseRecords$time),]
}

me2HandsOn4Sim <- function(p, bugsdata, subjectSpecific = FALSE){

  require(mvtnorm)
  require(odesolve)
  
  if(subjectSpecific){
    theta <- matrix(p[grep("thetaPD\\[", names(p))], ncol = bugsdata$nsub)
  }else{
    omega <- matrix(p[grep("omegaPK\\[", names(p))],
                    nrow = sqrt(length(p[grep("omegaPK\\[", names(p))])),
                    byrow = TRUE)
    theta <- exp(cbind(rmvnorm(bugsdata$nsub,
                               log(c(p["kaHat"], p["CLHat"], p["QHat"], p["V1Hat"], p["V2Hat"])),
                               omega),
                       rnorm(bugsdata$nsub, log(p["alphaHat"]), p["omegaAlpha"]),
                       rnorm(bugsdata$nsub, log(p["mttHat"]), p["omegaMtt"]),
                       rnorm(bugsdata$nsub, log(p["circ0Hat"]), p["omegaCirc0"])))
  }
    colnames(theta) <- c("ka", "CL", "Q", "V1", "V2", "alpha", "mtt", "circ0")
    weight <- bugsdata$weight[!duplicated(bugsdata$subject)]
    theta[,c("CL", "Q")] <- theta[,c("CL", "Q")] * (weight/70)^0.75
    theta[,c("V1", "V2")] <- theta[,c("V1", "V2")] * weight / 70

  data <-  as.data.frame(bugsdata[c("subject", "weight", "time", "amt", "ii",
                                    "evid", "addl")])
  subjects <- unique(data$subject)
  
  ## following assumes that subjects is an integer sequence starting at 1
  xHat <- lapply(subjects,
                function(subject, theta, gamma, data, dydx){
                   theta <- theta[subject,]
                   data <- data[data$subject == subject,]
                   x <- data$time
                   nx <- length(x)
                   doseRecords <- getDoseRecords(data)
                   doses <- doseRecords$amt
                   doseTimes <- doseRecords$time
                  nDoseTimes <- length(doseTimes)
                  y0 <- c(0, 0, 0, rep(theta["circ0"], 5))
                  y <- NULL
                  for(i in 1:nDoseTimes){
                    y0[1] <- y0[1] + doses[i]             
                    xstart <- doseTimes[i]
                    if(i == nDoseTimes) xend <- x[nx] else xend <- doseTimes[i+1]
                    xpred <- sort(unique(c(xstart, xend, x[x > xstart & x < xend])))
                    y <- rbind(y,
                               lsoda(y = y0, times = xpred, func = dydx,
                                     parms = c(theta, gamma = as.vector(gamma)),
                                     rtol = 1.0E-8, atol = 1.0E-10))
                    y0 <- y[nrow(y), -1]
                  }
                  y <- y[!duplicated(y[,1]),]
                  y[y[,1] %in% x, -1]
                },
                dydx = dydx2cptNeut, theta = theta, gamma = p["gamma"],
                data = data)

##  x2Hat <- unlist(lapply(xHat, function(x) x[,2]))
  conc <- unlist(lapply(1:length(xHat),
                        function(i, xHat, theta)
                        1000 * xHat[[i]][,2] / theta[i, "V1"],
                        xHat = xHat, theta = theta))
  neut <- unlist(lapply(xHat, function(x) x[,8]))
##  conc <- 1000 * t(t(x2Hat) / theta[,"V1"])

  cbind(ifelse(conc <= 0, NA, exp(rnorm(length(conc), log(conc), p["sigmaPK"]))),
        exp(rnorm(length(neut), log(neut), p["sigmaNeut"])))
}

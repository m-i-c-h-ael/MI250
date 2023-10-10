trapez = function(x,y){
	nx = length(x)
	0.5*sum((y[-1]+y[-nx])*(x[-1]-x[-nx]))
}

pkpar = function(x,y){
	nx = length(x)
	auc = sum((x[-1]-x[-nx])*(y[-1]+y[-nx]),na.rm=T)/2
	nmax = order(y,na.last=F)[nx]
	cmax = y[nmax]
	tmax = x[nmax]
	c(auc=auc,cmax=cmax,tmax=tmax)
}

frac = function(x) sum(x)/length(x)

swap = function(x, oldval, newval)
{
	if(length(oldval) != length(newval))
		stop("length(oldval)!=length(newval)")
	x1 = match(x, oldval, 0)
	x[as.logical(x1)] = newval[x1]
	x
}

mvdnorm = function(y,mu,sigma2){
# bivariate norma density
	d = length(mu)
	det.sigma2 = det(sigma2)
	inv.sigma2 = solve(sigma2)
	exp(-t(y-mu)%*%inv.sigma2%*%(y-mu)/2)/sqrt(((2*pi)^d)*det.sigma2)
}


# generalized logit and inverse logit functions
logit.inv = function(x, lower = 0, upper = 1)
{
	x1 <- exp(x)
	ifelse(is.infinite(x1),upper,lower + ((upper - lower) * x1)/(1 + x1))
}

logit = function(x, lower = 0, upper = 1)
{
	x <- (x - lower)/(upper - lower)
	log(x/(1 - x))
}

factor2char.data.frame = function(x){
	for(i in 1:ncol(x)){
		if(class(x[[i]])=="factor") x[[i]] = as.character(x[[i]])
	}
	x
}

strip.lead.blanks = function(x){
	w = regexpr(" +",x)
	as.character(ifelse(w==1,substring(x,attr(w,"match.length")+1),x))
}

symmat = function(a){
# generate symmetric matrix from its lower triangle given as
# (a[1,1],a[2,1],a[2,2],a[3,1],a[3,2],a[3,3],...)
	n = (-1+sqrt(1+8*length(a)))/2
	x = matrix(rep(0,n*n),ncol=n)
	x[t(lower.tri(x,T))] = a
	y = x + t(t(lower.tri(x))*x)
	y
}

"nonmem.plot.1" = function(x, xname = "TIME", y = "DV", pred = "PRED", ipred = "IPRE", res = "RES", 
	wres = "WRES", ires = "IRES", iwres = "IWRE", id = "ID", xlab = list(xlab = "time", cex = 1.2), 
	ylab = list(ylab = "concentration", cex = 1.2), layout = c(5, 4), main = NULL, 
	scales = list(cex = 1, relation = "same"), xlim = NULL, ylim = NULL, 
	par.strip.text = list(cex = 1.2), strip = function(...)
	strip.default(..., style = 1), id.plot = T, fit.plot = T)
{
	# diagnostic plots of NONMEM results
	# x = data.frame including columns named ID, TIME, DV, PRED, IPRE, RES, WRES, IRES, 
	# and IWRE. it is usually constructed by reading in a NONMEM table file generated
	# with an NMTRAN statement like (recall that DV, PRED, RES & WRES are defaults):
	#   $TABLE ID TIME IPRED IRES IWRES NOPRINT ONEHEADER FILE=nonmem.pred
	individual.fit = y.vs.pred = y.vs.ipred = res.vs.pred = res.vs.xname = wres.vs.pred = 
		wres.vs.xname = ires.vs.ipred = ires.vs.xname = iwres.vs.ipred = iwres.vs.xname = NULL
	# individual plots of observed and predicted time courses
	xnames = c(id, xname, y, pred, ipred, res, wres, ires, iwres)
	x = x[xnames[xnames %in% names(x)]]
	names(x) = c("id", "xname", "y", "pred", "ipred", "res", "wres", "ires", "iwres")[xnames %in% names(x)]
	x <- x[!is.na(x$ipred),  ]
	x = x[order(x$xname),  ]
	x1 <- x[c("id", "xname", "pred")]
	# x2 <- x[x$y != 0,  ][c("id", "xname", "ipred")]
	# x3 <- x[x$y != 0,  ][c("id", "xname", "y")]
	x2 <- x[c("id", "xname", "ipred")]
	x3 <- x[c("id", "xname", "y")]
	names(x1)[3] <- names(x2)[3] <- "y"
	x1$type <- rep("predicted-population", length(x1[, 1]))
	x2$type <- rep("predicted-individual", length(x2[, 1]))
	x3$type <- rep("observed", length(x3[, 1]))
	xnew <- rbind(x1, x2, x3)
	# x1 <- x3 <- x[x$y != 0,  ]
	x1 = x
	if(id.plot)
		individual.fit <- xyplot(y ~ xname | as.factor(id), groups = type, data = xnew, 
		panel = panel.superpose, lty = c(1, 1, 3), type = c("p", "l", "l"), lwd = 3, xlab = xlab, 
		ylab = ylab, layout = layout, main = main, col = c(1, 6, 8), scales = scales, xlim = xlim, 
		ylim = ylim, par.strip.text = par.strip.text, strip = strip)
	if(fit.plot) {
		# various diagnostic plots
		y.vs.pred <- xyplot(y ~ pred, x1, panel = function(x, y, ...)
		{
			panel.xyplot(x, y, col = 1, ...)
			panel.abline(0, 1, col = 8, lwd = 3, ...)
		}
		, scales = list(cex = 1), xlab = list(xlab = paste("predicted", ylab$ylab), cex = 1.2), 
		ylab = list(ylab = paste("observed", ylab$ylab), cex = 1.2), main = list(main = "population predictions", cex = 1.2))
		y.vs.ipred <- xyplot(y ~ ipred, x1, panel = function(x, y, ...)
		{
			panel.xyplot(x, y, col = 1, ...)
			panel.abline(0, 1, col = 8, lwd = 3, ...)
		}
		, scales = list(cex = 1), xlab = list(xlab = paste("predicted", ylab$ylab), cex = 1.2), 
		ylab = list(ylab = paste("observed", ylab$ylab), cex = 1.2), main = list(main = "individual predictions", cex = 1.2))
		res.vs.pred <- xyplot(res ~ pred, x1, panel = function(x, y, ...)
		{
			panel.xyplot(x, y, col = 1, ...)
			panel.abline(0, 0, col = 8, lwd = 3, ...)
		}
		, scales = list(cex = 1), xlab = list(xlab = paste("predicted", ylab$ylab), cex = 1.2), 
		ylab = list(ylab = "population residuals", cex = 1.2))
		res.vs.xname <- xyplot(res ~ xname, x1, panel = function(x, y, ...)
		{
			panel.xyplot(x, y, col = 1, ...)
			panel.abline(0, 0, col = 8, lwd = 3, ...)
		}
		, scales = list(cex = 1), xlab = list(xlab = xlab$xlab, cex = 1.2), 
		ylab = list(ylab = "population residuals", cex = 1.2))
		wres.vs.pred <- xyplot(wres ~ pred, x1, panel = function(x, y, ...)
		{
			panel.xyplot(x, y, col = 1, ...)
			panel.abline(0, 0, col = 8, lwd = 3, ...)
		}
		, scales = list(cex = 1), xlab = list(xlab = paste("predicted", ylab$ylab), cex = 1.2), 
		ylab = list(ylab = "weighted population residuals", cex = 1.2))
		wres.vs.xname <- xyplot(wres ~ xname, x1, panel = function(x, y, ...)
		{
			panel.xyplot(x, y, col = 1, ...)
			panel.abline(0, 0, col = 8, lwd = 3, ...)
		}
		, scales = list(cex = 1), xlab = list(xlab = xlab$xlab, cex = 1.2), 
		ylab = list(ylab = "weighted population residuals", cex = 1.2))
		if("ires" %in% names(x1)) {
			ires.vs.ipred <- xyplot(ires ~ ipred, x1, panel = function(x, y, ...)
			{
				panel.xyplot(x, y, col = 1, ...)
				panel.abline(0, 0, col = 8, lwd = 3, ...)
			}
			, scales = list(cex = 1), xlab = list(xlab = paste("predicted", ylab$ylab), cex = 1.2), 
			ylab = list(ylab = "individual residuals", cex = 1.2))
			ires.vs.xname <- xyplot(ires ~ xname, x1, panel = function(x, y, ...)
			{
				panel.xyplot(x, y, col = 1, ...)
				panel.abline(0, 0, col = 8, lwd = 3, ...)
			}
			, scales = list(cex = 1), xlab = list(xlab = xlab$xlab, cex = 1.2), 
			ylab = list(ylab = "individual residuals", cex = 1.2))
		}
		if("iwres" %in% names(x1)) {
			iwres.vs.ipred <- xyplot(iwres ~ ipred, x1, panel = function(x, y, ...)
			{
				panel.xyplot(x, y, col = 1, ...)
				panel.abline(0, 0, col = 8, lwd = 3, ...)
			}
			, scales = list(cex = 1), xlab = list(xlab = paste("predicted", ylab$ylab), cex = 1.2), 
			ylab = list(ylab = "weighted individual residuals", cex = 1.2))
			iwres.vs.xname <- xyplot(iwres ~ xname, x1, panel = function(x, y, ...)
			{
				panel.xyplot(x, y, col = 1, ...)
				panel.abline(0, 0, col = 8, lwd = 3, ...)
			}
			, scales = list(cex = 1), xlab = list(xlab = xlab$xlab, cex = 1.2), 
			ylab = list(ylab = "weighted individual residuals", cex = 1.2))
		}
	}
	list(individual.fit = individual.fit, y.vs.pred = y.vs.pred, y.vs.ipred = y.vs.ipred,
		res.vs.pred = res.vs.pred, res.vs.xname = res.vs.xname, wres.vs.pred = wres.vs.pred, 
		wres.vs.xname = wres.vs.xname, ires.vs.ipred = ires.vs.ipred, ires.vs.xname = ires.vs.xname, 
		iwres.vs.ipred = iwres.vs.ipred, iwres.vs.xname = iwres.vs.xname)
}


"print.nonmem.plot" = function(nmplots)
{
	print(nmplots$individual.fit)
	print(nmplots$y.vs.pred, position = c(0, 0, 0.525, 1), more = T)
	print(nmplots$y.vs.ipred, position = c(0.475, 0, 1, 1), more = F)
	print(nmplots$res.vs.pred, position = c(0, 0.475, 0.525, 1), more = T)
	print(nmplots$res.vs.xname, position = c(0.475, 0.475, 1, 1), more = T)
	print(nmplots$wres.vs.pred, position = c(0, 0, 0.525, 0.525), more = T)
	print(nmplots$wres.vs.xname, position = c(0.475, 0, 1, 0.525), more = F)
	print(nmplots$ires.vs.ipred, position = c(0, 0.475, 0.525, 1), more = T)
	print(nmplots$ires.vs.xname, position = c(0.475, 0.475, 1, 1), more = !is.null(nmplots$iwres.vs.ipred))
	print(nmplots$iwres.vs.ipred, position = c(0, 0, 0.525, 0.525), more = T)
	print(nmplots$iwres.vs.xname, position = c(0.475, 0, 1, 0.525), more = F)
	NULL
}

qqhist <- function(x, label = NULL, main=NULL)
{
# plot histogram and quantile-normal quantile plot for vector x

	plot1 <- histogram(~x,
		ylab=list(cex=1.2),	xlab=list(label=label,cex=1.2),
		scales=list(cex=1.2))
	plot2 <- qqmath(~x,prepanel = prepanel.qqmathline,
		panel=function(x){
			panel.qqmath(x,cex=1.2,pch=16)
			panel.qqmathline(x,distribution=qnorm,col=3,lwd=3)
		},
		ylab=list(label=label,cex=1.2),
		xlab=list(label="quantiles of standard normal",cex=1.2),
		scales=list(cex=1.2),ylim=range(x))
	list(plot1=plot1,plot2=plot2)
}

qnorm.trunc = function(p,mean=0,sd=1,lower=-Inf,upper=Inf)
	qnorm(p*pnorm(upper,mean,sd)+(1-p)*pnorm(lower,mean,sd),mean,sd)

rnorm.trunc = function(n,mean=0,sd=1,lower=-Inf,upper=Inf)
	qnorm.trunc(runif(n),mean,sd,lower,upper)

# Prepanel and panel plot functions for error bars taken from R help response by
# Deepayan Sarkar
prepanel.ci <- function(y, ly, uy, subscripts, ...)
{
    y <- as.numeric(y)
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    list(ylim = range(y, uy, ly, finite = TRUE))
}

panel.ci <- function(x, y, ly, uy, subscripts, pch = 16, col.line =
'black', arrow.length = 0.25, ...)
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    panel.arrows(x, ly, x, uy, col = col.line,
                 length = arrow.length, unit = "native",
                 angle = 90, code = 3)
    panel.xyplot(x, y, pch = pch, col.line = col.line, ...)
}


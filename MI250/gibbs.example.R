# simple Gibbs sampling example
# From Gelman, Carlin, Stern & Rubin. Bayesian Data Analysis, 2nd ed. pp. 288-289

# c(y1,y2) ~ bivariate normal(c(mu1,mu2),matrix(c(1,rho,rho,1),ncol=2)),
#            c(mu1,mu2) unknown, rho known
#            (improper) uniform prior on c(mu1,mu2)

y = c(0,0)
rho = 0.8
sigma2 = matrix(c(1,rho,rho,1),ncol=2)

rnorm.cond = function(n,i,mu,sigma2,y){
# random samples from distribution of ith element of a
# multivariate normal distribution (y[i]) conditioned on the
# values of all remaining elements (y[j], j != i)

	d = length(mu)
	det.sigma2 = det(sigma2)
	inv.sigma2 = solve(sigma2)

	cond.coef = diag(d) - solve(diag(diag(inv.sigma2)))%*%inv.sigma2
	mean.cond = mu[i] + sum(cond.coef[i,]*(y[-i]-mu[-i]))
	sd.cond = sqrt(1/inv.sigma2[i,i])
	rnorm(n,mean.cond,sd.cond)
}

gibbs1 = function(nsamp,y,sigma2,mu.init){
# Gibbs sampler for bivariate normal mean with known covariance matrix (sigma2),
# an improper uniform prior and given one data point (y)
	mu = matrix(double(2*nsamp),ncol=2)
	mu[1,] = mu.init
	for(i in 2:nsamp){
		mu[i,1] = rnorm.cond(1,1,y,sigma2,mu[i-1,])
		mu[i,2] = rnorm.cond(1,2,y,sigma2,mu[i,])
	}
	mu
}

# first 10 iterations
nsamp = 10
chain1 = gibbs1(nsamp,y,sigma2,c(2.5,2.5))
chain2 = gibbs1(nsamp,y,sigma2,c(2.5,-2.5))
chain3 = gibbs1(nsamp,y,sigma2,c(-2.5,2.5))
chain4 = gibbs1(nsamp,y,sigma2,c(-2.5,-2.5))

# open graphics device
pdf(file = file.path(getwd(), "gibbsSamplingExample%03d.pdf"), width = 6, height = 6, onefile = F)

plot(chain1[,1],chain1[,2],type="l",col=1,xlim=c(-4,4),ylim=c(-4,4),
     xlab=expression(mu[1]),ylab=expression(mu[2]), cex.lab = 1.5)
points(chain1[1,1],chain1[1,2],cex=2,pch=19,col=1)
lines(chain2[,1],chain2[,2],col=2)
points(chain2[1,1],chain2[1,2],cex=2,pch=19,col=2)
lines(chain3[,1],chain3[,2],col=3)
points(chain3[1,1],chain3[1,2],cex=2,pch=19,col=3)
lines(chain4[,1],chain4[,2],col=4)
points(chain4[1,1],chain4[1,2],cex=2,pch=19,col=4)

# 1000 iterations
nsamp = 990
chain1 = rbind(chain1, gibbs1(nsamp,y,sigma2, chain1[10,]))
chain2 = rbind(chain2, gibbs1(nsamp,y,sigma2, chain2[10,]))
chain3 = rbind(chain3, gibbs1(nsamp,y,sigma2, chain3[10,]))
chain4 = rbind(chain4, gibbs1(nsamp,y,sigma2, chain4[10,]))

# spaghetti plot
plot(chain1[,1],chain1[,2],type="l",col=1,xlim=c(-4,4),ylim=c(-4,4),
     xlab=expression(mu[1]),ylab=expression(mu[2]), cex.lab = 1.5)
points(chain1[1,1],chain1[1,2],cex=2,pch=19,col=1)
lines(chain2[,1],chain2[,2],col=2)
points(chain2[1,1],chain2[1,2],cex=2,pch=19,col=2)
lines(chain3[,1],chain3[,2],col=3)
points(chain3[1,1],chain3[1,2],cex=2,pch=19,col=3)
lines(chain4[,1],chain4[,2],col=4)
points(chain4[1,1],chain4[1,2],cex=2,pch=19,col=4)

# point cloud
plot(chain1[,1],chain1[,2],type="p",col=1,pch=19,cex=0.5,xlim=c(-4,4),ylim=c(-4,4),
     xlab=expression(mu[1]),ylab=expression(mu[2]), cex.lab = 1.5)
points(chain2[,1],chain2[,2],col=2,pch=19,cex=0.5)
points(chain3[,1],chain3[,2],col=3,pch=19,cex=0.5)
points(chain4[,1],chain4[,2],col=4,pch=19,cex=0.5)

dev.off()


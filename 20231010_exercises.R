library('rjags')

modStr= 'model{
  for(i in 1:Ntotal){
    y[i] ~ dnorm(mu[i],1/sig^2)
    mu[i]= a + b*x[i]
  }
  #priors
  a ~ dnorm(0, 1/100^2)
  b ~ dnorm(0, 1/100^2)
  sig ~ dunif(0, 10000)
}'

dat= data.frame(x=1:10, y= c(5.19,6.56,9.19,8.09,7.6,7.08,6.74,9.3,8.98,11.5))
datList= list(x=dat$x,y=dat$y,Ntotal=length(dat$x))

jagMod= jags.model(textConnection(modStr),n.chains = 3,data=datList)
cS= coda.samples(jagMod, n.iter= 10000,variable.names = c('a','b','sig'))

df= data.frame( cS[[1]] )
head(df)
hist(df$a,col='skyblue')

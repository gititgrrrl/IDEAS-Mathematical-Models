### R code from vignette source 'deterministic-models.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: deterministic-models.rnw:104-105
###################################################
require(deSolve)                          #deSolve library needed for this computing session


###################################################
### code chunk number 2: deterministic-models.rnw:109-124
###################################################
sir.model.closed <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
       as.list(params),                   #this argument to "with" lets us use the variable names
       {                                  #the system of rate equations
         dS <- -beta*S*I
         dI <- beta*S*I-gamma*I
         dR <- gamma*I
         dx <- c(dS,dI,dR)                #combine results into a single vector dx
         list(dx)                         #return result as a list
       }
       )
}


###################################################
### code chunk number 3: deterministic-models.rnw:133-136
###################################################
times <- seq(0,120,by=5)                    #function seq returns a sequence
params <- c(beta=0.3,gamma=1/7)             #function c "c"ombines values into a vector
xstart <- c(S=9999/10000,I=1/10000,R=0)     #initial conditions


###################################################
### code chunk number 4: deterministic-models.rnw:140-141
###################################################
out <- as.data.frame(ode(xstart,times,sir.model.closed,params))  #result stored in dataframe


###################################################
### code chunk number 5: deterministic-models.rnw:144-149
###################################################
op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))                  #set graphical parameters
plot(I~time,data=out,type='b')                              #plot the I variable against time
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                  #re-set graphical parameters
plot(I~S,data=out,type='b',yaxt='n',xlab='S')               #plot phase portrait
par(op)                                                     #re-set graphical parameters


###################################################
### exercise number 1
###################################################

# explore system dynamics at different parameter values

params2 <- c(beta=0.5,gamma=1/7)             #larger value of beta
params3 <- c(beta=0.3,gamma=1/8)             #smaller value of gamma
params4 <- c(beta=0.6,gamma=2/7)             #multiply both by 2

out2 <- as.data.frame(lsoda(xstart,times,sir.model.closed,params2))
out3 <- as.data.frame(lsoda(xstart,times,sir.model.closed,params3))
out4 <- as.data.frame(lsoda(xstart,times,sir.model.closed,params4))

op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I~time,data=out,type='l', lwd=2, ylim=c(0, 0.45))
lines(I~time,data=out2,type='l', lwd=2, col='#8566B1')
lines(I~time,data=out3,type='l',lwd=2, col='#498EA4')
lines(I~time,data=out4,type='l', lwd=2, col='#04A472')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I~S,data=out,type='l',lwd=2, log='x',yaxt='n',xlab='S', ylim=c(0.01, 0.45), xlim=c(0.01, 1))
lines(I~S,data=out2,type='l', lwd=2, col='#8566B1')
lines(I~S,data=out3,type='l', lwd=2, col='#498EA4')
lines(I~S,data=out4,type='l', lwd=2, col='#04A472')
legend('topleft', col=c('#8566B1','#498EA4','#04A472'), legend=c('Increase beta', 'Decrease gamma', 'Increase beta and gamma'), lty=1, bty='n')
par(op)

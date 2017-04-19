

f = function(theta){

     ### computes the negative log-likelihood
     ### theta = (mu, b, K0, c, p, a, d, q)
     #
     #mu = theta[1]; b = theta[2]; K0 = theta[3]; c = theta[4]; p = theta[5]; a = theta[6]; d = theta[7]; q = theta[8]
     #if (mu <= 0 || b <= 0 || K0 <= 0 || c <= 0 || p < 1 || a <= 0 || d <= 0 || q < 1) return(99999999)
     #
     #### This is for ETAS (2.3).
     #### lambda = mu b exp(-b(m-m0)) +
     ####    K0 b exp(-b(m-m0)) SUM [(t+c)^(-p)] exp[a(mi - m0)] (x^2+y^2+d)^-q
     #### We assume magnitudes are all relative to m0 already!
     #
     #lam1 = rep(mu * b * exp(-b*hm[1]) , n)
     #for(i in 2:n){
     #	r2	= (x[i]-x[1:(i-1)])^2+(y[i]-y[1:(i-1)])^2 # not using Great Circle Distance
     #   lam1[i] = b * exp(-b*(hm[i])) * (mu + K0 * sum((t[i]-t[1:(i-1)]+c)^(-p)*exp(a*hm[1:(i-1)])*(r2+d)^(-q)))
     #
     #}
     #if(min(lam1) < .00000000001) return(999999999)
     #
     #int1 = mu*x1*y1*t1 + K0 * c^(1-p)/(p-1) * pi * d^(1-q)/(q-1) * sum(exp(a*hm))
     #
     #### Note that this integral is only approximate!
     #### It is over infinite time and space.
     #
     #cat("\n Integral = ", int1, ", n=",n,", negative loglikelihood is ",int1-sum(log(lam1)),".")
     #cat("\n theta = ", theta, "\n")
     #return(int1 - sum(log(lam1)))


     ## theta = (a, b, p)
     #v = theta[4]; 
     a = theta[1]; b = theta[2]; p = theta[3]
     if (p <= 0 || p > 1) return(99999999)

     lam1 = rep((1-p)*muf(x[1],y[1])/t1,n)
     for(i in 2:n){
     	r2	= (x[i]-x[1:(i-1)])^2+(y[i]-y[1:(i-1)])^2 # not using Great Circle Distance
        lam1[i] = (1-p)*muf(x[i],y[i])/t1 + p*sum((a*b/pi)*exp(-a*(t[i]-t[1:(i-1)])-b*r2))
        #lam1[i] = (1-p)*muf(x[i],y[i]) + p*sum((a*b/pi)*exp(-a*(t[i]-t[1:(i-1)])-b*r2))
     }
     if(min(lam1) < .00000000001) return(999999999)

     #int1 = (1-p)*sum(muf(x,y)) + p*n
     int1 = (1-p)*x1*y1*mean(kern.est$z) + p*n

     cat("\n Integral = ", int1, ", n=",n,", negative loglikelihood is ",int1-sum(log(lam1)),".")
     cat("\n theta = ", theta, "\n")
     return(int1 - sum(log(lam1)))

}

## theta = (a, b, p)
#theta1 = c(1, 5, .5)

#theta1 = c(0.08483664, 0.03669618, 0.508506) 
theta1 = c(0.07735132, 0.03056830, 0.5664537) 
theta1 = c(0.1036544984, 0.0004771806, 0.7795091750)
theta1 = c(0.07605042, 0.02923889, 0.5749955) 
theta1 = c(0.07605042, 1.4, 0.5749955) 
theta1 = c(0.07604419, 0.02924009, 0.57702962)

b1 = optim(theta1,f,hessian=T,control=list(maxit=600))
## for a catalog of 1000 events, this line above can take 10-15 min.
## Note that, for suitable parameters, the integral should be roughly equal to n.



theta2 = b1$par
b2 = b1$hess
b3 = sqrt(diag(solve(b2)))  ## this gives the asymptotic standard errors
cat("\n The asymptotic standard errors are:\n",signif(b3,2),"\n")



theta1 <- b1$par



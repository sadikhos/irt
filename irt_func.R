##
# given theta, a, b, c, and the type of model ("3pl" or "norm"), this function returns the probability the examinee will get the item correct. By default the guessing parameter cc is set to 0.
irf <- function( theta = 0      # ability parameter
                ,a = 1          # item discrimination
                ,b = 0          # item difficulty
                ,cc = 0         # guessing parameter
                ,type = "3pl"){
    if (type == "3pl"){
        prob <- cc+(1-cc)/(1+exp(-a*(theta-b)))
    }
    if (type == "norm"){
        prob <- cc+(1-cc)*pnorm(a*(theta-b))
    }
    return(prob)
}


# irfplot - will give a plot of the specified item response function for abilities ranging from -3.5 to +3.5
irfplot <- function( a=1
                    ,b=0
                    ,cc=0
                    ,type="3pl"
                    ,label=""
                    ,dashed= FALSE){
    prob <- numeric(71)
    theta <- rep(0,71)
    for (i in 1:71){
        theta[i]<-(i-36)/10
        prob[i]<-irf(theta[i],a,b,cc,type)
    }
    if (!dashed){
        plot(theta, prob, xlab="theta", ylab="prob", type="l",
             xlim=c(-3.5,3.5), ylim=c(0,1), main=label)
    }
    if (dashed){
        plot(theta,prob,xlab="theta",ylab="prob",type="l",
             xlim=c(-3.5,3.5),ylim=c(0,1),main=label,lty=2)
    }
}


## # Example
## irfplot(a=1.8,b=1.0,cc=0.0,label="item 1 comparison")
## par(new=T)
## irfplot(a=1.8,b=1.0,cc=0.0,type="norm",dashed=TRUE)

##########################################################################################
#
#           Cross-sectional (single time point) IRT, multiple groups
#
##########################################################################################
# simulate data from an exam that follows either the normal-ogive or logistic models.
irtgen <- function( n.s = 10                # number of subjects
                   ,avec = c(1)             # discrimination parameter
                   ,bvec = c(0)             # difficulty parameter
                   ,cvec = c(0)             # guessing parameter
                   ,type = "3pl"
                   ,mtheta = c(-0.5, 0.5)   # ability parameter mean
                   ,sdtheta = 1){           # ability parameter sd
    n <- length(avec)
    data <- matrix(0, nrow = n.s, ncol = n)
    cut.p <- floor(n.s/2)
    if (length(mtheta) > 1) {
        theta1 <- rnorm(cut.p, mean = mtheta[1], sd = sdtheta)
        theta2 <- rnorm(n.s-cut.p, mean = mtheta[2], sd = sdtheta)
        theta <- c(theta1, theta2)
    }
        else {theta <- rnorm(n.s, mean = mtheta, sd = sdtheta)}
    p <- t(apply(as.array(theta), 1, FUN=irf, avec, bvec, cvec, type))
    try <- runif(n.s*n)
    data[try < p] <- 1
    return(data)
}


## # Generate abilities centered at 0 with an AR(1) covariance matrix
## # theta <- mvtnorm::rmvnorm(n=1000, mean = rep(0,5), sigma = V.theta)
## # for two groups with different means at post baseline
## theta1 <- mvtnorm::rmvnorm(n= 500, mean = c(0, rep(mu.vec[1], 4)), sigma = V.theta, method = "chol")
## theta2 <- mvtnorm::rmvnorm(n= 500, mean = c(0, rep(mu.vec[2], 4)), sigma = V.theta, method = "chol")
## theta <- rbind(theta1, theta2)

# simulate longitudinal IRT data that follows either the normal-ogive or logistic models.
# produces IRT data in a long format (time as rows)
lirtgen <- function( n.s = 10               # number of subjects
                   ,amat = c(1)             # discrimination parameter
                   ,bmat = c(0)             # difficulty parameter
                   ,cmat = c(0)             # guessing parameter
                   ,type = "3pl"
                   ,mtheta = c(-1.5, 1.5)   # ability parameter mean vector
                   ,vcov = diag(length(mtheta))  # ability parameter cov matrix
                   ,trend = c("linear", "quadratic")
                    ){
    n.times <- dim(vcov)[1]
    n <- nrow(amat)                         # number of items
    cut.n <- floor(n.s/2)
    p <- matrix(0, nrow = n.s, ncol = n)
    out <- vector("list", n.times)          # output data set
    if (trend=="linear") {
    mean.vec1 <- c(0, 1:(n.times-1)*mtheta[1]/(n.times-1))
    mean.vec2 <- c(0, 1:(n.times-1)*mtheta[2]/(n.times-1))
    }
    else if (trend == "quadratic") {
        x.mat <- rbind(rep(0, n.times),
                        0:(n.times-1),
                       (0:(n.times-1))^2)
    beta1.1 <- 2*mtheta[1]/(n.times-1)
    beta2.1 <- (mtheta[1]-(n.times-1)*beta1.1)/(n.times-1)^2
    mean.vec1 <- t(c(0, beta1.1, beta2.1))%*%x.mat
    beta1.2 <- 2*mtheta[2]/(n.times-1)
    beta2.2 <- (mtheta[2]-(n.times-1)*beta1.2)/(n.times-1)^2
    mean.vec2 <- t(c(0, beta1.2, beta2.2))%*%x.mat
    }
    if (length(mtheta) > 1) {

        # set the first visit (baseline) means equal to 0 - no trt effect
        theta1 <- mvtnorm::rmvnorm(n= cut.n, mean = mean.vec1, sigma = vcov)
        theta2 <- mvtnorm::rmvnorm(n= n.s - cut.n, mean = mean.vec2, sigma = vcov)
        theta <- rbind(theta1, theta2)
    }
        else {theta <- mvtnorm::rmvnorm(n= n.s, mean = rep(mtheta, n.times), sigma = vcov)
          }

    try <- matrix(runif(n.s*n), ncol = n)
    for (i in 1:n.times){
        p <- t(apply(as.array(theta[, i]), 1, FUN=irf, amat[, i], bmat[, i], cmat[, i], type))
        out[[i]] <- (try < p)*1
    }
    return(list(as.data.frame(unlist(do.call("rbind", out))),
                theta))
}


### Fast calculation of the A statistic
A <- function(x, y)
{
   nx <- length(x)
   ny <- length(y)
   rx <- sum(rank(c(x, y))[1:nx])
   return((rx / nx - (nx + 1) / 2) / ny)
}

### Calculate the SE and construct a CI for the A statistic using bootstrap methods
Bootstrap.SE.CI.A <- function(x, y, B = 1999, Conf.Level = .95, seed = 1)
{
   # initialize variables
   set.seed(seed)
   nx <- length(x)
   ny <- length(y)
   A.obs <- A(x, y)
   Alpha <- 1 - Conf.Level
   CI.Lower <- CI.Upper <- pi

   # perform bootstrap to generate B values of A
   BS.Values <- rep(0, B)
   for (i in 1:B)
      BS.Values[i] <- A(sample(x, replace = T), sample(y, replace = T))
   BS.Values <- sort(BS.Values)

   # if all bootstrap samples yield same value for A, use it for both ends of CI
   if (min(BS.Values) == max(BS.Values))
      CI.Lower <- CI.Upper <- BS.Values[1]

   # if sample value not within range of bootstrap values, revert to percentile CI
   if ((A.obs < min(BS.Values)) | (A.obs > max(BS.Values)))
   {
      CI.Lower <- BS.Values[round((Alpha / 2) * B)]
      CI.Upper <- BS.Values[round((1 - Alpha / 2) * B)]
   }

   # otherwise, use BCA CI
   if ((CI.Lower == pi) & (CI.Upper == pi))
   {
      # calculate bias-correction and acceleration parameters (z0 and a)
      z0 <- qnorm(mean(BS.Values < A.obs))

      jk <- rep(0, (nx + ny))
      for (i in 1:nx)
         jk[i] <- A(x[-i], y)
      for (i in 1:ny)
         jk[nx + i] <- A(x, y[-i])
      Diff <- mean(jk) - jk
      a <- sum(Diff ^ 3) / (6 * (sum(Diff ^ 2)) ^ 1.5)

      # adjust location of endpoints
      Alpha1 <- pnorm(z0 + (z0 + qnorm(Alpha/2)) / (1 - a * (z0 + qnorm(Alpha/2))))
      Alpha2 <- pnorm(z0 + (z0 - qnorm(Alpha/2)) / (1 - a * (z0 - qnorm(Alpha/2))))

      # if either endpoint undefined, replace it with value for percentile CI
      if (is.na(Alpha1)) Alpha1 <- Alpha / 2
      if (is.na(Alpha2)) Alpha2 <- 1 - Alpha / 2

      if (round(Alpha1 * B) < 1) CI.Lower <- BS.Values[1]
         else CI.Lower <- BS.Values[round(Alpha1 * B)]	
      CI.Upper <- BS.Values[round(Alpha2 * B)]	
   }

   # return A, SE of A, lower limit of CI, upper limit of CI
   return(list(A = A.obs, SE = sd(BS.Values), LL = CI.Lower, UL = CI.Upper))
}


### ==================================================
### This file contains the functions used by NorInAliS
### ==================================================


## Auxiliary functions
## -------------------


# Combines text strings
"%+%" <- function(string1, string2) paste0(string1, string2)


# Removes an element of a set
"%-%" <- function(arg1, arg2) sort(unique(arg1[which(!(arg1 %in% na.omit(arg2)))]))


# Calculates the union of two sets (vectors)
"%A%" <- function(set1, set2)
  if (is.null(set1)) logical(0) else as.vector(na.omit(set1[set1 %in% set2]))


# I just find this more intuitive...
"%contains%" <- function (textstring, searchtext) grepl(searchtext, textstring)


# Calculates a running mean
running.mean <- function(x, n) {
  z <- y <- x
  if (n > 1) {
    for (i in 2:n) {
      z <- c(z[-1], NA)
      y <- y + z
    }
  }
  return(y / n)
}



## Functions to calculate indicator S2
## -----------------------------------


asFreq  <- function(x) {
  # translates a frequency interval into a numerical value
  # Note that input is in events per decade,
  #  whereas output is in events per year!
  sapply(x, function(z) switch(z,
                               "<1"=0.07, "1-8"=0.5, "9-19"=1.4, ">19"=4.6, NA))
}


asAbund <- function(x) {
  # translates an abundance interval into a numerical value
  sapply(x, function(z) switch(z,
                               "1"=1, "2-10"=5, "11-100"=50, "101-1000"=500, ">1000"=5000, NA))
}


rincr <- function(n, min, max, r=TRUE) {
  # generates random numbers according to an increasing
  # triangular probability distribution
  min + sqrt(    if(r) runif(n) else 0.5) * (max - min)
}


rdecr <- function(n, min, max, r=TRUE) {
  # generates random numbers according to a decreasing
  # triangular probability distribution
  max - sqrt(1 - if(r) runif(n) else 0.5) * (max - min)
}


assignFrequencies <- function(n, pathways, pw, r=TRUE) {
  # assigns unknown frequencies according to the distribution of known frequencies,
  # separately for main pathway categories
  freq <- table(
    pathways$Freq[which(pathways$Freq != "unknown" & pathways$Cat == pw)]
  )[c("<1", "1-8", "9-19", ">19")]
  if (any(is.na(freq))) {
    freq[which(is.na(freq))] <- 0
  }
  freq <- freq / sum(freq)
  for (i in 4:2) {
    freq[i] <- sum(freq[1:i])
  }
  seq <- if (r) runif(n) else (1:n - 0.5) / n
  return(c("<1", "1-8", "9-19", ">19")[sapply(seq, function(x) which(freq > x)[1])])
}


S2a <- function(pathways) {
  # estimates S2(a) according to Equation ¤1
  return(length(unique(
    (pathways$Name %+% pathways$Subcat)[pathways$Introd & pathways$Time == "current"]
  )))
}


S2b <- function(pathways, nsim=0, CL=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  # estimates S2(b) according to Equation ¤2
  pathways <- pathways[which(pathways$Introd & pathways$Time == "current"),]
  results <- matrix(0, nrow(pathways), max(1, nsim))
  freq <- matrix(pathways$Freq, nrow(pathways), max(1, nsim))
  for (pw in unique(pathways$Cat)) {
    w <- which(pathways$Freq == "unknown" & pathways$Cat == pw)
    if (length(w)) {
      freq[w,] <- assignFrequencies(length(w) * max(1, nsim),
                                    pathways, pw, r = nsim >= 1)
    }
  }
  W <- freq == "<1"
  results[W] <-                rincr(sum(W), 0.01,  0.1, r = nsim >= 1)
  W <- freq == "1-8"
  results[W] <- if (nsim >= 1) runif(sum(W), 0.10,  0.9) else 0.5
  W <- freq == "9-19"
  results[W] <- if (nsim >= 1) runif(sum(W), 0.90,  1.9) else 1.4
  W <- freq == ">19"
  results[W] <-                rdecr(sum(W), 1.90, 10.0, r = nsim >= 1)
  if (nsim >=1 ) {
    freq <- c(mean(apply(results, 2, sum)),
              sd(apply(results, 2, sum)),
              quantile(apply(results, 2, sum), CL))
    names(freq) <- c("Average", "St.Dev.", (CL * 100) %+% "%CL")
    return(freq)
  } else {
    return(sum(asFreq(freq)))
  }
}


S2c <- function(pathways) {
  # estimates S2(c) according to Equation ¤3
  # (Note: this is a minimum implementation of indicator definition S2c!
  #  It illustrates the way S2c would be estimated, but it currently ignores
  #  uncertainty and simply assumes that the (utterly few) abundances reported
  #  are in fact representative of the unknown abundances, which is unlikely.)
  freq <- asFreq(pathways$Freq[which(pathways$Introd & pathways$Time == "current")])
  w <- which(is.na(freq))
  if (length(w)) {
    mod <- lm(asFreq(Freq) ~ Cat, data=pathways, subset=which(Freq != "unknown"))
    freq[w] <- predict(mod, data.frame(Cat = pathways$Cat[w]))
  }
  abund <- asAbund(pathways$Abund[pathways$Introd & pathways$Time == "current"])
  w <- which(is.na(abund))
  if (length(w)) {
    abund[w] <- mean(abund[-w])
  }
  return(sum(freq * abund))
}



## Functions to calculate indicator S3
## -----------------------------------


asFreq  <- function(x) {
  # translates a frequency interval into a numerical value
  # Note that input is in events per decade,
  #  whereas output is in events per year!
  sapply(x, function(z) switch(z,
                               "<1"=0.07, "1-8"=0.5, "9-19"=1.4, ">19"=4.6, NA))
}


asAbund <- function(x) {
  # translates an abundance interval into a numerical value
  sapply(x, function(z) switch(z,
                               "1"=1, "2-10"=5, "11-100"=50, "101-1000"=500, ">1000"=5000, NA))
}


rincr <- function(n, min, max, r=TRUE) {
  # generates random numbers according to an increasing
  # triangular probability distribution
  min + sqrt(    if(r) runif(n) else 0.5) * (max - min)
}


rdecr <- function(n, min, max, r=TRUE) {
  # generates random numbers according to a decreasing
  # triangular probability distribution
  max - sqrt(1 - if(r) runif(n) else 0.5) * (max - min)
}


assignFrequencies <- function(n, pathways, pw, r=TRUE) {
  # assigns unknown frequencies according to the distribution of known frequencies,
  # separately for main pathway categories
  freq <- table(
    pathways$Freq[which(pathways$Freq != "unknown" & pathways$Cat == pw)]
  )[c("<1", "1-8", "9-19", ">19")]
  if (any(is.na(freq))) {
    freq[which(is.na(freq))] <- 0
  }
  freq <- freq / sum(freq)
  for (i in 4:2) {
    freq[i] <- sum(freq[1:i])
  }
  seq <- if (r) runif(n) else (1:n - 0.5) / n
  return(c("<1", "1-8", "9-19", ">19")[sapply(seq, function(x) which(freq > x)[1])])
}


S3a <- function(pathways) {
  # estimates S3(a) according to Equation ¤1
  return(length(unique(
    (pathways$Name %+% pathways$Subcat)[!pathways$Introd & pathways$Time == "current"]
  )))
}


S3b <- function(pathways, nsim=0, CL=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  # estimates S3(b) according to Equation ¤2
  pathways <- pathways[which(!pathways$Introd & pathways$Time == "current"),]
  results <- matrix(0, nrow(pathways), max(1, nsim))
  freq <- matrix(pathways$Freq, nrow(pathways), max(1, nsim))
  for (pw in unique(pathways$Cat)) {
    w <- which(pathways$Freq == "unknown" & pathways$Cat == pw)
    if (length(w)) {
      freq[w,] <- assignFrequencies(length(w) * max(1, nsim),
                                    pathways, pw, r = nsim >= 1)
    }
  }
  W <- freq == "<1"
  results[W] <-                rincr(sum(W), 0.01,  0.1, r = nsim >= 1)
  W <- freq == "1-8"
  results[W] <- if (nsim >= 1) runif(sum(W), 0.10,  0.9) else 0.5
  W <- freq == "9-19"
  results[W] <- if (nsim >= 1) runif(sum(W), 0.90,  1.9) else 1.4
  W <- freq == ">19"
  results[W] <-                rdecr(sum(W), 1.90, 10.0, r = nsim >= 1)
  if (nsim >=1 ) {
    freq <- c(mean(apply(results, 2, sum)),
              sd(apply(results, 2, sum)),
              quantile(apply(results, 2, sum), CL))
    names(freq) <- c("Average", "St.Dev.", (CL * 100) %+% "%CL")
    return(freq)
  } else {
    return(sum(asFreq(freq)))
  }
}


S3c <- function(pathways) {
  # estimates S3(c) according to Equation ¤3
  # (Note: this is a minimum implementation of indicator definition S3c!
  #  It illustrates the way S3c would be estimated, but it currently ignores
  #  uncertainty and simply assumes that the (utterly few) abundances reported
  #  are in fact representative of the unknown abundances, which is unlikely.)
  freq <- asFreq(pathways$Freq[!pathways$Introd & pathways$Time == "current"])
  w <- which(is.na(freq))
  if (length(w)) {
    freq[w] <- mean(freq[-w])
  }
  abund <- asAbund(pathways$Abund[!pathways$Introd & pathways$Time == "current"])
  w <- which(is.na(abund))
  if (length(w)) {
    abund[w] <- mean(abund[-w])
  }
  return(sum(freq * abund))
}



## Functions to calculate indicator P1
## -----------------------------------


as.nr <- function(x)
  # assigns numerical values from 1 to 5 to the five ecological impact categories
  return(sapply(x, function(z) switch(z, "SE"=5,"HI"=4,"PH"=3,"LO"=2,"NK"=1,NA)))


P1a <- function(dataset, column = "Impact")
  # estimates P1(a)
  return(length(which(dataset[, column] %in% c("SE", "HI", "PH", "LO", "NK"))))


P1b <- function(dataset, column = c("minImp", "Impact", "maxImp")) {
  # estimates P1(b)
  q0 <-       sum(as.nr(dataset[, column[1]]))  # minimum and 2.5% confidence level
  q1 <- round(sum(as.nr(dataset[, column[2]]) * 0.75 +
                    as.nr(dataset[, column[1]]) * 0.25))  # 1st quartile
  q2 <-       sum(as.nr(dataset[, column[2]]))  # best estimate (mean and median)
  q3 <- round(sum(as.nr(dataset[, column[2]]) * 0.75 +
                    as.nr(dataset[, column[3]]) * 0.25))  # 3rd quartile
  q4 <-       sum(as.nr(dataset[, column[3]]))  # maximum and 97.5% confidence level
  SD <- round(sd(c(rep(q0,   floor(P1a(dataset, column[2]) / 4)),
                   rep(q2, ceiling(P1a(dataset, column[2]) / 2)),
                   rep(q4,   floor(P1a(dataset, column[2]) / 4)))))
  # approximation of the S.D.
  p1b <- c(q2, SD, q0, q1, q2, q3, q4)
  names(p1b) <- c("Average", "St.Dev.", c(2.5, 25, 50, 75, 97.5) %+% "%CL")
  return(p1b)
}



## Functions to calculate indicator P2
## -----------------------------------


# The first two functions use maximum-likelihood estimation to infer the 
# standard deviation of the AOOs, based on the best estimate (median), 
# low estimate (1st quartile) and high estimate (3rd quartile), 
# and assuming a log-normal distribution.


f <- function(s, mean, q1, q3) return(
  (q1 - qlnorm(0.25, log(mean) - exp(2*s)/2, exp(s)))^2 +
    (q3 - qlnorm(0.75, log(mean) - exp(2*s)/2, exp(s)))^2
)


findSD <- function(Ex, q1, q3) return(
  exp(optimise(f, c(-12, 12), mean=Ex, q1=q1, q3=q3)$min)
)


P2 <- function(dataset,
               column = c("known", "low", "best", "high"),
               nsim = 100000, 
               maxArea = 323800) {
  # simulates AOOs for all species, which is the basis of indicator P2
  N <- nsim     # random numbers per species
  M <- maxArea  # maximum possible area in km^2
  aoo <- dataset[, column]
  AOO <- matrix(0, N, nrow(aoo))
  for (i in 1:nrow(aoo)) {
    if (any(is.na(aoo[i, c("low", "best", "high")]))) {
      if (is.na(aoo$best[i])) {
        if (is.na(aoo$known[i])) {
          # if no AOO is provided, assume it is 0
          AOO[, i] <- 0
        } else {
          # of no total AOOs are provided, assume they equal the known AOO
          AOO[, i] <- aoo$known[i]
        }
      } else {
        # if no low and high estimates are provided, use the best estimate
        AOO[, i] <- aoo$best[i]
      }
    } else {
      if (aoo$best[i] == 0) {
        # if the best estimate is 0, keep it 0
        AOO[, i] <- 0
      } else {
        # if low, best and high estimates are provided, estimate the standard
        # deviation and generate log-normally distributed random numbers
        SD <- findSD(aoo$best[i], aoo$low[i], aoo$high[i])
        AOO[, i] <- qlnorm(runif(N), log(aoo$best[i]) - SD * SD / 2, SD)
        AOO[, i] <-  sapply(AOO[, i], min, M)  # constrain to area of Norway
        AOO[, i] <- ceiling(AOO[, i] / 4) * 4  # ensure multiples  of 4 km^2
      }
    }
  } # i 
  return(AOO)
} # P2


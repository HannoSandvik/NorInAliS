
# This code calculates indicator S3 of NorInAliS



# Load data from the Alien Species List 2018
# NB: These datasets are not part of NorInAliS and have to be downloaded separately!
# The "fab"  data are available from https://doi.org/10.5061/dryad.8sf7m0cjc
# The "path" data are available from https://doi.org/10.5061/dryad.4b8gthtg7
# After downloading you have to either place these datasets in the working directory
# or adjust the file name/path in the commands!
if (file.exists("assess.txt")) {
  fab  <- read.csv2("assess.txt", as.is=T)
} else {
  cat("Please download \"assess.txt\" from https://doi.org/10.5061/dryad.8sf7m0cjc\n")
}
if (file.exists("path.csv")) {
  path <- read.csv2("path.csv", as.is=T)
} else {
  cat("Plase download \"path.csv\" from https://doi.org/10.5061/dryad.4b8gthtg7\n")
}


# Restrict data to alien species that are reproducing unaidedly in mainland Norway
fab <- fab[which(fab$Status == "reproducing" & fab$Mainl),]


# Restrict data to current pathways of secondary spread
path <- path[which(!path$Introd & path$Time == "current" & path$Name %in% fab$Name), ]


# Load auxiliary functions
eval(parse(text=readLines("function.r")))


# Define additional auxiliary functions
asFreq  <- function(x) {
  # translates a frequency interval into a numerical value
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
  #  are in fact representatative of the unknown abundances, which is unlikely.)
  freq <- asFreq(pathways$Freq[which(!pathways$Introd & pathways$Time == "current")])
  w <- which(is.na(freq))
  if (length(w)) {
    freq[w] <- mean(freq[-w])
  }
  abund <- asAbund(pathways$Abund[which(!pathways$Introd & pathways$Time == "current")])
  w <- which(is.na(abund))
  if (length(w)) {
    abund[w] <- mean(abund[-w])
  }
  return(sum(freq * abund))
}



{# Summarise available data on frequency and abundance
  N <- length(which(!path$Introd & path$Time == "current"))
  cat("Distribution of frequencies of seconadry spread events (per decade):\n")
  print(table(path$Freq )[c(1, 3, 4, 2, 5   )] / N)
  cat("\n\nDistribution of abundances (individuals) per secondary spread event:\n")
  print(table(path$Abund)[c(4, 3, 2, 1, 6)] / N)
}


# Show how frequencies of secondary spread vary between main pathway categories
summary(lm(asFreq(Freq) ~ Cat, data=path, subset=which(Freq != "unknown")))


# Tabulate frequencies separately for each main pathway category
for (pw in unique(path$Cat)) {
  cat("\n\n" %+% pw %+% ":\n")
  print(table(path$Freq[which(path$Cat == pw)]))
}
# This shows that "corridor" has only 1 occurrence with unknown frequencies.
# Therefore, corridors are excluded from S3b below.


{# Estimate the three indicator definitions for the current Alien Species List
  cat("\nIndicator value for S3(a): " %+% S3a(path), "\n")
  cat("\nIndicator value for S3(b):\n")
  print(S3b(path[path$Cat != "corridor",], nsim=100000))  # this one takes a while!
  cat("\nIndicator value for S3(c): " %+% round(S3c(path)), "\n\n")
}


# Disaggregate indicator S3(b) for main pathway categories
for (pw in c("release", "contaminant", "stowaway", "unaided")) {
  cat("\n\n" %+% pw %+% ":\n")
  print(S3b(path[which(path$Cat == pw),], nsim=100000))
}



# Load data from the Alien Species List 2018
fab <- read.csv2("c:\\art\\publ.z\\ecosolev.20\\assess.txt", as.is=T) #¤
aoo <- read.csv2("c:\\art\\neobindc\\aoo.txt", as.is=T) #¤

# Load auxiliary functions
eval(parse(text=readLines("c:\\art\\neobindc\\function.r"))) #¤

# Restrict data to alien species that are reproducing unaidedly in mainland Norway
fab <- fab[which(fab$Status == "reproducing" & fab$Mainl),]

# Restrict the data to freshwater species
w <- which(fab$LifeSt %in% c("lim", "lim,mar", "lim,par", "lim,mar,par"))
fab <- fab[w,]
aoo <- aoo[w,]


# Make sure that the two data frames are compatible
if (all(fab$Name == aoo$Name)) {
  cat("Everything is fine.\n")
} else {
  cat("ERROR: For some reason, the two dataframes are not compatible!\n")
}

# Check for obvious errors in the original data:
# (1) Is any low estimate greater than the corresponding best estimate? 
w <- which(aoo$low > aoo$best)
if (length(w)) {
  print(aoo[w,])
} else {
  cat("Everything is fine.\n")
}

# (2) Is any high estimate less than the corresponding best estimate?
w <- which(aoo$high < aoo$best)
if (length(w)) {
  print(aoo[w,])
} else {
  cat("Everything is fine.\n")
}

# (3) By definition, AOOs are multiples of 4 square kilometres. Some figures are
# incompatible with this definition. We solve this by rounding upwards:
aoo[, 2:5] <- ceiling(aoo[, 2:5] / 4) * 4

# (4) Is any AOO greater than the area of mainland Norway?
w <- which(aoo$high > 323800)
if (length(w)) {
  print(aoo[w,])
} else {
  cat("Everything is fine.\n")
}

# Define additional auxiliary functions:
# These two function use maximum-likelihood estimation to infer the standard deviation
# of the AOOs, based on the best estimate (median), low estimate (1st quartile) and 
# high estimate (3rd quartile), and assuming a log-normal distribution.
f <- function(s, mean, q1, q3) return(
  (q1 - qlnorm(0.25, log(mean) - exp(2*s)/2, exp(s)))^2 +
  (q3 - qlnorm(0.75, log(mean) - exp(2*s)/2, exp(s)))^2
)
findSD <- function(Ex, q1, q3) exp(optimise(f, c(-12, 12), mean=Ex, q1=q1, q3=q3)$min)


# Simulations
# This version creates the overall indicator value for P3(e)
N <- 100000  # random numbers per species
M <- 323800  # maximum possible area
AOO <- matrix(0, N, nrow(aoo))
for (i in 1:nrow(aoo)) {
  if (any(is.na(aoo[i, 3:5]))) {
    if (is.na(aoo$best[i])) {
      if (is.na(aoo$known[i])) {
        # if no AOO is provided, assume it is 0
        AOO[, i] <- 0
      } else {
        # of not total AOOs are provided, assume they equal the known AOO
        AOO[, i] <- aoo$known[i]
      }
    } else {
      # if no low and high estimates are provided, use the best estimate
      AOO[, i] <- aoo$best[i]
    }
  } else {
    if (aoo$best[i] == 0) {
      # is the best estimate is 0, keep it 0
      AOO[, i] <- 0
    } else {
      # if low, best and high estimates are provided, estimate the standard
      # deviation and generate log-normally distributed random numbers
      SD <- findSD(aoo$best[i], aoo$low[i], aoo$high[i])
      AOO[, i] <- qlnorm(runif(N), log(aoo$best[i]) - SD * SD / 2, SD)
      AOO[, i] <-  sapply(AOO[, i], min, M)  # constrain to ara of Norway
      AOO[, i] <- ceiling(AOO[, i] / 4) * 4  # ensure multiples of 4 km^2
    }
  }
}


{ # Output of the results
  cat("Pe(e) is " %+% round(mean(apply(AOO, 1, sum))) %+% " ± " %+%
                   round(  sd(apply(AOO, 1, sum))) %+% " (mean ± SD)\n\n")
  cat("Confidence levels:\n")
  print(round(quantile(apply(AOO, 1, sum), c(0.025, 0.25, 0.5, 0.75, 0.975))))
}


################################ AOO-ene

#png("c:\\nina\\fremmede\\tiltaksp.lan\\opsjon\\figP2b.png", 1000, 1000, res=180)
par(mai=c(0.96, 0.96, 0.06, 0.06), lwd=1.8)
H <- hist(log10(as.numeric(aoo$best)), breaks=seq(0, 6, 1/3),
  xlab="Area of occupancy (km²)", ylab="Number of species", main="",
  xaxt="n", cex.axis=1.2, cex.lab=1.8, lwd=1.8, col=grey(0.84))
axis(1, 0:6, c(1, 10, 100, expression(10^3), expression(10^4), expression(10^5),
  expression(10^6)), cex.axis=1.2, lwd=1.8)
#dev.off()





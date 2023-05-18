
# This code calculates indicator P2 of NorInAliS



# Load data from the Alien Species List 2018
# The data are available from https://doi.org/10.5061/dryad.8sf7m0cjc
fab  <- read.csv2(url("https://datadryad.org/stash/downloads/file_stream/359484"),
                  as.is=T)
# The AOOs (areas of occupancy) provided in the above dataset are only the _best_
# estimates of the _total_ AOO for each species (i.e. the _median_ expert judgement
# of the _real_ AOO, including "dark figures" or unreported occurrences). However,
# this indicator needs more than that, viz. the _low_ and _high_ estimates of the
# total AOO (i.e. lower and upper _quartiles_) as well as the _known_ AOO 
# (i.e. excluding dark figures). These values are here read from a separate file.
# Their source is https://artsdatabanken.no/fremmedartslista2018
aoo <- read.csv2("aoo.txt", as.is=T)


# Restrict data to alien species that are reproducing unaidedly in mainland Norway
fab <- fab[which(fab$Status == "reproducing" & fab$Mainl),]


# Restrict the data to terrestrial species
w <- which(fab$LifeSt %in%
           c("lim", "lim,mar", "lim,par", "lim,mar,par", "mar", "mar,par"))
fab <- fab[-w,]
aoo <- aoo[-w,]


# Make sure that the two data frames are compatible
if (all(fab$Name == aoo$Name)) {
  cat("Everything is fine.\n")
} else {
  cat("ERROR: For some reason, the two dataframes are not compatible!\n")
}


# Load auxiliary functions
eval(parse(text=readLines("function.r")))


# Define additional functions:
# These two function use maximum-likelihood estimation to infer the std. deviation
# of the AOOs, based on the best estimate (median), low estimate (1st quartile)  
# and high estimate (3rd quartile), and assuming a log-normal distribution.
f <- function(s, mean, q1, q3) return(
  (q1 - qlnorm(0.25, log(mean) - exp(2*s)/2, exp(s)))^2 +
  (q3 - qlnorm(0.75, log(mean) - exp(2*s)/2, exp(s)))^2
)
findSD <- function(Ex, q1, q3) return(
  exp(optimise(f, c(-12, 12), mean=Ex, q1=q1, q3=q3)$min)
)


# Check for obvious errors in the original data:
# (1) Is any low estimate greater than the corresponding best estimate? 
w <- which(aoo$low > aoo$best)
if (length(w)) {
  print(aoo[w,])
}
# This looks very much like a punching error in line 28. Most likely, the best
# estimate should have been in the middle between the low and high estimate
# (implying a dark figure of 1000 rather than 10). The error is corrected manually:
aoo$best[w] <- 80000


# (2) Is any high estimate less than the corresponding best estimate?
w <- which(aoo$high < aoo$best)
if (length(w)) {
  print(aoo[w,])
}
# Here, it seems that the low and high estimate simply haven't been provided,
# and that, as a default, the known AOO was listed instead. In the absence of
# other information, we have to assume the low and high estimates to be equal
# to the best estimate:
aoo$low[w] <- aoo$high[w] <- aoo$best[w]


# (3) By definition, AOOs are multiples of 4 square kilometres. Some figures
# are incompatible with this definition. We solve this by rounding upwards:
aoo[, 2:5] <- ceiling(aoo[, 2:5] / 4) * 4


# (4) Is any AOO greater than the area of mainland Norway?
w <- which(aoo$high > 323800)
if (length(w)) {
  print(aoo[w,])
}
# Indeed. That is clearly an overestimate. The high estimate is reduced to
# (a bit less than) the area of Norway. All other estimates are possible 
# in principle, and are thus left unchanged.
aoo$high[w] <- 323000


# Simulations
# This version creates the overall indicator value for P2.
# For high- and severe-impact species, do this first:
# aoo <- aoo[which(fab$Impact %in% c("HI", "SE")),]
# For the remaining species, do this first:
# aoo <- aoo[which(fab$Impact %in% c("NK", "LO", "PH")),]
N <- 100000  # random numbers per species
M <- 323800  # maximum possible area in km^2
AOO <- matrix(0, N, nrow(aoo))
for (i in 1:nrow(aoo)) {
  if (any(is.na(aoo[i, 3:5]))) {
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
}


{ # Output of the results
  cat("P2 is " %+% round(mean(apply(AOO, 1, sum)), -2) %+% " ± " %+%
                   round(  sd(apply(AOO, 1, sum)), -2) %+% " (mean ± SD)\n\n")
  cat("Confidence levels:\n")
  print(round(quantile(apply(AOO, 1, sum), c(0.025, 0.25, 0.5, 0.75, 0.975)), -2))
}


# Show a graph of the distribution of the best estimates of total AOOs
par(mai=c(0.96, 0.96, 0.06, 0.06), lwd=1.8)
H <- hist(log10(as.numeric(aoo$best)), breaks=seq(0, 6, 1/3),
  xlab="Area of occupancy (km²)", ylab="Number of species", main="",
  xaxt="n", cex.axis=1.2, cex.lab=1.8, lwd=1.8, col=grey(0.84))
axis(1, 0:6, c(1, 10, 100, expression(10^3), expression(10^4), expression(10^5),
  expression(10^6)), cex.axis=1.2, lwd=1.8)


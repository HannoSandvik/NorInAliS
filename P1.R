
# Load data from the Alien Species List 2018
fab  <- read.csv2("c:\\art\\publ.z\\ecosolev.20\\assess.txt", as.is=T) #¤

# Restrict data to alien species that are reproducing unaidedly in mainland Norway
fab <- fab[which(fab$Status == "reproducing" & fab$Mainl),]

# Load auxiliary functions
eval(parse(text=readLines("c:\\art\\neobindc\\function.r"))) #¤

# Define an additional auxiliary function
# This one assign numerical values from 1 to 5 to the five ecological impact catgeories
as.nr <- function(x) sapply(x, function(z) switch(z, "SE"=5,"HI"=4,"PH"=3,"LO"=2,"NK"=1,NA))

{ # P1(a) is simply the number of reproducing alien species in mainland Norway
  P1a <- length(which(fab$Impact != "NR"))
  cat("The total number of alien species recorded as reproducing unaidedly in Norway: " %+%
      P1a %+% "\n\n")
}

{ # The disaggretation for ecological impact categories looks like this:
  cat("Alien species in Norway by impact category:\n")
  print(table(fab$Impact)[c("NK", "LO", "PH", "HI", "SE")])
}

{ # P1(b) is a weighted sum of impact categories
  # It takes the uncertainty in impact categories into account, quantifiable as quartiles:
  q0 <-       sum(as.nr(fab$minImp))  # minimum and 2.5% confidence level
  q1 <- round(sum(as.nr(fab$Impact) * 0.75 + as.nr(fab$minImp) * 0.25))
  q2 <-       sum(as.nr(fab$Impact))  # best estimate (mean and median)
  q3 <- round(sum(as.nr(fab$Impact) * 0.75 + as.nr(fab$maxImp) * 0.25))
  q4 <-       sum(as.nr(fab$maxImp))  # maximum and 97.5% confidence level
  SD <- round(sd(c(rep(q0,   floor(S1a / 4)),
                   rep(q2, ceiling(S1a / 2)),
                   rep(q4,   floor(S1a / 4))))) # approximation of the S.D.
  conf <- c(q0, q1, q2, q3, q4)
  names(conf) <- c("2.5", "25", "50", "75", "97.5") %+% "% CL"
  cat("P1(b) is " %+% q2 %+% " ± " %+% SD %+% " (mean ± SD)\n\n")
  cat("Confidence levels:\n")
  print(conf)
}


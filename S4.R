
# This code calculates indicator S4 of NorInAliS



# Load auxiliary functions
eval(parse(text=readLines("function.r")))


# Input the data for 2016 to 2020
# For the time being, this has to be done manually.
# Source: https://view.nina.no/planteimport/
# Data were read from the graph at the tab "Cumulative graphs"
# with the following parameter choices:
# - one year at a time (2016 to 2020)
# - one taxon at a time (vascular plants and arthropods)
# - "Imported item" = "All"
# - "Export country" = "All"
# - "Data to plot" = "Taxon
# - "Exclude juveniles" was _not_ ticked
# - "Show only alien species" _was_ ticked
Arthropods <- matrix(
  c(58, 38, 33, 23, 21,
    17, 12, 14, 11, 12),
  2, 5, TRUE,
  list(c("species", "containers"), "yr" %+% 2016:2020))
Plants <- matrix(
  c(33, 26, 22, 16, 11,
    17, 12, 14, 09, 10),
  2, 5, TRUE,
  list(c("species", "containers"), "yr" %+% 2016:2020))


{ # Output indicator values
  # Results are standardised for the contents of 10 containers,
  # but are expressed per container
  pla <- mean(Plants[1,] / Plants[2,])
  art <- mean(Arthropods[1,] / Arthropods[2,])
  cat("Plants:\n"     %+%
       pla %+% " ± "  %+% (sqrt(pla * 10) / 10) %+% " (mean ± SD)\n\n")
  cat("Arthropods:\n" %+%
       art %+% " ± "  %+% (sqrt(art * 10) / 10) %+% " (mean ± SD)\n\n")
  cat("Totals:\n"     %+% (pla + art) %+% " ± " %+%
                  (sqrt((art + pla) * 10) / 10) %+% " (mean ± SD)\n\n")
  cat("Confidence levels:\n")
  conf <- qpois(c(0.025, 0.25, 0.5, 0.75, 0.975), (art + pla) * 10) / 10
  names(conf) <- c("2.5", "25", "50", "75", "97.5") %+% "% CL"
  print(conf)
}


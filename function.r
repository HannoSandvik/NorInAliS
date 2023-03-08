### Auxiliary functions

# Combines text strings
"%+%" <- function(string1, string2) paste0(string1, string2)

# Removes an element of a set
"%-%" <- function(arg1, arg2) sort(unique(arg1[which(!(arg1 %in% na.omit(arg2)))]))

# Calculates the union of two sets (vectors)
"%A%" <- function(set1, set2)
  if (is.null(set1)) logical(0) else as.vector(na.omit(set1[set1 %in% set2]))

# I just find this more intuitive...
"%contains%" <- function (textstring, searchtext) grepl(searchtext, textstring)

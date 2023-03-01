### Auxiliary functions

# Combines text strings
"%+%" <- function(string1, string2) paste0(string1, string2)

# Removes an element of a set
"%-%" <- function(arg1, arg2) sort(unique(arg1[which(!(arg1 %in% na.omit(arg2)))]))

# I just find this more intuitive...
"%contains%" <- function (textstring, searchtext) grepl(searchtext, textstring)


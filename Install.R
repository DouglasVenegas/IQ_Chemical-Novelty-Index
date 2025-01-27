install.packages(c(
  "repr", "IRdisplay", "evaluate", "crayon", "pbdZMQ", "devtools", "uuid", 
  "digest", "RColorBrewer", "tidyverse", "readxl", "rvest", "dplyr", 
  "tidyr", "igraph", "visNetwork", "readr", "ggplot2", "ggraph", 
  "graphTweets", "writexl", "httr", "proxy"
))
install.packages("devtools")
devtools::install_github("IRkernel/IRkernel")
IRkernel::installspec(user = FALSE)

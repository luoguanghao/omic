library(tidyverse)



get_rda <- function(rda){
    l <- load(rda)
    data <- eval(parse(text = l))
    return(data)
}




















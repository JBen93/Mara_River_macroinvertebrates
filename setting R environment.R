install.packages("renv")
renv::init()
install.packages("tidyverse")

renv::snapshot()
renv::restore()
summary(tidyverse)

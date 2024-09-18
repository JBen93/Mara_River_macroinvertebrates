install.packages("renv")
renv::init()
install.packages("tidyverse")
install.packages("psych")
install.packages("vegan")

renv::snapshot()
renv::restore()
summary(tidyverse)

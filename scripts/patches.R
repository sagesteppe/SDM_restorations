library(terra)
library(tidyverse)
r <- rast(nrows=18, ncols=36, xmin=0)
r[1:2, 5:8] <- 1.1
r[5:8, 2:6] <- 1.2
r[7:12, 22:36] <- 1.3
r[15:16, 18:29] <- 1.4
p <- patches(r)

plot(p)

tibble(
  Patch = values(p)[!is.na(values(p))],
  Suitability = values(r)[!is.na(values(p))]
) %>% 
  sample_n(15)


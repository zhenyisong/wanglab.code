library(tidyverse)

ino80.exprs.df <- data.frame( reads   = c(269, 193, 1451, 1224),
                              spieces = c('human','human','rat','rat'),
                              treat   = c('OE','GFP','OE','GFP'))

ggplot(ino80.exprs.df, aes(x = treat, y = reads, group = spieces, col = spieces )) + 
      geom_path(size = 2)
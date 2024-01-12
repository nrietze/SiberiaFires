library(tidyverse)
library(cowplot)

cci_data <- read.csv2('C:/Users/nrietze/Documents/1_PhD/8_CHAPTER2/output/cavm_vs_burned_area/Arctic_BA_FireCCI51.csv',
                      dec = '.')

cci_data$total <- rowSums(cci_data[,2:8])

ggplot(cci_data) +
  geom_bar(aes(x = Year, y = total), stat = 'identity') +
  xlab("") +
  ylab("Burned area in tundra (ha)") +
  theme_cowplot()
ggsave('C:/Users/nrietze/Documents/1_PhD/8_CHAPTER2/output/cavm_vs_burned_area/total_burned_area.png',
       width = 10, height = 5, dpi = 200)

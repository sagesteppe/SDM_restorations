setwd('/media/steppe/hdd/SDM_restorations/scripts')

library(hrbrthemes)
library(GGally)
library(viridis)
library(tidyverse)

p2dat <- '../results/summary'
f <- file.path(p2dat, list.files(p2dat)[ grep('1k.csv', list.files(p2dat))])

dat <- do.call("rbind", sapply(f, read.csv, simplify = FALSE)) |>
  rownames_to_column('Taxon') |>
  mutate(Taxon = gsub('1k.*', '', basename(Taxon))) |>
  filter(Metric %in% c('Independent_Variables', 'AUC', 'Sensitivity', 'Specificity',
                       'Kappa', 'Balanced Accuracy')) |>
  pivot_wider(names_from = Metric, values_from = Value) |>
  mutate(AUC_copy = AUC)

rm(p2dat, f)


# Plot
ggparcoord(dat,
           columns = c( 4:5, 2, 6:7), groupColumn = 8,
           showPoints = TRUE, scale = "globalminmax",
           title = "Common evaluation metrics for each model",
           alphaLines = 0.3 
) + 
  labs(y = 'Value', x = 'Measurement') + 
  scale_color_gradient(high = '#426b69', low = '#8BB174') +

  ylim(0, 1) + 
  theme(
    aspect.ratio = 4/3,
    text = element_text(colour = "white"),
    axis.text.x=element_text(colour="white"),
    axis.text.y=element_text(colour="white"),
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

ggsave('../plots/Eval_metrics.png', bg = 'transparent')

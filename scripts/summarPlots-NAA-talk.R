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
           title = "Evaluation metrics for models",
           alphaLines = 0.3 
) + 
  labs(y = 'Value', x = 'Metric') + 
  scale_color_gradient(high = '#426b69', low = '#8BB174') +
  theme(
    aspect.ratio = 4/3,
  #  text = element_text(colour = "white"),
   # axis.text.x=element_text(colour="white"),
  #  axis.text.y=element_text(colour="white"),
    axis.text.x = element_text(angle = 25, vjust = 0.65),
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

ggsave('../plots/Eval_metrics.png', bg = 'transparent')



#################################################################################
####   TABLE OF HOW MANY RECORDS WERE USED AND WHAT THEIR SOURCES WERE    ######

source('functions.R')

base_p <- '../data/raw/occurrence'

targets <- read.csv('../data/raw/2024sos_species.csv', na.strings = "") %>% 
  ## we will not model a few species which are widely abundant, hence do not require modelling, and would utilize enormous computational resources.
  filter(!species %in% c('tridentata', 'millefolium', 'dioica', 'negundo')) %>%  # purshia,
  unite('taxon', genus:infraspecies, na.rm = T, sep = ' ') %>% # artemisia, and larrea !
  mutate(taxon = str_squish(taxon)) |>
  pull(taxon) |>
  unique()

  
bien <- read.csv(file.path(base_p, 'BIEN_records.csv')) %>% 
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) %>% 
  mutate(date = as.Date(date_collected),
         datasource = case_when(
           datasource == 'VegBank' ~ 'VegBank', 
           datasource == 'FIA' ~ 'FIA', 
           .default = as.character('Herbarium')
         )) |>
  select(taxon = scrubbed_species_binomial, date, datasource)


idig <- read.csv(file.path(base_p, 'IDIG.csv')) %>% 
  st_as_sf(coords = c('geopoint.lon', 'geopoint.lat'), crs = 4326) %>% 
  mutate(date = as.Date(datecollected), datasource = 'Herbarium') %>% 
  select(taxon = scientificname, date, datasource)

aim <- st_read(file.path(base_p, 'AIM', 'AIM_records.shp'), quiet = TRUE) %>% 
  mutate(date = as.Date(DatVstd)) %>% 
  select(taxon, date, PrmryKy) %>% 
  st_transform(4326)

gbif1 <- read.delim(file.path(base_p, 'gbif1', 'occurrence.txt'), na.strings = "") %>% 
  drop_na(decimalLongitude, decimalLatitude) %>% 
  filter(str_detect(decimalLatitude, '[A-Z]|[a-z]', negate = T)) %>% 
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>% 
  mutate(date = as.Date(eventDate), .before = geometry, 
         coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
         datasource = if_else(publisher == 'iNaturalist.org', 'iNaturalist', 'Herbarium') ) %>% 
  filter(coordinateUncertaintyInMeters <= 250) %>% 
  select(taxon = species, date, datasource)

gbif2 <- read.delim(file.path(base_p, 'gbif2', 'occurrence.txt'), na.strings = '')%>% 
  drop_na(decimalLongitude, decimalLatitude) %>% 
  filter(str_detect(decimalLatitude, '[A-Z]|[a-z]', negate = T)) %>% 
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>% 
  mutate(date = as.Date(eventDate), .before = geometry, 
         coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
         datasource = if_else(publisher == 'iNaturalist.org', 'iNaturalist', 'Herbarium') ) %>% 
  filter(coordinateUncertaintyInMeters <= 250) %>% 
  select(taxon = species, date, datasource)

records <- bind_rows(aim, bien, gbif1, gbif2, idig) %>% 
  mutate(taxon = str_to_sentence(taxon), 
         taxon = str_squish(taxon)) %>% 
  filter(taxon %in% targets) %>% 
  arrange(taxon) # recs need be > 90 m from each other

rm(aim, bien, gbif1, gbif2, idig, base_p, targets)

splicies <- split(records, records$taxon)

thinned <- lapply(splicies, ONEperCELL)
thinned <- bind_rows(thinned)

thinned <- mutate(thinned, datasource = if_else( !is.na(PrmryKy), 'AIM', datasource ))|>
  select(-PrmryKy)

any(is.na(thinned$datasource))
table(thinned$datasource)



################################################################################
#########    WHICH VARIABLES HAD MOST OF THE INFLUENCE ON THE MODELS?   ########
returnGini <- function(x){
  
  f <- paste0(x, list.files(x, pattern = '1k-'))
  mods <- lapply(f, readRDS)
  mods <- parallel::mclapply(mods, randomForest::importance)
  
  names(mods) <- gsub('1k-.*$', '', basename(f))
  mods <- purrr::map(mods, data.frame) |>
    purrr::map(tibble::rownames_to_column, var = 'Variable') |>
    dplyr::bind_rows(.id = 'Taxon')
  
  return(mods)
}

giniDecrease <- returnGini(x = '../results/rf_models/')

gd_summ <- giniDecrease |>
  group_by(Taxon) |>
  mutate(
    sumGini = sum(MeanDecreaseGini),
    propDecreaseGini = MeanDecreaseGini/sumGini
  ) |>
  ungroup() |>
  group_by(Variable) |>
  summarize(
    n = n(), 
    mean = mean(propDecreaseGini)
  ) |>
  arrange(-n) |>
  mutate(
    Theme = case_when(
      Variable %in% c(
        'bio1', 'bio10', 'bio11', 'bio12', 'bio18', 'bio19', 'bio4', 'bio5',
        'bio6', 'gdd5', 'gdgfgd5', 'hurs', 'ngd5', 'scd', 'swe', 'vpd') ~ 'Climate',
      Variable %in% c('cfvo_30_60', 'depth2bedrock', 'phh2o_30_60', 'salinity', 
                      'sand_5_15', 'silt_5_15', 'soc_0_5') ~ 'Soil',
      Variable %in% c('aspect', 'cti', 'dem', 'slope', 'tpi', 'tri') ~ 'Topography',
      Variable %in% c('herb', 'lastwild', 'shrub', 'tree') ~ 'Land Cover/Landuse'
    )
  )


model <- lm(gd_summ$mean ~ gd_summ$n)

ggplot(data = gd_summ, aes(x = n, y = mean, color = Theme, label = Variable)) + 
  geom_smooth(method = 'glm', formula = 'y ~ x', aes(x = n, y = mean),
              inherit.aes = F, color = 'black') +
  annotate(
    geom = 'text', x = 300, y = 0.025, parse = TRUE, 
    label = paste("~R^2==~",  round(broom::glance(model)[,2], 3)))  + 
  scale_color_manual(
    breaks = c("Climate", "Soil", "Land Cover/Landuse", "Topography"),
    values=c("#773344", "#20A4F3", "#63A375", "#EA7317")
  ) + 
  geom_jitter() + 
  ggrepel::geom_label_repel(max.overlaps = 20, alpha = 0.9) + 
  labs(
    title = 'Relationship between the number of models featuring\n a variable and the variables mean contribution to all models', 
    x = 'Number of models including variable', 
    y = ' Proportional mean decrease in Gini', 
    caption = expression(~italic('Variables become more important as they increase the mean decrease in Gini'))
  ) + 
  theme_classic() + 
  theme(
    axis.line.y = element_blank(), 
    axis.line.x = element_blank(),
    plot.title = element_text(hjust = 0.5), 
    legend.position = 'bottom'
  ) 


rm(ob, giniDecrease, gd_summ, model, returnGini)

ggsave('../plots/VariableImportance.png', bg = 'transparent', 
       width = 8, height = 4, units = 'in')

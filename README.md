## Variables and Sources for Environmental Niche Models

Two major models are created, one which focuses on identifying areas of fundamental environmental niches in the wild, and a second which focuses on identifying farms sites areas which the species has considerable niche overlap with. 


| Layer |                       Description                       |              Source                   |     Model    |      
| :---: | :-----------------------------------------------------: | :-----------------------------------: | :----------: | 
|  1.   |              Mean Annual Air Temperature (BIO1)         |        PRISM / CLIMATENA/ DISMO       |     both     |
|  2.   |               Temperature seasonality (BIO4)            |        PRISM / CLIMATENA / DISMO      |     both     |
|  3.   |         Max Temperature of Warmest Month (BIO5)         |        PRISM / CLIMATENA / DISMO      |     both     |
|  4.   |         Min Temperature of Coldest Month (BIO6)         |        PRISM / CLIMATENA / DISMO      |     both     |
|  5.   |        Mean Temperature of Warmest Quarter (BIO10)      |        PRISM / CLIMATENA / DISMO      |     both     |
|  6.   |        Mean Temperature of Coldest Quarter (BIO11)      |        PRISM / CLIMATENA / DISMO      |     both     |
|  7.   |              Mean annual precipitation (BIO12)          |        PRISM / CLIMATENA / DISMO      |     both     |
|  8.   |         Precipitation of Warmest Quarter (BIO18)        |        PRISM / CLIMATENA / DISMO      |     both     |
|  9.   |        Precipitation of Colest Quarter (BIO19)          |        PRISM / CLIMATENA / DISMO      |     both     |
| 10.   |                Mean annual cloudiness - MODIS           |          Wilson et al. 2016           |     both     |
| 12.   |         Beginning of the frost-free period (gdgfgd0)    |              Wang et al.              |     both     |
| 12.   |                   Climatic moisture deficit             |              Wang et al.              |     both     |
| 13.   |                  Degree-days above 5C (gdd5)            |              Wang et al.              |     both     |
| 14.   |               Mean annual precipitation as snow         |              Wang et al.              |     both     |
| 15.   |                 Percent Herbaceous Vegetation           |               EarthEnv                |     wild     |
| 16.   |                     Percent Shrub Cover                 |               EarthEnv                |     wild     |
| 17.   |                      Percent Tree Cover                 |               EarthEnv                |     wild     |
| 18.   |          Soil probability of bedrock (R Horizon)        |              SoilGrids                |     wild     |
| 19.   |                Soil organic carbon (Tonnes / ha)        |              SoilGrids                |     wild     |
| 20.   |                Surface (0-5 cm) soil pH in H~2~O        |              SoilGrids                |     both     |  
| 21.   |                   30-60 cm soil pH in H~2~O             |              SoilGrids                |     both     |
| 22.   |                Surface (0-5 cm) soil % sand             |              SoilGrids                |     both     |
| 23.   |                    5-15 cm  soil % sand                 |              SoilGrids                |     both     |
| 24.   |                    15-30 cm soil % sand                 |              SoilGrids                |     both     |
| 25.   |                Surface (0-5 cm) soil % clay             |              SoilGrids                |     both     |
| 26.   |                    5-15 cm  soil % clay                 |              SoilGrids                |     both     |
| 27.   |                    15-30 cm soil % clay                 |              SoilGrids                |     both     |
| 28.   |              Surface (0-5 cm) coarse fragments          |              SoilGrids                |     wild     | 
| 29.   |                      Soil USDA class                    |              SoilGrids                |     both     |
| 30.   |                         Elevation                       |             Geomorpho90               |     both     |
| 31.   |                          Slope                          |             Geomorpho90               |     wild     |
| 32.   |                          Aspect                         |             Geomorpho90               |     wild     |
| 33.   |                 Topographic Wetness Index               |             Geomorpho90               |     wild     |
| 34.   |                  Terrain ruggedness Index               |             Geomorpho90               |     wild     |
| 35.   |                       Geomorphon                        |             Geomoprho90               |     wild     |
| 36.   |        Estimated actual (w/-cloud) solar radiation      |      r.sun / Wilson et al. 2016       |     both     |
| 37.   |        Log-transformed distance to surface water        |     Global Surface Water Explorer     |     both     |


Some workflows for dealing with temporal mismatch between environmental variables and species occurrence records. 

In order to reduce the number of records which have questionable geographic integrity, no records before 1950 were acquired. 

1) Reducing Historic Records  
a) remove records with noted accuracy > 250 m  
b) if multiple records exist per raster cell, reduce to most recent record.
c) remove records on county centroids  
d) flag all records which are located within a burn scar, which occurred after the date of observation.   
d1) for that species, use plot based occurrence data, and years since burn interacting with cover of invasive annual grasses to model probability of an extant population at the site.   
d2) randomly sample flagged records with probabilities defined in d1.   


Arid and semi-arid lands constitute the plurality of earths terrestrial surface. 
These areas have been used both historically and currently for raising animals, such as cattle and sheep, without wholesale conversion of the land to agronomic settings, i.e. pastoral or rangelands.
These agricultural related uses of land have often lead to degradation of these areas.
Historic actions and phenomena such as fire suppression, shrub encroachment, over grazing, in combination with climate change and droughts are resulting in wildfires and there subsequent burns scars which are often dominated by a handful of invasive species.
These invasive species decrease the ability of the land to serve the natural wildlife, and human uses of the land, and in extreme cases preclude all uses of the land. 

In order to prevent the degradation of burned areas into stands of invasive species, restoration plantings have long been used. 
More recently a focus on utilizing native species, and now locally adapted genotypes, has been adopted as these species and populations have been shown to provide ecosystem services to a variety of wildlife and stakeholders, and they are often more capable of establishing at sites. 
However, the acquisition of native plant materials in the volumes required for restoration, especially in fire prone years, is not yet feasible. 
In an attempt to provide for the demand of locally adapted native plant materials partnerships with farmers who are planting wild harvested seeds, and harvesting their seed for sale to restoration performing agencies has increased. 

This solution still faces three challenging steps. 
1) Locating populations of adequate size that they can have enough viable seed collected from that that it is  economically viable for farmers to grow out from the seed lots while not adversely affecting the demographic vital rates of the population. 
2) Minimizing the adverse affects of artificial selection of these seed lines while being increased in agricultural settings. 
3) Ensuring that farmed seeds are applied to the areas in burn scars which they have the highest probability of success. 

The United States Bureau of Land Management (BLM) is currently addressing solutions to these three challenges. 




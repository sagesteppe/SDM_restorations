## Variables and Sources for Environmental Niche Models

Two major models are created, one which focuses on identifying areas of fundamental environmental niches in the wild, and a second which focuses on identifying farms which the species has considerable niche overlap with. The full 'wild' models feature up to 38 variables, while the 'farm' models features ca. 26 variables. Both endeavors run step wise. 


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
|  9.   |        Precipitation of Coldest Quarter (BIO19)         |        PRISM / CLIMATENA / DISMO      |     both     |
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
| 38.   |                   Human Influence Index                 | TIGRIS/ BTS NARN/ HDX/ NLCD/ NASA     |     wild     |


## Workflow for Presence/Absence generation and cleaning  

Some considerations for dealing with temporal mismatch between environmental variables and species occurrence records. 

In order to reduce the number of records which have questionable geographic integrity, no records before 1950 were acquired. 

### Thinning Historic Records  
a) remove records with noted coordinate uncertainty of > 270m  
b) remove records on county centroids  
c) if multiple records exist per raster cell, reduce to that with the lowest coordinate uncertainty, or the most recent

### Processing Post Burn Records
a) flag all historic records which are located within a burn scar, which occurred after the date of observation.    
b) for that species, use plot based occurrence data, and years since burn interacting with cover of invasive annual grasses to model probability of an extant population at the site.  $\text{occurrence ~ year since fire * modelled invasive annual grass cover}$    
c) randomly sample flagged records with probabilities defined in *b*  

### Processing Records in Sites Invaded by Invasive Annual Grasses
a) flag all historic records which are located within an area with greater >10% invasive annual grass cover    
b) for that species, use plot based occurrence data, and cover of invasive annual grasses to model probability of an extant population at the site.  $\text{occurrence ~ modelled invasive annual grass cover}$     
c) randomly sample flagged records with probabilities defined in *b*  

### Identifying Infraspecies Records

Many collectors of field data do not match taxa to infraspecies due to sampling at times when the population lacks the appropriate character states. Some of these records may be recoverable by identifying their 1) nearest neighbors 2) whethery they lay within the n-dimensional convex hulls of hypervolumes generated from unambigiously identified records.

## A Better Test Set of Ensembled Predictions



## Better matching the Fundamental Niche to Population Occurrence via Cost Surfaces  

Cost surfaces reflect historic connectivity of locations prior to European settlement.

|  Layer  |                       Description                       |              Source                   |    Effect      |
| :-----: | :-----------------------------------------------------: | :-----------------------------------: | :------------: |
|   1.    |                 Terrain Ruggedness Index                |             Geomorpho90               |    negative    |
|   2.    |        Log-transformed distance to surface water        |     Global Surface Water Explorer     |    negative    |
|   3.    |      Log-transformed distance to perennial streams      |     National Hydrography Dataset      |    negative    |
|   4.    |                      Percent Barren                     |              EarthEnv                 |    negative    |
|   5.    |          Max Temperature of Warmest Month (BIO5)        |        PRISM / CLIMATENA / DISMO      |    negative    |
|   6.    |                     Wind Direction*                     |              rWind/GFS                |   anisotropic  |
|   7.    |                       Wind Speed*                       |              rWind/GFS                |    positive    |

Variables indicated by a '*' exist only for taxa which have evident adaptions to wind dispersal, notably those of the Asteraceae.

## Correlation between patch habitat suitable and abundance

Assess, Inventory, and Monitor data can be used to test the relationship between modelled habitat suitability and field measured abundance. Rather than simple relating the value of modelled suitability at a single cell, we can sum the modelled suitability of multiple cells, which together theoretically form a population. Patches may be identified using terras 'get_patches' function. 

Thereafter, 
The suitability of any discrete patch of cells can be defined as: 

$$PS =  \sum_{\text{i = 1}}^{n} {x_{i}}$$
Where *n* represents each cell in a patch. Hence it is a summation of total suitability across the patch to create *PS*, patch suitability. I anticipate this value best be transformed to reduce the variation between large and small events, by either $\log(PS)$ or $\sqrt{PS}$. Subsequently, each species should be tested as an effect, wherein under the maximal model they are treated as a fixed effect. 

These steps would develop this set of candidate models

<center>
  
% Cover ~ PS | species   
% Cover ~ sqrt(PS) | species   
% Cover ~ log(PS) | species   

% Cover ~ PS * species   
% Cover ~ sqrt(PS) * species   
% Cover ~ log(PS) * species   

% Cover ~ PS + species   
% Cover ~ sqrt(PS) + species   
% Cover ~ log(PS) + species   

</center> 

The model with the highest correlation coefficient to the occurrence data would be the final model indicating the strength of correlation between the SDM's and in field measured occurrence.

In instances where multiple field observations exist per patch, patches will be decomposed so that the cells closest to the plots become the unit at which PS is calculated.

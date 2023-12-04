## Variables and Sources for Environmental Niche Models

Two major models are created, one which focuses on identifying areas of fundamental environmental niches in the wild, and a second which focuses on identifying farms which the species has considerable niche overlap with. The full 'wild' models feature up to 38 variables, while the 'farm' models features ca. 26 variables. Both endeavors run step wise. 


| Layer |                       Description                       |              Source                   |     Model    |      
| :---: | :-----------------------------------------------------: | :-----------------------------------: | :----------: | 
|  1.   |              Mean Annual Air Temperature (BIO1)         |              Chelsa                   |     both     |
|  2.   |               Temperature seasonality (BIO4)            |              Chelsa                   |     both     |
|  3.   |         Max Temperature of Warmest Month (BIO5)         |              Chelsa                   |     both     |
|  4.   |         Min Temperature of Coldest Month (BIO6)         |              Chelsa                   |     both     |
|  5.   |        Mean Temperature of Warmest Quarter (BIO10)      |              Chelsa                   |     both     |
|  6.   |        Mean Temperature of Coldest Quarter (BIO11)      |              Chelsa                   |     both     |
|  7.   |              Mean annual precipitation (BIO12)          |              Chelsa                   |     both     |
|  8.   |         Precipitation of Warmest Quarter (BIO18)        |              Chelsa                   |     both     |
|  9.   |        Precipitation of Coldest Quarter (BIO19)         |              Chelsa                   |     both     |
| 10.   |         Beginning of the frost-free period (gdgfgd0)    |              Chelsa                   |     both     |
| 11.   |      Mean Monthly vapour pressure deficit (vpd_mean)    |              Chelsa                   |     both     |
| 12.   |     Heat accumulation of  Degree-days above 5C (gdd5)   |              Chelsa                   |     both     |
| 13.   |           Number of Degree-days above 5C (ngd5)         |              Chelsa                   |     both     |
| 14.   |   Heat accumulation of  Degree-days above 10C (gdd10)   |              Chelsa                   |     both     |
| 15.   |         Number of Degree-days above 10C (ngd10)         |              Chelsa                   |     both     |
| 16.   |           Mean Monthly cloudiness (tcc mean)            |              Chelsa                   |     both     |
| 17.   |        Mean monthly near surface humidity (hurs_mean)   |              Chelsa                   |     both     |
| 18.   |            Number of Days with Snow Cover(scd)          |              Chelsa                   |     both     |
| 19.   |            Annual Snow Water Equivalent (swe)           |              Chelsa                   |     both     |
| 20.   |          Number of days where temps cross 0C (fcf)      |              Chelsa                   |     both     |
| 21.   |           Mean annual precipitation as snow             |              Chelsa                   |     both     |
| 22.   |                Percent Herbaceous Vegetation            |              EarthEnv                 |     wild     |
| 23.   |                   Percent Shrub Cover                   |              EarthEnv                 |     wild     |
| 24.   |                   Percent Tree Cover                    |              EarthEnv                 |     wild     |
| 25.   |            Soil probability of bedrock (R Horizon)      |              SoilGrids                |     wild     |
| 26.   |                Soil organic carbon (Tonnes / ha)        |              SoilGrids                |     wild     |
| 27.   |                Surface (0-5 cm) soil pH in H~2~O        |              SoilGrids                |     both     |
| 28.   |                   30-60 cm soil pH in H~2~O             |              SoilGrids                |     both     |
| 29.   |                Surface (0-5 cm) soil % sand             |              SoilGrids                |     both     |
| 30.   |                    5-15 cm  soil % sand                 |              SoilGrids                |     both     |
| 31.   |                    15-30 cm soil % sand                 |              SoilGrids                |     both     |
| 32.   |                Surface (0-5 cm) soil % clay             |              SoilGrids                |     both     |
| 33.   |                    5-15 cm  soil % clay                 |              SoilGrids                |     both     |
| 34.   |                    15-30 cm soil % clay                 |              SoilGrids                |     both     |
| 35.   |              Surface (0-5 cm) coarse fragments          |              SoilGrids                |     wild     | 
| 36.   |                      Soil USDA class                    |              SoilGrids                |     both     |
| 37.   |                   Global Soil Salinity                  |              SoilGrids                |     both     |
| 38.   |                         Elevation                       |             Geomorpho90               |     both     |
| 39.   |                          Slope                          |             Geomorpho90               |     wild     |
| 40.   |                          Aspect                         |             Geomorpho90               |     wild     |
| 41.   |                 Topographic Wetness Index               |             Geomorpho90               |     wild     |
| 42.   |                  Terrain ruggedness Index               |             Geomorpho90               |     wild     |
| 43.   |                       Geomorphon                        |             Geomoprho90               |     wild     |
| 44.   |        Log-transformed distance to surface water        |     Global Surface Water Explorer     |     both     |
| 45.   |                   Human Influence Index                 |           NASA Earth Data             |     wild     |


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

Many collectors of field data do not match taxa to infraspecies due to sampling at times when the population lacks the appropriate character states. Some of these records may be recoverable by identifying their:
1) nearest neighbors  
2) whether they lay within the n-dimensional convex hulls of hyper volumes generated from unambiguously identified records.

## Modelling

### Ensemble learning methods

Boosted Regression Trees, and Random Forests models will be generated for each species. \
The Boruta algorithm will be applied to the full stack of predictor variables to remove un-informative predictors, to reduce the costs associated with overfitting, but largely in an effort to speed up projecting models onto gridded surfaces. 
The projection of models onto raster surfaces has been identified as the rate limiting step in performing species distribution modelling using my setup. 
Subsequent to the thinning of un-informative variables via the boruta algorithm, groups of variables with high degrees of multicollinearity (e.g. % soil textural components at varying depths), will be thinned by selection of a single variable per group which *appeared* to have the highest variable importance via boruta analysis. 

## A Better Test Set of Ensembled Predictions

In order to generate test sets which may most accurately evaluate ensembled predictions of SDM's stratitied sampling of a reduced set hypervolume is used to identify 100 occurrence and 100 absence records. 

## Better matching the Fundamental Niche to Population Occurrence via Cost Surfaces  

Discrepancies exist between the niche which a species is capable of establishing a poopulation in, and it's existence there. These mismatches are due to dispersal limitations. In order to identify areas which have fundamental niche, but may not have extent populations the distance between the nearest known population and suitable patches of habitat. 

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

% Cover ~ PS  
% Cover ~ sqrt(PS)  
% Cover ~ log(PS)  

</center> 

The model with the highest correlation coefficient to the occurrence data would be the final model indicating the strength of correlation between the SDM's and in field measured occurrence.

In instances where multiple field observations exist per patch, patches will be decomposed so that the cells closest to the plots become the unit at which PS is calculated.

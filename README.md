## Variables and Sources for Environmental Niche Models

Two major models are created, one which focuses on identifying areas of fundamental environmental niches in the wild, and a second which focuses on identifying farms which the species has considerable niche overlap with. The full 'wild' models feature up to 33 variables, while the 'farm' models features ca. 26 variables. Both endeavors run step wise. 


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
| 10.   |      Mean Monthly vapour pressure deficit (vpd)         |              Chelsa                   |     both     |
| 11.   |     Heat accumulation of  Degree-days above 5C (gdd5)   |              Chelsa                   |     both     |
| 12.   |           Number of Degree-days above 5C (ngd5)         |              Chelsa                   |     both     |
| 13.   |   Heat accumulation of  Degree-days above 10C (gdd10)   |              Chelsa                   |     both     |
| 14.   |         Number of Degree-days above 10C (ngd10)         |              Chelsa                   |     both     |
| 15.   |        Mean monthly near surface humidity (hurs)        |              Chelsa                   |     both     |
| 16.   |            Number of Days with Snow Cover(scd)          |              Chelsa                   |     both     |
| 17.   |            Annual Snow Water Equivalent (swe)           |              Chelsa                   |     both     |
| 18.   |                Percent Herbaceous Vegetation            |              EarthEnv                 |     wild     |
| 19.   |                   Percent Shrub Cover                   |              EarthEnv                 |     wild     |
| 20.   |                   Percent Tree Cover                    |              EarthEnv                 |     wild     |
| 21.   |                 Depth to Bedrock (R Horizon)            |              SoilGrids                |     wild     |
| 22.   |                Soil organic carbon (Tonnes / ha)        |              SoilGrids                |     wild     |
| 23.   |                Surface (0-5 cm) soil pH in H~2~O        |              SoilGrids                |     both     |
| 24.   |                   30-60 cm soil pH in H~2~O             |              SoilGrids                |     both     |
| 25.   |                Surface (0-5 cm) soil % sand             |              SoilGrids                |     both     |
| 26.   |                    5-15 cm  soil % sand                 |              SoilGrids                |     both     |
| 27.   |                    15-30 cm soil % sand                 |              SoilGrids                |     both     |
| 28.   |                Surface (0-5 cm) soil % clay             |              SoilGrids                |     both     |
| 29.   |                    5-15 cm  soil % clay                 |              SoilGrids                |     both     |
| 30.   |                    15-30 cm soil % clay                 |              SoilGrids                |     both     |
| 31.   |              Surface (0-5 cm) coarse fragments          |              SoilGrids                |     wild     | 
| 32.   |                       Soil Salinity                     |              SoilGrids                |     both     |
| 33.   |                         Elevation                       |             MERIT - DEM               |     both     |
| 34.   |                          Slope                          |             Geomorpho90               |     wild     |
| 35.   |                          Aspect                         |             Geomorpho90               |     wild     |
| 36.   |                  Topographic Position Index             |             Geomorpho90               |     wild     |
| 37.   |        Log-transformed distance to surface water        |     Global Surface Water Explorer     |     both     |
| 38.   |                   Human Influence Index                 |           NASA Earth Data             |     wild     |


## Workflow for Presence/Absence generation and cleaning  

Some considerations for dealing with temporal mismatch between environmental variables and species occurrence records. 

In order to reduce the number of records which have questionable geographic integrity, no records before 1950 were acquired. 

### Thinning Historic Records  
a) remove records with noted coordinate uncertainty of > 270m  
b) if multiple records exist per raster cell, reduce to that with the lowest coordinate uncertainty, or the most recent

### Manual Review of Presence Records

### Acquisition & Generation of Absences

True absences are acquired from AIM plots... True absences are generated around a buffered convex hull of the species. 

### Pseudo-absence generation and reduction

PA's are generated within the range of the species... 

## Modelling

### Fitting Models

Random Forests models will be generated for each species. \
The Boruta algorithm will be applied to the full stack of predictor variables to remove un-informative predictors, to reduce the costs associated with overfitting, but largely in an effort to speed up projecting models onto gridded surfaces. 
The projection of models onto raster surfaces has been identified as the rate limiting step in performing species distribution modelling using my setup. 
Subsequent to the thinning of un-informative variables via the boruta algorithm, groups of variables with high degrees of multicollinearity (e.g. % soil textural components at varying depths), will be thinned by selection of a single variable per group which *appeared* to have the highest variable importance via boruta analysis. 
Following the dropping of collinear variables, Recursive Feature Elimination (rfe), is performed using the Caret package. 
RFE identifies variables which are uninformative using cross-fold validation and stepwise combinations of all predictors. 
The algorithm will suggest variables to remove, but given that the penalty for extraneous variables is low in random forest, it generally suggests utilizing all features. 
In order to better reduce the number of variables to speed up prediction of models onto rasters any number of predictors which result in a less than 1.5% decrease in mean accuracy, relative to the full model, are retained and the terms contributing to the model with the fewest variable are selected. 

## A Better Test Set of Ensembled Predictions



## Better matching the Fundamental Niche to Population Occurrence via Cost Surfaces  

Discrepancies exist between areas with suitable environmental niches which a species is capable of establishing a population in, and it's existence there, due most notably to dispersal limitation. 
In order to identify areas which have fundamental niche, but may not have extent populations, the distance between the nearest known population and suitable patches of habitat are analyzed using connectivity analysis. 
Cost surfaces are generated via the inversion of hypothesized habitat suitability maps, and circuit theory is used to score predicted suitable habitats from known populations. 


# this content goes somewhere.... 

### Processing Post Burn Records
a) flag all historic records which are located within a burn scar, which occurred after the date of observation.    
b) for that species, use plot based occurrence data, and years since burn interacting with cover of invasive annual grasses to model probability of an extant population at the site.  $\text{occurrence ~ year since fire * modelled invasive annual grass cover}$    
c) randomly sample flagged records with probabilities defined in *b*  

### Processing Records in Sites Invaded by Invasive Annual Grasses
a) flag all historic records which are located within an area with greater >10% invasive annual grass cover    
b) for that species, use plot based occurrence data, and cover of invasive annual grasses to model probability of an extant population at the site.  $\text{occurrence ~ modelled invasive annual grass cover}$     
c) randomly sample flagged records with probabilities defined in *b*  



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

# Introduction 

A primary challenge to restoring Earths terrestrial ecosystems is the lack of available plant germplasm (@national2023assessment, @merritt2011restoration).
Given the scale of our restoration needs, in most scenarios, the only sustainable source of seed is from grow-outs of wild harvested seed in agricultural settings (@pedrini2020collection, @broadhurst2015seeding,  @national2023assessment).
Enormous efforts are now underway to increase the number of species, the number of populations within these species, and the genetic diversity of these populations, available for restoration @national2023assessment.
However, difficulties exist in both the wild harvest and increase of seed which are limiting our ability to develop adequate amounts of germplasm.   

While most species desired in restorations have historically had relatively large geographic ranges, numbers of populations, and number of individuals per populations, the development of native germplasm remains behind targets (@national2023assessment). 
We posit that in part this is due to the difficulty of finding populations with the appropriate number of individuals, which are experiencing climatic conditions conducive to producing enough viable seed to being agricultural increase, a complications borne of widespread habitat degradation and unnatural wildfires *(cite on degradation and wildfires)*. 
Tools which are capable of predicting a species geographic range, the presence and size of populations across the range and in seed collection target units (such as empirical or provisional seed transfer zones or ecoregions), such as Species Distribution Models offer promise to increase the rate at which native germplasm can be developed.

However, while SDM's generate hypothesis of whether areas have environmental conditions similar to the observed environmental niche of the species they do not consider the probabilities of colonization, nor the populations census sizes.
A possible tool to associate a probability of occurrence of a species in a suitable habitat patch, is connectivity analysis.
Utilizing the predicted unsuitability of habitat, with extreme barriers to dispersal, to create cost-resistance surfaces between patches with known populations and predicted populations allows for simulating the probability of occurrence. 
While previous correlations between habitat suitability and population size have often been low, we posit that patch level parameters offer a more useful prediction of population size (@weber2017there, @waldock2022quantitative)
Patch metrics will more appropriately reflect parameters of a potential population, e.g. geographic extent, which we posit is a more informative proxy of population size than environmental niche correlation.  
Utilizing this approach will allow field crews to associate probabilities with previously unground truthed patches. 

Here we showcase the utility of SDM's to predict suitable habitat for common species in natural settings and use those data to test two specific hypotheses. 
Firstly that SDM's can identify more populations than exist in the occurrence based records they were derived from. 
And secondly patch metrics derived from SDM's correlate with observed population census size. 
In order to determine whether SDM's are useful in detecting new populations, over ~10~ field crews were given 50 putative populations to survey for the presence and absence of the modelled species. 
To determine if patch metrics can be used to predict population census sizes, these crews noted whether the populations was large enough to support 

# 2 MATERIALS AND METHODS  
## 2.1 Study system  

The 355 modelled species were selected by land managers in arid & semi-arid areas Western North America and are species they already or immediately intend to use in restoration. 
The spatial domain of modelling broadly encompasses the Western United States being bounded by the Pacific Ocean, the 50th parallel at north, -100 degrees at East and Mexico to the south. 
This domain has tremendous variation in amounts of elevation, temperature, and precipitation. 

## Species Occurrence Records  
Species presence records were collected from Natural History Museums, citizen science initiatives, and standardized ecological monitoring programmes using R (@chamberlain2024rgbif, @maitner2023bien, @michonneau2023idigbio). 
These records were filtered to only those collected after 1950, with coordinate uncertainty less than 250 meters, and only the most recent record per 90m cell was retained. 
All of these records were then manually reviewed by species, where records with more 2 or more variables in the 2.5% quantiles of several environmental variables (bio1, 4, 10, 12, & 19; TABLE X) and distance were flagged.
During subjective manual review all digitized herbarium specimens which were suspect were reviewed, while other records were dropped based on the analysts review of other occurrence records.   

Species absences were generated using three processes, and we sought to have a 1:1 ratio of presences to absences. 
True absences were acquired from a massive ecological monitoring program Assess, Inventory, and Monitor (AIM), these absences accounted for a percent equivalent to land managed by BLM within the species range (e.g. if 20% of land ownership across a species range was BLM, than the number of presences * 0.2 defined the number of these absences). 
Likely absences, representing 15% of records, were generated outside the known range of the species, but bounded to be within 50 miles of the extent extent of occurrences. 
Pseudo-absences (PA) were randomly selected from areas in, or within XX km, of the species range but greater than 10km from an occurrence, these records accounted for ((1 - (% BLM Land + 15 % LA records)) * 1.25).
Because most of these species are highly abundant in order to reduce the probability that a PA was drawn within a species the environmental niche linear discriminant analysis (LDA) was used. The LDA utilized the presences, true absences and likely absences as the classes and several environmental variables (@venables2002mass), identified by vifstep with theta = 10 as displaying at most moderate collinearity, as independent variables (@naimi2014usdm). 
As many records, classified by the LDA as originating from the Presence data set, from the pseudo-absences could be removed as to achieve a class balance between presences and absences. 

## Environmental Variables  
XX variables were used in generating the species distribution models (table X, @karger2017climatologies, @hengl2017soilgrids250m, @ivushkin2019global,  @tuanmu2014global, @yamazaki2017merit, @amatulli2020geomorpho90m, @sanderson2002human).
These variables were selected based on authors previous work, and represent variables commonly associated with the empirical distributions of taxa in arid and semi-arid regions. 
All variables were downloaded from source and re-sampled to 90m resolution using terra (@hijman2023terra).

## Species Distribution Modelling  
Random Forests models were generated for each species following several steps to reduce the number of features.
While superfluous features generally do not notably decrease the performance of random forests, they increase the amount of time required to predict models onto gridded surfaces. 
The Boruta algorithm was used on the full stack of 44?? independent variables to remove un-informative variables @kursa2010boruta), and the variable importance factor (VIF) scores were then used to subset the most informative features in several auto-correlated pairs, reducing the independent variables to at most 33.  
Recursive Feature Elimination (rfe) was then used to identify the fewest number of independent variables which could either increase model performance (measured using accuracy), or only decrease it by 1.5% relative to the full set of variables, these variables were subset and used as features for random forest modelling with all default values except for optimized mtrys (@kuhn2008caret, @liaw2002randomForest). 
Models were then predicted onto gridded surfaces which exceeded the species range by 50 miles (@terra). 

## Patch identification 
In order to identify patches which may support populations raster values less than 0.8 were masked as NA's and a rooks case search was conducted using @LANDSCAPEMETRICS. 
Multiple landscape metrics were calculated the Class metrics Euclidean Nearest Neighbor of both mean and coefficient of variation, and the Patch metrics: Euclidean Nearest Neighbor, Core Area Index, PARA, Fractal Index LANDSCAPE METRICS.
Because the > 0.8 threshold used above did not capture all colonized patches, all patches with predicted suitability scores of >0.55 which contained species observations were then identified using Terra's patches to create a data set of known populations. 

## Predicting Species Occurrence in Patches

All patches identified above were used to identify patches nearest neighbors. The 5 nearest neighbors of each patch were detected using @SPDEP. 


## Variables and Sources for Environmental Niche Models

Two major models are created, one which focuses on identifying areas of fundamental environmental niches in the wild, and a second which focuses on identifying farms which the species has considerable niche overlap with. 
The full 'wild' models feature up to 33 variables, while the 'farm' models features ca. 26 variables. Both endeavors run step wise. 


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
| 37.   |                   Human Influence Index                 |           NASA Earth Data             |     wild     |


data sources (@karger2017climatologies, @hengl2017soilgrids250m, @ivushkin2019global,  @tuanmu2014global, @yamazaki2017merit, @amatulli2020geomorpho90m, @sanderson2002human)


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



To determine whether habitat suitability could predict species abundance mean weighted patch suitability $PS =  \sum_{\text{i = 1}}^{n} {x_{i}}$ was regressed against abundance using GLM. 
All AIM plots were queried for all occasions which a target species was observed via either the Species Richness, or Line-Point Intercept method (150 points along three equiangular transects over a 25m radius outside a 5m radius buffer zone). 
Occurrence of a species only along SR was considered 1% cover, while the canopy cover from LPI was used. 

To determine whether connectivity could predict the presence of a population in a patch, connectivity was calculated between every patch, and jacknife-resampled for logistic regression where connectivity scaled from 0 to 100 served as a the independent variable. 


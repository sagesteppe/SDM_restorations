# Introduction 

A primary challenge to restoring Earths terrestrial ecosystems is the lack of available plant germplasm (@national2023assessment, @merritt2011restoration). 
Given the scale of our restoration needs, in most scenarios, the only sustainable source of seed is from grow-outs of wild harvested seed in agricultural settings (@pedrini2020collection, @broadhurst2015seeding,  @national2023assessment).
Enormous efforts are now underway to increase the number of species, the number of populations within these species, and the genetic diversity of the seed available for restoration @national2023assessment. 
However, numerous difficulties exist in both the wild harvest seed and it's increase which are limiting our ability to develop adequate amounts of germplasm.   

While most species desired in restorations have historically had relatively large geographic ranges, numbers of populations, and number of individuals per populations, the development of native germplasm remains behind targets (@national2023assessment). 
We posit that in part this is due to the difficulty of finding populations with the appropriate number of individuals, which are experiencing climatic conditions conducive to producing enough viable seed to being agricultural increase, a complications borne of widespread habitat degradation and unnatural wildfires *(@abatzoglou2021projected)*. 
Tools which are capable of predicting a species geographic range, the presence and size of populations across the range and in seed collection target units (such as empirical or provisional seed transfer zones or ecoregions), such as Species Distribution Models offer promise to increase the rate at which native germplasm can be developed.

However, while SDM's generate hypothesis of whether areas have environmental conditions similar to the observed environmental niche of the species they do not consider the probabilities of colonization, nor the populations census sizes. 
A possible tool to associate a probability of occurrence of a species in a suitable habitat patch, is connectivity analysis. 
Utilizing the predicted unsuitability of habitat, with extreme barriers to dispersal, to create cost-resistance surfaces between patches with known populations and predicted populations allows for simulating the probability of occurrence. 
While previous correlations between habitat suitability and population size have often been low, we posit that patch level parameters offer a more useful prediction of population size (@weber2017there, @waldock2022quantitative). 
Patch metrics will more appropriately reflect parameters of a potential population, e.g. geographic extent, which we posit is a more informative proxy of population size than environmental niche correlation.  
Utilizing this approach will allow field crews to associate probabilities with previously unground-truthed patches. 

Here we showcase the utility of SDM's to predict suitable habitat for common species in natural settings and use those data to test two specific hypotheses. 
Firstly that SDM's can identify more populations than exist in the occurrence based records they were derived from. 
And secondly patch metrics derived from SDM's correlate with observed population census size. 
In order to determine whether SDM's are useful in detecting new populations, over ~10~ field crews were given 50 putative populations to survey for the presence and absence of the modeled species. 
To determine if patch metrics can be used to predict population census sizes, these crews noted whether the populations was large enough to support 

# 2 MATERIALS AND METHODS  
## 2.1 Study system  

The 353 modelled species were selected by land managers in arid & semi-arid areas Western North America and are species they already or immediately intend to use in restoration. 
The spatial domain of modelling broadly encompasses the Western United States being bounded by the Pacific Ocean, the 50th parallel at north, -100 degrees at East and Mexico to the south. 
This domain has tremendous variation in amounts of elevation, temperature, and precipitation. 

## Species Occurrence Records  
Species presence records were collected from Natural History Museums, citizen science initiatives, and standardized ecological monitoring programmes using R (@chamberlain2024rgbif, @maitner2023bien, @michonneau2023idigbio). 
These records were filtered to only those collected after 1950, with coordinate uncertainty less than 250 meters, and only the most recent record per 90m cell was retained. 
All of these records were then manually reviewed by species, where records with more 2 or more variables in the 2.5% quantiles of several environmental variables (bio1, 4, 10, 12, & 19; see TABLE X) and distance were flagged.
During subjective manual review all digitized herbarium specimens which were suspect were reviewed, while other records were dropped based on the analysts review of other occurrence records.   

Species absences were generated using three processes, and we sought to have a 1:1 ratio of presences to absences. 
True absences were acquired from a massive ecological monitoring program Assess, Inventory, and Monitor (AIM), these absences accounted for a percent equivalent to land managed by BLM within the species range (e.g. if 20% of land ownership across a species range was BLM, than the number of presences * 0.2 defined the number of these absences). 
Likely absences, representing 15% of records, were generated outside the known range of the species, but bounded to be within 50 miles of the extent extent of occurrences. 
Pseudo-absences (PA) were randomly selected from areas in, or within XX km, of the species range but greater than 10km from an occurrence, these records accounted for ((1 - (% BLM Land + 15 % LA records)) * 1.25).
Because most of these species are highly abundant in order to reduce the probability that a PA was drawn within a species the environmental niche linear discriminant analysis (LDA) was used. The LDA utilized the presences, true absences and likely absences as the classes and several environmental variables (@venables2002mass), identified by vifstep with theta = 10 as displaying at most moderate collinearity, as independent variables (@naimi2014usdm). 
As many records, classified by the LDA as originating from the Presence data set, from the pseudo-absences could be removed as to achieve a class balance between presences and absences. 

## Environmental Variables  
Up to 44 variables were used in generating the species distribution models (table X, @karger2017climatologies, @hengl2017soilgrids250m, @ivushkin2019global,  @tuanmu2014global, @yamazaki2017merit, @amatulli2020geomorpho90m, @sanderson2002human). 
These variables were selected based on authors previous work, and represent variables commonly associated with the empirical distributions of taxa in arid and semi-arid regions. 
All variables were downloaded from source and re-sampled to 90m resolution using terra (@hijman2023terra). 

## Species Distribution Modelling  
Random Forests models were generated for each species following several steps to reduce the number of features.
While superfluous features generally do not notably decrease the performance of random forests, they increase the amount of time required to predict models onto gridded surfaces. 
The Boruta algorithm was used on the full stack of 44 independent variables to remove un-informative variables @kursa2010boruta), and the variable importance factor (VIF) scores were then used to subset the most informative features in several auto-correlated pairs, reducing the independent variables to at most 33.  
Recursive Feature Elimination (rfe) was then used to identify the fewest number of independent variables which could either increase model performance (measured using accuracy), or only decrease it by 1.5% relative to the full set of variables, these variables were subset and used as features for random forest modelling with all default values except for optimized mtrys (@kuhn2008caret, @liaw2002randomForest). 
Models were then predicted onto gridded surfaces which exceeded the species range by 50 miles (@terra). 

## Patch identification 
To identify putative metapopulations raster cells predicted as having less than 0.8 probability of suitable habitat were masked as NA's, and then areas which were crossed by streams with Strahler orders of three of more (@usgs2004nhd), and the divides of HU10 watersheds (@usgs2023wbd) were 'burnt' away from the raster. 
The resulting rasters were aggregated by a factor of 2 to 5, depending of their sizes (< 300 MiB 2, < 500 MiB 3, < 700 MiB 4, > 700 MiB 5) to accommodate, a rooks case search using terra in a practical time period. 
Resultant patches < 5 acres were then discarded, and the rasters was resampled to it's input resolution.  

Putative populations were identified using the patches generated above. These patches followed the same processing, except that HU12 watersheds were used to delineate populations. 

Because the > 0.8 threshold used above did not capture all colonized patches, all patches with predicted suitability scores of >0.55 which contained species observations were then identified using Terra's patches to create a data set of known populations. 

## Predicting Species Occurrence in Patches

All patches identified above were used to identify patches within 5 contiguous neighbors, or 5 kilometers, of a patch known to be occupied.
Occupied patches were determined using both the training and test occurrences data set used to generate the SDMs. 
To identify each patch within 5 orders of contiguous neighbors `nblag`, and to determine all patches within 5k of a populated patch `dnearneigh`, were used (@bivand2018spdep). 
For each of these non-occupied patches the number of occupied patches at different lags were counted.  

Each non-occupied patch was assigned an arbitrary rank based upon whether they were contiguous with an occupied patch, and if contiguous than their lag number to the nearest occupied patch, and the number of occupied patches connected (TABLE XX). 
The arbitrarily assigned numbers increase from '1', for an occupied patch, to '7' for a patch which has fewer than 3 second-order contiguous neighbors to an occupied patch(es). 

|  Connection   |  No. Connections   |   Rank   | 
| :-----------: | :----------------: | :------: | 
|     1^st^     |        >= 2        |    2     | 
|     1^st^     |         1          |    3     | 
|     2^nd^     |        >= 3        |    4     | 
|     2^nd^     |        <= 2        |    5     | 
|     3^rd^     |        >= 4        |    6     | 
|     3^rd^     |        <= 3        |    7     | 

Patches which have no contiguous neighbors, but which have neighbors within 5km were also assigned rank values based in this system (TABLE XX). 

|  No. Occupied Patches   |   Rank   | 
| :---------------------: | :------: | 
|           >= 3          |    5     |  
|           <= 2          |    6     | 

## Patch Metrics 

Multiple landscape metrics were calculated; the Class metrics Euclidean Nearest Neighbor of both mean and coefficient of variation, and the Patch metrics: Euclidean Nearest Neighbor, Core Area Index, PARA, Fractal Index (@hesselbarth2019lsm).


# Update

Unfortunately someone leaked some data (not me or anyone associated with our crews; actually a government employee) and we no longer have access to all of the data for evaluation. 
If you want to see my final report on this project please look here https://sagesteppe.github.io/Modelling-Natural-History-to-Increase-Native-Seed-Collection-Efficiency/ModellingNaturalHistory.html#/title-slide . 


## Appendix 2

Correlation between patch habitat suitable and abundance

Assess, Inventory, and Monitor data can be used to test the relationship between modeled habitat suitability and field measured abundance. Rather than simple relating the value of modeled suitability at a single cell, we can sum the modeled suitability of multiple cells, which together theoretically form a population. Patches may be identified using terras 'get_patches' function. 

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


## Appendix 1 

Variables and Sources for Environmental Niche Models

| Layer |                       Description                       |              Source                   |     Abbrev   |      
| :---: | :-----------------------------------------------------: | :-----------------------------------: | :----------: | 
|  1.   |              Mean Annual Air Temperature (BIO1)         |              Chelsa                   |     bio1     |
|  2.   |               Temperature seasonality (BIO4)            |              Chelsa                   |     bio4     |
|  3.   |         Max Temperature of Warmest Month (BIO5)         |              Chelsa                   |     bio5     |
|  4.   |         Min Temperature of Coldest Month (BIO6)         |              Chelsa                   |     bio6     |
|  5.   |        Mean Temperature of Warmest Quarter (BIO10)      |              Chelsa                   |     bio10    |
|  6.   |        Mean Temperature of Coldest Quarter (BIO11)      |              Chelsa                   |     bio11    |
|  7.   |              Mean annual precipitation (BIO12)          |              Chelsa                   |     bio12    |
|  8.   |         Precipitation of Warmest Quarter (BIO18)        |              Chelsa                   |     bio18    |
|  9.   |        Precipitation of Coldest Quarter (BIO19)         |              Chelsa                   |     bio19    |
| 10.   |      Mean Monthly vapour pressure deficit (vpd)         |              Chelsa                   |     vpd      |
| 11.   |     Heat accumulation of  Degree-days above 5C (gdd5)   |              Chelsa                   |     gdd5     |
| 12.   |         First growing degree day above 5C (gdgfgd5)     |              Chelsa                   |    gdgfgd5   |
| 13.   |           Number of Degree-days above 5C (ngd5)         |              Chelsa                   |     ngd5     |
| 14.   |   Heat accumulation of  Degree-days above 10C (gdd10)   |              Chelsa                   |     gdd10    |
| 15.   |         First growing degree day above 10C (gdgfgd10)   |              Chelsa                   |   gdgfgd10   |
| 16.   |         Number of Degree-days above 10C (ngd10)         |              Chelsa                   |     ngd10    |
| 17.   |        Mean monthly near surface humidity (hurs)        |              Chelsa                   |     hurs     |
| 18.   |            Number of Days with Snow Cover(scd)          |              Chelsa                   |     scd      |
| 19.   |            Annual Snow Water Equivalent (swe)           |              Chelsa                   |     swe      |
| 20.   |                Percent Herbaceous Vegetation            |              EarthEnv                 |     herb     |
| 21.   |                   Percent Shrub Cover                   |              EarthEnv                 |    shrub     |
| 22.   |                   Percent Tree Cover                    |              EarthEnv                 |     tree     |
| 23.   |               Soil Depth to Bedrock (R Horizon)         |              SoilGrids                | depth2bedrock|
| 24.   |                Soil organic carbon (Tonnes / ha)        |              SoilGrids                |    soc_0_5   |
| 25.   |                Soil Surface (0-5 cm) pH in H~2~O        |              SoilGrids                |   phh2o_0_5  |
| 26.   |                   Soil 30-60 cm pH in H~2~O             |              SoilGrids                |  phh2o_30_60 |
| 27.   |                Soil Surface (0-5 cm) % clay             |              SoilGrids                |   clay_0_5   |
| 28.   |                     Soil 5-15 cm % clay                 |              SoilGrids                |   clay_5_15  |
| 29.   |                Soil Surface (0-5 cm) % silt             |              SoilGrids                |   silt_0_5   |
| 30.   |                    Soil 5-15 cm  % silt                 |              SoilGrids                |   silt_5_15  |
| 31.   |                    Soil 15-30 cm % silt                 |              SoilGrids                |  silt_15_30  |
| 32.   |                Soil Surface (0-5 cm) % sand             |              SoilGrids                |   sand_0_5   |
| 33.   |                    Soil 5-15 cm % sand                  |              SoilGrids                |   sand_5_15  |
| 34.   |                    Soil 15-30 cm % sand                 |              SoilGrids                |  sand_15_30  |
| 35.   |           Soil Surface (0-5 cm) coarse fragments        |              SoilGrids                |    cfvo_0_5  | 
| 36.   |              Soil (30-60 cm) coarse fragments           |              SoilGrids                |   cfvo_30_60 | 
| 37.   |                       Soil Salinity                     |              SoilGrids                |   salinity   |
| 38.   |                         Elevation                       |             MERIT - DEM               |      dem     |
| 39.   |                          Slope                          |             Geomorpho90               |     slope    |
| 40.   |                          Aspect                         |             Geomorpho90               |    aspect    |
| 41.   |                  Topographic Position Index             |             Geomorpho90               |      tpi     |
| 42.   |                  Compound Topographic Index             |               terra                   |      cti     |
| 43.   |                  Topographic Roughness Index            |               terra                   |      tri     |
| 44.   |                   Human Influence Index                 |           NASA Earth Data             |   lastwild   |


data sources (@karger2017climatologies, @hengl2017soilgrids250m, @ivushkin2019global,  @tuanmu2014global, @yamazaki2017merit, @amatulli2020geomorpho90m, @sanderson2002human)

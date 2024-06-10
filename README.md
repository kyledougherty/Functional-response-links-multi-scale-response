# Functional-response-links-multi-scale-response
'Functional response links multi-scale response.R' contains code to reproduce results from the manuscript 'A Functional Response in Resource Selection Links Multi-Scale Responses of a Large Carnivore to Human Mortality Risk'. 

# Sections of 'Functional response links multi-scale response.R'

## Depredation Mortality Risk
This section contains code to:
* Download the following files:
  * Depredation_Mortality_Risk_Input_Data.csv
      * Each row represents a used or available location for use in fitting depredation mortality risk models. Used locations are locations where mountain lions were legally killed following livestock depredation, and available locations are locations systematically distributed throughout known or approximated home ranges.
      * Columns include:
        * ID: ID code for the mountain lion (DM#).
        * County: County code for the county where the mortality occured.
        * Used: 1 = used, 0 = available
        * GPS_Tracked: Logical column indicating whether the mountain lion was tracked with GPS telemetry.
        * Dist_Cover: Euclidean distance to the nearest occurence of forest or shrub cover (derived from NLCD Landcover).
        * Dist_Developed_Low: Euclidean distance to the nearest occurence of low-intensity development (derived from NLCD Landcover).
        * Dist_Developed_Open: Euclidean distance to the nearest developed open area (derived from NLCD Landcover).
        * Dist_Local_Roads: Euclidean distance to the nearest local road (derived from USGS National Transportation Dataset).
        * Distance Building: Euclidean distance to the nearest building (derived from Microsoft Buildings Dataset).
        * elevation: elevation as measured by Farr et al. 2007.
        * slope: local gradient computed using the 4-connected neighbors (30m pixels) around each location (derived from elevation layer).
        * HR_Sample_#: There are 100 columsn, each eith names ranging from HR_Sample_1 through HR_Sample_100. Each is a logical column indicating whether a location intersects the approximated home range in the given sample #.  
  * Depredation_Mortality_Risk_Layer.tif
      * Raster with 30x30m pixels that each contain a value (0-1) for predicted depredation mortality risk (derived from depredation mortality risk models).
  * Study_Area_Boundaries.gpkg
      * Geometry column for study area boundaries. 
  * State_Boundaries.gpkg
      * Geometry column for U.S. state bounadries.
  * California_Boundaries.gpkg
      * Geometry column for California county boundaries.
* Reproduce the following:
    * Depredation mortality risk models used to predict spatially varying depredation mortality risk.
    * Depredation mortality risk model using only data from individuals that were tracked with GPS telemetry (reported in Table S3)
    * Figure 2

## Coarse-Scale Resource Selection
This section contains code to: 
* Download the file:
  * Second_Order_RSF_Input_Data.csv
    * Each row represent a used or available location for use in fitting second-order resource selection functions. Used locations are locations systematically distributed throughout mountain lions home ranges, and available locations are locations systematically distributed throughout the study area.
    * Columns include:
      * ID: ID code for the mountain lion
      * Male: 1 = male, 0 = female
      * Study: Study area code for the mountain lion
      * Used: 1 = used, 0 = available
      * Depredation_Risk: Predicted depredation mortality risk
      * Dist_Herbaceous: Euclidean distance to the nearest occurence of herbaceous cover (derived from NLCD Landcover).
      * Dist_Shrub: Euclidean distance to the nearest occurence of shrub cover (derived from NLCD Landcover).
      * Dist_Forest: Euclidean distance to the nearest occurence of forest cover (derived from NLCD Landcover).
      * elevation: elevation as measured by Farr et al. 2007.
      * slope: local gradient computed using the 4-connected neighbors (30m pixels) around each location (derived from elevation layer).
* Reproduce the following:
    * Second-order resource selection function reported in the main text (shown in Figure 3A)
    * AICc table for second-order resource selection functions (reported in Table S4)
    * Second-order resource selection function with random slope to derive individual-level responses to depredation mortality risk.
    * Sex-specific second-order resource selection function (reported in Table S5). 

## Fine-Scale Resource Selection
This section contains code to: 
* Download the file:
  * SSF_Input_Data.csv
    * Each row represent a used or available location for use in fitting step-selection functions. Used locations are locations from GPS collars, and available locations are locations of potential steps generated from individual- and state-specific distributions of step lengths and turning angles.  
    * Columns include:
      * ID: ID code for the mountain lion
      * STUDY: Study area code for the mountain lion
      * Male: 1 = male, 0 = female
      * TOD: Time of day (day, night, crepuscular)
      * Night: 1 = Night, 0 = day/crepuscular
      * N_Available: Number of available steps per used location (1-10, 10, 15)
      * Interval_Length: Time between GPS fixes (1, 2, 3, or 4 hours)
      * State: Movement state derived from hidden Markov models ("Encamped", "Meandering", "Directed")
      * Directed: 1 = Directed movement state, 0 = Encamped/Meandering movement state
      * Step_Length: Distance in meters from previous used location
      * Turning_Angle: Turning angle (degrees) from previous bearing
      * Used: 1 = used, 0 = available
      * Depredation_Risk: Predicted depredation mortality risk
      * Dist_Herbaceous: Euclidean distance to the nearest occurence of herbaceous cover (derived from NLCD Landcover).
      * Dist_Shrub: Euclidean distance to the nearest occurence of shrub cover (derived from NLCD Landcover).
      * Dist_Forest: Euclidean distance to the nearest occurence of forest cover (derived from NLCD Landcover).
      * Elevation: elevation as measured by Farr et al. 2007.
      * Slope: local gradient computed using the 4-connected neighbors (30m pixels) around each location (derived from elevation layer).
* Reproduce the following:
  *  Step-selection function reported in the main text and shown in Figure 3B
  *  Individual-level functional response model with local availability (Reported in Table S11 and shown in Figure 5)
  *  Step-selection functions with interactions between resource variables and Male, Directed, and Night (reported in Tables S7-9)
  *  QIC table for second-order resource selection functions (reported in Table S6)
    
## Functional Responses in Resource Selection
This section contains code to reproduce the following: 
* Individual step-selection functions from data with 2 hour interval length
* Fit functional response models for depredation mortality risk, forest, shrub, and herbaceous areas (Reported in main text)
* Evaluate support for sex-specific depredation mortality risk functional response (Reported in main text).
* Evaluate if observed functional response to depredation mortality risk robust to uncertainty around individual-level beta coefficients (reported in main text and shown in Figure 4)
* Evaluate correlation between coarse-scale response to depredation mortality risk and mean depredation mortality risk within home range (reported in Table S10)

  ## Supporting Figures
  This section contains code to reproduce the following:
  * Figure S1: Figure showing beta coefficients and 95% confidence intervals from step selection functions evaluating movement-based resource selection of radio collared mountain lions (n = 144, 2 hour intervals) in California, USA (2004-2020) with varying numbers of available steps per used step. For distance-based variables, values less than 0 indicate selection, whereas values greater than 0 indicate avoidance. For non-distance-based variables, values less than 0 indicate avoidance, whereas values greater than 0 indicate selection. 
  * Figure S2: Figure beta coefficients and 95% confidence intervals from models evaluating landscape features influencing spatial variation in risk of mountain lions being killed on depredation permits (n = 466) in California, USA (2016-2020). We fit 100 sets of models, each with different approximations of home ranges for animals whose true home ranges were unknown. For distance-based variables, values less than 0 indicate increased mortality risk, whereas values greater than 0 indicate reduced mortality risk. For non-distance-based variables, values less than 0 indicate reduced mortality risk, whereas values greater than 0 indicate increased mortality risk.
  * Figure S3: Figure showing beta coefficients and 95% confidence intervals from step selection functions evaluating movement-based resource selection of radio collared mountain lions (n = 244) in California, USA (2004-2020). For distance-based variables, values less than 0 indicate selection, whereas values greater than 0 indicate avoidance. For non-distance-based variables, values less than 0 indicate avoidance, whereas values greater than 0 indicate selection.

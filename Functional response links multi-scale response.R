library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggsn)
library(ggrepel)
library(patchwork)
library(glmmTMB)
library(survival)
library(broom)
library(parallel)
library(MuMIn)
library(bfp)

# Data to reproduce results is provided via onedrive links and code throughout
# this file will download and store in relevant directories. 
# Note that some files are rather large, the download timeout below may need to
# be adjusted based on download speeds. 
options(timeout = max(7200, getOption("timeout")))
dir.create("Data")
# Depredation mortality risk ---------------------------------------------------

## Download depredation mortality risk data using onedrive link and store locally
download.file("https://uofnelincoln-my.sharepoint.com/:x:/g/personal/kdougherty8_unl_edu/EXumiNHNNDFGr7JODDgbVeEBmx6DchPLfiTLFOQK5GDLQA?e=mc8Xlw&download=1", 
              "Data/Depredation_Mortality_Risk_Data/Depredation_Mortality_Risk_Input_Data.csv")
## Create directory to store results
dir.create("Results/Depredation_Mortality_Risk_Models", recursive = TRUE)

## Fit depredation mortality risk models reported in the main text. 
## Each of the 100 model sets uses different approximations of home ranges for 
## individuals that were not tracked with GPS telemetry. 
map(paste0("HR_Sample_", 1:100),
    function(hr_sample){
      
      # browser()
      # Import depredation mortality risk data and retain only columns containing 
      # covariates, and the relevant HR_Sample_#  
      #   - The HR_Sample_# is a logical column indicating whether a 
      #     location intersected the approximated home range. 
      depredation_mortality_risk_data = read_csv("Data/Depredation_Mortality_Risk_Data/Depredation_Mortality_Risk_Input_Data.csv", 
                                                 col_select = c(1:15, hr_sample), 
                                                 lazy = TRUE) %>%
        # Remove available locations that do not intersect the current 
        # approximation of a home range for individuals not tracked with 
        # GPS telemetry
        filter(eval(parse(text = hr_sample)) == TRUE | 
                 GPS_Tracked == TRUE | 
                 Used == 1) %>%
        # Rescale each continuous variable by
        # subtracting its mean and dividing by two
        # standard deviations
        mutate(across(c(Dist_Cover:slope),
                      function(x){
                        (x - mean(x))/(2*sd(x))
                      }))
      
      # Create a correlation matrix and use it to create a subset expression that
      # will prevent covariates with Pearson's correlation coefficients greater
      # than or equal to |0.6| from being included in the same model.
      correlation_matrix = cor(depredation_mortality_risk_data %>%
                                 select(Dist_Developed_Open,
                                        Dist_Developed_Low,
                                        Distance_Developed,
                                        Dist_Cover,
                                        Dist_Local_Roads,
                                        Distance_Building,
                                        slope))
      
      corr_matrix_include = abs(correlation_matrix) <= .6
      corr_matrix_include[!lower.tri(corr_matrix_include)] = NA
      i = as.vector(corr_matrix_include == FALSE & !is.na(corr_matrix_include))
      nm = colnames(corr_matrix_include)
      
      subset_expression = parse(text = paste("!(", paste("(",
                                                         nm[col(corr_matrix_include)[i]], "&&",
                                                         nm[row(corr_matrix_include)[i]], ")",
                                                         sep = "",
                                                         collapse = " || "),
                                             ")"))
      
      # Fit global model with strata(ID) and cluster(County)
      global = clogit(Used ~ Dist_Developed_Open + Dist_Developed_Low + Distance_Developed + 
                        Distance_Building + Dist_Local_Roads + Dist_Cover +
                        slope + strata(ID) + cluster(County),
                      data = depredation_mortality_risk_data,
                      method = "efron",
                      na.action = na.fail)
      
      # Generate formulas with all possible combinations of covariates
      model_formulas = dredge(global,
                              fixed = "strata(ID)",
                              m.lim = c(3, NA),
                              subset = subset_expression,
                              evaluate = FALSE)
      
      rm(list = setdiff(ls(), c("depredation_mortality_risk_data", "model_formulas", "correlation_matrix", "hr_sample", "start")),
         envir = environment())
      gc()
      
      # Fit models in parallel using 8 cores
      cl = makeCluster(8)
      clusterExport(cl, list("model_formulas", "depredation_mortality_risk_data"),
                    envir=environment())
      clusterEvalQ(cl, {
        library(tidyverse)
        library(broom)
        library(survival)
        library(MuMIn)})
      
      model_summaries = parLapply(cl,
                                  model_formulas,
                                  function(formula){
                                    
                                    # Fit Model
                                    model = eval(formula)
                                    
                                    # Sumarise Model
                                    model_summary = tibble(Formula = paste(names(model$coefficients), collapse = "+"),
                                                           QIC = MuMIn::QIC(model), 
                                                           Summary = list(broom::tidy(model) %>%
                                                                            select(-statistic) %>%
                                                                            left_join(as_tibble(summary(model)$conf.int, 
                                                                                                rownames = "term") %>%
                                                                                        mutate(across(contains("95"), 
                                                                                                      log)) %>%
                                                                                        select(term, contains("95")), 
                                                                                      by = "term")))
                                    
                                    model_summary
                                    
                                  }) %>%
        bind_rows() %>%
        # Bind into a single QIC table
        arrange(QIC)
      
      stopCluster(cl)
      
      # Save model summaries 
      saveRDS(model_summaries,
              file = paste0("Results/Depredation_Mortality_Risk_Models/Depredation_Mortality_Risk_Models_",
                            hr_sample,
                            ".RDS"))
      gc()
    })

## Summarize most strongly supported depredation mortality risk models
Depredation_Mortality_Risk_Models <- map(paste0("HR_Sample_", seq(1, 100)), 
                                            function(hr_sample){
                                              
                                              model_summaries = readRDS(paste0("Results/Depredation_Mortality_Risk_Models/Depredation_Mortality_Risk_Models_",
                                                                               hr_sample,
                                                                               ".RDS"))
                                              
                                              model_summaries %>% 
                                                mutate(HR_Sample = hr_sample, 
                                                       Delta_QIC = QIC - min(QIC)) %>% 
                                                select(HR_Sample, Formula:QIC, Delta_QIC, Summary)
                                              
                                            })

## Summarize most strongly supported depredation mortality risk models 
## with mean/sd of beta coefficients for the 100 sets (Reported in results of main text) 
Depredation_Mortality_Risk_Model_Summary  <- Depredation_Mortality_Risk_Models %>% 
  bind_rows() %>%
  group_by(HR_Sample) %>% 
  filter(QIC == min(QIC)) %>% 
  unnest(Summary) %>% 
  group_by(term) %>% 
  summarize(Mean = round(mean(estimate), 2),
            SD = round(sd(estimate), 2))

## Fit depredation mortality risk model with only individuals that were tracked 
## with GPS telemetry. 
GPS_Depredation_Mortality_Risk_Data = read_csv("Data/Depredation_Mortality_Risk_Data/Depredation_Mortality_Risk_Input_Data.csv", 
                                               col_select = c(1:15), 
                                               lazy = TRUE) %>% 
  filter(GPS_Tracked == TRUE) %>%
  # Rescale each continuous variable by
  # subtracting its mean and dividing by two
  # standard deviations
  mutate(across(c(Dist_Cover:slope),
                function(x){
                  (x - mean(x))/(2*sd(x))
                }))

## Create a correlation matrix and use it to create a subset expression that
## will prevent covariates with Pearson's correlation coefficients greater
## than or equal to |0.6| from being included in the same model.
correlation_matrix = cor(GPS_Depredation_Mortality_Risk_Data %>%
                           select(Dist_Developed_Open,
                                  Dist_Developed_Low,
                                  Distance_Developed,
                                  Dist_Cover,
                                  Dist_Local_Roads,
                                  Distance_Building,
                                  slope))

corr_matrix_include = abs(correlation_matrix) <= .6
corr_matrix_include[!lower.tri(corr_matrix_include)] = NA
i = as.vector(corr_matrix_include == FALSE & !is.na(corr_matrix_include))
nm = colnames(corr_matrix_include)
subset_expression = parse(text = paste("!(", paste("(",
                                                   nm[col(corr_matrix_include)[i]], "&&",
                                                   nm[row(corr_matrix_include)[i]], ")",
                                                   sep = "",
                                                   collapse = " || "),
                                       ")"))

## Fit global model with strata(ID) and cluster(County)
Global_GPS_Depredation_Risk_Model = clogit(Used ~ Dist_Developed_Open + Dist_Developed_Low + Distance_Developed +
                  Distance_Building + Dist_Local_Roads + Dist_Cover +
                  slope + strata(ID) + cluster(County),
                data = GPS_Depredation_Mortality_Risk_Data,
                method = "efron",
                na.action = na.fail)

## Generate formulas with all possible combinations of covariates
GPS_Depredation_Risk_Formulas = dredge(Global_GPS_Depredation_Risk_Model,
                                       fixed = "strata(ID)",
                                       m.lim = c(3, NA),
                                       subset = subset_expression,
                                       evaluate = FALSE)

## Fit depredation mortality risk models using only data from individuals 
## that were GPS tracked
GPS_Depredation_Risk_Models <- map_df(GPS_Depredation_Risk_Formulas, 
                                      function(formula){
                                        
                                        model = eval(formula)
                                        
                                        model_summary = tibble(Formula = paste(names(model$coefficients), collapse = "+"),
                                                               QIC = MuMIn::QIC(model), 
                                                               Summary = list(broom::tidy(model) %>%
                                                                                select(-statistic) %>%
                                                                                left_join(as_tibble(summary(model)$conf.int, 
                                                                                                    rownames = "term") %>%
                                                                                            mutate(across(contains("95"), 
                                                                                                          log)) %>%
                                                                                            select(term, contains("95")), 
                                                                                          by = "term")))
                                        
                                        model_summary
                                        
                                      })

Table_S3 <- GPS_Depredation_Risk_Models  %>% 
  filter(QIC == min(QIC)) %>% 
  unnest(Summary) %>% 
  mutate(estimate = round(estimate, 2), 
         `lower .95` = round(`lower .95`, 2), 
         `upper .95` = round(`upper .95`, 2), 
         p.value = round(p.value, 3)) %>% 
  select(term, estimate, `lower .95`, `upper .95`, p.value)

## Code to produce Figure 2
##  Download data needed to reproduce figure 2: 
##   -
##   -
##   -
##   -
getwd()
dir.create("Data/Figure_2_Data")
download.file("https://uofnelincoln-my.sharepoint.com/:u:/g/personal/kdougherty8_unl_edu/ETvsVLhIw7JMm-9dMrE0uTUB2ncUDe34EInQdMVefPqHEQ?e=xursKc&download=1",
              mode = "wb",
              "Data/Figure_2_Data/State_Boundaries.gpkg")
download.file("https://uofnelincoln-my.sharepoint.com/:u:/g/personal/kdougherty8_unl_edu/EQCaVwSOtkxKgFI1oqqvtmABL9ALt3P5TKc23-mEdoerrQ?e=ko30yy&download=1", 
              mode = "wb",
              "Data/Figure_2_Data/California_Boundaries.gpkg")
download.file("https://uofnelincoln-my.sharepoint.com/:u:/g/personal/kdougherty8_unl_edu/EWBlRhgCmfFIsLt6rACsE5wB1xCPWsnUbscyP6lfgSNLJg?e=SGJJh4&download=1", 
              mode = "wb",
              "Data/Figure_2_Data/Study_Area_Boundaries.gpkg")
download.file("https://uofnelincoln-my.sharepoint.com/:i:/g/personal/kdougherty8_unl_edu/Ec1YM4vsxXpKguvIImxZH8QBHtgSY-l_5PIGDMJV6Dnasg?e=lQXmtw&download=1", 
              mode = "wb",
              "Data/Figure_2_Data/Depredation_Mortality_Risk_Layer.tif")

## Import depredation mortality risk layer 
Depredation_Risk <- rast("Data/Figure_2_Data/Depredation_Mortality_Risk_Layer.tif")

## Import state boundaries shapefile
States <- st_read("Data/Figure_2_Data/State_Boundaries.gpkg") %>%
  filter(!NAME %in% c("Alaska", "Hawaii", "Guam", "United States Virgin Islands", 
                      "Puerto Rico",
                      "American Samoa", "Commonwealth of the Northern Mariana Islands")) %>%
  st_union()

## Import California County Boundaries, remove islands, and transform to CRS fo the depredation
## mortality risk layer
CA <- st_read("Data/Figure_2_Data/California_Boundaries.gpkg") %>% 
  st_cast("POLYGON") %>%
  mutate(Area = as.numeric(st_area(.))) %>% 
  filter(Area > 1e9) %>% 
  st_union() %>%
  st_transform(st_crs(Depredation_Risk)) 

Study_Areas <- st_read("Data/Figure_2_Data/Study_Area_Boundaries.gpkg")

## Crop depredation mortality risk layer to California boundaries
Depredation_Risk <- Depredation_Risk %>% 
  mask(CA %>% vect())

## Create dataframe with City Coordinates and Labels
Cities <- tibble(City = c("Los Angeles", "San Francisco", "San Jose", "San Diego"), 
                 X = c(-118.24, -122.43, -121.89, -117.16), 
                 Y = c(34.05, 37.67, 37.33, 32.71)) %>%
  st_as_sf(coords = c("X", "Y"), 
           crs = 4326, 
           remove = FALSE)

Main_Panel <- ggplot() + 
  geom_spatraster(data = Depredation_Risk, 
                  maxcell = 1e6) + 
  geom_sf(data = CA, 
          fill = NA, 
          color = "black", 
          lwd = 1.5) + 
  scale_fill_gradientn(colors = viridis::viridis(100),# colorRamps::blue2red(100), 
                       limits = c(0, 0.962),
                       breaks = c(0, 0.962),
                       labels = c("      Low", "High       "),
                       na.value = "white") + 
  geom_sf(data = Study_Areas %>% 
            mutate(Label = "Study Areas"),
          aes(color = Label), 
          fill = NA, 
          # color = "red", 
          lwd = 3) + 
  scale_color_manual(values = c("Study Areas" = "red")) +
  guides(color = guide_legend(title = "", 
                              # title.position = "top",
                              label.position = "top",
                              keywidth = unit(12, "cm"),
                              order = 1), 
         fill = guide_colorbar(ticks.colour = NA,
                               title = "Depredation Mortality Risk",
                               title.position = "top",
                               order = 2)) + 
  # ggsn::north(data = ext(Depredation_Risk) %>% 
  #                  as.polygons() %>%
  #                  st_as_sf() %>% 
  #                  st_set_crs(st_crs(Depredation_Risk)), 
  #             location = "bottomleft", 
  #             anchor = c(x = -119.5, y = 40.75)) +
  ggsn::scalebar(data = ext(Depredation_Risk) %>% 
                   as.polygons() %>%
                   st_as_sf() %>% 
                   st_set_crs(st_crs(Depredation_Risk)), 
                 transform = TRUE, 
                 dist = 200,
                 location = "bottomleft",
                 anchor = c(x = -119.25, y = 40.25),
                 dist_unit = "km", 
                 st.size = 10) + 
  # geom_rect(aes(xmin = -126,
  #               xmax = -121,
  #               ymin = 32,
  #               ymax = 34),
  #           fill = NA,
  #           color = "black") +
  # annotate("text",
  #          label = "Study Areas and Depredation Mortality Risk",
  #          fontface = "bold",
  #          size = 7,
  #          x = -116.75,
#          y = 41.75) +
geom_label_repel(data = Cities %>%
                   filter(City %in% c("San Francisco",
                                      "San Jose")),
                 aes(x = X,
                     y = Y,
                     label = City),
                 size = 10,
                 nudge_x = -3) +
  geom_label_repel(data = Cities %>%
                     filter(!City %in% c("San Francisco",
                                         "San Jose")),
                   aes(x = X,
                       y = Y,
                       label = City),
                   size = 10,
                   nudge_x = -2) +
  geom_label(data = tibble(Label = "Sierra Nevada\nMountains",
                           X = -120.25,
                           Y = 38),
             aes(x = X,
                 y = Y,
                 label = Label),
             size = 9) +
  geom_label(data = tibble(Label = "Modoc Plateau",
                           X = -120.822888,
                           Y = 40.75),
             aes(x = X,
                 y = Y,
                 label = Label),
             size = 10) +
  geom_label(data = tibble(Label = "North Coast",
                           X = -123.5,
                           Y = 40.4),
             aes(x = X,
                 y = Y,
                 label = Label),
             size = 10) +
  theme_void() + 
  theme(legend.position = c(0.5425, 0.87), 
        legend.direction = "horizontal",
        legend.key.width = unit(2.5, "cm"),
        legend.justification = "left",
        legend.title = element_text(face = "bold", 
                                    size = 20,
                                    vjust = 1, 
                                    hjust = 0.5),
        legend.text = element_text(size = 20, 
                                   face = "bold"))

CA_Inset <- ggplot() + 
  geom_sf(data = States, 
          fill = "grey50", 
          color = "grey50") + 
  geom_sf(data = CA, 
          fill = "black",
          color = "black") + 
  theme_void()

Main_Panel + 
  inset_element(CA_Inset, 
                left = -0.05, 
                right = 0.425, 
                bottom = -0.05,
                top = 0.425)

# Coarse-scale resource selection ----------------------------------------------

## Download second-order RSF data using onedrive link and store locally
dir.create("Data/Second_Order_RSF_Data", recursive = TRUE)
download.file("https://uofnelincoln-my.sharepoint.com/:x:/g/personal/kdougherty8_unl_edu/EeLR36Jh23BKshuDrsrkZtYBtnypRwjf3jL2d-jFMWvcfA?e=Dg69Ti&download=1", 
              "Data/Second_Order_RSF_Data/Second_Order_RSF_Input_Data.csv")

## Create directory to store results
dir.create("Results/Second_Order_RSF_Models", recursive = TRUE)

## Import data for second-order RSF
Second_Order_RSF_Data <- read_csv("Data/Second_Order_RSF_Data/Second_Order_RSF_Input_Data.csv")

## Fit Second Order RSF Reported in Main Text ----
Main_Second_Order_RSF <- glmmTMB(Used ~ Depredation_Risk + Dist_Forest + Dist_Shrub + Dist_Herbaceous + 
                                   elevation + slope + (1 | Study/ID), 
                                 family = binomial(link = "logit"), 
                                 data = Second_Order_RSF_Data %>%
                                   # Rescale each continuous variable by
                                   # subtracting its mean and dividing by two 
                                   # standard deviations
                                   mutate(across(c(Depredation_Risk:slope),
                                                 function(x){
                                                   (x - mean(x))/(2*sd(x))
                                                 })))

saveRDS(Main_Second_Order_RSF, 
        "Results/Second_Order_RSF_Models/Main_Second_Order_RSF.RDS")

Second_Order_Model_Summary <- Main_Second_Order_RSF %>% 
  mutate(Var = case_when(Var == "Depredation_Risk" ~ "Mortality \nRisk", 
                         Var == "Dist_Forest" ~ "Distance to \n Forest", 
                         Var == "Dist_Herbaceous" ~ "Distance to \nHerbaceous",
                         Var == "Dist_Shrub" ~ "Distance to \nShrub",
                         Var == "elevation" ~ "Elevation", 
                         Var == "slope" ~ "Slope"), 
         Var = factor(Var, 
                      levels = c("Mortality \nRisk", 
                                 "Distance to \n Forest",
                                 "Distance to \nShrub",
                                 "Distance to \nHerbaceous",
                                 "Elevation", 
                                 "Slope")))

## Create Figure 3a
ggplot(data = Second_Order_Model_Summary, 
                              aes(x = Var, 
                                  y = Estimate.x, 
                                  ymin = `2.5 %`, 
                                  ymax = `97.5 %`)) + 
  geom_hline(yintercept = 0, 
             col = "red", 
             lty = 2) + 
  geom_errorbar(position = position_dodge2(width = 0.25), 
                width = 0.05) + 
  geom_point(size = 1.5, 
             position = position_dodge2(width = 0.25)) + 
  labs(title = "A) Coarse-Scale Resource Selection", 
       y = "Beta Coefficient") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0,
                                   vjust = 0.5,
                                   color = "black", 
                                   size = 14, 
                                   face = "bold"),
        axis.text.y = element_text(color = "black", 
                                   size = 14, 
                                   face = "bold"),
        axis.title.y = element_text(color = "black", 
                                    size = 16, 
                                    face = "bold"),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "black", 
                                  size = 16, 
                                  face = "bold", 
                                  vjust = -6, 
                                  hjust = 0.02), 
        plot.tag = element_text(face = "bold", 
                                size = 16), 
        plot.tag.position = c(0.1, 0.98))

## Fit Second Order RSF Null Model ----
Second_Order_RSF_Null_Model <- glmmTMB(Used ~ (1|Study/ID), 
                                       family = binomial(link = "logit"), 
                                       data = Second_Order_RSF_Data)

saveRDS(Second_Order_RSF_Null_Model, 
        "Results/Second_Order_RSF_Models/Second_Order_RSF_Null_Model.RDS")

## Fit Sex-Specific Second Order RSF ----
Sex_Specific_Second_Order_RSF <- glmmTMB(Used ~ Depredation_Risk*Male + Dist_Forest*Male + 
                                           Dist_Shrub*Male + Dist_Herbaceous*Male + 
                                           elevation*Male + slope*Male + (1 | Study/ID), 
                                         family = binomial(link = "logit"), 
                                         data = RSF_Data %>%
                                           # Rescale each continuous variable by
                                           # subtracting its mean and dividing by two 
                                           # standard deviations
                                           mutate(across(c(Depredation_Risk:slope),
                                                         function(x){
                                                           (x - mean(x))/(2*sd(x))
                                                         })))

saveRDS(Sex_Specific_Second_Order_RSF, 
        "Results/Second_Order_RSF_Models/Sex_Specific_Second_Order_RSF.RDS")

## Fit Random Slope Second Order RSF ----
Random_Slope_Second_Order_RSF <- glmmTMB(Used ~ Depredation_Risk + Dist_Forest + Dist_Herbaceous + 
                                           Dist_Shrub + elevation + slope + (1 + Depredation_Risk|ID), 
                                         family = binomial(link = "logit"), 
                                         data = Second_Order_RSF_Data %>%
                                           # Rescale each continuous variable by
                                           # subtracting its mean and dividing by two 
                                           # standard deviations
                                           mutate(across(c(Depredation_Risk:slope),
                                                         function(x){
                                                           (x - mean(x))/(2*sd(x))
                                                         })))

saveRDS(Random_Slope_Second_Order_RSF, 
        "Results/Second_Order_RSF_Models/Random_Slope_Second_Order_RSF.RDS")

## Import Models Fit Above
Main_Second_Order_RSF <- readRDS("Results/Second_Order_RSF_Models/Main_Second_Order_RSF.RDS")
Sex_Specific_Second_Order_RSF <- readRDS("Results/Second_Order_RSF_Models/Sex_Specific_Second_Order_RSF.RDS")
Second_Order_RSF_Null_Model <- readRDS("Results/Second_Order_RSF_Models/Second_Order_RSF_Null_Model.RDS")

## AICc reported in Supporting Information
Table_S4 <- list(`Null Model` = Second_Order_RSF_Null_Model, 
                `No Interactions` = Main_Second_Order_RSF, 
                `Male x Resource Variables` = Sex_Specific_Second_Order_RSF) %>% 
  map_df(~tibble(AICc = MuMIn::AICc(.x)), 
         .id = "Model") %>% 
  arrange(AICc) %>% 
  mutate(Delta_AICc = AICc - min(AICc))

## Sex-specific second-order RSF model summary reported in Supporting Information
Table_S5 <- summary(Sex_Specific_Second_Order_RSF)$coeff$cond %>%
  as_tibble(rownames = "term") %>%
  left_join(confint(Sex_Specific_Second_Order_RSF, 
                    method = "wald", 
                    level = 0.95, 
                    estimate = FALSE,
                    component = "cond")%>%
              as_tibble(rownames = "term"), 
            by = "term")

# Fine-scale resource selection ------------------------------------------------

## Download SSF data using onedrive link and store locally
dir.create("Data/SSF_Data", recursive = TRUE)
download.file("https://uofnelincoln-my.sharepoint.com/:x:/g/personal/kdougherty8_unl_edu/Ec6IkgLuVJpGq9nvkXWsP4UBVSjOXUx9bImJlbVRombMyQ?e=SiOkgm&download=1", 
              "Data/SSF_Data/SSF_Input_Data.csv")

## Import data for SSFs and split into separate datasets for: 
##   - Number of available locations: (1-10, 10, 15)
##   - Interval Length: 1-4 hours
SSF_Data <- read_csv("Data/SSF_Data/SSF_Input_Data.csv") %>% 
  group_by(Interval_Length) %>% 
  group_split() %>% 
  set_names(map(., ~unique(.x$Interval_Length))) %>%
  map(function(ssf_data){
    ssf_data %>% 
      group_by(N_Available) %>% 
      group_split() %>% 
      set_names(map(., ~unique(.x$N_Available)))
  })

## Fit SSF Reported in Main Text ----
##    - Interval Length = 2 Hours, 
##    - n available steps = 5
Main_SSF <- clogit(Used ~ Depredation_Risk + Dist_Herbaceous + 
                     Dist_Shrub + Dist_Forest + Elevation + Slope + 
                     strata(StepID) + cluster(ID), 
                   # Retain only steps in the 'Meandering' or 'Directed' 
                   # movement state and rescale each continuous variable by
                   # subtracting its mean and dividing by two standard deviations
                   data = SSF_Data$`2`$`5` %>%
                     filter(State %in% c("Meandering", "Directed")) %>% 
                     mutate(across(c(Depredation_Risk:Slope), 
                                   function(x){
                                     (x - mean(x))/(2*sd(x))
                                   })), 
                   method = "approximate")

## Summarise SSF Reported in Main Text: 
##   - Beta Coefficients
##   - 95% Confidence Intervals
Main_SSF_Summary <- tidy(Main_SSF) %>%
  left_join(confint(Main_SSF) %>% 
              as_tibble(rownames = "term")) %>% 
  mutate(term = factor(term, 
                       levels = c("Depredation_Risk", 
                                  "Dist_Forest", 
                                  "Dist_Shrub", 
                                  "Dist_Herbaceous", 
                                  "Elevation", 
                                  "Slope"), 
                       labels = c("Mortality \nRisk", 
                                  "Distance to \n Forest", 
                                  "Distance to \n Shrub", 
                                  "Distance to \n Herbaceous", 
                                  "Elevation", 
                                  "Slope")),
         across(c(estimate, `2.5 %`, `97.5 %`), ~round(.x, 2))) %>%
  select(term, estimate, `2.5 %`, `97.5 %`, p.value)

## Create Figure 3b
ggplot(Main_SSF_Summary, 
       aes(x = term, 
           y = estimate, 
           ymin = `2.5 %`, 
           ymax = `97.5 %`)) + 
  geom_hline(yintercept = 0, 
             col = "red", 
             lty = 2) + 
  geom_errorbar(position = position_dodge2(width = 0.25), 
                width = 0.25) + 
  geom_point(size = 2, 
             position = position_dodge2(width = 0.25)) + 
  labs(y = "Beta Coefficients") + 
  theme_classic() + 
  theme(axis.text = element_text(angle = 0, 
                                 vjust = 0.5, 
                                 size = 11,
                                 face = "bold",
                                 color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14,
                                    face = "bold", 
                                    color = "black"))

## Fit individual-level functional response model with interactions between 
## depredation mortality risk, herbaceous cover, shrub cover, and mean local 
## depredation mortality risk
Individual_Level_FR_SSF_Data <- SSF_Data$`2`$`5` %>%
  filter(State %in% c("Meandering", "Directed")) %>% 
  # Add mean local depredation mortality risk at available locations
  left_join(SSF_Data$`2`$`5` %>% 
              group_by(ID, StepID, Used) %>%
              summarise(Mean_Local_Depredation_Risk = mean(Depredation_Risk)) %>%
              filter(Used == 0) %>%
              select(-Used)) %>%
  ungroup()

Individual_Level_FR_SFF <- clogit(Used ~ Depredation_Risk + Depredation_Risk:Mean_Local_Depredation_Risk + 
                                    Dist_Herbaceous + Dist_Herbaceous:Mean_Local_Depredation_Risk + 
                                    Dist_Shrub + Dist_Shrub:Mean_Local_Depredation_Risk + 
                                    Dist_Forest + Dist_Forest:Mean_Local_Depredation_Risk + 
                                    Elevation + 
                                    Slope + 
                                    strata(StepID) + cluster(ID), 
                                  data = Individual_Level_FR_SSF_Data %>%
                                    # Retain only steps in the 'Meandering' or 'Directed' 
                                    # movement state and rescale each continuous variable by
                                    # subtracting its mean and dividing by two standard deviations
                                    mutate(across(c(Depredation_Risk:Slope, Mean_Local_Depredation_Risk), 
                                                  function(x){
                                                    (x - mean(x))/(2*sd(x))
                                                  })), 
                                  method = "approximate", 
                                  na.action = na.fail)

# Summarize individual-level functional response model
Table_S11 <- tidy(Individual_Level_FR_SFF) %>%
  left_join(confint(Individual_Level_FR_SFF) %>% 
              as_tibble(rownames = "term")) %>% 
  mutate(term = factor(term, 
                       levels = c("Depredation_Risk", 
                                  "Depredation_Risk:Mean_Local_Depredation_Risk", 
                                  "Dist_Forest", 
                                  "Dist_Forest:Mean_Local_Depredation_Risk", 
                                  "Dist_Herbaceous", 
                                  "Dist_Herbaceous:Mean_Local_Depredation_Risk", 
                                  "Dist_Shrub", 
                                  "Dist_Shrub:Mean_Local_Depredation_Risk", 
                                  "Elevation", 
                                  "Slope"),
                       labels = c("Mortality Risk", 
                                  "Mortality Risk x Mean Local Mortality Risk", 
                                  "Distance to Forest", 
                                  "Distance to Forest x Mean Local Mortality Risk", 
                                  "Distance to Herbaceous", 
                                  "Distance to Herbaceous x Mean Local Mortality Risk", 
                                  "Distance to Shrub", 
                                  "Distance to Shrub x Mean Local Mortality Risk", 
                                  "Elevation",  
                                  "Slope")), 
         across(c(estimate, `2.5 %`, `97.5 %`), ~round(.x, 2))) %>%
  select(term, estimate, `2.5 %`, `97.5 %`, p.value) %>% 
  arrange(term)

# Create data for predicting conditional probability of use shown in Figure 5
#   - Low (0.25), Moderate (0.5), and High (0.75) Mortality Risk
#   - Mean Local risk ranging from 0-1
expand.grid(Depredation_Risk = c(0.25, 0.5, 0.75), 
            Mean_Local_Depredation_Risk = seq(0, 
                                              1, 
                                              length.out = 100)) %>% 
  # Set all other values to 0 (mean value after rescaling)
  bind_cols(tibble(Elevation = 0, 
                   Slope = 0, 
                   Dist_Forest = 0, 
                   Dist_Shrub = 0, 
                   Dist_Herbaceous = 0, 
                   StepID = unique(SSF_Data$`2`$`5`$StepID)[1])) %>% 
  # Add columns for depredation mortality risk and mean local depredation mortality
  # risk on original scale (OS) and rescaled using mean/sd from the SSF data
  mutate(Depredation_Risk_OS = factor(Depredation_Risk, labels = c("Low (0.25)", "Moderate (0.50)", "High (0.75)")), 
         Mean_Local_Depredation_Risk_OS = Mean_Local_Depredation_Risk, 
         Mean_Local_Depredation_Risk = (Mean_Local_Depredation_Risk - mean(Individual_Level_FR_SSF_Data$Mean_Local_Depredation_Risk))/
           (2 * sd(Individual_Level_FR_SSF_Data$Mean_Local_Depredation_Risk)),
         Depredation_Risk = (Depredation_Risk - mean(Individual_Level_FR_SSF_Data$Depredation_Risk))/
           (2 * sd(Individual_Level_FR_SSF_Data$Depredation_Risk))) %>% 
  # Predict conditional probability of use and bind columns
  bind_cols(predict(Individual_Level_FR_SFF, 
                    newdata = ., 
                    se.fit = TRUE,
                    type = "lp", 
                    reference = "sample")) %>% 
  mutate(Conditional_Probability = 1/(1 + exp(-fit)), 
         LCL = 1/(1 + exp(-(fit - 1.959964 * se.fit))), 
         UCL = 1/(1 + exp(-(fit + 1.959964 * se.fit)))) %>% 
  ggplot(aes(y = Conditional_Probability, 
             ymin = LCL, 
             ymax = UCL,
             x = Mean_Local_Depredation_Risk_OS, 
             group = Depredation_Risk)) + 
  lims(y = c(0, 1)) + 
  labs(x = "Mean Mortality Risk at Available Steps", 
       y = "Conditional Probability of Use") +
  geom_ribbon(aes(fill = Depredation_Risk_OS), 
              alpha = 0.5) + 
  geom_line(aes(col = Depredation_Risk_OS), 
            size = 1.5) + 
  scale_fill_manual(values = c("Low (0.25)" = viridis::viridis(3)[[1]], 
                               "Moderate (0.50)" = viridis::viridis(3)[[2]], 
                               "High (0.75)" = viridis::viridis(3)[[3]])) + 
  scale_color_manual(values = c("Low (0.25)" = viridis::viridis(3)[[1]], 
                                "Moderate (0.50)" = viridis::viridis(3)[[2]], 
                                "High (0.75)" = viridis::viridis(3)[[3]])) + 
  labs(col = "Mortality Risk", 
       fill = "Mortality Risk") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0,
                                   vjust = 0.5,
                                   color = "black", 
                                   size = 14, 
                                   face = "bold"),
        axis.text.y = element_text(color = "black", 
                                   size = 14, 
                                   face = "bold"),
        axis.title.x = element_text(color = "black", 
                                    size = 16, 
                                    face = "bold"),
        axis.title.y = element_text(color = "black", 
                                    size = 16, 
                                    face = "bold"),
        plot.tag = element_text(face = "bold", 
                                size = 16), 
        legend.position = c(0.85, 0.85),
        legend.title = element_text(face = "bold", 
                                    color = "black", 
                                    size = 14), 
        legend.text = element_text(face = "bold", 
                                   color = "black", 
                                   size = 12))

## Fit SSFs with Interactions Reported in Supporting Information ----
##    - Resource Variables x Male (1 = Male, 0 = Female)
##    - Resource Variables x Directed (1 = Directed, 0 = Meandering)
##    - Resource Variables x Night (1 = Night, 0 = all other times)
Sex_Specific_2_Hour_SSF <- clogit(Used ~ Depredation_Risk + Depredation_Risk:Male + 
                                    Dist_Herbaceous + Dist_Herbaceous:Male + 
                                    Dist_Shrub + Dist_Shrub:Male + 
                                    Dist_Forest + Dist_Forest:Male + 
                                    Elevation + Elevation:Male + 
                                    Slope + Slope:Male + 
                                    strata(StepID) + cluster(ID), 
                                  # Retain only steps in the 'Meandering' or 'Directed' 
                                  # movement state and rescale each continuous variable by
                                  # subtracting its mean and dividing by two standard deviations
                                  data = SSF_Data$`2`$`5` %>%
                                    filter(State %in% c("Meandering", "Directed")) %>% 
                                    ungroup() %>%
                                    mutate(across(c(Depredation_Risk:Slope), 
                                                  function(x){
                                                    (x - mean(x))/(2*sd(x))
                                                  })), 
                                  method = "approximate")

State_Specific_2_Hour_SSF <- clogit(Used ~ Depredation_Risk + Depredation_Risk:Directed + 
                                      Dist_Herbaceous + Dist_Herbaceous:Directed + 
                                      Dist_Shrub + Dist_Shrub:Directed + 
                                      Dist_Forest + Dist_Forest:Directed + 
                                      Elevation + Elevation:Directed + 
                                      Slope + Slope:Directed + 
                                      strata(StepID) + cluster(ID), 
                                    # Retain only steps in the 'Meandering' or 'Directed' 
                                    # movement state and rescale each continuous variable by
                                    # subtracting its mean and dividing by two standard deviations
                                    data = SSF_Data$`2`$`5` %>%
                                      filter(State %in% c("Meandering", "Directed")) %>% 
                                      mutate(across(c(Depredation_Risk:Slope), 
                                                    function(x){
                                                      (x - mean(x))/(2*sd(x))
                                                    })), 
                                    method = "approximate")

TOD_Specific_2_Hour_SSF <- clogit(Used ~ Depredation_Risk + Depredation_Risk:Night + 
                                    Dist_Herbaceous + Dist_Herbaceous:Night + 
                                    Dist_Shrub + Dist_Shrub:Night + 
                                    Dist_Forest + Dist_Forest:Night + 
                                    Elevation + Elevation:Night + 
                                    Slope + Slope:Night + 
                                    strata(StepID) + cluster(ID), 
                                  # Retain only steps in the 'Meandering' or 'Directed' 
                                  # movement state and rescale each continuous variable by
                                  # subtracting its mean and dividing by two standard deviations
                                  data = SSF_Data$`2`$`5` %>%
                                    filter(State %in% c("Meandering", "Directed")) %>% 
                                    mutate(across(c(Depredation_Risk:Slope), 
                                                  function(x){
                                                    (x - mean(x))/(2*sd(x))
                                                  })), 
                                  method = "approximate")

Table_S6 <- list(`No Interactions` = Main_SSF,
                `Time of Day x Resource Variables` = TOD_Specific_2_Hour_SSF, 
                `Sex x Resource Variables` = Sex_Specific_2_Hour_SSF, 
                `Movement State x Resource Variables` = State_Specific_2_Hour_SSF) %>% 
  map_df(~tibble(QIC = MuMIn::QIC(.x)), 
         .id = "Model") %>% 
  arrange(QIC) %>% 
  mutate(Delta_QIC = round(QIC - min(QIC), 2))


Table_S7 <- tidy(TOD_Specific_2_Hour_SSF) %>%
  left_join(confint(TOD_Specific_2_Hour_SSF) %>% 
              as_tibble(rownames = "term")) %>% 
  mutate(term = factor(term, 
                       levels = c("Depredation_Risk", 
                                  "Depredation_Risk:Night", 
                                  "Dist_Forest", 
                                  "Dist_Forest:Night", 
                                  "Dist_Herbaceous", 
                                  "Dist_Herbaceous:Night", 
                                  "Dist_Shrub", 
                                  "Dist_Shrub:Night", 
                                  "Elevation", 
                                  "Elevation:Night", 
                                  "Slope", 
                                  "Slope:Night"),
                       labels = c("Mortality Risk", 
                                  "Mortality Risk x Night", 
                                  "Distance to Forest", 
                                  "Distance to Forest x Night", 
                                  "Distance to Herbaceous", 
                                  "Distance to Herbaceous x Night", 
                                  "Distance to Shrub", 
                                  "Distance to Shrub x Night", 
                                  "Elevation", 
                                  "Elevation x Night", 
                                  "Slope", 
                                  "Slope x Night")), 
         across(c(estimate, `2.5 %`, `97.5 %`), ~round(.x, 2))) %>%
  select(term, estimate, `2.5 %`, `97.5 %`, p.value) %>% 
  arrange(term)

Table_S8 <- tidy(State_Specific_2_Hour_SSF) %>%
  left_join(confint(State_Specific_2_Hour_SSF) %>% 
              as_tibble(rownames = "term")) %>% 
  mutate(term = factor(term, 
                       levels = c("Depredation_Risk", 
                                  "Depredation_Risk:Directed", 
                                  "Dist_Forest", 
                                  "Dist_Forest:Directed", 
                                  "Dist_Herbaceous", 
                                  "Dist_Herbaceous:Directed", 
                                  "Dist_Shrub", 
                                  "Dist_Shrub:Directed", 
                                  "Elevation", 
                                  "Elevation:Directed", 
                                  "Slope", 
                                  "Slope:Directed"),
                       labels = c("Mortality Risk", 
                                  "Mortality Risk x Directed", 
                                  "Distance to Forest", 
                                  "Distance to Forest x Directed", 
                                  "Distance to Herbaceous", 
                                  "Distance to Herbaceous x Directed", 
                                  "Distance to Shrub", 
                                  "Distance to Shrub x Directed", 
                                  "Elevation", 
                                  "Elevation x Directed", 
                                  "Slope", 
                                  "Slope x Directed")), 
         across(c(estimate, `2.5 %`, `97.5 %`), ~round(.x, 2))) %>%
  select(term, estimate, `2.5 %`, `97.5 %`, p.value) %>% 
  arrange(term)

Table_S9 <- tidy(Sex_Specific_2_Hour_SSF) %>%
  left_join(confint(Sex_Specific_2_Hour_SSF) %>% 
              as_tibble(rownames = "term")) %>% 
  mutate(term = factor(term, 
                       levels = c("Depredation_Risk", 
                                  "Depredation_Risk:Male", 
                                  "Dist_Forest", 
                                  "Dist_Forest:Male", 
                                  "Dist_Herbaceous", 
                                  "Dist_Herbaceous:Male", 
                                  "Dist_Shrub", 
                                  "Dist_Shrub:Male", 
                                  "Elevation", 
                                  "Elevation:Male", 
                                  "Slope", 
                                  "Slope:Male"),
                       labels = c("Mortality Risk", 
                                  "Mortality Risk x Male", 
                                  "Distance to Forest", 
                                  "Distance to Forest x Male", 
                                  "Distance to Herbaceous", 
                                  "Distance to Herbaceous x Male", 
                                  "Distance to Shrub", 
                                  "Distance to Shrub x Male", 
                                  "Elevation", 
                                  "Elevation x Male", 
                                  "Slope", 
                                  "Slope x Male")), 
         across(c(estimate, `2.5 %`, `97.5 %`), ~round(.x, 2))) %>%
  select(term, estimate, `2.5 %`, `97.5 %`, p.value) %>% 
  arrange(term)



# Functional responses in resource selection -----------------------------------

## Create a list where each element is SSF data from an individual mountain lion with: 
##    - 2 Hour Interval Length
##    - 5 Available Steps 
Individual_SSF_Data <- SSF_Data$`2`$`5` %>% 
  group_by(ID) %>% 
  group_split() %>% 
  set_names(map(., ~unique(.x$ID)))

## Fit Individual SSFs 
Individual_SSF_Models <- map(Individual_SSF_Data, 
                             function(individual_ssf_data){
                               
                               model = clogit(Used ~ Depredation_Risk +
                                                Dist_Forest + Dist_Shrub + Dist_Herbaceous +
                                                Elevation + Slope + strata(StepID),
                                              data = individual_ssf_data %>% 
                                                filter(State %in% c("Meandering", "Directed")) %>% 
                                                mutate(across(c(Depredation_Risk:Slope), 
                                                              function(x){
                                                                (x - mean(x))/(2*sd(x))
                                                              })))
                               
                               model_summary = summary(model)$coefficients %>%
                                 as_tibble(rownames = "term") %>% 
                                 left_join(confint(model) %>% 
                                             as_tibble(rownames = "term"))
                               
                               model_summary
                               
                             })

## Get mean values of depredation mortality risk from each home range (used locations
## in the second order RSF data)
Mean_HR_Depredation_Risk <- Second_Order_RSF_Data %>% 
  filter(Used == 1) %>% 
  group_by(ID, Study, Male) %>% 
  summarise(Depredation_Risk_Mean = mean(Depredation_Risk))

## Create functional response data by adding values for mean depredation mortality 
## risk to each row
Functional_Response_Data <- Individual_SSF_Models %>% 
  bind_rows(.id = "ID") %>% 
  left_join(Mean_HR_Depredation_Risk) %>% 
  select(ID, Study, Male, term, Fine_Scale_Response = coef, SE = `se(coef)`, Depredation_Risk_Mean) %>% 
  filter(term %in% c("Depredation_Risk", "Dist_Forest", "Dist_Shrub", "Dist_Herbaceous"))

## Fit Functional Response Models for: 
##   - Depredation Mortality Risk
##   - Distance Forest
##   - Distance Shrub
##   - Distance Herbaceous
Depredation_Mortality_Risk_FR <- lm(Fine_Scale_Response ~ Depredation_Risk_Mean, 
                                    data = Functional_Response_Data %>% 
                                      filter(term == "Depredation_Risk"))

Forest_FR <- lm(Fine_Scale_Response ~ Depredation_Risk_Mean, 
                data = Functional_Response_Data %>% 
                  filter(term == "Dist_Forest"))

Shrub_FR <- lm(Fine_Scale_Response ~ Depredation_Risk_Mean, 
               data = Functional_Response_Data %>% 
                 filter(term == "Dist_Shrub"))

Herbaceous_FR <- lm(Fine_Scale_Response ~ Depredation_Risk_Mean, 
                    data = Functional_Response_Data %>% 
                      filter(term == "Dist_Herbaceous"))

## Summarize functional response models
summary(Depredation_Mortality_Risk_FR)
summary(Forest_FR)
summary(Shrub_FR)
summary(Herbaceous_FR)

## Fit Sex-Specific Functional Response Models
Sex_Specific_Depredation_Mortality_Risk_FR <- lm(Fine_Scale_Response ~ Depredation_Risk_Mean*Male, 
                                                 data = Functional_Response_Data %>% 
                                                   filter(term == "Depredation_Risk"))

Sex_Specific_Forest_FR <- lm(Fine_Scale_Response ~ Depredation_Risk_Mean*Male, 
                             data = Functional_Response_Data %>% 
                               filter(term == "Dist_Forest"))

Sex_Specific_Shrub_FR <- lm(Fine_Scale_Response ~ Depredation_Risk_Mean*Male, 
                            data = Functional_Response_Data %>% 
                              filter(term == "Dist_Shrub"))

Sex_Specific_Herbaceous_FR <- lm(Fine_Scale_Response ~ Depredation_Risk_Mean*Male, 
                                 data = Functional_Response_Data %>% 
                                   filter(term == "Dist_Herbaceous"))

AICc(Depredation_Mortality_Risk_FR, 
     Sex_Specific_Depredation_Mortality_Risk_FR) %>% 
  mutate(Delta_AICc = round(AICc - min(AICc), 2))

## Evaluate if observed functional response is robust to uncertainty around 
## individual-level beta coefficients
set.seed(1)

## Generate 10,000 samples from a normal distribution centered on the value 
## of each individual's beta coefficient
Depredation_Risk_FR_Evaluation_Data <- map(1:10000, 
                                           function(x){
                                             
                                             Functional_Response_Data %>% 
                                               filter(term == "Depredation_Risk") %>% 
                                               rowwise() %>% 
                                               mutate(Sampled_Estimate = rnorm(1, 
                                                                               Fine_Scale_Response, 
                                                                               SE))
                                             
                                           }) %>% 
  set_names(1:10000)

Depredation_Risk_FR_Evaluation_Data_Summary <- Depredation_Risk_FR_Evaluation_Data %>%
  bind_rows(.id = "Sample") %>%
  group_by(ID, Depredation_Risk_Mean) %>%
  summarise(Mean_Fine_Scale_Response = c(mean(Sampled_Estimate), mean(Sampled_Estimate)),
            HPD_Level = c("90%", "95%"),
            L_HPD_Fine_Scale_Response = c(empiricalHpd(Sampled_Estimate, 0.9)[1], empiricalHpd(Sampled_Estimate, 0.95)[1]),
            U_HPD_Fine_Scale_Response = c(empiricalHpd(Sampled_Estimate, 0.9)[2], empiricalHpd(Sampled_Estimate, 0.95)[1]))

## Fit the observed functional response with each of the 10,000 datasets and create 
## a dataframe that retains values for the intercept and depredation mortality risk 
## beta coefficient for each sample
Depredation_Risk_FR_Evaluation_Coefficients <- map(Depredation_Risk_FR_Evaluation_Data, 
                                                   function(functional_response_data){
                                                     
                                                     fr_model = lm(Sampled_Estimate ~ Depredation_Risk_Mean,
                                                                   data = functional_response_data)
                                                     
                                                     summary(fr_model)$coefficients %>% 
                                                       as_tibble(rownames = "term")
                                                     
                                                   }) %>% 
  bind_rows(.id = "Sample") %>% 
  select(Sample, term, Estimate) %>% 
  pivot_wider(names_from = term, 
              values_from = Estimate)

## Evaluates whether the observed functiona response is robust to uncertainty by 
## calculating 90% and 95% emperical highest posteriod density (HPD) intervals 
Depredation_Risk_FR_Evaluation_Summary <- Depredation_Risk_FR_Evaluation_Coefficients %>% 
  summarise(Mean_Intercept = c(mean(`(Intercept)`), mean(`(Intercept)`)), 
            HPD_Level = c("90%", "95%"),
            Intercept_L_HPD = c(empiricalHpd(`(Intercept)`, 0.90)[1], empiricalHpd(`(Intercept)`, 0.95)[1]),
            Intercept_U_HPD = c(empiricalHpd(`(Intercept)`, 0.90)[2], empiricalHpd(`(Intercept)`, 0.95)[2]),
            Depredation_Risk_Mean_Coeff = c(mean(Depredation_Risk_Mean), mean(Depredation_Risk_Mean)),
            Depredation_Risk_Mean_L_HPD = c(empiricalHpd(Depredation_Risk_Mean, 0.90)[1], empiricalHpd(Depredation_Risk_Mean, 0.95)[1]),
            Depredation_Risk_Mean_U_HPD = c(empiricalHpd(Depredation_Risk_Mean, 0.90)[2], empiricalHpd(Depredation_Risk_Mean, 0.95)[2])) %>% 
  mutate(across(where(is.numeric), ~round(.x, 2)))

## Figure 4
Ribbon = tibble(x = seq(0, 0.5, length.out = 10000), 
                Lower = Depredation_Risk_FR_Evaluation_Summary$Intercept_L_HPD[2] + 
                  (x * Depredation_Risk_FR_Evaluation_Summary$Depredation_Risk_Mean_L_HPD[2]),
                Mean = Depredation_Risk_FR_Evaluation_Summary$Mean_Intercept[2] + 
                  (x * Depredation_Risk_FR_Evaluation_Summary$Depredation_Risk_Mean_Coeff[2]),
                Upper = Depredation_Risk_FR_Evaluation_Summary$Intercept_U_HPD[2] + 
                  (x * Depredation_Risk_FR_Evaluation_Summary$Depredation_Risk_Mean_U_HPD[2]))

ggplot(Depredation_Risk_FR_Evaluation_Data_Summary, 
       aes(x = Depredation_Risk_Mean, 
           y = Mean_Fine_Scale_Response)) + 
  geom_hline(yintercept = 0,
             lty = 2,
             color = "red") + 
  geom_ribbon(data = Ribbon, 
              aes(x = x, 
                  ymin = Lower, 
                  ymax = Upper),
              inherit.aes = FALSE, 
              fill = "blue",
              color = NA,
              alpha = 0.5) + 
  geom_line(data = Ribbon, 
            aes(x = x, 
                y = Mean), 
            inherit.aes = FALSE, 
            color = "blue") + 
  geom_errorbar(aes(ymin = L_HPD_Fine_Scale_Response,
                    ymax = U_HPD_Fine_Scale_Response),
                color = "black",
                alpha = 0.5, 
                width = 2) + 
  geom_point(color = "black") +
  scale_x_continuous(limits = c(0, 0.5),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3, 1.5),
                     expand = c(0, 0)) +
  labs(x = "Mean Mortality Risk",
       y = "Beta Coefficient") +#"Fine-Scale Response to Mortality Risk") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0,
                                   vjust = 0.5,
                                   color = "black",
                                   size = 10,
                                   face = "bold"),
        axis.text.y = element_text(color = "black",
                                   size = 10,
                                   face = "bold"),
        axis.title.y = element_text(color = "black",
                                    size = 14,
                                    face = "bold"),
        axis.title.x = element_text(color = "black",
                                    size = 14,
                                    face = "bold"),
        plot.tag = element_text(face = "bold",
                                size = 14),
        plot.tag.position = c(0.1, 0.98), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"))

### Evaluate correlation between second order random slope coefficients 
### and mean values of depredation mortality risk within an individuals home range
Second_Order_Random_Slope_Coefficients <- readRDS("Results/Second_Order_RSF_Models/Random_Slope_Second_Order_RSF.RDS") %>%
  coef()

Table_S10 <- Second_Order_Random_Slope_Coefficients$cond[[1]] %>%
  as_tibble(rownames = "ID") %>% 
  left_join(Mean_HR_Depredation_Risk) %>% 
  select(ID, Study, Coarse_Scale_Response = Depredation_Risk, Depredation_Risk_Mean) %>% 
  group_by(Study) %>% 
  filter(n() > 3) %>% 
  summarise(Corr = round(cor(Coarse_Scale_Response, Depredation_Risk_Mean), 2),
            P = round(cor.test(Coarse_Scale_Response, Depredation_Risk_Mean)$p.value, 3),
            N = n()) %>% 
  mutate(P = case_when(P <= 0.000 ~ "<0.001", 
                       TRUE ~ as.character(P))) %>% 
  arrange(desc(N))


# Supporting Figures ------------------------------------------------------

## Figure S1 ----
### Fit SSF Reported in Figure S1 
###    - Interval Length = 2 Hours, 
###    - n available steps = 1-10, 15, and 20
Variable_N_Available_SSFs <- map(SSF_Data$`2`, 
                                 function(ssf_data_n){
                                   
                                   ssf_model = clogit(Used ~ Depredation_Risk + Dist_Herbaceous + 
                                                        Dist_Shrub + Dist_Forest + Elevation + Slope + 
                                                        strata(StepID) + cluster(ID), 
                                                      # Retain only steps in the 'Meandering' or 'Directed' 
                                                      # movement state and rescale each continuous variable by
                                                      # subtracting its mean and dividing by two standard deviations
                                                      data = ssf_data_n %>%
                                                        filter(State %in% c("Meandering", "Directed")) %>% 
                                                        mutate(across(c(Depredation_Risk:Slope), 
                                                                      function(x){
                                                                        (x - mean(x))/(2*sd(x))
                                                                      })), 
                                                      method = "approximate")
                                   
                                   tidy(ssf_model) %>%
                                     left_join(confint(ssf_model) %>% 
                                                 as_tibble(rownames = "term")) %>% 
                                     mutate(across(c(estimate, `2.5 %`, `97.5 %`), ~round(.x, 2)), 
                                            N_Available = unique(ssf_data_n$N_Available)) %>%
                                     select(term, estimate, `2.5 %`, `97.5 %`, p.value, N_Available)
                                   
                                 })

ggplot(data = Variable_N_Available_SSFs %>% 
         bind_rows() %>%
         mutate(term = factor(term, 
                              levels = c("Depredation_Risk", 
                                         "Dist_Forest", 
                                         "Dist_Shrub", 
                                         "Dist_Herbaceous", 
                                         "Elevation", 
                                         "Slope"),
                              labels = c("Mortality Risk", 
                                         "Distance to Forest", 
                                         "Distance to Shrub", 
                                         "Distance to Herbaceous", 
                                         "Elevation", 
                                         "Slope"))), 
       aes(x = N_Available, 
           y = estimate, 
           ymin = `2.5 %`, 
           ymax = `97.5 %`)) + 
  geom_hline(yintercept = 0, 
             color = "red", 
             lty = 2) + 
  labs(x = "Number of Available Steps", 
       y = "Beta Coefficient") + 
  geom_errorbar() + 
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0,
                                   vjust = 0.5,
                                   color = "black", 
                                   size = 10, 
                                   face = "bold"),
        axis.text.y = element_text(color = "black", 
                                   size = 10, 
                                   face = "bold"),
        axis.title.y = element_text(color = "black", 
                                    size = 14, 
                                    face = "bold"),
        axis.title.x = element_text(color = "black", 
                                    size = 14, 
                                    face = "bold"),
        plot.tag = element_text(face = "bold", 
                                size = 14), 
        plot.tag.position = c(0.075, 0.98), 
        strip.text.x = element_text(size = 10, 
                                    face = "bold")) +
  facet_wrap(~term)

## Figure S2 ----
Depredation_Mortality_Risk_Models %>% 
  bind_rows() %>%
  group_by(HR_Sample) %>% 
  filter(QIC == min(QIC)) %>% 
  unnest(Summary) %>%
  mutate(term = factor(term, 
                       levels = c("Dist_Cover", 
                                  "Dist_Developed_Low", 
                                  "Dist_Local_Roads", 
                                  "Distance_Building", 
                                  "slope"), 
                       labels = c("Distance to \nCover", 
                                  "Distance to \nLow-Intensity \nDevelopment", 
                                  "Distance to \nLocal Roads", 
                                  "Distance to \nBuildings", 
                                  "Slope")), 
         HR_Sample = as.numeric(str_remove(HR_Sample, "HR_Sample_"))) %>%
  arrange(HR_Sample) %>%
  ggplot(aes(x = term, 
             y = estimate, 
             ymin = `lower .95`, 
             ymax = `upper .95`)) +
  geom_hline(yintercept = 0, 
             col = "red", 
             lty = 2) + 
  geom_errorbar(position = position_dodge2(width = 0.5), 
                width = 0.5,
                alpha = 0.25) + 
  geom_point(aes(color = HR_Sample), 
             size = 2, 
             position = position_dodge2(width = 0.5)) + 
  labs(y = "Beta Coefficients",
       color = "Availability Sample") +
  theme_classic() + 
  theme(axis.text = element_text(angle = 0, 
                                 vjust = 0.5, 
                                 size = 11,
                                 face = "bold",
                                 color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14,
                                    face = "bold", 
                                    color = "black"),
        legend.position = c(0.15, 0.25),
        legend.title = element_text(face = "bold"))

## Figure S3 ----
### Fit SSF Reported in Figure S3
###    - Interval Lengths = 1-4 Hours, 
###    - n available steps = 5
Variable_Interval_Length_SSFs <- map(SSF_Data, 
                                     function(ssf_data){
                                       
                                       ssf_model = clogit(Used ~ Depredation_Risk + Dist_Herbaceous + 
                                                            Dist_Shrub + Dist_Forest + Elevation + Slope + 
                                                            strata(StepID) + cluster(ID), 
                                                          # Retain only steps in the 'Meandering' or 'Directed' 
                                                          # movement state and rescale each continuous variable by
                                                          # subtracting its mean and dividing by two standard deviations
                                                          data = ssf_data$`5` %>%
                                                            filter(State %in% c("Meandering", "Directed")) %>% 
                                                            mutate(across(c(Depredation_Risk:Slope), 
                                                                          function(x){
                                                                            (x - mean(x))/(2*sd(x))
                                                                          })), 
                                                          method = "approximate")
                                       
                                       tidy(ssf_model) %>%
                                         left_join(confint(ssf_model) %>% 
                                                     as_tibble(rownames = "term")) %>% 
                                         mutate(across(c(estimate, `2.5 %`, `97.5 %`), ~round(.x, 2)), 
                                                Interval_Length = unique(ssf_data$`5`$Interval_Length)) %>%
                                         select(term, estimate, `2.5 %`, `97.5 %`, p.value, Interval_Length)
                                       
                                     })

ggplot(data = Variable_Interval_Length_SSFs %>% 
         bind_rows() %>%
         mutate(Interval_Length = factor(Interval_Length), 
                term = factor(term, 
                              levels = c("Depredation_Risk", 
                                         "Dist_Forest", 
                                         "Dist_Shrub", 
                                         "Dist_Herbaceous", 
                                         "Elevation", 
                                         "Slope"),
                              labels = c("Mortality Risk", 
                                         "Distance to Forest", 
                                         "Distance to Shrub", 
                                         "Distance to Herbaceous", 
                                         "Elevation", 
                                         "Slope"))), 
       aes(x = term, 
           y = estimate, 
           ymin = `2.5 %`, 
           ymax = `97.5 %`, 
           color = Interval_Length)) + 
  geom_hline(yintercept = 0, 
             col = "red", 
             lty = 2) + 
  geom_errorbar(position = position_dodge2(width = 0.5), 
                width = 0.5) + 
  geom_point(size = 2, 
             position = position_dodge2(width = 0.5)) + 
  scale_color_manual(values = c("1" = "blue", 
                                "2" = "black", 
                                "3" = "darkgreen", 
                                "4" = "red")) + 
  labs(y = "Beta Coefficients", 
       color = "Interval Length (hours)") + 
  theme_classic() + 
  theme(axis.text = element_text(angle = 0, 
                                 vjust = 0.5, 
                                 size = 11,
                                 face = "bold",
                                 color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14,
                                    face = "bold", 
                                    color = "black"),
        legend.position = c(0.15, 0.2),
        legend.title = element_text(face = "bold"))

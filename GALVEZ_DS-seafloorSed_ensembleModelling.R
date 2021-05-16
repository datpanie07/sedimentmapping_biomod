##############################################################################################
###### This script is a supplementary material of the research article entitled: Ensemble mapping
#### and change analysis of sediment distributions in the Sylt Outer Reef, German North Sea from 2016-2018
#### Galvez, D.S. et al., 2021#####
##### script author:  Galvez,D.S.; written in May 2021
#############################################################################################
## set working directory
setwd("E:/H3/2018")

## load libraries.use install.packages() to install the packages 
## if they are not yet available in your library
library(biomod2) 
library(sdm) 
library(rgdal) 
library(ggplot2)
library(raster) 
library(rasterVis)
library(rpart)
library(gridExtra)
library(inlmisc)
library(dplyr)
library(dismo)
library(tidyr)
library(usdm)

##########################################################################################
########################### DATA PREPARATION #############################################
##########################################################################################

### Load response variables--sediment classes or species that you would like to model
### Here the files are in shapefile format, where the attribute tables contain:
#####the coordinates of the samples and a column stating the presence/absence (1-0) of the response variable

lag18<- readOGR(dsn = "E:/H3/2018/shp_utm", layer = "H32018LagSed_presence_utm2")
sand18 <- readOGR(dsn = "E:/H3/2018/shp_utm", layer = "H32018Sand_presence_utm2")

### Load the predictor variable--textural features or environmental layers (e.g. bathy, slope, etc.)

#### Textural features from GLCM

tex5 <-list.files(path='E:/H3/2018/glcm1',pattern='.tif$', full.names=T)
tex5 <-stack(tex5)

tex10 <-list.files(path='E:/H3/2018/glcm2',pattern='.tif$', full.names=T)
tex10 <-stack(tex10)

### Load geophysical/environmental predictor variables
gp <- list.files(path='E:/H3/2018/geophysical',pattern='.tif$', full.names=T)
gp <- stack(gp)

### Resample tex resolution to the resolution of gp. Only do this if the resolution
#### or spatial extents of your variables do not match
tex5 <- resample(tex5, gp)
tex10 <- resample(tex10, gp)

### Merge the tex and gp variables into a RasterStack
pred <- stack(tex5, tex10, gp)

### check if the layers overlap
plot(pred$sss)
plot(lag18, add=T)

##############################################################################
########### FEATURE SELECTION: VIF TEST FOR MULTICOLLINEARITY ################
##############################################################################
## The VIF test is to evaluate if there are multi-collinearity or high correlation between your predictors
### the following steps will identify predictors with high correlation value and 
#### remove them from your set of predictor variables
#####

### Extract the raster values of your predictors in the location of your sampling points/response variables
lag.ex <- raster::extract(pred,lag18)
sand.ex <- raster::extract(pred, sand18)

### Assess the variable for multi-collinearity; you can use either 'vifcor' or 'vifstep'
vif(pred)

v1 <- vifstep(lag.ex, th = 5) #81 out of 88 variables have collinearity prob
v1

v2 <- vifstep(sand.ex, th = 5)#78 out of 88 variables have collinearity prob
v2

### Exclude the variables from your set of predictors that were suggested by vif

pred.fsp <- exclude(pred, v1)
pred.lsp <- exclude(pred, v2)

# #### Optional: You can export the final predictors into your working directory, so you will not
# #### need to repeat the feature selection step when you re-run the model
# ExportRasterStack(pred.fsp, path="E:/H3/2018/pred_fsp")
# ExportRasterStack(pred.lsp, path="E:/H3/2018/pred_lsp")


###################################################################################
##################### BIOMOD MODELLING STEP 1: FORMAT DATA ########################
###################################################################################

### Load the predictor variables from the VIF test; Please note that each sediment class 
##### usually require different sets of predictor variables, therefore you cannot use the same predictors
###### for two different sediment classes. 

pred.lsp <- list.files(path='E:/H3/2018/pred_lsp',pattern='.tif$', full.names=T)
pred.fsp <- list.files(path='E:/H3/2018/pred_fsp',pattern='.tif$', full.names=T)

pred.lsp <- stack(pred.lsp)
pred.fsp <- stack(pred.fsp)

# ### Optional: rename the predictor variables--for easier reference
# 
# names(pred.lsp) <- c("glcm59", "glcm64", "glcm68",
#                      "bathy","eastness", "max.friction", "northness",
#                      "res.velV", "slope","SSS")
# 
# 
# names(pred.fsp) <- c( "glcm71", "glcm72", "glcm80", "eastness", 
#                       "northness", "res.velU" , "SSS")

###### Subset variables that have high variable importance score-- This step can only be done 
###### after the first model run, because you will need to get the variable importance score

# pred.lsp1 <- stack(pred.lsp$eastness,
#                    pred.lsp$glcm68,
#                    pred.lsp$glcm64,
#                    pred.lsp$SSS)
#                    
# 
# pred.fsp1 <- stack(pred.fsp$SSS,
#                    pred.fsp$glcm71)
                 

##### Format  your response and explanatory variables as BIOMOD data
set.seed(555) ##set seed to generate the same datasets in every model runs
data_lag18<- BIOMOD_FormatingData(resp.var = lag18['LagSed'],##name of your response variable in the column
                                   expl.var = pred.lsp, ## predictor variables--rasterStack
                                   resp.name = 'LagSed',
                                   PA.nb.rep = 3, ## pseudo-absences will be generated 3x in 3 different random locations
                                   PA.nb.absences = 200,## the no. of pseudo-absences will depend on the no. of your sample pts
                                   PA.strategy = 'random') ## pseudo-absences strategy 

set.seed(555)
data_sand18<- BIOMOD_FormatingData(resp.var = sand18['fSa'],
                                  expl.var = pred.fsp1,
                                  resp.name = 'Sand',
                                  PA.nb.rep = 3,
                                  PA.nb.absences = 200,
                                  PA.strategy = 'random')



###################################################################################
##################### BIOMOD MODELLING STEP 2: SET INDIVIDUAL MODEL SETTINGS ######
###################################################################################

##### For each sediment class, you have to set the parameters of the models that you would like
###### to use for the ensemble models. Settings need to be set accordingly during model calibration/fitting

### Model settings for Lag sediment models

set.seed(555)
Optmod_lag <- BIOMOD_ModelingOptions(
  GBM = list(distribution= "bernoulli", ntrees =5000, shrinkage= 0.001, 
             interaction.depth = 5, bag.fraction = 0.5, cv.folds=10),
  CTA= list(method="class", control=rpart.control(minsplit= 1, xval=10, cp = 0.01),
   RF  = list(ntree = 1000, type = "classification", importance = TRUE)),
  ANN = NULL)

###  Model settings for sand models
set.seed(555)
Optmod_sand <- BIOMOD_ModelingOptions(
  GBM = list(distribution= "bernoulli", ntrees =5000, shrinkage= 0.001, 
             interaction.depth= 5,bag.fraction = 0.5, cv.folds=10),
  CTA= list(method= "class", control= rpart.control(cp = 0.05,minsplit= 5)),
  RF  = list(ntree = 2500, do.classif= TRUE, importance = TRUE),
  ANN = NULL)
# ANN = list (size = 2, rang= 0.5, decay= 5e-04, maxit = 300))


##################################################################################
##################### BIOMOD MODELLING STEP 3: MODEL CALIBRATION #################
##################################################################################
### In this step, you need to run and re-run the models until you achieve a desirable
#### model score performance (e.g. TSS > 0.7)
####### Here, each sediment class are modelled separately instead of running in loop
######## because each sediment class require a different set of predictors and model parameter settings


## LagSed single model runs
set.seed(555)
mod.lag <- BIOMOD_Modeling( data = data_lag18,## data in biomod format
                              models = c( 'CTA','GBM', 'RF', 'ANN'), ## model choice
                              models.options = Optmod_lag, ## individual model settings
                              NbRunEval = 20, ## no. of evaluation runs
                              DataSplit = 70, ## data split for training and testing; 70-30%
                              Yweights = NULL,
                              VarImport = 10, ## no. of permutations to test variable importance
                              models.eval.meth = c('TSS', 'ROC', 'KAPPA'),
                              SaveObj = TRUE,
                              rescal.all.models = FALSE,
                              do.full.models = TRUE,
                              modeling.id = "lag_h32018")


## Sand single model runs

set.seed(555)
mod.sand <- BIOMOD_Modeling( data = data_fsp18,
                              models = c('GBM', 'RF', 'ANN','CTA'),
                              models.options = Optmod_sand,
                              NbRunEval = 20,
                              DataSplit = 70,
                              Yweights = NULL,
                              VarImport = 10,
                              models.eval.meth = c('TSS', 'ROC','KAPPA'),
                              SaveObj = TRUE,
                              rescal.all.models = FALSE,
                              do.full.models = FALSE,
                              modeling.id = "sand_h32018")    


###################################################################################
##################### BIOMOD MODELLING STEP 4: SINGLE MODEL EVALUATION ############
###################################################################################
### Here, you can assess your model performance using different graphs

######## Plot model scores according to model algorithms,
##########number of cross-validation runs, and datasets. 
########### Models scores are ploted based on their TSS and ROC values

lag_bymod <- models_scores_graph(mod.lag,
                    by = "models",
                    metrics = c("ROC", "TSS"),
                    xlim= c(0,1), ylim= c(0,1))
lag_byCVrun <-models_scores_graph(mod.lag, 
                    by = "cv_run",
                    metrics = c("ROC", "TSS"),
                    xlim= c(0,1), ylim= c(0,1))
lag_byDataset <- models_scores_graph(mod.lag,
                    by = "data_set",
                    metrics = c("ROC", "TSS"),
                    xlim= c(0,1), ylim= c(0,1))

grid.arrange(lsp18_bymod, lsp18_byCVrun, lsp18_byDataset) ## arrange scores into one image


####### Plot variable responses to each model. 
####### Here you can evaluate which predictors influence your model performance
 
lag_cta <- BIOMOD_LoadModels(mod.lag, models= 'CTA')
lag_gbm <- BIOMOD_LoadModels(mod.lag, models= 'GBM')
lag_rf <- BIOMOD_LoadModels(mod.lag, models= 'RF')
lag_ann <- BIOMOD_LoadModels(mod.lag, models= 'ANN')

lag_cta_rp <- biomod2::response.plot2(
  models = lag_cta,
  Data = get_formal_data(mod.lag, 'expl.var'),
  show.variables = get_formal_data(mod18_lsp, 'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title =TRUE,
  plot = TRUE)

lag_gbm_rp <- biomod2::response.plot2(
  models = lag_gbm,
  Data = get_formal_data(mod.lag, 'expl.var'),
  show.variables = get_formal_data(mod18_lsp, 'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title =TRUE,
  plot = TRUE)

lag_rf_rp <- biomod2::response.plot2(
  models = lag_rf,
  Data = get_formal_data(mod.lag, 'expl.var'),
  show.variables = get_formal_data(mod.lag, 'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title =TRUE,
  plot = TRUE)

lag_ann_rp <- biomod2::response.plot2(
  models = lag_ann,
  Data = get_formal_data(mod.lag, 'expl.var'),
  show.variables = get_formal_data(mod.lag, 'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title =TRUE,
  plot = TRUE)


####### Check the variable importance of your predictors ######
#### This step is important to identify which predictors you need to remove from your
###### original set of predictors; using predictors with variable importance score of >0.1 can improve your models

lsp18_VI <- get_variables_importance(mod18_lsp)
fsp18_VI <- get_variables_importance(mod_fsp18)

apply(lsp18_VI, c(1,2), mean) ## make the mean of variable importance by algorithm
apply(fsp18_VI, c(1,2), mean)


###################################################################################
##################### BIOMOD MODELLING STEP 5: ENSEMBLE MODELLING #################
###################################################################################

#### After you achieved a satisfying performance score of your single models, 
##### you can now build your ensemble models

### Ensemble model for lag sediment
lag.EM <- BIOMOD_EnsembleModeling(
              modeling.output = mod.lag,
              chosen.models = 'all',
              em.by = 'all', 
              eval.metric = c('TSS'),  ### select which metric to use to filter the models that you want to include in the EM  
              eval.metric.quality.threshold = c(0.7), ## only models with TSS score of at least 0.7 will be included
              models.eval.meth = c('TSS', 'ROC', 'KAPPA'), 
              prob.mean = TRUE, 
              prob.cv = TRUE,
              prob.ci = FALSE, 
              prob.ci.alpha = 0.05, 
              prob.median = FALSE, 
              committee.averaging = TRUE,
              prob.mean.weight = TRUE,
              prob.mean.weight.decay = 'proportional')

### Ensemble model for sand
sand.EM <- BIOMOD_EnsembleModeling(
                    modeling.output = mod.sand,
                    chosen.models = 'all',
                    em.by = 'all', 
                    eval.metric = c('TSS'),
                    eval.metric.quality.threshold = c(0.7),
                    models.eval.meth = c('TSS', 'ROC', 'KAPPA'), 
                    prob.mean = TRUE, 
                    prob.cv = TRUE,
                    prob.ci = FALSE, 
                    prob.ci.alpha = 0.05, 
                    prob.median = FALSE, 
                    committee.averaging = TRUE,
                    prob.mean.weight = TRUE,
                    prob.mean.weight.decay = 'proportional')

###############################################
##### Evaluate the ensemble models ############
###############################################

#### Load score for each ensemble models

lagEM_scores <- get_evaluations(lag.EM)
sandEM_scores <- get_evaluations(sand.EM)

lagEM_scores
sandEM_scores


### Get the name of models that were kept in the final ensemble

lagEM_kept <- get_kept_models(lag.EM,.__C__.name)
sandEM_kept <-get_kept_models(sand.EM,.__C__.name)

summary(lagEM_kept)
summary(sandEM_kept)


###################################################################################
##################### BIOMOD MODELLING STEP 6: SINGLE MODELS PROJECTION ###########
###################################################################################

###### Project the sediment distribution using the single models

sp_lag <- BIOMOD_Projection(
                      modeling.output = mod.lag,
                      new.env = pred.lsp1, ## make sure to select the same predictors that you used in step 2
                      proj.name = "lag2018", 
                      selected.models = 'all', 
                      binary.meth = "TSS",
                      filtered.meth = NULL, 
                      compress = TRUE,
                      build.clamping.mask = TRUE,
                      do.stack = FALSE,
                      output.format = ".img")


sp_sand <- BIOMOD_Projection(
                      modeling.output = mod.sand,
                      new.env = pred.fsp1,
                      proj.name = "sand2018", 
                      selected.models = 'all', 
                      binary.meth = "TSS",
                      filtered.meth = NULL, 
                      compress = TRUE,
                      build.clamping.mask = TRUE,
                      do.stack = FALSE,
                      output.format = ".img")

### get the name of projected models
lag_sp <- get_projected_models(sp_lag)
sand_sp <- get_projected_models(sp_sand)

summary(lag_sp)
summary(sand_sp)

####################################################################################
##################### BIOMOD MODELLING STEP 7: ENSEMBLE MODEL PROJECTION ###########
####################################################################################

###### Project the sediment distribution using the ensemble model
lag.EMproj <- BIOMOD_EnsembleForecasting(
                            EM.output = lag.EM,
                            projection.output = sp_lag, 
                            binary.meth = "TSS",
                            compress = TRUE, 
                            do.stack = FALSE,
                            output.format = ".img")

sand.EMproj <- BIOMOD_EnsembleForecasting(
                            EM.output = sand.EM,
                            projection.output = sp_sand, 
                            binary.meth = "TSS",
                            compress = TRUE, 
                            do.stack = FALSE,
                            output.format = ".img")



###################################################################################
##################### BIOMOD MODELLING STEP 8: PLOT PROJECTIONS ###################
###################################################################################

### The output from the single and ensemble projections area saved in your working directory as .img or rasters
##### You can view the results here or in your GIS software


## Plot lag sediment projection
lag.plot <- lag.EMproj
lag.plot@models.projected <- c("mean", "cv", "ca","Wmean")
    plot(lsp18_Proj.plot, str.grep = "ca|Wmean")
    plot(lsp18_Proj.plot, str.grep = "mean|cv")

lag_levplot <- get_predictions(lag.EMproj)
names(lag_levplot) <- c("mean", "cv", "ca","Wmean")
levelplot(lsp18_levplot, par.settings = BuRdTheme, 
          main = "Lag Sediment ensemble projections\nin 2018")

## Plot sand projection
sand.plot <- sand.EMproj
sand.plot@models.projected <- c("mean", "cv", "ca","Wmean")
plot(sand.plot, str.grep = "ca|Wmean")


sand.levplot <- get_predictions(sand.EMproj)
names(sand.levplot) <- c("mean", "cv", "ca","Wmean")
levelplot(sand.levplot, main = "Fine sand ensemble projections\nin 2018",
          par.settings = BuRdTheme)

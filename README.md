# Interpretability of Machine Learning Models

Package in R to interpret predictions from machine learning models in a local and agnostic way.

The package is based on __Shimizu, G. Y., Izbicki, R., e Carvalho, A.C.: Model interpretation using improved local regression with variable importance (2022)__. 

Two interpretability methods are available: VarImp and SupClus. While VarImp can be used in more complex datasets (eg with non-linear relationships), SupClus can be used in simpler datasets generating explanations for clusters of instances with similar interpretations.

## Installation

Use the devtools package to install interpretability directly from github.
```
install_github("gilsonshimizu/interpretability")
```

## Simple example

```
#Install package
library(devtools)
install_github("gilsonshimizu/interpretability")
library(interpretability)

#Load artificial dataset
data("artificial_data")

#VarImp explanations 
x_coef <- varimp(x = artificial_data$x_test,
                 y = artificial_data$y_pred,
                 h = 0.1,
                 dmax = ncol(artificial_data$x_test))
#SupClus explanations
supcluscoef <- supclus(k = 2,
                       x = artificial_data$x_test,
                       y = artificial_data$y_pred,
                       h = 0.1,
                       dmax = ncol(artificial_data$x_test))

#Plot VarImp local slopes of some points for one feature.
plotslopes(sampsize = 50,
           varnumber = 1,
           x = artificial_data$x_test,
           y = artificial_data$y_pred,
           xcoef = x_coef)

#Plot VarImp effects of one instance.
ploteffect(instancen = 1,
           x = artificial_data$x_test,
           y = artificial_data$y_pred,
           xcoef = x_coef)

#Plot SupClus effects of one instance.
ploteffectclus(instancen = 1,
               x = artificial_data$x_test,
               y = artificial_data$y_pred, 
               clustercoef = supcluscoef$coef,  
               cluster = supcluscoef$cluster)
```

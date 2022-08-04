#' VarImp explanations
#' 
#' This function calculates the local coefficients using the VarImp method for a sample.
#' 
#' @param x Matrix \eqn{n \times d} with covariates for each sample unit.
#' @param y \eqn{n}-dimensional vector with predictions from the original model.
#' @param h Bandwidth of gaussian kernel.
#' @param dmax Maximum number of covariates (\eqn{\le d}).
#' @export
#' @return Matrix \eqn{n \times dmax} with local coefficients using the VarImp method.
#' @examples
#' library(interpretability)
#' 
#' data("artificial_data")
#' 
#' x_coef <- varimp(x = artificial_data$x_test,
#'                  y = artificial_data$y_pred,
#'                  h = 0.1,
#'                  dmax = ncol(artificial_data$x_test))

varimp <- function(x, y, h, dmax){
  varimportance <- varimpcal(x = x, y = y)$localimp
  cvarimp <- t(apply(x, 1, function(arg) locstep(arg, x = x, y = y, h = h, dmax = dmax, varimp = varimportance)))
  return(cvarimp)
}

#' SupClus explanations
#' 
#' This function calculates the local coefficients using the SubClus method for a sample.
#' 
#' @param k Number of interpretable clusters.
#' @param x Matrix \eqn{n \times d} with covariates for each sample unit.
#' @param y \eqn{n}-dimensional vector with predictions from the original model.
#' @param h Bandwidth of gaussian kernel.
#' @param dmax Maximum number of covariates (\eqn{\le d}).
#' @export
#' @return Matrix \eqn{n \times dmax} with local coefficients and \eqn{n}-dimensional vector with cluster of each instance. 
#' @examples
#' library(interpretability)
#' 
#' data("artificial_data")
#' 
#' supcluscoef <- supclus(k = 2,
#'                        x = artificial_data$x_test,
#'                        y = artificial_data$y_pred,
#'                        h = 0.1,
#'                        dmax = ncol(artificial_data$x_test))


supclus <- function(k, x, y, h, dmax){
  x_coef <- varimp(x = x, y = y, h = h, dmax = dmax)
  cluster_coef <- regcluster(k = k, xcoef = x_coef, x = x, y = y, dmax = dmax)
  csupclus <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  csupclus <- cluster_coef$coefs[cluster_coef$cluster, -ncol(cluster_coef$coefs)]
  return(list(coef = csupclus, cluster = cluster_coef$cluster))
}

#' Plot slopes
#'
#' This function plot local slopes of some points for one feature.
#' 
#' @param sampsize Sample size for plotting slopes.
#' @param varnumber Covariate index.
#' @param x Matrix \eqn{n \times d} with covariates for each sample unit.
#' @param y \eqn{n}-dimensional vector with predictions from the original model.
#' @param xcoef Matrix \eqn{n x d} with local slopes for each sample unit.
#' @export
#' @return Graph with local slopes of some points for one feature.
#' @examples
#' library(interpretability)
#' 
#' data("artificial_data")
#' 
#' x_coef <- varimp(x = artificial_data$x_test,
#'                  y = artificial_data$y_pred,
#'                  h = 0.1,
#'                  dmax = ncol(artificial_data$x_test))
#' 
#' plotslopes(sampsize = 50,
#'            varnumber = 1,
#'            x = artificial_data$x_test,
#'            y = artificial_data$y_pred,
#'            xcoef = x_coef)
#'            


plotslopes <- function(sampsize, varnumber, x, y, xcoef){ 
  xcoef <- as.matrix(xcoef[ , -ncol(xcoef)]) #remove intercept
  i <- sample(1:nrow(x), sampsize)
  slopes <- data.frame(intercept = y[i] - xcoef[i, varnumber] * x[i, varnumber], 
                       slope = xcoef[i, varnumber],
                       x = x[i, varnumber],
                       y = y[i],
                       xstart = x[i, varnumber] - sd(x[, varnumber])/5,
                       xend = x[i, varnumber] + sd(x[, varnumber])/5,
                       ystart = y[i] - xcoef[i, varnumber] * sd(x[, varnumber])/5,
                       yend = y[i] + xcoef[i, varnumber] * sd(x[, varnumber])/5)
  g <- ggplot(data.frame(x=x[,varnumber], y), aes(x = x, y = y)) + 
    geom_point(color="royalblue", alpha = 0.4)+
    labs(x = colnames(x)[varnumber], y = expression(f(x)) ) +
    theme_minimal() +
    geom_segment(data = slopes, aes(x = xstart, xend = xend, y = ystart, yend = yend), col = "green", size = 2, alpha = 0.4) +
    geom_point(data = slopes, aes(x = x, y = y), color="red", size = 2, alpha = 0.5)
  g
}

#' Plot VarImp effect
#'
#' This function plot VarImp effect of one instance.
#' 
#' @param instancen Instance index.
#' @param x Matrix \eqn{n \times d} with covariates for each sample unit.
#' @param y \eqn{n}-dimensional vector with predictions from the original model.
#' @param xcoef Matrix \eqn{n \times d} with local slopes for each sample unit.
#' @export
#' @return Graph with effect of one instance.
#' @examples
#' library(interpretability)
#' 
#' data("artificial_data")
#' 
#' x_coef <- varimp(x = artificial_data$x_test,
#'                  y = artificial_data$y_pred,
#'                  h = 0.1,
#'                  dmax = ncol(artificial_data$x_test))
#'
#'ploteffect(instancen = 1,
#'           x = artificial_data$x_test,
#'           y = artificial_data$y_pred,
#'           xcoef = x_coef)
 


ploteffect <- function(instancen, x, y, xcoef){
  effect <- xcoef[, -ncol(xcoef)] * x
  df <- data.frame(Feature = factor(names(x[instancen,]), levels = unique(colnames(x)), ordered = TRUE), 
                   Effect = effect[instancen,], 
                   Value = factor(paste(names(x[instancen,]), round(x[instancen, ], digits = 4), sep = "="),
                                  levels = unique(paste(names(x[instancen,]), round(x[instancen, ], digits = 4), sep = "=")), ordered = TRUE) 
  )
  g <- ggplot2::ggplot(data = df, aes(y = Effect, x = Value, fill = Effect > 0)) +
    scale_fill_manual(name = 'Effect > 0', values = setNames(c('green', 'red'), c(T, F))) +
    geom_col(show.legend = FALSE, alpha = 0.6) +
    coord_flip() +
    geom_hline(yintercept=0, color = "gray") +
    labs(y="Effect", x = "Feature", 
         caption = paste("Instance number:", instancen),
         title = paste(" Actual prediction:", round(y[instancen], digits = 4), 
                       "\n Average prediction:", round(mean(y), digits = 4)
         )
    )
  print(g)
}

#' Plot SupClus effect
#'
#' This function plot SupClus effect of one instance.
#' 
#' @param instancen Instance index.
#' @param x Matrix \eqn{n \times d} with covariates for each sample unit.
#' @param y \eqn{n}-dimensional vector with predictions from the original model.
#' @param clustercoef Matrix \eqn{k x d} with local slopes for each cluster.
#' @param cluster \eqn{n}-dimensional vector with cluster of each instance.
#' @export
#' @return Graph with effect of one instance.
#' @examples
#' library(interpretability)
#' 
#' data("artificial_data")
#' 
#' supcluscoef <- supclus(k = 2,
#'                        x = artificial_data$x_test,
#'                        y = artificial_data$y_pred,
#'                        h = 0.1,
#'                        dmax = ncol(artificial_data$x_test))
#'                        
#'ploteffectclus(instancen = 1,
#'               x = artificial_data$x_test,
#'               y = artificial_data$y_pred, 
#'               clustercoef = supcluscoef$coef,  
#'               cluster = supcluscoef$cluster)

ploteffectclus <- function(instancen, x, y, clustercoef, cluster){
  effect <- clustercoef[cluster, ] * x
  df <- data.frame(Feature = factor(names(effect[instancen,]), levels = unique(colnames(x)), ordered = TRUE), 
                   Effect = effect[instancen,], 
                   Value = factor(paste(names(effect[instancen,]), round(x[instancen, ], digits = 4), sep = "="),
                                  levels = unique(paste(names(effect[instancen,]), round(x[instancen, ], digits = 4), sep = "=")), ordered = TRUE) 
  )
  g <- ggplot2::ggplot(data = df, aes(y = Effect, x = Value, fill = Effect > 0)) +
    scale_fill_manual(name = 'Effect > 0', values = setNames(c('green', 'red'), c(T, F))) +
    geom_col(show.legend = FALSE, alpha = 0.6) +
    coord_flip() +
    geom_hline(yintercept=0, color = "gray") +
    labs(y="Effect", x = "Feature", 
         caption = paste("Instance number:", instancen),
         title = paste(" Actual prediction:", round(y[instancen], digits = 4), 
                       "\n Average prediction:", round(mean(y), digits = 4)
         )
    )
  print(g)
}


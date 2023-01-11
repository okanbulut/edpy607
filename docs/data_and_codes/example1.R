
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                      EDPY 607 - MEASUREMENT THEORY II                    ----
##                          FACTOR ANALYTIC METHODS                           ~~
##                                  EXAMPLE 1                                 ~~
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Install the packages
install.packages(c("DataExplorer", "ggcorrplot", "psych", "lavaan", "semPlot", 
                   "ggplot2", "mirt", "MASS"))

# Activate the required packages
library("DataExplorer") # for data summarization
library("ggcorrplot") # for creating correlation plots
library("psych") # for exploratory factor analysis
library("lavaan") # for confirmatory factor analysis
library("semPlot") # for creating path diagrams
library("ggplot2") # for general visualization tasks
library("mirt") # for simulating response data
library("MASS") # for simulating correlated multivariate data

# Read the data in
sapa <- read.csv("sapa_data.csv", header = TRUE)

# View the first 6 rows
head(sapa)

#....................Exploratory Data Analysis...................

DataExplorer::introduce(sapa)

DataExplorer::plot_intro(sapa)

DataExplorer::plot_missing(sapa)

# Save the correlation matrix
cormat <- psych::tetrachoric(x = sapa)$rho

# Print the correlation matrix
print(cormat)

# Correlation matrix plot
ggcorrplot::ggcorrplot(corr = cormat, # correlation matrix
                       type = "lower", # print only the lower part of the correlation matrix
                       hc.order = TRUE, # hierarchical clustering
                       show.diag = TRUE, # show the diagonal values of 1
                       lab = TRUE, # add correlation values as labels
                       lab_size = 3) # Size of the labels


#..................Exploratory Factor Analysis...................

# Try one-factor EFA model --> nfactors = 1
efa.model1 <- psych::fa(r = sapa, nfactors = 1, fm = "ml", cor = "tet")

# Print the results 
print(efa.model1, sort = TRUE) # Show the factor loadings sorted by absolute value

# Try two-factor EFA model --> nfactors=2
efa.model2 <- psych::fa(sapa, nfactors = 2, rotate = "oblimin", fm = "ml", cor = "tet")

# Print the results 
print(efa.model2, sort = TRUE)

# Visualize the model
plot(efa.model2)

# Try four-factor EFA model --> nfactors=4
efa.model4 <- psych::fa(sapa, nfactors = 4, rotate = "oblimin", fm = "ml", cor = "tet")

# Print the results 
print(efa.model4, sort = TRUE)

# Screeplot
data_eigen <- data.frame(ev = efa.model4$values, 
                         # 16 eigenvalues due to having 16 items in SAPA
                         factor_number=factor(1:16)) 

ggplot2::ggplot(data_eigen, 
                aes(x = factor_number, y = ev, group = 1)) + 
  geom_point() + 
  geom_line() +
  theme_bw(base_size=15) +
  xlab("Factor number") + 
  ylab("Eigenvalue")

# Bi-factor EFA model
omega.model <- psych::omega(sapa, nfactors = 4, fm = "ml", poly = TRUE)

# Print the results 
print(omega.model)

#............................Exercises...........................

# 1. Run the four-factor EFA model using rotate = "varimax" for orthogonal rotation and 
# check the model fit (i.e., whether forcing the factors to be uncorrelated improved the 
# model fit).


# 2. The **psych** package help page for the `fa` function has the following information:
  
# A very strong argument against using MLE is found in the chapter by MacCallum, Brown and 
# Cai (2007) who show that OLS approaches produce equivalent solutions most of the time, 
# and better solutions some of the time. This particularly in the case of models with some 
# unmodeled small factors. (See sim.minor to generate such data.) Principal axes may be used 
# in cases when maximum likelihood solutions fail to converge, although fm="minres" will also 
# do that and tends to produce better (smaller RMSEA) solutions.

# Run the four-factor EFA model using fm = "pa" and fm = "minres" and check whether the 
# resulting models produce smaller residuals and RMSEA values. 


  
# 3. Now, let's simulate some ordinal data (e.g., Likert scales) using the Graded Response Model 
# and check which estimator would yield better results. We will use the mirt and MASS packages for 
# data simulation (check out ?psych::sim for built-in simulation functions in the psych package). 
# The following function simulates 12 items based on a two-factor model (6 items per factor). 
# We can determine the sample size using `sample.size` and the correlation between the two factors 
# using `cor`. The final `seed` argument allows to fix the seed in the simulation so that we can 
# create the same dataset (or different datasets) in the future. The first six items are mostly 
# loaded on the first dimension, while the second set of 6 items are mostly loaded on the second 
# dimension. The function returns a list consisting of the factor scores, item parameters, and 
# response data. 

# Let's define a simulation function called simGRMdata
simGRMdata <- function(sample.size, cor, seed) {
  
  require("mirt")
  require("MASS")
  
  # Seed will allow us to generate the same data again later on
  if(!is.null(seed)) {set.seed(seed)}
  
  # Define multidimensional abilities (i.e., factor scores)
  theta <- MASS::mvrnorm(n = sample.size, 
                         # mean for factor scores
                         mu = rep(0, 2),
                         # variance-covariance matrix of matrix scores
                         Sigma = matrix(c(1,cor,cor,1),2,2)) 
  
  # Generate slope (i.e., discrimination) parameters
  a1 <- c(runif(n = 6, min = 0.9, max = 2.4), # First 6 items load heavily on the first factor
          runif(n = 6, min = 0.1, max = 0.4)) # But they load low on the second factor
  
  a2 <- c(runif(n = 6, min = 0.1, max = 0.4), # Second six items load low on the first factor
          runif(n = 6, min = 0.9, max = 2.4)) # But they load heavily on the second factor
  
  a <- as.matrix(cbind(a1, a2), ncol = 2)
  
  # Generate intercept (i.e., difficulty) parameters
  b1 <- runif(n = 12, min = 0.67, max = 2)
  b2 <- b1 - runif(n = 12, min = 0.67, max = 1.34)
  b3 <- b2 - runif(n = 12, min = 0.67, max = 1.34)
  b <- as.matrix(cbind(b1, b2, b3), ncol = 3)
  
  # Generate item responses based on GRM
  resp <- mirt::simdata(a = a, d = b, itemtype = 'graded', Theta = theta)
  
  # Return all the parameters and data
  data <- list(theta = theta,
               parameters = cbind(a, b),
               response = resp)
  
  return(data)
}

# Example dataset
mydata <- simGRMdata(sample.size = 1000, cor = 0.8, seed = 2023)

# View the factor scores
head(mydata$theta)

# View the parameters
head(mydata$parameters)

# View the response data
head(mydata$response)

# Using the function we defined above, generate datasets with different sample sizes and inter-factor 
# correlations and then apply EFA to the response dataset to evaluate the dimensionality. You can use 
# fm = "pa" and fm = "ml" in the estimation. Compare the performance of the two estimation methods.










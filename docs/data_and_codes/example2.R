
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                      EDPY 607 - MEASUREMENT THEORY II                    ----
##                          FACTOR ANALYTIC METHODS                           ~~
##                                  EXAMPLE 2                                 ~~
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
finance <- read.csv("finance.csv", header = TRUE)

# View the first six rows
head(finance)

# The variables in the datasets are as follows:
# * respondent_id: Respondent ID
# * age: Age category
# * gender: Male or Female
# * education: Education level
# * employment: Employment status
# * marital_status: Marital status
# * fwb1_1 to fwb1_6: First six items with 1=Not at all to 5=Completely
# * fwb2_1 to fwb2_4: The last four items with 1=Never to 5=Always
# * raise2000: Confidence in ability to raise $2000 in 30 days
# * financial_knowledge: Overall financial knowledge (self-reported)
# * debt_collector: Contacted by a debt collector in past 12 months
# * hardship_food: Food didn’t last and didn’t have money to get more
# * hardship_doctor: Any household member couldn't afford to see doctor or go to hospital

#....................Exploratory Data Analysis...................

# Describe the dataset
response <- finance %>%
  dplyr::select(dplyr::starts_with("fwb"))

psych::describe(response)

# Recode missing items
response <- apply(response, # data to apply the function
                  2, # 1 to apply to each row; 2 to apply to each column
                  function(x) ifelse(x %in% c(-1, -4), NA, x)) %>%
  as.data.frame()

# Check response utilization
apply(response, 2, table)

# Visualize the data
DataExplorer::introduce(response)

DataExplorer::plot_intro(response)

DataExplorer::plot_histogram(response) 

# Save the correlation matrix
cormat <- psych::polychoric(x = response)$rho

# Print the correlation matrix
print(cormat)

# First, let's save the original variable names
names_fwb <- colnames(response)

# Reverse code negatively-worded items 3,5,6,7,9, and 10
keys <- c(1,1,-1,1,-1,-1,-1,1,-1,-1)
response_final <- psych::reverse.code(keys, 
                                      items = response,
                                      mini = rep(1, 10),
                                      maxi = rep(5, 10))

# View the recoded dataset
head(response_final)

# Now let's change the names to the original names to avoid "-" added to the names
colnames(response_final) <- names_fwb

# Finally, recalculate the correlation matrix
cormat <- psych::polychoric(x = response_final)$rho

# Create a correlation matrix plot with the new data
ggcorrplot::ggcorrplot(corr = cormat, # correlation matrix
                       type = "lower", # print only the lower part of the correlation matrix
                       hc.order = TRUE, # hierarchical clustering
                       show.diag = TRUE, # show the diagonal values of 1
                       lab = TRUE, # add correlation values as labels
                       lab_size = 3) # Size of the labels


#........................Parallel Analysis.......................

parallel.test <- psych::fa.parallel(x = response_final, 
                                    fm="ml", # maximum likelihood
                                    fa="fa", # principal axis factor analysis
                                    cor="poly", # use polychoric correlations
                                    n.iter = 20, # number of iterations
                                    error.bars = TRUE)

print(parallel.test)


#..................Confirmatory Factor Analysis..................

# Define a two-factor model
model.two <- '# Define factors
              pw =~ fwb1_3 + fwb2_4 + fwb2_1 + fwb1_5 + fwb1_6 + fwb2_3
              nw =~ fwb1_2 + fwb1_4 + fwb1_1 + fwb2_2
          
              # Variances and covariances
              pw ~~ pw
              nw ~~ nw
              pw ~~ nw'

# model.two <- '# Define factors
#               pw =~ NA*fwb1_3 + fwb2_4 + fwb2_1 + fwb1_5 + fwb1_6 + fwb2_3
#               nw =~ NA*fwb1_2 + fwb1_4 + fwb1_1 + fwb2_2
#           
#               # Variances and covariances
#               pw ~~ 1*pw
#               nw ~~ 1*nw
#               pw ~~ nw'

# Fit the model - WLSMV
cfa.wlsmv <- lavaan::cfa(model.two, data = response_final, estimator = "WLSMV")

# Model summary
summary(cfa.wlsmv, fit.measures=TRUE, standardized=TRUE)

# Model path diagram
semPlot::semPaths(cfa.wlsmv, "std") # "std" gives standardized loadings

# Fit the model - MLR
cfa.mlr <- lavaan::cfa(model.two, data = response_final, estimator = "MLR")

# Model summary
summary(cfa.mlr, fit.measures=TRUE, standardized=TRUE)

# Bi-factor model
bi.model <- '# Define factors
             g =~ fwb1_3 + fwb2_4 + fwb2_1 + fwb1_5 + fwb1_6 + fwb2_3 + fwb1_2 + fwb1_4 + fwb1_1 + fwb2_2
             pw =~ fwb1_3 + fwb2_4 + fwb2_1 + fwb1_5 + fwb1_6 + fwb2_3
             nw =~ fwb1_2 + fwb1_4 + fwb1_1 + fwb2_2'

bifactor.wlsmv <- cfa(bi.model, data=response_final, estimator = "WLSMV", 
                      orthogonal=TRUE, std.lv=TRUE)
summary(bifactor.wlsmv, fit.measures=TRUE, standardized=TRUE)

semPlot::semPaths(bifactor.wlsmv, "std") # "std" gives standardized loadings

# Model comparison
anova(cfa.wlsmv, bifactor.wlsmv)


#..................Principal Component Analysis..................

# Select the variables of interest
data <- finance %>%
  dplyr::select(raise2000, financial_knowledge,	debt_collector,	
                hardship_food, hardship_doctor)

# Recode missing items
data <- apply(data, # data to apply the function
              2, # 1 to apply to each row; 2 to apply to each column
              function(x) ifelse(x %in% c(-1, 8), NA, x)) %>%
  as.data.frame()

# PCA - 1 Component
pca_results1 <- psych::principal(
  r = data, 
  nfactors = 1, # Number of components to extract
  rotate="varimax",
  cor = "mixed" # for a mixture of tetrachorics, polychorics, and so on
)

print(pca_results1)

# PCA - 2 Components (Uncorrelated)
pca_results2a <- psych::principal(
  r = data, 
  nfactors = 2, 
  rotate="varimax",
  cor = "mixed"
)

print(pca_results2a)

biplot(pca_results2a,
       main = "PCA with Financial Variables: Varimax Rotation")

# PCA - 2 Components (Correlated)
pca_results2b <- psych::principal(
  r = data, 
  nfactors = 2, 
  rotate="promax",
  cor = "mixed"
)

print(pca_results2b)

biplot(pca_results2b,
       main = "PCA with Financial Variables: Promax Rotation")


#............................Exercises...........................

# 1. Fit the two-factor CFA model to the data by forcing the correlation between the factors 
# to be zero and compare the model fit to that of the correlated model we estimated earlier 
# using the `anova` function.  

# 2. Use the following model to fit a second-order (i.e., higher-order) model where the pw and nw 
# factors define a higher-order factor of g (i.e., general factor). Note that there is no perfect 
# way to specify a second-order factor when there are only two first-order factors. However, 
# for the sake of our demonstration, we will give it a try and check its model it. 

so.model <- '# Define factors
             pw =~ fwb1_3 + fwb2_4 + fwb2_1 + fwb1_5 + fwb1_6 + fwb2_3
             nw =~ fwb1_2 + fwb1_4 + fwb1_1 + fwb2_2
             g =~ 1*pw + 1*nw
             g ~~ g'


# Nieto et al. (2021) show that there is an alternative way to model item wording effects in 
# psychological instruments--the random intercept item factor analysis (RIIFA) model. In the 
# RIIFA model introduced by @maydeu2006random, the researcher defines a factor associated with 
# all items and an additional "wording" factor. This additional factor is also associated with 
# all items. Its factor loadings are fixed to 1 but its variance is estimated (see Figure 1, Step 3 
# in Nieto et al. (2021)). Fit the RIIFA model to the FWB scale and evaluate the model fit. 
# For this model, you will use the response dataset instead of response_final. 


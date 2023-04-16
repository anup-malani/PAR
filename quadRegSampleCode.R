# Note that the R code below assumes that the haven library is installed to read Stata files using the read_dta() function. The code first reads the data from the Stata file into a dataframe named data. Then, it generates quadratic terms for the e0 and ed variables. It fits a quadratic regression model using the lm() function, with y as the dependent variable and e0, ed, e02, ed2, and e0d as the independent variables. Finally, the code calculates and tests for the presence of the DC and PAR effects using the t.test() function and the estimated coefficients from the regression model.

# Load data.  If the data are in R format, you can just direcly load the data.
library(haven)
data <- read_dta("filename.dta")

# Clean data
# We assume the data are labeled y, e0, and ed using the definitions in the background section.
# We also assume there are no additional controls. These include group or individual fixed effects.
# I.e., we are assuming the analysis is cross-sectional.

# Generate quadratic terms
data$e02 <- data$e0^2
data$ed2 <- data$ed^2
data$e0d <- data$e0 * data$ed

# Basic quadratic regression
model <- lm(y ~ e0 + ed + e02 + ed2 + e0d, data = data)

# Test for DC and PAR at mean of data

# Gather means of each input variable
m0 <- mean(data$e0)
md <- mean(data$ed)

# Test for DC
# Check if the estimate below is positive and the p-value is below your desired critical value
dc_est <- coef(model)[2] + (2 * coef(model)[4] * m0) + (coef(model)[6] * md)
dc_test <- t.test(dc_est, alternative = "greater")
dc_pval <- dc_test$p.value

# Test for PAR
# Check if the estimate below is negative and the p-value is below your desired critical value
par_est <- coef(model)[3] + (2 * coef(model)[5] * md) + (coef(model)[6] * m0)
par_test <- t.test(par_est, alternative = "less")
par_pval <- par_test$p.value
# Load data.  If the data are in R format, you can just direcly load the data.
library(haven)
data <- read_dta("filename.dta")

# Clean data
# We assume the data are labeled y, e0, and ed using the definitions in the background section.
# We also assume there are no additional controls. These include group or individual fixed effects.
# I.e., we are assuming the analysis is cross-sectional.

# Generate quadratic terms
data$e02 <- data$e0^2
data$ed2 <- data$ed^2
data$e0d <- data$e0 * data$ed

# Basic quadratic regression
model <- lm(y ~ e0 + ed + e02 + ed2 + e0d, data = data)

# Test for DC and PAR at mean of data

# Gather means of each input variable
m0 <- mean(data$e0)
md <- mean(data$ed)

# Test for DC
# Check if the estimate below is positive and the p-value is below your desired critical value
dc_est <- coef(model)[2] + (2 * coef(model)[4] * m0) + (coef(model)[6] * md)
dc_test <- t.test(dc_est, alternative = "greater")
dc_pval <- dc_test$p.value

# Test for PAR
# Check if the estimate below is negative and the p-value is below your desired critical value
par_est <- coef(model)[3] + (2 * coef(model)[5] * md) + (coef(model)[6] * m0)
par_test <- t.test(par_est, alternative = "less")
par_pval <- par_test$p.value
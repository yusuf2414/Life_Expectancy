#### Predicting Life expectancy at Birth from Mortality rates using Lee carter model with Arima 
# -----------------------------------
# Step 1: Load Required Libraries
# -----------------------------------

library("StMoMo")
library("demography")
library("forecast")
library("lifecontingencies")
library(ggplot2)
library(tseries)

# -----------------------------------
# Step 2: Load and Prepare Data
# -----------------------------------

# Fetch mortality data from Human Mortality Database (HMD) ##This example is for Australian Data 
#https://www.mortality.org/Data/DataAvailability  #You will need to create an account 
data <- hmd.mx(country = "AUS", username = "your_email@example.com",
               password = "your_password", label = "Australia")

# View structure of the data
str(data)

# -----------------------------------
# Step 3: Plot Mortality Rates
# -----------------------------------

# Plot mortality rates for males, females, and the total population
plot(data, series = "male", datatype = "rate", main = "Male Rates")
plot(data, series = "female", datatype = "rate", main = "Female Rates")
plot(data, series = "total", datatype = "rate", main = "Total Rates")

# -----------------------------------
# Step 4: Apply Lee-Carter Model
# -----------------------------------

# Lee-Carter model for male, female, and total populations (1991-2020)
# lengthing the years beyond 30 years may cause the 
LcaM <- lca(data, series = "male", max.age = 100, years = 1991:2020)
LcaF <- lca(data, series = "female", max.age = 100, years = 1991:2020)
LcaT <- lca(data, series = "total", max.age = 100, years = 1991:2020)

# Summary of the Lee-Carter Model results
summary(LcaM)
summary(LcaF)
summary(LcaT)

# -----------------------------------
# Step 5: Plot Lee-Carter Model Components
# -----------------------------------

# Plot ax (baseline age-specific mortality rates), bx (age-specific sensitivity to time), and kt (time-varying mortality trends)
par(mfrow = c(1, 3))

# Plot ax for Total, Female, and Male populations
plot(LcaT$ax, main = "ax (Total)", xlab = "Age", ylab = "ax", type = "l")
lines(x = LcaF$age, y = LcaF$ax, col = "red")
lines(x = LcaM$age, y = LcaM$ax, col = "blue")

# Plot bx for Total, Female, and Male populations
plot(LcaT$bx, main = "bx (Total)", xlab = "Age", ylab = "bx", type = "l", col = "black")
lines(x = LcaF$age, y = LcaF$bx, col = "red")
lines(x = LcaM$age, y = LcaM$bx, col = "blue")

# Add legend to bx plot
legend("topleft", legend = c("Male", "Female", "Total"), cex = 0.8, col = c("blue", "red", "black"), lty = 1)

# Plot kt for Total, Female, and Male populations
plot(LcaT$kt, main = "kt (Total)", xlab = "Year", ylab = "kt", type = "l", col = "black")
lines(x = LcaF$year, y = LcaF$kt, col = "red")
lines(x = LcaM$year, y = LcaM$kt, col = "blue")

# Add legend to kt plot
legend("topright", legend = c("Male", "Female", "Total"), cex = 0.8, col = c("blue", "red", "black"), lty = 1)

# -----------------------------------
# Step 6: ADF Test for Stationarity
# -----------------------------------

# Perform ADF test for stationarity on kt for Male, Female, and Total populations
adf_result_T <- adf.test(LcaT$kt)
print("ADF Test for Total Population:")
print(adf_result_T)

adf_result_F <- adf.test(LcaF$kt)
print("ADF Test for Female Population:")
print(adf_result_F)

adf_result_M <- adf.test(LcaM$kt)
print("ADF Test for Male Population:")
print(adf_result_M)

# -----------------------------------
# Step 7: Differencing for Non-Stationary Series
# -----------------------------------

# Difference the series for non-stationary results
if (adf_result_M$p.value > 0.05) {
  LcaM_diff <- diff(LcaM$kt, differences = 1)
  print("Differenced Male Population Series:")
  print(head(LcaM_diff))
} else {
  LcaM_diff <- LcaM$kt  # If already stationary, no differencing
}

if (adf_result_F$p.value > 0.05) {
  LcaF_diff <- diff(LcaF$kt, differences = 1)
  print("Differenced Female Population Series:")
  print(head(LcaF_diff))
} else {
  LcaF_diff <- LcaF$kt  # If already stationary, no differencing
}

if (adf_result_T$p.value > 0.05) {
  LcaT_diff <- diff(LcaT$kt, differences = 1)
  print("Differenced Total Population Series:")
  print(head(LcaT_diff))
} else {
  LcaT_diff <- LcaT$kt  # If already stationary, no differencing
}

# -----------------------------------
# Step 8: Fit ARIMA Model
# -----------------------------------

# Fit ARIMA model for Male, Female, and Total populations
arima_M <- auto.arima(LcaM_diff)
print("ARIMA Model for Male Population:")
print(arima_M)

arima_F <- auto.arima(LcaF_diff)
print("ARIMA Model for Female Population:")
print(arima_F)

arima_T <- auto.arima(LcaT_diff)
print("ARIMA Model for Total Population:")
print(arima_T)

# -----------------------------------
# Step 9: Residuals Analysis
# -----------------------------------

# Extract residuals from the ARIMA models
residuals_M <- residuals(arima_M)
residuals_F <- residuals(arima_F)
residuals_T <- residuals(arima_T)

# Plotting Residuals and ACF for Male Population
par(mfrow = c(2, 2))  # Plot in 2x2 grid
plot(residuals_M, main = "Residuals for Male Population", ylab = "Residuals")
acf(residuals_M, main = "ACF of Residuals for Male Population")

# Plotting Residuals and ACF for Female Population
plot(residuals_F, main = "Residuals for Female Population", ylab = "Residuals")
acf(residuals_F, main = "ACF of Residuals for Female Population")

# Shapiro-Wilk Normality Test for Residuals
shapiro_test_M <- shapiro.test(residuals_M)
shapiro_test_F <- shapiro.test(residuals_F)
shapiro_test_T <- shapiro.test(residuals_T)

# Print Shapiro-Wilk test results
cat("\nShapiro-Wilk Test for Male Population Residuals:\n")
print(shapiro_test_M)

cat("\nShapiro-Wilk Test for Female Population Residuals:\n")
print(shapiro_test_F)

cat("\nShapiro-Wilk Test for Total Population Residuals:\n")
print(shapiro_test_T)

# -----------------------------------
# Step 10: Forecasting with ARIMA Models
# -----------------------------------

# Define forecast years (next 20 years)
forecast_years <- 2021:2041

# Forecast Male Population for the next 20 periods
forecast_M <- forecast(arima_M, h = 20)

# Forecast Female Population for the next 20 periods
forecast_F <- forecast(arima_F, h = 20)

# Forecast Total Population for the next 20 periods
forecast_T <- forecast(arima_T, h = 20)

# -----------------------------------
# Step 11: Reverse Differencing for Forecasts
# -----------------------------------

# Reverse differencing for Male and Female populations
forecast_M_rev <- cumsum(c(LcaM$kt[length(LcaM$kt)], forecast_M$mean))
forecast_F_rev <- cumsum(c(LcaF$kt[length(LcaF$kt)], forecast_F$mean))
forecast_T_rev <- cumsum(c(LcaT$kt[length(LcaT$kt)], forecast_T$mean))

# Remove the first value for correct alignment
forecast_M_rev <- forecast_M_rev[-1]
forecast_F_rev <- forecast_F_rev[-1]
forecast_T_rev <- forecast_T_rev[-1]

# -----------------------------------
# Step 12: Combine Original and Forecasted Series
# -----------------------------------

# For Male Population
forecast_years_M <- seq(max(LcaM$year) + 1, max(LcaM$year) + 20)
combined_years_M <- c(LcaM$year, forecast_years_M)
combined_kt_M <- c(LcaM$kt, forecast_M_rev)

# For Female Population
forecast_years_F <- seq(max(LcaF$year) + 1, max(LcaF$year) + 20)
combined_years_F <- c(LcaF$year, forecast_years_F)
combined_kt_F <- c(LcaF$kt, forecast_F_rev)

# For Total Population
forecast_years_T <- seq(max(LcaT$year) + 1, max(LcaT$year) + 20)
combined_years_T <- c(LcaT$year, forecast_years_T)
combined_kt_T <- c(LcaT$kt, forecast_T_rev)

# -----------------------------------
# Step 13: Plot Combined Time Series
# -----------------------------------

par(mfrow = c(1, 3))

# Plot Male Population (forecasted and original)
plot(combined_years_M, combined_kt_M, type = "l", col = "blue", 
     xlab = "Year", ylab = "Death Rates (kt)", main = "Male Population Forecast")
lines(LcaM$year, LcaM$kt, col = "black")

# Plot Female Population (forecasted and original)
plot(combined_years_F, combined_kt_F, type = "l", col = "red", 
     xlab = "Year", ylab = "Death Rates (kt)", main = "Female Population Forecast")
lines(LcaF$year, LcaF$kt, col = "black")

# Plot Total Population (forecasted and original)
plot(combined_years_T, combined_kt_T, type = "l", col = "black", 
     xlab = "Year", ylab = "Death Rates (kt)", main = "Total Population Forecast")
lines(LcaT$year, LcaT$kt, col = "black")

# Add legends
legend("topright", legend = c("Forecasted Male", "Original Male"), 
       cex = 0.8, col = c("blue", "black"), lty = 1)
legend("topright", legend = c("Forecasted Female", "Original Female"), 
       cex = 0.8, col = c("red", "black"), lty = 1)
legend("topright", legend = c("Forecasted Total", "Original Total"), 
       cex = 0.8, col = c("black", "black"), lty = 1)

# -----------------------------------
# Step 20: Predict on the Test Set
# -----------------------------------

# Define the forecast period (e.g., next 10 years)
forecast_years <- 2021:2031

# Calculate the predicted death rates (mx) for each age group for each forecast year (2022-2032)
mx_predicted_10_years <- matrix(nrow = length(LcaM$age), ncol = length(forecast_years))

for (i in 1:length(forecast_years)) {
  # Calculate predicted mx for each age group at each forecast year
  mx_predicted_10_years[,i] <- exp(LcaM$ax + LcaM$bx * combined_kt_M[i])
}

# Initialize lx for each forecast year (starting with a radix of 100,000)
lx_10_years <- matrix(nrow = length(LcaM$age), ncol = length(forecast_years))
lx_10_years[1,] <- 100000  # Starting with a radix of 100,000

# Calculate lx for each forecast year using the predicted mx values
for (i in 2:length(LcaM$age)) {
  for (j in 1:length(forecast_years)) {
    lx_10_years[i,j] <- lx_10_years[i-1,j] * (1 - mx_predicted_10_years[i-1,j])
  }
}

# Calculate life expectancy (e_x) for each forecast year
e_x_10_years <- matrix(nrow = length(LcaM$age), ncol = length(forecast_years))

for (j in 1:length(forecast_years)) {
  for (i in 1:length(LcaM$age)) {
    e_x_10_years[i,j] <- sum(lx_10_years[i:length(LcaM$age), j]) / lx_10_years[i,j]
  }
}

# Summarize life expectancy for the entire population (at age 0)
life_expectancy_0_10_years <- e_x_10_years[1,]  # Life expectancy at age 0 for each forecast year

# Create a data frame with the life expectancy predictions for forecast years
life_expectancy_df <- data.frame(
  Year = forecast_years,
  LifeExpectancy_0 = life_expectancy_0_10_years
)

# Display and write the forecasted life expectancy to an Excel file
write.xlsx(life_expectancy_df, "forecasted_life_expectancy_new.xlsx", rowNames = FALSE)

# -----------------------------------
# Step 21: Plot Life Expectancy Forecast
# -----------------------------------

# Plot the Life Expectancy at Birth for the next 10 years
ggplot(life_expectancy_df, aes(x = Year, y = LifeExpectancy_0)) +
  geom_line(color = "blue", size = 1) +   # Line plot for life expectancy
  geom_point(color = "red", size = 3) +   # Points to highlight each year
  labs(title = "Forecasted Life Expectancy at Age 0 in 10 years",
       x = "Year",
       y = "Life Expectancy (Years)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center title
  scale_x_continuous(breaks = seq(min(life_expectancy_df$Year), max(life_expectancy_df$Year), by = 1))  # Ensure x-axis shows years correctly

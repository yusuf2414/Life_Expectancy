library("StMoMo")
library("demography")
library("forecast")
library("lifecontingencies")
library(ggplot2)
library(tseries)

data <- hmd.mx(country = "AUS", username = "kiwanuka2414@gmail.com",
               password = "7cXSreMTTkSv?pM", label = "Australia")
#Email that helps use get country code and label 
#https://www.mortality.org/Data/DataAvailability
#You need to visit the HMD data base , get an account and that will assist you to download
#This prediction for the country of australia 
#


str(data)

plot(data,series="male",datatype="rate", main="Male rates")
plot(data,series="female",datatype="rate", main="Female rates") 
plot(data,"total",datatype="rate", main="Total rates")

##### Lee Carter model from Stmomo package , I specify the years to 
#### about 30 years because if you extend it there is alot that could have happened
#### like wars Civil wars or plagues that affect mortality rate trends 
LcaM<-lca(data,series="male",max.age=100,years = 1991:2020) 
LcaF<-lca(data,series="female",max.age=100,years = 1991:2020) 
LcaT<-lca(data,series="total",max.age=100,years = 1991:2020)

#### Summary of the Lee Carter Model
summary(LcaM)
summary(LcaF)
summary(LcaT)

# Plot ax, bx, and kt for Male, Female, and Total populations
#ax is the baseline age-specific mortality rate (ax), capturing mortality at each age.
#bx is the age-specific sensitivity to time trends (bx), showing how each age group's mortality responds to changes over time.
#kt is the time-varying component (kt), representing general mortality trends (e.g., improvements or increases in life expectancy).


par(mfrow=c(1,3)) 

# Plot ax (age-specific mortality rates) for Total, Female, and Male
plot(LcaT$ax, main="ax (Total)", xlab="Age", ylab="ax", type="l") 
lines(x=LcaF$age, y=LcaF$ax, col="red") 
lines(x=LcaM$age, y=LcaM$ax, col="blue")

# Plot bx (age-specific mortality rates) for Total, Female, and Male
plot(LcaT$bx, main="bx (Total)", xlab="Age", ylab="bx", type="l", col="black")
lines(x=LcaF$age, y=LcaF$bx, col="red")
lines(x=LcaM$age, y=LcaM$bx, col="blue")

# Add legend to the bx plot
legend("topleft", legend=c("Male", "Female", "Total"), cex=0.8, col=c("blue", "red", "black"), lty=1)

# Plot kt (time-specific mortality rates) for Total, Female, and Male
plot(LcaT$kt, main="kt (Total)", xlab="Year", ylab="kt", type="l", col="black")
lines(x=LcaF$year, y=LcaF$kt, col="red")
lines(x=LcaM$year, y=LcaM$kt, col="blue")

##### plot shows that mortality rate is reducing with time -- mortality is kt
# Add legend to the kt plot
legend("topright", legend=c("Male", "Female", "Total"), cex=0.8, col=c("blue", "red", "black"), lty=1)


#Perform ADF Test for Stationarity on kt

adf_result_T <- adf.test(LcaT$kt)
print("ADF Test for Total Population:")
print(adf_result_T)

adf_result_F <- adf.test(LcaF$kt)
print("ADF Test for Female Population:")
print(adf_result_F)

adf_result_M <- adf.test(LcaM$kt)
print("ADF Test for Male Population:")
print(adf_result_M)

#Differencing if necessary (for non-stationary series)
### I am differencing for Male and Female Populations 
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
  print("Difference total Population Series:")
  print(head(LcaT_diff))
} else {
  LcaT_diff <- LcaT$kt  # If already stationary, no differencing
}

#Fit ARIMA model for Male and Female populations
# ARIMA model for Male Population (using differenced data if necessary)
arima_M <- auto.arima(LcaM_diff)
print("ARIMA Model for Male Population:")
print(arima_M)

# ARIMA model for Female Population (using differenced data if necessary)
arima_F <- auto.arima(LcaF_diff)
print("ARIMA Model for Female Population:")
print(arima_F)

# ARIMA model for Total Population (using differenced data if necessary)
arima_T <- auto.arima(LcaT_diff)
print("ARIMA Model for Female Population:")
print(arima_T)

# Extracting Residuals for Male and Female ARIMA Models
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


#Forecasting with ARIMA models

# Forecast Male Population for the next 50 periods
forecast_M <- forecast(arima_M, h=20)
# Forecast Female Population for the next 50 periods
forecast_F <- forecast(arima_F, h=20)

# Forecast Total Population for the next 50 periods (adjust if needed)
forecast_T <- forecast(arima_T, h=20)

#Reverse differencing for Male and Female populations to get forecasted values on original scale
forecast_M_rev <- cumsum(c(LcaM$kt[length(LcaM$kt)], forecast_M$mean))  # Reversing differencing for Male
forecast_F_rev <- cumsum(c(LcaF$kt[length(LcaF$kt)], forecast_F$mean))  # Reversing differencing for Female
forecast_T_rev <- cumsum(c(LcaT$kt[length(LcaT$kt)], forecast_T$mean))

# Remove the first value for correct alignment
forecast_M_rev <- forecast_M_rev[-1]
forecast_F_rev <- forecast_F_rev[-1]
forecast_T_rev <- forecast_T_rev[-1]


#Combine original and forecasted series for Male, Female, and Total populations

# For Male Population
forecast_years_M <- seq(max(LcaM$year) + 1, max(LcaM$year) + 20)  # Extend years for forecasted data
combined_years_M <- c(LcaM$year, forecast_years_M)  # Combine original years and forecast years
combined_kt_M <- c(LcaM$kt, forecast_M_rev)  # Combine original kt and forecasted kt

# For Female Population
forecast_years_F <- seq(max(LcaF$year) + 1, max(LcaF$year) + 20)
combined_years_F <- c(LcaF$year, forecast_years_F)
combined_kt_F <- c(LcaF$kt, forecast_F_rev)

# For Total Population
forecast_years_T <- seq(max(LcaT$year) + 1, max(LcaT$year) + 20)
combined_years_T <- c(LcaT$year, forecast_years_T)
combined_kt_T <- c(LcaT$kt, forecast_T_rev)


# Step 8: Plot the combined time series for Male, Female, and Total populations

par(mfrow = c(1, 3))  # 3 plots in 1 row

# Plot Male Population (forecasted and original)
plot(combined_years_M, combined_kt_M, type = "l", col = "blue", 
     xlab = "Year", ylab = "Death Rates (kt)", main = "Male Population Forecast")
lines(LcaM$year, LcaM$kt, col = "black")  # Original data in black

# Plot Female Population (forecasted and original)
plot(combined_years_F, combined_kt_F, type = "l", col = "red", 
     xlab = "Year", ylab = "Death Rates (kt)", main = "Female Population Forecast")
lines(LcaF$year, LcaF$kt, col = "black")  # Original data in black

# Plot Total Population (forecasted and original)
plot(combined_years_T, combined_kt_T, type = "l", col = "black", 
     xlab = "Year", ylab = "Death Rates (kt)", main = "Total Population Forecast")
lines(LcaT$year, LcaT$kt, col = "black")  # Original data in black

# Step 9: Add legends to the plots
legend("topright", legend = c("Forecasted Male", "Original Male"), 
       cex = 0.8, col = c("blue", "black"), lty = 1)
legend("topright", legend = c("Forecasted Female", "Original Female"), 
       cex = 0.8, col = c("red", "black"), lty = 1)
legend("topright", legend = c("Forecasted Total", "Original Total"), 
       cex = 0.8, col = c("black", "black"), lty = 1)

# Forecast years
forecast_years <- 2021:2031

# 1. Calculate the predicted mx for each age group and for each forecast year (2022-2032)
mx_predicted_10_years <- matrix(nrow = length(LcaM$age), ncol = length(forecast_years))

for (i in 1:length(forecast_years)) {
  # Calculate mx for each age group at each forecast year
  mx_predicted_10_years[,i] <- exp(LcaM$ax + LcaM$bx * combined_kt_M[i])
}

# 2. Initialize lx for each forecast year (starting with a radix of 100,000)
lx_10_years <- matrix(nrow = length(LcaM$age), ncol = length(forecast_years))
lx_10_years[1,] <- 100000  # Starting with a radix of 100,000

# 3. Calculate lx for each forecast year (using the predicted mx)
for (i in 2:length(LcaM$age)) {
  for (j in 1:length(forecast_years)) {
    lx_10_years[i,j] <- lx_10_years[i-1,j] * (1 - mx_predicted_10_years[i-1,j])
  }
}

# 4. Calculate life expectancy (e_x) for each forecast year
e_x_10_years <- matrix(nrow = length(LcaM$age), ncol = length(forecast_years))

for (j in 1:length(forecast_years)) {
  for (i in 1:length(LcaM$age)) {
    e_x_10_years[i,j] <- sum(lx_10_years[i:length(LcaM$age), j]) / lx_10_years[i,j]
  }
}

# 5. Summarize life expectancy for the entire population (at age 0)
life_expectancy_0_10_years <- e_x_10_years[1,]  # Life expectancy at age 0 for each forecast year

# 6. Create a data frame with the life expectancy predictions
life_expectancy_df <- data.frame(
  Year = forecast_years,
  LifeExpectancy_0 = life_expectancy_0_10_years
)

# Display the life expectancy predictions for the next 10 years

#Write the life expectancy forecast to an Excel file
write.xlsx(life_expectancy_df, "forecasted_life_expectancy_new.xlsx", rowNames = FALSE)

# Life Expectancy at Birth Plot
ggplot(life_expectancy_df, aes(x = Year, y = LifeExpectancy_0)) +
  geom_line(color = "blue", size = 1) +   # Line plot for life expectancy
  geom_point(color = "red", size = 3) +   # Points to highlight each year
  labs(title = "Forecasted Life Expectancy at Age 0 in 10 years",
       x = "Year",
       y = "Life Expectancy (Years)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center title
  scale_x_continuous(breaks = seq(min(life_expectancy_df$Year), max(life_expectancy_df$Year), by = 1))  # Ensure x-axis shows years correctly



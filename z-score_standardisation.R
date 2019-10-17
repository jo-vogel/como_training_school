# Z-score standardisation
# Author: Johannes Vogel

load('./Processed data/Seasonal_monthly_variables.RData') # Matrix with all pretictands
Season_month_variables_stand <- apply(Data,2,scale) # z-score standardisation

# Check out the z-score standardisation visually
mean_original <- apply(Data,2,mean)
sd_original <- apply(Data,2,sd)

mean_stand <- apply(Season_month_variables_stand,2,mean)
sd_stand <- apply(Season_month_variables_stand,2,sd)

plot(mean_original)
plot(mean_stand)
plot(sd_original)
plot(sd_stand)
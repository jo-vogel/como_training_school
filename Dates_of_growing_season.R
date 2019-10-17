# Create dates of growing season
Dates_of_growing_season_CN <- seq(as.Date("2035-09-29"), as.Date("2036-07-18"), by="days")
Dates_of_growing_season_FR <- seq(as.Date("2035-10-31"), as.Date("2036-09-04"), by="days")
Dates_of_growing_season_US <- seq(as.Date("2035-10-03"), as.Date("2036-08-12"), by="days")

# Save data
save(Dates_of_growing_season_CN,Dates_of_growing_season_FR,Dates_of_growing_season_US, file="Dates_of_growing_season.RData")
write.csv(Dates_of_growing_season_CN,file='Dates_of_growing_season_CN.csv')
write.csv(Dates_of_growing_season_CN,file='Dates_of_growing_season_FR.csv')
write.csv(Dates_of_growing_season_CN,file='Dates_of_growing_season_US.csv')

# Load data
load('Dates_of_growing_season.RData')
Dates_of_growing_season_CN <- read.csv('Dates_of_growing_season_CN.csv',header=T)
Dates_of_growing_season_FR <- read.csv('Dates_of_growing_season_FR.csv',header=T)
Dates_of_growing_season_US <- read.csv('Dates_of_growing_season_US.csv',header=T)

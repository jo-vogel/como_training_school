#create index matrix for months and seasons
# file by Elisabeth and Leila

load('cube_FR.RData')
load('Dates_of_growing_season.RData')

growing_season<-cube_FR
Dates_of_growing_season_1year<- read.csv("Dates_of_growing_season_1year.csv")

# months <- as.numeric(format(as.Date(time_gs_1year),
                            # "%m"))
months <- as.numeric(format(as.Date(Dates_of_growing_season_1year[[1]]),
                            "%m"))
index_matrix <- matrix(NA, ncol = 310, nrow=16)
for(m in 1:12){
  index_matrix[m,] <- months==m
}
index_matrix[13,] <- months %in% c(12,1,2)
index_matrix[14,] <- months %in% c(3,4,5)
index_matrix[15,] <- months %in% c(6,7,8)
index_matrix[16,] <- months %in% c(9,10,11)

save(index_matrix, file="index.RData")

#monthly means (starting with november)
monthly_seasonal_mean <- array(NA, c(1600,16,6)) #the last four columns are seasonal means
for(n in 1:6){
  if(n==2){
    for(m in 1:12){
      if(m!=10){
        M<-(m+2)%%12
        monthly_seasonal_mean[,M,n] <-
          apply(growing_season[,index_matrix[m,],n], 1, sum)
      }
    }
    
  }else{
    for(m in 1:12){
      if(m!=10){
        M<-(m+2)%%12
        monthly_seasonal_mean[,M,n] <-
          apply(growing_season[,index_matrix[m,],n], 1, mean)
      }
    }
  }
}
#seasonal means (starting with winter)
for(n in 1:6){
  if(n==2){
    for(m in 13:16){
      monthly_seasonal_mean[,m,n] <-
        apply(growing_season[,index_matrix[m,],n], 1, sum)
    }
  }else{
    for(m in 13:16){
      monthly_seasonal_mean[,m,n] <-
        apply(growing_season[,index_matrix[m,],n], 1, mean)
    }
  }
}


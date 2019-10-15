# Binomial logist regression model
setwd("C:/Users/vogel/boxup/Promotion/Damocles COST Action/Training school/Group1_Project")
load('./Processed data/cube_FR.RData')
load('./Processed data/Dates_of_growing_season.RData')
# time_gs_1year2 <- read.csv("Dates_of_growing_season_1year.csv") # does not work, object has wrong class
library(ROCR)
library(car)


# Calculate dew point depression (Tmax-Dewpoint)
dpd <- cube_FR[,,5]-cube_FR[,,1]

# select pretictands ####


# First approach: start with monthly temperature maximum / minimum, precipitation sum

# Nov-Aug: 10 months; 3 variables; 10*3=30 pretictands
months <- unlist(strsplit(as.character(time_gs_1year),split="-"))
vec <- seq(2,929,3)
months <- months[vec]

# Aggregate months ####
months2 <- c("01","02","03","04","05","06","07","08","11","12")

tmax_list <- vector("list", 10)
j <- 1
for (i in c(01,02,03,04,05,06,07,08,11,12)){
  # nov_tmax <- cube_FR[,which(months==11),5]
  tmax_list[[j]] <- cube_FR[,which(months==months2[j]),5]
  j <- j+1
}

tmin_list <- vector("list", 10)
j <- 1
for (i in c(01,02,03,04,05,06,07,08,11,12)){
  # nov_tmax <- cube_FR[,which(months==11),5]
  tmin_list[[j]] <- cube_FR[,which(months==months2[j]),6]
  j <- j+1
}

prec_list <- vector("list", 10)
j <- 1
for (i in c(01,02,03,04,05,06,07,08,11,12)){
  # nov_tmax <- cube_FR[,which(months==11),5]
  prec_list[[j]] <- cube_FR[,which(months==months2[j]),2]
  j <- j+1
}

dps_list <- vector("list", 10)
j <- 1
for (i in c(01,02,03,04,05,06,07,08,11,12)){
  # nov_tmax <- cube_FR[,which(months==11),5]
  dps_list[[j]] <- cube_FR[,which(months==months2[j]),1]
  j <- j+1
}

dpd_list <- vector("list", 10)
j <- 1
for (i in c(01,02,03,04,05,06,07,08,11,12)){
  # nov_tmax <- cube_FR[,which(months==11),5]
  dpd_list[[j]] <- dpd[,which(months==months2[j])]
  j <- j+1
}

message('attribute month names')


# Monthly means and sums ####
tmax_mean <- vector("list", 10)
for (i in 1:10){
  tmax_mean [[i]] <- apply(tmax_list[[i]],1,mean)
}

tmin_mean <- vector("list", 10)
for (i in 1:10){
  tmin_mean [[i]] <- apply(tmin_list[[i]],1,mean)
}

prec_sums <- vector("list", 10)
for (i in 1:10){
  prec_sums [[i]] <- apply(prec_list[[i]],1,sum)
}

dps_mean <- vector("list", 10)
for (i in 1:10){
  dps_mean [[i]] <- apply(dps_list[[i]],1,mean)
}

dpd_mean <- vector("list", 10)
for (i in 1:10){
  dpd_mean [[i]] <- apply(dpd_list[[i]],1,mean)
}

# Annual means and sums for  max, min temp, prec. sum ####
tmax_annual_mean <- apply(cube_FR[,,5],1,mean)
tmin_annual_mean <- apply(cube_FR[,,6],1,mean)
prec_annual_sum <- apply(cube_FR[,,2],1,sum)
prec_spr_sum_sum <- apply(cube_FR[,131:310,2],1,sum) # spring & summer
dps_annual_mean <- apply(cube_FR[,,2],1,mean)
dpd_annual_mean <- apply(dpd,1,mean)

# transform yield data to binary data: severe loss, no severe loss
low_yield <- quantile(cy_gsl_FR[,1],0.05)
# ind <- cy_gsl_FR[which(cy_gsl_FR[,1]<low_yield),1]
ind <- which(cy_gsl_FR[,1]<low_yield)
cy <- cy_gsl_FR[,1]
cy[ind] <- 0
cy[-ind] <- 1
low_yield_FR <- cy_gsl_FR[ind,1]




# build models #####

annual_tmax_min_pr_cy <- cbind(tmax_annual_mean,tmin_annual_mean,prec_annual_sum,cy)
annual_tmax_min_pr_cy_df <- as.data.frame(annual_tmax_min_pr_cy) 

training <- annual_tmax_min_pr_cy_df [1:1000,1:4]
testing <- annual_tmax_min_pr_cy_df [1001:1600,1:4]


model1<- glm(cy ~.,family=binomial(link='logit'),data=training)
model2<- glm(cy ~ tmax_annual_mean*prec_annual_sum,family=binomial(link='logit'),data=training)
# model <- glm(annual_tmax_min_pr_cy_df[,4] ~.,family=binomial(link='logit'),data=train)
summary(model1)
# model$data
# model$fitted.values
# model$model

ano_mod <- anova(model1, model2)
summary(ano_mod)
step(model1,model2)
step(model2)

# Test data set
fitted.results <- predict(model1,testing[,1:3],type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError <- mean(fitted.results != testing$cy)
print(paste('Accuracy',1-misClasificError))


# contingency table
obs_pred <- cbind(testing[,4],fitted.results)
tt <- sum(rowSums(obs_pred)==2)
ff <- sum(rowSums(obs_pred)==0)
tf <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
ft <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)

tt+ff+tf+ft

con_tab <- matrix(c(tt,ft,tf,ff),nrow=2,ncol=2)
con_tab2 <- con_tab
colnames(con_tab) <- c('TRUE','FALSE')
rownames(con_tab) <- c('TRUE','FALSE')
con_tab[,1] <- con_tab[,1]/540
con_tab[,2] <-con_tab[,2]/60


# Most relevant variables:
# dew point depression
# (Jun),July/Aug
# Spring, Summer prec.



# Monthly predictors ####
# starting set: july / aug dew point depression


monthly_dpd_jul_aug_cy <- cbind(dpd_mean[[7]],dpd_mean[[8]],cy)
monthly_dpd_jul_aug_cy_df <- as.data.frame(monthly_dpd_jul_aug_cy) 


training <- monthly_dpd_jul_aug_cy_df [1:1000,1:3]
testing <- monthly_dpd_jul_aug_cy_df [1001:1600,1:3]


model1<- glm(cy ~.,family=binomial(link='logit'),data=training)
summary(model1)


# Combination 1 - 3 ###

# Comb. 1
# tmax in aug, prec. spring/summer, dpd in aug
monthly_tmax_prec_dpd_aug_sp.su_aug_cy <- cbind(tmax_mean[[8]],dpd_mean[[8]],prec_spr_sum_sum,cy)
monthly_tmax_prec_dpd_aug_sp.su_aug_cy_df <- as.data.frame(monthly_tmax_prec_dpd_aug_sp.su_aug_cy) 
colnames(monthly_tmax_prec_dpd_aug_sp.su_aug_cy_df) <- c('T_Max_Aug','Dpd_Aug','P_Spr_Sum','Crop_yield')

training <- monthly_tmax_prec_dpd_aug_sp.su_aug_cy_df [1:1000,1:4]
testing <- monthly_tmax_prec_dpd_aug_sp.su_aug_cy_df [1001:1600,1:4]

model1<- glm(Crop_yield ~.,family=binomial(link='logit'),data=training)
summary(model1)

fitted.results <- predict(model1,testing[,1:3],type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

# contingency table
obs_pred <- cbind(testing[,4],fitted.results)
tt <- sum(rowSums(obs_pred)==2)
ff <- sum(rowSums(obs_pred)==0)
tf <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
ft <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
con_tab1 <- matrix(c(tt,ft,tf,ff),nrow=2,ncol=2)
con_tab1b <- con_tab1
colnames(con_tab1) <- c('TRUE','FALSE');rownames(con_tab1) <- c('TRUE','FALSE')
con_tab1[,1] <- con_tab1[,1]/540
con_tab1[,2] <-con_tab1[,2]/60

# Comb. 2
# monthly_tmax_prec_dpd_juag_sp.su_juag_cy <- cbind(tmax_mean[[7]]+tmax_mean[[8]],dpd_mean[[7]]+dpd_mean[[8]],prec_spr_sum_sum,cy)
monthly_tmax_prec_dpd_juag_sp.su_juag_cy <- cbind(tmax_mean[[7]],tmax_mean[[8]],dpd_mean[[7]],dpd_mean[[8]],prec_spr_sum_sum,cy)
monthly_tmax_prec_dpd_juag_sp.su_juag_cy_df <- as.data.frame(monthly_tmax_prec_dpd_juag_sp.su_juag_cy) 
colnames(monthly_tmax_prec_dpd_juag_sp.su_juag_cy_df) <- c('T_Max_Jul','T_Max_Aug','Dpd_Jul','Dpd_Aug','P_Spr_Sum','Crop_yield')


# training <- monthly_tmax_prec_dpd_juag_sp.su_juag_cy_df [1:1000,1:4]
# testing <- monthly_tmax_prec_dpd_juag_sp.su_juag_cy [1001:1600,1:4]
training <- monthly_tmax_prec_dpd_juag_sp.su_juag_cy_df [1:1000,1:6]
# testing <- as.data.frame(monthly_tmax_prec_dpd_juag_sp.su_juag_cy [1001:1600,1:6])
testing <- monthly_tmax_prec_dpd_juag_sp.su_juag_cy_df [1001:1600,1:6]

model2<- glm(Crop_yield ~.,family=binomial(link='logit'),data=training)
summary(model2)

fitted.results <- predict(model2,testing[,1:5],type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

# contingency table
obs_pred <- cbind(testing[,6],fitted.results)
tt <- sum(rowSums(obs_pred)==2)
ff <- sum(rowSums(obs_pred)==0)
tf <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
ft <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
con_tab2 <- matrix(c(tt,ft,tf,ff),nrow=2,ncol=2)
con_tab2b <- con_tab2
colnames(con_tab2) <- c('TRUE','FALSE');rownames(con_tab2) <- c('TRUE','FALSE')
con_tab2[,1] <- con_tab2[,1]/540
con_tab2[,2] <-con_tab2[,2]/60
  
# Comb. 3
# monthly_tmax_prec_dpd_JJA_sp.su_JJA_cy <- cbind(tmax_mean[[6]]+tmax_mean[[7]]+tmax_mean[[8]],dpd_mean[[6]]+dpd_mean[[7]]+dpd_mean[[8]],prec_spr_sum_sum,cy)
monthly_tmax_prec_dpd_JJA_sp.su_JJA_cy <- cbind(tmax_mean[[6]],tmax_mean[[7]],tmax_mean[[8]],dpd_mean[[6]],dpd_mean[[7]],dpd_mean[[8]],prec_spr_sum_sum,cy)
monthly_tmax_prec_dpd_JJA_sp.su_JJA_cy_df <- as.data.frame(monthly_tmax_prec_dpd_JJA_sp.su_JJA_cy) 
colnames(monthly_tmax_prec_dpd_JJA_sp.su_JJA_cy_df) <- c('T_Max_Jun','T_Max_Jul','T_Max_Aug','Dpd_Jun','Dpd_Jul','Dpd_Aug','P_Spr_Sum','Crop_yield')

training <- monthly_tmax_prec_dpd_JJA_sp.su_JJA_cy_df [1:1000,1:8]
testing <- monthly_tmax_prec_dpd_JJA_sp.su_JJA_cy_df [1001:1600,1:8]

model3<- glm(Crop_yield ~.,family=binomial(link='logit'),data=training)
summary(model3)

fitted.results <- predict(model3,testing[,1:7],type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

# contingency table
obs_pred <- cbind(testing[,8],fitted.results)
tt <- sum(rowSums(obs_pred)==2)
ff <- sum(rowSums(obs_pred)==0)
tf <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
ft <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
con_tab3 <- matrix(c(tt,ft,tf,ff),nrow=2,ncol=2)
con_tab3b <- con_tab3
colnames(con_tab3) <- c('TRUE','FALSE');rownames(con_tab3) <- c('TRUE','FALSE')
con_tab3[,1] <- con_tab3[,1]/540
con_tab3[,2] <-con_tab3[,2]/60

# Comb. 4: Interactions
monthly_tmax_prec_dpd_aug_sp.su_aug_cy_int <- cbind(tmax_mean[[8]],dpd_mean[[8]],prec_spr_sum_sum,cy)
monthly_tmax_prec_dpd_aug_sp.su_aug_cy_df_int <- as.data.frame(monthly_tmax_prec_dpd_aug_sp.su_aug_cy_int ) 
colnames(monthly_tmax_prec_dpd_aug_sp.su_aug_cy_df_int) <- c('T_Max_Aug','Dpd_Aug','P_Spr_Sum','Crop_yield')

training <- monthly_tmax_prec_dpd_aug_sp.su_aug_cy_df_int  [1:1000,1:4]
testing <- monthly_tmax_prec_dpd_aug_sp.su_aug_cy_df_int  [1001:1600,1:4]

model4<- glm(Crop_yield ~T_Max_Aug*Dpd_Aug*P_Spr_Sum,family=binomial(link='logit'),data=training)
summary(model4)

fitted.results <- predict(model4,testing[,1:3],type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

# contingency table
obs_pred <- cbind(testing[,4],fitted.results)
tt <- sum(rowSums(obs_pred)==2)
ff <- sum(rowSums(obs_pred)==0)
tf <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
ft <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
con_tab4 <- matrix(c(tt,ft,tf,ff),nrow=2,ncol=2)
con_tab4b <- con_tab4
colnames(con_tab4) <- c('TRUE','FALSE');rownames(con_tab4) <- c('TRUE','FALSE')
con_tab4[,1] <- con_tab4[,1]/540
con_tab4[,2] <-con_tab4[,2]/60

# ROC ####

p <- predict(model1, newdata=subset(testing,select=c(1,2,3)), type="response")
pr <- prediction(p, testing$Crop_yield)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc@y.values[[1]]

# jpeg('Combination1.jpg')
# plot(prf)
# dev.off()




# Evaluate model performance: select best model ####
summary(anova(model3,model2))
summary(anova(model2,model1))
summary(anova(model4,model1))

step(model1)
# step(model2) # does not work for unknown reason
# step(model3) # does not work for unknown reason
step(model4)

confint(model1)
exp(coef(model1))


vif(model1) # multicollinearity (https://rpubs.com/ranvirkumarsah/LR)
1/vif(model1)
mean(vif(model1))

# stepwise regression (forward / backward selection)

# - 1) Select models: account for interactions in the glm
# - 2) Run models
# - 3) contigency matrix
# - 4) ROC
# - 5) Compare models

save.image(file='./Code/Workspaces/Binomial logit regression model.RData')

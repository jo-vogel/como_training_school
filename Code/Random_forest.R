# Random forest appproach

library(randomForest)
library(pbapply)
library(tidyverse)
library(viridis)

source('./Code/Lasso_interact_global_preparation_incl_ext_ind.R') # monthly data including extreme indices


# Create a Random Forest model with default parameters

# model_rf_reg <- randomForest(x=x1_train_list[[1]],y=y1_train_list[[1]],importance = TRUE)
# model_rf_reg
model_rf <- randomForest(x=x1_train_list[[1]],y=as.factor(y1_train_list[[1]]),importance = TRUE) # mtry = 5 in this case
model_rf

# Predicting on train set
pred_set <- matrix(unlist(x1_train_list[[1]]), nrow=dim(x1_train_list[[1]])[1], byrow=F)
# mytrainset <- cbind(mymat,unlist(y1_train_list[[1]]))
predTrain <- predict(model_rf, pred_set, type = "class")
# Checking classification accuracy
mean(predTrain == y1_train_list[[1]])
table(predTrain, y1_train_list[[1]])


# model_rf$predicted and predTrain are different, reason not known so far

# Predicting on Validation set
validset <- matrix(unlist(x1_test_list[[1]]), nrow=dim(x1_test_list[[1]])[1], byrow=F)
predValid <- predict(model_rf, validset, type = "class")
# Checking classification accuracy
mean(predValid == y1_test_list[[1]])                    
table(predValid,y1_test_list[[1]])



# To check important variables
importance(model_rf)        
varImpPlot(model_rf) 



# Using For loop to identify the right mtry for model
a=c()
i=5
for (i in 3:8) {
  model_rf3 <- randomForest(x=x1_train_list[[1]],y=as.factor(y1_train_list[[1]]), ntree = 500, mtry = i,importance = TRUE)
  predValid <- predict(model_rf3, validset, type = "class")
  a[i-2] = mean(predValid == y1_test_list[[1]])
}

a

plot(3:8,a)





# all pixels ####
#################

model_rf <- pblapply(1:length(x1_train_list), function(x) randomForest(x=x1_train_list[[x]],y=as.factor(y1_train_list[[x]]),importance = TRUE)) # mtry = 5 in this case

# Predicting on train set
pred_set <- pblapply (1:length(x1_train_list),function(x) matrix(unlist(x1_train_list[[x]]), nrow=dim(x1_train_list[[x]])[1], byrow=F))
# mytrainset <- cbind(mymat,unlist(y1_train_list[[1]]))
# predTrain <- predict(model_rf, mypred, type = "class")
predTrain <- pblapply (1:length(x1_train_list),function(x) predict(model_rf[[x]], pred_set[[x]], type = "class"))
# Checking classification accuracy
mean_acc_train <- sapply (1:length(x1_train_list),function(x) mean(predTrain[[x]] == y1_train_list[[x]]))
plot(mean_acc_train)
table(predTrain[[100]], y1_train_list[[100]])


# model_rf$predicted and predTrain are different, reason not known so far, see contingency_table_problem.r

# Predicting on Validation set
valid_set <- pblapply (1:length(x1_train_list),function(x) matrix(unlist(x1_test_list[[x]]), nrow=dim(x1_test_list[[x]])[1], byrow=F))
# predValid <- predict(model_rf, myvalid, type = "class")
predValid <- pblapply (1:length(x1_train_list),function(x) predict(model_rf[[x]], valid_set[[x]], type = "class"))
# Checking classification accuracy
mean_acc_test <- sapply (1:length(x1_train_list),function(x) mean(predValid[[x]] == y1_test_list[[x]]))   
plot(mean_acc_test)
table(predValid[[1]],y1_test_list[[1]])


# To check important variables for examplary pixel
importance(model_rf[[100]])      
varImpPlot(model_rf[[100]]) 


# Using For loop to identify the right mtry for model
mean_acc_test=c()
speci <- c()
model_mtry_all <- vector("list",length=5)
predValid_all <- vector("list",length=5)

# for (i in 3:8) {
# for (i in 9:14) {  
# for (i in 15:20) {  
# for (i in 21:25) {
for (i in 17) { # you get 2 warnings, because there are 2 pixel with only 16 predictors, which is lower than the number of trees mtrees
# for (i in 26:30) {
# for (i in 31:35) { # speci does not increase anymore
  k <- i - 16
  model_mtry <-  pblapply(1:length(x1_train_list), function(x) randomForest(x=x1_train_list[[x]],y=as.factor(y1_train_list[[x]]), mtry=i, importance = TRUE))
  predValid <- pblapply (1:length(x1_train_list),function(x) predict(model_mtry[[x]], valid_set[[x]], type = "class"))
  model_mtry_all[[k]] <- model_mtry
  predValid_all[[k]] <- predValid
  mean_acc_test[k] <- mean(sapply (1:length(x1_train_list),function(x) mean(predValid_all[[k]][[x]] == y1_test_list[[x]])))
  plot(mean_acc_test)
  speci[k] <- mean(sapply(1:length(x1_train_list), function(x){InformationValue::specificity(y1_test_list[[x]],(as.numeric(predValid_all[[k]][[x]])-1))}))
}
save.image(file="D:/user/vogelj/Group_project/Code/Workspaces/rf_21_25.RData")

mean_acc_test
speci
plot(mean_acc_test)
plot(speci)


# Further performance assessments
work_pix <- 1:965
obs_pred <- lapply(seq_along(work_pix), function(x){cbind(y1_test_list[[x]],as.numeric(predValid_all[[6]][[x]])-1)})
# the conversion of predValid from factor to numeric sets (0,1) to (1,2), therefore this correction by subtracting 1
tp <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==2)}) # Correct rejections
tn <- sapply(seq_along(work_pix), function(x){sum(rowSums(obs_pred[[x]])==0)}) # Hits
fp <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[x]][,1]==0 & obs_pred[[x]][,2]==1)}) # Misses
fn <- sapply(seq_along(work_pix), function(x){sum(obs_pred[[x]][,1]==1 & obs_pred[[x]][,2]==0)}) # False alarm


# Calculate specificity
speci <- rep(NA,965)
speci[work_pix] <- sapply(seq_along(work_pix), function(x){InformationValue::specificity(y1_test_list[[x]],(as.numeric(predValid[[x]])-1))})
mean(speci)
plot(speci)

# Calculate sensitivity
sensi <- rep(NA,965)
sensi[work_pix] <- sapply(seq_along(work_pix), function(x) {InformationValue::sensitivity(y1_test_list[[x]],(as.numeric(predValid[[x]])-1))})
mean(sensi)

# Calculate CSI 
csi_rf <- rep(NA,965)
csi_rf[work_pix] <- sapply(seq_along(work_pix), function(x){tn[x]/(tn[x]+fn[x]+fp[x])})
mean(csi_rf)
plot(csi_rf)



# Spatial plots
model_name <- "Random forest"
world <- map_data("world")
DF_speci <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], specificity = speci)

ggplot(data = DF_speci, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=speci)) +
  scale_fill_viridis(na.value="grey50")+
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Specif.",
       title = paste("Specificity, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/Specificity_rand_forest_map.png")


DF_csi <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], Critical_success_index = csi_rf)

ggplot(data = DF_csi, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=csi_rf)) +
  scale_fill_viridis(na.value="grey50")+
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="CSI",
       title = paste("Critical succes index, simple",model_name,"regression"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))  +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/CSI_rand_forest_map.png")



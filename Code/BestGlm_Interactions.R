##### CODE BESTGLM FR #####
#PRESELECTION: Linear Correlation #
#Author: Cristina Deidda #


setwd("C:/Users/39349/Documenti/DAMOCLES")
load('cube_FR.RData')
# time_gs_1year2 <- read.csv("Dates_of_growing_season_1year.csv") # does not work, wrong class
load('Dates_of_growing_season.RData')
#library(ROCR)

library(lattice)
library(bestglm)
library(ROCR) 
library(tictoc)
library(car)
#create index matrix for months and seasons

growing_season<-cube_FR
Dates_of_growing_season_1year<- read.csv("Dates_of_growing_season_1year.csv")


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

###### MATRIX with month mean and seasonal mean ###########
### Matrix MOnth starts from November to October and concatenate the 4 season##### 

monthly_season<-monthly_seasonal_mean
monthly_seasonal_mean2<-monthly_seasonal_mean[,-11,]   #Delete September#
monthly_seasonal_mean2<-monthly_seasonal_mean2[,-11,]   #Delete October#
monthly_seasonal_mean2<-monthly_seasonal_mean2[,-14,]  #Delete Autumn#

Data_matrix<-matrix(monthly_seasonal_mean2,1600,13*6)           #

#colnames(monthly_seasonal_mean)<-c("November","December","January","February","March","April","May","June","July","August","September","October","Winter","Spring","Summer","Autumn")

#Variable_name<-c("cy","gsl","dps","pr","sfcwind","rsds","tasmax","tasmin")


Season_name<-c("November","December","January","February","March","April","May","June","July","August","Winter","Spring","Summer")

Variable_name<-c("dps","pr","sfcwind","rsds","tasmax","tasmin")

names_v<-c()
for (i in 1:6)
    {  for(j in 1:13)
  {   names_v<-c(names_v,paste(Variable_name[i], Season_name[j], sep="_"))}
 }


Data_All<-cbind(cy_gsl_FR[,1],Data_matrix)

names_v<-c("Yields",names_v)
colnames(Data_All)<-names_v

Data_All<-data.frame(Data_All)
   
   
   ################################
#############SEASONAL_MATRIX#######################àà

x2<-matrix()
a<-12
b<-14
for(i in 1:6)
{
 x1<-Data_All[,a:b]
 a<-b+11
 b<-a+2
 x2<-cbind(x2,x1)  
    
}

x2[,1]<-Data_All[,1]
Season_matrix<-x2

#save(Season_matrix,file='Season_matrix.RData')


########## ONLY MONTHLY DATA###########

x3<-matrix()
a<-2
b<-11
for(i in 1:6)
{
 x1<-Data_All[,a:b]
 a<-b+4
 b<-a+9
 x3<-cbind(x3,x1)  
    
}

x3[,1]<-Data_All[,1]
Monthly_matrix<-x3


############### CHOOSE MATRIX ###########################

 #CH<-1
#for ( CH in 1){
   
#if(CH==1) {Data<-Monthly_matrix }
#if(CH==2) {Data<-Season_matrix }
#if(CH==3) {Data<-Data_All }

Data0<-Season_matrix
 

Summary_m<-matrix(,7,3)
colnames(Summary_m) <-c("0.025","0.05","0.1")
rownames(Summary_m) <-c("Lowyields","Time","MisClassification","AUC","AIC","Spec","Sens")


################### PRESELECTION_ Linear correlation ###########################


                     #####  Standardized data ############
                     
Data <- apply(Data0,2,scale)

corr<-cor(Data)[,1]
Correlation<-sort(abs(corr),decreasing=T)
Final_variable<-Correlation[which(Correlation>0.30)]
Final_variable<- Final_variable [-1]
if(length(Final_variable)>14){Final_variable<-Final_variable[1:14]}



#######################     Predictors matrix    ###########################


Predictors<-names(Final_variable)

j<-1
x<-c(NA)
for(i in 1:length(Predictors))
{
x[i]<-which(colnames(Data)==Predictors[i])
}
 
Data2<-Data[,x]


######################
for (gg in 1:3)
{

quant<- c(0.025,0.05,0.1)
# transform yield data to binary data: severe loss, no severe loss
low_yield <- quantile(Data[,1],quant[gg])
# ind <- cy_gsl_FR[which(cy_gsl_FR[,1]<low_yield),1]
ind <- which(Data[,1]<low_yield)
cy <- Data[,1]
cy[ind] <- 0
cy[-ind] <- 1
#low_yield_FR <- cy_gsl_FR[ind,1]
Summary_m[1,gg]<- quantile(Data0[,1],quant[gg])

Data_test<-cbind(Data2,cy)
#Xy<-data.frame(Data2,y=Data[,1])

###################
################### BEST GLM ####################################

#########################


# Split data into training and testing data set
set.seed(1994)
training_indices <- sample(1:1600, size = floor(1600*0.6))
testing_indices <- (1:1600)[-sort(training_indices)]
Training_Data <- Data_test[training_indices,]
Testing_Data <- Data_test[testing_indices,]

Training_Data<-data.frame(Training_Data)
Testing_Data<-data.frame(Testing_Data)
###################     MODEL 0   ######################################


XY<-Training_Data
#X <- XY[,1:(ncol(XY)-1)]

tim<-system.time(
{
model_CV<-bestglm(Xy = XY, IC = "CV", CVArgs = list(Method = "HTF", K = 10, REP =1000), 
        family=binomial)
 
})

Summary_m[2,gg]<- tim[3]
       
#model_BIC<-bestglm(Xy = XY, IC = "BIC", CVArgs = list(Method = "HTF", K = 10, REP =10), 
        #family=binomial)

#model_AIC<-bestglm(Xy = XY, IC = "AIC", CVArgs = list(Method = "HTF", K = 10, REP =10), 
        #family=binomial)

if(gg==1){Summary_CV_025<-summary(model_CV$BestModel)
          Correlation_025<-Correlation
          Final_variable_025<-Predictors
          tim_025<-tim}

if(gg==2){Summary_CV_05<-summary(model_CV$BestModel)
          Correlation_05<-Correlation 
          Final_variable_05<-Predictors
          tim_05<-tim}

if(gg==3){Summary_CV_01<-summary(model_CV$BestModel)
          Correlation_01<-Correlation 
          Final_variable_01<-Predictors
          tim_01<-tim}

#Summary_BIC<-summary(model_BIC$BestModel)
#Summary_AIC<-summary(model_AIC$BestModel)


######

# Test data set

NR<-ncol(Testing_Data)-1
fitted.results <- predict(model_CV$BestModel,Testing_Data[1:NR],type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError <- mean(fitted.results != Testing_Data$cy)
print(paste('Accuracy',1-misClasificError))
mean(fitted.results != Testing_Data$cy)




### ROC ####

#p <- predict(model1, newdata=subset(testing,select=c(1,2,3)), type="response")
pr <- prediction(fitted.results, Testing_Data$cy)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)


auc <- performance(pr, measure = "auc")
auc@y.values[[1]]


# Confusion matrix ####
# a) Confusion matrix manually calculated

obs_pred <- cbind(Testing_Data$cy,fitted.results)
tp <- sum(rowSums(obs_pred)==2)
tn <- sum(rowSums(obs_pred)==0)
fp <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
fn <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
con_tab1<- matrix(c(tp,fn,fp,tn),nrow=2,ncol=2)
con_tab1b <- con_tab1
colnames(con_tab1) <- c('Actual TRUE','Actual FALSE');rownames(con_tab1) <- c('Predicted TRUE','Predicted FALSE')
con_tab1[,1] <- con_tab1[,1]/ sum(Testing_Data$cy==1)
con_tab1[,2] <-con_tab1[,2]/sum(Testing_Data$cy==0)
# b) Confusion matrix from package InformationValue
con_tab <- InformationValue::confusionMatrix(Testing_Data$cy,fitted.results)

# Sensitivity and Specificity
# a) Manually calculated
spec1 <- tn/(tn+fp) 
sens1 <- tp/(tp+fn) 


sens <-InformationValue::sensitivity(Testing_Data$cy,fitted.results)
spec<-InformationValue::specificity(Testing_Data$cy,fitted.results)


          ####### fill Summary ##########

Summary_m[3,gg]<- misClasificError
Summary_m[4,gg]<- auc@y.values[[1]]
Summary_m[5,gg]<- model_CV$BestModel$aic
Summary_m[6,gg]<- spec
Summary_m[7,gg]<- sens


################## MULTI COLLINEARITY ############################################
# multicollinearity (https://rpubs.com/ranvirkumarsah/LR), https://rpubs.com/andreasme/mlw
                      
M0_VIF<-vif(model_CV$BestModel)                # VIF must be <10 
M0_VIFtolerance<-1/vif(model_CV$BestModel)     #tolerance < 0,1 grave problem ; <0.2 potential problem
M0_VIF_mean<-mean(vif(model_CV$BestModel))     # MEAN must be circa 1, if it is much more than 1 -> problems


VIF_matrix<-matrix(M0_VIF,M0_VIFtolerance,M0_VIF_mean)
 
if(gg==1){M0_VIF_025<-M0_VIF
          M0_VIFtolerance_025<-M0_VIFtolerance
          M0_VIF_mean_025<- M0_VIF_mean
          con_tab_025<-con_tab 
          prf_025<-prf}
          
if(gg==2){M0_VIF_05<-M0_VIF
          M0_VIFtolerance_05<-M0_VIFtolerance
          M0_VIF_mean_05<- M0_VIF_mean 
          con_tab_05<-con_tab
          prf_05<-prf}

if(gg==3){M0_VIF_1<-M0_VIF
          M0_VIFtolerance_1<-M0_VIFtolerance
          M0_VIF_mean_1<- M0_VIF_mean
          con_tab_1<-con_tab
          prf_1<-prf}
 
 

#write.table(Summary_m,"Summary_monthly_TH.csv",sep=" ")


#save(Summary_CV_01,Summary_CV_025,Summary_CV_05,Correlation,Final_variable_01,Final_variable_025,Final_variable_05, file= "GLM_Linear_corr_Monthly.RData")

#save(Summary_CV_01,Summary_CV_025,Summary_CV_05,Correlation,Final_variable_01,Final_variable_025,Final_variable_05, con_tab_1, con_tab_05,con_tab_025,prf_1,prf_05,prf_025,tim_01,tim_025,tim_05,Predictors,file= "GLM_Linear_corr_Season_REP1000.RData")

#setwd('C:/Users/39349/Documents/DAMOCLES/MODEL bestglm/Monthly_data')
#plot(prf_025, main="Monthtly data REP 1000")
#plot(prf_05,add=T,col="red")
#plot(prf_1,add=T,col="green")
#legend("bottomright",legend=c("TH 05","TH 025","TH 1"), col=c("red","black","green"), lty=1)
 
      
      
        ################## END MODEL 0 ##################
                          
                          
                           ##########
      ################# Interaction PART ###############################

Summary_INT<-matrix(,7,3)
colnames(Summary_INT) <-c("0.025","0.05","0.1")
rownames(Summary_INT) <-c("Lowyields","Time","MisClassification","AUC","AIC","Spec","Sens")



#if(gg==1){ model_CV


#########   Selecting the 4 main drivers from the model  #####################

Main_drivers<-  model_CV$BestModel$coefficients
Names_main_d<-names(Main_drivers[2:5])

j<-1
x1<-c(NA)
for(i in 1:length(Names_main_d))
{
x1[i]<-which(colnames(Data)==Names_main_d[i])
}
 
Data_maind<-Data[,x1]


############## ALL COMBINATIONS OF MAIN DRIVERS #####################

df1<-Data_maind
n=ncol(df1)
combb=combn(n,2)
#combb=cbind(combb, sapply(1:n, function(i) rep(i,2)))
res=apply(df1, 1, function(x) { apply(combb, 2, function(y) prod(x[y])) })
Interaction<-t(res)
v_name<-colnames(df1)
Name_interact<-apply(combb, 2, function(y) (v_name[y]))


nc<-paste(Name_interact[1,], Name_interact[2,], sep="-")

colnames(Interaction)<-nc

New_data<-cbind(Data_maind,Interaction,cy)
New_data<-data.frame(New_data)

Training_Data_INT <- New_data[training_indices,]
Testing_Data_INT <- New_data[testing_indices,]

Training_Data_INT<-data.frame(Training_Data_INT)

###################     MODEL INTERACTION   ####################################

#NCT<-ncol(Training_Data_INT)-1
XY<-Training_Data_INT
#X <- XY[,1:(ncol(XY)-1)]


tim_INT<-system.time(
{
model_INT<-bestglm(Xy = XY, IC = "CV", CVArgs = list(Method = "HTF", K = 10, REP =1000), 
        family=binomial)
}) 


Summary_INT[2,gg]<- tim_INT[3]
       
#model_BIC<-bestglm(Xy = XY, IC = "BIC", CVArgs = list(Method = "HTF", K = 10, REP =10), 
        #family=binomial)

#model_AIC<-bestglm(Xy = XY, IC = "AIC", CVArgs = list(Method = "HTF", K = 10, REP =10), 
        #family=binomial)

if(gg==1){Summary_INT_CV_025<-summary(model_INT$BestModel)
          Correlation_INT_025<-Correlation_INT
          Final_variable_INT_025<-Predictors_INT
          tim_INT_025<-tim_INT}

if(gg==2){Summary_INT_CV_05<-summary(model_INT$BestModel)
          Correlation_INT_05<-Correlation_INT 
          Final_variable_INT_05<-Predictors_INT
          tim_INT_05<-tim_INT}

if(gg==3){Summary_INT_CV_01<-summary(model_INT$BestModel)
          Correlation_INT_01<-Correlation_INT 
          Final_variable_INT_01<-Predictors_INT
          tim_INT_01<-tim_INT}

#Summary_BIC<-summary(model_BIC$BestModel)
#Summary_AIC<-summary(model_AIC$BestModel)


######

# Test data set

NR<-ncol(Testing_Data_INT)-1
fitted.results <- predict(model_INT$BestModel,Testing_Data_INT[1:NR],type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError_INT <- mean(fitted.results != Testing_Data_INT$cy)
print(paste('Accuracy',1-misClasificError_INT))
mean(fitted.results != Testing_Data_INT$cy)




### ROC ####

#p <- predict(model1, newdata=subset(testing,select=c(1,2,3)), type="response")
pr_INT <- prediction(fitted.results, Testing_Data_INT$cy)
prf_INT <- performance(pr_INT, measure = "tpr", x.measure = "fpr")
#plot(prf_INT)

figure()
plot(prf_INT, main=paste0("Interaction ",quant[gg]))
end()
auc_INT <- performance(pr_INT, measure = "auc")
auc_INT@y.values[[1]]


# Confusion matrix ####
# a) Confusion matrix manually calculated

obs_pred <- cbind(Testing_Data_INT$cy,fitted.results)
tp <- sum(rowSums(obs_pred)==2)
tn <- sum(rowSums(obs_pred)==0)
fp <- sum(obs_pred[,1]==0 & obs_pred[,2]==1)
fn <- sum(obs_pred[,1]==1 & obs_pred[,2]==0)
con_tab1<- matrix(c(tp,fn,fp,tn),nrow=2,ncol=2)
con_tab1b <- con_tab1
colnames(con_tab1) <- c('Actual TRUE','Actual FALSE');rownames(con_tab1) <- c('Predicted TRUE','Predicted FALSE')
con_tab1[,1] <- con_tab1[,1]/ sum(Testing_Data_INT$cy==1)
con_tab1[,2] <-con_tab1[,2]/sum(Testing_Data_INT$cy==0)
# b) Confusion matrix from package InformationValue
con_tab <- InformationValue::confusionMatrix(Testing_Data_INT$cy,fitted.results)

# Sensitivity and Specificity
# a) Manually calculated
spec1 <- tn/(tn+fp) 
sens1 <- tp/(tp+fn) 


sens_INT <-InformationValue::sensitivity(Testing_Data_INT$cy,fitted.results)
spec_INT<-InformationValue::specificity(Testing_Data_INT$cy,fitted.results)


          ####### fill Summary ##########

Summary_INT[3,gg]<- misClasificError_INT
Summary_INT[4,gg]<- auc_INT@y.values[[1]]
Summary_INT[5,gg]<- model_INT$BestModel$aic
Summary_INT[6,gg]<- spec_INT
Summary_INT[7,gg]<- sens_INT


################## MULTI COLLINEARITY ############################################
# multicollinearity (https://rpubs.com/ranvirkumarsah/LR), https://rpubs.com/andreasme/mlw
                      
M0_VIF_INT<-vif(model_INT$BestModel)                # VIF must be <10 
M0_VIFtolerance_INT<-1/vif(model_INT$BestModel)     #tolerance < 0,1 grave problem ; <0.2 potential problem
M0_VIF_mean_INT<-mean(vif(model_INT$BestModel))     # MEAN must be circa 1, if it is much more than 1 -> problems

M0_INT_VIF_1
M0_INT_VIFtolerance_1
M0_INT_VIF_mean_1

#VIF_matrix<-matrix(M0_VIF,M0_VIFtolerance,M0_VIF_mean)
 
if(gg==1){M0_INT_VIF_025<-M0_VIF_INT
          M0_INT_VIFtolerance_025<-M0_VIFtolerance_INT
          M0_INT_VIF_mean_025<- M0_VIF_mean_INT
          con_tab_INT_025<-con_tab_INT 
          prf_INT_025<-prf_INT}
          
if(gg==2){M0_INT_VIF_05<-M0_VIF_INT
          M0_INT_VIFtolerance_05<-M0_VIFtolerance_INT
          M0_INT_VIF_mean_05<- M0_VIF_mean_INT 
          con_tab_INT_05<-con_tab_INT
          prf_INT_05<-prf_INT}

if(gg==3){M0_INT_VIF_1<-M0_VIF_INT
          M0_INT_VIFtolerance_1<-M0_VIFtolerance_INT
          M0_INT_VIF_mean_1<- M0_VIF_mean_INT
          con_tab_INT_1<-con_tab_INT
          prf_INT_1<-prf}
 
}



###############END OF THE CODE ##################################
#####################################################################
     ########## ALL SEASONAL INTERACTIONS ######################


df1<-Season_matrix[-1]
n=ncol(df1)
combb=combn(n,2)
#combb=cbind(combb, sapply(1:n, function(i) rep(i,2)))
res=apply(df1, 1, function(x) { apply(combb, 2, function(y) prod(x[y])) })
Interaction_seas<-t(res)
v_name<-colnames(df1)
Name_interact<-apply(combb, 2, function(y) (v_name[y]))


nc<-paste(Name_interact[1,], Name_interact[2,], sep="-")
colnames(Interaction_seas)<-nc


Data_seas<-Interaction_seas

#cc<-c(0.45,0.3)
Corr_seas<-cor(Data_seas)[,1]
Correlation_seas<-sort(abs(Corr_seas),decreasing=T)

#ALL_corr<- c(Correlation,Correlation_INT_seas)


Final_variable_seas<-Correlation_seas[which(Correlation_seas>0.45)]
Final_variable_seas<- Final_variable_seas [-1]


Final_variable_INT<- names(Final_variable_seas[1:10])


#Names_main_d<-names(Main_drivers[2:5])

j<-1
xx1<-c(NA)
for(i in 1:length(Final_variable_INT))
{
xx1[i]<-which(colnames(Data_seas)==Final_variable_INT[i])
}
 
Interaction<-Data_seas[,xx1]

New_data<-cbind(Interaction,Data_maind,cy)

Training_Data_int <- New_data[training_indices,]
Testing_Data_int <- New_data[testing_indices,]


###################     MODEL INTERAC   ########################################

XY<-Training_Data_int
#X <- XY[,1:(ncol(XY)-1)]

tic("model")
model_INT<-bestglm(Xy = XY, IC = "CV", CVArgs = list(Method = "HTF", K = 10, REP =10), 
        family=binomial)
toc() 


   
if(CH==1) {Model_monthly_CV<-model_CV
           Model_monthly_BIC<-model_BIC
           Model_monthly_AIC<-model_AIC
           Predictors_monthly<- Predictors
           Corr_monthly<-Correlation
          }
if(CH==2) {Model_season<-model
          Predictors_season<- Predictors
          Corr_season<-Correlation}
if(CH==3) {Model_All<-model
          Predictors_All<- Predictors
          Corr_All<-Correlation}



save(Model_All,Model_monthly,Model_season,Predictors_All,Predictors_monthly,Predictors_season, Corr_All,Corr_monthly, Corr_season, file= "GLM_Linear_corr_0.3.RData")



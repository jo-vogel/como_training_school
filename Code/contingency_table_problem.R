# Contigency table problem

# Actual values
# Transform categories to numbers
TrainSet[TrainSet[,7]=="acc",8] <- 1
TrainSet[TrainSet[,7]=="good",8] <- 2
TrainSet[TrainSet[,7]=="unacc",8] <- 3
TrainSet[TrainSet[,7]=="vgood",8] <- 4

TrainSet[,8] # actuals

# Predictions
# predicted values are assigned to a category by column
pred_val <- apply(model2[["votes"]],1,which.max) # model predictions, contingency table in model is in line with this
diffs <- TrainSet[,8]-pred_val # 
sum(diffs!=0) # 31 misclassifications, as in the contingency table of the model
model2$confusion # contingency table

# Second set of predictions
predTrain # also model predictions
# predTrain_ext <- matrix(NA,nrow=length(predTrain),ncol=2)
# predTrain_ext[predTrain=="acc",2] <- 1
# predTrain_ext[predTrain=="good",2] <- 2
# predTrain_ext[predTrain=="unacc",2] <- 3
# predTrain_ext[predTrain=="vgood",2] <- 4
# predTrain_ext[,1] <- predTrain
predTrain_vec <- vector("numeric",length=length(predTrain)) # easier way
predTrain_vec[1:1208] <- predTrain

diffs2 <- TrainSet[,8] - predTrain_vec
plot(diffs2)
sum(diffs2!=0) # predTrain is identical with Trainset[,8], no misclassifications


# Why are there two different contigency tables? The second one seems wrong, so can we trust the validation table, which is calculated in the same way?

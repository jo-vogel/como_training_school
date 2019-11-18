#it is not working#

tic()
no_cores <- detectCores() / 2 - 1
cl<-makeCluster(no_cores)
clusterEvalQ(cl, {
  library(glinternet)
  library(dplyr)
}) # parallelisation has own environment, therefore some packages and variables need be loaded again
registerDoParallel(cl)

res<-vector(mode = "list", length = pix_num)
res[[1]]<-try(exp(12),silent = T)  

   #Pauline: i cannot go up to (pix_num+1) if you're calling Training_Data [[i]] (list of size pix_num)
   
model_CV<- foreach (i=2:(pix_num+1),.multicombine=TRUE,.packages="bestglm") %dopar% {
  
  if( class(res[[i-1]])!="try-error") #Pauline: res[[i-1]] cannot be called in a paralellized loop: maybe res[[i-1]] was not created yet. In a paralellized loop we don't know in what order the results are productes
  {        
  res<-try(bestglm(Xy = data.frame(Training_Data [[i]]), IC = "CV", CVArgs = list(Method = "HTF", K = 10, REP =1000), 
        family=binomial),silent=T)} else{0} 
   
  } 
  
      
stopCluster(cl)
toc()





####### Modif from Pauline #######

tic()
no_cores <- detectCores() / 2 - 1
cl<-makeCluster(no_cores)
clusterEvalQ(cl, {
  library(bestglm)
  library(glinternet)
  library(dplyr)
}) # parallelisation has own environment, therefore some packages and variables need be loaded again
registerDoParallel(cl)

res <- list()

   
model_CV<- foreach (i=1:pix_num,.multicombine=TRUE) %dopar% {
  
  TRY <- try(bestglm(Xy = data.frame(Training_Data[[i]]), IC = "CV", CVArgs = list(Method = "HTF", K = 10, REP =1000), family=binomial), silent=T)}
  
  if( class(TRY)!="try-error"){#there is no error at this pixel, thus TRY is the result of bestglm
  
  res[[i]] <- TRY
  
  } else { #bestglm doesn't give any result, put a messsage instead of a result, so that the loop doesn't crash
  
  res[[i]] <- "problem at this pixel"
  
  }#end ifelse
  
}#end foreach
      
stopCluster(cl)
toc()
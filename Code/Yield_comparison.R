# Creator: Johannes Vogel

# to be run with Problematic_pixels.r

# Comparison of yield in areas with problematic pixels compared to whole northern hemisphere

load("D:/user/vogelj/Data/Group project Como/extremeindices_and_monthlymeteovar.Rdata")
yield <- Data_non_standardized$yield

yield_prob <- yield[no_neg,] # Yield of problematic pixels (no negative predictions)
yield_non_prob <- yield[!no_neg,] # Yield of problematic pixels (no negative predictions)

yield_sd <- apply(yield,1,sd,na.rm=T)
yield_prob_sd <- apply(yield_prob,1,sd,na.rm=T)
yield_non_prob_sd <- apply(yield_non_prob,1,sd,na.rm=T)

yield_median <- apply(yield,1,median,na.rm=T)
yield_prob_median <- apply(yield_prob,1,median,na.rm=T)
yield_non_prob_median <- apply(yield_non_prob,1,median,na.rm=T)

# hist(yield_median)
# hist(yield_prob_median)
# summary(yield_median)
# summary(yield_prob_median)
par(mar=c(5, 4, 4, 2) + 0.1)
# boxplot(yield_median,yield_prob_median, names=c("Median yield (all pixels)","Median yield (problematic pixels)"), col="lightblue",main="Glinternet incl. extreme indicators")
# boxplot(yield_median,yield_prob_median, names=c("Median yield (all pixels)","Median yield (problematic pixels)"), col="lightblue",main="Glmnet incl. extreme indicators")
# make a data frame: one column with all observations another column with yes/no (bad or not), then proceed as in the example
yield_median_df <- data.frame(yield_median,no_neg)
names(yield_median_df)[2] <- "No negative predictions"
ggplot(data=yield_median_df)+
  geom_violin(aes(`No negative predictions`,yield_median),fill="lightblue",title="violin")+
  labs( x = "No negative predictions", y = "Median yield",
        title ="Median yield for pixels with (left) and without (right) \n negative prections for glinternet")
      # title ="Median yield for pixels with (left) and without (right) \n negative prections for glmnet")

yield_sd_df <- data.frame(yield_sd,no_neg)
ggplot(data=yield_sd_df)+
  geom_violin(aes(no_neg,yield_sd),fill="lightblue",title="violin")+
  labs( x = "No negative predictions", y = "Standard deviation of yield",
        title ="Standard deviation of yield for pixels with (left) and without (right) \n negative prections for glinternet")
      # title ="Standard deviation of yield for pixels with (left) and without (right) \n negative prections for glmnet")
# boxplot(yield_sd,yield_prob_sd, col="lightblue",names=c("Standard deviation of yield \n (all pixels)","Standard deviation of yield \n (problematic pixels)"),main="Glinternet incl. extreme indicators")
# boxplot(yield_sd,yield_prob_sd, col="lightblue",names=c("Standard deviation of yield \n (all pixels)","Standard deviation of yield \n (problematic pixels)"),main="Glmnet incl. extreme indicators")


# compare how far 5th percentile deviates from the median

yield_5perc <- apply(yield,1,quantile,0.05,na.rm=T) # same as "low_yield", but non-standardised
yield_5perc_prob <- apply(yield_prob,1,quantile,0.05,na.rm=T)
yield_5perc_non_prob <- apply(yield_non_prob,1,quantile,0.05,na.rm=T)

# boxplot(yield_5perc,yield_5perc_prob, col="lightblue",names=c("5th perc. yield \n (all pixels)","5th perc. yield \n (problematic pixels)"),main="Glinternet incl. extreme indicators")
# boxplot(yield_5perc,yield_5perc_prob, col="lightblue",names=c("5th perc. yield \n (all pixels)","5th perc. yield \n (problematic pixels)"),main="Glmnet incl. extreme indicators")
yield_5perc_df <- data.frame(yield_5perc,no_neg)
ggplot(data=yield_5perc_df)+
  geom_violin(aes(no_neg,yield_5perc),fill="lightblue",title="violin")+
  labs( x = "No negative predictions", y = "5th percentile yield",
        title ="5th percentile yield for pixels with (left) and without (right) \n negative prections for glinternet")
      # title ="5th percentile yield for pixels with (left) and without (right) \n negative prections for glmnet")
  

par(mar=c(5,10,3,2))
barplot(c(median(yield_prob_median),median(yield_non_prob_median),median(yield_5perc_prob),median(yield_5perc_non_prob)), col="lightblue"
        ,horiz=T,names=(c("Median yield \n (problematic pixels)","Median yield \n (Non-problematic pixels)","5th perc. yield \n (problematic pixels)","5th perc. yield \n (Non-problematic pixels)")),las=1,main="Median for glinternet incl. extreme indicators")
          # ,horiz=T,names=(c("Median yield \n (problematic pixels)","Median yield \n (Non-problematic pixels)","5th perc. yield \n (problematic pixels)","5th perc. yield \n (Non-problematic pixels)")),las=1,main="Median for glmnet incl. extreme indicators")





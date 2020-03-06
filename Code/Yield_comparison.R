# Creator: Johannes Vogel

# to be run with Problematic_pixels.r

# Comparison of yield in areas with problematic pixels compared to whole northern hemisphere


yield_prob <- yield[no_neg,] # Yield of problematic pixels (no negative predictions)

yield_sd <- apply(yield,1,sd,na.rm=T)
yield_prob_sd <- apply(yield_prob,1,sd,na.rm=T)

yield_median <- apply(yield,1,median,na.rm=T)
yield_prob_median <- apply(yield_prob,1,median,na.rm=T)

# hist(yield_median)
# hist(yield_prob_median)
# summary(yield_median)
# summary(yield_prob_median)
par(mar=c(5, 4, 4, 2) + 0.1)
# boxplot(yield_median,yield_prob_median, names=c("Median yield (all pixels)","Median yield (problematic pixels)"), col="lightblue",main="Glinternet incl. extreme indicators")
# boxplot(yield_sd,yield_prob_sd, col="lightblue",names=c("Standard deviation of yield \n (all pixels)","Standard deviation of yield \n (problematic pixels)"),main="Glinternet incl. extreme indicators")
boxplot(yield_median,yield_prob_median, names=c("Median yield (all pixels)","Median yield (problematic pixels)"), col="lightblue",main="Glmnet incl. extreme indicators")
boxplot(yield_sd,yield_prob_sd, col="lightblue",names=c("Standard deviation of yield \n (all pixels)","Standard deviation of yield \n (problematic pixels)"),main="Glmnet incl. extreme indicators")



# compare how far 5th percentile deviates from the median

# yield_5perc <- apply(yield,1,quantile,0.05,na.rm=T) # same as "low_yield"
low_yield_prob <- apply(yield_prob,1,quantile,0.05,na.rm=T)

# boxplot(low_yield,low_yield_prob, col="lightblue",names=c("5th perc. yield \n (all pixels)","5th perc. yield \n (problematic pixels)"),main="Glinternet incl. extreme indicators")
boxplot(low_yield,low_yield_prob, col="lightblue",names=c("5th perc. yield \n (all pixels)","5th perc. yield \n (problematic pixels)"),main="Glmnet incl. extreme indicators")

par(mar=c(5,9,3,2))
barplot(c(median(yield_prob_median),median(yield_median),median(low_yield_prob),median(low_yield)), col="lightblue"
        # ,horiz=T,names=(c("Median yield \n (problematic pixels)","Median yield \n (all pixels)","5th perc. yield \n (problematic pixels)","5th perc. yield \n (all pixels)")),las=1,main="Glinternet incl. extreme indicators")
          ,horiz=T,names=(c("Median yield \n (problematic pixels)","Median yield \n (all pixels)","5th perc. yield \n (problematic pixels)","5th perc. yield \n (all pixels)")),las=1,main="Glmnet incl. extreme indicators")





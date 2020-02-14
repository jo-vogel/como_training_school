library(ggplot2)

load("D:/user/vogelj/Group_project/Data/lambda_values_Lasso_glmnet_monthly_xtrm_indices.Rdata")
lambda_gl <- read.csv("D:/user/vogelj/Group_project/Data/glinternet_lasso_without_interactions_lambdaHat1Std.csv")
# in folder Models on gdrive

# make sure the coordinates are the same




# Load model output ####
########################

path_to_NH_files <- "D:/user/vogelj/Data/Group project Como"
source('./Code/Lasso_interact_global_preparation.R') # monthly data


# without interactions
# Models/Lasso (glinternet)/LASSO_without_interactions/cv_fit_no_int.RData
load("D:/user/vogelj/Group_project/Code/Workspaces/cv_fit_no_int.RData") # monthly model without interactions




# make a comparison plot
world <- map_data("world")

# lambda_curr <- lambda_gl$x # lambda from glinternet
lambda_curr <- LAMBDAS[[2]] # lambda from glmnet

DF_lambda <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], lamb = lambda_curr)

ggplot(data = DF_lambda, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.1) +
  geom_point(shape=15, aes(color=lamb),size=0.7) +
  scale_color_gradient2(limits=c(min(lambda_curr,na.rm=T),max(lambda_curr,na.rm=T)),midpoint=min(lambda_curr,na.rm=T)+(max(lambda_curr,na.rm=T)-min(lambda_curr,na.rm=T))/2,
                        low = "yellow", mid = "red3", high = "black") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  ylab("Lat (°N)") +
  xlab("Lon (°E)") +
  coord_fixed(xlim = c(-120, 135),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1.3)+
  labs(color="Lambda",
       title = paste("Lambda comparison"),
       subtitle = paste("Bad yield threshold=", threshold,
                        ", segregation threshold=", segreg_th, sep = ""))+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
  X11(width = 20, height = 7)
ggsave(file="D:/user/vogelj/Group_project/Output/Plots/lambda_map.png")

boxplot(lambda_gl$x,LAMBDAS[[2]],names=c("glinternet","glmnet"),ylab="Lambda (1 Std)")

boxplot(lambda_gl$x,LAMBDAS[[2]],LAMBDAS[[1]],names=c("glinternet (1Std)","glmnet (1Std)", "glmnet (Min)"))


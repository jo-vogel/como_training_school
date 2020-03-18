get_firstcoeffs <- function(coeff, nb_of_coeff=1){
  N <- nb_of_coeff
  coeff_name_date <- matrix(data = NA, nrow = N, ncol = 2)
  row.names(coeff_name_date) <- as.character(1:N)
  colnames(coeff_name_date) <- c("var", "month")
  
  for (ind in 1:N) {
    coeff_names <- rownames(coeff)[sort(abs(as.numeric(coeff)), decreasing = T,
                                        index.return=T)$ix[-which(sort(abs(as.numeric(coeff)), 
                                                                       decreasing = T,
                                                                       index.return=T)$ix==1)]][ind]
    if(substr(coeff_names, start = 1, stop = 3)=="dtr"){
      coeff_name_date[ind,1] <- "dtr"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="frs"){
      coeff_name_date[ind,1] <- "frs"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="txx"){
      coeff_name_date[ind,1] <- "txx"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="tnn"){
      coeff_name_date[ind,1] <- "tnn"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="rx5"){
      coeff_name_date[ind,1] <- "rx5"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="tx90p"){
      coeff_name_date[ind,1] <- "tx90p"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="tn10p"){
      coeff_name_date[ind,1] <- "tn10p"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 4)=="tmax"){
      coeff_name_date[ind,1] <- "tmax"
      coeff_name_date[ind,2] <- substr(coeff_names, start = 6, stop = 11)
    }
    
    if(substr(coeff_names, start = 1, stop = 1)=="v"){
      coeff_name_date[ind,1] <- "vpd"
      coeff_name_date[ind,2] <- substr(coeff_names, start = 5, stop = 10)
    }
    
    if(substr(coeff_names, start = 1, stop = 1)=="p"){
      coeff_name_date[ind,1] <- "pr"
      coeff_name_date[ind,2] <- substr(coeff_names, start = 4, stop = 9)
    }
    
  }#end for ind
  
  return(coeff_name_date)
  
}#end fct get_firstcoeffs
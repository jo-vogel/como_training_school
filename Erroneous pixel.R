# Erroneous pixel


for (i in 1:930){
   vect[i] <- sum(is.nan(Model_data_stand[i,,]))
}
plot(vect) # pixel 10 is erroneous

pix_in[10]

plot(nh_variables[[5]][180,12,5,]) # For this pixel precipitation in May is always 0, gives NaN (and therefore error in the model) if you standardise it
plot(monthly_data3dim[3700,5,]) # same pixel

# Conversion of dimensions 
# example: 3700 corresponds to 180,12
pixel <- 3700
numrow <- 3700%%320
numcol <- 3700%/%320
numcol*320+numrow
newpixel <- c(numrow,numcol+1)

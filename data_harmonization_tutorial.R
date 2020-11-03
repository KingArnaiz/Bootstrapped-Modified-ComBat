#*************************************************************************************************************************************
#**************************************** TESTING THE BOOTSTRAPPED MODIFIED COMBAT ***************************************************
#*************************************************************************************************************************************
#@author  Ronrick Da-ano <ronrickarnaiz@gmail.com>
#@version 2.0, 03/11/2020
#@since   R version (3.4.4).
#_____________________________________________________________________________________________________________________________________
#  LIBRARY DECLARATION
#-------------------------------------------------------------------------------------------------------------------------------------
# install sva package from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("sva")
library(sva)
library(boot)
source("M_ComBat.R")
set.seed(2020)
start.time <- Sys.time()
#-------------------------------------------------------------------------------------------------------------------------------------
#load the sample data
#_____________________________________________________________________________________________________________________________________
data_1 <- read.csv("sample_data.csv")

#-------------------------------------------------------------------------------------------------------------------------------------
######################################################## Modified_ComBat ############################################################
#_____________________________________________________________________________________________________________________________________

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Separate the data into clinical and radiomic feature!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data_1_sub <- data_1[,c(1:4)]  #subset with the clinical features only
data_2_sub <- data_1[,-c(1:4)] #subset of the data with radiomics fetures only
newdata_1 <- as.data.frame(data_2_sub) #make the radiomics data a data frame format
newdata_1 <- as.matrix(newdata_1) #convert the data frame into a matrix
newdata_1 <- t(newdata_1) #transpose the matrix format

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declare the scanner / center value !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
batch_1 <- data_1_sub$Institution #get the scanner/instution value to refer the harmonization
batch_new_1 <- as.integer(batch_1) #convert the column as integer/numeric type
mod_1 <- matrix(rep(1,length(batch_new_1)),length(batch_new_1),1) #create a baseline model for the batches

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Harmonize the data based on the center of your reference !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data_harmonized_1 <- M_ComBat(newdata_1 ,batch_new_1 , center=1, mod_1)  # perform  M-ComBat centered at batch 1 (center=1)

data_harmonized_2 <- M_ComBat(newdata_1, batch_new_1 , center=2, mod_1)  # perform  M-ComBat centered at batch 2 (center=2)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Put back the data into its orginal format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
finalized_data_1 <- t(data_harmonized_1$dat.combat) #obtain only the feture and transpose back
combined_data_1 <- cbind(data_1_sub,finalized_data_1) #combine the original clinical data and the harmonized features

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Save the harmonized overall data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write.csv(combined_data_1, "m_combat.csv", row.names = FALSE)
#_____________________________________________________________________________________________________________________________________


#-------------------------------------------------------------------------------------------------------------------------------------
################################################## Bootstrapped Modified_ComBat ######################################################
#_____________________________________________________________________________________________________________________________________

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! the bootstrap function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
proposed_bootstrap = function(data,formula){  
  return(M.COMBAT(dat=data_harmonized_1$dat.combat, batch_new_1 , center=1, mod_1)$dat.combat)  
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Perform the bootstrap based on the last harmonization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
boot_data <- boot(data=data_harmonized_1, statistic=proposed_bootstrap, R=1000) 

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Put back the data into its orginal format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
finalized_boot_data <- boot_data$t0
finalized_boot_data <- t(finalized_boot_data) #obtain only the feture and transpose back
combined_boot_data <- cbind(data_1_sub,finalized_boot_data) #combine the original clinical data and the harmonized features

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Save the harmonized overall data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write.csv(combined_boot_data, "bm_combat_data.csv", row.names = FALSE)

#----------------------------------------------------------------------------------------------------------------------------------------------
#__________________________________________________________________________________________________________________________________________________________________________________________ 
end.time <- Sys.time()
Algo_time<-end.time - start.time
print("ALGO TIME:")
print(Algo_time)


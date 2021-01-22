## code to prepare `SPOT_6_7` spectral response
# 1- define sensor name
SensorName <-'SPOT_6_7'
Path_SRF <- file.path('data-raw',paste(SensorName,'_Spectral_Response.csv',sep=''))
# 2- read file containing spectral response
message('_____ reading spectral response corresponding to ______')
print(SensorName)
SRFraw <- read.csv(Path_SRF,header = TRUE,sep = '\t')

# 3- identify initial spectral sampling
OriginalBands <- SRFraw[,1]
# 4- identify name of spectral bands for sensor
Spectral_Bands <- colnames(SRFraw)
Spectral_Bands <- Spectral_Bands[-1]
# 5- check if conversion of spctral bands into numeric values
SensorRadiometry <-SRFraw[,-1]
SensorRadiometry <- t(SensorRadiometry)
SPOT_6_7 <- list("Spectral_Response"=SensorRadiometry, "Spectral_Bands"=Spectral_Bands,'OriginalBands'=OriginalBands)

## code to prepare `SPOT_6_7_Spectral_Response` dataset goes here
usethis::use_data(SPOT_6_7,compress = 'xz',overwrite = TRUE)
save(SPOT_6_7,file =  file.path('data',paste(SensorName,'.RData',sep = '')))




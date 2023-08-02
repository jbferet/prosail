## code to prepare `Venus` spectral response
# 1- define sensor name
SensorName <-'Venus'
Path_SRF <- file.path('data-raw',paste(SensorName,'_Spectral_Response.csv',sep=''))
# 2- read file containing spectral response
message('_____ reading spectral response corresponding to ______')
print(SensorName)
SRFraw <- read.csv(Path_SRF,header = TRUE,sep = '\t')

# 3- identify initial spectral sampling
Original_Bands <- SRFraw[,1]
# 4- identify name of spectral bands for sensor
Spectral_Bands <- colnames(SRFraw)
Spectral_Bands <- Spectral_Bands[-1]
# 5- check if conversion of spctral bands into numeric values
SensorRadiometry <-SRFraw[,-1]
SensorRadiometry <- t(SensorRadiometry)
Central_WL <- c(420, 443, 490, 555, 620, 620, 667, 702, 742, 782, 865, 910)
Venus <- list('Spectral_Response' = SensorRadiometry,
              'Spectral_Bands' = Spectral_Bands,
              'Original_Bands' = Original_Bands,
              'Central_WL' = Central_WL)

## code to prepare `Venus_Spectral_Response` dataset goes here
usethis::use_data(Venus,compress = 'xz',overwrite = TRUE)
save(Venus,file =  file.path('data',paste(SensorName,'.RData',sep = '')))




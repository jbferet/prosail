## code to prepare `MODIS` spectral response
# 1- define sensor name
SensorName <-'MODIS'
Path_SRF <- file.path('data-raw',paste0(SensorName,'_Spectral_Response.csv'))
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
Central_WL <- c(659, 865, 470, 555, 1240, 1640, 2130, 415, 443, 490, 531,
                565, 653, 681, 750, 865, 905, 936, 940)
MODIS <- list('Spectral_Response' = SensorRadiometry,
              'Spectral_Bands' = Spectral_Bands,
              'Original_Bands' = Original_Bands,
              'Central_WL' = Central_WL)

## code to prepare `MODIS_Spectral_Response` dataset goes here
usethis::use_data(MODIS,compress = 'xz',overwrite = TRUE)
save(MODIS,file =  file.path('data',paste0(SensorName,'.RData')))

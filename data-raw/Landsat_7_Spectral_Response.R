## code to prepare `Landsat_7` spectral response
# 1- define sensor name
SensorName <-'Landsat_7'
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
Central_WL <- c(477.5, 560, 661.5, 835, 1648, 2204.5, 710)
Landsat_7 <- list('Spectral_Response' = SensorRadiometry,
                  'Spectral_Bands' = Spectral_Bands,
                  'Original_Bands' = Original_Bands,
                  'Central_WL' = Central_WL,
                  'spectral_response' = SensorRadiometry,
                  'spectral_bands' = Spectral_Bands,
                  'original_bands' = Original_Bands,
                  'central_wl' = Central_WL)

## code to prepare `Sentinel_2_Spectral_Response` dataset goes here
usethis::use_data(Landsat_7,compress = 'xz',overwrite = TRUE)
save(Landsat_7,file =  file.path('data',paste0(SensorName,'.RData')))

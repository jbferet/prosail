## code to prepare `Sentinel_2C` spectral response
# 1- define sensor name
SensorName <-'Sentinel_2C'
Path_SRF <- file.path('data-raw', paste0(SensorName,'_Spectral_Response.csv'))
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
Central_WL <- c(489.0, 560.6, 666.5, 707.1, 741.1, 784.7, 834.6, 865.6,
                1612.0, 2191.3)
Sentinel_2C <- list('Spectral_Response' = SensorRadiometry,
                    'Spectral_Bands' = Spectral_Bands,
                    'Original_Bands' = Original_Bands,
                    'Central_WL' = Central_WL,
                    'spectral_response' = SensorRadiometry,
                    'spectral_bands' = Spectral_Bands,
                    'original_bands' = Original_Bands,
                    'central_wl' = Central_WL)

## code to prepare `Sentinel_2C_Spectral_Response` dataset goes here
usethis::use_data(Sentinel_2C,compress = 'xz',overwrite = TRUE)
save(Sentinel_2C,file =  file.path('data', paste0(SensorName, '.RData')))

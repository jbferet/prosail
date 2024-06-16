## code to prepare `Sentinel_2B` spectral response
# 1- define sensor name
SensorName <-'Sentinel_2B'
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
Central_WL <- c(492.1, 559, 664.9, 703.8, 739.1, 779.7, 832.9, 864,
                1610.4, 2185.7)
Sentinel_2B <- list('Spectral_Response' = SensorRadiometry,
                    'Spectral_Bands' = Spectral_Bands,
                    'Original_Bands' = Original_Bands,
                    'Central_WL' = Central_WL)

## code to prepare `Sentinel_2B_Spectral_Response` dataset goes here
usethis::use_data(Sentinel_2B,compress = 'xz',overwrite = TRUE)
save(Sentinel_2B,file =  file.path('data', paste0(SensorName, '.RData')))

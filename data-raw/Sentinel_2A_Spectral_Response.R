## code to prepare `Sentinel_2A` spectral response
# 1- define sensor name
SensorName <-'Sentinel_2A'
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
Central_WL <- c(492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7,
                1613.7, 2202.4)
Sentinel_2A <- list('Spectral_Response' = SensorRadiometry,
                    'Spectral_Bands' = Spectral_Bands,
                    'Original_Bands' = Original_Bands,
                    'Central_WL' = Central_WL)

## code to prepare `Sentinel_2A_Spectral_Response` dataset goes here
usethis::use_data(Sentinel_2A,compress = 'xz',overwrite = TRUE)
save(Sentinel_2A,file =  file.path('data', paste0(SensorName,'.RData')))

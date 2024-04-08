## code to prepare `Pleiades_1A` spectral response
# 1- define sensor name
SensorName <-'Pleiades_1A'
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
Central_WL <- c(488, 551, 645, 832, 635)
Pleiades_1A <- list('Spectral_Response' = SensorRadiometry,
                    'Spectral_Bands' = Spectral_Bands,
                    'Original_Bands' = Original_Bands,
                    'Central_WL' = Central_WL)

## code to prepare `Pleiades_1A_Spectral_Response` dataset goes here
usethis::use_data(Pleiades_1A,compress = 'xz',overwrite = TRUE)
save(Pleiades_1A,file =  file.path('data',paste0(SensorName,'.RData')))

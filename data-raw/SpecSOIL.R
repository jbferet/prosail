# read dataspec file and save it in a dataframe
dataspec    = read.csv2(file = 'data-raw/dataSpec_SOIL_ATM.txt',header = FALSE,sep = '\t')

lambda        = as.numeric(as.character(dataspec[,1]))
Dry_Soil      = as.numeric(as.character(dataspec[,4]))
Wet_Soil      = as.numeric(as.character(dataspec[,5]))
SpecSOIL      =  data.frame('lambda'=lambda,'Dry_Soil'=Dry_Soil,'Wet_Soil'=Wet_Soil)

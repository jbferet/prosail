# read dataspec file and save it in a dataframe
dataspec    = read.csv2(file = 'data-raw/dataSpec_SOIL_ATM.txt',header = FALSE,sep = '\t')

lambda        = as.numeric(as.character(dataspec[,1]))
Direct_Light  = as.numeric(as.character(dataspec[,2]))
Diffuse_Light = as.numeric(as.character(dataspec[,3]))

SpecATM       =  data.frame('lambda'=lambda,'Direct_Light'=Direct_Light,'Diffuse_Light'=Diffuse_Light)

# read dataspec file and save it in a dataframe
dataspec    = read.csv2(file = 'data-raw/SpecPROSPECT.txt',header = FALSE,sep = '\t')
lambda      = as.numeric(as.character(dataspec[,1]))
nrefrac     = as.numeric(as.character(dataspec[,2]))
SAC_CHL     = as.numeric(as.character(dataspec[,3]))
SAC_CAR     = as.numeric(as.character(dataspec[,4]))
SAC_ANT     = as.numeric(as.character(dataspec[,5]))
SAC_BROWN   = as.numeric(as.character(dataspec[,6]))
SAC_EWT     = as.numeric(as.character(dataspec[,7]))
SAC_LMA     = as.numeric(as.character(dataspec[,8]))
SAC_PROT    = as.numeric(as.character(dataspec[,9]))
SAC_CBC     = as.numeric(as.character(dataspec[,10]))

SpecPROSPECT  =  data.frame('lambda'=lambda,'nrefrac'=nrefrac,'SAC_CHL'=SAC_CHL,'SAC_CAR'=SAC_CAR,
                            'SAC_ANT'=SAC_ANT,'SAC_BROWN'=SAC_BROWN,'SAC_EWT'=SAC_EWT,'SAC_LMA'=SAC_LMA,
                            'SAC_PROT'=SAC_PROT,'SAC_CBC'=SAC_CBC)

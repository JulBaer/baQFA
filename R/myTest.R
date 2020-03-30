#setwd to the folder where the .dat files are

#devtools::load_all()
#devtools::document()
#
#setwd("C:/Users/julia/OneDrive/Dokumente/QFA2wip")
#
#setwd("C:/Users/julia/OneDrive/Dokumente/QFA_data/Colonyzer_modified_new/Images_normal_color/Output_Data")
#setwd("C:/Users/julia/Google Drive/QFA Images/2019-04-30-Saureus_Staphman_CI1995+2002+2007+2012")


#control=try(colonyzer.read(experiment="SAU1ExptDescription.txt", ORF2gene="ORF2GENE.txt",libraries="LibraryDescription.txt" , screenID="SAUtest1"))

#Strip non-experimental edge cultures
#control=control[(control$Row!=1)&(control$Col!=1)&(control$Row!=8)&(control$Col!=12),]


# Define which measure of cell density to use
#control$Growth=control$Intensity

#QFA run
#control.fit = NA
#control.fit<-qfa.fit2(control,inocguess=NULL,detectThresh=0,globalOpt=F,AUCLim=NA,STP=2,nrate=T,TimeFormat="h", Model="Glog",lowK=NA,upK=NA,lowr=NA,upr=NA,lowg=NA,upg=NA,lowv=NA,upv=NA,lowb=NA,upb=NA)

#control.fit$r[1] = NA
#
#control.fit=makeFitness2(control.fit,AUCLim=NA, plotFitness="MDR", filename=NA)
#
#ts = normalisePlates(control.fit, "b", "ORF")

#qfa.plot2("Gmp_GC2D2.pdf",control.fit,control)
#
#pl1 <- qfa.AgrregateCurves(control.fit, control)
#pl1




#setwd("C:/Users/julia/OneDrive/Dokumente/QFA_data/PlateReaderImportTest")
#
#GC24_1564 <- read.table("1564MIC_test_24h_GC.csv",sep=";",header=T)
#mta24_1564 <- read.table("1564MIC_test_24h_meta.csv",sep=";",header=T)
#control <- GC2QFA(GC24_1564,mta24_1564)
#
#control$Growth <- control$OD_calib
#
#control=control[control$Selected==T,]
#
#mta=mta24_1564

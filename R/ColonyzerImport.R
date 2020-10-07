#### Read in Colonyzer output and return a rod.read-like object ####
# This function essentially replaces the database functionalities of ROD and is therefore quite slow.
# I recommend that you run it once, and write the returned table to file for future use.
# It is *extremely* fast and flexible compared to using ROD itself however...
####


#' Read raw cell density timecourse data from Colonyzer output files
#'
#' Reads in and binds together all of the Colonyzer output files in a
#' directory and puts together a data.frame suitable for qfa.fit2 input.
#' Colonyzer is an open source image analysis tool for quantifying photographic
#' images by calculating cell densities of bacterial or eukaryotic spots and
#' colonies on agar plates. We recommend using the easy-to install and
#' easy to use version called BaColonyzer found on https://github.com/judithbergada/bacolonyzer .
#' The original version and some documentation about it can be found on http://research.ncl.ac.uk/colonyzer/
#' Required meta-data file processing is reduced compared to the read.colonyzer function of the
#' original QFAr package.
#'
#' The .dat files should contain the following parameters (they are automatically there if you used BaColonyzer or Colonyzer):
#'  \itemize{
#'
#' \item Image.Name - Full name at image capture (includes barcode and
#' date-time) of image from which data are derived \item Row - Row number
#' (counting from top of image) of culture in rectangular gridded array \item
#' Col - Column number (counting from left of image) of culture in rectangular
#' gridded array \item X.Offset - x-coordinate of top left corner of
#' rectangular tile bounding culture (number of pixels from left of image)
#' \item Y.Offset - y-coordinate of top left corner of rectangular tile
#' bounding culture (number of pixels from top of image) \item Area - Culture
#' area (pixels) \item Trimmed - Integrated Optical Density, sum of pixel
#' intensities within culture area \item Threshold - Global pixel intensity
#' threshold used for image segmentation (after lighting correction) \item
#' Intensity - Total pixel intensity for square tile containing culture \item
#' Edge Pixels - Number of pixels classified as culture on edge of square tile
#' \item Colony.Color.R - Culture red channel intensity \item Colony.Color.G -
#' Culture green channel intensity \item Colony.Color.B - Culture blue channel
#' intensity \item Background.Color.R - Background red channel intensity (for
#' current tile) \item Background.Color.G - Background green channel intensity
#' (for current tile) \item Background.Color.B - Background blue channel
#' intensity (for current tile) \item Edge.length - Number of culture pixels
#' classified as being microcolony edge pixels (useful for classifying
#' contaminants in cultures grown from dilute inoculum) \item Tile.Dimensions.X
#' - Culture tile width (pixels) \item Tile.Dimensions.Y - Culture tile height
#' (pixels) \item Growth - Default measure of cell density (direct copy of one
#' of Trimmed, Threshold or Intensity) \item Barcode - Unique plate identifier
#' \item Date.Time - Timestamp of image capture (extracted from image filename)
#' \item Inoc.Time - User specified date and time of inoculation (specified in
#' ExptDescription.txt file) \item Treatments - Conditions applied externally
#' to plates (e.g. temperature(s) at which cultures were grown, UV irradiation
#' applied, etc.) \item Medium - Nutrients/drugs in plate agar \item
#' Screen.Name - Name of screen (identifies biological repeats, and experiment)
#' \item RepQuad - Integer identifying which of the quadrants of a 1536 plate
#' were used to inoculate the current 384 plate (set equal to 1 for all
#' cultures for 1536 format for example) \item MasterPlate Number - Library
#' plate identifier \item Timeseries order - Sequential photograph number \item
#' Library.Name - Name of library, specifying particular culture location \item
#' ORF - Systematic, unique identifier for genotype in this position in arrayed
#' library \item Gene - Standard gene name for genotype in this position in
#' arrayed library.  Note that this can be set equal to ORF for example \item
#' ScreenID - Unique identifier for this QFA screen \item Client - Client for
#' whom screen was carried out \item ExptDate - A representative/approximate
#' date for the experiment (note that genome-wide QFA screens typically take
#' weeks to complete) \item User - Person who actually carried out screen \item
#' PI - Principal investigator leading project that screen is part of \item
#' Condition - The most important defining characteristic of screen, as
#' specified by user (e.g. the temperature screen was carried out at if screen
#' is part of multi-temperature set of screens, or the query mutation if part
#' of a set of screens comparing query mutations, or the drugs present in the
#' medium if part of a set of drug screens) \item Inoc - Qualitative identifier
#' of inoculation type (e.g. "DIL" for dilute inoculum, "CONC" for
#' concentrated).  Used to distinguish between experiments carried out with
#' different methods of inoculation. \item Expt.Time - Time (days) since
#' user-specified inoculation date (Inoc.Time) that current image was captured
#'
#' }
#'
#'
#' @param path String. The path to the folder containing the Colonyzer .dat files to be
#' read. Set to working directory by default.
#' @param files String vector. (Optional). Vector giving locations of Colonyzer .dat files to be
#' read (overrides path specified).
#' @param experiment String. (Optional). Name of text file describing the inoculation times,
#' library and plate number for unique plates. If this file is not specified,
#' the variables are taken if possible from the data in the .dat files and if not possible an arbitraryy "UNKNOWN"
#' value is given (see information in brackets below). Filename is taken relative to path if path is.
#' specified. File must be a tab-delimited text file with no header containing
#' the following columns: \itemize{
#' \item Barcode - unique identifier for each plate (no file specified: Barcodes from libraries file)
#' \item Start.time - Time of inoculation of the plate in format YYYY-MM-DD_hh_mm_ss (no file specified: first image date)
#' \item Treatment - Whatever treatment you applied to the plate (no file specified: UNKNOWN)
#' \item Medium - Whatver medium you used (no file specified: UNKNOWN)
#' \item Screen - Screening ID (no file specified: use screenID)
#' \item Library - Library ID (no file specified: Library values from libraries file)
#' \item Plate - ID replicates (no file specified: Plate values from libraries file)
#' \item RepQuad - Which quadrant was used for inoculation (no file specified: UNKNOWN)
#' }
#' @param ORF2gene String. (Optional). Filename of a tab-delimited text file containing two columns
#' (with no headers) associating unique, systematic strain identifiers (e.g.
#' yeast ORF Y-numbers) with human readable gene names (e.g. standard names
#' from SGD). If not specified, Gene and ORF are set to the same value
#' @param libraries String. The only necesarry meta-data in the new QFA iteration.
#' Tab-delimited text file describing each well as row-column coordinate of each plate in a series of
#' rectangular arrayed libraries. Header row format is: "Library ORF Plate Row
#' Column Notes".  Columns are: \itemize{
#' \item Library - Library identifier
#' (e.g. SAU1)
#' \item ORF - Systematic strain identifier (e.g. Cowan)
#' \item Plate - Plate number
#' \item Row - Row number
#' \item Column - Column number
#' \item Notes - Optional strain notes
#' }
#' @param screenID String. Unique experiment identifier (not the same as the Barcode which is a unique plate identifier).
#' For example, we use it as an identifier of the grid layout for competition assays (Grid vs.
#' Single strain vs Block).
#' @param Growth String. (Optional). Rowname of the .out datafiles to use as
#' the target Growth parameter for following analysis. Set to Intensity by default.
#'
#' @return An R data.frame where each row corresponds to a single observation
#' on a single colony, with the value of the growth measurement in 'Growth',
#' and the date and time of the measurement in 'Date.Time'. Other information
#' about the observation is stored in the other columns.  Several columns
#' returned are direct copies of Colonyzer output and mapped as follows:
#' \itemize{ \item Image.Name - Image Name \item Row - Spot Row \item Col -
#' Spot Column \item X.Offset - X Offset \item Y.Offset - Y Offset \item Area -
#' Area \item Trimmed - Trimmed Area \item Threshold - Threshold \item
#' Intensity - Intensity \item Edge.Pixels - Edge Pixels \item Colony.Color.R -
#' Colony Color R \item Colony.Color.G - Colony Color G \item Colony.Color.B -
#' Colony Color B \item Background.Color.R - Background Color R \item
#' Background.Color.G - Background Color G \item Background.Color.B -
#' Background Color B \item Edge.length - Edge length \item Tile.Dimensions.X -
#' Tile Dimensions X \item Tile.Dimensions.Y - Tile Dimensions Y }
#'
#' Extra columns are automatically added as follows.  Some of this information
#' is derived from auxiliary files passed to the function such as the
#' experimental description file, the orf-gene dictionary and the library
#' description file: \itemize{ \item Growth - A cell density surrogate built
#' from trimmed Area normalised by tile area and maximum achievable pixel
#' intensity: Trimmed/(Tile.Dimensions.X*Tile.Dimensions.Y*255) \item Barcode -
#' Plate identifier, essentially image name with date time and file extension
#' stripped \item Date.Time - Date time of image capture in YYYY-MM-DD_hh-mm-ss
#' format \item Inoc.Time - Date time that plate was inoculated.  If plate is
#' grown at a high temperature, date time at which plate was moved into high
#' temperature incubator.  The assumption in this case being that negligible
#' growth occurred before plate temperature was shifted the the target
#' temperature. \item Treatments - Treatments applied to plate (e.g.
#' temperature) \item Medium - Medium contained in agar (e.g. nutrients or
#' drugs added to agar) \item Screen.Name - Unique identifier for experiment
#' (usually identifies repeat number also if multiple repeats carried out).
#' \item RepQuad - Identifier for experiments scaling down from 1536 format
#' plates to 384, indicating which quadrant on the original 1536 source plate
#' the current 384 format plate belongs to. \item MasterPlate.Number -
#' Identifies which plate in the source library (as described in the library
#' description file) corresponds to the current plate \item Timeseries.order -
#' Ordinal describing which photograph captured \item Library.Name - Identifies
#' which of the libraries identified in the library description file was used
#' to construct this plate \item ORF - Unique systematic identifier for the
#' genotype of the strain at this location (e.g. yeast Y-number), as defined by
#' library description file \item Gene - Standard, human readable genotype
#' identifier for the strain at this location, as defined by the ORF-Gene
#' dictionary \item Background - Tag identifying experiment, typically used to
#' construct file names and axes titles in plots \item Expt.Time - Number of
#' days passed between inoculation (start of experiment) and current time }
#'
#' Finally, as well as returning the object above, this function prints a small
#' report to screen, summarising the data returned.  This includes number of
#' unique barcodes read, number of photos read, number of genotypes in
#' experiment, number of unique culture observations made, a list of treatments
#' applied, a list of media used, a list of unique screen names (e.g.
#' replicates carried out), the plate dimensions (e.g. 1536, 384 or 96 format)
#' and a list of unique inoculation dates.
#'
#' @examples
#' #qfa.testdata was generated with the call in the Not run section
#' #(files not included in package)
#' data(qfa.testdata)
#' #Strip non-experimental edge cultures
#' qfa.testdata = qfa.testdata[(qfa.testdata$Row!=1) & (qfa.testdata$Col!=1) & (qfa.testdata$Row!=8) & (qfa.testdata$Col!=12),]
#' # Define which measure of cell density to use
#' qfa.testdata$Growth = qfa.testdata$Intensity
#' GmpFit = qfa.fit2(qfa.testdata, inocguess=NULL, detectThresh=0, globalOpt=F, AUCLim=NA, TimeFormat="h", Model="Gmp")
#' # Construct fitness measures
#' GmpFit = makeFitness2(GmpFit, AUCLim=NA, plotFitness="All", filename="Example_Gmp_fitness.pdf")
#' # Create plot
#' qfa.plot2("Example_Gmp_GrowthCurves.pdf", GmpFit, qfa.testdata, maxt=30)
#'
#' \dontrun{
#' qfa.testdata = try(colonyzer.read(experiment = "SAU1ExptDescription.txt", ORF2gene = "ORF2GENE.txt",
#' libraries = "LibraryDescription1.txt" , screenID = "SAUtest1"))
#' }
#' @keywords qfa
colonyzer.read<-function(path=".", libraries="LibraryDescriptions.txt", Growth="Intensity",files=c(),experiment=NA,ORF2gene=NA,screenID=""){
  # Are we reading all files in a folder?
  if (length(strsplit(path,".")[[1]])>1){pathT=TRUE}else{pathT=FALSE}
  # If no Colonyzer files specified, use all .dat in working directory
  if (length(files)==0){
    fs<-list.files(path=path,pattern="\\.out$")
    # For historical datasets, there are no .out files, only headerless .dat files
    if(length(fs)==0) fs<-list.files(path=path,pattern="\\.dat$")
  } else {fs<-files}
  # If we're using a path, then add path to filenames
  if (pathT==TRUE) fs<-paste(path,fs,sep="/")
  print("List of data files to read:")
  print(fs)

  # Open the library descriptions
  libs=read.delim(libraries,sep="\t",header=TRUE,stringsAsFactors=FALSE)
  libs$ORF=toupper(libs$ORF)

  # What library are we using? e.g. library="AndrewsOE_v2" library="SDL_v2"
  liblist=unique(libs$Library)
  NLIB=length(liblist)
  libdict=1:NLIB
  names(libdict)=liblist
  #libs=libs[libs$Library==library,]
  # Make an array for storing info about spots
  NROW=max(libs$Row,na.rm=TRUE)
  NCOL=max(libs$Column,na.rm=TRUE)
  NPLATE=max(libs$Plate,na.rm=TRUE)

  SPOTARRAY=array("missing",dim=c(NLIB,NPLATE,NROW,NCOL))
  # Fill the spot array object
  for (x in 1:length(libs$ORF)) SPOTARRAY[libdict[[libs[x,"Library"]]],libs[x,"Plate"],libs[x,"Row"],libs[x,"Column"]]=libs[x,"ORF"]
  getORF<-function(lib,plate,row,col) SPOTARRAY[libdict[[lib]],plate,row,col]




  # Read in the image analysis output
  # Have to define colClasses here since from R3.1-10, automatic conversion for string representation of numbers
  if(grepl("\\.out",fs[1])){hasHeader=TRUE}else{hasHeader=FALSE}
  if("data.table" %in% rownames(installed.packages())){
    library(data.table) # Faster version of read.delim...
    iman=do.call(rbind, lapply(fs, data.table::fread,header=hasHeader,sep="\t",stringsAsFactors=FALSE))
    iman=data.frame(iman)
  }else{
    iman=do.call(rbind, lapply(fs, read.delim,header=hasHeader,sep="\t",stringsAsFactors=FALSE))
  }
  if(!hasHeader){
    colnames(iman)=c("Filename","Row","Column","X.Offset","Y.Offset","Area","Trimmed","Threshold","Intensity","Edge.Pixels","redMean","greenMean","blueMean","redMeanBack","greenMeanBack","blueMeanBack","Edge.Length","Tile.Dimensions.X","Tile.Dimensions.Y")
    iman$x=iman$X.Offset+iman$Tile.Dimensions.X/2.0
    iman$y=iman$Y.Offset+iman$Tile.Dimensions.Y/2.0
    iman$Diameter=(iman$Tile.Dimensions.X+iman$Tile.Dimensions.Y)/2.0
    if(max(iman$Intensity)>1.0) iman$Intensity=iman$Intensity/(iman$Tile.Dimensions.X*iman$Tile.Dimensions.Y*255)
    if(max(iman$Trimmed)>1.0) iman$Trimmed=iman$Trimmed/(iman$Tile.Dimensions.X*iman$Tile.Dimensions.Y*255)
  }

  # Create extra columns
  if(nchar(iman$Filename[1])==31){
    if(!"Barcode"%in%colnames(iman)) iman$Barcode=substr(iman$Filename,1,11)
    iman$Date.Time=substr(iman$Filename,13,31)
  }else{
    if(!"Barcode"%in%colnames(iman)) iman$Barcode=substring(iman$Filename, 1, regexpr("_", iman$Filename) -1)
    iman$Date.Time=substr(iman$Filename,nchar(iman$Barcode)+2,nchar(iman$Filename))
  }


  # Open the experimental description
  if (!is.na(experiment)){
    expt=read.delim(experiment,sep="\t",header=TRUE,stringsAsFactors=FALSE)
  }else{
    bcsi=unique(iman$Barcode)
    expt=data.frame(Barcode=bcsi)
    for (i in 1:length(bcsi)){
      expt$Start.Time[expt$Barcode==bcsi[i]]=min(iman$Date.Time[iman$Barcode==bcsi[i]])
    }

    expt$Treatment="UNKNOWN"
    expt$Medium="UNKNOWN"
    expt$Screen=screenID
    expt$Library=unique(libs$Library)
    expt$Plate=unique(libs$Plate)
    expt$RepQuad="UNKNOWN"
  }




  # Watch the mismatch between library names (e.g. BooneSDLV2 and SDL_v2)
  # barcode-> plate dictionaries for all treatments and repeats
  barcStart=expt$Start.Time
  names(barcStart)=expt$Barcode

  barcTreat=expt$Treatment
  names(barcTreat)=expt$Barcode

  barcMed=expt$Medium
  names(barcMed)=expt$Barcode

  barcScreen=expt$Screen
  names(barcScreen)=expt$Barcode

  if(!"RepQuad"%in%colnames(expt)) stop("Woah there!  We need a RepQuad column in expt. description file")
  barcQuad=expt$RepQuad
  names(barcQuad)=expt$Barcode

  barcPlate=expt$Plate
  names(barcPlate)=expt$Barcode

  barcLib=expt$Library
  names(barcLib)=expt$Barcode

  if("Client"%in%colnames(expt)){
    barcClient=expt$Client
    names(barcClient)=expt$Barcode
  }

  if("ExptDate"%in%colnames(expt)){
    barcExptDate=expt$ExptDate
    names(barcExptDate)=expt$Barcode
  }

  if("User"%in%colnames(expt)){
    barcUser=expt$User
    names(barcUser)=expt$Barcode
  }

  if("PI"%in%colnames(expt)){
    barcPI=expt$PI
    names(barcPI)=expt$Barcode
  }

  if("Condition"%in%colnames(expt)){
    barcCond=expt$Condition
    names(barcCond)=expt$Barcode
  }

  if("Inoc"%in%colnames(expt)){
    barcInoc=expt$Inoc
    names(barcInoc)=expt$Barcode
  }

  if (!is.na(ORF2gene)){
    # Open the ORF2GENE file
    orf2gene=read.delim(ORF2gene,sep="\t",header=FALSE,stringsAsFactors=FALSE)
    colnames(orf2gene)=c("orf","gene")
    orf2gene$orf=toupper(orf2gene$orf)
    # Add a "missing" row
    orf2gene=rbind(orf2gene,c("missing","missing"))
    orf2gene=rbind(orf2gene,c("MISSING","MISSING"))
    # Create an ORF2Gene dictionary
    getGene=orf2gene$gene
    names(getGene)=orf2gene$orf
  }else{

  }



  # Dump any images which are not in the experimental description file
  iman=iman[iman$Barcode%in%expt$Barcode,]
  smalliman=iman[(iman$Row==1)&(iman$Column==1),]

  # Make sure that there are not multiple copies of files present
  counts=c(table(smalliman$Filename))
  probs=counts[counts>1]
  if(length(probs)>0){
    print("ERROR: The following images were found in multiple locations")
    print(names(probs))
    print("You should probably look up their locations in the relevant folder, decide which location is correct and delete the files at other locations before reading data in again")
    return(NULL)
  }

  fmt="%Y-%m-%d_%H-%M-%S"

  # Create a dictionary for filename->photo number
  getPhotoNum<-function(filename,fmt="%Y-%m-%d_%H-%M-%S"){
    # Get plate name from filename
    #nlst=strsplit(filename,"_")[[1]]
    #lenlst=length(nlst)
    #platename=paste(nlst[1:(lenlst-2)],collapse="_")
    lenFname=nchar(filename)
    lenDate=nchar(format(Sys.time(), fmt))
    platename=substring(filename,1,lenFname-lenDate-1)
    # Filter iman data frame by filename
    #tmp=na.omit(smalliman[(smalliman$Barcode==platename),])
    tmp=smalliman[(smalliman$Barcode==platename),]
    tmp=tmp[order(tmp$Filename),]
    tmp$PhotoNum=1:length(tmp$Filename)
    return(as.numeric(tmp$PhotoNum[tmp$Filename==filename]))
  }

  fnames=unique(iman$Filename)
  photoNum=sapply(fnames,getPhotoNum,fmt)
  names(photoNum)=fnames

  iman$Inoc.Time=barcStart[iman$Barcode]
  iman$Treatments=barcTreat[iman$Barcode]
  iman$Medium=barcMed[iman$Barcode]
  iman$Screen.Name=barcScreen[iman$Barcode]
  iman$RepQuad=barcQuad[iman$Barcode]
  iman$MasterPlate.Number=barcPlate[iman$Barcode]
  iman$Timeseries.order=as.numeric(photoNum[iman$Filename])
  iman$Library.Name=barcLib[iman$Barcode]
  iman$ORF=mapply(getORF, iman$Library.Name, iman$MasterPlate.Number, iman$Row, iman$Column)
  if (!is.na(ORF2gene)){
    iman$Gene=getGene[toupper(iman$ORF)]
  }else{
    iman$Gene=toupper(iman$ORF)
  }

  iman$ScreenID=rep(screenID,length(iman$Filename))
  if("Client"%in%colnames(expt)){iman$Client=barcClient[iman$Barcode]}
  if("ExptDate"%in%colnames(expt)){iman$ExptDate=barcExptDate[iman$Barcode]}
  if("User"%in%colnames(expt)){iman$User=barcUser[iman$Barcode]}
  if("PI"%in%colnames(expt)){iman$PI=barcPI[iman$Barcode]}
  if("Condition"%in%colnames(expt)){iman$Condition=barcCond[iman$Barcode]}
  if("Inoc"%in%colnames(expt)){iman$Inoc=barcInoc[iman$Barcode]}

  t0<-as.POSIXlt(as.character(iman$Inoc.Time),format=fmt)
  t1<-as.POSIXlt(as.character(iman$Date.Time),format=fmt)
  iman$Expt.Time=as.numeric(difftime(t1,t0,units="days"))

  # Print checks so people are sure they're using the correct data #
  print(paste("Number of Barcodes :",length(unique(iman$Barcode))))
  print(paste("Number of Plate Photos :",length(unique(iman$Filename))))
  print(paste("Number of ORFs :",length(unique(iman$ORF))))
  print(paste("Number of Culture Images :",length(iman$Filename)))
  print("Treatments :")
  print(unique(iman$Treatments))
  print("Media:")
  print(unique(iman$Medium))
  print("Screens:")
  print(unique(iman$Screen.Name))
  if("Client"%in%colnames(iman)){print("Client :"); print(unique(iman$Client))}
  if("ExptDate"%in%colnames(iman)){print("Experiment date :"); print(unique(iman$ExptDate))}
  if("User"%in%colnames(iman)){print("User :"); print(unique(iman$User))}
  if("PI"%in%colnames(iman)){print("PI :"); print(unique(iman$PI))}
  if("Condition"%in%colnames(iman)){print("Condition :"); print(unique(iman$Condition))}
  if("Inoc"%in%colnames(iman)){print("Inoculation type :"); print(unique(iman$Inoc))}
  platesize<-max(as.numeric(iman$Row))*max(as.numeric(iman$Column))
  if (length(iman$Date.Time)%%platesize!=0){
    warning("Number of cultures not multiple of plate size")}
  print("Inoculation DateTimes:")
  print(unique(iman$Inoc.Time))

  if (Growth %in% names(iman)){
    grw=iman[names(iman)==Growth]
    iman$Growth=grw[,]
  }else{
    warning("The specified Growth variable is no valid rowname. Please specify Growth manually")
  }



  return(iman)
}


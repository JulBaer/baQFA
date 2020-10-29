##### This is for the qfa.plot2 pdf generation function----

#' Creating plots of growth curves including models for QFA dataframe in a pdf
#'
#' Produces a multipage pdf of growth curves including model fits if present.
#' Each page corresponds to a
#' single plate and growth curves are arrayed on the page according to their
#' position on the plate. Both observations (red crosses) and the fitted
#' growth curve from \code{\link{qfa.fit2}} (solid black curves) are shown for each culture.
#' The time at which the maximum slope of the observed growth curve on the log
#' scale occurs (solid blue line) and the time at which the maximum slope of
#' the raw observed data occurs (dashed blue line) are also indicated. Where
#' available, various fitness estimates are displayed for each culture,
#' together with culture genotype. These pdfs are useful for visually
#' checking quality of model fit & data.
#'
#'
#' @param file String. The pdf name to store the plots
#' @param results data.frame. The output of qfa.fit2 which contains the model fit data. If curves = F, results
#' should be the same as d for plotting of data without model fits
#' @param d data.frame. The original data.frame fed to qfa.fit2 containing all of the
#' timecourse data. Output of either Colonyzer.Read or GC2QFA
#' @param fmt String. The format in which Exp.Time of measurement and inoculation time
#' are stored. Default set to standard Colonyzer.Read or GC2QFA output
#' @param barcodes Vector of strings. Plot only for the plates with barcodes in this character
#' vector; all by default
#' @param master.plates Vector of strings. Plot only for the plates from master.plates in this
#' character vector; all by default
#' @param treatments Vector of strings. Plot only for the plates with treatments in this character
#' vector; all by default
#' @param screen.names Vector of strings. Plot only for the plates with screen.names in this
#' character vector; all by default
#' @param screenIDs Vector of strings. Plot only for the plates with screenIDs in this character
#' vector; all by default
#' @param maxg Numeric. Upper cell density (y-axis limit) for all growth curve plots.
#' Default value is NA.  If maxg=NA, the y-limit is chosen automatically to show all datapoints
#' @param maxt Numeric. Growth curve is plotted from time t = 0 to maxt. If set to NA,
#' maxt is set from the data as the latest timepoint found in the data.frame d
#' @param logify logical. Indicates whether growth curve plots should be on a
#' semilog scale. Not recommended.
#' @param densityCol String. Name of column in data frame d which contains cell density
#' estimate.  Note that the image analysis software Colonyzer provides several
#' possible alternatives. Default set to Growth which needs to be manually specified
#' before calling qfa.fit2
#' @param curves logical. Indicates whether fitted model curves should be
#' drawn. Useful to set this to false when generating diagnostic, data-only
#' growth curves. To do this, set curves = FALSE and pass the same (data)
#' object twice, as both the d and results arguments.  For example:
#' qfa.plot2("test.pdf",df,df,curves=FALSE)
#' @param ylabel String. y-axis label on growth curve plots.
#' @param ptype String. Plot type for data: "p": points, "l": lines, "b": both. Default: "p"
#' @param ming Numeric. Lower y-axis limit. If set to NA, minimal Growth value is used to set ming. Default set to NA
#' @param TimeFormat String. Either "h" for hours or "d" for days. This defines the
#' format of the x-axis. Although, it only needs to be specified if curves=F
#' as it is read otherwise from the results data.frame
#' @param mainTitle String. A string for the main title of the pdf plot. If set to NA
#' (default), a combination of Barcode, Treatment, Medium and Master.Plate columns
#' of the results data.frame will be used for the main title
#'
#' @return This function creates a pdf file stored in the working directory
#' named acording to the user input. On the pdf file, growth curves of the observation
#' data.frame d and, if curves=T and the specified results data.frame contains
#' results from the qfa.fit2 function, curves for the model fits are displayed
#' per position.
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
qfa.plot2<-function(file="QFA_GrowthCurves.pdf", results, d, curves=TRUE, maxg=NA, maxt=NA, ming=NA,mainTitle=NA, TimeFormat="d",
                    ylabel="Cell density (AU)",ptype="p", fmt="%Y-%m-%d_%H-%M-%S", logify=FALSE, densityCol="Growth",
                    barcodes=c(), master.plates=c(), treatments=c(), screen.names=c(), screenIDs=c()){
  # creates a pdf plot of the input data showing the specified growth measurements, fitted model and estimated parameters.

  #if results equals d, no curves need to be plotted, only the data:
  if (sum(dim(results)==dim(d))==2 && sum(results!=d, na.rm=T)==0) curves=F
  if (curves) {
    #check which TimeFormat was used, adjust max time for plotting if not set by user manually
    TimeFormat=unique(results$TimeFormat)
    #check which model was used in qfa.fit2
    Model <- unique(results$Model)
    if (Model=="Gmp") {Gmp=T;fixG=F; print(paste("Gompertz model plot, TimeFormat set to ", TimeFormat, sep=""))
    }else if (Model=="Glog") {Gmp=F;glog=T; print(paste("Generalized logistic model plot, TimeFormat set to ", TimeFormat, sep=""))
    }else if (Model=="Slog") {Gmp=F;glog=F; print(paste("Standard logistic model plot, TimeFormat set to ", TimeFormat, sep=""))
    }else {stop("Your dataset was generated with an older version. Please re-run qfa.fit2 with a Model specified")}


  }else{
    Model=NA
    print("Data only plotting")
  }


  #define x-axis maxima
if (is.na(maxt)) {
  if (TimeFormat=="h") {
    maxt=ceiling(max(d$Expt.Time, na.rm = T))*24
  }else{
    maxt=ceiling(max(d$Expt.Time, na.rm = T))
  }
}

  #turn off all existing pdf devices
  while (!is.null(dev.list())) dev.off()


  #remnant of old coding
  if(!"Column"%in%colnames(d)) d$Column=d$Col
  if(!"Col"%in%colnames(results)) results$Col=results$Column

  # Sort the data to be plotted sensibly, allowing easy comparison between repeats
  results=results[order(results$MasterPlate.Number,results$Treatment,results$Screen.Name),]
  # Get character vectors of requested barcodes, treatements,etc.; all if none specified
  #otherwise subset for these inputs
  if (!exists("master.plates")) master.plates=vector()
  if (length(master.plates)==0){master.plates<-unique(results$MasterPlate.Number)}
  results<-results[results$MasterPlate.Number%in%master.plates,]
  d<-d[d$MasterPlate.Number%in%master.plates,]
  if (length(treatments)==0){treatments<-unique(results$Treatment)}
  results<-results[results$Treatment%in%treatments,]
  d<-d[d$Treatments%in%treatments,]
  if (!exists("screen.names")) screen.names=vector()
  if (length(screen.names)==0){screen.names<-unique(results$Screen.Name)}
  results<-results[results$Screen.Name%in%screen.names,]
  d<-d[d$Screen.Name%in%screen.names,]
  if (!exists("screenIDs")) screenIDs=vector()
  if (length(screenIDs)==0){screenIDs<-unique(results$ScreenID)}
  results<-results[results$ScreenID%in%screenIDs,]
  d<-d[d$ScreenID%in%screenIDs,]
  if (length(barcodes)==0){barcodes<-sort(unique(results$Barcode))}
  results<-results[results$Barcode%in%barcodes,]
  d<-d[d$Barcode%in%barcodes,]

  if (curves){
    print(paste("Plotting for ",length(results[,1])," colonies", sep=""))
  }else{
    print(paste("Plotting for ",length(unique(paste(results$Col,results$Row,results$Barcode,sep="_")))," colonies", sep=""))
  }

  # Produce PDF. Modify some of the paramters depending on how many spots there are
  colmin<-min(results$Col); rowmin<-min(results$Col)
  colmax<-max(results$Col); rowmax<-max(results$Row)
  #add .pdf if not done by user:
  if (substr(file, nchar(file)-4+1, nchar(file))!=".pdf") {
    file = paste(file,".pdf", sep="")
  }
  #pdf size:
  pdf(file,4*(colmax-colmin+1),4*(rowmax-rowmin+1))
  cexfctr<-(rowmax*colmax)/384; marge=c(2,1,2,0.75)

  #if(rowmax*colmax==96) {cexfctr=1.0;marge=c(2,1,2,0.75)}
  if(rowmax*colmax>384) {cexfctr<-(rowmax*colmax)/1536; marge=c(2.0,1.5,2.0,1.5)}
  #cexfctr=cexfctr
  #cexfctr=1.0

  #go over all barcodes (identifiers of individual plates) for plotting:
  bcount<-0; nbc<-length(barcodes)
  a=as.logical(is.finite(unlist(d$Growth)))
  if (is.na(maxg)) {
    maxg=try(round(max(d$Growth[a], na.rm=T),3))
        if (is.na(maxg) | is.infinite(maxg)){
          maxg=1
          warning("maxg was not possible to determine programmatically, set to 1. Check your input data!")
        }
  }
  if (is.na(ming)) {
      ming=try(round(min(d$Growth[a], na.rm=T),3))
      if (is.na(ming) | is.infinite(ming)){
        ming=0
        warning("ming was not possible to determine programmatically, set to 0. Check your input data!")
      }
      if (logify & ming <= 0) ming = 0.00001
  }


  for (bcode in barcodes){
    bcode<-as.character(bcode)
    # Say what you're doing
    bcount<-bcount+1
    print(paste("Plotting for Plate ",bcount,"/",nbc,": ",bcode,". Please wait...", sep=""))
    # Restricts results and data to this plate
    rbc<-results[as.character(results$Barcode)==bcode,]
    dbc<-d[as.character(d$Barcode)==bcode,]
    # Get starting time for this plate
    inoctime<-rbc$Inoc.Time[1]
    #as.POSIXlt(as.character(bigd$Date.Time),format=fmt)
    # Find number of rows and columns on this plate
    #nrow<-max(rbc$Row)-min(rbc$Row)+1; ncol<-max(rbc$Column)-min(rbc$Column)+1
    nrow<-rowmax-rowmin+1; ncol=colmax-colmin+1
    #nrow=length(unique(rbc$Row)); ncol=length(unique(rbc$Col))

    # Find max growth to set y-axis for this plate if not specified by user

    # Set graphics parameters for each plate
    op<-par(mfrow=c(nrow,ncol),oma=c(13,15,22,1),
            mar=marge,mgp=c(3,1,0),cex=cexfctr*2.5)
    ## Plot for each row of results for that bcode ##
    for(rno in rowmin:rowmax){
      for(cno in colmin:colmax){
        params=rbc[(rbc$Row==rno)&(rbc$Col==cno),]
        if(nrow(params)>0){
          #actual call for the plot:
          rowplot2(params,dbc,inoctime,maxg,fmt,maxt,logify,densityCol=densityCol,curves=curves,ptype=ptype,ming=ming,TimeFormat=TimeFormat,Gmp=Gmp)
        }else{
          frame()
        }
      }
    }
    #z<-apply(rbc,1,rowplot,dbc,inoctime,maxg,fmt,maxt,logify,densityCol=densityCol,curves=curves,ptype=ptype)
    # Title for the plate. if not user specified, use this default combination of barcode, treatment, medium and plate nr
    if (is.na(mainTitle)){
      maintit<-paste(rbc$Barcode[1],"Treatment:",rbc$Treatment[1],
                     "Medium:",rbc$Medium[1],"Plate:",rbc$MasterPlate.Number[1],sep=" ")
    }else{
      if (length(barcodes>1)) {
        maintit=paste(as.character(mainTitle), rbc$Barcode[1], sep=", ")
      }else{
        maintit=as.character(mainTitle)
      }

    }

    #this is the x-axis label (in correct timeformat) at the bottom of the page
    cextit<-500/nchar(maintit)
    if (TimeFormat=="h") {
      title(main=maintit,xlab="Time since inoculation (hours)",line=7,
            ylab=ylabel,cex.main=cextit,cex.lab=8,outer=TRUE)
    }else{
      title(main=maintit,xlab="Time since inoculation (days)",line=7,
            ylab=ylabel,cex.main=cextit,cex.lab=8,outer=TRUE)
    }

    print(paste("Plotting for Plate ",bcount,"/",nbc,": ",bcode," finished.", sep=""))
    par(op)} #bcode
  dev.off()
}

#### Plot a colony's timecourse from a row of the results
rowplot2<-function(resrow,dbc,inoctime,maxg,fmt,maxt,logify,densityCol="Growth",curves=TRUE,ptype="p",ming=0.00001, TimeFormat="d",Gmp=Gmp){
  #get which row and col we are in
  if (curves){
    row<-as.numeric(resrow['Row']); col<-as.numeric(resrow['Col'])
  }else{
    row<-as.numeric(unique(resrow['Row'])); col<-as.numeric(unique(resrow['Col']))
    }

  #remnant of old coding:
  if ('Gene'%in%names(resrow)){gene<-unique(resrow['Gene'])} else {gene<-unique(resrow['ORF'])}
  # Get data for that colony:
  dcol<-dbc[(dbc$Row==row)&(dbc$Column==col),]
  #nozero replaces all values below zero with a very low number close to zero (mainly for log plots useful)
  growth<-sapply(dcol[[densityCol]],nozero)
  if (is.list(growth)) {growth=unlist(growth)}
  #define timevector
  if (TimeFormat=="h") {
    tim<-dcol$Expt.Time*24
  } else {
    tim<-dcol$Expt.Time
  }

  #get the timeshift. needed for model curve plot of Glog and Slog
  if("tshift"%in%names(resrow)){tshift<-unique(as.numeric(resrow[["tshift"]]))[1]}else{tshift=0}
  # Draw the curves and data. Finally :)
  logdraw2(row,col,resrow,tim,growth,gene,maxg,maxt=maxt,logify=logify,densityCol=densityCol,curves=curves,ptype=ptype,tshift=tshift,logmin=ming,TimeFormat,Gmp=Gmp)
}

### Converts row no. to position vector ###
index2pos<-function(index,dbc) c(dbc[[index,'Row']],dbc[[index,'Column']])

### Do individual timecourse plot given parameters & data ###
logdraw2<-function(row,col,resrow,tim,growth,gene,maxg,fitfunct,maxt=0,scaleT=1.0,logify=FALSE,densityCol="Growth",curves=TRUE,ptype="p",tshift=0,logmin=0.00001,TimeFormat="d",Gmp=F){
  #define ylim:
  if(logify) {
    ylog="y"
    ylim=c(logmin,1.1*maxg)
  }else{
    ylog=""
    ylim=c(logmin,1.1*maxg)
  }
  #create empty plot with all the necessary infos
  plot(NULL,type="n",xlim=c(0,maxt),ylim=ylim,log=ylog,xlab="",ylab="",main=gene,frame.plot=0,cex.main=3*scaleT,cex.axis=1*scaleT)
  # Add data points in red
  points(tim,growth,col="red",cex=1.25*scaleT,pch=4,lwd=2,type=ptype,xlim=c(0,maxt),ylim=c(logmin,1.2*maxg))
  #and if the model plot is requested:
  if(curves){
    if (Gmp==T) {
      if("b"%in%names(resrow)) abline(v=resrow['b'],col="darkgreen")
    }
    
    #draw blue lines for model-free estimation of time to max growth from linear and log scale
    if("nr_t"%in%names(resrow)) abline(v=resrow['nr_t'],col="blue")
    if("maxslp_t"%in%names(resrow)) abline(v=resrow['maxslp_t'],col="blue",lty=2)
    

    # Get the following model-free parameters. if not specfied, replace with NA as placeholder
    if ("maxslp" %in% names(resrow)) {
      maxslp <- as.numeric(resrow['maxslp']);
      maxslp_t <- as.numeric(resrow['maxslp_t']);
      Logmaxslp <- as.numeric(resrow['nr']);
      Logmaxslp_t <- as.numeric(resrow['nr_t']);
    }else{
      maxslp <- NA
      maxslp_t <- NA
      Logmaxslp <- NA
      Logmaxslp_t <- NA
    }

    #same here. Get values of the model estimated parameters, if not present, take NA as placeholder
    if ("MDR" %in% names(resrow)) {
      MDR<-as.numeric(resrow['MDR']);
      MDP<-as.numeric(resrow['MDP']);
      AUC<-as.numeric(resrow['AUC']);
      DT<-as.numeric(resrow['DT']);
    }else{
      MDR<-NA
      MDP<-NA
      AUC<-NA
      DT<-NA
    }



    if(logify) {ylog="y"}else{ylog=""}
    # Add logistic curve
    if(maxt==0) maxt=ceiling(max(tim))
    x=0
    #get the corresponding model parameters for the different  models, plot the curve
    if (!Gmp) {
      K<-as.numeric(resrow['K']); r<-as.numeric(resrow['r']); g<-as.numeric(resrow['g']); v<-as.numeric(resrow['v']);
      curve(Glogist(K,r,g,v,x-tshift),n=31,lwd=2.5,add=TRUE,from=0,to=maxt,xlim=c(0,maxt),ylim=c(0.00001,1.2*maxg))
    }else{
      K<-as.numeric(resrow['K']); r<-as.numeric(resrow['r']); b<-as.numeric(resrow['b']); g<-as.numeric(resrow['g']);
      plShift=0
      curve((Gmprtz(K,r,b,x)+plShift),n=31,lwd=2.5,add=TRUE,from=0,to=maxt,xlim=c(0,maxt),ylim=c(0.00001,1.2*maxg))
    }



  }
  #add text with model parameters
  if(curves){
    # Add legend
    if (!Gmp) {
      legt1<-paste(c("K=","r=","g=","v=","MDR=","MDP=","AUC=","DT=","maxslp=","maxslp t=", "Lmaxslp=","Lmaxslp t="),
                   c(signif(K,3),signif(r,3),signif(g,3),signif(v,3),signif(MDR,3),signif(MDP,3),signif(AUC,3),signif(DT,3),
                     signif(maxslp,3),signif(maxslp_t,3),signif(Logmaxslp,3),signif(Logmaxslp_t,3)),sep="")
    }else{
      legt1<-paste(c("K=","r=","b=","g=","MDR=","MDP=","AUC=","DT=","maxslp=","maxslp t=", "Lmaxslp=","Lmaxslp t="),
                   c(signif(K,3),signif(r,3),signif(b,3),signif(g,3),signif(MDR,3),signif(MDP,3),signif(AUC,3),signif(DT,3),
                     signif(maxslp,3),signif(maxslp_t,3),signif(Logmaxslp,3),signif(Logmaxslp_t,3)),sep="")
    }

    if(logify){legend("bottomright",legt1,box.lty=0,cex=scaleT*1.5)}else{legend("topleft",legt1,box.lty=0,cex=scaleT*1.5)}
  }else{
    legt1 = " "
    if(logify){legend("bottomright",legt1,box.lty=0,cex=scaleT*1.5)}else{legend("topleft",legt1,box.lty=0,cex=scaleT*1.5)}
  }
  #legend("topright",sprintf("R%02dC%02d",row,col),box.lty=0,cex=scaleT)
}

### Get the mode of a vector (why isn't there such a function in base?)
getMode <- function(x) {
  ux = unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}





##### qfa.AggregateCurves----


#' Plot all data measurement points including model curves grouped by userinput in a ggplot object
#'
#' This functions creates a ggplot object with the input data.frames from either Read.Colonyzer
#' or GC2QFA functions (data.frame GC) and the model fits from a subsequent call to qfa.fit2
#' (data.frame results) of the same data.frame. The data can also be displayed without the model-fits.
#' For this, the user inputs the same data.frame for results and GC which was generated by either Read.Colonyzer
#' or GC2QFA. The qfa.AggregateCurves functions lets the user specify by which factors the plot should be facetted and colored.
#' These factors need to be present in both results and GC data.frame as a column.
#'
#' @param results data.frame. Usually, the output data.frame of qfa.fit2 to create model curves. It can also be the output data.frame of Read.Colonyzer
#' or GC2QFA if the data without model curves should be displayed
#' @param GC data.frame. The data.frame used as an input for qfa.fit2 to generate the results data.frame. If the results data.frame is the output
#' of either Read.Colonyzer or GC2QFA, then GC needs to be the same data.frame
#' @param plotCol String. The name of the column header which should be used for coloring of the curves
#' @param plotFacet String. The name of the column header which should be used for facetting of the plots
#' @param yname String. Label of the y-axis
#' @param legName String. Label for the color legend
#' @param PlotName String. Title for the plot
#' @param markb logical. Indicator if the b parameter of the Gompertz model should be marked as a
#' colored vertical line. Ignored if another model was used to create the results data.frame
#' @param TimeFormat String. Either "d" for days or "h" for hours. If set to NA and results contains model-fits, the TimeFormat
#' specified in the call to qfa.fit2 is used
#'
#' @return This function just plots the ggplot object if not used with an assignment. If it is
#' called as an assignment, the return is a ggplot object which can further be modified by
#' the user if desired.
#'
#' @examples
#' #qfa.testdata was generated with the call in the Not run section
#' #(these files are not included in the package)
#' data(qfa.testdata)
#' #Strip non-experimental edge cultures
#' qfa.testdata = qfa.testdata[(qfa.testdata$Row!=1) & (qfa.testdata$Col!=1) & (qfa.testdata$Row!=8) & (qfa.testdata$Col!=12),]
#' # Define which measure of cell density to use
#' qfa.testdata$Growth = qfa.testdata$Intensity
#' GmpFit = qfa.fit2(qfa.testdata, inocguess=NULL, detectThresh=0, globalOpt=F, AUCLim=NA, TimeFormat="h", Model="Gmp")
#' # Construct fitness measures
#' GmpFit = makeFitness2(GmpFit, AUCLim=NA, plotFitness="All", filename="Example_Gmp_fitness.pdf")
#' # Create aggregated curves plot
#' a = qfa.AggregateCurves(GmpFit, qfa.testdata, plotFacet="ORF", yname="Intesity [AU]", legName="Strain", markb=F)
#' # Demonstrating subsequent ggplot object manipulations:
#' a = a + ggplot2::ggtitle("Alternate title")
#' a
#'
#' @keywords qfa
#'
qfa.AggregateCurves <- function(results, GC, plotCol="ORF", plotFacet="ScreenID", yname="Growth",legName=NA, PlotName=NA, markb=F, TimeFormat=NA) {
  #plot the model for each position along data measurement points. Let user specify color and facetting variable

  #input variales:
  #results: output of qfa.fit2 and optionally makeFitness2
  #GC: the growth curves, either from GC2QFA or colonyzer.read with optional qfa.normalize
  #timeFormat: h or d for plotting. Read it out of the
  #plotCol: Factor column that is used for color in plot. need to be a column name of GC and results
  #plotFacet: Factor column that is used for facetting in plot. need to be a column name of GC and results
  #yname: y-axis label
  #legName: color legend title
  #PlotName: main title of plot


  qfa_check_packages()




  #get which model
  Model <- unique(results$Model)
  if (sum(dim(results)==dim(GC))==2 && sum(results!=GC, na.rm=T)==0) Model="none"
  if (Model=="Gmp") {Gmp=T;fixG=F; print("Gompertz model used"); mdl <- "Gompertz model"
  }else if (Model=="Glog") {Gmp=F;glog=T;markb=F; print("Generalized logistic model used"); mdl <- "Generalized logistic model"
  }else if (Model=="Slog") {Gmp=F;glog=F;markb=F; print("Standard logistic model used"); mdl <- "Standard logistic model"
  }else if (Model=="none") {
    print("Data without a model fit is going to be plotted"); mdl = "no model"; Gmp=F;glog=F;markb=F;

  }else {stop("Please make sure the input data.frame match and are created according to the description in the help-page")}

  #get timeformat, modify GC time column
  if (is.na(TimeFormat) & Model!="none") {
    TimeFormat=unique(results$TimeFormat)
  }else{
    stop("Please specify TimeFormat for plotting without model-fit")
  }

  if (TimeFormat=="h") {
      GC$Time=GC$Expt.Time*24
      xname="Time [h]"
    }else{
      GC$Time=GC$Expt.Time
      xname="Time [d]"
  }


  #paste row and col for unique identifier
  results$RowCol=as.factor(paste(results$Row,results$Column,results$Barcode, sep="_"))
  GC$RowCol=as.factor(paste(GC$Row,GC$Column,GC$Barcode, sep="_"))
  id <- unique(GC$RowCol)
  cnt=1

  if (length(plotFacet)==2) {
    plotFacet1 = plotFacet[1]
    plotFacet2 = plotFacet[2]
    if (plotFacet1==plotFacet2) {
      stop("Please do not specify twice the same factor name for facetting")
    }
    EqNames = (sum(plotCol!=plotFacet)==2)
  }else if (length(plotFacet)>2) {
    stop("Please specify at max 2 faceting variables")
  }else{
    EqNames = (plotCol==plotFacet)
  }


  #modify the colnames of results and GC file for easy ggplotting
  #names(results)[names(results)==plotX] <- "plotX"
  if (!EqNames) {
    names(GC)[names(GC)==plotCol] <- "plotCol"
    names(results)[names(results)==plotCol] <- "plotCol"
    if (length(plotFacet)==1){
      names(results)[names(results)==plotFacet] <- "plotFacet"
      names(GC)[names(GC)==plotFacet] <- "plotFacet"
    }else{
      names(results)[names(results)==plotFacet1] <- "plotFacet1"
      names(GC)[names(GC)==plotFacet1] <- "plotFacet1"
      names(results)[names(results)==plotFacet2] <- "plotFacet2"
      names(GC)[names(GC)==plotFacet2] <- "plotFacet2"
    }

  }else{
    GCbck = GC
    resultsbck = results
    names(GC)[names(GC)==plotCol] <- "plotCol"
    names(results)[names(results)==plotCol] <- "plotCol"

    if (length(plotFacet)==1){
    names(resultsbck)[names(resultsbck)==plotFacet] <- "plotFacet"
    names(GCbck)[names(GCbck)==plotFacet] <- "plotFacet"
    GC$plotFacet = GCbck$plotFacet
    results$plotFacet = resultsbck$plotFacet
    }else{
      names(resultsbck)[names(resultsbck)==plotFacet1] <- "plotFacet1"
      names(GCbck)[names(GCbck)==plotFacet1] <- "plotFacet1"
      GC$plotFacet1 = GCbck$plotFacet1
      results$plotFacet1 = resultsbck$plotFacet1

      names(resultsbck)[names(resultsbck)==plotFacet2] <- "plotFacet2"
      names(GCbck)[names(GCbck)==plotFacet2] <- "plotFacet2"
      GC$plotFacet2 = GCbck$plotFacet2
      results$plotFacet2 = resultsbck$plotFacet2
    }

  }

  if (is.na(legName)) {legName=plotCol}
  if (is.na(PlotName)) {PlotName=paste("Growth Curve", mdl, sep=", ")}




  #go over all positions
  for (i in 1:length(id)) {
    #get parameters for model plot, paste these to the temporary Gmp column
    if (Model=="Gmp") {
      K <- results$K[results$RowCol==id[i]]
      r <- results$r[results$RowCol==id[i]]
      b <- results$b[results$RowCol==id[i]]

      tm = GC$Time[GC$RowCol==id[i]]
      for (idT in 1:length(tm)) {
        GC$Gmp[GC$RowCol==id[i]&GC$Time==tm[idT]] = Gmprtz(K,r,b,tm[idT])
        GC$b[GC$RowCol==id[i]&GC$Time==tm[idT]] = b
      }

      #GC$Gmp[cnt:(cnt+length(GC$Time[GC$RowCol==id[1]])-1)] <- Gmprtz(K,r,b,)
    }else{
      K <- results$K[results$RowCol==id[i]]
      r <- results$r[results$RowCol==id[i]]
      v <- results$v[results$RowCol==id[i]]
      g <- results$g[results$RowCol==id[i]]
      tshift <- results$tshift[results$RowCol==id[i]]

      tm = GC$Time[GC$RowCol==id[i]]
      for (idT in 1:length(tm)) {
        GC$Gmp[GC$RowCol==id[i]&GC$Time==tm[idT]] = Glogist(K,r,g,v,tm[idT]-tshift)
      }

      #GC$Gmp[cnt:(cnt+length(GC$Time[GC$RowCol==id[1]])-1)] <- Glogist(K,r,g,v,GC$Time[GC$RowCol==id[i]]-tshift)
    }


    cnt=cnt+length(GC$Time[GC$RowCol==id[1]])
  }

  #time vector for plotting
  tvect=GC$Time[GC$RowCol==id[1]]

  GC <- droplevels(GC)
  NrCol <- length(unique(GC$plotCol))
  if (NrCol<3) {NrCol <- 3}
  if (NrCol<8) {pal <- RColorBrewer::brewer.pal(NrCol, "Set2")}else{pal <- RColorBrewer::brewer.pal(8, "Set2")}
  colVect <- colorRampPalette(pal)(NrCol)

  #actual plot:
  pl = ggplot2::ggplot(data=GC) +
    ggplot2::geom_point(ggplot2::aes(y=Growth,x=Time, color=plotCol),shape=4, size=3)+ #the measured point
    #geom_point(aes(y=Growth,x=Time),color="black",fill="black",shape=16, size=0.5)+ #the measured point
    ggplot2::geom_line(ggplot2::aes(x=Time, y=Gmp, color=plotCol, group=RowCol),size=1) + #model line, group them by rowcol
    ggplot2::scale_color_manual(values=colVect, name=legName) +
    #geom_point(aes(y=0.25, x=Time25), col="black", shape=3, size=2) +
    ggplot2::ylab(yname) + ggplot2::xlab(xname) + ggplot2::ggtitle(PlotName)

  if (length(plotFacet)==1) {
    pl = pl + ggplot2::facet_grid(~plotFacet) #faceting
  }else{
    pl = pl + ggplot2::facet_grid(plotFacet1~plotFacet2) #faceting
  }



  if (markb&Model=="Gmp") {
    pl = pl + ggplot2::geom_vline(data=GC,ggplot2::aes(xintercept=b, color=plotCol,group=RowCol),lwd=0.5)
  }

  pl
  return(pl)
}



#FitnessPlot----
#
#'Boxplots displying specified fitness measurements, colored and facetted by user-input
#'
#'This function creates ggplot boxplots visualizing fitness measurements. It is integrated into the makeFitness2 function as well.
#'There are three methods used to define fitness:
#'1) By multiplying maximal doubling potential (MDP) and maximal doubling rate (MDR) accoring to the original QFA work.
#'2) By divinding the maximal absolute growth rate (r) of the Gompertz model with the time to reach r (b).
#'3) By dividing the model-free estimation of r with the model free estimation of b. The user can specify which should be plotted by setting
#'the plotFitness parameters to either "All", "MDR" (1) "Gmp" (2) or "ModelFree" (3). Factor used on the x-axis is defined with plotX.
#'Additionally, coloring and facetting factor names can be specified with plotFill and plotFacet. If filename is specified,
#'the plots will be saved as a pdf with the specified name.
#'
#'@param results This is the output data.frame of the qfa.fit2 function
#'@param filename If this is anything else than NA, a pdf file will be saved into the current working directory named filename
#'@param plotX String. Name of the column of the results data.frame used as a factor seperating boxplots on the x-axis
#'@param plotoFacet String. Name of the column of the results data.frame used as a factor creaeting facets
#'@param plotFill String. Name of the column of the results data.frame used as a factor coloring the boxplots
#'@param plotFitness String. Either set to "All" to plot all three types of fitness measurements if the Gompertz model was used, otherwise
#'only MDR-based and modelfree fitness estimates are plotted. If set to "MDR" only MDR*MDP based fitness estimates are displayed.
#'If set to "Gmp" only Gompertz model based fitness estimates are displayed. If set to "ModelFree" only model free fitness estimates
#'are displayed.
#'
#' @return This function plots the ggplot object if not used with an assignment. If it is
#' called as an assignment, the return is a ggplot object which can further be modified by
#' the user if desired.
#'
#' @examples
#' #qfa.testdata was generated with the call in the Not run section
#' #(these files are not included in the package)
#' data(qfa.testdata)
#' #Strip non-experimental edge cultures
#' qfa.testdata = qfa.testdata[(qfa.testdata$Row!=1) & (qfa.testdata$Col!=1) & (qfa.testdata$Row!=8) & (qfa.testdata$Col!=12),]
#' # Define which measure of cell density to use
#' qfa.testdata$Growth = qfa.testdata$Intensity
#' GmpFit = qfa.fit2(qfa.testdata, inocguess=NULL, detectThresh=0, globalOpt=F, AUCLim=NA, TimeFormat="h", Model="Gmp")
#' # Construct fitness measures
#' GmpFit = makeFitness2(GmpFit, AUCLim=NA, plotFitness="no")
#' qfaFitnessPlot(GmpFit, plotFitness="All", filename="Example_Gmp_fitness.pdf")
#'
#' @keywords qfa
#'
qfaFitnessPlot <- function(results, plotFitness="no", filename = NA, plotX = "ORF", plotFacet = "ScreenID", plotFill = "Treatment") {
  #check if ggplot and cowplot are installed:
  qfa_check_packages()
  #to use same number of digits for all legends
  scaleFUN <- function(x) sprintf("%.4f", x)
  #for some functionality, create a backup of the results
  resBck <- results
  resBck1 <- results
  #get the model used
  Model <- unique(results$Model)

  if (is.na(plotX) | is.na(plotFacet) | is.na(plotFill)) {
    stop("Please specify all plotting variables")
  }

  if (!(plotX %in% names(results)) | !(plotFacet%in% names(results)) | !(plotFill %in% names(results))) {
    stop("Please specify all plotting variables")
  }



  #if the user wants a pdf, tell him to wait until all is done:
  if (!is.na(filename)) {
    print("Please wait until plots are finished exporting...")
    if (substr(filename, nchar(filename)-4+1, nchar(filename))!=".pdf") {
      filename = paste(filename,".pdf", sep="")
    }
  }


  #there is probably a better way to tell R which columns I want to plot... But that works as well
  if (plotX!=plotFill & plotX!=plotFacet & plotFacet!=plotFill) {
    names(results)[names(results) == plotX] <- "plotX"
    names(results)[names(results) == plotFill] <- "plotFill"
    names(results)[names(results) == plotFacet] <- "plotFacet"
  }else if (plotX==plotFill) {
    names(results)[names(results) == plotX] <- "plotX"
    names(resBck1)[names(resBck1)==plotFill] <- "plotFill"
    results$plotFill = resBck1$plotFill
    names(results)[names(results) == plotFacet] <- "plotFacet"
  }else if (plotX==plotFacet) {
    names(results)[names(results) == plotX] <- "plotX"
    names(resBck1)[names(resBck1)==plotFacet] <- "plotFacet"
    results$plotFacet = resBck1$plotFacet
    names(results)[names(results) == plotFacet] <- "plotFacet"
  }else if (plotFacet==plotFill) {
    names(results)[names(results) == plotFill] <- "plotFill"
    names(resBck1)[names(resBck1)==plotFacet] <- "plotFacet"
    results$plotFacet = resBck1$plotFacet
    names(results)[names(results) == plotX] <- "plotX"
  }else{
    resBck2=resBck1
    names(results)[names(results) == plotX] <- "plotX"
    names(resBck1)[names(resBck1)==plotFill] <- "plotFill"
    results$plotFill = resBck1$plotFill
    names(resBck2)[names(resBck2)==plotFill] <- "plotFacet"
    results$plotFacet = resBck2$plotFacet
  }



  results <- droplevels(results)
  NrCol <- length(unique(results$plotFill))
  if (NrCol<3) {NrCol <- 3}
  if (NrCol<8) {pal <- RColorBrewer::brewer.pal(NrCol, "Set2")}else{pal <- RColorBrewer::brewer.pal(8, "Set2")}
  colVect <- colorRampPalette(pal)(NrCol)

  pl1 <- ggplot2::ggplot(results) +
    ggplot2::geom_boxplot(ggplot2::aes(y = MDRMDP, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::ggtitle("MDR*MDP Fitness") +
    ggplot2::xlab("") +
    ggplot2::ylab("MDR * MDP") +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_grid(plotFacet ~ ., scales = "free_x", space = "free_x", shrink = T) +
    ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  pl2 <- ggplot2::ggplot(results) +
    ggplot2::geom_boxplot(ggplot2::aes(y = MDR, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::ggtitle("MDR") +
    ggplot2::xlab("") +
    ggplot2::ylab("MDR") +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_grid(plotFacet ~ ., scales = "free_x", space = "free_x", shrink = T)+
    ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  pl3 <- ggplot2::ggplot(results) +
    ggplot2::geom_boxplot(ggplot2::aes(y = MDP, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::ggtitle("MDP") +
    ggplot2::xlab("") +
    ggplot2::ylab("MDP") +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_grid(plotFacet ~ ., scales = "free_x", space = "free_x", shrink = T)+
    ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (Model == "Gmp") {
    pl4 <- ggplot2::ggplot(results) +
      ggplot2::geom_boxplot(ggplot2::aes(y = r/b, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
      ggplot2::ggtitle("GMP fitness") +
      ggplot2::xlab("") +
      ggplot2::ylab(paste("fitness: r/b [AU/",unique(results$TimeFormat), "^2]", sep = "")) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::facet_grid(plotFacet ~ ., scales = "free_x", space = "free_x", shrink = T) +
      ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    pl5 <- ggplot2::ggplot(results) +
      ggplot2::geom_boxplot(ggplot2::aes(y = r, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
      ggplot2::ggtitle("GMP max r") +
      ggplot2::xlab("") +
      ggplot2::ylab(paste("max r [AU/", unique(results$TimeFormat), "]", sep = "")) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::facet_grid(plotFacet ~ ., scales = "free_x", space = "free_x", shrink = T) +
      ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    pl6 <- ggplot2::ggplot(results) +
      ggplot2::geom_boxplot(ggplot2::aes(y = b, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
      ggplot2::ggtitle("GMP time to max r") +
      ggplot2::xlab("") +
      ggplot2::ylab(paste("time to max r: b [", unique(results$TimeFormat), "]", sep = "")) +
      ggplot2::scale_y_continuous(labels = scaleFUN) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::facet_grid(plotFacet ~ ., scales = "free_x", space = "free_x", shrink = T) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  }


  pl7 <- ggplot2::ggplot(results) +
    ggplot2::geom_boxplot(ggplot2::aes(y = maxslp/maxslp_t, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::ggtitle("model-free \nfitness") +
    ggplot2::xlab("") +
    ggplot2::ylab(paste("r*/b* [AU/",unique(results$TimeFormat), "^2]", sep = "")) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_grid(plotFacet ~ ., scales = "free_x",space = "free_x", shrink = T) +
    ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  pl8 <- ggplot2::ggplot(results) +
    ggplot2::geom_boxplot(ggplot2::aes(y = maxslp, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::ggtitle("model-free \nmax slope r*") +
    ggplot2::xlab("") +
    ggplot2::ylab(paste("r* [AU/",unique(results$TimeFormat), "]", sep = "")) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_grid(plotFacet ~ ., scales = "free_x",space = "free_x", shrink = T) +
    ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  pl9 <- ggplot2::ggplot(results) +
    ggplot2::geom_boxplot(ggplot2::aes(y = maxslp_t, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::ggtitle("model-free \ntime to max slope b*") +
    ggplot2::xlab("") +
    ggplot2::ylab(paste("b* [", unique(results$TimeFormat), "]", sep = "")) +
    ggplot2::scale_y_continuous(labels = scaleFUN) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_grid(plotFacet ~., scales = "free_x", space = "free_x", shrink = T) +
    ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))




  plLeg <- ggplot2::ggplot(results) +
    ggplot2::geom_boxplot(ggplot2::aes(y = maxslp_t, x = plotX, fill = plotFill), color = "black", position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::ggtitle("MDR*MDP Fitness") +
    ggplot2::xlab("") +
    ggplot2::ylab("MDR * MDP") +
    ggplot2::scale_y_continuous(labels = scaleFUN) + ggplot2::labs(fill = plotFill) +
    ggplot2::facet_grid(plotFacet ~ ., scales = "free_x", space = "free_x", shrink = T) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = ceiling(length(unique(results$plotFill))/6))) +
    ggplot2::scale_fill_manual(values=colVect, name=plotFill) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  while (!is.null(dev.list())) dev.off()
  leg <- cowplot::get_legend(plLeg)


  hVar = length(unique(results$plotFacet))
  if (hVar>=1) hVar = 2
  if (plotFitness == "All") hVar=hVar*3
  hVar2 = hVar * 2
  if (hVar2 < 8) hVar2 = 8
  results = resBck



  if (!is.na(filename)) {
    pdf(filename, width = 10, height = 3.5*hVar)
  }

  if (Model == "Gmp") {
    if (plotFitness == "All") {
      # grid.arrange(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,pl9, ncol=3)
      grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
    } else if (plotFitness == "MDR") {
      grd1 = gridExtra::grid.arrange(pl1, pl2, pl3, ncol = 3)
    } else if (plotFitness == "Gmp") {
      # grid.arrange(pl4,pl5,pl6, ncol=3)
      grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl4, pl5, pl6, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
    } else if (plotFitness == "ModelFree") {
      grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl7, pl8, pl9, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
      # grid.arrange(pl7,pl8,pl9, ncol=3)
    } else {
      print("Wrong fitness type defined for plotting")
    }

  } else {
    if (plotFitness == "All") {
      # grid.arrange(pl1,pl2,pl3,pl7,pl8,pl9, ncol=3)
      grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl1, pl2, pl3, pl7, pl8, pl9, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
    } else if (plotFitness == "MDR") {
      # grid.arrange(pl1,pl2,pl3, ncol=3)
      grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl1, pl2, pl3, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
    } else if (plotFitness == "ModelFree") {
      # grid.arrange(pl7,pl8,pl9, ncol=3)
      grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl7, pl8, pl9, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
    } else {
      print("Wrong fitness type defined for plotting")
    }
  }


  if (!is.na(filename)) {
    dev.off()
    if (Model == "Gmp") {
      if (plotFitness == "All") {
        # grid.arrange(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,pl9, ncol=3)
        grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
      } else if (plotFitness == "MDR") {
        grd1 = gridExtra::grid.arrange(pl1, pl2, pl3, ncol = 3)
      } else if (plotFitness == "Gmp") {
        # grid.arrange(pl4,pl5,pl6, ncol=3)
        grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl4, pl5, pl6, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
      } else if (plotFitness == "ModelFree") {
        grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl7, pl8, pl9, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
        # grid.arrange(pl7,pl8,pl9, ncol=3)
      } else {
        print("Wrong fitness type defined for plotting")
      }

    } else {
      if (plotFitness == "All") {
        # grid.arrange(pl1,pl2,pl3,pl7,pl8,pl9, ncol=3)
        grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl1, pl2, pl3, pl7, pl8, pl9, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
      } else if (plotFitness == "MDR") {
        # grid.arrange(pl1,pl2,pl3, ncol=3)
        grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl1, pl2, pl3, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
      } else if (plotFitness == "ModelFree") {
        # grid.arrange(pl7,pl8,pl9, ncol=3)
        grd1 = gridExtra::grid.arrange(gridExtra::arrangeGrob(pl7, pl8, pl9, ncol = 3), leg, nrow = 2, heights = c(hVar2, 1))
      } else {
        print("Wrong fitness type defined for plotting")
      }
    }
  }
  if (!is.na(filename)) print("Plot exporting finished!")
  return(grd1)

}




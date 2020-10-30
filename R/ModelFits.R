#' Growth curve modelling
#'
#' This function will fit the specified growth model to timecourse
#' observations by least squares using either the L-BFGS-B
#' algorithm in R's optim function, or the differential evolution, stochastic
#' global optimisation package DEoptim. The input needs to be in the correct format
#' preferably created with either the colonyzer.read or GC2QFA function.
#' The user can specify to use one of the three following growth models with the argument Model:
#' Standard logistic model (Slog), Generalized linear model (Glog) or Gompertz
#' growth model (Gmp).  qfa.fit2 will also calculate a numerical
#' Area Under Curve (nAUC) fitness measure by integrating under a loess
#' smooothed version of the dataset if there are sufficient observations or
#' under a linear interpolation between observations if observations are too
#' infrequent.
#'
#'
#' @param d The data.frame containing the timecourse data for each colony
#' (returned from colonyzer.read or GC2QFA).
#' @param Model Either "Slog" for standard logistic model, "Glog" for generalized
#' logistic model or "Gmp" for Gompertz model. The parameters returned by qfa.fit2
#' will change corresponding to the chosen model. The Gmp variation used is formula (22)
#' from Tjorve and Tjorve 2017 (PLOS ONE) which is a unified version with an
#' absolute growth rate
#' @param inocguess Only relevant for Slog and Glog.
#' Should be either numerical or NULL.
#' The best guess for starting density of viable cells in each
#' colony.  This is the g parameter in the Slog and Glog.
#' Typically, for dilute inoculum 384 format spotted cultures, this value
#' cannot be observed directly by photography.  inocguess should be in the same
#' units as the values in the Growth column in d. If fixG=TRUE, only values of
#' g within the range 0.9*inocguess and 1.1*inocguess will be assessed during
#' optimisation.  Otherwise values within 1e-10*inocguess and 1e+10*inocguess
#' will be tried.  Without a sensible independent estimate for inoculum
#' density, the best we can do is to estimate it based on observed data.
#' Estimating inocguess happens if inocguess is set to NULL.
#' Estimating inoculum density will only work well if the inoculum density is
#' high enough to be measurable (e.g. pinned cultures or conc. spotted) and is
#' clearly observed.  Clearly observed means: no condensation on plates
#' immediately after they are placed in incubator for example.  If we are
#' making an independent estimate of inoculum density, then we should also
#' reset the time at which the experiment "begins".  This experiment start time
#' should be the time at which the inoculum density is observed.
#' @param fmt The date.time format that the inoculation time (Inoc.Time) and
#' measurement times (Date.Time) are stored in. Default to "%Y-%m-%d_%H-%M-%S"
#' which is the output format of Colonyzer.
#' @param TimeFormat Either set to hours (h) or days (d). Only defines the format
#' of the outputs in growth rates and for subsequent plotting.
#' @param minK The minimum value of K above which a strain is said to be alive.
#' Strains with K optimised to lie below this value will be classified as dead,
#' by setting r to be zero.
#' @param detectThresh The minimum detectable cell density (or Growth value)
#' which reliably identifies the presence of cells. Cell densities below this
#' value are classified as noise and repalced with the detectThresh value.
#' Can also be set =0. Then, the software is trying to estimate a threshold as
#' the mean between the two smallest Growth values per position.
#' @param globalOpt logical. Indicates whether qfa.fit2 should use the slower, but
#' more robust DEoptim global optimisation functions to fit the growth
#' model to the data, or the quicker optim function.
#' @param logTransform logical. Indicating if data should be log-transformed before
#' model fit. Not recommended.
#' @param fixG logical. Only relevant for Slog and Glog.
#' Indicates whether to allow g parameter to vary over a wide 1e-10*inocguess to 1e+10*inocguess
#' or narrow range 0.9*inocguess to 1.1*inocguess during optimisation.
#' fixG=TRUE corresponds to narrow constraints on g.
#' @param AUCLim Numerical AUC (nAUC) is calculated as the integral of an
#' approximation of the growth curve between time 0 and AUCLim. If set to NA (default),
#' AUClim will be set to the maximum time in the dataset.
#' @param STP Time to use for "Single Timepoint" fitness estimate.
#' Defaults to 20 days (very late in growth curve) which is like carrying
#' capacity. Untested functionality of the first QFA package.
#' @param nCores Can attempt to split model fitting load across multiple
#' parallel cores.  Experimental, probably best to leave this value set to
#' default (1).
#' @param modelFit logical. Specifies whether to carry out any
#' model fitting at all.  When set to FALSE, only numerical fitness estimates
#' such as nr, nMDP, nAUC are generated
#' @param checkSlow logical. Specifies whether to re-optimise
#' curve-fitting for slow-growing strains.  If TRUE, slow-growing or dead
#' strains are identified heuristically and a second round of curve fitting
#' using global (but slower) optimisation is carried out.  Heuristic
#' identification of slow-growing strains is currently experimental, it seems
#' we have over-tuned these to datasets we capture at Newcastle.  If you notice
#' a banding pattern in your MDR or r fitness distributions, please set
#' checkSlow to FALSE.
#' @param nrate Boolean specifiying whether to include numerical fitness estimates
#' like maximum growth rate and time to reach this maximum growth rate. These
#' estimates are derived from model-free numerical integration of the data
#' @param lowK,upK,etc Set the lower and upper boundaries of the optimisation
#' algorithm for the corresponding model parameters. Parameter possibilites:
#' K, r, g, b, v
#' @param ... Extra arguments passed to optim
#'
#' @return R data.frame, similar to that returned by the colonyzer.read
#' function.  The major difference is that instead of a row for every cell
#' density observation for every culture, this object summarises all timecourse
#' density observations for each culture with fitted grwoth model
#' parameters and numerical fitness estimates.
#'
#' \itemize{ \item Barcode - Unique plate identifier \item Row - Row number
#' (counting from top of image) of culture in rectangular gridded array \item
#' Col - Column number (counting from left of image) of culture in rectangular
#' gridded array \item ScreenID - Unique identifier for this QFA screen \item
#' Treatment - Conditions applied externally to plates (e.g. temperature(s) at
#' which cultures were grown, UV irradiation applied, etc.) \item Medium -
#' Nutrients/drugs in plate agar \item ORF - Systematic, unique identifier for
#' genotype in this position in arrayed library \item Screen.Name - Name of
#' screen (identifies biological repeats, and experiment) \item Library.Name -
#' Name of library, specifying particular culture location \item MasterPlate
#' Number - Library plate identifier \item Timeseries order - Sequential
#' photograph number \item Inoc.Time - User specified date and time of
#' inoculation (specified in ExptDescription.txt file) \item TileX - Culture
#' tile width (pixels) \item TileY - Culture tile height (pixels) \item XOffset
#' - x-coordinate of top left corner of rectangular tile bounding culture
#' (number of pixels from left of image) \item YOffset - y-coordinate of top
#' left corner of rectangular tile bounding culture (number of pixels from top
#' of image) \item Threshold - Global pixel intensity threshold used for image
#' segmentation (after lighting correction) \item EdgeLength - Number of
#' culture pixels classified as being microcolony edge pixels (useful for
#' classifying contaminants in cultures grown from dilute inoculum) \item
#' EdgePixels - Number of pixels classified as culture on edge of square tile
#' \item RepQuad - Integer identifying which of the quadrants of a 1536 plate
#' were used to inoculate the current 384 plate (set equal to 1 for all
#' cultures for 1536 format for example) \item K - carrying capacity (upper asymptote) for all models
#' \item r - Rate parameters (Slog and Glog), maximum absolute growth rate (Gmp)
#' \item g - For Slog and Glog: inoculum density (lower asymptote) (referred to in vignette as
#' g_0). For Gmp: Calculated based on the three Gmp parameters
#' \item v - Only for Slog and Glog: Generalised logistic model shape parameter (=1 for Slog)
#' \item b Only for Gmp: Time to reach max growth rate (r).
#' \item yshift Only for Gmp: Min growth (data is shifted down by that amount for Gmp fit)
#' \item objval - Objective function value (sum of squares) at selected optimum.
#'  \item tshift - Only for Slog and Glog: Shift applied to observation times before fitting
#' logistic model (need to apply same shift before overlaying curve on expt.
#' obs.).  Set to first timepoint at which growth value is equal or above to detectThresh
#' \item t0 - Time of first detectable cell
#' density observation (i.e. above detectThresh) \item d0 - Normalised cell
#' density of first observation (be careful about condensation on plates when
#' using this).  Note this is not necessarily the density at t0. \item nAUC -
#' Numerical Area Under Curve.  This is a model-free fitness estimate. \item
#' nSTP - Single Time Point fitness.  Cell density at time STP, as estimated
#' with approximating function.  This is a model-free fitness estimate. \item
#' nr - Numerical estimate of intrinsic growth rate.  Growth rate estimated by
#' fitting smoothing function to log of data, calculating numerical slope
#' estimate across range of data and selecting the maximum estimate (should
#' occur during exponential phase). \item nr_t - Time at which maximum slope of
#' log observations occurs \item maxslp - Numerical estimate of maximum slope
#' of growth curve.  Slope estimated by fitting smoothing function to
#' untransformed data and calculating numerical slope estimate of smoothed
#' version of data and selecting the maximum estimate (should occur
#' approximately half way through growth).  This fitness measure will be
#' affected by both rate of growth and final colony size.  Final colony size is
#' expected to be strongly affected by competition between cultures. \item
#' maxslp_t - Time at which maximum slope of observations occurs \item ExptDate - A
#' representative/approximate date for the experiment (note that genome-wide
#' QFA screens typically take weeks to complete) \item User - Person who
#' actually carried out screen \item PI - Principal investigator leading
#' project that screen is part of \item Condition - The most important defining
#' characteristic of screen, as specified by user (e.g. the temperature screen
#' was carried out at if screen is part of multi-temperature set of screens, or
#' the query mutation if part of a set of screens comparing query mutations, or
#' the drugs present in the medium if part of a set of drug screens) \item Inoc
#' - Qualitative identifier of inoculation type (e.g. "DIL" for dilute
#' inoculum, "CONC" for concentrated).  Used to distinguish between experiments
#' carried out with different methods of inoculation. \item Gene - Identifier
#' for genotype at a particular location on an agar plate.  Typically prefer
#' unambiguous, systematic gene names here. \item TrtMed - Combination of
#' treatment and medium identifiers, specifying the environment in which the
#' cells have grown }
#' @keywords qfa.fit2
#' @examples
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
qfa.fit2 <- function(d,  Model = "Gmp", TimeFormat = "d", inocguess = NULL, fmt = "%Y-%m-%d_%H-%M-%S", minK = 0.025, detectThresh = 5e-04, globalOpt = FALSE, logTransform = FALSE, fixG = TRUE, AUCLim = NA,
    STP = 20, nCores = 1, modelFit = TRUE, checkSlow = F, nrate = T, lowK = NA, upK = NA, lowr = NA, upr = NA, lowg = NA, upg = NA,
    lowv = NA, upv = NA, lowb = NA, upb = NA, ...) {
    # this is the main function of qfa. It estimates a growth model according to the user input and outputs a dataframe consisiting of relevant model parameters and all
    # the necessary auxilliary information the user may want to know. It's column called Growth is going to be used as targetvariable input parameters:
    # d: the output of
    # either colonyzer.read or GC2QFA or any other dataframe that is formatted the correct way to work with the qfa.fit2
    # fmt: specify the timeformat of Inoc.Time and
    # Exp.Time
    # minK:min value to consider as growing and not dead.
    # detectThresh: minimum detectable Growth value. Everything below is considered noise and is going to be
    # repalced with detectThresh Value. Can also be set =0. Then, the software is trying to estimate a threshold as the mean between the two smallest Growth values per
    # position.
    # globalOpt: =F: use optim function. =T: use DEoptim (more robust but slower)
    # logTransform: =F: do not transform data before model optimization. =T: do log
    # transformation. Not recommended.
    # modelFit: =T: do modelfits. =F: do not.
    # Model: ='Gmp': use gompertz model with parameters K (max Growth, carrying capacity, upper
    # asymptote), r (max absolute growth rate), b (time to reach max growth rate) ='Slog': Standard logistic model with parameters K(max Growth, carrying capacity, upper
    # asymptote), r (growth rate), g (inoculation Growth, lower asymptote) ='Glog': generalized logistic model with parameters K(max Growth, carrying capacity, upper
    # asymptote), r (growth rate), g (inoculation Growth, lower asymptote), and v (shape parameter, describes asymmetry of growth curve) inocguess: needed for Glog and
    # Slog models to estimate the starting value of g for optimization. If =NULL: try to estimate it from given data fixG: =T: narrow spectrum of varying G parameter for
    # Glog model. =F: wider spectrum AUCLim: manual specification where the limit for area under curve calculation is. If =NA, the max time of the experiment is taken as
    # AUCLim STP: Time for single image calculation. Never tested. nCores: number of cores if multithreading is possible. Avoid problems by leaving this set to 1
    # checkSlow: over optimised call for slow growing colonies re-fitting. Best to leave it =F nrate: =T: calculate numerical fitnes values (like max slope and time to
    # max slope) TimeFormat: specify if the output values (like growth rate and time to max growth rate) should be coded in h or d all lowX, upX variables: manually set
    # boundaries for the optimization functions. Lower upr to values around 0.01 if TimeFormat='h' can sometimes improve fit

  if (! "Growth" %in% names(d)){
    stop("Please indicate the target growth variable by a naming a column 'Growth'")
  }

    # first check if the model is specified correctly, set some control variables
    if (Model == "Gmp") {
        Gmp = T
        fixG = F
        glog = F
        print("Gompertz model used")
    } else if (Model == "Glog") {
        Gmp = F
        glog = T
        print("Generalized logistic model used")
    } else if (Model == "Slog") {
        Gmp = F
        glog = F
        print("Standard logistic model used")
    } else {
        stop("Please specify the model correctly with Model= one of the following: Gmp (Gompertz), Glog (Generalized logistic), Slog (standard logistic).
              Do not forget quotation marks!")
    }

    # check if timeformat is specfied
    if (TimeFormat != "d" & TimeFormat != "h") {
        stop("Please specify the TimeFormat as either d or h.")
    }


    if (!"Column" %in% colnames(d))
        d$Column = d$Col


    # check if an inocguess is specified
    if (!is.null(inocguess)) {
        if (length(inocguess) == 0) {
            print("ERROR: must specify an inoculum density guess (or NULL)")
            return()
        }
    }

    if (nCores > 1) {
        cl = makeCluster(nCores)
        clusterCall(cl, function() library(qfa))
    } else {
        cl = NULL
    }

    if ((!is.na(lowK) & !is.na(upK)) && upK<=lowK) {
      stop("The boundaries for K are incorrect. Please correct")
      return()
    }

    if ((!is.na(lowr) & !is.na(upr)) && upr<=lowr) {
      stop("The boundaries for r are incorrect. Please correct")
      return()
    }

    if ((!is.na(lowg) & !is.na(upg)) && upg<=lowg) {
      stop("The boundaries for r are incorrect. Please correct")
      return()
    }

    if ((!is.na(lowv) & !is.na(upv)) && upv<=lowv) {
      stop("The boundaries for r are incorrect. Please correct")
      return()
    }

    if ((!is.na(lowg) & !is.na(upb)) && upb<=lowb) {
      stop("The boundaries for r are incorrect. Please correct")
      return()
    }

    # Rename columns if necessary
    if (!"Tile.Dimensions.X" %in% colnames(d))
        d[, "Tile.Dimensions.X"] = d$Diameter
    if (!"Tile.Dimensions.Y" %in% colnames(d))
        d[, "Tile.Dimensions.Y"] = d$Diameter
    if (!"X.Offset" %in% colnames(d))
        d[, "X.Offset"] = round(d$x - d$Diameter/2)
    if (!"Y.Offset" %in% colnames(d))
        d[, "Y.Offset"] = round(d$y - d$Diameter/2)
    if (!"Edge.length" %in% colnames(d))
        d[, "Edge.length"] = d$Perimeter
    if (!"Edge.Pixels" %in% colnames(d))
        d[, "Edge.Pixels"] = d$Perimeter

    # check if there is a normalization or not
    if ("Normalization" %in% names(d)) {
        Nrm = unique(d$Normalization)
    } else {
        Nrm = "Raw"
    }


    if (is.na(AUCLim)) {
        if (TimeFormat == "d") {
            AUCLim = max(d$Expt.Time, na.rm=T)
        } else {
            AUCLim = max(d$Expt.Time, na.rm=T) * 24
        }
    }

    if (is.na(AUCLim)) {
      AUCauto = T
    }else{
      AUCauto = F
    }

    # Vector of barcodes
    barcodes <- unique(d$Barcode)
    nbc <- length(barcodes)
    # Get big data frame ready for results
    results <- data.frame()
    # For each barcode, optimize
    bcount <- 0
    for (bcode in barcodes) {

      if (AUCauto) {
        if (TimeFormat == "d") {
          AUCLim = max(d$Expt.Time)
        } else {
          AUCLim = max(d$Expt.Time) * 24
        }
      }


        bcount <- bcount + 1
        print(paste("Optimizing Plate", bcount, "/", nbc, ":", bcode))
        # Restrict data to just this barcode
        dbc <- d[d$Barcode == bcode, ]

        if (AUCauto) {
          if (TimeFormat == "d") {
            AUCLim = max(dbc$Expt.Time)
          } else {
            AUCLim = max(dbc$Expt.Time) * 24
          }
        }
        inoctime = dbc$Inoc.Time[1]
        # Get list of unique colony positions
        positions <- lapply(1:length(dbc[, 1]), index2pos, dbc)
        positions <- unique(positions)
        # Fit logistic model to each colony
        if (nCores > 1) {
            print(as.vector(lsf.str(envir = .GlobalEnv)))
            ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
            clusterExport(cl, ex)
            clusterExport(cl, as.vector(lsf.str(envir = .GlobalEnv)))
            bcfit <- t(parSapply(cl, positions, colony.fit2, dbc, inocguess, fixG, globalOpt, detectThresh, minK, logTransform, AUCLim, STP, glog=glog, modelFit, checkSlow,
                nrate, TimeFormat = TimeFormat, Gmp = Gmp, lowK = lowK, upK = upK, lowr = lowr, upr = upr, lowg = lowg, upg = upg, lowv = lowv, upv = upv, lowb = lowb,
                upb = upb, Nrm = Nrm))
        } else {
            # this is the call for the fiting of growth model, pass all arguments down to that
            bcfit <- t(sapply(positions, colony.fit2, dbc, inocguess, fixG, globalOpt, detectThresh, minK, logTransform, AUCLim, STP, glog=glog, modelFit, checkSlow, nrate,
                TimeFormat = TimeFormat, Gmp = Gmp, lowK = lowK, upK = upK, lowr = lowr, upr = upr, lowg = lowg, upg = upg, lowv = lowv, upv = upv, lowb = lowb, upb = upb,
                Nrm = Nrm))
        }
        # extract all the meta data things
        info <- data.frame(t(sapply(positions, colony.info, dbc)))
        rows <- sapply(positions, rcget, "row")
        cols <- sapply(positions, rcget, "col")
        # Bind Data frame of barcode results to overall results
        barcMetadata <- data.frame(Barcode = as.character(info$Barcode), Row = rows, Column = cols, Col = cols, ScreenID = as.character(info$ScreenID), Treatment = as.character(info$Treatments),
            Medium = as.character(info$Medium), ORF = as.character(info$ORF), Screen.Name = as.character(info$Screen.Name), Library.Name = as.character(info$Library.Name),
            MasterPlate.Number = as.numeric(info$MasterPlate.Number), Timeseries.order = as.numeric(info$Timeseries.order), Inoc.Time = inoctime, TileX = as.numeric(info$Tile.Dimensions.X),
            TileY = as.numeric(info$Tile.Dimensions.Y), XOffset = as.numeric(info$X.Offset), YOffset = as.numeric(info$Y.Offset), Threshold = as.numeric(info$Threshold),
            EdgeLength = as.numeric(info$Edge.length), EdgePixels = as.numeric(info$Edge.Pixels), RepQuad = as.numeric(info$RepQuad), stringsAsFactors = FALSE)
        barcFitness <- data.frame(bcfit)
        barcResults = cbind(barcMetadata, barcFitness)
        if ("Client" %in% colnames(info)) {
            barcResults$Client = as.character(info$Client)
        }
        if ("ExptDate" %in% colnames(info)) {
            barcResults$ExptDate = as.character(info$ExptDate)
        }
        if ("User" %in% colnames(info)) {
            barcResults$User = as.character(info$User)
        }
        if ("PI" %in% colnames(info)) {
            barcResults$PI = as.character(info$PI)
        }
        if ("Condition" %in% colnames(info)) {
            barcResults$Condition = as.character(info$Condition)
            barcResults$Condition[is.na(barcResults$Condition)] = ""
        }
        if ("Inoc" %in% colnames(info)) {
            barcResults$Inoc = as.character(info$Inoc)
        }
        if ("Gene" %in% colnames(info)) {
            barcResults$Gene = as.character(info$Gene)
        }

        # add some more things to the output and done!

        barcResults$Model <- Model
        barcResults$TimeFormat <- TimeFormat
        barcResults$Normalization <- Nrm
        barcResults$AUCLim = AUCLim

        results = rbind(results, barcResults)
    }  #bcode
    if (nCores > 1) {
        stopCluster(cl)
    }
    print("QFA model fit finished")
    return(results)

}

#' Internal QFA function called by qfa.fit2. Do not call directly
colony.fit2 <- function(position, bcdata, inocguess, fixG = TRUE, globalOpt = FALSE, detectThresh = 0, minK = 0, logTransform = FALSE, AUCLim = 5, STP = 10, glog = TRUE,
    modelFit = TRUE, checkSlow = TRUE, nrate = TRUE, TimeFormat = "d", Gmp = F, lowK = NA, upK = NA, lowr = NA, upr = NA, lowg = NA, upg = NA, lowv = NA, upv = NA, lowb = NA,
    upb = NA, Nrm = "Raw", ...) {
    # next step in calling model fit stuff

    # Get row & column to restrict data
    row <- position[1]
    col <- position[2]

    # subset the data for one position
    do <- bcdata[(bcdata$Row == row) & (bcdata$Column == col), ]
    # if the first value is bigger (sometimes happens due to our setup), set it to the value of the second
    if (!is.na(do$Growth[1]) && !is.na(do$Growth[2])) {
      if (do$Growth[1] > do$Growth[2]) {
        do$Growth[1] = do$Growth[2]
    }
    }


    # now adjust the timeformat
    if (TimeFormat == "h") {
        obsdat = data.frame(Expt.Time = as.numeric(do$Expt.Time * 24), Growth = as.numeric(do$Growth))
    } else {
        obsdat = data.frame(Expt.Time = as.numeric(do$Expt.Time), Growth = as.numeric(do$Growth))
    }


    # call for the model fit:
    pars = makefits2(obsdat, inocguess, fixG, globalOpt, detectThresh, minK, logTransform, AUCLim, STP, glog=glog, modelFit, checkSlow, nrate, TimeFormat = TimeFormat, Gmp = Gmp,
        lowK = lowK, upK = upK, lowr = lowr, upr = upr, lowg = lowg, upg = upg, lowv = lowv, upv = upv, lowb = lowb, upb = upb, Nrm = Nrm)
    return(pars)
}


#' Internal QFA function called by qfa.fit2. Do not call directly
makefits2 <- function(obsdat, inocguess, fixG = TRUE, globalOpt = FALSE, detectThresh = 0, minK = 0, logTransform = FALSE, AUCLim = 5, STP = 10, glog = TRUE, modelFit = TRUE,
    checkSlow = TRUE, nrate = TRUE, TimeFormat = "d", Gmp = F, lowK = NA, upK = NA, lowr = NA, upr = NA, lowg = NA, upg = NA, lowv = NA, upv = NA, lowb = NA, upb = NA,
    Nrm = "Raw", ...) {
    # next call for model fit things


    # get numerical estimates of max slope & time to max slope:
    numfit = numericalfitness2(obsdat, AUCLim, STP, nrate = nrate)


    if (modelFit) {
        # can also use the function without... next call to the model fit things
        pars = growthcurve2(obsdat, inocguess, fixG, globalOpt, detectThresh, minK, logTransform, glog=glog, checkSlow, TimeFormat = TimeFormat, Gmp = Gmp, lowK = lowK, upK = upK,
            lowr = lowr, upr = upr, lowg = lowg, upg = upg, lowv = lowv, upv = upv, lowb = lowb, upb = upb, Nrm = Nrm)
    } else {
        # just replace all with NA
        pars = c(max(obsdat$Growth), NA, NA, NA, NA, NA)
        names(pars) = c("K", "r", "g", "v", "objval", "tshift")
    }
    # Add on time of first valid obs. and numerical AUC, STP
    valid = obsdat[obsdat$Growth >= detectThresh, ]
    if (dim(valid)[1] > 0) {
        t0 = min(valid$Expt.Time)
        # Calculate mean here in case multiple data for same timepoint...
        d0 = mean(valid$Growth[valid$Expt.Time == t0])
    } else {
        t0 = Inf
        d0 = NA
    }
    parsadd = c(t0, d0)
    names(parsadd) = c("t0", "d0")
    pars = c(pars, parsadd, numfit)
    # if(length(pars)!=9) {dput(pars); dput(obsdat)}
    return(pars)
}


#' Internal QFA function called by qfa.fit2. Do not call directly
growthcurve2 <- function(obsdat, iguess, fixG = TRUE, globalOpt = FALSE, detectThresh = 0, minK = 0, logTransform = FALSE, glog = TRUE, checkSlow = TRUE, TimeFormat = "d",
    Gmp = F, lowK = NA, upK = NA, lowr = NA, upr = NA, lowg = NA, upg = NA, lowv = NA, upv = NA, lowb = NA, upb = NA, Nrm = "Raw", ...) {

    # next step in growth model calculation.
    len1 = length(obsdat$Growth)

    # only for Glog and Slog:
    if (!Gmp) {
        # Throw away observations below the detectable threshold
        d = obsdat[obsdat$Growth >= detectThresh, ]
    } else {
        d = obsdat
    }

    # Throw away observations occurring before inoculation date (negative times...) ?from qfa V1. Do not know if that is necessary or why that should happen...
    d = d[d$Expt.Time >= 0, ]
    d = d[!is.na(d$Growth), ]
    len2 = length(d$Growth)
    # If this has left us with too few points, return 'dead colony'
    if ((len2/len1 < 0.25) | (len2 < 3)) {
        if (!is.null(iguess)) {
            pars = c(iguess, 0, iguess, 1, Inf, 0)
        } else {
            pars = c(minK, 0, minK, 1, Inf, 0)
        }
        names(pars) = c("K", "r", "g", "v", "objval", "tshift")
        return(pars)
    }


    tshift = 0
    # if detectThresh not =0, set all measurements smaller or equal to detectThresh to detectThresh. also needs to define the time shift for Glog and Slog models. The
    # fit should only start from the time where we actually have values otherwise the Glog and Slog give strange results as they start from t=0 with a paramter defining
    # the shift on x-axis (time)
    if (detectThresh != 0) {
        d$Growth[d$Growth <= detectThresh] = detectThresh
        # also, the first measurement is sometimes bigger than the rest, set that to detectThresh even if it is bigger
        if (d$Growth[1] > detectThresh) {
            d$Growth[1] = detectThresh
        }
        # calculate the timeshift for Glog and Slog models. Start the fit from there. First value bigger than detectThresh
        tshift = d$Expt.Time[min(which(d$Growth > detectThresh))]
    } else {
        # if detectThresh=0, take the first value bigger than 10% of the smallest value (ignoring first entry) as the value to detect tshift
        tst <- d$Growth[2:length(d$Growth)]  #remove first value
        candidate = NA
        if (sum(tst > 1.1 * min(tst)) == 0) {
            candidate = min(tst)
        } else {
            candidate = min(tst[tst > 1.1 * min(tst)], na.rm = T)
        }


        if (is.na(candidate))
            candidate = tst[1]
        ind <- match(candidate, d$Growth)  #get index

        # if that did not return the very first entry, take one before. Because, as soon as measurements are above that value, we consider them as real data
        if (ind != 1) {
            ind = ind - 1
        }
        tshift = d$Expt.Time[ind]
    }


    if (!Gmp) {
        if (is.null(iguess)) {
            # Without a sensible independent estimate for inoculum density, the best we can do is to estimate it based on observed data.  This strategy will only work well if
            # the inoculum density is high enough to be measurable (e.g. pinned cultures or conc. spotted) and is clearly observed.  Clearly observed means: no condensation on
            # plates immediately after they are placed in incubator for example.  If we are making an independent estimate of inoculum density, then we should also reset the
            # time at which the experiment 'begins'.  This experiment start time should be the time at which the inoculum density is observed.

            # do we have some data:
            if (length(d$Growth) > 0) {
                # candidate for iguess is detectThresh
                if (detectThresh != 0) {
                  candidate = detectThresh
                } else {
                  # if detectThresh=0, take the first value bigger than 10% of the smallest value (ignoring first entry) as the value to detect tshift
                  tst <- d$Growth[2:length(d$Growth)]
                  candidate = min(tst[tst > 1.1 * min(tst)])
                  # calculate the mean between the candidate and the min value as a iguess
                  candidate = (min(tst) + candidate)/2
                }

            } else {
                candidate = 0
            }
            # take the bigger of candidate and 1e-08
            iguess = max(1e-08, candidate)
        }
    }



    # for Gmp we do not need a tshift
    if (Gmp) {
        tshift = 0
        if (detectThresh != 0) {
            # but want to set all values above detectThresh to =detectThresh
            d$Growth[d$Growth <= detectThresh] = detectThresh
            if (d$Growth[1] > detectThresh) {
                d$Growth[1] = detectThresh
            }
        }
        iguess = tshift  #just because there is an error somewhere downstream if iguess is empty. But not needed for Gmp model
    }

    # Here happens the actual time shift:
    d$Expt.Time = d$Expt.Time - tshift
    # In the unlikely event that the smallest cell density did not occur at the earliest timepoint, need to delete earlier timepoints
    d = d[d$Expt.Time >= 0, ]

    if (Nrm != "ZeroShift") {
      yShift = min(d$Growth)
    } else {
      yShift = 0
    }
    d$Growth = d$Growth - yShift
    #print(min(d$Growth))


    # define bounds for the optim functions here:
    bounds = makeBoundsQFA2(iguess, d, minK, fixG, globalOpt, glog=glog, TimeFormat, Gmp = Gmp, lowK = lowK, upK = upK, lowr = lowr, upr = upr, lowg = lowg, upg = upg, lowv = lowv,
        upv = upv, lowb = lowb, upb = upb, Nrm = Nrm)

    if (Gmp) {
        bounds$b <- bounds$g  #just because I was lazy coding the makeBounds function, just taking the g for b
    }
    inocguess = bounds$inocguess  #copy that over
    # bounds$inocguess=NULL
    xybounds = bounds

    # First *detectable* observation at time t0
    t0 = min(d$Expt.Time) + tshift
    # xybounds$K=c(0.9*max(d$Growth),1.1*max(d$Growth))

    # two options: globaloptim or not:
    if (globalOpt) {
        # actual model fit call:
        pars = de.fit2(d$Expt.Time, d$Growth, inocguess, xybounds, initPop = TRUE, logTransform = logTransform, Gmp = Gmp, Nrm = Nrm, TimeFormat)
        # Check for high fraction of K at end of modelled experiment
        if (!Gmp) {
            GEnd = Glogist(pars[1], pars[2], pars[3], pars[4], max(d$Expt.Time))
            if (GEnd/pars[1] < 0.75) {
                # If experiment not quite finished... print('Modelled growth at end of experiment is less than 0.75 of K estimate...') Put tight bounds on K and optimise again
                Kmin = max(0.95 * 1.5 * GEnd, minK)
                Kmax = max(1.05 * 1.5 * GEnd, minK)
                xybounds$K = c(Kmin, Kmax)
                pars = de.fit2(d$Expt.Time, d$Growth, inocguess, xybounds, initPop = TRUE, logTransform = logTransform, Gmp = Gmp, Nrm = Nrm, TimeFormat = TimeFormat)
            }
        }

    } else {
        # use regular optim function actual model fit call:
        pars = data.fit2(d$Expt.Time, d$Growth, inocguess, xybounds, logTransform = logTransform, TimeFormat = TimeFormat, Gmp = Gmp, Nrm = Nrm)
        # Check for high fraction of K at end of modelled experiment
        if (!Gmp) {
            GEnd = Glogist(pars[1], pars[2], pars[3], pars[4], max(d$Expt.Time))

            if (GEnd/pars[1] < 0.75) {
                # If experiment not quite finished... print('Modelled growth at end of experiment is less than 0.75 of K estimate...') Put tight bounds on K and optimise again
                Kmin = max(0.95 * 1.5 * GEnd, minK)
                Kmax = max(1.05 * 1.5 * GEnd, minK)
                xybounds$K = c(Kmin, Kmax)
                pars = data.fit2(d$Expt.Time, d$Growth, inocguess, xybounds, logTransform = logTransform, TimeFormat = TimeFormat, Gmp = Gmp, Nrm = Nrm)
            }

        }
    }

    # Check for spurious combination of relatively high r, low K, spend more time optimising...
    if (!Gmp) {
        # Try optimising with sick colony as guess:
        if (checkSlow & (mdr(pars[1], pars[2], pars[3], pars[4]) > 1) & ((pars[1] < 0.05) | (max(d$Growth) < 0.05) | (tail(d$Growth, 1) < 0.05))) {
            # print('Attempting to do global fit alternative sick colony growth curve')
            Kmin = max(0.9 * inocguess, minK)
            Kmax = 1.5 * max(pars[1], minK)
            xybounds$K = c(Kmin, Kmax)
            xybounds$r = c(0, 3)  # Slow growth
            if (glog) {
                xybounds$v = c(0.75, 1.5)  # More logistic growth
            } else {
                xybounds$v = c(1, 1)
            }
            inits = list(K = pars[1], r = 0.6, g = inocguess, v = 1)
            newpars = de.fit2(d$Expt.Time, d$Growth, inocguess, xybounds, inits = inits, initPop = TRUE, widenr = FALSE, logTransform = logTransform, Gmp = Gmp, Nrm = Nrm, TimeFormat = TimeFormat)
            if (newpars[5] <= pars[5])
                pars = newpars
        }
        opt = sumsq(pars[1], pars[2], pars[3], pars[4], obsdat$Growth, obsdat$Expt.Time, logTransform = logTransform)  # Use all data (not just data below detection thresh)
        dead = sumsq(inocguess, 0, inocguess, 1, obsdat$Growth, obsdat$Expt.Time, logTransform = logTransform)
        if (dead <= opt) {
            pars = c(inocguess, 0, inocguess, 1, dead)
            names(pars) = c("K", "r", "g", "v", "objval")  # Try dead colony
        }
        if (!is.finite(pars[["K"]]))
            pars[["K"]] = minK
        if (!is.finite(pars[["r"]]))
            pars[["r"]] = 0
        if (!is.finite(pars[["g"]]))
            pars[["g"]] = 0
        if ("v" %in% names(pars))
            if (!is.finite(pars[["v"]]))
                pars[["v"]] = 1
        pars[["tshift"]] = tshift
    }else{
      if (Nrm != "ZeroShift") {
          pars[["yShift"]] = yShift
      } else {
          pars[["yShift"]] = 0
      }
    }




    # if(length(pars)!=5) dput(pars)
    return(pars)
}


#' Internal QFA function called by qfa.fit2. Do not call directly
makeBoundsQFA2 <- function(inocguess, d, minK = 0, fixG = FALSE, globalOpt = FALSE, glog = TRUE, TimeFormat = "d", Gmp = F, lowK = NA, upK = NA, lowr = NA, upr = NA,
    lowg = NA, upg = NA, lowv = NA, upv = NA, lowb = NA, upb = NA, Nrm = "Raw", ...) {
    # Find optimization bounds if not user-defined (values are not NA)


    # first for r
    if (is.na(lowr))
        lowr = 0
    if (Gmp) {
        # different values for Gmp vs Glog/Slog
        if (TimeFormat == "h") {
            if (is.na(upr))
                upr = max(d$Growth) * 10/24
        } else {
            if (is.na(upr))
                upr = max(d$Growth) * 10
        }
    } else {
        # Glog/Slog
        if (TimeFormat == "h") {
            if (is.na(upr))
                upr = 20
        } else {
            if (is.na(upr))
                upr = 200
        }
    }



    # then, v. does only matter for Glog and Slog. Set boundaries to (1,1) if Slog
    if (glog) {
        # lowv<-0.1; upv<-10.0
        if (is.na(lowv))
            lowv = 0.001
        if (is.na(upv))
            upv = 15
    } else {
            lowv = 1
            upv = 1
    }


    # now K: defined by max Growth values
    if (is.na(lowK))
        lowK <- max(0.85 * max(d$Growth), minK)
    if (is.na(upK))
        upK <- 1.15 * (max(d$Growth))

    # then the g (Glog and Slog). Output lowg and upg are used as b for Gmp
    if (!Gmp) {
        # We often fix inoculation density, but users might prefer to infer it from growth curve vary it only by +- 10% of inocguess
        if (fixG) {
            if (is.na(lowg))
                lowg = 0.9 * inocguess
            if (is.na(upg))
                upg = 1.1 * inocguess
        } else {
            # very high variablility. over 20 logs
            if (is.na(lowg))
                lowg = 1e-10 * inocguess
            if (is.na(upg))
                upg = 1e+10 * inocguess
        }
        lowg <- max(0, lowg)
        upg <- min(upK, upg)
    } else {
        # g is now the b=Tlag

        # g is b for Gmperz. Set boundaries as max time and min time
        if (is.na(lowb)) {
            lowg = min(d$Expt.Time)
        } else {
            lowg = lowb
        }
        if (is.na(upb)) {
            upg = max(d$Expt.Time)
        } else {
            upg = upb
        }
        lowv = 0
        upv = 0
    }


    # print(list(K=c(lowK,upK),r=c(lowr,upr),g=c(lowg,upg),v=c(lowv,upv),inocguess=inocguess))

    return(list(K = c(lowK, upK), r = c(lowr, upr), g = c(lowg, upg), v = c(lowv, upv), inocguess = inocguess))
}


#' Internal QFA function called by qfa.fit2. Do not call directly
data.fit2 <- function(tim, growth, inocguess, xybounds, inits = list(), logTransform = FALSE, verbose = FALSE, TimeFormat = "d", Gmp = F, Nrm = "Raw") {
    # regular optim model fit function


    if (length(inits) == 0) {
        # Get initial guess for parameters
        init <- guess2(tim, growth, inocguess, xybounds, TimeFormat = TimeFormat, Gmp = Gmp)
        if (!Gmp) {
            # for Slog, set v to NULL (don't need that parameter), and if g bounds are equal, that as well (Glog and Slog)
            if (xybounds$v[1] == xybounds$v[2])
                init$v = NULL
            if (xybounds$g[1] == xybounds$g[2])
                init$g = NULL
        } else {
            # don't need v in Gmp.
            init$v = NULL
        }
    } else {
        init = inits
    }

    # if the initial r guess is 10% smaller/bigger as the lower/upper boundary, set init r to mean of lower/upper boundaries
    if (init$r < 1.1 * xybounds$r[1] | init$r > 0.9 * xybounds$r[2]) {
        init$r <- mean(xybounds$r)
    }


    # Try to stop L-BFGS-B errors by moving away from minK xybounds$K[1]=1.01*xybounds$K[1] Function to be optimized

    # define parameter scale for the model fit
    if (!Gmp) {
        if (TimeFormat == "h") {
            paramscale = c(0.1, 0.1, 0.0015, 10)
        } else {
            paramscale = c(1, 500, inocguess, 15)
        }
        names(paramscale) = c("K", "r", "g", "v")
    } else {
        paramscale = c(2, 10, 1)
        names(paramscale) = c("K", "r", "b")
    }



    # depending on which model is used and how the boundaries are set, different function parameters for sum of square function are defined
    if (!Gmp) {
        if (xybounds$v[1] == xybounds$v[2]) {
            if (xybounds$g[1] == xybounds$g[2]) {
                objf <- function(modpars) {
                  return(sumsq(modpars[["K"]], modpars[["r"]], xybounds$g[1], xybounds$v[1], growth, tim, logTransform = logTransform))
                }
                lbounds = c(xybounds$K[1], xybounds$r[1])
                ubounds = c(xybounds$K[2], xybounds$r[2])
                pscl = as.numeric(paramscale[c("K", "r")])
            } else {
                objf <- function(modpars) {
                  return(sumsq(modpars[["K"]], modpars[["r"]], modpars[["g"]], xybounds$v[1], growth, tim, logTransform = logTransform))
                }
                lbounds = c(xybounds$K[1], xybounds$r[1], xybounds$g[1])
                ubounds = c(xybounds$K[2], xybounds$r[2], xybounds$g[2])
                pscl = as.numeric(paramscale[c("K", "r", "g")])
            }
        } else {
            if (xybounds$g[1] == xybounds$g[2]) {
                objf <- function(modpars) {
                  return(sumsq(modpars[["K"]], modpars[["r"]], xybounds$g[1], modpars[["v"]], growth, tim, logTransform = logTransform))
                }
                lbounds = c(xybounds$K[1], xybounds$r[1], xybounds$v[1])
                ubounds = c(xybounds$K[2], xybounds$r[2], xybounds$v[2])
                pscl = as.numeric(paramscale[c("K", "r", "v")])
            } else {
                objf <- function(modpars) {
                  return(sumsq(modpars[["K"]], modpars[["r"]], modpars[["g"]], modpars[["v"]], growth, tim, logTransform = logTransform))
                }
                lbounds = c(xybounds$K[1], xybounds$r[1], xybounds$g[1], xybounds$v[1])
                ubounds = c(xybounds$K[2], xybounds$r[2], xybounds$g[2], xybounds$v[2])
                pscl = as.numeric(paramscale[c("K", "r", "g", "v")])
            }
        }
    } else {
        objf <- function(modpars) {
            return(sumsq2(modpars[["K"]], modpars[["r"]], modpars[["b"]], growth, tim, logTransform = logTransform, Gmp = Gmp))
        }
        lbounds = c(xybounds$K[1], xybounds$r[1], xybounds$b[1])
        ubounds = c(xybounds$K[2], xybounds$r[2], xybounds$b[2])
        pscl = as.numeric(paramscale[c("K", "r", "b")])
    }


    # Perform optimization with given parameters
    optsol <- optim(par = unlist(init), fn = objf, gr = NULL, method = "L-BFGS-B", lower = lbounds, upper = ubounds, control = list(maxit = 1e+05, factr = 1e+08, trace = 0,
        parscale = pscl))
    pars = abs(optsol$par)
    objval = objf(pars)


    # a bit of name plumbing and output the values
    pars[["objval"]] = objval
    if (!Gmp) {
        pars = pars[c("K", "r", "g", "v", "objval")]
        if (!"g" %in% names(pars))
            pars[["g"]] = xybounds$g[1]
        if (!"v" %in% names(pars))
            pars[["v"]] = xybounds$v[1]
        # Sanity check for fitted parameters (no negative growth)
        if (pars[["K"]] <= pars[["g"]]) {
            pars[["K"]] = pars[["g"]]
            pars[["r"]] = 0
            pars[["v"]] = 1
            objval = objf(pars)
        }
        pars = pars[c("K", "r", "g", "v", "objval")]
    } else {
        if (objval > 5e-04) {


            if (TimeFormat == "h") {
                newr = init$K/10
            } else {
                newr = (init$K/10) * 24
            }
            # print(newr); print(TimeFormat)
            ubounds = c(xybounds$K[2], newr, xybounds$b[2])
            optsol <- optim(par = unlist(init), fn = objf, gr = NULL, method = "L-BFGS-B", lower = lbounds, upper = ubounds, control = list(maxit = 1e+05, factr = 1e+08,
                trace = 0, parscale = pscl))

            objval2 = objf(pars)
            if (objval2 < objval) {
                pars = abs(optsol$par)
            }

        }
        pars = pars[c("K", "r", "b", "objval")]


        pars[5] <- pars[1] * exp(-exp(exp(1) * pars[2] * pars[3]))
        names(pars)[5] <- "g"
        if (pars[["K"]] <= pars[["g"]]) {
            pars[["K"]] = pars[["g"]]
            pars[["r"]] = 0
            pars[["b"]] = NA
        }

    }

    if (verbose) {
        if (optsol$message != "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH")
            print(optsol$message)
    }
    return(pars)
}


#' Internal QFA function called by qfa.fit2. Do not call directly
de.fit2 <- function(tim, growth, inocguess, xybounds, inits = list(), initPop = FALSE, widenr = F, logTransform = FALSE, mxit = 2000, Gmp = F, Nrm="Raw", TimeFormat = TimeFormat) {
    ### Function that fits model to a timecourse with genetic optimisation (DEoptim) algorithm

    # Fit to growth curve with differential evolution Get initial guess for parameters
    if (length(inits) == 0) {
        # Get initial guess for parameters
        init <- guess2(tim, growth, inocguess, xybounds, TimeFormat = TimeFormat, Gmp = Gmp)
        if (!is.finite(init$r) | is.na(init$r)) {
            init$r = max(xybounds$K)/max(tim)
        }
    } else {
        init = inits
    }



    # Function to be optimized (minimze sum of squares)
    if (!Gmp) {
        if (widenr)
            xybounds$r = c(0.01 * init$r, 10 * init$r)
        objf <- function(modpars) {
            K <- modpars[1]
            r <- modpars[2]
            g <- modpars[3]
            v <- abs(modpars[4])
            res = sumsq(K, r, g, v, growth, tim, logTransform = logTransform)
            return(res)
        }

        low = c(xybounds$K[1], xybounds$r[1], xybounds$g[1], xybounds$v[1])
        up = c(xybounds$K[2], xybounds$r[2], xybounds$g[2], xybounds$v[2])


    } else {
        objf <- function(modpars) {
            K <- modpars[1]
            r <- modpars[2]
            b <- modpars[3]
            res = sumsq2(K, r, b, growth, tim, logTransform = logTransform, Gmp = Gmp, yShift = 0)
            return(res)
        }

        low = c(xybounds$K[1], xybounds$r[1], xybounds$b[1])
        up = c(xybounds$K[2], xybounds$r[2], xybounds$b[2])


        if (low[1] >= up[1]) low[1] = 0

    }



    # Format bounds
    NumParticles = 11 * length(low)

    # as a default, sample uniformly, random in the given boundaries as an initial population
    if (initPop) {
        # Randomly sample from within bounds for initial population
        if (!Gmp) {
            Klist = runif(NumParticles, min = low[1], max = up[1])
            rlist = runif(NumParticles, min = low[2], max = up[2])
            glist = runif(NumParticles, min = low[3], max = up[3])
            vlist = runif(NumParticles, min = low[4], max = up[4])
            pop = data.frame(K = Klist, r = rlist, g = glist, v = vlist)
            pop = as.matrix(pop)
            # Set first particle equal to initial guess
            pop[1, ] = as.numeric(t(init))
            # Set second particle equal to a dead colony
            pop[2, ] = c(min(growth), 0, low[3], 1)
        } else {
            Klist = runif(NumParticles, min = low[1], max = up[1])
            rlist = runif(NumParticles, min = low[2], max = up[2])
            blist = runif(NumParticles, min = low[3], max = up[3])
            pop = data.frame(K = Klist, r = rlist, b = blist)
            pop = as.matrix(pop)
            # Set first particle equal to initial guess
            pop[1, ] = as.numeric(t(init))
            # Set second particle equal to a dead colony
            pop[2, ] = c(min(growth), 0, up[3])
        }

    } else {
        pop = NULL
    }

    # this is the optimization call
    optsol = DEoptim::DEoptim(objf, lower = low, upper = up, DEoptim::DEoptim.control(trace = F, NP = NumParticles, initialpop = pop, itermax = mxit))

    # a bit of name plumbing and finally output the data
    pars = abs(as.numeric(optsol$optim$bestmem))
    objval = objf(pars)
    if (!Gmp) {
        # Sanity check for fitted parameters (no negative growth)
        if (pars[1] < pars[3]) {
            pars[1] = pars[3]
            pars[2] = 0
            pars[4] = 1
            objval = objf(pars)
        }
        pars = c(pars, objval)
        names(pars) = c("K", "r", "g", "v", "objval")
    } else {
        pars = c(pars, objval)
        names(pars) = c("K", "r", "b", "objval")

        pars[5] <- pars[1] * exp(-exp(exp(1) * pars[2] * pars[3]))
        names(pars)[5] <- "g"
        if (pars[["K"]] <= pars[["g"]]) {
            pars[["K"]] = pars[["g"]]
            pars[["r"]] = 0
            pars[["b"]] = NA
        }
    }

    return(pars)
}


#' Internal QFA function called by qfa.fit2. Do not call directly
guess2 <- function(tim, growth, inocguess, xybounds, minK = 0.025, TimeFormat = "d", Gmp = F) {
    # Sort time and growth
    growth = growth[order(tim)]
    tim = tim[order(tim)]
    n = sum((!is.na(tim)) & (!is.na(growth)))


    if (!Gmp) {
        # first the guesses for Glog and Slog model parameters define min g via as inocguess or if not specified by user via min growth value
        if (is.null(inocguess)) {
            G0g <- max(min(growth, na.rm = TRUE), 1e-07)
        } else {
            G0g <- inocguess
        }

        # guess for K depending on Growth values
        Kg <- max(max(growth, na.rm = TRUE), minK)
        vg = 0.2  # Start lower than logistic (which is v=1), experience showed that v is often around 0.2 for Colonyzer data

        rg = 0
        if (n > 3) {
            # Guess for r is a bit more complicated: calculate it from a smoothed approximation of the growth curve etc.
            targ <- min(growth) + (max(growth) - min(growth))/2
            approxspline = smooth.spline(tim, growth)
            approxcurve = approxfun(approxspline$x, approxspline$y, rule = 2)
            solvefn <- function(x) approxcurve(x) - targ
            if ((approxcurve(max(tim)) >= approxcurve(min(tim))) & (sign(solvefn(min(tim))) != sign(solvefn(tim[which.max(growth)])))) {
                sol = uniroot(solvefn, c(min(tim), tim[which.max(growth)]))
                tmrate = sol$root
            } else {
                # Too few points to fit spline curve, still do something sensible...
                tmrate = (min(tim) + max(tim))/2
            }
            if (Kg > G0g) {rg = log((Kg - G0g)/G0g)/tmrate}
        }
        # Sanity check for guessed parameter values If the data have low correlation, then set r=0 If all elements of growth are equal (e.g. zero) then correlation function
        # throws error...
        if ((length(unique(growth)) == 1) | (cor(tim, growth) < 0.1)) {
            rg = 0
            Kg = G0g
        }
        if (rg <= 0) {
            if (TimeFormat == "h") {
                rg = 2.5
            } else {
                rg = 20
            }
        }
        if (TimeFormat == "h" & rg > 5) {
            rg = 2.5
        }
        return(list(K = Kg, r = rg, g = G0g, v = vg))


    } else {
        # estimation for Gmp model K as max growth, b as time where Growth is 36.79% of max Growth (K). (model definition of where b is)
        Kg = max(growth)
        bg = tim[which(growth >= 0.3679 * max(growth), 1)[1]]

        # r is similarly estimated as for Glog/Slog
        targ <- min(growth) + (max(growth) - min(growth))/2
        approxspline = smooth.spline(tim, growth)
        approxcurve = approxfun(approxspline$x, approxspline$y, rule = 2)
        solvefn <- function(x) approxcurve(x) - targ
        if ((approxcurve(max(tim)) >= approxcurve(min(tim))) & (sign(solvefn(min(tim))) != sign(solvefn(tim[which.max(growth)])))) {
            sol = uniroot(solvefn, c(min(tim), tim[which.max(growth)]))
            tmrate = sol$root
        } else {
            # Too few points to fit spline curve, still do something sensible...
            tmrate = (min(tim) + max(tim))/2
        }
        if (Kg > min(growth)){
          rg = log((Kg - min(growth[growth > 0]))/min(growth[growth > 0])+0.00001)/tmrate
        }else{
          rg=0
        }
        return(list(K = Kg, r = rg, b = bg))

    }
}  #guess


#' Internal QFA function called by qfa.fit2. Do not call directly
numerical_r2 = function(obsdat, mkPlots = F, span = 0.3, nBrute = 10000, cDiffDelta = 1e-04, mlab = "") {
    # Generate numerical (model-free) estimate for nr_t intrinsic growth rate
    tims = obsdat$Expt.Time
    gdat = obsdat$Growth
    tmax = max(tims)

    # Smooth data, find slope as function of time and nr_t slope
    lgdat = log(gdat + abs(min(gdat[gdat != 0])))

    ltims = tims[!is.na(lgdat)]
    ltims = tims[!is.infinite(lgdat)]
    lgdat = lgdat[!is.na(lgdat)]
    lgdat = lgdat[!is.infinite(lgdat)]
    la = NA

    try(a <- loapproxfun(tims, gdat, span = span), silent = F)
    try(la <- loapproxfun(ltims, lgdat, span = span), silent = F)
    problem = list(nr = 0, nr_t = NA, mslp = 0, mslp_t = NA)

    if (!exists("a"))
        return(problem)
    if (!is.function(a))
        return(problem)
    if (!is.function(la))
        return(problem)

    centralDiff = function(f, delta) return(function(x) (f(x + delta/2) - f(x - delta/2))/delta)
    lslp = centralDiff(la, cDiffDelta)
    slp = centralDiff(a, cDiffDelta)
    # Brute force optimization
    stimes = seq(min(ltims), max(ltims), length.out = nBrute)
    vals = a(stimes)
    slps = slp(stimes)
    lvals = la(stimes)
    lslps = lslp(stimes)

    # Discard points too close to t=0 to avoid artificially high slopes
    opt = which.max(slps)
    lopt = which.max(lslps)
    res = list(nr = lslps[lopt], nr_t = stimes[lopt], mslp = slps[opt], mslp_t = stimes[opt])
    # maxlslp=optimize(lslp,lower=min(ltims),upper=max(ltims),nr_t=TRUE) res$nr=maxlslp$objective res$nr_t=maxlslp$maximum
    maxslope = res$nr


    # and here is a part of the function where the vizualisation of this part can be possible. I don't even know how to activate that. Old part, unused...
    if (mkPlots) {
        mainlab = paste("Max. intrinsic growth rate =", formatC(res$nr, 3), "(estimated on log scale)", mlab)
        # Plot synthetic data, Loess approximation and estimated slope
        op = par(mar = c(5, 4, 4, 5) + 0.1, mfrow = c(1, 2))
        plot(NULL, xlab = "Time (d)", ylab = "", xlim = c(-0.1 * tmax, tmax), ylim = range(lgdat), main = mainlab)
        abline(v = res$nr_t, lwd = 3, col = "green", lty = 2)
        slope = res$nr
        intercept = la(res$nr_t) - slope * res$nr_t
        abline(a = intercept, b = slope, lwd = 3, col = "green")
        points(ltims, lgdat)
        mtext("Log Cell density (AU)", side = 2, line = 3, col = "red")
        curve(la, from = -0.1 * tmax, to = tmax, add = TRUE, col = "red", lwd = 2, xlim = c(-0.1 * tmax, tmax))
        par(new = TRUE)

        curve(lslp, from = -0.1 * tmax, to = tmax, col = "blue", xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = 2, xlim = c(-0.1 * tmax, tmax))
        axis(4)
        mtext("Slope of Log Cell Density", side = 4, line = 3, col = "blue")

        mainlab = paste("Max. slope =", formatC(res$mslp, 3), "(estimated on linear scale)", mlab)
        # Plot synthetic data, Loess approximation and estimated slope
        op = par(mar = c(5, 4, 4, 5) + 0.1)
        plot(NULL, xlab = "Time (d)", ylab = "", xlim = c(-0.1 * tmax, tmax), ylim = range(obsdat$Growth), main = mainlab)
        abline(v = res$mslp_t, lwd = 3, col = "green", lty = 2)
        slope = slp(res$mslp_t)
        intercept = a(res$mslp_t) - slope * res$mslp_t
        abline(a = intercept, b = slope, lwd = 3, col = "green", untf = FALSE)
        points(obsdat$Expt.Time, obsdat$Growth)
        mtext("Cell density (AU)", side = 2, line = 3, col = "red")
        curve(a, from = -0.1 * tmax, to = tmax, add = TRUE, col = "red", lwd = 2, xlim = c(-0.1 * tmax, tmax))
        par(new = TRUE)

        curve(slp, from = -0.1 * tmax, to = tmax, col = "blue", xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = 2, xlim = c(-0.1 * tmax, tmax))
        axis(4)
        mtext("Slope of Cell Density", side = 4, line = 3, col = "blue")

        par(op)
    }
    return(res)
}

#' Internal QFA function called by qfa.fit2. Do not call directly
sumsq2 <- function(K, r, b, growth, tim, logTransform = FALSE, Gmp = F, yShift = 0) {
    # this is the sumsq function for Gmp. For Glog and Slog, sumsq is called

    ss = sqrt(sum((growth - Gmprtz(K, r, b, tim) + yShift)^2))/length(growth)
    if (is.na(ss)) {
        return(Inf)
    } else {
        return(ss)
    }
}


#' Internal QFA function. Do not call directly
#'
#' Closure returning a loess-based approximating function for a timeseries
loapproxfun = function(t, g, span = 0.2) {
  lo = loess(g ~ t, span = span)
  # If loess smoothing with specified span does not work (e.g. too few points), revert to linear interpolation
  if (is.na(lo$fitted[1])) {
    loex = approxfun(t, g, method = "linear", rule = 2)
  } else {
    loex = function(t) {
      cdens = predict(lo, t)
      cdens[t < min(lo$x)] = predict(lo, min(lo$x))
      cdens[t > max(lo$x)] = predict(lo, max(lo$x))
      return(cdens)
    }
  }
  return(loex)
}

# Sum of squared error
sumsq <- function(K, r, g, v, growth, tim, logTransform = FALSE) {
  # For generalised logistic function, K/g must be positive to avoid complex cell density estimates
  if ((K/g) < 0)
    return(Inf)
  if (logTransform) {
    glog = growth[growth > 0]
    tlog = tim[growth > 0]
    ss = sqrt(sum((log(glog) - log(Glogist(K, r, g, v, tlog)))^2))/length(glog)
  } else {
    ss = sqrt(sum((growth - Glogist(K, r, g, v, tim))^2))/length(growth)
  }
  if (is.na(ss)) {
    return(Inf)
  } else {
    return(ss)
  }
}


# Extract row & col from list of position vectors
rcget <- function(posvec, rc) posvec[match(rc, c("row", "col"))]

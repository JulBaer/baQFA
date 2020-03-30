#General use functions----

#' Internal QFA function. Do not call directly
#'
#' Prevent zero or negative growth from occurring
nozero<-function(growth){
  if (!is.na(growth)) {
    if (growth<=0){growth<-0.0000000000001} else {growth<-growth}
  }else{
    growth=NA
  }
  }


#' Internal QFA function. Do not call directly
#'
#' Remnant of old coding, does nothing
na2zero <- function(num) {
  return(num)
}


#' Internal QFA function. Do not call directly
#'
#' Create inoculation time environment
bctimes <- function(bctimef) {
  startframe = read.delim(bctimef, colClasses = c("character"))
  starts <- new.env(hash = TRUE)
  for (k in 1:length(startframe$Barcode)) {
    st = startframe$Start.Time[k]
    assign(startframe$Barcode[k], st, envir = starts)
  }
  return(starts)
}  # bctimes


#' Internal QFA function. Do not call directly
#'
#' Convert barcode to start time
bc2st <- function(bc, inocenv) get(as.character(bc), envir = inocenv)


#' Internal QFA function. Do not call directly
#'
#' Create ORF2Gene dictionary
orf2gdict <- function(dictionary) {
  cache = new.env(hash = TRUE)
  # Read in dictionary txt file
  orf2g = read.delim(dictionary, header = FALSE, colClasses = c("character"), col.names = c("ORF", "Gene"))
  orf2g$ORF = toupper(orf2g$ORF)
  z = apply(orf2g, 1, orfun, cache)
  return(cache)
}

#' Internal QFA function. Do not call directly
#'
#' orf2gdict subfunction
orfun <- function(row, environ) {
  orf <- as.character(row[1])
  gene <- as.character(row[2])
  return(assign(orf, gene, envir = environ))
}  #orfun

#' Internal QFA function. Do not call directly
#'
#' orf to gene translation
orf2g <- function(orf, dictenv) {
  get(orf, envir = dictenv)
}

#' Internal QFA function. Do not call directly
#'
#' standard error of mean
sterr <- function(x) sqrt(var(x)/length(x))

#' Internal QFA function. Do not call directly
#'
#' custom aggregate call
summarise <- function(df, fdef = "fit", summarystat = mean) {
  summ = aggregate(formula(paste(fdef, "ORF", sep = "~")), df, summarystat)
  return(setNames(summ[[fdef]], summ$ORF))
}

#' Internal QFA function. Do not call directly
#'
#' Grabs posteriors for each variable for an ORF
posmake <- function(orf, orfs, varpos) {
  orfn <- match(orf, orfs)
  norfs <- length(orfs)
  return(lapply(varpos, varposget, orfn, norfs))
}


#' Internal QFA function. Do not call directly
#'
#' Returns posterior if variable repeated for all ORFs
varposget <- function(var, orfn, norfs) {
  if (!is.null(dim(var))) {
    z <- var[, orfn]
  } else {
    z <- NULL
  }
  return(z)
}

#' Internal QFA function. Do not call directly
#'
#' Extracting meta data of a position from the barcode metadata
colony.info <- function(position, bcdata) {
    # Get row & column to restrict data
    row <- position[1]
    col <- position[2]
    d <- bcdata[(bcdata$Row == row) & (bcdata$Column == col), ]
    # Want to use the LAST datapoint Most relevant for edge length (possibly most reliable for position)
    d <- d[order(d$Date.Time, decreasing = TRUE), ]
    d <- d[1, ]
    # Character vector of colony info
    cvec = list(Barcode = as.character(d$Barcode), Treatments = as.character(d$Treatments), Medium = as.character(d$Medium), ORF = as.character(d$ORF), Screen.Name = as.character(d$Screen.Name),
        Library.Name = as.character(d$Library.Name), MasterPlate.Number = as.numeric(d$MasterPlate.Number), Timeseries.order = as.numeric(d$Timeseries.order), ScreenID = as.character(d$ScreenID),
        Tile.Dimensions.X = as.numeric(d$Tile.Dimensions.X), Tile.Dimensions.Y = as.numeric(d$Tile.Dimensions.Y), X.Offset = as.numeric(d$X.Offset), Y.Offset = as.numeric(d$Y.Offset),
        Threshold = as.numeric(d$Threshold), Edge.length = as.numeric(d$Edge.length), Edge.Pixels = as.numeric(d$Edge.Pixels), RepQuad = as.numeric(d$RepQuad))


    if ("Client" %in% colnames(d)) {
        cvec$Client = as.character(d$Client)
    }
    if ("ExptDate" %in% colnames(d)) {
        cvec$ExptDate = as.character(d$ExptDate)
    }
    if ("User" %in% colnames(d)) {
        cvec$User = as.character(d$User)
    }
    if ("PI" %in% colnames(d)) {
        cvec$PI = as.character(d$PI)
    }
    if ("Condition" %in% colnames(d)) {
        cvec$Condition = as.character(d$Condition)
    }
    if ("Inoc" %in% colnames(d)) {
        cvec$Inoc = as.character(d$Inoc)
    }
    if ("Gene" %in% colnames(d)) {
        cvec$Gene = as.character(d$Gene)
    }

    if ("QFAtype" %in% names(d)) {
        if ("QFAtype" == "liquid") {
            ind1 <- which(names(control) == "Expt.Time")
            Variables = names(control)[1:ind1]

        }
    }
    return(cvec)
}


#'Internal QFA function. Do not call directly.
#'
#'calculates timedifference in days from input string
tconv <- function(tstring, startt, fmt) {
  t <- as.POSIXlt(as.character(tstring), format = fmt)
  return(as.numeric(difftime(t, startt, units = "days")))
}

#'Internal QFA function. Do not call directly.
#'
#'Timing for time to calculate stuff
trep <- function(thing, tobj) {
  print(paste("Finished", thing, "in", round(tobj[3], 2), "seconds"))
  print(paste(round(tobj[3]/60, 2), "minutes"))
}



#Model functions and small calculation functions (fitness etc)----

#'Internal QFA function. Do not call directly.
#'
#' Called by qfa.fit2 for logistic model calls. Do not call directly
Glogist <- function(K, r, g, v, t) {
    K/(1 + (-1 + (K/g)^v) * exp(-r * v * t))^(1/v)
}

#' Internal QFA function. Do not call directly
#'
#' Maximum slope of generalised logistic growth curve (on linear scale)
gen_logistic_maxslp <- function(K, r, g, v) r * v * K/(v + 1)^((v + 1)/v)

#'Internal QFA function. Do not call directly.
#'
#' Called by qfa.fit2 for gompertz model calculation
Gmprtz <- function(K, r, b, t) {
    K * exp(-(exp(-(exp(1) * r * (t - b))/K)))
}  #the absolute growth rate version

#'Internal QFA function. Do not call directly.
#'
# 'Calculates the derivative of a generalised logistic model
GlogistDer <- function(K, r, g, v, t) {
    a = K * r * ((K/g)^v - 1)
    b = -1 * r * v
    c1 = (K/g)^v - 1
    d = -1 * r
    e = (-1 * (v + 1)/v)
    a * exp(b * t) * (c1 * exp(d * t) + 1)^e
    K * r * ((K/g)^v - 1) * exp(-1 * r * t * v) * (((K/g)^v - 1) * exp(-1 * r * t) + 1)^(-1 * (v + 1)/v)
}

#'Internal QFA function. Do not call directly.
#'
#' Called by qfa.fit2 for MDR calculation.
mdr <- function(K, r, g, v, Gmp = F) {
  if (!Gmp) {
        sapply((r * v)/log(1 - (2^v - 1)/((2^v) * (g/K)^v - 1)), na2zero)
    } else {
        for (i in 1:length(g)) {
            if (g[i] <= 0) {
                g[i] <- 1e-25
            }
        }
        sapply(r/log((2 * (K - g))/(K - 2 * g)), na2zero)
    }
}

#'Internal QFA function. Do not call directly.
#'
#' Called by qfa.fit2 for MDP calculation.
mdp <- function(K, r, g = 0, v = 1, Gmp = F) {
    for (i in 1:length(g)) {
        if (g[i] <= 0) {
            g[i] <- 1e-25
        }
    }
    sapply(log(K/g)/log(2), na2zero)
}

#'Internal QFA function. Do not call directly.
#'
#' Calculate mdr*mdp
mdrmdp <- function(K, r, g = 0, v = 1, Gmp = F) mdr(K, r, g, v) * mdp(K, r, g = 0, v = 1)

#'Internal QFA function. Do not call directly.
#'
#'Doubling time for generalised logistic function, as a function of time t
dtl <- function(K, r, g, v, t) sapply((-(r * t * v) + log((2^v * (1 - (K/g)^v) * ((1 + (-1 + (K/g)^v)/exp(r * t * v))^(-1/v))^v)
                                                          /(-1 + 2^v * ((1 + (-1 + (K/g)^v)/exp(r * t * v))^(-1/v))^v)))/(r * v), na2zero)

#'Internal QFA function. Do not call directly.
#'
#'Logistic growth model calculation
logist <- function(K, r, g, t) (K * g * exp(r * t))/(K + g * (exp(r * t) - 1))




#

#check if plotting packages are installed----

#'Internal QFA function. Do not call directly.
#'
#' Internal QFA function to check if pachages for plotting are installed.
  qfa_check_packages <- function() {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

    if (packageVersion("ggplot2") < 3) {
      stop("Package \"ggplot2\" version 3.0.0 or higher for this function to work. Please update it.")
    }

    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Package \"gridExtra\" needed for this function to work. Please install it.",
           call. = FALSE)
    }


    if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package \"cowplot\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("Package \"RColorBrewer\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }




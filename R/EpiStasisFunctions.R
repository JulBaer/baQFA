################################################## Epistasis Function ###########################################################

#' Finds genetic interaction strengths and p-values
#'
#' This function is from the original QFA v1 release and was not tested by the authors
#' of the re-release of the bacterial adaptation v2. Use with caution.
#' Fits a genetic independence model between control strains and double mutant
#' strains, either using rjags and a Bayesian linear regression model, or lm
#' and maximum likelihood. For each ORF, the probability that it is a false
#' discovery of a suppressor or enhancer is calculated. These probabilities are
#' then fdr corrected and returned along with genetic interaction scores.
#'
#'
#' @param double Either a qfa.posterior or the results of qfa.fit for the
#' double mutants
#' @param control Either a qfa.posterior or the results of qfa.fit for the
#' control strains
#' @param qthresh The FDR corrected cut off
#' @param orfdict Location of file giving a column of ORFs first and a column
#' of corresponding gene names second - so gene names can be plotted
#' @param GISthresh When returning interaction hitlists, this variable
#' determines the cutoff for strength of genetic interaction.
#' @param plot If TRUE, then a 2-way fitness plot is made.
#' @param modcheck If TRUE then diagnostic residual plots are output to
#' \dQuote{ModelCheck.pdf}
#' @param wctest If TRUE, then use the Wilcoxon test for differences in medians
#' as a measure of statistical significance of genetic interaction.  This is
#' the default.  If FALSE, then use a t-test for difference in mean fitnesses
#' instead.
#' @param bootstrap If TRUE, then use bootstrapping procedure to check if
#' genetic interactions are significant.  If false, then use linear regression
#' and t-test or wilcoxon test.
#' @param Nboot Number of bootstrap samples to generate if using bootstrapping
#' procedure
#' @param subSamp Number of subsamples of available replicates to sample when
#' bootstrapping (default, Inf, uses all available replicates, i.e. each
#' summary (each bootstrap sample) is based on sampling subSamp from N with
#' replacement.  If subSamp==Inf, then subSamp is set equal to N.
#' @param reg String specifying what type of regression to use.  Default is
#' least squares regression as implemented in lm function: "lmreg".
#' Alternatives include "quantreg", "splitreg" and "perpreg".  See lm.epi
#' function help for further details.
#' @param fdef String specifying what fitness definition to use.  Must be the
#' name of a column common to double and control.  Typical options include:
#' "nAUC", "r", "MDRMDP".  The default "fit" is included for backwards
#' compatibility with earlier versions of this function which relied on users
#' manually creating a "fit" column that includes their required fitness
#' definition values.  This was usually achieved by copying an existing column
#' (e.g. "MDRMDP").
#' @return Returns an R list containing three data frames: Results, Enhancers
#' and Suppressors.  Each data frame has the following columns: \itemize{ \item
#' ORF - Unique strain genotype identifier (e.g. Y-number for yeast strains)
#' \item Gene - Human readable genotype identifier \item P - p-value for
#' significance of difference between control and query strain fitnesses \item
#' Q - q-value for significance of difference between control and query strain
#' fitnesses.  This is FDR corrected p-value \item GIS - Genetic interaction
#' strength.  Deviation of (mean or median, depending on value of wctest)
#' observed query strain fitness from expected fitness given control query
#' strain fitness and a multiplicative model of genetic interaction.  \item
#' QueryFitnessSummary - Summary statistic for all available replicate
#' observations of query strain fitness (mean or median, depending on value of
#' wctest).  \item ControlFitnessSummary - Summary statistic for all available
#' replicate observations of control strain fitness (mean or median, depending
#' on value of wctest).  \item QuerySE - Standard error on mean of query strain
#' fitness observations \item ControlSE - Standard error on mean of control
#' strain fitness observations \item TestType - Type of statistical test for
#' significant difference carried out (i.e. Wilcoxon or t-test) \item
#' SummaryType - Type of summary statistic used for fitnesses (i.e. mean or
#' median) \item cTreat - Treatment applied to control plates \item cMed -
#' Medium added to agar in control plates \item cBack - Control plate
#' background tag (experiment identifier) \item qTreat - Treatment applied to
#' query plates \item qMed - Medium added to agar in query plates \item qBack -
#' Query plate background tag (experiment identifier) \item Type - Type of
#' genetic interaction observed (suppressor, enhancer, positive, negative).
#' This is assigned for strains with abs(GIS)>GISthresh and by comparing
#' q-value with qthresh. }
#' @keywords qfa
qfa.epi <- function(double, control, fdef = "fit", qthresh = 0.05, orfdict = "ORF2GENE.txt", GISthresh = 0, plot = TRUE, modcheck = FALSE, wctest = TRUE, bootstrap = NULL, Nboot = 5000,
                    subSamp = Inf, reg = "lmreg") {
  ###### Get ORF median fitnesses for control & double #######
  print("Calculating median (or mean) fitness for each ORF")
  ## LIK ## Get orfs in question
  dorfs <- unique(as.character(double$ORF))
  corfs <- unique(as.character(control$ORF))
  orfs <- sort(intersect(dorfs, corfs))
  control = control[control$ORF %in% orfs, ]
  double = double[double$ORF %in% orfs, ]

  # Get lists with fitnesses for each repeat
  cFstats <- split(control[[fdef]], control$ORF)
  dFstats <- split(double[[fdef]], double$ORF)
  # Get means or medians for each ORF
  if (!is.null(bootstrap)) {
    cFms <- sapply(cFstats, bootstrap)
    dFms <- sapply(dFstats, bootstrap)
  } else {
    if (wctest) {
      cFms <- summarise(control, fdef, median)
      dFms <- summarise(double, fdef, median)
    } else {
      cFms <- summarise(control, fdef, mean)
      dFms <- summarise(double, fdef, mean)
    }
  }
  cSe <- summarise(control, fdef, sterr)
  dSe <- summarise(double, fdef, sterr)
  cCount <- summarise(control, fdef, length)
  dCount <- summarise(double, fdef, length)

  cFitDef <- findFit(control)
  dFitDef <- findFit(double)

  diffFrac = function(GIS, direction = `>`) {
    return(length(GIS[direction(GIS, 0)])/length(GIS))
  }

  ### Fit genetic independence model ###
  m <- lm.epi(dFms, cFms, modcheck, reg = reg)
  print(paste("Ratio of background mutant fitness to wildtype fitness =", round(m, 4)))
  ###### Estimate probability of interaction #######
  print("Calculating interaction probabilities")
  if (!is.null(bootstrap)) {
    GISlist <- sapply(orfs, pgis_bootstrap, cFstats, dFstats, sampSumm = bootstrap, Nreps = Nboot, subSamp = subSamp)
    g = apply(GISlist, 2, bootstrap)
    greaterFrac = apply(GISlist, 2, diffFrac, `>`)
    lessFrac = apply(GISlist, 2, diffFrac, `<`)
    p = ifelse(g < 0, 1 - lessFrac, 1 - greaterFrac)
    pg = rbind(p, g)
  } else {
    pg <- sapply(orfs, pgis, m, cFstats, dFstats, wilcoxon = wctest)
  }

  pg <- as.data.frame(t(pg))
  colnames(pg) = c("p", "gis")
  p <- pg$p

  # Adjust for multiple comparisons
  q <- p.adjust(p, "fdr")
  # Get gene names if needed lik
  genes <- as.character(double$Gene[match(orfs, double$ORF)])
  # genes<-genes[genes%in%as.character(control$Gene)] If gene names missing, find them with orf2g Create orf-gene dictionary
  if (is.null(genes)) {
    orfdict <- orf2gdict(orfdict)
    genes <- sapply(orfs, orf2g, orfdict)
  }

  # Get genetic interaction scores
  gis <- as.numeric(pg$gis)
  # Put into data.frame
  nObs = length(p)
  if (wctest) {
    testType = "wilcoxon"
    sumType = "median"
  } else {
    testType = "t-test"
    sumType = "mean"
  }
  results <- data.frame(ORF = orfs, Gene = genes, P = p, Q = q, GIS = gis, QueryFitnessSummary = dFms, ControlFitnessSummary = cFms, QuerySE = dSe, ControlSE = cSe,
                        QueryCount = dCount, ControlCount = cCount, QueryFit = dFitDef, ControlFit = cFitDef, TestType = rep(testType, nObs), SummaryType = rep(sumType, nObs), cTreat = rep(control$Treatment[1],
                                                                                                                                                                                             nObs), cMed = rep(control$Medium[1], nObs), cScrID = rep(control$ScreenID[1], nObs), qTreat = rep(double$Treatment[1], nObs), qMed = rep(double$Medium[1],
                                                                                                                                                                                                                                                                                                                                      nObs), qScrID = rep(double$ScreenID[1], nObs), qGen = rep(double$Screen.Name[1], nObs), cGen = rep(control$Screen.Name[1], nObs), cLib = rep(paste(unique(control$Library.Name),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         collapse = "_"), nObs), qLib = rep(paste(unique(double$Library.Name), collapse = "_"), nObs))

  if ("Client" %in% colnames(control))
    results$cClient = rep(control$Client[1], nObs)
  if ("Client" %in% colnames(double))
    results$qClient = rep(double$Client[1], nObs)
  if ("ExptDate" %in% colnames(control))
    results$cDate = rep(gsub("'", "", control$ExptDate[1]), nObs)
  if ("ExptDate" %in% colnames(double))
    results$qDate = rep(gsub("'", "", double$ExptDate[1]), nObs)
  if ("User" %in% colnames(control))
    results$cUser = rep(control$User[1], nObs)
  if ("User" %in% colnames(double))
    results$qUser = rep(double$User[1], nObs)
  if ("PI" %in% colnames(control))
    results$cPI = rep(control$PI[1], nObs)
  if ("PI" %in% colnames(double))
    results$qPI = rep(double$PI[1], nObs)
  if ("Condition" %in% colnames(control))
    results$cCond = rep(control$Condition[1], nObs)
  if ("Condition" %in% colnames(double))
    results$qCond = rep(double$Condition[1], nObs)
  if ("Inoc" %in% colnames(control))
    results$cInoc = rep(control$Inoc[1], nObs)
  if ("Inoc" %in% colnames(double))
    results$qInoc = rep(double$Inoc[1], nObs)

  # results$Type<-apply(results,1,typemake,m)
  results$Type = ifelse(results$GIS > 0, "S", "E")
  results <- results[order(results$GIS, results$Q, results$Type), ]
  # Get rid of duplicate entries in results
  orflist <- unique(as.character(results$ORF))
  results <- results[match(orflist, results$ORF), ]
  # Plot results
  if (plot == TRUE) {
    qfa.epiplot(results, qthresh, m)
  }
  final <- list(Results = results, Enhancers = gethits(results, qthresh, type = "E", GISthresh = GISthresh), Suppressors = gethits(results, qthresh, type = "S", GISthresh = GISthresh),
                slope = m, reg = reg)
  if (!is.null(bootstrap))
    final$GISsummary = GISlist
  return(final)
}


############### Epistasis Functions ################## Makes epistasis plot for a given fdr level #
#' Makes an epistasis plot from the full results of qfa.epi
#'
#' This function is from the original QFA v1 release and was not tested by the authors
#' of the re-release of the bacterial adaptation v2. Use with caution.
#' Creates a scatterplot of control fitnesses on the x-axis and query fitnesses
#' on the y-axis, with those deemed to be hits (by FDR adjusted p-value)
#' coloured.  Essentially, this function assumes that the experiment consists
#' of a series of paired fitness observations for a collection (typically
#' genome-wide) of deletion mutations either single mutations (x-axis) or the
#' same single mutation in combination with a common background or query
#' mutation (y-axis).  Fitting a linear regression to all observations (forced
#' through the origin) and searching for significant deviations from that
#' regression is equivalent to searching for mutations which show significant
#' deviation from a multiplicative model of genetic interaction.  Genes whose
#' deletions deviate from this model significantly can be said to interact with
#' the query gene.
#'
#'
#' @param results The results of interaction analysis returned by the qfa.epi
#' function.
#' @param qthresh The fdr adjusted cutoff point for determining hits.
#' @param fitratio The ratio of background mutant fitness to wildtype fitness,
#' from the genetic indepdence model. If FALSE, this is estimated from results
#' using linear regression.
#' @param ref.orf ORF for a reference strain (typically wild-type or a
#' surrogate), whose fitness will be marked on the control and query axes of
#' the interaction plot (horizontal and vertical blue lines). HIS3 is the
#' default strain.
#' @param xxlab x axis label
#' @param yylab y axis label
#' @param mmain Plot label
#' @param fmax Maximum fitness range for both x-axis (control axis) and y-axis
#' (query axis).  If 0, axis ranges are automatically chosen to include all
#' data points.
#' @keywords qfa
qfa.epiplot <- function(results, qthresh, fitratio = FALSE, ref.orf = "YOR202W", xxlab = "Control Fitness", yylab = "Query Fitness", mmain = "Epistasis Plot", fmax = 0) {
  enhancers <- gethits(results$Results, qthresh, type = "E")
  suppressors <- gethits(results$Results, qthresh, type = "S")
  others <- results$Results[(!results$Results$ORF %in% enhancers$ORF) & (!results$Results$ORF %in% suppressors$ORF), ]
  if (fmax == 0) {
    # Get plot parameters
    ymax = 1.1 * max(results$Results$QueryFitnessSummary)
    ymin = 0
    xmax = 1.1 * max(results$Results$ControlFitnessSummary)
    xmin = 0
  } else {
    ymin = 0
    ymax = fmax
    xmin = 0
    xmax = fmax
  }
  plot(NULL, type = "n", xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = xxlab, ylab = yylab, main = mmain, col = 8, pch = 19, cex = 0.5)
  # Add line for genetic independence
  if (fitratio != FALSE) {
    abline(0, fitratio, lwd = 2, col = 8)
  } else {
    slope = results$slope
    abline(0, slope, lwd = 2, col = 8)
  }
  # Add 1:1 fitness line
  abline(0, 1, lwd = 2, lty = 4, col = 8)
  # Add reference ORF fitnesses lines
  if (ref.orf != FALSE) {
    reforf <- results$Results[results$Results$ORF == ref.orf, ]
    abline(v = reforf$ControlFitnessSummary, col = "lightblue", lwd = 2)
    abline(h = reforf$QueryFitnessSummary, col = "lightblue", lwd = 2)
  }
  if (length(others$ORF) > 0) {
    # Add points for non-suppressors & non-enhancers
    points(others$ControlFitnessSummary, others$QueryFitnessSummary, col = "grey", cex = 0.5, pch = 19)
    text(others$ControlFitnessSummary, others$QueryFitnessSummary, others$Gene, col = "grey", pos = 4, offset = 0.1, cex = 0.4)
  }
  # Add suppressors & enhancers
  if (length(enhancers$ORF) > 0) {
    points(enhancers$ControlFitnessSummary, enhancers$QueryFitnessSummary, col = "green", pch = 19, cex = 0.5)
    text(enhancers$ControlFitnessSummary, enhancers$QueryFitnessSummary, enhancers$Gene, col = 1, pos = 4, offset = 0.1, cex = 0.4)
  }
  if (length(suppressors$ORF) > 0) {
    points(suppressors$ControlFitnessSummary, suppressors$QueryFitnessSummary, col = "red", pch = 19, cex = 0.5)
    text(suppressors$ControlFitnessSummary, suppressors$QueryFitnessSummary, suppressors$Gene, col = 1, pos = 4, offset = 0.1, cex = 0.4)
  }
  Corr = signif(cor(results$Results$QueryFitnessSummary, results$Results$ControlFitnessSummary), 3)
  Slope = signif(slope, 3)
  legend("topleft", c(paste("Correlation: ", Corr), paste("Slope: ", Slope)), bty = "n")
}

## Extract hits from epistasis results object ##
gethits <- function(results, qthresh, type = "S", all.types = FALSE, GISthresh = 0) {
  results <- results[results$Q < qthresh, ]
  if (all.types == FALSE) {
    results <- results[results$Type == type, ]
  }
  results <- results[abs(results$GIS) >= GISthresh, ]
  results <- results[!is.na(results$ORF), ]
  return(results[order(results$GIS, results$Q, results$Type), ])
}

# Function to calculate fitnesses for repeats
orfstat <- function(orf, fitframe, fitfunct) {
  orfd <- fitframe[fitframe$ORF == orf, ]
  return(fitfunct(orfd$K, orfd$r, orfd$g, orfd$v))
}

# Regression through origin and point that splits dataset in two
splitRegression = function(x, y) {
  # Find line that goes through one of the points (& origin) and has about half of data above and about half of data below line
  fdat = data.frame(x = x, y = y)
  fdat = fdat[(fdat$x > 0) & (fdat$y > 0), ]
  slopes = fdat$y/fdat$x
  pred = slopes %*% t(fdat$x)
  vert = t(pred) - fdat$y
  vert[abs(vert) < 1e-08] = 0
  above = sign(vert)

  fracs = apply(above > 0, 1, sum)/(length(slopes) - 1)
  best = which.min(abs(fracs - 0.5))
  slopes[best]
}
# Regression minimising perpendicular distance between points and line
perpRegression = function(x, y) {

  perpDist = function(x, y, m) {
    sqrt(((x + m * y)/(m^2 + 1) - x)^2 + (m * (x + m * y)/(m^2 + 1) - y)^2)
  }

  obf = function(m) {
    return(sum(perpDist(x, y, m)))
  }

  perpRes = optim(50, obf, lower = 0, upper = Inf)$par

  return(perpRes)
}

# Linear regression genetic ind. model fit
lm.epi <- function(doubles, controls, modcheck = FALSE, reg = "splitreg") {
  if (reg == "lmreg") {
    indmod = lm(doubles ~ 0 + controls)
    m = as.numeric(indmod$coefficients)
  }
  if (reg == "quantreg") {
    indmod = rq(doubles ~ 0 + controls)
    m = as.numeric(indmod$coefficients)
  }
  if (reg == "splitreg") {
    m = splitRegression(controls, doubles)
  }
  if (reg == "perpreg") {
    m = perpRegression(controls, doubles)
  }

  # Check if linear regression OK
  if ((modcheck == TRUE) & (reg %in% c("quantreg", "lmreg"))) {
    pdf("ModelCheck.pdf")
    resids <- indmod$residuals
    hist(resids, xlab = "Fitness Residuals", main = "Histogram of Residuals")
    qqnorm(resids/sd(resids), main = "Normal QQ Plot")
    abline(0, 1, lwd = 2, col = "red")
    dev.off()
  }
  return(m)
}

# Estimates p-value and estimated strength of interaction
pgis <- function(orf, m, cFs, dFs, wilcoxon = TRUE) {
  # If this orf is not present in both lists, return appropriate p,gis
  if ((length(dFs[[orf]]) == 0) | (length(cFs[[orf]]) == 0)) {
    return(c(1, 0))
  }
  # If fitnesses are not unique in one condition (e.g. all dead) and only one repeat in another (e.g. after stripping) This would cause wilcoxon test to fail due to
  # ties, and t.test to fail due to insufficient y?
  ldFS = length(dFs[[orf]])
  lcFS = length(cFs[[orf]])
  ludFS = length(unique(dFs[[orf]]))
  lucFS = length(unique(cFs[[orf]]))
  if (((ludFS == 1) & (ldFS > 1) & (lcFS == 1)) | ((lucFS == 1) & (lcFS > 1) & (ldFS == 1))) {
    return(c(1, median(dFs[[orf]], na.rm = TRUE) - m * median(cFs[[orf]], na.rm = TRUE)))
  }
  if ((ludFS == 1) & (lucFS == 1)) {
    return(c(1, median(dFs[[orf]], na.rm = TRUE) - m * median(cFs[[orf]], na.rm = TRUE)))
  }
  # If both sets of cultures are dead (and fitnesses therefore equal) return appropriate p,gis
  if (sum(dFs[[orf]] == m * cFs[[orf]]) == length(dFs[[orf]] == m * cFs[[orf]])) {
    return(c(1, 0))
  } else {
    if (wilcoxon) {
      valdiff = 0
      p = 1
      # Returns p-value for significance of difference, and estimate of difference between medians
      tryCatch({
        ctest <- wilcox.test(dFs[[orf]], m * cFs[[orf]], alternative = "two.sided", conf.int = TRUE)
        valdiff <- as.numeric(ctest$estimate)
        p <- as.numeric(ctest$p.value)
      }, error = function(e) {
        print(paste("Wilcox test failed for", orf, ", using t-test instead for this gene"))
        ctest <- t.test(dFs[[orf]], m * cFs[[orf]], alternative = "two.sided", conf.int = TRUE)
        valdiff <- as.numeric(ctest$estimate)
        valdiff <<- valdiff[1] - valdiff[2]
        p <<- as.numeric(ctest$p.value)
      })

      return(c(p, valdiff))
    } else {
      # t-test fails if only one element in either list, do one sample test
      if ((ldFS <= 1) | (lcFS <= 1)) {
        ctest <- t.test(dFs[[orf]] - m * cFs[[orf]])
        p <- as.numeric(ctest$p.value)
        valdiff <- as.numeric(ctest$estimate)
        return(c(p, valdiff))
      } else {
        # Returns p-value for significance of difference, and the difference between the means
        ctest <- t.test(dFs[[orf]], m * cFs[[orf]], alternative = "two.sided", conf.int = TRUE)
        p <- as.numeric(ctest$p.value)
        valdiff <- as.numeric(ctest$estimate)
        valdiff <- valdiff[1] - valdiff[2]
        return(c(p, valdiff))
      }
    }
  }
}

# Estimates p-value and estimated strength of interaction
pgis_bootstrap <- function(orf, cFs, dFs, sampSumm = mean, wt = "YOR202W", Nreps = 10000, subSamp = Inf) {
  # Need to come up with a good estimate for a minumum non-zero fitness in order to avoid problems with division by zero later
  fitMin = median(c(unlist(cFs, use.names = FALSE), unlist(dFs, use.names = FALSE)))/200

  # If this orf is not present in both lists, return appropriate p,gis
  if ((length(dFs[[orf]]) == 0) | (length(cFs[[orf]]) == 0)) {
    return(c(1, 0))
  }

  # Bootstrap estimates of fitness summary distribution for relevant genotypes
  obsDoubleMut = bsSamp(dFs[[orf]], Nreps, sampSumm, subSamp)
  arrayMut = bsSamp(cFs[[orf]], Nreps, sampSumm, subSamp)
  backMut = bsSamp(dFs[[wt]], Nreps, sampSumm, subSamp)
  wtMut = bsSamp(cFs[[wt]], Nreps, sampSumm, subSamp)

  # Uncertainty about summary of predicted double mutant fitness
  predDoubleMut = backMut * arrayMut/pmax(wtMut, fitMin)
  GIS = obsDoubleMut - predDoubleMut
  return(GIS)
}

bsSamp = function(A, Nrep = 1e+05, sampSumm = mean, subSamp = Inf) {
  bsreps = replicate(Nrep, sampSumm(sample(as.numeric(A), min(subSamp, length(A)), replace = TRUE)))
  return(bsreps)
}

# Get type of interaction (DEPRECATED)
typemake <- function(row, m) {
  summd <- as.numeric(row["QueryFitnessSummary"])
  summc <- as.numeric(row["ControlFitnessSummary"])
  if (m < 1) {
    if (summd > m * summc) {
      type <- "S"
    } else {
      type <- "E"
    }
  } else {
    if (summd > m * summc) {
      type <- "E"
    } else {
      type <- "S"
    }
  }
  return(type)
}

#' A screen of three Staphylococcus aureus strains arranged in a grid
#'
#' This is the result of a call of Colonyzer.read of an experiment in which
#' three different Staphylococcus aureus strains have been grown in a grid.
#' The three strains are Cowan and JE2 lab strains and a clinical isolate
#' named CI1149.
#'
#' @docType data
#'
#' @usage data(qfa.testdata)
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references unpublished
#'
#'
#' @examples
#' data(qfa.testdata)
#' #Strip non-experimental edge cultures
#' qfa.testdata=qfa.testdata[(qfa.testdata$Row!=1) & (qfa.testdata$Col!=1) & (qfa.testdata$Row!=8) & (qfa.testdata$Col!=12),]
#' # Define which measure of cell density to use
#' qfa.testdata$Growth = qfa.testdata$Intensity
#' GmpFit = qfa.fit2(qfa.testdata, inocguess=NULL, detectThresh=0, globalOpt=F, AUCLim=NA, TimeFormat="h", Model="Gmp")
#' # Construct fitness measures
#' GmpFit = makeFitness2(GmpFit, AUCLim=NA, plotFitness="All", filename="Example_Gmp_fitness.pdf")
#' # Create plot
#' qfa.plot2("Example_Gmp_GrowthCurves.pdf", GmpFit, qfa.testdata, maxt=30)


"qfa.testdata"

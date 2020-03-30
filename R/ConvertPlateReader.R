#' Converting plate reader growth curve output for QFA usage
#'
#' This function transforms the output of plate reader growth curves saved in 'Time' format to a data.frame suitable for calls to qfa.fit2.
#' Together with the growth curve data (GC) as a data.frame, a second data.frame containing metadata information (info) about the wells of the GC
#' data.frame. There are specific formatting conditions that need to be met by the input data.frames. There needs to be a column named 'Time'
#' or 'time' indicating time for the GC data.frame. This column should either be in hours, min or days and this information is specified with
#' 'timeFormat'. The other column names of the GC data.frame should indicate the 'Well' name. Each of these wells need to be represented with one
#' row in the info data.frame. The info data.frame needs to have a calumn called 'Well' with these matching well names. The second mandatory
#' column should be named 'Condition' and describes the content of each well. Here, the user can also specify which wells are blanks to be used
#' as a correction mean. Additionally, the user needs to specify the numbers of rows and cols. The default values are for 96 well formats.
#' The output data.frame contains aggregated data with one row per well and time combination together with all info data.frame data
#' and is ready to be used with qfa.fit2. The output data is stored in columns named 'OD_raw' and 'OD_calib'.
#'
#' Different calibration methods can be used. The default method is 'blank'. A mean of all wells indicated as 'blank' within the
#' info data.frame column 'Condition' for each measurement timepoint is calculated. This mean is then subtracted from all other wells.
#' Another options is 'individual': the min value of each well is subtracted from itself. With 'none', no calibration is calculated.
#'
#' @param GC data.frame. This contains the growth curve data from a plate reader stored in 'Time' format. Column names should represent well names.
#' One column needs to be named 'Time' and indicate the timecourse of the growth curve either in minutes, hours or days. No other columns should be present
#' @param info data.frame. This contains the metadata of the GC data.frame. Two columns are necessary: One named 'Well' corresponding to the column names
#' of the GC data.frame and a second named 'Condition' describing the well contents. All other columns are added to the final output of GC2QFA
#' @param rows numerical. Indicating number of rows of the growth curve plate. Default set to 8 for 96 well formats
#' @param cols numerical. Indicating number of cols of the growth curve plate. Default set to 12 for 96 well formats. Please note that the length of the
#' info data.frame needs to match rows times cols and that the number of columns of the GC data.frame needs to match rows times cols + 1 (for the time column)
#' @param correction String. Indicating the type of correction method used. Either 'blank', 'individual' or 'none'
#' @param timeFormat String. Indicating the format of the 'Time' column of the GC data.frame. Either 'min', 'h' or 'd'
#' @param shiftZero logical. Indicating if each wells growth curve should be normalised to start from 0 by subtracting from each well its minimum
#' @param ExpName String. Necesarry information for qfa.fit2 to describe the experiment (Barcode). Use this as identifier for each plate/run to ensure
#' smooth operation of qfa.fit2 and following functions
#' @param Inoc.Time String. Inoculation time start in specific format 'YYYY-MM-DD hh:mm:ss'. Necesarry for qfa.fit2 to calculate timecourse
#' @param Library.Name String. Necesarry information for qfa.fit2 to describe the experiment
#' @param MasterPlate.Number Numerical. Necesarry information for qfa.fit2 to describe the experiment
#' @param ScreenID String. Necesarry information for qfa.fit2 to describe the experiment
#' @param Screen.Name String. Necesarry information for qfa.fit2 to describe the experiment
#' @param Treatments Numerical or String. Necesarry information for qfa.fit2 to describe the experiment
GC2QFA <- function(GC, info, rows = 8, cols = 12, correction = "blank", timeFormat = "min", ExpName = "Default", Inoc.Time = "2019-01-20 03:21:00",
                   Library.Name = "Default", MasterPlate.Number = 1, ScreenID = "SaureusTest", Screen.Name = "SauresTest", Treatments = 1,
                   shiftZero = F) {
    # This function transforms the output of plate reader growth curves which have been saved in the 'Time' format to a format that can be read by qfa.fit2 input
    # variables: GC: the plate reader output file. Needs to have one column named Time in which the time is indicated in minutes (or whatever Format is specified by
    # timeFormat). Each other column should contain the measurements for one well ordered by the column Time with the column Name Well (eg 'A3') info: metadata info for
    # each well. It number of Rows need to match the number of samples in GC (eg 96 for 96 well plate). It needs to have a column named 'Well' which corresponds to the
    # colnames of the GC file. Minimum requirement is a column named Condition apart from the 'Well' column. If there are blank wells, they can be specified in
    # Condition. The mean of these will be calculated and subtracted from the other wells if correction='blank' correction: either ='blank' for subtracting from each
    # well*time combination the mean of all wells marked as 'blank' in the metadata file column called Condition. ='individual': the min value of each well is subtracted
    # from each timepoint of that well. ='none': no correction done timeFormat: indicates in which format the Time column in the GC file is. either ='min', ='h' or ='d'
    # shiftZero: logical input specifying if the values for each well should be shifted to start at 0 by subtracting from each value of it its minimum Inoc.Time: string
    # indicating the inoculation time start as needed for the qfa.fit2 in the following format: 'yyyy-mm-dd hh:minmin:secsec' ExpName: string input for specifing the
    # column named ExpName needed for qfa.fit2 Library.Name: string input for specifing the column named Library.Name needed for qfa.fit2 ScreenID: string input for
    # specifing the column named ScreenID needed for qfa.fit2 Screen.Name: string input for specifing the column named Screen.Name needed for qfa.fit2
    # MasterPlate.Number: number input for specifing the column named MasterPlate.Number needed for qfa.fit2 Treatments: number input for specifing the column named
    # Treatments needed for qfa.fit


  mta <- info  # meta data for the wells
  # extract the timevector, modify it to be in min
  if (!"Time" %in% names(GC)) {
    if (!"time" %in% names(GC)) {
      stop("Be sure there is a column called 'Time' in the GC data.frame. It should indicate time in specified timeFormat")
    }else{
      names(GC)[names(GC)=="time"] = "Time"
    }
  }

  if (!"Condition" %in% names(mta)) {
    stop("Be sure there is a column called 'Condition' in the info data.frame. It should describe the content of each well")
  }


  if (timeFormat == "h") {
    Time <- GC$Time * 60
  } else if (timeFormat == "d") {
    Time <- Gc$Time * 24 * 60
  } else if (timeFormat == "min") {
    Time <- GC$Time
  } else {
    stop("Please specify timeFormat as either min, h or d.")
  }




    colV <- rep(c(1:cols), rows)  #vector for the cols
    rowV <- rep(1:rows, each = cols)  # and rows

    if (length(colV) != dim(mta)[1]) {
        stop("the dimension of your specified cols and rows input does not match the dimension of the metadata file")
    }

    if (!"Well" %in% names(mta)) {
      stop("Be sure there is a column called 'Well' in the metadata file. Its content should correspond to the column-names of the growth curve dataframe")
    }

    if (sum(!(mta$Well %in% names(GC)[names(GC)!="Time"]))!=0){
      stop("The content of the Well (names) column in the metadata file does not match the column-names ofo the growth curve dataframe")
    }


    # change all factors to character vectors
    for (i in 1:length(mta)) {
        if (is.factor(mta[, i])) {
            mta[, i] <- as.character(mta[, i])
        }
    }


    #Order the metadata by row and column name (aka by well name)
    mta$Row = substring(mta$Well, 1, 1)
    C1 = as.numeric(substring(mta$Well, 2))
    mta$Column = match(C1, LETTERS)
    if (sum(is.na(mta$Column)) != 0) mta$Column = match(C1, letters)
    if (sum(is.na(mta$Column)) != 0) mta$Column = C1
    mta = mta[order(mta$Row,mta$Column),]








    # initialize the dataframe for export
    GC.dt <- data.frame(matrix(ncol = length(names(mta)) + 3, nrow = length(Time) * length(mta$Well)))
    # take column names from meta and add the new info stuff
    names(GC.dt) <- c(names(mta), "Time", "OD_raw", "OD_calib")

    rwN <- 1  # count for which row to add data



    # over all wells, paste the data in correct format
    for (i in 1:length(mta$Well)) {
        indx <- rwN:(rwN + length(Time) - 1)  # the rows to add for current well

        # first gather all the meta data and add these for each row of measurement (length(Time))
        for (i2 in 1:length(mta)) {
            GC.dt[indx, i2] <- rep(mta[i, i2], times = length(Time))
        }

        # then add the time for the current well, OD, Column and Row
        GC.dt$Time[indx] <- Time
        GC.dt$OD_raw[indx] <- GC[, names(GC) == mta$Well[i]]
        GC.dt$Column[indx] <- as.numeric(rep(colV[i], times = length(Time)))
        GC.dt$Row[indx] <- as.numeric(rep(rowV[i], times = length(Time)))

        # increase row indexing
        rwN <- rwN + length(Time)
    }


    # calculate the correction vector

    if (correction == "blank") {
        # take the average of all blank wells at each time, subtract that from each well
        calib <- vector(length = length(Time), mode = "numeric")
        if (sum(GC.dt$Condition == "blank") == 0) {
            stop("No wells are marked as blanks. Please correct that or choose different correction method")
        }
        # calculate the calibration part from the blanks as the mean of all the blanks at each time
        for (i in 1:length(Time)) {
            calib[i] <- mean(GC.dt$OD_raw[GC.dt$Condition == "blank" & GC.dt$Time == Time[i]], na.rm = T)
        }

        # calculate the calibrated OD by subtracting the calib vector
        for (i in 1:length(mta$Well)) {
            GC.dt$OD_calib[GC.dt$Well == mta$Well[i]] <- GC.dt$OD_raw[GC.dt$Well == mta$Well[i]] - calib
        }

    } else if (correction == "individual") {
        # subtract from each well its min value
        for (i in 1:length(mta$Well)) {
            GC.dt$OD_calib[GC.dt$Well == mta$Well[i]] <- GC.dt$OD_raw[GC.dt$Well == mta$Well[i]] - min(GC.dt$OD_raw[GC.dt$Well == mta$Well[i]])
        }

    } else if (correction == "none") {
        # na calibration done
        GC.dt$OD_calib = NA
    } else {
        GC.dt$OD_calib = NA
        warning("An incorrect correction method was specifed. Correction omited, OD_calib output left empty")
    }





    # add the Expt.Time in days for QFA
    GC.dt$Expt.Time <- (GC.dt$Time/60)/24

    # add some arbitrary stuff for QFA
    GC.dt$Barcode <- ExpName
    GC.dt$ORF <- GC.dt$Condition
    GC.dt$Gene <- GC.dt$Condition
    GC.dt$Library.Name <- Library.Name
    GC.dt$MasterPlate.Number = MasterPlate.Number
    GC.dt$ScreenID <- ScreenID
    GC.dt$Screen.Name <- Screen.Name
    GC.dt$Treatments <- Treatments


    # also, we need the time of experiment start and measurement in a very specific format for QFA to work...
    init <- as.POSIXlt(Inoc.Time)
    expDate <- init + 60 * Time
    GC.dt$Inoc.Time <- gsub(":", "-", gsub(" ", "_", as.character(Inoc.Time)))

    # add this to the final output as well. Also do the downshift to zero here
    rwN <- 1
    for (i in 1:length(mta$Well)) {
        indx <- rwN:(rwN + length(Time) - 1)
        GC.dt$Date.Time[indx] <- gsub(":", "-", gsub(" ", "_", as.character(expDate)))
        rwN <- rwN + length(Time)

        if (shiftZero) {
            GC.dt$OD_calib[indx] <- GC.dt$OD_calib[indx] - min(GC.dt$OD_calib[indx])
        }
    }

    # indicate here that this output file comes from a liquid plate reader growth curve
    GC.dt$QFAtype = "liquid"
    GC.dt$Row = as.numeric(GC.dt$Row)
    # done :)
    return(GC.dt)
}




qfa.mergeMeta <- function(dbc, mta1, rows = NA, cols = NA, Barcode=NA) {
    # this merges a qfa.fit2 output file (dbc) with a metadata file (mta1). If the input for qfa.fit2 was from a GC2QFA call, you can use the same metadata file here.
    # Both files need to have a 'Row' and 'Column' column which are corresponding to each other for copying files input variables: dbc: the output of qfa.fit2
    # data.frame. mta1: metadata file. Either needs to contain Row and Column columns that match the QFA.fit2 output or rows and cols need to be specified rows: number
    # of rows of the input files (eg. 8 for 96well plate) cols: number of cols of the input files (eg. 12 for 96well plate) both rows and cols are only necessary if no
    # Row and Col specified in the metadata

  if (length(unique(dbc$Barcode))>1 && is.na(Barcode)){
    stop("There are more than one plate in the dataset. Please indicate for which barcode the metadata is.")
  }

  if (!is.na(Barcode)){
    dbcBck = dbc[dbc$Barcode!=Barcode,]
    dbc = dbc[dbc$Barcode==Barcode,]
    mltp=T
  }
    # first, check if Row and Column are colnames of both input files
    if ((!("Row" %in% names(dbc)) | !("Row" %in% names(mta1)) | !("Column" %in% names(mta1)) | !("Column" %in% names(mta1)))) {
        if ((!is.na(rows) & !is.na(cols))) {
            # and if not and rows and cols are specified, create input variable stop if the dimensions missmatch
            if (cols * rows > dim(mta1)[1]) {
                stop("The dimension of the rows and cols you specfied does not match the dimension of the metadata you specfied. Please correct.")
            }
            colV <- rep(c(1:cols), rows)  #vector for the cols
            rowV <- rep(1:rows, each = cols)  # and rows

            #Order the metadata by row and column name (aka by well name)
            mta1$Row = substring(mta1$Well, 1, 1)
            mta1$Column = as.numeric(substring(mta1$Well, 2))
            mta1 = mta1[order(mta1$Row,mta1$Column),]

            mta1$Column <- colV
            mta1$Row <- rowV
        } else {
            stop("Make sure both Row and Column are specified in both the QFA output and the metadata or specifiy number of rows and columns in the metadata (make sure the meta is ordered)!")
        }
    }

    # get the names of the metadata columns
    vrs <- names(mta1)

    if (is.na(Barcode)){
      for (ind1 in 1:length(vrs)) {
        # rename the colnames if they are the same as the output from the qfa.fit2
        if ((vrs[ind1] %in% names(dbc)) && (vrs[ind1] != "Row" & vrs[ind1] != "Column")) {
          if (vrs[ind1] != "Condition") {
            names(mta1)[names(mta1) == vrs[ind1]] = paste(vrs[ind1], "1", sep = "_")
            print(paste("variable name", vrs[ind1], "renamed to", paste(vrs[ind1], "1", sep = "_"), sep = " "))
          }else{
            mta1 = mta1[ , !(names(mta1) %in% c("Condition"))]
          }

        }
    }
    }else{
      mta1 = mta1[ , !(names(mta1) %in% c("Condition"))]
      nms = names(mta1)
      nms = nms[!(nms %in% c("Row","Column"))]
      dbc = dbc[, !(names(dbc) %in% nms)]
    }


    # and now, we finally can merge it

    barcMetadata = merge.data.frame(dbc, mta1, by.x = c("Row", "Column"), by.y = c("Row", "Column"))

    if (mltp) {
      oldDim = dim(dbcBck)[2]
      newDim = dim(barcMetadata)[2]
      if (oldDim!=newDim){
        oldDim=oldDim+1
      dbcBck[, oldDim:newDim] = NA
      names(dbcBck)[oldDim:newDim] = names(barcMetadata)[oldDim:newDim]
      }

      barcMetadata = rbind(barcMetadata, dbcBck)
    }

    return(barcMetadata)
}

# # ██▄      ▄▄▄▄▄   ▄████  ██   █    ▀▄    ▄ ▄▄▄▄▄▄   ▄███▄   █▄▄▄▄
# # █  █    █     ▀▄ █▀   ▀ █ █  █      █  █ ▀   ▄▄▀   █▀   ▀  █  ▄▀
# # █   █ ▄  ▀▀▀▀▄   █▀▀    █▄▄█ █       ▀█   ▄▀▀   ▄▀ ██▄▄    █▀▀▌
# # █  █   ▀▄▄▄▄▀    █      █  █ ███▄    █    ▀▀▀▀▀▀   █▄   ▄▀ █  █
# # ███▀              █        █     ▀ ▄▀              ▀███▀     █
# # ▀      █                                 ▀

library(tidyverse)
library(ggplot2)
library(drc)


#' importData
#'
#' Importer for the Roche LightCycler exported data
#' Need to specify if these are raw curves or the first derivative export
#'
#' @param file
#' df path
#' @param type
#' 'Raw': Raw data output
#' 'FirstDeriv': First derivative data output
#' @return
#' df containing cleaned Roche Thermal Cycler data.
#' @export
#'
#' @examples
#' importData(file = './firstderiv.txt', type = 'FirstDeriv')
importData <-
  function(filePath = '',
           type = '') {
    # Bit of rough error-checking for the function.
    if (filePath == '') {
      stop('No file path provided for function.')
    }
    if (type == '') {
      stop('No file type provided for function.')
    }
    if (type != 'Raw' && type != 'FirstDeriv') {
      stop('Incorrect file type provided for function. Choose \'Raw\' or \'FirstDeriv\'')
    }
    
    # Read in tab-delimited file and extract values and a single temp column.
    df <- read.delim(file = filePath,
                     sep = "\t")
    df <-
      tibble(select(df, matches("X(?!.)", perl = TRUE)), select(df, contains('Sample')))
    colnames(df)[1] <- 'Temp'
    colnames(df) <- gsub('\\.\\.(.*)', '', colnames(df), perl = TRUE)
  }

#' dfDerivative
#' 
#' x'[i] = (x[i+1] + x[i]) / 2
#' y' at x'[i] = (y[i+1] - y[i]) / (x[i+1] - x[i])
#'
#' @param df df containing the x,y array to have first derivative taken
#' column 'Temp' should be passed
#' @param plot do we plot the results
#'
#' @return df containing the y' values
#' @export
#'
#' @examples
dfDerivative <- function(df, plot = FALSE) {
  deriv <- df
  deriv[paste(colnames(deriv)[1],'.Deriv',sep='')] <- 0
  deriv[paste(colnames(deriv)[2],'.Deriv',sep='')] <- 0
  for(i in 1:nrow(deriv)){
    deriv[i,3] <- (deriv[i+1,1] + deriv[i,1]) / 2
    deriv[i,4] <- (deriv[i+1,2] - deriv[i,2]) / (deriv[i+1,1] + deriv[i,1])
  }
}
  
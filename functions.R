# # ██▄      ▄▄▄▄▄   ▄████  ██   █    ▀▄    ▄ ▄▄▄▄▄▄   ▄███▄   █▄▄▄▄
# # █  █    █     ▀▄ █▀   ▀ █ █  █      █  █ ▀   ▄▄▀   █▀   ▀  █  ▄▀
# # █   █ ▄  ▀▀▀▀▄   █▀▀    █▄▄█ █       ▀█   ▄▀▀   ▄▀ ██▄▄    █▀▀▌
# # █  █   ▀▄▄▄▄▀    █      █  █ ███▄    █    ▀▀▀▀▀▀   █▄   ▄▀ █  █
# # ███▀              █        █     ▀ ▄▀              ▀███▀     █
# # ▀      █                                 ▀
# #
# # DSFalyzer
# #   Software to import, clean, and analyze DSF data from Roche LightCycler by
# #   second derivative Tm calling and isothermal analysis.
# # Michael Ronzetti, 2021

library(tidyverse)
library(ggplot2)
library(drc)
library(spatialEco)
library(readxl)

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
      tibble(dplyr::select(df, matches("X(?!.)", perl = TRUE)), dplyr::select(df, contains('Sample')))
    colnames(df)[1] <- 'Temp'
    colnames(df) <-
      gsub('\\.\\.(.*)', '', colnames(df), perl = TRUE)
    
    return(df)
  }

#' add_rowcol
#' Add row and column to a tidy dataframe (columns are each temperatures, rows are wells/conditions)
#'
#' @param df 
#' @param well_num 
#'
#' @return
#' @export
#'
#' @examples
add_rowcol <- function(df, well_num) {
  if (well_num == 96) {
    col_by_row <-
      expand.grid(row = sprintf('%.2d', 1:8), col = sprintf('%.2d', 1:12)) %>%
      arrange(., row)
  }
  else if (well_num == 384) {
    col_by_row <-
      expand.grid(row = sprintf('%.2d', 1:16),
                  col = sprintf('%.2d', 1:24)) %>%
      arrange(., row)
  }
  message('Row + Column assignments created for ',
          well_num,
          '-well plate')
  df <- cbind(col_by_row, df)
  return(df)
}

#' well_assignment
#' Add well assignmnets for each plate
#' 
#' @param df 
#' @param well_num 
#'
#' @return
#' @export
#'
#' @examples
well_assignment <- function(df, well_num) {
  if (well_num == 96) {
    letter <- LETTERS[1:8]
    number <- c(1:12)
    number <- str_pad(number, 2, pad = '0')
    tracker <- 1
    temp_df <- tibble(well = c(1:384))
    for (val in letter) {
      for (num in number) {
        temp_df$well[tracker] <- paste(val, num, sep = '')
        tracker <- tracker + 1
      }
    }
  }
  else if (well_num == 384) {
    letter <- LETTERS[1:16]
    number <- c(1:24)
    number <- str_pad(number, 2, pad = '0')
    tracker <- 1
    temp_df <- tibble(well = c(1:384))
    for (val in letter) {
      for (num in number) {
        temp_df$well[tracker] <- paste(val, num, sep = '')
        tracker <- tracker + 1
      }
    }
  }
  message('Well assignments created for ', well_num, '-well plate.')
  df <- cbind(temp_df, df)
  return(df)
}

#' plate_assignment
#' Assign compound ids and concentration from platemap
#'
#' @param df 
#' @param platemap_file 
#'
#' @return
#' @export
#'
#' @examples
plate_assignment <- function(df, platemap_file) {
  id_df <- read_excel(platemap_file, sheet = 'sample') %>%
    dplyr::select(-1) %>%
    pivot_longer(., cols = 1:ncol(.)) %>%
    rename(ligand = value) %>%
    dplyr::select(-c('name'))
  id_df$ligand <- gsub('empty', 'vehicle', id_df$ligand)
  conc_df <- read_excel(platemap_file, sheet = 'conc') %>%
    dplyr::select(-1) %>%
    pivot_longer(., cols = 1:ncol(.)) %>%
    rename(conc = value) %>%
    dplyr::select(-c('name'))
  df <- cbind(id_df, conc_df, df)
  message('Plate assignment attached to dataframe.')
  df$row <- as.numeric(df$row)
  df$col <- as.numeric(df$col)
  return(df)
}

add_tempheaders <- function(df,
                            start_temp = 37,
                            end_temp = 90) {
  temperature_df <-
    seq(start_temp, end_temp, by = ((end_temp - start_temp) / (ncol(df) - 1))) %>%
    round(., digits = 1)
  for (i in 1:ncol(df)) {
    colnames(df)[i] <- paste('t_', temperature_df[i], sep = '')
  }
  message('Temperature assignments changed for ',
          ncol(df),
          ' points.')
  return(df)
}

#' targetMinMax
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
targetMinMax <- function(df, plotMinMax = FALSE) {
  deriv <- df
  deriv[paste(colnames(deriv)[1], '.Deriv', sep = '')] <- 0
  deriv[paste(colnames(deriv)[2], '.Deriv', sep = '')] <- 0
  for (i in 1:nrow(deriv)) {
    deriv[i, 3] <- (deriv[i + 1, 1] + deriv[i, 1]) / 2
    deriv[i, 4] <-
      (deriv[i + 1, 2] - deriv[i, 2]) / (deriv[i + 1, 1] + deriv[i, 1])
  }
  
  # Making an extra df for the Minima/Maxima while we sort out if we can get away with removing NAs.
  df.minMax <- deriv %>%
    filter(!is.na(A1.Deriv))
  minMax.Spatial <- spatialEco::local.min.max(x = df.minMax$A1, plot = plotMixMax)
  return(minMax.Spatial)
  
  # Figure out the global minimum and maximum, and their respective temperatures.
  
  
  # Figure out which of these global values is the furthest, from the avaerage of the points.
  # How do we use a statistical test to determine this?
}

#' determineTagg
#'
#' @param df 
#' Need a spatialEco df input into this function.
#' minima, maxima, devmin, devmax
#' @param plotVisual
#' Do we plot the minima and maxima plot and show what point we choose?
#' @param forceDirection
#' Do we force the function to look in a certain direction?
#' 'neg' for global minimum
#' 'pos'for global maximum
#'
#' @return
#' Return a df with the point of maximum difference
#' @export
#'
#' @examples
determineTagg <- function(df, plotVisual = FALSE, forceDirection = '') {
  Tagg.df <- tibble()
}
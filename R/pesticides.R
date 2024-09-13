#' Pesticides data
#'
#' A synthetic proprietary pesticide usage survey data from GfK Kynetec in 2020.
#' 
#' The original data is contaminated using by adding noise and creating missing values and imputation.
#'
#' @name acresCRD
#' @docType data
#' @format A data frame with 1197 rows on the following 30 variables:\describe{
#' \item{Corn}{Corn 10, 20, 30, 40, 50, 60, 70}
#' \item{Soybean}{Soybean 10, 20, 30, 40, 50, 60, 70, 90}
#' \item{Alfalfa}{Alfalfa 10, 30, 40, 50, 70, 80}
#' \item{Pasture}{Pasture 10, 20, 30, 40, 50, 60, 70, 80, 90}}
#'
#' @keywords datasets
#' @examples
#' data(pesticides)
#'
#' calibration <- GECal::GEcalib(~ 0, dweight = d_S, data = acresCRD,
#'                               const = numeric(0),
#'                               entropy = "EL", method = "DS")
#' GECal::estimate(pesticideusage ~ 1, calibration = calibration, pimat = pimat)$estimate
#' 
#' 
#' calibration <- GECal::GEcalib(~ 0 + ., dweight = d_S, data = acresCRD,
#'                               const = total,
#'                               entropy = "SL", method = "DS")
#' GECal::estimate(pesticideusage ~ 1, calibration = calibration, pimat = pimat)$estimate
#' 
#' calibration <- GECal::GEcalib(~ 0 + ., dweight = d_S, data = acresCRD,
#'                               const = c(total),
#'                               entropy = "ET", method = "DS")
#' GECal::estimate(pesticideusage ~ 1, calibration = calibration, pimat = pimat)$estimate
#' 
#' 
#' calibration <- GECal::GEcalib(~ 0 + . + g(d_S), dweight = d_S, data = acresCRD,
#'                               const = c(total, NA),
#'                               entropy = "HD", method = "GEC")
#' GECal::estimate(pesticideusage ~ 1, calibration = calibration, pimat = pimat)$estimate
#' 
#' calibration <- GECal::GEcalib(~ 0 + ., dweight = d_S, data = acresCRD,
#'                               const = total,
#'                               entropy = "HD", method = "GEC0")
#' GECal::estimate(pesticideusage ~ 1, calibration = calibration, pimat = pimat)$estimate
NULL
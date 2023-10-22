#' ecostatscale: Statistical Scaling Functions for Ecological Systems
#' 
#' Implementation of the scaling functions presented in
#' "General statistical scaling laws for stability in
#' ecological systems" by Clark et al. 2021 in Ecology Letters
#' (DOI: 10.1111/ele.13760). Includes functions 
#' for extrapolating variability, resistance, and resilience
#' across spatial and ecological scales, as well as a basic
#' simulation function for producing time series, and a
#' regression routine for generating unbiased parameter
#' estimates. See the main text of the paper for more details.
#' 
#' @section Author:
#' Adam Clark
#'
#'
#' @docType package
#' @source Clark et al. (2021). General statistical scaling laws for stability in ecological systems. Ecology Letters. DOI:10.1111/ele.13760.
#' @name ecostatscale
#' @keywords internal 
#' "_PACKAGE"
#' @examples
#' # simulate a time series:
#' ?symdyn  # one species one patch
#' ?symdynN # multiple species or patches
#'
#' # get unbiased parameter estimates:
#' ?xt2fun
#'
#' #variance scaling function:
#' ?var_scale
#'
#' # resistance scaling function:
#' ?sd_scale
#'
#' # resilience scaling function:
#' ?res_scale
NULL

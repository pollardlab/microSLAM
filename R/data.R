#' Example Binary Gene Presence Absence Data
#'
#' This dataset contains input data used in microSLAM analysis.
#'
#' @format A binary sample-by-gene data frame: with 1 as present, 0 as absent
#' \describe{
#'   \item{sample_name}{Character. The sample identifier.}
#'   \item{gene1}{Character. Gene ID 1.}
#'   \item{gene2}{Character. Gene ID 2.}
#'   \item{...}{Additional gene IDs.}
#' }
#' @source Generated  from microSLAM
"exp_genedata"


#' Example Metadata
#'
#' This dataset contains example metadata about samples.
#'
#' @format A data frame with X rows and Y columns:
#' \describe{
#'   \item{sample_name}{Character. The sample identifier.}
#'   \item{y}{Numeric. Disease condition (e.g., "healthy", "disease").}
#'   \item{age}{Numeric. Age of the sample}
#'   \item{strain}{Numeric. Disease condition (e.g., "healthy", "disease").}
#'   \item{uncorrelated_strain}{Numeric. Disease condition (e.g., "healthy", "disease").}
#' }
#' @source Generated internally from microSLAM
"exp_metadata"

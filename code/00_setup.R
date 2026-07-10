#' ======================================================================
#' 00_setup.R
#'
#' DRGSNMI main processing setup script. Loads libraries, helper funcs, dependencies, etc.
#'
#' Inputs:  (none — setup only)
#' Outputs: (none — setup only)
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' Last modified: 2024-11-29
#' ======================================================================

### SETWD ###
# Set working directory to script location (portable: Rscript or RStudio)
local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) setwd(dirname(normalizePath(sub("^--file=", "", f))))
  else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
})


### IMPORTS ###

# Load libraries
source('utils/libraries.R')

# Load helper functions
source('utils/helper_functions.R')

# Load gene lists
source('utils/gene_lists.R')

# Load common variables
source('utils/variables.R')


### RUN ONLY ONCE ###

# Make project directories
#dir.create('../r_objects')
#dir.create('../r_objects/single_objects')
#dir.create('../r_objects/batch_objects')
#dir.create('../r_objects/study_objects')
#dir.create('../dge')
#dir.create('../dge_GO')
#dir.create('../markers')
#dir.create('../metadata')
#dir.create('../plots')
#dir.create('../figures')
#dir.create('../solo_results')

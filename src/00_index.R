# Header start ==============================================================================
# 00_index.R
#
# Author: Hennch Cornelius (cornelius.hennch@charite.de)
#
# Description: Template for analysis index script that integrates multiple 
# scripts/reports
#
# Input: 
# Output:
#
# Code written according to Hadley Wickhams "tidyverse style guide"

# Packages & dependencies ------------------------------------------------------
starttime <- Sys.time()
work_dir <- getwd()

# packages and functions
library(tidyverse)

# custom functions
source("src/00_functions.R")

# today's date
today <- format(Sys.Date(), "%Y%m%d")

# Header end ===================================================================

# print intermediate figures/tables/reports?
save_output <- FALSE
render_reports <- TRUE

# 1. Import and tidy -----------------------------------------------------------

# 2. Transform -----------------------------------------------------------------

# 3. Explore  ------------------------------------------------------------------

# 4. Visualize -----------------------------------------------------------------


# 5. Communicate  --------------------------------------------------------------
# render .Rmd report files
# they can take objects from the global environment that have been created by
# the previous analysis steps

if (render_reports) {
  # powerpoint presentations
  rmarkdown::render(input = "markdown/presentation.Rmd",
                    output_dir = "output/reports/presentations/",
                    knit_root_dir = work_dir,
                    intermediates_dir = work_dir)

  # html report
  rmarkdown::render(input = "markdown/report.Rmd",
                    output_dir = "output/reports/",
                    output_format = "html_document",
                    knit_root_dir = work_dir,
                    intermediates_dir = work_dir)

  # pdf report
  rmarkdown::render(input = "markdown/report.Rmd",
                    output_dir = "output/reports/",
                    output_format = "pdf_document",
                    knit_root_dir = work_dir,
                    intermediates_dir = work_dir)

}

# runtime
runtime <- Sys.time() - starttime
runtime




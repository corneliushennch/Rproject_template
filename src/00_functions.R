# Header start =================================================================
# 00_functions.R
#
# Author: Hennch Cornelius (cornelius.hennch@charite.de)
#
# Description: Custom functions for the current project
#
# Code written according to Hadley Wickhams "tidyverse style guide"
# Header end ===================================================================

# 1. Data wrangling functions --------------------------------------------------
## 1.1 as.numeric.factor   -----------------------------------------------------
# convert factor back to numeric
as.numeric.factor <- function(x) {
  as.numeric(levels(x))[x]
}

# 2. Model functions -----------------------------------------------------------


# 3. Plotting functions --------------------------------------------------------

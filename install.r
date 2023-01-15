# Package names
packages <- c(
  "pracma",
  "contfrac",
  "stats",
  "R.matlab",
  "writexl",
  "dplyr",
  "rbenchmark",
  "data.table",
  "tictoc",
  "NCmisc",
  "sysid",
  "minpack.lm"
)

#devtools::install_github("jessevent/crypto")
#devtools::install_github("deanfantazzini/bitcoinFinance")
#devtools::install_github("cvarrichio/rowr", force=TRUE)

# for (i in 1:length(c)) {
#  install.packages(packages[i])
# }
# suppressPackageStartupMessages(library(plyr))



# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE, quietly = TRUE))


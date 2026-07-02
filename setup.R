# Run this once after creating the mamba environment to install packages
# that are not available on conda-forge.
#
# Usage:
#   mamba env create -f environment.yml
#   mamba run -n pleiodbr Rscript setup.R

options(repos = c(
  coolbutuseless = "https://coolbutuseless.r-universe.dev",
  CRAN           = "https://cloud.r-project.org"
))

install.packages("zstdlite")

# Install pleiodbr itself (from source in the current directory)
devtools::install(".", dependencies = FALSE)

# Optional: TwoSampleMR for the MR vignette
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  install.packages('TwoSampleMR', repos = c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))
}

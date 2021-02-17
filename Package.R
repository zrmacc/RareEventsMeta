setwd("RareEventsMeta/")

# Documents the package.
# Run from within package directory.
devtools::document()

# Install package.
# Run from outside the package directory.
setwd("..")
devtools::install(pkg = "RareEventsMeta", reload = TRUE)

# Check package.
# Run from inside the package directory.
setwd("RareEventsMeta")
devtools::check()
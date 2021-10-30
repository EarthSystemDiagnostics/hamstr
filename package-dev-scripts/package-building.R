# devtools::document()
# pkgbuild::compile_dll(force = TRUE)
# devtools::load_all()

# devtools::install(quick = FALSE)
# devtools::install(quick = FALSE, dependencies = FALSE)
# devtools::install(quick = FALSE, build = FALSE, dependencies = FALSE)
# devtools::install(quick = TRUE, build = TRUE, dependencies = FALSE)

# devtools::install(quick = FALSE, build = FALSE, dependencies = FALSE)

# using ../rstanlm because already inside the rstanlm directory
install.packages("../hamstr", repos = NULL, type = "source", INSTALL_opts = "--no-multiarch", clean = TRUE)


remotes::install_github("earthsystemdiagnostics/hamstr", args = "--preclean", build_vignettes = FALSE)

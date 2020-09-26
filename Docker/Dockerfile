FROM mattocci/rstan

## Do not upgrade dependencies to avoid reinstalling rstan and rstanheaders
## Should be quite new based on base image
RUN Rscript -e 'remotes::install_github(c("earthsystemdiagnostics/hamstr"), dependencies=NA, upgrade = "never", build_vignettes=FALSE)'

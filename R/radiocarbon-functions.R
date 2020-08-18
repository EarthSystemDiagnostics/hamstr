#' @title Calibrate Radiocarbon Dates with Bchron::BchronCalibrate
#' @description Calibrates a set of 14C ages using BchronCalibrate and
#'   optionally summarises the empirical PDFs of calendar age to mean and
#'   standard deviation and appends these to the input dataframe
#' @param dat A dataframe containing the radiocarbon dates and uncertainties
#' @param age.14C Name of column with 14C ages, Default: 'age.14C'
#' @param age.14C.se Name of column with 1se 14C age uncertainty, Default:
#'   'age.14C.se'
#' @param cal_curve Calibration curve, Default: 'intcal13', see
#'   \code{\link[Bchron]{BchronCalibrate}}
#' @param return.type Return the ammended dataframe or additionally the list of
#'   PDFs, Default: 'dat'
#' @param offset Optional offset applied to all 14C ages, Default: 0
#' @return A dataframe or list
#' @details A wrapper for Bchron::Bchroncalibrate
#' @examples
#' # With defaults
#' dat <- data.frame(age.14C = c(2000, 20000),
#'                   age.14C.se = c(100, 200))
#' calibrate_14C_age(dat)
#'
#' # Change the calibration
#' calibrate_14C_age(dat, curve = "marine13")
#'
#' # Return the PDFs
#' cal.lst <- calibrate_14C_age(dat, curve = "marine13", return = "lst")
#' with(cal.lst[[2]][[1]][[1]], {plot(ageGrid, densities)})
#'
#' # Use different column names
#' dat <- data.frame(radiocarbon.age = c(2000, 20000),
#'                  se = c(100, 200))
#' calibrate_14C_age(dat, age.14C = "radiocarbon.age", age.14C.se = "se")
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso \code{\link[Bchron]{BchronCalibrate}}
#' @rdname calibrate_14C_age
#' @export
#' @importFrom Bchron BchronCalibrate
calibrate_14C_age <- function(dat, age.14C = "age.14C",
                              age.14C.se = "age.14C.se",
                              cal_curve = "intcal13",
                              return.type = "dat", offset = 0){
  
  return.type <- match.arg(return.type, choices = c("data.frame", "list"))
  cal_curve <- 
    match.arg(cal_curve,
              choices = c("intcal13", "shcal13", "marine13", "normal"))
  
  cal.ages <- lapply(1:nrow(dat), function(x) {
    tryCatch(Bchron::BchronCalibrate(
      ages = dat[[age.14C]][x] + offset,
      ageSds = dat[[age.14C.se]][x],
      calCurves = cal_curve,
      ids = x),
      error = function(i){
        cat(strsplit(as.character(i), " : ", fixed = TRUE)[[1]][2])
        
        NA
      })
  })
  
  
  # Use mean and sd of empirical PDFs as point estimates of calendar ages
  dat$age.14C.cal <- sapply(cal.ages, function(x){
    if (is.na(x) == FALSE)
    {SummariseEmpiricalPDF(x[[1]]$ageGrid, x[[1]]$densities)["mean"]} else {NA}
  })
  
  dat$age.14C.cal.se <- sapply(cal.ages, function(x){
    if (is.na(x) == FALSE)
    {SummariseEmpiricalPDF(x[[1]]$ageGrid, x[[1]]$densities)["sd"]} else {NA}
  })
  
  if (return.type == "data.frame"){
    out <- dat
  }
  
  if (return.type == "list"){
    out <- list(dat = dat, cal.ages = cal.ages)
  }
  
  return(out)
}


#' Summarise an Empirical Probability Distribution Function.
#'
#' @param x A vector of values of empirical PDF
#' @param p A vector of probabilities
#' @details Calculation of the mode is naive. For a multimodal distribution only
#' the highest is returned, in the case of 2 or more modes with exactly the same
#'  probability, the first is returned.
#' @return Returns a named vector with the mean, median, mode, and standard
#' deviation of the empirical PDF
#' @keywords internal
#' @examples
#' df <- data.frame(x = 1:10)
#' df$p <- dnorm(df$x, 5, 2)
#' SummariseEmpiricalPDF(df$x, df$p)
SummariseEmpiricalPDF <- function(x, p){
  
  # Ensure x and p are sorted
  p <- p[order(x)]
  x <- sort(x)
  
  # Ensure p sum to 1
  p <- p / sum(p)
  
  # Mean
  w.mean <- sum(x * p)
  
  # SD
  M <- sum(p > 0)
  w.sd <- sqrt(sum(p * (x-w.mean)^2) / ((M-1)/M * sum(p)))
  
  # Median
  csum.p <- cumsum(p)
  med.ind <- which.min(abs(csum.p - 0.5))
  w.median <- x[med.ind]
  
  # Mode
  max.wt <- max(p)
  n.max <- sum(p == max.wt)
  if (n.max > 1)
    warning(paste0(n.max,
                   " x with equal maximum probability. Returning the first"))
  mode <- x[which.max(p)]
  
  return(c("mean" = w.mean, "median" = w.median, "mode" = mode, "sd" = w.sd))
}


#' Compare the full empirical calendar age PDF of a radiocarbon date with a
#' t-distribution approximation
#'
#' @param age.14C vector of radiocarbon dates in years BP
#' @param age.14C.se vector of radiocarbon date uncertainties
#' @param cal_curve calibration curve
#' @param t_df degrees of freedom of the t-distribution
#' @param return.type return a ggplot object or a list containing the ggplot
#'   object and two data frames with the empirical and t-distributions
#' @return
#' @export
#'
#' @examples
#' compare_14C_PDF(age.14C = c(1000, 4000), age.14C.se = c(100, 150),
#'  cal_curve = "intcal13", return.type = "plot")

compare_14C_PDF <- function(age.14C, age.14C.se, cal_curve = "intcal13", t_df = 6,
                             return.type = c("plot", "list")){
  
  dt_ls <- function(x, dat=1, mu=0, sigma=1) 1/sigma * dt((x - mu)/sigma, dat)
  
  stopifnot(length(age.14C) == length(age.14C.se))
  
  cal.dat <- data.frame(age.14C = age.14C, age.14C.se = age.14C.se)
  
  return.type <- match.arg(return.type, choices = c("plot", "list"))
  
  cal_curve <-
    match.arg(cal_curve,
              choices = c("intcal13", "shcal13", "marine13", "normal"))
  
  calib <- calibrate_14C_age(cal.dat,
                             return.type = "list",
                             offset = 0, cal_curve = cal_curve)
  
  # The summarised calendar ages are appended to the input data
  C14 <- calib$dat
  C14$.id <- 1:nrow(cal.dat)
  
  # These are the full PDFs
  cal.ages <- calib$cal.ages
  
  cali.pdf.dat <- plyr::ldply(1:length(cal.ages), function(i){
    x <- cal.ages[[i]]
    if (is.na(x)==FALSE){
      data.frame(age = x[[1]]$ageGrid, density = x[[1]]$densities, .id = i)
      }else{
      data.frame(age = 0, density = 0, .id = i)
    }
  })
  
  t_df <- C14 %>%
    group_by(.id) %>%
    do({
      rng <- .$age.14C.cal.se * 5
      age = seq(.$age.14C.cal - rng, .$age.14C.cal + rng, length.out = 100)
      data.frame(
        age = age,
        density = dt_ls(age, dat = t_df,
                        mu = .$age.14C.cal,
                        sigma = .$age.14C.cal.se)
      )
    })
  
  gg <- cali.pdf.dat %>%
    ggplot(aes(x = age/1000, y = density, group = .id)) +
    geom_line(aes(colour = "Empirical PDF")) +
    geom_line(data = t_df, aes(y = density, colour = "t-distribution")) +
    labs(colour = "", x = "Calendar age [ka BP]", y = "Density") +
    facet_wrap(~.id, scales = "free") +
    theme_bw()
  
  if (return.type == "list"){
    return(list(plot = gg, cal.age.pdf = cali.pdf.dat, t.dist.age = t_df))
  } else if (return.type == "plot") {
    return(gg)
  }
}
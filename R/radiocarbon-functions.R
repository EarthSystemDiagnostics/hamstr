#' @title Calibrate Radiocarbon Dates with rcarbon::calibrate
#' @description Calibrates a set of 14C ages using rcarbon::calibrate and
#'   optionally summarises the empirical PDFs of calendar age to mean and
#'   standard deviation and appends these to the input dataframe
#' @param dat A dataframe containing the radiocarbon dates and uncertainties
#' @param age.14C Name of column with 14C ages, Default: 'age.14C'
#' @param age.14C.se Name of column with 1se 14C age uncertainty, Default:
#'   'age.14C.se'
#' @param cal_curve Calibration curve, Default: 'intcal20', see
#'   \code{\link[rcarbon]{calibrate}}
#' @param return.type Return the ammended dataframe or additionally the list of
#'   PDFs, Default: 'dat'
#' @param offset Name of offset column, e.g. reservoir age. If column does not 
#' exist, no offset is applied. 
#' @param offset.se Name of offset uncertainty column, e.g. sigmaDelatR. If
#'  column does not exist, no offset uncertainty is applied.
#' @importFrom rcarbon calibrate
#' @return A dataframe or list
#' @details A wrapper for rcarbon::calibrate
#' @examples
#' # With defaults
#' dat <- data.frame(age.14C = c(2000, 20000),
#'                   age.14C.se = c(100, 200))
#' calibrate_14C_age(dat)
#'
#' # Change the calibration
#' calibrate_14C_age(dat, cal_curve = "marine13")
#'
#' # Return the PDFs
#' cal.lst <- calibrate_14C_age(dat, cal_curve = "marine13", return = "list")
#' cal.lst

#' # Use different column names
#' dat <- data.frame(radiocarbon.age = c(2000, 20000),
#'                  se = c(100, 200))
#' calibrate_14C_age(dat, age.14C = "radiocarbon.age", age.14C.se = "se")
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso \code{\link[rcarbon]{calibrate}}
#' @rdname calibrate_14C_age
#' @export
calibrate_14C_age <- function(dat,
                                age.14C = "age.14C",
                                age.14C.se = "age.14C.se",
                                cal_curve = "intcal20",
                                return.type = "dat",
                                offset = "offset", 
                                offset.se = "offset.se") {
  return.type <-
    match.arg(return.type, choices = c("data.frame", "list"))
  cal_curve <-
    match.arg(
      cal_curve,
      choices = c(
        "intcal20",
        "marine20",
        "shcal20",
        "intcal13",
        "marine13",
        "shcal13",
        "normal"
      )
    )
  
  if (is.null(dat[[offset]])) {
    dat$offset <- rep(0, nrow(dat))
  } else{
    dat$offset <- dat[[offset]]
  }
  
  if (is.null(dat[[offset.se]])) {
    dat$offset.se <- rep(0, nrow(dat))
  } else{
    dat$offset.se <- dat[[offset]]
  }
  
  
  cal2tib <- function(cal){
    a <- dplyr::as_tibble(cal[[1]])
    b <- dplyr::bind_rows(cal[[2]], .id = "DateID") 
    dplyr::left_join(a, b, by = "DateID")
  }
  
  ## need to catch errors on a per point basis
  cal.ages <- lapply(1:nrow(dat), function(x) {
    tryCatch(rcarbon::calibrate(dat[[age.14C]][x], dat[[age.14C.se]][x],
                                resOffsets = -dat[["offset"]][x],
                                resErrors = dat[["offset.se"]][x],
                                calCurves = cal_curve, ids = x, verbose = FALSE),
             error = function(i){
               cat(strsplit(as.character(i), " : ", fixed = TRUE)[[1]][2])
               NA
             })
  })
  
  cal.ages <- lapply(cal.ages, cal2tib) %>%
    dplyr::bind_rows(.) %>%
    dplyr::mutate(DateID = as.numeric(DateID))

  cal.ages.sum <- cal.ages %>%
    dplyr::group_by(DateID) %>%
    dplyr::summarise(dplyr::as_tibble(as.list(
      suppressWarnings(hamstr:::SummariseEmpiricalPDF(calBP, PrDens))
    ))) %>%
    dplyr::select(-mode,-mean)
  
  dat <- cal.ages.sum %>%
    dplyr::mutate(DateID = as.numeric(DateID)) %>%
    dplyr::arrange(DateID) %>%
    dplyr::bind_cols(dat, .) %>%
    dplyr::mutate(age.14C.cal = median, age.14C.cal.se = sd) %>%
    dplyr::select(-median,-sd)
  
  if (return.type == "data.frame") {
    out <- dat
  }
  
  if (return.type == "list") {
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
#' \dontrun{
#' df <- data.frame(x = 1:10)
#' df$p <- dnorm(df$x, 5, 2)
#' hamstr:::SummariseEmpiricalPDF(df$x, df$p)
#'
#' x <- 1:100
#' y <- dnorm(x, 50, 10)
#' plot(x, y)
#' SummariseEmpiricalPDF(x, y)
#'
#'
#' x2 <- x[c(1:50, seq(51, 70, 3), seq(71, 100, 1))]
#' y2 <- y[x2]
#'
#' plot(x2, y2)
#' SummariseEmpiricalPDF(x2, y2)
#' }
SummariseEmpiricalPDF <- function(x, p){

  # Ensure x and p are sorted
  p <- p[order(x)]
  x <- sort(x)


  # Allow for changes in resolution of x
  dx <- diff(x)

  if (max(abs(dx - mean(dx))) > stats::median(dx) / 100){
    warning("Empirical PDF has variable resolution - this is accounted for but the results may be less reliable.")

  }

  dx <- c(dx[1], dx)

  p <- p * dx

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
  if (n.max > 1) warning(
    paste0("In mode calculation, these values have equal probability: ",
           paste0(x[p == max.wt], collapse = ", "), ". returning first"))
  mode <- x[which.max(p)]

  return(c("mean" = w.mean, "median" = w.median, "mode" = mode, "sd" = w.sd))
}



#' Compare Empirical and t-distribution Approximated Calendar Age PDFs
#' @description Compare the full empirical calendar age PDFs of calibrated
#'   radiocarbon dates with the t-distribution approximations use by hamstr
#' @param age.14C vector of radiocarbon dates in years BP
#' @param age.14C.se vector of radiocarbon date uncertainties
#' @param offset vector of offsets, e.g. reservoir ages. 
#' @param offset.se vector of offset uncertainties, e.g. sigmaDelatR.
#' @param cal_curve calibration curve
#' @param nu degrees of freedom of the t-distribution approximation, default in
#'   hamstr is 6
#' @param return.type return a ggplot object or a list containing the ggplot
#'   object and two data frames with the empirical and t-distributions
#' @return A ggplot2 object or list with data and ggplot2 object
#' @export
#' @importFrom rlang .data
#' @examples
#' compare_14C_PDF(age.14C = c(1000, 4000), age.14C.se = c(100, 150),
#'  cal_curve = "intcal20", return.type = "plot")
compare_14C_PDF <- function(age.14C, age.14C.se,
                            offset = 0, offset.se = 0,
                            cal_curve = "intcal20", nu = 6,
                            return.type = c("plot", "list")
                            ){

  dt_ls <- function(x, dat=1, mu=0, sigma=1) {
    1/sigma * stats::dt((x - mu)/sigma, dat)
    }

  stopifnot(length(age.14C) == length(age.14C.se))

  cal.dat <- data.frame(age.14C = round(age.14C),
                        age.14C.se = round(age.14C.se),
                        offset = offset,
                        offset.se = offset.se)

  return.type <- match.arg(return.type, choices = c("plot", "list"))

  cal_curve <-
    match.arg(cal_curve,
              choices = c("intcal20", "marine20", "shcal20",
                          "intcal13", "marine13", "shcal13",
                          "normal"))

  calib <- calibrate_14C_age(cal.dat,
                             return.type = "list",
                             cal_curve = cal_curve)

  # The summarised calendar ages are appended to the input data
  C14 <- calib$dat
  C14$id <- 1:nrow(cal.dat)

  # These are the full PDFs
  cal.ages <- calib$cal.ages

  cal.ages.df <- cal.ages %>% 
    dplyr::mutate(id = DateID)

  cali.pdf.dat <- cal.ages.df %>%
    dplyr::group_by(.data$id) %>%
    dplyr::reframe(
      age = calBP,
      density = PrDens
        )


  tdist_df <- C14 %>%
    dplyr::filter(is.na(.data$age.14C.cal) == FALSE) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::do({
      rng <- .data$age.14C.cal.se * 5
      age <- seq(.data$age.14C.cal - rng, .data$age.14C.cal + rng, length.out = 100)
      data.frame(
        age = age,
        density = dt_ls(age, dat = nu,
                        mu = .data$age.14C.cal,
                        sigma = .data$age.14C.cal.se)
      )
    })

  gg <- cali.pdf.dat %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$age/1000, y = .data$density, group = .data$id)) +
    ggplot2::geom_line(ggplot2::aes(colour = cal_curve)) +
    ggplot2::geom_line(data = tdist_df, ggplot2::aes(y = .data$density, colour = paste0("t-distribution: nu = ", nu))) +
    ggplot2::labs(colour = "", x = "Calendar age [ka BP]", y = "Density") +
    ggplot2::facet_wrap(~.data$id, scales = "free") +
    ggplot2::theme_bw()

  if (return.type == "list"){
    return(list(plot = gg, cal.age.pdf = cali.pdf.dat, t.dist.age = tdist_df))
  } else if (return.type == "plot") {
    return(gg)
  }
}

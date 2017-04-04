## Import Example Data

## Example core from Bacon v2.2
MSB2K <- read.csv("data-raw/MSB2K.csv") %>% 
  dplyr::tbl_df()
devtools::use_data(MSB2K, overwrite = TRUE)

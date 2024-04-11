source("util.R")

## NOTE: We are keeping this task as it is part of the original structure of 
##       our workflow for running the analysis as published. Please do note 
##       that this task handles identifiable data, which under GDPR
##       regulations CANNOT be made publicly available. Hence, we are solely
##       keeping the aggregated data that is necessary to reproduce all our
##       analysis. Please also note that a number of columns are left blank
##       (NA_real_), such as age-specific PCR cases and deaths, as these are
##       expected by the model given its inference machinery structure, but 
##       we did not use any of such data in our analysis.
##       
dat <- read_csv("data_timeseries.csv")
write_csv(dat, "data_timeseries.csv")

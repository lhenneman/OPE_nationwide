library( ncdf4)
library( data.table)
library( searchAQ)

colrow <- expand.grid( col.n = 1:152,
                       row.n = 1:116)
yrs <- 1990:2010

starttime <- Sys.time()

## == ## == ## == ## == ## == ## == ## == ## == ## 
##       Define import function functions       ##
## == ## == ## == ## == ## == ## == ## == ## == ## 

read.ioapi <- function( yr,
                        string_start,
                        col.n,
                        row.n){
  
  ## define file name and open ncdf file
  file_in <- paste0( string_start, yr, '.ioapi')
  ncin_in <- nc_open( file_in)
  
  ## read variable names and import time
  vnames <- names(ncin_in$var)[-1]
  time <- ncvar_get( ncin_in,'TFLAG', count = c(1,1,-1))
  time.date <-  as.IDate( as.character( time), format = '%Y%j')
  print( paste0( 'Length of time variable for year ', yr, ' is ', length( time)))
  
  ## define start and count strings that define which data to import
  start <- c( col.n, row.n, 1, 1)
  count <- c( 1,     1,     1, length( time))
  
  ## read variables
  daily_vals <- sapply( vnames,
                        function( variable)
                          ncvar_get( ncin_in, variable, start = start, count = count))
  
  ## combine with dates and export
  daily_vals.dt <- data.table( date = time.date, 
                               daily_vals)
  return( daily_vals.dt)
}

## == ## == ## == ## == ## == ## == ## == ## == ## 
##              Define worker function          ##
## == ## == ## == ## == ## == ## == ## == ## == ## 
calc.ope <- function( n,
                      colrow,
                      yrs = yrs) {
  
  colrow <- data.table( colrow)
  coln <- unlist( colrow[ n, .(col.n)])
  rown <- unlist( colrow[ n, .(row.n)])
  
  ## ========== import 11am-3pm file ============= ##
  string_11am3pm <- "~/Dropbox/Harvard/Meetings_and_People/ChristianEPA/Analysis/CMAQdata/hr2day_11am3pm_doe36km_"
  
  met_11am3pm <- rbindlist( lapply( yrs,
                                    read.ioapi,
                                    string_start = string_11am3pm,
                                    col.n = coln,
                                    row.n = rown))
  
  ## ========== import 2pm-3pm file ============= ##
  string_2pm3pm <- "~/Dropbox/Harvard/Meetings_and_People/ChristianEPA/Analysis/CMAQdata/hr2day_2pm3pm_doe36km_"
  
  conc_2pm3pm <- rbindlist( lapply( yrs,
                                    read.ioapi,
                                    string_start = string_2pm3pm,
                                    col.n = coln,
                                    row.n = rown))
  
  ## ========== import 24hours file ============= ##
  string_24hours <- "~/Dropbox/Harvard/Meetings_and_People/ChristianEPA/Analysis/CMAQdata/hr2day_24hours_doe36km_"
  
  conc_24hours <- rbindlist( lapply( yrs,
                                     read.ioapi,
                                     string_start = string_24hours,
                                     col.n = coln,
                                     row.n = rown))
  
  ## ========== complete the dataset ============= ##
  date.complete <- data.table( date = seq.Date( as.Date( paste0( yrs[1], '-01-01')),
                                                as.Date( paste0( yrs[length(yrs)], '-12-31')),
                                                by = 'day'))
  
  data.dt <- merge( merge( merge( date.complete,
                                  met_11am3pm, by = 'date', all = T), 
                           conc_2pm3pm, by = 'date', all = T), 
                    conc_24hours, by = 'date', all = T)
  data.dt[data.dt < -1e30] <- NA
  
  data.dt[, `:=` (NOX_AVG = NO2_AVG + NO_AVG,
                  NOZ_AVG = NOY_AVG - NO2_AVG - NO_AVG)]
  
  ## ========== detrend the data ============= ##
  
  gas <- data.dt[, .(date, NOZ_AVG, NOY_AVG, MDA8O3)]
  met  <- data.dt[, .(date, RH_AVG, TE_AVG,  WS_AVG,  TEDMAX)]
  rf   <- data.dt[, .(date, RAINTOT)]
  
  setnames( met, 
            c( 'RH_AVG',    'TE_AVG',      'WS_AVG',    'TEDMAX'),
            c( 'rh_midday', 'temp_midday', 'ws_midday', 'temp_max'))
  
  detrend <- sapply( names( gas)[-1],
                     detrending_SEARCH_gm,
                     data.frame( gas),
                     data.frame( met),
                     data.frame( rf),
                     simplify = F,
                     USE.NAMES = T)
  
  detrend.dt <- data.table( date = detrend$MDA8O3$data$date,
                            o3obs = detrend$MDA8O3$data$obs,
                            ps = detrend$MDA8O3$data$PS,
                            noz = detrend$NOZ_AVG$data$obs)
  
  
  ## ========== calculate OPE ============= ##
  # subset data
  detrend.dt.highps <- detrend.dt[ps > quantile( ps,
                                                 probs = .8,
                                                 na.rm = T)]
  
  #calculate OPE, extract correct data
  ope <- make_splineOPE( y = detrend.dt.highps$o3obs,
                         x = detrend.dt.highps$noz,
                         nk = 3)
  detrend.dt.highps[, ope := ope$ope]
  
  out <- list( detrenddata = detrend.dt.highps,
               intercept = ope$intercept)
  
}


## == ## == ## == ## == ## == ## == ## == ## == ## 
##              Run worker function             ##
## == ## == ## == ## == ## == ## == ## == ## == ## 
ope1000 <- calc.ope( 10000,
                     colrow,
                     yrs = yrs)


## ========== Check time ============= ##
print( paste( 'Runtime is', Sys.time() - starttime))


plot( ope1000$detrenddata$date,
      ope1000$detrenddata$ope)
plot( ope1000$detrenddata$date,
      ope1000$detrenddata$ps)
plot( ope1000$detrenddata$ps,
      ope1000$detrenddata$ope)
plot( ope1000$detrenddata$ps,
      ope1000$detrenddata$o3obs)
plot( ope1000$detrenddata$noz,
      ope1000$detrenddata$ope)

library(rvest)
#Function adjusted from Trinh Nguyen's code to pull salvage datasets from SacPAS
pull_salvage <- function(salvageURL = "http://www.cbr.washington.edu/sacramento/data/query_loss_detail.html") {
  startingSession <- session(salvageURL)
  startingForm <- html_form(startingSession)[[1]]

  df <- lapply(startingForm$fields$year$options, function(x) {
    filledForm <- set_values(startingForm,
                             year = x,
                             species = "26:all") # This part tells it what species

    submittedFormURL <- suppressMessages(submit_form(session = startingSession,
                                                     form = filledForm, POST = salvageURL)$url)

    csvLink <- submittedFormURL

    if (length(csvLink) == 0) {
      return(NULL)
    } else {
      csvDownload <- csvLink
    }

    df <- csvDownload %>%
      read_csv()  }) %>%
    bind_rows()
  df
}

# Get NWIS Flow -------------------------------------------------------------------
# Default is Vernalis
f_get_NWIS_flow <- function(siteNumbers=11303500, parameterCd = c('00060'), startDate = start, endDate = end, tz = "Etc/GMT+8"){

  # get data
  print("Downloading data...")
  data <- dataRetrieval::readNWISuv(siteNumbers, parameterCd, startDate, endDate, tz)

  # fix names

  data2 <- dataRetrieval::renameNWISColumns(data)

  print("Data downloaded!")

  # clean names
  data2 <- janitor::clean_names(data2)

  # write out
  saveRDS(data2, paste0("data_raw/USGS_NWIS_", siteNumbers, "_flow.rds"))

  # print message!
  print("Data saved in data_raw")
}

# Clean CDEC -------------------------------------------------------------------

clean_cdec <- function(cdec_df)

cdec_df %>%
  bind_rows() %>%
  mutate(datetime = ymd_hms(datetime),
         date = date(datetime),
         month = as.numeric(month)) %>%
  select(datetime,date,month, wy = water_year, station = station_id, parameter=value) %>%
  filter(!is.na(station),
         !is.na(date))

library(deltafish)
library(dplyr)
library(readr)


# Chipps --------------------------------------------------------
devtools::install_github("Delta-Stewardship-Council/deltafish")
create_fish_db()

# open our two data files
surv <- open_survey()
fish <- open_fish()

# filter for sources and taxa of interest
surv_DJFMP <- surv %>%
  filter(Source == "DJFMP")


fish_salmonids <- fish %>%
  filter(Taxa %in% c("Oncorhynchus mykiss", "Oncorhynchus tshawytscha"))


# do a join and collect the resulting data frame
# collect executes the sql query and gives you a table
df <- left_join(surv_DJFMP, fish_salmonids) %>%
  collect() %>%
  filter(Method != "Beach seine")



# From EDI -----------------------------------------------------




# Salvage ------------------------------------------------------
library(rvest)
#Function adjusted from Trinh Nguyen's code to pull salvage datasets from SacPAS
pull_salvage <- function(salvageURL = "http://www.cbr.washington.edu/sacramento/data/query_loss_detail.html") {
  startingSession <- session(salvageURL)
  startingForm <- html_form(startingSession)[[1]]

  df <- lapply(startingForm$fields$year$options, function(x) {
    filledForm <- set_values(startingForm,
                             year = x,
                             species = "1:all") # This part tells it what species

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

salmon <- pull_salvage()

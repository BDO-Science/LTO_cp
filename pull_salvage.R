smelt <- pull_salvage() %>%
  filter(!is.na(Species))

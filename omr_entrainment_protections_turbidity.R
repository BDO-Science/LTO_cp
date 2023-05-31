library(tidyverse)
library(rvest)
library(lubridate)
library(rgdal)
library(janitor)
library(dataRetrieval)
library(sharpshootR)
library(padr)
library(readr)
library(deltafish)


###########
# LTO work to look at delta smelt entrainment protection dates (hypothetical) - turbidity bridge and larval/juvenile protections
# Catarina Pien, come code from Brian Mahardja (adult turbidity downloads)
# Last updated 5/30/2023
############


# Adult --------------------------------------------------------------
ff_dates <- read_csv(here::here("data_clean", "FirstFlushFirstDateExceeded.csv")) %>%
  rename(wy = WY) %>%
  rename(date.ff = date) %>%
  mutate(date.end.ff = date.ff + 14) %>%
  select(wy, date.ff, date.end.ff)


## Download turbidity --------------------------------------------
start_date<-"2011-12-01"
end_date<-"2023-03-15"

sd_turbidity <- readRDS("data_raw/turbidity_cdec_southdelta_2000-2022.rds")

## Get turbidity data
OH4_turb <- CDECquery(id='OH4', sensor=221, interval='E', start=start_date, end=end_date)
OBI_turb <- CDECquery(id='OBI', sensor=221, interval='E', start=start_date, end=end_date)
HOL_turb <- CDECquery(id='HOL', sensor=221, interval='E', start=start_date, end=end_date)
# ORQ_turb <- CDECquery(id='ORQ', sensor=221, interval='E', start=start_date, end=end_date)

turb_data<-bind_rows(OH4_turb,OBI_turb,HOL_turb) %>%
  filter(value <200 & value > 0)

# saveRDS(turb_data, "data_clean/turbidity_data_OH4_OBI_HOL.csv")

## Calculate daily average turbidity -----------------------
turb_mean<- turb_data %>% mutate(date=date(datetime)) %>%
  group_by(date,station_id) %>%
  summarise(turbidity=mean(value,na.rm=T)) %>%
  ungroup() %>%
  mutate(year = lubridate:: year(date),
         month = lubridate::month(date),
         wy = ifelse(month>9, year+1, year)) %>%
  pad %>%
  arrange(station_id, date) %>%
  filter(!is.na(date)) %>%
  filter(!is.na(station_id))

## Turbidity bridge met 2020 ROD -------------------
turb_old <- turb_mean %>%
  filter(station_id == "OBI") %>%
  mutate(thresh = if_else(turbidity>=12, 1L, 0L))

turb_instances_old <- turb_old %>%
  group_by(x1 = cumsum(replace(thresh, is.na(thresh), 0) == 0)) %>%
  mutate(counter = (row_number() -1 * thresh)) %>%
  ungroup%>%
  mutate(counter = replace(counter, thresh == 0, 0))%>%
  mutate(Reg = "2020 ROD")

turb_old_ff <- turb_instances_old %>%
  arrange(date) %>%
  left_join(ff_dates, by = "wy") %>%
  mutate(Feb1 = paste0(year, "-02-01"),
         April1 = paste0(year, "-04-01"),
         afterFF = if_else(thresh == 1L & is.na(date.end.ff) & date>=Feb1, 1L,
                           if_else(thresh == 1L & date>date.end.ff , 1L,
                                   if_else(thresh == 1L & date>Feb1, 1L,
                                           0L))))


turb_old_sum <- turb_old_ff %>%
  filter(afterFF == 1L,
         thresh == 1L,
         date<April1) %>%
  group_by(wy) %>%
  summarize(instance_ROD = length(unique(x1)))



## Turbidity bridge met PA ---------------------
turb_new <- turb_mean %>%
  pivot_wider(names_from = "station_id", values_from = "turbidity") %>%
  mutate(thresh = if_else(HOL >=12 & OBI >=12 & OH4 >=12, 1L, 0L)) %>%
  ungroup()

turb_instances_new <- turb_new %>%
  group_by(x1 = cumsum(replace(thresh, is.na(thresh), 0) == 0)) %>%
  mutate(counter = (row_number() -1 * thresh)) %>%
  ungroup%>%
  mutate(counter = replace(counter, thresh == 0, 0))

turb_new_ff <- turb_instances_new %>%
  arrange(date) %>%
  left_join(ff_dates, by = "wy") %>%
  mutate(Feb1 = paste0(year, "-02-01"),
         April1 = paste0(year, "-04-01"),
    afterFF = if_else(thresh == 1L & is.na(date.end.ff) & date>=Feb1, 1L,
                      if_else(thresh == 1L & date>date.end.ff , 1L,
                              if_else(thresh == 1L & date>Feb1, 1L,
                                      0L))))

turb_new_sum <- turb_new_ff %>%
  filter(thresh == 1L,
         afterFF == 1L,
         date<April1) %>%
  group_by(wy) %>%
  summarize(instance_PA = length(unique(x1)))

turb_instances_long <- turb_new_ff %>%
  pivot_longer(cols = c("HOL", "OBI", "OH4"), names_to = "station_id", values_to = "turbidity")%>%
  mutate(Reg = "Proposed Action")

# turbidity summary table
turb_sum <- left_join(turb_old_sum, turb_new_sum)

# turbidity table for plotting
turb_all <- rbind(turb_instances_long, turb_old_ff) %>%
  rename(action_triggered = afterFF) %>%
  mutate(action_triggered = factor(action_triggered)) %>%
  filter(!is.na(action_triggered))

## Comparison plots ----------------------
(plot_turb_all <- ggplot(turb_all) +
  geom_point(aes(date, turbidity, shape = action_triggered, size = action_triggered, color = station_id)) +
  geom_hline(yintercept = 12) +
  facet_wrap(~Reg, nrow = 2, scales = "free_y") +
  scale_shape_manual(values = c(1,17)) +
    scale_size_manual(values = c(1, 3)) +
  scale_color_manual(values = c("pink", "steelblue4", "orange")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%b", expand = c(0,0))+
  labs(y = "Mean Daily Turbidity (FNU)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.title.x = element_blank(),
        legend.position = "top",
        strip.text = element_text(size=11)))

## Write -------------------------
png(filename = here::here("figures", "tba_old_new.png"), width = 8, height = 6, units = "in", family = "sans", res = 300)
plot_turb_all
dev.off()

# write_csv(turb_sum, "data_clean/turbidity_bridge_avoidance_instances.csv")

# Juvenile  ---------------------------------------------------------

## Pull Dayflow Data ------------------------------

dayflow_1984_1996<- read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/cb04e626-9729-4105-af81-f6e5a37f116a/download/dayflow-results-1984-1996.csv")%>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))
dayflow_1997_2020<- read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/21c377fe-53b8-4bd6-9e1f-2025221be095/download/dayflow-results-1997-2020.csv")%>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))
dayflow_2021<-read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/83122ce7-e7f5-4ad1-b5e9-6a7032cb1117/download/dayflowcalculations2021.csv")

dayflow_combined <- bind_rows(dayflow_1984_1996,dayflow_1997_2020,dayflow_2021) %>%
  # using the newer col names
  mutate(EXPORTS = coalesce(EXPORT, EXPORTS),
         DIVER = coalesce(DIVE, DIVER),
         EFFEC = coalesce(EFFECT, EFFEC),
         EFFDIV = coalesce(EFFD, EFFDIV), Month=month(Date)) %>% select(-c(Y, EXPORT, DIVE, EFFECT, EFFD))
# remove(dayflow_1984_1996,dayflow_1997_2020,dayflow_2021)
#filter just relevant columns
qwest <- dayflow_combined %>% select(Date, WEST) %>% rename(QWEST = WEST) %>% filter(QWEST > -20000)


## Pull deltafish secchi depth data -------------------------
create_fish_db()
surv <- open_survey()
fish <- open_fish()

# filter for sources and taxa of interest
surveys <- surv %>%
  filter(Source %in% c("EDSM", "20mm", "SLS", "DJFMP"))

# do a join and collect the resulting data frame
# collect executes the sql query and gives you a table
df <- left_join(surveys, fish) %>%
  collect() %>%
  select(-Taxa, -Notes_catch, -Count, -Length) %>%
  distinct()


# filter stations
library(sf)
library(deltamapr)

df_sf <- df %>%
  ungroup() %>%
  filter(!is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
df_sf <- st_transform(df_sf, crs = st_crs(R_EDSM_Regions_1718P1))

df_south <- df_sf %>%
  sf::st_join(R_EDSM_Regions_1718P1, join = st_within) %>%
  filter(Region == "South")

WW_Delta <- st_transform(WW_Delta, crs = st_crs(R_EDSM_Regions_1718P1)) %>%
  filter(HNAME!= "SAN FRANCISCO BAY")

bbox(R_EDSM_Regions_1718P1)
map_data <- ggplot() +
  geom_sf(data = R_EDSM_Regions_1718P1, aes(fill = Region), alpha = 0.5, inherit.aes = FALSE)+
  geom_sf(data = WW_Delta, color = "white", alpha = 0.5, inherit.aes = FALSE)+
geom_sf(data = df_south, color = "black", shape = 1, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(560000, 650000)) +
  scale_y_continuous(limits = c(4180000, 4270000)) +
  theme_bw()

png(filename = here::here("figures", "larval_smelt_secchi_map.png"), width = 7, height = 7, units = "in", family = "sans", res = 300)
map_data
dev.off()

# saveRDS(df, "data_raw/secchi_depths.rds")
# saveRDS(df_south, "data_raw/secchi_depths_south.rds")

## Join secchi depth and qwest data -----------------------------------------
larval0 <- left_join(df_south, qwest, by = "Date") %>%
  select(Date, Station, Secchi, QWEST, Region) %>%
  unique() %>%
  filter(!is.na(Secchi))

larval <- larval0 %>%
  mutate(Date2 = ymd(paste0("1980-", month(Date), "-", day(Date))),
         WY = ifelse(month(Date)>=10, year(Date) + 1, year(Date) )) %>%
  mutate(week = week(Date)) %>%
  group_by(WY, week) %>%
  mutate(meanSecchi = mean(Secchi, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(thresh_Secchi = if_else(meanSecchi < 100, 1L, 0L),
         thresh = if_else(QWEST < 0 & meanSecchi < 100, 1L, 0L))

secchi_thresh <- larval %>%
  pivot_longer(cols = c("Secchi", "QWEST"), values_to = "value", names_to = "parameter") %>%
  filter(!is.na(thresh)) %>%

  filter(WY > 2009 & WY<2020,
         month(Date) >=3 & month(Date)<7) %>%
  rename(larjuv_action_triggered = thresh) %>%
  mutate(larjuv_action_triggered = factor(larjuv_action_triggered))


## Plot and save data ---------------------------------
secchi_thresh_only <- secchi_thresh %>% filter(parameter == "Secchi")

(plot_secchi_south <- ggplot(secchi_thresh_only) +
  geom_point(aes(Date2, value, color = larjuv_action_triggered, size = larjuv_action_triggered, shape = larjuv_action_triggered)) +
  facet_wrap(~WY, scales= "free_y", nrow = 12, strip.position = "right") +
    geom_hline(yintercept = 100, linetype = "dotted") +
  scale_shape_manual(values = c(1,17)) +
  scale_size_manual(values = c(1, 3)) +
  scale_color_manual(values = c( "gray80", "steelblue4")) +
  scale_x_date(date_breaks = "1 months", date_labels = "%b", expand = c(0,0))+
  labs(y = "Secchi Depth (cm)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.title.x = element_blank(),
        legend.position = "top",
        strip.text = element_text(size=11)))

png(filename = here::here("figures", "larval_smelt_protections.png"), width = 8, height = 7, units = "in", family = "sans", res = 300)
plot_secchi_south
dev.off()

secchi_thresh_dates <- secchi_thresh %>%
  filter(thresh == 1)

 # write_csv(secchi_thresh_dates, "data_clean/secchi_thresh_dates.csv")


## New -----------------
secchi_thresh_new <- larval %>%
  mutate(Date2 = ymd(paste0("1980-", month(Date), "-", day(Date))),
         WY = ifelse(month(Date)>=10, year(Date) + 1, year(Date) )) %>%
  filter(WY > 2009 & WY<2020,
         month(Date) >=3 & month(Date)<7) %>%
  rename(larjuv_action_triggered = thresh_Secchi) %>%
  mutate(larjuv_action_triggered = factor(larjuv_action_triggered)) %>%
  arrange(Date)

(plot_secchi_new <- ggplot(secchi_thresh_new) +
    geom_point(aes(Date2, Secchi, color = larjuv_action_triggered, size = larjuv_action_triggered, shape = larjuv_action_triggered)) +
    geom_hline(yintercept = 100, linetype = "dotted") +
    facet_wrap(~WY, scales= "free_y", nrow = 12, strip.position = "right") +
    scale_shape_manual(values = c(1,17)) +
    scale_size_manual(values = c(1, 3)) +
    scale_color_manual(values = c( "gray80", "steelblue4")) +
    scale_x_date(date_breaks = "1 months", date_labels = "%b", expand = c(0,0))+
    labs(y = "Secchi Depth (cm)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 10),
          axis.title.x = element_blank(),
          legend.position = "top",
          strip.text = element_text(size=11)))

png(filename = here::here("figures", "larval_smelt_protections_new.png"), width = 8, height = 7, units = "in", family = "sans", res = 300)
plot_secchi_new
dev.off()

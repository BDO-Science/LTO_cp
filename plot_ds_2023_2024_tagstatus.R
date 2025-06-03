
adults_2023 <- read_csv("data_clean/adult_dsm_2023.csv") %>%
  mutate(Tagged = if_else(Mark == "None", "Not tagged", "Tagged")) %>%
  group_by(Tagged) %>%
  summarize(n = n())%>%
  mutate(Year = "2023")

adults_2024 <- read_csv("data_clean/dsm_adult_2024.csv") %>%
  mutate(Tagged = if_else(Mark == "None", "Not tagged", "Tagged")) %>%
  group_by(Tagged) %>%
  summarize(n = n()) %>%
  mutate(Year = "2024")

adults <- rbind(adults_2023, adults_2024)

ggplot(adults, aes(x = Tagged, y = n)) + geom_col(aes(fill = Tagged)) +facet_wrap(~Year, nrow = 2) +
  labs(y = "Count") +
  scale_fill_manual(values= c("goldenrod", "steelblue4"))  +
                      theme_bw()

yearplot <- ggplot(adults, aes(x = Year, y = n, fill = Tagged)) + geom_col(aes(fill = Tagged)) +
  labs(y = "Count", x = "Water Year") +
  scale_fill_manual(values= c("goldenrod", "steelblue4"))  +
  theme_bw()

png(filename = here("figures", "dsm_tagproportion_20232024.png"), width = 4, height = 4, units = "in", pointsize = 12, family = "sans", res = 300)
yearplot
dev.off()

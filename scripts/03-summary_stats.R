library(tidyverse)
library(here)
library(metafor)
library(meta)
library(glue)
library(openxlsx)
library(paletteer)
library(ggrepel)
library(extrafont)
library(rcartocolor)
library(patchwork)
library(bernr)

data <- read_rds(here("processed", "01-hbv_fibrosis.rds"))

data$n_chb %>% sum()
data$country %>% n_distinct()

n = data %>% 
  group_by(country, studypop2) %>% 
  summarise(n_cohorts = n(), 
            n_patients = sum(n_chb))


data %>% 
  group_by(studypop2) %>% 
  summarise(n_cohorts = n(),
            n_patients = sum(n_chb)) %>% 
  ungroup() %>% 
  mutate(p = round(n_patients / sum(n_patients)*100, 1))

data %>% 
  count(country) %>% 
  mutate(sum = sum(n))

data %>% 
  select(group, country, n_chb) %>% 
  group_by(country) %>% 
  summarise(n_cohort = n(), 
            n_patients = sum(n_chb)) %>% 
  ungroup() %>% 
  mutate(prop_patients = paste0(round(n_patients / sum(n_patients) * 100, 1), "%")) %>% 
  arrange(desc(n_patients))

data %>% 
  select(group, country, n_chb) %>% 
  group_by(group) %>% 
  summarise(n_cohort = n(), 
            n_patients = sum(n_chb)) %>% 
  ungroup() %>% 
  mutate(prop_patients = paste0(round(n_patients / sum(n_patients) * 100, 1), "%")) %>% 
  arrange(desc(n_patients))

data$median_age_hbvpatients %>% range(na.rm = TRUE)
data$prop_hbe_ag_tested_positive %>% range(na.rm = TRUE)

data %>% 
  select(n_chb, n_chb_female) %>% 
  mutate(p_female = n_chb_female / n_chb * 100) %>% 
  pull(p_female) %>% 
  range(na.rm = TRUE)

n %>% ungroup() %>% 
  mutate(n = row_number()) %>% 
  mutate(cat = paste(country, "//", studypop2)) %>% 
  ggplot(aes(x = n, y = 1)) + 
  geom_point(aes(size = n_patients, color = cat)) + 
  scale_size(range = c(2, 10))



n %>% 
  ungroup() %>% 
  group_by(studypop2) %>% 
  summarise(n_cohorts = sum(n_cohorts),
            n_patients = sum(n_patients)) %>%
  ungroup() %>% 
  mutate(sum_patients = sum(n_patients),
         prop_patients = round(n_patients/sum_patients*100,1),
         sum_cohorts = sum(n_cohorts)) %>% 
  select(studypop2,n_patients, prop_patients, n_cohorts, sum_cohorts, sum_patients)

data %>% 
  count(group)

data %>% 
  group_by(group) %>% 
  summarise(n_cohorts  = n(),
            n_patients = sum(n_chb)) %>% 
  ungroup() %>% 
  mutate(prop = round(100*n_patients/sum(n_patients), 1))


fibtest <- data %>% 
  gather(n_sig_fibrosis_fibroscan,n_sig_fibrosis_apri, 
         n_sig_fibrosis_fib4, n_sig_fibrosis_fibrotest,
         key = fib_test, value = fibrosis) %>% 
  mutate(fib_test = str_remove(fib_test, "n_sig_fibrosis_")) %>% 
  mutate(fib_test = factor(fib_test, levels = c(
    "fibroscan", "fib4", "apri", "fibrotest"
  ))) %>% 
  filter(!is.na(fibrosis)) %>% 
  group_by(group, author, year, country) %>% 
  arrange(group, author, year, country, fib_test) %>% 
  slice(1) %>% 
  ungroup()

fibtest %>% 
  count(fib_test)

# Create an Africa Map

theme_set(theme_bw())

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)


world <- ne_countries(scale = "medium", 
                      returnclass = "sf")

africa <- world %>% 
  filter(continent == "Africa") %>% 
  mutate(name = if_else(name == "Côte d'Ivoire",
                        "Ivory Coast", name))


population <- data %>% 
  transmute(name = country, studypop2, n_chb) %>% 
  add_count(name, name = "cohorts") %>% 
  group_by(name, studypop2, cohorts) %>% 
  summarise(n_chb = sum(n_chb)) %>% 
  group_by(name) %>% 
  mutate(tot = sum(n_chb)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = studypop2, 
              values_from = n_chb) %>% 
  janitor::clean_names() %>% 
  rename(n_chb = tot)

names <- population %>% 
  count(name) %>% pull(name)

pop_africa <- africa %>% 
  full_join(population)

pop_africa %>% 
  filter(name %in% names) %>% 
  select(name, n_chb, geometry)

map <- st_centroid(pop_africa)
map_points <- cbind(pop_africa, st_coordinates(st_centroid(pop_africa$geometry)))
map_points$label <- !is.na(map_points$n_chb)

breaks <- c(900, 600, 300, 100)

overall <- map_points %>% 
  ggplot() + 
  geom_sf(aes(fill = n_chb), na.rm = TRUE, size = 0.2,
          color = "grey50") + 
  scale_fill_carto_c(palette = "Emrld") +
  # scale_fill_viridis_c(option = "magma", trans = "sqrt", na.value = "#FFFCF1", 
  #                      direction = -1, end = 0.9, breaks = breaks,
  #                      begin = 0.2, limits = c(20, 1190)) +
  coord_sf(ylim = c(-40, 40), 
           xlim = c(-20, 50)) + 
  labs(fill = "N", 
       subtitle = "All patients",
       x = NULL, y = NULL) + 
  guides(fill = guide_colorbar(title.position = "top",
                               direction = "vertical",
                               title.hjust = 0.15)) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        legend.position = "right",
        legend.justification = "center",
        legend.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(c(0, 10, 0, 0)))

overall
ggsave("graph/03-africa_overall.pdf", width = 6, height = 8,
       device = cairo_pdf)


primary <- map_points %>% 
  ggplot() + 
  geom_sf(aes(fill = primary_care_general_population), na.rm = TRUE, size = 0.1,
          color = "grey50") + 
  scale_fill_viridis_c(option = "viridis", trans = "sqrt", na.value = "#FFFCF1", 
                       direction = -1, end = 0.95, breaks = breaks,
                       begin = 0.2) +
  coord_sf(ylim = c(-40, 40), 
           xlim = c(-20, 50)) + 
  labs(fill = "N", 
       subtitle = "Primary care and general population",
       x = NULL, y = NULL) + 
  guides(fill = guide_colorbar(title.position = "top",
                               direction = "vertical",
                               title.hjust = 0.15)) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        legend.position = "right",
        legend.justification = "center",
        legend.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(c(0, 10, 0, 0)))

tertiary <- map_points %>% 
  ggplot() + 
  geom_sf(aes(fill = referral_teaching_hospital), na.rm = TRUE, size = 0.1,
          color = "grey50") + 
  scale_fill_viridis_c(option = "cividis", trans = "sqrt", na.value = "#FFFCF1", 
                       direction = -1, end = 0.95, breaks = breaks,
                       begin = 0.2) +
  coord_sf(ylim = c(-40, 40), 
           xlim = c(-20, 50)) + 
  labs(fill = "N", 
       subtitle = "Referral or teaching hospitals",
       x = NULL, y = NULL) + 
  guides(fill = guide_colorbar(title.position = "top",
                               direction = "vertical",
                               title.hjust = 0.15)) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        legend.position = "right",
        legend.justification = "center",
        legend.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        plot.margin = margin(c(10, 10, 0, 0)))

library(patchwork)

patchwork <- overall | (primary / tertiary + 
                          plot_layout(heights = c(1, 1)))
patchwork + 
  plot_layout(width = c(1, 0.45), height = c(1, 1))+
  plot_annotation(title = "Distribution of all patients included in the review\n") &
  theme(text = element_text(family = "Lato"),
        plot.title = element_text(face = "bold"))

ggsave("graph/03-africa_combined.pdf", width = 7, height = 5, device = cairo_pdf)
ggsave("graph/03-africa_combined.png", width = 7, height = 5, dpi = 300)
embed_fonts("graph/03-africa_combined.pdf")






# Another Version of the map with labels and a barchart next to it

africamap <- map_points %>% 
  ggplot() + 
  geom_sf(aes(fill = n_chb), na.rm = TRUE, size = 0.2,
          color = "grey50") + 
  geom_point(aes(X, Y), color = "black", size = 1, 
             filter(map_points, !is.na(n_chb))) +
  geom_label_repel(aes(label = name,
                       x = X, 
                       y = Y), 
                   filter(map_points, !is.na(n_chb)),
                   family = "Roboto Condensed",
                   box.padding = 0.5,
                   size = 4,
                   direction = "both",
                   alpha = 0.7) + 
  scale_fill_distiller(palette = "YlOrRd", na.value = "white", 
                       trans = "log10", limits = c(20, 1190), direction = 1) +
  coord_sf(ylim = c(-40, 40), 
           xlim = c(-20, 50)) + 
  labs(fill = "Patients included",
       x = NULL, y = NULL) + 
  guides(fill = guide_colorbar(title.position = "top",
                               direction = "horizontal",
                               title.hjust = 0.5,
                               barheight = 0.7,
                               title.theme = element_text(family = "Roboto Condensed",
                                                          size = 12),
                               label.theme = element_text(family = "Roboto Condensed", 
                                                          size = 9))) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        
        legend.position = c(0.15, 0.08),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(c(0, 10, 0, 0)))


bardata <- data %>% 
  select(country, n_chb, group, studypop2) %>% 
  group_by(country, studypop2) %>% 
  summarise(n = sum(n_chb)) %>% 
  group_by(country) %>% 
  mutate(n_country = sum(n)) %>% 
  ungroup() %>% 
  mutate(country = fct_reorder(country, n_country)) %>%
  arrange(n) 


barchart <- bardata %>% 
  mutate(studypop2 = fct_rev(studypop2)) %>% 
  ggplot(aes(y = country, x = n)) +
  geom_col(aes(fill = studypop2), width = 0.5, alpha = 1) + 
  geom_text(aes(x = n_country, label = n_country), 
            nudge_x = 75, 
            family = "Roboto Condensed") +
  labs(x = "Number of patients",
       y = NULL, 
       fill = "Type of Cohort Setting") + 
  scale_fill_brewer(palette = 5) + 
  scale_x_continuous(position = "bottom") + 
  theme_minimal(base_family = "Roboto Condensed") + 
  theme(legend.position = c(0.45, 0.1),
        legend.justification = "left",
        panel.grid = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white", color = "NA", size = 0.2),
        legend.text = element_text(margin = margin(2, 0, 2, 0, unit = "mm"),
                                   size = 10),
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size = 13),
        panel.grid.major.x = element_line(color = "grey60", size = 0.2, linetype = 1),
        panel.grid.minor.x = element_line(color = "grey80", size = 0.1, linetype = 1),
        plot.margin = margin(10, 10, 10, 10, unit = "mm"))

barchart
library(patchwork)  

africamap + barchart +
  plot_layout(widths = c(2, 1)) + 
  plot_annotation(tag_levels = "A")

ggsave("graph/03-africa_overall_studypop2.png", width = 12, height = 9, dpi = 300)
ggsave("graph/03-africa_overall_studypop2.pdf", width = 12, height = 9,
       device = cairo_pdf)
embed_fonts("graph/03-africa_overall_studypop2.pdf")





## Another version in BW

a <- map_points %>% 
  ggplot() + 
  geom_sf(aes(fill = n_chb), na.rm = TRUE, size = 0.3,
          color = "grey50") + 
  scale_fill_distiller(palette = "YlOrRd", na.value = "white", 
                       trans = "log10", limits = c(20, 1190), direction = 1) +
  coord_sf(ylim = c(-40, 40), 
           xlim = c(-20, 50)) + 
  labs(fill = "HBV Patients",
       x = NULL, y = NULL) + 
  guides(fill = guide_colorbar(title.position = "top",
                               direction = "horizontal",
                               title.hjust = 0.5,
                               barheight = 0.7,
                               title.theme = element_text(family = "Roboto Condensed",
                                                          size = 12),
                               label.theme = element_text(family = "Roboto Condensed", 
                                                          size = 9))) +
  theme(panel.background = element_rect(fill = "grey97"),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.15, 0.08),
        legend.background = element_rect(fill = "white", color = NA, size = 0.2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(c(0, 10, 0, 0)))

a_label <- a +
  geom_label_repel(aes(label = name,
                       x = X, 
                       y = Y), 
                   filter(map_points, !is.na(n_chb)),
                   family = "Roboto Condensed",
                   size = 4,
                   direction = "both",
                   alpha = 0.8, label.size = NA,min.segment.length = 0)


b <- bardata %>% 
  mutate(studypop2 = fct_rev(studypop2)) %>% 
  ggplot(aes(y = country, x = n)) +
  geom_col(aes(fill = studypop2), width = 0.5, color = "white", size = 0.5) + 
  geom_text(aes(x = n_country, label = n_country), 
            nudge_x = 75, 
            family = "Roboto Condensed") +
  labs(x = "Number of HBV patients",
       y = NULL, 
       fill = "Cohort Setting") + 
  scale_fill_manual(values = c("grey85", "grey65"))+
  theme_minimal(base_family = "Roboto Condensed") + 
  theme(legend.position = c(0.45, 0.1),
        legend.justification = "left",
        panel.grid = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white", color = NA, size = 0.2),
        legend.text = element_text(margin = margin(2, 0, 2, 0, unit = "mm"),
                                   size = 10),
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size = 13),
        plot.margin = margin(10, 10, 10, 10, unit = "mm"))

a_label+b+
  plot_layout(widths = c(2.5, 1)) + 
  plot_annotation(tag_levels = "A")


ggsave("graph/03-africa_overall_studypop2_bw_label.png", width = 13, height = 9, dpi = 300)
ggsave("graph/03-africa_overall_studypop2_bw_label.pdf", width = 13, height = 9,
       device = cairo_pdf)
embed_fonts("graph/03-africa_overall_studypop2_bw_label.pdf")

# Remove labels
a+b+
  plot_layout(widths = c(2.5, 1)) + 
  plot_annotation(tag_levels = "A")


ggsave("graph/03-africa_overall_studypop2_bw.png", width = 13, height = 9, dpi = 300)
ggsave("graph/03-africa_overall_studypop2_bw.pdf", width = 13, height = 9,
       device = cairo_pdf)
embed_fonts("graph/03-africa_overall_studypop2_bw.pdf")


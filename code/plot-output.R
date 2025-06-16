# Plot simulation output

source(here("code", "utils.R"))

### Set up ###

date <- "2025-06-10"
figure_path <- here("figures", "simulation", date)
output_path <- here("output", "simulation", date)
ifelse(!dir.exists(figure_path), dir.create(figure_path, recursive = TRUE), FALSE)

## Load in inhib values

inhib_table <- qread(here("data", "wmel-inhib.qs")) |>
  mutate(mosquito = "wmel") |> 
  bind_rows(qread(here("data", "walbb-inhib.qs")) |> mutate(mosquito = "walbb"))

summary_inhib <- inhib_table |> 
  group_by(mosquito) |> 
  summarise(mean_inhib = mean(dissemination_norm, na.rm = TRUE), 
            median_inhib = median(dissemination_norm, na.rm = TRUE),
            q_0.025 = quantile(dissemination_norm, 0.025, na.rm = TRUE),
            q_0.975 = quantile(dissemination_norm, 0.975, na.rm = TRUE)) 


brazil_demog <- load_demography()

## Load in results 

files <- list.files(here("output", "simulation", date), pattern = "^time-to-emergence.*\\.csv$", full.names = TRUE)
reemergence_times <- map_dfr(files, read_csv)


reemergence_quantile <- reemergence_times |>
  group_by(inhib_level, r0) |> 
  summarise(median = median(time_to_reemergence/365, na.rm = TRUE),
            q_0.025 = quantile(time_to_reemergence/365, 0.025, na.rm = TRUE),
            q_0.975 = quantile(time_to_reemergence/365, 0.975, na.rm = TRUE), 
            n_reemerged = sum(!is.na(time_to_reemergence)),
            p_20 = sum(time_to_reemergence <= 20*365)/n(),
            p_20 = case_when(is.na(p_20) ~ 0, T ~ p_20)) 


#files <- list.files(here("output", "simulation", date), pattern = "^strain-proportions.*\\.csv$", full.names = TRUE)
#strain_proportions <- map_dfr(files, read_csv)

strain_proportions <- read.csv(here("output", "simulation", date, "strain-proportions_wmel.csv")) |> 
  mutate(mosquito = "wmel") |>
  bind_rows(read.csv(here("output", "simulation", date, "strain-proportions_walbb.csv")) |> 
              mutate(mosquito = "walbb")) |> 
  left_join(inhib_table |>  select(virus, dissemination_norm, mosquito), by = c("virus", "mosquito")) |> 
  mutate(year_sim = floor(time/365),
        date = time/365 + as.Date("2025-01-01"), 
        year = year(date)) |> 
  left_join(brazil_demog |> group_by(year) |> summarise(population = sum(population)), by = "year") |> 
  mutate(inc = value/population * 100000, 
         is_circulating = case_when(inc >= 0.1 ~ 1,
                                    T ~ 0)) |>
  group_by(time, mosquito, draw, run) |> 
  mutate(total_circulating = sum(is_circulating)) |> 
  ungroup() |> 
  group_by(time, mosquito) |> 
  mutate(median_circulating = median(total_circulating)) |>
  ungroup() 

example_dynamics <- qread(here("output", "strain-simulation", date, "strain-quantile-out_2.qs")) |> 
  filter(name == "inf")

# Mean and probability of reemergence < 20 years

wol_statistics <- tibble(r0 = NA, mosquito = NA, p_20_dist = list(), p_20_mean = NA)

for(mosquito in c("wmel", "walbb")){ 
for(r0_input in unique(reemergence_quantile$r0)){
  temp <- tibble(r0 = NA, mosquito = NA, p_20_dist = list(), p_20_mean = NA)
  df_reemergence <- reemergence_quantile |> 
    filter(r0 == as.double(r0_input))
  
  # mean time to reemegence
  
  temp[1, "r0"] <- r0_input
  temp[1, "mosquito"] <- mosquito
  
  # probability of reemergence < 20 years
  
  prob_fun <- cinterpolate::interpolation_function(df_reemergence$inhib_level, df_reemergence$p_20, type = "linear")
  temp$p_20_dist[[1]] <- prob_fun(inhib_table$dissemination_norm[inhib_table$mosquito == mosquito])
  temp[1, "p_20_mean"] <- mean(temp$p_20_dist[[1]])
  wol_statistics <- rbind(wol_statistics, temp)
  rm(temp)
  }
}

### Plotting ###

## Panel A - time to reemergence for different R0 and inhibition combinations

## Probability reemergence column

wol_long <- wol_statistics |>  
  bind_rows(add_0) |> 
  unnest_longer(p_20_dist) |> 
  mutate(mosquito = factor(mosquito , levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel")))

 probability_reemergence <- wol_long |> 
   filter(r0 >= 1) |> 
   mutate(r0 = factor(r0, levels = seq(5, 1, by = -0.5))) |> 
  ggplot() +
  geom_density(data = filter(wol_long, r0 != "0.5" & r0 != "0"), aes(x = p_20_dist, group = mosquito, fill = mosquito), alpha = 0.5) + 
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_fill_manual(values = c("#FED789FF", "#72874EFF")) +
  facet_grid(rows = vars(r0), switch = "y", scale = "free_y") +
  #scale_x_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.25)) +
  labs(x = "Probability of\nreemergence ≤ 20 years", y = NULL, fill = "Mosquito strain") +
   theme(strip.text.y.left = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.line.y.left = element_blank(),
         legend.position = "none") +
   theme(plot.margin = margin(10, 10, 20, 0))
 
  
 ggsave(here(figure_path, "prob-reemergence.png"), plot = probability_reemergence, width = 180, height = 150, unit = "mm", dpi = 300)
 

inhib_density <- inhib_table |> 
  mutate(mosquito = factor(mosquito, levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel"))) |> 
  ggplot() +
  geom_density(aes(x = dissemination_norm, fill = mosquito), alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,1, by = 0.1), limits = c(0,1), expand = c(0,0)) +
  scale_fill_manual(values = c("#FED789FF", "#72874EFF")) +
  theme(legend.position = "right")  +
  labs(x = NULL, fill = "Mosquito strain", y = "Density") +
  theme(plot.margin = margin(10, 0, 10, 15),
        axis.text.y = element_blank())


reemergence_plot <- reemergence_quantile |> 
  filter(r0 >= 1) |> 
  mutate(median = case_when(is.na(median) ~ 75, T~ median)) |> 
  ggplot() +
  geom_tile(aes(x = inhib_level, y = r0, fill = median)) +
  scale_y_continuous(breaks = seq(1, 5, by = 0.5), limits = c(0.75, 5.25)) +
  labs(x = "Relative transmissibility\n", y = "R0", fill = "Time to\nreemergence\n(years)") +
  scale_fill_gradient(low = "#023743FF",
                      high = "#fff7fb",
                      na.value = "#969696", breaks = c(0, 15, 30, 45, 60, 75), 
                      labels = c("0", "15", "30", "45","60","75+")) +
  scale_x_continuous(breaks = seq(0,1, by = 0.1), limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar(title.position = "left", title.hjust = 1, title.vjust = 0.8)) +
  theme(plot.margin = margin(10, 0, 20, 10), 
        legend.justification = "left",
        legend.spacing.x = unit(-0.3, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.spacing = unit(0, "pt")) 

legend <- ggpubr::get_legend(reemergence_plot)
panel_a_lower <- plot_grid(reemergence_plot + theme(legend.position = "none"), probability_reemergence, rel_widths = c(4,1))
panel_a_lower <- plot_grid(panel_a_lower, legend, ncol = 1, rel_heights = c(1, 0.1), labels = c("B", ""), label_size = 10, label_fontfamily = plot_font)
panel_a_upper <- plot_grid(inhib_density, legend_inhib, rel_widths = c(40,1), ncol = 2, nrow = 1, labels = c("A", ""), label_size = 10, label_fontfamily = plot_font)
panel_a <- plot_grid(panel_a_upper, panel_a_lower, rel_heights = c(1,5), ncol = 1, label_size = 10, label_fontfamily = plot_font)

ggsave(here(figure_path, "panel_a.png"), plot = panel_a, width = 180, height = 150, unit = "mm", dpi = 300)

## Panel C - selection over time

#pal <- c("#435F90", "#7c4b73", "#ed968c", "#e78429") # from Archaumbault in met brewer palette

panel_labels <- strain_proportions |> 
  group_by(year_sim, mosquito, median_circulating) |> 
  summarise(median_circulating = unique(median_circulating)) |> 
  ungroup() |> 
  mutate(year_sim = factor(year_sim, levels = c(0,5,20,50,75), 
                           labels = c("Pre-release", "5 years\npost-release", "20 years\npost-release", 
                                      "50 years\npost-release", "75 years\npost-release"))) 

panel_c <- strain_proportions |> 
  group_by(name, time, draw, mosquito, virus, dissemination_norm) |> 
  mutate(median = median(proportions)) |> 
  mutate(year_sim = factor(year_sim, levels = c(0,5,20,50,75), 
         labels = c("Pre-release", "5 years\npost-release", "20 years\npost-release", 
                    "50 years\npost-release", "75 years\npost-release"))) |>
  ggplot() +
 geom_violin(aes(factor(year_sim), proportions, fill = mosquito), position = position_dodge(width = 0.5)) +
  geom_point(aes(factor(year_sim), median, col = dissemination_norm, group = mosquito), position = position_dodge(width = 0.5), size = 0.5, alpha = 0.3, width = 0.01) +
  geom_text(data = panel_labels, aes(x = factor(year_sim), y = 1.0, group = mosquito, label = median_circulating), 
            position = position_dodge(width = 0.5),
            vjust = -0.5) +
  scale_color_gradient(low = "#7640A9FF", 
                      high = "#F8B150FF") +
  scale_y_continuous(breaks = seq(0,1, by = 0.25), limits = c(0,1), labels = c("Absent", "0.25", "0.5", "0.75", "1.00")) +
  labs(x = "", y = "Strain prevalence", col = "Relative transmissibility") +
  ggtitle("wMel releases") +
  theme(plot.title = element_text(family = plot_font, face = "bold", size = 9))

ggsave(here(figure_path, "panel_c.png"), plot = panel_b, width = 180, height = 150, unit = "mm", dpi = 300)

figure <- plot_grid(panel_a, panel_c, nrow = 1, rel_widths = c(1, 0.6), labels = c("", "C"), label_size = 10, label_fontfamily = plot_font)
ggsave(here(figure_path, "modelling-fig.png"), plot = figure, width = 240, height = 120, unit = "mm", dpi = 300, bg = "white")

## Panel C - example dynamics

dynamics_fig <- example_dynamics |> 
  mutate(date = time + as.Date("2025-01-01") - 1) |> 
  filter(date <= "2100-01-01") |>
  ggplot() +
  geom_line(aes(x = date, y = median, col = as.factor(level), group = as.factor(level))) +
  geom_ribbon(aes(x = date, ymin = q_0.025, ymax = q_0.975, fill = as.factor(level), group = as.factor(level)), alpha = 0.2) +
  labs(y = "Infections", x = "Date", fill = "Strain", col = "Strain")

ggsave(here(figure_path, "example-dynamics-75-years.png"), plot = dynamics_fig, width = 180, height = 150, unit = "mm", dpi = 300)

# Plot simulation output

source(here("code", "utils.R"))

### Set up ###

date <- "2025-06-25"

figure_path <- here("figures")
output_path <- here("output", "simulation", date)
ifelse(!dir.exists(figure_path), dir.create(figure_path, recursive = TRUE), FALSE)

## Load in inhib values

inhib_table <- qread(here("data", "wmel-inhib.qs")) |>
  mutate(mosquito = "wmel") |> 
  bind_rows(qread(here("data", "walbb-inhib.qs")) |> mutate(mosquito = "walbb")) |> 
  mutate(inhib_level = round(dissemination_norm,3))

summary_inhib <- inhib_table |> 
  group_by(mosquito) |> 
  summarise(mean_inhib = mean(dissemination_norm, na.rm = TRUE), 
            median_inhib = median(dissemination_norm, na.rm = TRUE),
            q_0.025 = quantile(dissemination_norm, 0.025, na.rm = TRUE),
            q_0.975 = quantile(dissemination_norm, 0.975, na.rm = TRUE)) 


brazil_demog <- load_demography()

## Load in results 

# Reemergence results
files <- list.files(here("output", "simulation", date), pattern = "^time-to-emergence.*\\.csv$", full.names = TRUE)
files_strain <- list.files(here("output", "simulation", date, "strain-specific"), pattern = "^time-to-emergence.*\\.csv$", full.names = TRUE)

reemergence_times <- map_dfr(c(files), read_csv) |> 
  mutate(type = "grid") |> 
  mutate(diff = time_to_reemergence - time_to_extinction) 

reemergence_quantile <- reemergence_times |>
  group_by(inhib_level, r0, serotype_varying, type) |> 
  summarise(median = median(time_to_reemergence/365, na.rm = TRUE),
            q_0.025 = quantile(time_to_reemergence/365, 0.025, na.rm = TRUE),
            q_0.975 = quantile(time_to_reemergence/365, 0.975, na.rm = TRUE))

reemergence_probs <- map_dfr(c(files_strain), read_csv) |> mutate(type = "strain") |> 
  group_by(inhib_level, r0, serotype_varying, type) |> 
  summarise(n_reemerged_sub_20 = sum(time_to_reemergence <= 20*365, na.rm = TRUE),
            n_reemerged_sub_40 = sum(time_to_reemergence <= 40*365, na.rm = TRUE),
            n_total = n(), 
            p_20 = n_reemerged_sub_20/n_total,
            p_40 = n_reemerged_sub_40/n_total, 
            p_20 = case_when(is.na(p_20) ~ 0, T ~ p_20),
            p_40 = case_when(is.na(p_40) ~ 0, T ~ p_40))

# Create data frame for each R0 and mosquito

reemergence_strains <- map_dfr(seq(1,5, by = 0.5), ~ inhib_table |>  mutate(r0 = .x) |> rename(serotype_varying = serotype)) |> 
  left_join(reemergence_probs, by = c("inhib_level", "r0", "serotype_varying")) |> 
  mutate(r0 = factor(r0, levels = seq(5, 1, by = -0.5)),
        mosquito = factor(mosquito, levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel")),
        serotype_varying = factor(serotype_varying, levels = c(1, 2, 3, 4), 
                                         labels = c("DENV-1", "DENV-2", "DENV-3", "DENV-4")))

# Strain selection results

#files <- list.files(here("output", "simulation", date), pattern = "^strain-out.*\\.csv$", full.names = TRUE)
#strain_out <- map_dfr(files, read_csv)

strain_out <- read.csv(here("output", "simulation", date, "strain-out_wmel.csv")) |> 
  mutate(mosquito = "wmel") |>
  bind_rows(read.csv(here("output", "simulation", date, "strain-out_walbb.csv")) |> 
              mutate(mosquito = "walbb")) |> 
  left_join(inhib_table |>  select(virus, dissemination_norm, mosquito), by = c("virus", "mosquito")) |> 
  mutate(year_sim = floor(time/365),
        date = time + as.Date("2025-01-01"), 
        year = year(date)) |> 
  left_join(brazil_demog |> group_by(year) |> summarise(population = sum(population)), by = "year") |> 
  mutate(inc = value/population * 100000, 
         is_circulating = case_when(inc >= 0.1 ~ 1,
                                    T ~ 0),
         weight_dissemination = proportions * dissemination_norm) |> 
  group_by(time, date, mosquito, draw, run) |> 
  mutate(total_circulating = sum(is_circulating),
         abundance_weighted_dissem = sum(weight_dissemination)) |> 
  ungroup() |> 
  group_by(time, date, year_sim, mosquito) |>
  mutate(median_circulating = median(total_circulating))

strain_plot <- strain_out |>
  group_by(time, date, mosquito, draw, run) |> 
  mutate(total_circulating = sum(is_circulating),
         abundance_weighted_dissem = sum(weight_dissemination)) |> 
  ungroup() |> 
  group_by(time, date, year_sim, mosquito) |>
  summarise(median_circulating = median(total_circulating), 
            q_0.025_circulating = quantile(total_circulating, 0.025), 
            q_0.975_circulating = quantile(total_circulating, 0.975),
            median_transmissibility = median(abundance_weighted_dissem),
            q_0.025_transmissibility = quantile(abundance_weighted_dissem, 0.025),
            q_0.975_transmissibility = quantile(abundance_weighted_dissem, 0.975)) |>
  ungroup() 

example_dynamics <- qread(here("output", "simulation", date, "quantiles-out_inhib-level-0.75.qs"))

### Plotting ###

## Panel A - time to reemergence for different R0 and inhibition combinations

## Probability reemergence column
    
reemergence_counts <- reemergence_strains |> 
  group_by(serotype_varying, mosquito, r0, p_40) |> 
  summarise(count = n())
    
probability_reemergence <- reemergence_strains |> 
    group_by(serotype_varying, mosquito,  r0) |> 
    summarise(
      mean = mean(p_40), 
      median = median(p_40),
      q_0.25 = quantile(p_40, 0.25), 
      q_0.75 = quantile(p_40, 0.75),
      q_0.025 = quantile(p_40, 0.025), 
      q_0.975 = quantile(p_40, 0.975),
      .groups = "drop")  |> 
    ggplot(aes(x = mosquito, y = median, group = mosquito)) +
    geom_rect(aes(xmin = as.numeric(factor(mosquito)) - 0.3,
                  xmax = as.numeric(factor(mosquito)) + 0.3,
                  ymin = q_0.25, ymax = q_0.75, fill = mosquito),
              alpha = 0.5) +
    geom_point(data = reemergence_counts, aes(x = mosquito, y = p_40, color = mosquito, size = count), shape = 20, alpha = 0.9) +
    geom_segment(aes(x = as.numeric(factor(mosquito)) - 0.3,
                     xend = as.numeric(factor(mosquito)) + 0.3,
                     y = median, yend = median),
                 color = "black", linewidth = 0.5) +
    facet_grid(r0 ~ serotype_varying, scales = "free_y") +
    coord_flip() +
    labs(y = "Probability of\nreemergence ≤ 40 years", x = NULL, fill = "Mosquito strain", size = "Number of\nDENV strains") +
    scale_fill_manual(values = c("#C5692D", "#132F5B"), name = "Mosquito strain") +
    scale_color_manual(values = c("#C5692D", "#132F5B"), name = "Mosquito strain") +
    scale_size_continuous(range = c(0.5,2.2)) +
    scale_y_continuous(breaks = c(0,0.5,1),labels = c("0", "0.5", "1")) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y.left = element_blank()) +
    guides(size = guide_legend(title.position = "left", title.hjust = 1, title.vjust = 0.8), color = "none", fill = "none") +
    theme(plot.margin = margin(15, 10, 20, 5),
          legend.justification = "right") 
  
    ggsave(here(figure_path, "panel_c_probability-reemergence_40-years.png"), plot = probability_reemergence, width = 180, height = 120, unit = "mm", dpi = 300)
    
panel_a <- inhib_table |> 
  mutate(mosquito = factor(mosquito, levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel"))) |> 
  ggplot() +
  geom_histogram(aes(x = dissemination_norm, fill = mosquito),position = "identity", alpha = 0.5, binwidth = 0.05) +
  scale_fill_manual(values = c("#C5692D", "#132F5B")) +
  scale_x_continuous(breaks = seq(0,1, by = 0.1), limits = c(-0.04,1.02), expand = c(0,0), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")) +
  theme(legend.position = "right")  +
  labs(x = NULL, fill = "Mosquito strain", y = "Frequency") +
  guides(fill = guide_legend(direction = "vertical")) +
  theme(plot.margin = margin(10, 5, 10, 15),
              legend.key.width = unit(0.6, "cm"),
              legend.key.height = unit(0.6, "cm"),
              legend.title = element_text(size = 9, family = plot_font),
              legend.position = c(1, 1),
              legend.justification = c(1, 1))

ggsave(here(figure_path, "panel_a.png"), plot = panel_a, width = 180, height = 120, unit = "mm", dpi = 300)

legend_inhib <- ggpubr::get_legend(inhib_density)

reemergence_plot <- reemergence_quantile |> 
  filter(r0 >= 1) |> 
  mutate(median = case_when(is.na(median) ~ 75, T~ median)) |> 
  ggplot() +
  geom_tile(aes(x = inhib_level, y = r0, fill = median)) +
  labs(x = "Relative dissemination\n", y = "R0", fill = "Time to\nreemergence\n(years)") +
  scale_fill_gradient(low = "#3A5F44",
                      high = "#fff7fb",
                      na.value = "#969696", breaks = c(0, 15, 30, 45, 60, 75), 
                      labels = c("0", "15", "30", "45","60","75+")) +
  scale_x_continuous(breaks = seq(0,1, by = 0.1), limits = c(-0.04,1.02), expand = c(0,0), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")) +
scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar(title.position = "left", title.hjust = 1, title.vjust = 0.8)) +
  theme(plot.margin = margin(30, 5, 20, 17.5), 
        legend.justification = "left",
        legend.spacing.x = unit(-0.3, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.spacing = unit(0, "pt"),
        legend.title = element_text(size = 8, family = plot_font))

legend <- ggpubr::get_legend(reemergence_plot)
legend_prob <- ggpubr::get_legend(probability_reemergence)
ggsave(here(figure_path, "panel-c_reemergence-tile.png"), plot = reemergence_plot, width = 240, height = 180, unit = "mm", dpi = 300, bg = "white")

## Panel C - selection over time

#pal <- c("#435F90", "#7c4b73", "#ed968c", "#e78429") # from Archaumbault in met brewer palette

plot_times <- c(1,182, 2007, 3832, 7482, 18432, 27192)

panel_labels <- strain_out |> 
  group_by(time, date, mosquito, median_circulating) |> 
  summarise(median_circulating = unique(median_circulating)) |> 
  filter(time %in% plot_times) |> 
  mutate(time_label = factor(time, levels = plot_times, labels = c("Pre-release", "6\nmonths", "5\nyears", "10\nyears", "20\nyears", "50\nyears", "75\nyears"))) |> 
  ungroup() 


# Panel B

  panel_b <- strain_plot |> 
  mutate(mosquito = factor(mosquito, levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel"))) |> 
  filter(time %% 365 != 0 ) |> 
  ggplot() +
  geom_line(aes(x = date, y = median_transmissibility, col = mosquito, group = mosquito), lty = "dashed") +
  geom_ribbon(aes(x = date, ymin = q_0.025_transmissibility, ymax = q_0.975_transmissibility, fill = mosquito, group = mosquito), alpha = 0.2) +
    scale_fill_manual(values = c("#C5692D", "#132F5B"), name = "Mosquito strain") +
    scale_color_manual(values = c("#C5692D", "#132F5B"), name = "Mosquito strain") +
  scale_y_continuous(name = "Relative dissemination\n of circulating strains") +
  scale_x_date(labels = label_date_short(), breaks = seq(as.Date("2025-01-01"), as.Date("2100-02-01"), by = "10 years")) +
  labs(x = "Date", col = NULL, fill = NULL) + 
  theme(legend.position = "none")
  
  ggsave(here(figure_path, "panel_b.png"), plot = panel_b, width = 180, height = 120, unit = "mm", dpi = 300)
  
  top_row <- plot_grid(panel_a, panel_b, labels = c("A", "B"), rel_widths = c(2,1),label_size = 10, label_fontfamily = plot_font)
  panel_c <- plot_grid(reemergence_plot + theme(legend.position = "none"), probability_reemergence + theme(legend.position = "none"), rel_widths = c(2,1), labels = c("C", "D"), label_size = 10, label_fontfamily = plot_font)
  legends <- plot_grid(legend, legend_prob, ncol = 2)
  panel_c <- plot_grid(panel_c, legends, ncol = 1, rel_heights = c(1, 0.1))

  
  figure <- plot_grid(top_row, panel_c, nrow = 2, rel_heights = c(1,2.5), labels = c("", "C"), label_size = 10, label_fontfamily = plot_font)
  ggsave(here(figure_path, "modelling-fig.png"), plot = figure, width = 240, height = 180, unit = "mm", dpi = 300, bg = "white")
  
## Panel D - example dynamics

dynamics_fig <- example_dynamics |> 
  filter(name == "inf" & serotype_varying == 4 & r0 == 1) |> 
  mutate(date = time + as.Date("2025-01-01") - 1) |> 
  filter(date <= as.Date("2034-01-01")) |> View()
  ggplot() +
  geom_line(aes(x = date, y = median/population*100000, col = as.factor(level), group = as.factor(level))) +
  geom_ribbon(aes(x = date, ymin = q_0.025/population*100000, ymax = q_0.975/population*100000, fill = as.factor(level), group = as.factor(level)), alpha = 0.2) +
 # coord_cartesian(ylim = c(0,10000)) +
  labs(y = "Infections", x = "Date", fill = "Strain", col = "Strain")

ggsave(here(figure_path, "example-dynamics.png"), plot = dynamics_fig, width = 180, height = 120, unit = "mm", dpi = 300)

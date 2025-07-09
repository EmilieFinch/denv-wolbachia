# Plot simulation output

source(here("code", "utils.R"))

### Set up ###

date <- "2025-06-19"

figure_path <- here("figures")
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

# Serotype results
files <- list.files(here("output", "simulation", date), pattern = "^time-to-emergence.*\\.csv$", full.names = TRUE)
reemergence_times <- map_dfr(files, read_csv)


reemergence_quantile <- reemergence_times |>
  group_by(inhib_level, r0) |> 
  summarise(median = median(time_to_reemergence/365, na.rm = TRUE),
            q_0.025 = quantile(time_to_reemergence/365, 0.025, na.rm = TRUE),
            q_0.975 = quantile(time_to_reemergence/365, 0.975, na.rm = TRUE), 
            n_reemerged = sum(!is.na(time_to_reemergence)),
            n_reemerged_sub_20 = sum(time_to_reemergence <= 20*365, na.rm = TRUE),
            n_total = n(),
            p_20 = n_reemerged_sub_20/n_total,
            p_20 = case_when(is.na(p_20) ~ 0, T ~ p_20)) |> 
  ungroup()

# Mean and probability of reemergence < 20 years

wol_statistics <- tibble(r0 = NA, mosquito = NA,p_20_dist = list(), p_20_mean = NA)

for(mosquito_input in c("wmel", "walbb")){ 
  for(r0_input in unique(reemergence_quantile$r0)){
      temp <- tibble(r0 = NA, mosquito = NA, p_20_dist = list(), p_20_mean = NA)
      df_reemergence <- reemergence_quantile |> 
        filter(r0 == as.double(r0_input)) 

      inhib_in <- inhib_table |> 
        filter(mosquito == mosquito_input) |> 
        mutate(dissemination_norm = case_when(dissemination_norm == 0 ~ dissemination_norm + 1e-10, 
                                              T ~ dissemination_norm))
      
      # mean time to reemegence
      
      temp[1, "r0"] <- r0_input
      temp[1, "mosquito"] <- mosquito_input

      # probability of reemergence < 20 years
      
      prob_fun <- cinterpolate::interpolation_function(df_reemergence$inhib_level, df_reemergence$p_20, type = "linear")
      temp$p_20_dist[[1]] <- prob_fun(inhib_in$dissemination_norm)
      temp[1, "p_20_mean"] <- mean(temp$p_20_dist[[1]])
      wol_statistics <- rbind(wol_statistics, temp)
      rm(temp)
    }
  }

# Multi-strain results
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

#example_dynamics <- qread(here("output", "simulation", date, "strain-quantile-out_wmel_2.qs")) |> 
# filter(name == "inf")


### Plotting ###

## Panel A - time to reemergence for different R0 and inhibition combinations

## Probability reemergence column

wol_long <- wol_statistics |>  
  unnest_longer(p_20_dist) |> 
  mutate(mosquito = factor(mosquito , levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel"))) |> 
  group_by(r0, mosquito) |> 
  mutate(
    has_variation = var(p_20_dist, na.rm = TRUE) > 0) 

probability_reemergence <- wol_long |> 
  mutate(r0 = factor(r0, levels = seq(5, 1, by = -0.5))) |> 
  ggplot() +
  geom_density(data = wol_long |> filter(has_variation), aes(x = p_20_dist, group = mosquito, fill = mosquito), alpha = 0.5) + 
  geom_vline(aes(xintercept = p_20_mean, group = mosquito, col = mosquito), lty = "longdash") +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_fill_manual(values = c("#FED789FF", "#72874EFF")) +
  scale_color_manual(values = c("#FED789FF", "#72874EFF")) +
  facet_grid(rows = vars(r0), switch = "y", scale = "free_y") +
  #scale_x_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.25)) +
  labs(x = "Probability of\nreemergence ≤ 20 years", y = NULL, fill = "Mosquito strain") +
  theme(strip.text.y.left = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.left = element_blank(),
        legend.position = "none") +
  theme(plot.margin = margin(10, 10, 20, 0)) 

ggsave(here(figure_path, "prob-reemergence_20-years.png"), plot = probability_reemergence, width = 180, height = 150, unit = "mm", dpi = 300)

box_plots <- wol_long |> 
  ggplot() +
  geom_boxplot(aes(y = p_20_dist, group = interaction(serotype_varying,mosquito), fill = mosquito)) +
  facet_grid(r0 ~ serotype_varying, scales = "free_y")  
#ggsave(here(figure_path, "prob-reemergence.png"), plot = probability_reemergence, width = 180, height = 150, unit = "mm", dpi = 300)


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
ggsave(here(figure_path, "inhib-density-serotype.png"), plot = inhib_density, width = 180, height = 150, unit = "mm", dpi = 300)


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
panel_a_upper <- plot_grid(inhib_density, legend, rel_widths = c(40,1), ncol = 2, nrow = 1, labels = c("A", ""), label_size = 10, label_fontfamily = plot_font)
panel_a <- plot_grid(panel_a_upper, panel_a_lower, rel_heights = c(1,5), ncol = 1, label_size = 10, label_fontfamily = plot_font)

ggsave(here(figure_path, "panel_a.png"), plot = panel_a, width = 180, height = 150, unit = "mm", dpi = 300)

## Panel C - selection over time

#pal <- c("#435F90", "#7c4b73", "#ed968c", "#e78429") # from Archaumbault in met brewer palette

plot_times <- c(1,182, 2007, 3832, 7482, 18432, 27192)

panel_labels <- strain_out |> 
  group_by(time, date, mosquito, median_circulating) |> 
  summarise(median_circulating = unique(median_circulating)) |> 
  filter(time %in% plot_times) |> 
  mutate(time_label = factor(time, levels = plot_times, labels = c("Pre-release", "6\nmonths", "5\nyears", "10\nyears", "20\nyears", "50\nyears", "75\nyears"))) |> 
  ungroup() 

panel_c <- strain_out |> 
  filter(time %in% plot_times) |>  
  mutate(mosquito = factor(mosquito, levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel")),
         time_label = factor(time, levels = plot_times, labels = c("Pre-release", "6\nmonths", "5\nyears", "10\nyears", "20\nyears", "50\nyears", "75\nyears"))) |> 
  group_by(name, time_label, date, draw, mosquito, virus, dissemination_norm) |> 
  mutate(median = median(proportions)) |> 
  ggplot() +
  geom_point(aes(x = time_label,y = median, col = dissemination_norm, group = mosquito, shape = mosquito), position = position_dodge(width = 0.5), size = 1, alpha = 0.3) +
  geom_text(data = panel_labels, aes(x = time_label, y = 1.0, group = mosquito, label = median_circulating), 
            position = position_dodge(width = 0.5),
            vjust = -0.5,
            family = plot_font,
            size = 2) +
  scale_color_gradient(low = "#7640A9FF", 
                       high = "#F8B150FF") +
  scale_y_continuous(breaks = seq(0,1, by = 0.25), limits = c(0,1), labels = c("Absent", "0.25", "0.5", "0.75", "1.00")) +
  labs(x = "", y = "Strain prevalence", col = "Relative transmissibility", shape = "Mosquito") +
  theme(plot.title = element_text(family = plot_font, face = "bold", size = 9),
        legend.box = "vertical")

ggsave(here(figure_path, "panel_c.png"), plot = panel_c, width = 180, height = 150, unit = "mm", dpi = 300)

# Alternate for Panel C

panel_c_alt <- strain_plot |> 
  mutate(mosquito = factor(mosquito, levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel"))) |> 
  filter(time %% 365 != 0 ) |> 
  ggplot() +
  geom_line(aes(x = date, y = median_circulating, col = mosquito, group = mosquito)) +
  geom_ribbon(aes(x = date, ymin = q_0.025_circulating, ymax = q_0.975_circulating, fill = mosquito, group = mosquito), alpha = 0.2) +
  geom_line(aes(x = date, y = median_transmissibility*25, col = mosquito, group = mosquito), lty = "dashed") +
  geom_ribbon(aes(x = date, ymin = q_0.025_transmissibility*25, ymax = q_0.975_transmissibility*25, fill = mosquito, group = mosquito), alpha = 0.2) +
  scale_color_manual(values = c("#FED789FF", "#72874EFF")) +
  scale_fill_manual(values = c("#FED789FF", "#72874EFF")) +
  scale_y_continuous(name = "Number of circulating strains", 
                     sec.axis = sec_axis(~./25, name = "Relative transmissibility")) +
  scale_x_date(labels = label_date_short(), breaks = seq(as.Date("2025-01-01"), as.Date("2100-02-01"), by = "10 years")) +
  labs(x = "Date", col = NULL, fill = NULL) + 
  theme(legend.position = "none")

ggsave(here(figure_path, "panel_c_alt.png"), plot = panel_c_alt, width = 180, height = 150, unit = "mm", dpi = 300)

figure <- plot_grid(panel_a, panel_c, nrow = 1, rel_widths = c(1, 0.6), labels = c("", "C"), label_size = 10, label_fontfamily = plot_font)
ggsave(here(figure_path, "modelling-fig.png"), plot = figure, width = 240, height = 120, unit = "mm", dpi = 300, bg = "white")

figure <- plot_grid(panel_a, panel_c_alt, nrow = 1, rel_widths = c(1, 0.6), labels = c("", "C"), label_size = 10, label_fontfamily = plot_font)
ggsave(here(figure_path, "modelling-fig_alt.png"), plot = figure, width = 240, height = 120, unit = "mm", dpi = 300, bg = "white")

## Panel D - example dynamics

dynamics_fig <- example_dynamics |> 
  mutate(date = time + as.Date("2025-01-01") - 1) |> 
  filter(date <= "2031-01-01" & date >= "2030-01-01") |>
  ggplot() +
  geom_line(aes(x = date, y = median, col = as.factor(level), group = as.factor(level))) +
  geom_ribbon(aes(x = date, ymin = q_0.025, ymax = q_0.975, fill = as.factor(level), group = as.factor(level)), alpha = 0.2) +
  labs(y = "Infections", x = "Date", fill = "Strain", col = "Strain")

ggsave(here(figure_path, "example-dynamics-75-years.png"), plot = dynamics_fig, width = 180, height = 150, unit = "mm", dpi = 300)

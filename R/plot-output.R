# Plot simulation output

## Set up

date <- "2025-06-10"
figure_path <- here("figures", "simulation", date)
output_path <- here("output", "simulation", date)
ifelse(!dir.exists(figure_path), dir.create(figure_path, recursive = TRUE), FALSE)

## Load in inhib values

inhib_table <- qread(here("data", "wmel-inhib.qs")) |>
  mutate(mosquito = "wmel") |> 
  bind_rows(qread(here("data", "walbb-inhib.qs")) |> mutate(mosquito = "walbb"))

## Load in results 

files <- list.files(here("output", "simulation", date), pattern = "^time-to-emergence.*\\.csv$", full.names = TRUE)
reemergence_times <- map_dfr(files, read_csv)

reemergence_quantile <- reemergence_times |>
  group_by(inhib_level, r0) |> 
  summarise(median = median(time_to_reemergence, na.rm = TRUE),
            q_0.025 = quantile(time_to_reemergence, 0.025, na.rm = TRUE),
            q_0.975 = quantile(time_to_reemergence, 0.975, na.rm = TRUE), 
            n_reemerged = sum(!is.na(time_to_reemergence)))


## Panel A - time to reemergence for different R0 and inhibition combinations

inhib_density <- inhib_table |> 
  mutate(mosquito = factor(mosquito, levels = c("walbb", "wmel"), labels = c("wAlbB", "wMel"))) |> 
  ggplot() +
  geom_density(aes(x = dissemination_norm, fill = mosquito), alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,1, by = 0.1), limits = c(0,1)) +
  scale_fill_manual(values = c("#FED789FF", "#72874EFF")) +
  theme(
    axis.text.y = element_text(color = "white"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.key = element_rect(color = NA)) +
  labs(x = NULL, fill = "Mosquito strain", y = "Density") 
  
reemergence_plot <- ggplot(reemergence_quantile) +
  geom_tile(aes(x = inhib_level, y = r0, fill = median/365)) +
  labs(x = "Level of Wolbachia inhibition", y = "r0", fill = "Time to reemergence (years)") +
  scale_fill_gradient(low = "#023743FF",
                      high = "#D6E1F2FF",
                      na.value = "#969696", breaks = c(0, 10, 20, 30, 40, 50, 60, 70)) +
  scale_x_continuous(breaks = seq(0,1, by = 0.1), limits = c(0,1)) +
  theme(plot.margin = margin(10, 10, 10, 10))

panel_a <- plot_grid(inhib_density, reemergence_plot, ncol = 1, rel_heights = c(1, 3), labels = c("A", "B"), label_size = 10, label_fontfamily = plot_font)

ggsave(here(figure_path, "panel_a.png"), plot = panel_a, width = 180, height = 150, unit = "mm", dpi = 300)

## Panel B - selection over time

pal <- c("#435F90", "#7c4b73", "#ed968c", "#e78429") # from Archaumbault in met brewer palette
panel_b <- strain_proportions |> 
  mutate(year = floor(time/365) + 1) |> 
  ggplot() +
  geom_violin(aes(factor(year), proportions)) +
  geom_point(aes(factor(year), proportions, col = factor(serotype)), size = 1.25, alpha = 0.5) +
  scale_color_manual(values = pal)

ggsave(here(figure_path, "panel_b.png"), plot = panel_b, width = 180, height = 150, unit = "mm", dpi = 300)

figure <- plot_grid(panel_a, panel_b, nrow = 2, rel_heights = c(1.5,1), labels = c("", "C", label_size = 8, label_fontfamily = plot_font))
ggsave(here(figure_path, "modelling-fig.png"), plot = figure, width = 180, height = 150, unit = "mm", dpi = 300)



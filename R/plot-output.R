# Plot simulation output

## Set up

date <- "2025-06-02"
figure_path <- here("figures", "simulation", date)
ifelse(!dir.exists(figure_path), dir.create(figure_path, recursive = TRUE), FALSE)

## Load in results 

quantiles_out <- NULL
strain_div_out <- NULL

r0s <- c(2,4)
mosquitos <- c("wmel", "walbb", "wt")
for(r0 in r0s){
  for(mos in mosquitos){
    temp <- qread(here("output", "simulation", date, paste0("quantile-out_r0-", r0, "_", mos, ".qs")))
    temp2 <- qread(here("output", "simulation", date, paste0("species-div-out_r0-", r0, "_", mos, ".qs")))
    
    quantiles_out <- rbind(quantiles_out, temp)
    strain_div_out <- rbind(strain_div_out, temp2)
  }
}
rm(temp, temp2)

# Plotting

## Infection plot

strain_names <- c("DH002", "DH006", "DH009", "DH014", "DH019", "DH021", "DH022", "DH026", "DH035", "DH040", "DH048",
               "DH131", "DH141", "DH142", "DH143", "DH144", "DH145", "DH146", "DH147", "DH148", "DH060", "DH061", 
               "DH064", "DH069A", "DH071", "DH073", "DH074", "DH077", "DH081", "DH087", "DH091", "DH099", "DH101", 
               "DH134", "DH135", "DH154", "DH156", "DH158", "DH159", "DH161", "DH102", "DH104", "DH105", "DH110", 
               "DH162", "DH164", "DH165", "DH167", "DH168", "DH169", "DH102", "DH104", "DH105", "DH110", "DH162", 
               "DH164", "DH165", "DH167", "DH168", "DH169", "DH114", "DH116", "DH118", "DH119", "DH122", "DH139", 
               "DH140", "DH149", "DH157", "DH170", "DH114", "DH116", "DH118", "DH119", "DH122", "DH139", "DH140", 
               "DH149", "DH157", "DH170")

infection_plot <- quantiles_out |> 
  filter(name == "inf") |> 
  mutate(r0 = factor(r0, levels = c(2,4), labels = c("R0 = 2", "R0 = 4")), 
         mosquito = factor(mosquito, levels = c("wmel", "walbb", "wt"), labels = c("wMel", "wAlbB", "WT")),
         year = time + as.Date("2025-01-01"),
         strain = as.numeric(level), strain = factor(strain, levels = 1:80, labels = strain_names)) |>
  ggplot() + 
  geom_line(aes(x = year, y = median, col = strain, group = strain)) +
  geom_ribbon(aes(x = year, ymin = q_0.25, ymax = q_0.75, fill = strain, group = strain), alpha = 0.4) +
  geom_ribbon(aes(x = year, ymin = q_0.025, ymax = q_0.975, fill = strain, group = strain), alpha = 0.2) +
  facet_grid(mosquito ~ r0, scale = "free_y") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Year", y = "Infections", fill = "Strain", color = "Strain") +
  theme(legend.position = "bottom")

ggsave(here(figure_path, "strain-infection-plot.png"), width = 180, height = 240, unit = "mm", dpi = 300, plot = infection_plot)


infection_plot_zoom <- quantiles_out |> 
  filter(name == "inf") |> 
  mutate(r0 = factor(r0, levels = c(2,4), labels = c("R0 = 2", "R0 = 4")), 
         mosquito = factor(mosquito, levels = c("wmel", "walbb", "wt"), labels = c("wMel", "wAlbB", "WT")),
         year = time + as.Date("2025-01-01"),
         strain = as.numeric(level), strain = factor(strain, levels = 1:80, labels = strain_names)) |>
  ggplot() + 
  geom_line(aes(x = year, y = median, col = strain, group = strain)) +
  geom_ribbon(aes(x = year, ymin = q_0.25, ymax = q_0.75, fill = strain, group = strain), alpha = 0.4) +
  geom_ribbon(aes(x = year, ymin = q_0.025, ymax = q_0.975, fill = strain, group = strain), alpha = 0.2) +
  facet_grid(mosquito ~ r0, scale = "free_y") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(ylim = c(0,150000)) +
  labs(x = "Year", y = "Infections", fill = "Strain", color = "Strain") +
  theme(legend.position = "bottom")

ggsave(here(figure_path, "strain-infection-zoom-plot.png"), width = 180, height = 240, unit = "mm", dpi = 300, plot = infection_plot_zoom)


## Strain diversity plot

strain_diversity_plot <- strain_div_out |> 
  mutate(r0 = factor(r0, levels = c(2,4), labels = c("R0 = 2", "R0 = 4")), 
         mosquito = factor(mosquito, levels = c("wmel", "walbb", "wt"), labels = c("wMel", "wAlbB", "WT")),
         year = time + as.Date("2025-01-01")) |> 
  group_by(time, mosquito, r0, year) |> 
  summarise(median = median(simpson_index, na.rm = TRUE), q_0.025 = quantile(simpson_index, 0.025, na.rm = TRUE), q_0.975 = quantile(simpson_index, 0.975, na.rm = TRUE),
            q_0.25 = quantile(simpson_index, 0.25, na.rm = TRUE), q_0.75 = quantile(simpson_index, 0.75, na.rm = TRUE)) |> 
  ggplot() +
  geom_line(aes(x = year, y = median), col = "#a63603") +
  geom_ribbon(aes(x = year, ymin = q_0.25, ymax = q_0.75), alpha = 0.4, fill = "#a63603") +
  geom_ribbon(aes(x = year, ymin = q_0.025, ymax = q_0.975), alpha = 0.2, fill = "#a63603") +
  facet_grid(mosquito ~ r0, scale = "free_y") +
  labs(x = "Year", y = "Strain diversity") 
ggsave(here(figure_path, "strain-diversity-plot.png"), width = 180, height = 150, unit = "mm", dpi = 300, plot = strain_diversity_plot)


strain_number_plot <- strain_div_out |> 
  mutate(r0 = factor(r0, levels = c(2,4), labels = c("R0 = 2", "R0 = 4")), 
         mosquito = factor(mosquito, levels = c("wmel", "walbb", "wt"), labels = c("wMel", "wAlbB", "WT")),
         year = time + as.Date("2025-01-01")) |> 
  group_by(time, mosquito, r0, year) |> 
  summarise(median = median(abundance, na.rm = TRUE), q_0.025 = quantile(abundance, 0.025, na.rm = TRUE), q_0.975 = quantile(abundance, 0.975, na.rm = TRUE),
            q_0.25 = quantile(abundance, 0.25, na.rm = TRUE), q_0.75 = quantile(abundance, 0.75, na.rm = TRUE)) |> 
  ggplot() +
  geom_line(aes(x = year, y = median), col = "#980043") +
  geom_ribbon(aes(x = year, ymin = q_0.25, ymax = q_0.75), alpha = 0.4, fill = "#980043") +
  geom_ribbon(aes(x = year, ymin = q_0.025, ymax = q_0.975), alpha = 0.2, fill = "#a63603") +
  facet_grid(mosquito ~ r0, scale = "free_y") +
  labs(x = "Year", y = "Number of circulating strains") 
ggsave(here(figure_path, "strain-number-plot.png"), width = 180, height = 150, unit = "mm", dpi = 300, plot = strain_number_plot)


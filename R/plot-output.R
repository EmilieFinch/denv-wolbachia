# Plot simulation output

## Set up

date <- "2025-06-04"
figure_path <- here("figures", "simulation", date)
ifelse(!dir.exists(figure_path), dir.create(figure_path, recursive = TRUE), FALSE)

## Load in results 

quantiles_out <- NULL
strain_div_out <- NULL

r0s <- c(1.5,2,4)
mosquitos <- c("wmel", "walbb", "wt")
for(r0 in r0s){
  for(mos in mosquitos){
    quant <- qread(here("output", "simulation", date, paste0("quantile-out_r0-", r0, "_", mos, ".qs")))
    strain <- qread(here("output", "simulation", date, paste0("strain-div-out_r0-", r0, "_", mos, ".qs")))
    
    quantiles_out <- rbind(quantiles_out, quant)
    strain_div_out <- rbind(strain_div_out, strain)
  }
}
rm(quant, strain)

# Plotting

## Infection plot

strain_names <- c("DH002", "DH006", "DH009", "DH014", "DH019", "DH021", "DH022", "DH026", "DH035", "DH040", "DH040.1",
               "DH131", "DH141", "DH142", "DH143", "DH144", "DH145", "DH146", "DH147", "DH148", "DH060", "DH061", 
               "DH064", "DH069A", "DH071", "DH073", "DH074", "DH077", "DH081", "DH087", "DH091", "DH099", "DH101", 
               "DH134", "DH134.1", "DH154", "DH156", "DH158", "DH159", "DH161", "DH102", "DH104", "DH105", "DH110", 
               "DH162", "DH164", "DH165", "DH167", "DH168", "DH169", "DH102.1", "DH104.1", "DH105.1", "DH110.1", "DH162.1", 
               "DH164.1", "DH165.1", "DH167.1", "DH168.1", "DH169.1", "DH114", "DH116", "DH118", "DH119", "DH122", "DH139", 
               "DH140", "DH149", "DH157", "DH170", "DH114.1", "DH116.1", "DH118.1", "DH119.1", "DH122.1", "DH139.1", "DH140.1", 
               "DH149.1", "DH157.1", "DH170.1")

quantiles_plot <- quantiles_out |> 
  filter(name == "inf") |> 
  mutate(r0 = factor(r0, levels = c(1.5, 2,4), labels = c("R0 = 1.5", "R0 = 2", "R0 = 4")), 
         mosquito = factor(mosquito, levels = c("wmel", "walbb", "wt"), labels = c("wMel", "wAlbB", "WT")),
         year = time + as.Date("2025-01-01"),
         strain = as.numeric(level), strain = factor(strain, levels = 1:80, labels = strain_names))

infection_plot <-  quantiles_plot |>
  ggplot() + 
  geom_ribbon(aes(x = year, ymin = q_0.25, ymax = q_0.75, fill = strain, group = strain), alpha = 0.4) +
  geom_ribbon(aes(x = year, ymin = q_0.025, ymax = q_0.975, fill = strain, group = strain), alpha = 0.2) +
    geom_line(aes(x = year, y = median, col = strain, group = strain)) +
  facet_grid(mosquito ~ r0, scale = "free_y") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Year", y = "Infections", fill = "Strain", color = "Strain") +
  theme(legend.position = "bottom")

ggsave(here(figure_path, "strain-infection-plot.png"), width = 180, height = 240, unit = "mm", dpi = 300, plot = infection_plot)


infection_plot_zoom <- quantiles_plot |> 
  ggplot() + 
  geom_ribbon(aes(x = year, ymin = q_0.25, ymax = q_0.75, fill = strain, group = strain), alpha = 0.4) +
  geom_ribbon(aes(x = year, ymin = q_0.025, ymax = q_0.975, fill = strain, group = strain), alpha = 0.2) +
  geom_line(aes(x = year, y = median, col = strain, group = strain)) +
  facet_grid(mosquito ~ r0, scale = "free_y") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(ylim = c(0,1000)) +
  labs(x = "Year", y = "Infections", fill = "Strain", color = "Strain") +
  theme(legend.position = "bottom")

ggsave(here(figure_path, "strain-infection-zoom-plot.png"), width = 180, height = 240, unit = "mm", dpi = 300, plot = infection_plot_zoom)


## Strain diversity plot

strain_diversity_plot <- strain_div_out |> 
  mutate(r0 = factor(r0, levels = c(1.5, 2,4), labels = c("R0 = 1.5", "R0 = 2", "R0 = 4")), 
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
  mutate(r0 = factor(r0, levels = c(1.5, 2,4), labels = c("R0 = 1.5", "R0 = 2", "R0 = 4")), 
         mosquito = factor(mosquito, levels = c("wmel", "walbb", "wt"), labels = c("wMel", "wAlbB", "WT")),
         year = time + as.Date("2025-01-01")) |> 
  group_by(time, mosquito, r0, year) |> 
  summarise(median = median(number, na.rm = TRUE), q_0.025 = quantile(number, 0.025, na.rm = TRUE), q_0.975 = quantile(number, 0.975, na.rm = TRUE),
            q_0.25 = quantile(number, 0.25, na.rm = TRUE), q_0.75 = quantile(number, 0.75, na.rm = TRUE)) |> 
  ggplot() +
  geom_line(aes(x = year, y = median), col = "#980043") +
  geom_ribbon(aes(x = year, ymin = q_0.25, ymax = q_0.75), alpha = 0.4, fill = "#980043") +
  geom_ribbon(aes(x = year, ymin = q_0.025, ymax = q_0.975), alpha = 0.2, fill = "#a63603") +
  facet_grid(mosquito ~ r0, scale = "free_y") +
  labs(x = "Year", y = "Number of circulating strains") 
ggsave(here(figure_path, "strain-number-plot.png"), width = 180, height = 150, unit = "mm", dpi = 300, plot = strain_number_plot)


# Script to process experimental data
set.seed(84)
yale_results <- read_xlsx(here("data", "DENV1-4_Salje_format.xlsx")) |> 
  clean_names() |> 
  arrange(serotype)

dissemination_rates <- yale_results |> 
  group_by(mosquito_strain, virus) |> 
  mutate(group_id = cur_group_id()) |> 
  group_by(mosquito_strain, virus, serotype, genotype) |> 
  summarise(dissemination_raw = sum(leg_dissemination)/n()) |> 
  group_by(virus) |> 
  mutate(dissemination_wt = dissemination_raw[mosquito_strain == "WT"]) |> 
  ungroup() |> 
  mutate(dissemination_norm = dissemination_raw/dissemination_wt) |> 
  filter(mosquito_strain != "WT")

# Stratified nonparametric boostrapping (with replacement)

percent_diff <- function(data, indices) {
  d <- data[indices,]
  data_out <- d |> 
    group_by(mosquito_strain, virus, serotype, genotype) |> 
    summarise(dissemination_raw = sum(leg_dissemination)/n(), .groups = "keep") |> 
    group_by(virus) |> 
    mutate(dissemination_wt = dissemination_raw[mosquito_strain == "WT"]) |> 
    ungroup() |> 
    mutate(dissemination_norm = dissemination_raw/dissemination_wt,
           dissemination_norm = case_when(dissemination_wt == 0 ~ 1, 
                                          dissemination_norm > 1 ~ 1,
                                          T ~ dissemination_norm)) |> 
    filter(mosquito_strain != "WT") |> 
    pull(dissemination_norm)
  
  return(data_out)
}

mosquito_set <- c("wMel", "wAlbB")
strain_set <- unique(yale_results$virus)

# wMel estimates

wmel_mat <- matrix(nrow = 100, ncol = 0)
for(strain in strain_set){
subset <- yale_results |> 
  filter(virus == strain & (mosquito_strain == "wMel" | mosquito_strain == "WT")) |> 
  mutate(mosquito_strain = as.factor(mosquito_strain))
boot <- boot(data = subset, statistic = percent_diff, R = 100, strata = subset$mosquito_strain)
names <- c(colnames(wmel_mat), strain)
wmel_mat <- cbind(wmel_mat, boot$t)
colnames(wmel_mat) <- names
}
qsave(wmel_mat, here("data", "wmel-inhib-values.qs"))

# wAlbB estimates

walbb_mat <- matrix(nrow = 100, ncol = 0)
for(strain in strain_set){
  subset <- yale_results |> 
    filter(virus == strain & (mosquito_strain == "wAlbB" | mosquito_strain == "WT")) |> 
    mutate(mosquito_strain = as.factor(mosquito_strain))
  boot <- boot(data = subset, statistic = percent_diff, R = 100, strata = subset$mosquito_strain)
  names <- c(colnames(walbb_mat), strain)
  walbb_mat <- cbind(walbb_mat, boot$t)
  colnames(walbb_mat) <- names
}

qsave(walbb_mat, here("data", "walbb-inhib-values.qs"))
# Plot values

inhib_values <- as.data.frame(wmel_mat) |> 
  pivot_longer(everything()) |> 
  mutate(mosquito = "WMel") |> 
  rbind(as.data.frame(walbb_mat) |> 
              pivot_longer(everything()) |> 
              mutate(mosquito = "wAlbB")) |> 
  mutate(value = as.numeric(value)) 

  inhib_plot <- ggplot(inhib_values) +
  geom_histogram(aes(x = value, group = mosquito, fill = mosquito), alpha = 0.5, position = "identity") +
  facet_wrap(~name) +
  scale_fill_manual(values = c("#762a83", "#1b7837"))+
  scale_x_continuous(breaks = c(0,0.5,1)) +
    labs(y = "Frequency", x = "Inhibition", fill = NULL) +
    theme(legend.key.height = unit(0.3, "cm"), legend.position = "bottom")

  ggsave(here("figures", "wolbachia-inhib-values.png"),inhib_plot, width = 180, height = 150, units = "mm", dpi = 300)

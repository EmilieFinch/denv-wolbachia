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
           dissemination_norm = case_when(dissemination_wt == 0 ~ NA_integer_, 
                                          dissemination_norm > 1 ~ 1,
                                          T ~ dissemination_norm)) |> 
    filter(mosquito_strain != "WT") |> 
    pull(dissemination_norm)
  
  return(data_out)
}

strain_set1 <- unique(yale_results$virus[yale_results$serotype == 1])
strain_set1[11] <- "DH040" # replace strains with no dissemination in WT
strain_set2 <- unique(yale_results$virus[yale_results$serotype == 2])
strain_set2[15] <- "DH134" # replace strain with no dissemination in WT
# double serotypes 3 and 4 so total number of strains per serotype is consistent
strain_set3 <- rep(unique(yale_results$virus[yale_results$serotype == 3]),2) 
strain_set4 <- rep(unique(yale_results$virus[yale_results$serotype == 4]),2)

strain_set <- c(strain_set1, strain_set2, strain_set3, strain_set4)
rm(strain_set1, strain_set2, strain_set3, strain_set4)
mosquito_set <- c("wMel", "wAlbB")


# wMel estimates

wmel_mat <- matrix(nrow = 100, ncol = 0)
for(strain in strain_set){
subset <- yale_results |> 
  filter(virus == strain & (mosquito_strain == "wMel" | mosquito_strain == "WT")) |> 
  mutate(mosquito_strain = as.factor(mosquito_strain))
boot <- boot(data = subset, statistic = percent_diff, R = 1000, strata = subset$mosquito_strain)
valid_samples <- boot$t[!is.na(boot$t)]
selected_samples <- sample(valid_samples, 100, replace = FALSE)
names <- c(colnames(wmel_mat), strain)
wmel_mat <- cbind(wmel_mat, selected_samples)
colnames(wmel_mat) <- names
}
qsave(wmel_mat, here("data", "wmel-inhib-values.qs"))

# wAlbB estimates

walbb_mat <- matrix(nrow = 100, ncol = 0)
for(strain in strain_set){
  subset <- yale_results |> 
    filter(virus == strain & (mosquito_strain == "wAlbB" | mosquito_strain == "WT")) |> 
    mutate(mosquito_strain = as.factor(mosquito_strain))
  boot <- boot(data = subset, statistic = percent_diff, R = 1000, strata = subset$mosquito_strain)
  valid_samples <- boot$t[!is.na(boot$t)]
  selected_samples <- sample(valid_samples, 100, replace = FALSE)
  names <- c(colnames(walbb_mat), strain)
  walbb_mat <- cbind(walbb_mat, selected_samples)
  colnames(walbb_mat) <- names
}
qsave(walbb_mat, here("data", "walbb-inhib-values.qs"))

# Plot values
wmel_mat <- as.data.frame(wmel_mat)
names(wmel_mat) <- make.unique(names(wmel_mat))
walbb_mat <- as.data.frame(walbb_mat)
names(walbb_mat) <- make.unique(names(walbb_mat))

inhib_values <- wmel_mat |> 
  pivot_longer(everything()) |>
  mutate(mosquito = "WMel") |> 
  rbind(walbb_mat |> 
              pivot_longer(everything()) |> 
              mutate(mosquito = "wAlbB")) |> 
  mutate(value = as.numeric(value)) 

  inhib_plot <- ggplot(inhib_values) +
  geom_histogram(aes(x = value, group = mosquito, fill = mosquito), alpha = 0.5, position = "identity") +
  facet_wrap(~name) +
  scale_fill_manual(values = c("#762a83", "#1b7837"))+
  scale_x_continuous(breaks = c(0,0.5,1)) +
  scale_y_continuous(breaks = c(0,50,100)) +
    labs(y = "Frequency", x = "Inhibition", fill = NULL) +
    theme(legend.key.height = unit(0.3, "cm"), legend.position = "bottom")

  ggsave(here("figures", "wolbachia-inhib-values.png"),inhib_plot, width = 180, height = 150, units = "mm", dpi = 300)

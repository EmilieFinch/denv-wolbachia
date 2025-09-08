# Script to process experimental data
set.seed(84)
yale_results <- read_xlsx(here("data", "DENV1-4_Salje_format_08.26.2025.xlsx")) |> 
  clean_names() |> 
  arrange(serotype)

dissemination_rates_new <- yale_results |> 
  group_by(mosquito_strain, virus) |> 
  mutate(group_id = cur_group_id()) |> 
  group_by(mosquito_strain, virus, serotype, genotype) |> 
  summarise(dissemination_raw = sum(leg_dissemination)/n()) |> 
  group_by(virus) |> 
  mutate(dissemination_wt = dissemination_raw[mosquito_strain == "WT"]) |> 
  ungroup() |> 
  mutate(dissemination_norm = dissemination_raw/dissemination_wt) |> 
  filter(mosquito_strain != "WT")

wmel_out <- dissemination_rates |> 
  filter(!is.na(dissemination_norm)) |> 
  filter(mosquito_strain == "wMel") |> 
  select(virus, serotype, dissemination_norm) |> 
  arrange(serotype) 
qsave(wmel_out, here("data", "wmel-inhib.qs"))

walbb_out <- dissemination_rates |> 
  filter(!is.na(dissemination_norm)) |> 
  filter(mosquito_strain == "wAlbB") |> 
  select(virus, serotype, dissemination_norm) |> 
  arrange(serotype) 
qsave(walbb_out, here("data", "walbb-inhib.qs"))


  dissemination_out <- dissemination_rates |> 
    select(dissemination_norm) |> 
    mutate(dissemination_norm = round(dissemination_norm, 3)) |> 
    rename(dissemination = dissemination_norm) |> 
    filter(!is.na(dissemination)) |> 
    unique() |> 
    arrange(dissemination)
    
write.csv(dissemination_out, here("data", "strain-specific-input-values.csv"))
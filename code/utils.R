# Utils fns

# Function to load demography

load_demography <- function(){
wpp_brazil <- read.csv(here("data", "wpp-population-single-age_brazil.csv")) |> 
  clean_names() |> 
  select(c(year, starts_with("x"))) |> 
  rename_with(~ str_replace(., "^x", ""), starts_with("x")) |> 
  pivot_longer(cols = -year, names_to = "age", values_to = "population") |>
  mutate(age = as.numeric(age),
         age_group = case_when(age < 80 ~ as.character(age),
                               T ~ "80+"),
         age_group = factor(age_group, levels = c(0:79, "80+"))) |> 
  group_by(year, age_group) |> 
  summarise(population = sum(population*1000), .groups = "keep") |>  # wpp data is presented in thousands
  ungroup()

wpp_brazil_proj <- read.csv(here("data", "wpp-population-projections-single-age_brazil.csv")) |> 
  clean_names() |> 
  select(c(year, starts_with("x"))) |> 
  rename_with(~ str_replace(., "^x", ""), starts_with("x")) |> 
  mutate(across(-year, ~as.numeric(gsub(" ", "", .x)))) |> 
  pivot_longer(cols = -year, names_to = "age", values_to = "population") |>
  mutate(age = as.numeric(age),
         age_group = case_when(age < 80 ~ as.character(age),
                               T ~ "80+"),
         age_group = factor(age_group, levels = c(0:79, "80+"))) |> 
  group_by(year, age_group) |> 
  summarise(population = sum(population*1000), .groups = "keep") |>  # wpp data is presented in thousands
  ungroup()
wpp_all <- rbind(wpp_brazil, wpp_brazil_proj)

return(wpp_all)
}

wrangle_demography <- function(demog, year_start, year_end, pad_left = 0){
  demog <- brazil_demog |> 
    filter(year >= year_start & year <= year_end) |> 
    pivot_wider(id_cols = "age_group", names_from = "year", values_from = "population") |> 
    select(- age_group) |> 
    round() |> 
  as.matrix()
  
  if(pad_left != 0){ 
  demog <- cbind(matrix(demog[,1], nrow = nrow(demog), ncol = pad_left), demog)}
  
  return(demog)
}

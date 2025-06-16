if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here, dplyr, janitor, tidyr, readxl, ggplot2,
  showtext, stringr, ggplot2, patchwork, ggpubr,
  tibble, yaml,qs, lubridate, abind, data.table, 
  spdep, igraph, boot, scales, purrr, readr, cowplot
)

if(!require("odin2")) install.packages(
  "odin2",
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))

if(!require("dust2")) install.packages(
  "dust2",
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))

# Set up plot font
plot_font <- "Open Sans"
font_add_google(plot_font)
showtext_opts(dpi = 300)
showtext_auto()

# Set theme for plotting
theme_set(theme_classic() +
            theme(plot.title = element_text(size = 8, family = plot_font),
                  axis.title = element_text(size = 8, family = plot_font),
                  axis.text = element_text(size = 7, family = plot_font),
                  legend.title = element_text(size = 8, family = plot_font), 
                  legend.text = element_text(size = 7, family = plot_font),
                  legend.key.height = unit(0.2, "cm"), legend.position = "bottom",
                  strip.background = element_rect(fill = "#082544"), 
                  strip.text = element_text(color = "white", size = 6, family = plot_font),
                  axis.line = element_line(color = "#082544"),
                  axis.ticks = element_line(color = "#082544")))


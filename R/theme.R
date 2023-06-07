library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(patchwork)
library(paletteer)
library(ggsci)
library(scico)
library(colorspace)
library(drlib)
library(pantone35)
library(cowplot)
library(nationalparkcolors)
library(fishualize)
library(BuenColors)
library(ggstance)

# Define theme
clean_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
        axis.text.y = element_text(size = size_sm),
        axis.title.x = element_text(size = size_lg),
        axis.title.y = element_text(size = size_lg),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        strip.text = element_text(size = size_sm),
        # strip.background = element_rect(fill = "grey90", color = "grey90",
        #                                 size = 0),
        strip.background = element_blank(),
        axis.line.y = element_line(colour = "grey50"),
        axis.line.x = element_line(colour = "grey50"), 
        axis.ticks = element_line(colour = "grey50"),
        legend.position = "top",
        legend.text = element_text(size = size_sm),
        legend.title = element_text(size = size_sm),
        legend.key.size = unit(0.6, "lines"),
        legend.margin = margin(rep(0, 4)),
        # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
        # legend.box.spacing = ggplot2::margin(rep(0, 4)),
        legend.background = element_blank(),
        plot.title = element_text(size = size_lg, hjust = 0.5))
}

grid_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_text(size = size_lg),
          strip.text = element_text(size = size_sm),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_line(colour = "grey50"),
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
          # legend.box.spacing = ggplot2::margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5))
}

boxed_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_text(size = size_lg),
          panel.grid = element_blank(),
          strip.text = element_text(size = size_sm),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_line(colour = "grey50"),
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
          # legend.box.spacing = ggplot2::margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5))
}

darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # remove zero
  l <- gsub("0e\\+00", "0", l)
  # remove one
  l <- gsub("^1e\\+00", "1", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + from exponent
  l <- gsub("e\\+" ,"e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove 1 x 10^ (replace with 10^)
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

plot_pal = function(pal) {
  grid::grid.raster(pal, interpolate=F)
}

cubehelix = function(n_colors) {
  colours = c("#000000", "#1A1935", "#15474E", "#2B6F39", "#767B33", "#C17A6F",
              "#D490C6", "#C3C0F2")
  idxs = 0.3
  if (n_colors > 1)
    idxs = seq(0, 1, 1 / (n_colors - 1))
  colour_ramp(colours)(idxs)
}

kinney6 = c(
  "#c6c3bf",
  "#119e87",
  "#53bad3",
  "#559ed2",
  "#3b5687",
  "#e34d3b")

Gpal = c("#E30F17", "#0296E1", "#F49203", "#E5E5E5", "#E0A4D1") %>%
  # extended version
  ## https://www.sciencedirect.com/science/article/pii/S0896627316000106
  c('#2C8942', '#E8B820')

colours.cafe447 = c('#ffb838', '#fee5a5', '#f7f6fee', '#486d87')
colours.cafe433 = c('#077893', '#e3deca', '#fcfaf1', '#ff9465')
colours.cafe425 = c('#2B5B6C', '#C7CFAC', '#FCFAF1', '#E34F33', '#FFC87E')
colours.cafe322 = c("#7bbaea", "#d46363", "#fbdaa7", "#fcfaf2", "#30598c")

winsorize = function(vec, limits) {
  lower_limit = limits[1]
  upper_limit = limits[2]
  if (!is.na(upper_limit))
    vec[vec > upper_limit] = upper_limit
  if (!is.na(lower_limit))
    vec[vec < lower_limit] = lower_limit
  return(vec)
}

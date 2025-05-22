
#https://r-graph-gallery.com/330-bubble-map-with-ggplot2.html
# Libraries
library(ggplot2)
library(dplyr)
# Get the world polygon and extract UK
#get giscoR package
if(!require(giscoR)){
  #https://ropengov.github.io/giscoR/
  #remotes::install_github("rOpenGov/giscoR")
  install.packages("giscoR",repos = c("https://ropengov.r-universe.dev", 
                                      "https://cloud.r-project.org")
  )
}  
library(giscoR)
UK <- gisco_get_countries(country = "UK", resolution = 1)
# Get a data frame with longitude, latitude, and size of bubbles (a bubble = a city)
library(maps)
data <- world.cities %>% filter(country.etc == "UK")
# Left chart
ggplot() +
  geom_sf(data = UK, fill = "grey", alpha = 0.3) +
  geom_point(data = data, aes(x = long, y = lat)) +
  theme_void() +
  ylim(50, 59)

# Second graphic with names of the 10 biggest cities
library(ggrepel)
ggplot() +
  geom_sf(data = UK, fill = "grey", alpha = 0.3) +
  geom_point(data = data, aes(x = long, y = lat, alpha = pop)) +
  geom_text_repel(
    data = data %>% arrange(pop) %>% tail(10),
    aes(x = long, y = lat, label = name), size = 5
  ) +
  geom_point(
    data = data %>% arrange(pop) %>% tail(10), aes(x = long, y = lat),
    color = "red", size = 3
  ) +
  theme_void() +
  ylim(50, 59) +
  theme(legend.position = "none")
# Left: use size and color
ggplot() +
  geom_sf(data = UK, fill = "grey", alpha = 0.3) +
  geom_point(data = data, aes(x = long, y = lat, size = pop, color = pop)) +
  scale_size_continuous(range = c(1, 12)) +
  scale_color_viridis_c(trans = "log") +
  theme_void() +
  ylim(50, 59)

# Center: reorder your dataset first! Big cities appear later = on top
data %>%
  arrange(pop) %>%
  mutate(name = factor(name, unique(name))) %>%
  ggplot() +
  geom_sf(data = UK, fill = "grey", alpha = 0.3) +
  geom_point(aes(x = long, y = lat, size = pop, color = pop), alpha = 0.9) +
  scale_size_continuous(range = c(1, 12)) +
  scale_color_viridis_c(trans = "log") +
  theme_void() +
  ylim(50, 59) +
  theme(legend.position = "none")

# Right: just use arrange(desc(pop)) instead
data %>%
  arrange(desc(pop)) %>%
  mutate(name = factor(name, unique(name))) %>%
  ggplot() +
  geom_sf(data = UK, fill = "grey", alpha = 0.3) +
  geom_point(aes(x = long, y = lat, size = pop, color = pop), alpha = 0.9) +
  scale_size_continuous(range = c(1, 12)) +
  scale_color_viridis_c(trans = "log") +
  theme_void() +
  ylim(50, 59) +
  theme(legend.position = "none")
# Create breaks for the color scale
mybreaks <- c(0.02, 0.04, 0.08, 1, 7)

# Reorder data to show biggest cities on top
data <- data %>%
  arrange(pop) %>%
  mutate(name = factor(name, unique(name))) %>%
  mutate(pop = pop / 1000000)

# Build the map
data %>%
  ggplot() +
  geom_sf(data = UK, fill = "grey", alpha = 0.3) +
  geom_point(aes(x = long, y = lat, size = pop, 
                 color = pop, alpha = pop),
             shape = 20, stroke = FALSE
  ) +
  scale_size_continuous(
    name = "Population (in M)", trans = "log",
    range = c(1, 12), breaks = mybreaks
  ) +
  scale_alpha_continuous(
    name = "Population (in M)", trans = "log",
    range = c(0.1, .9), breaks = mybreaks
  ) +
  scale_color_viridis_c(
    option = "magma", trans = "log",
    breaks = mybreaks, name = "Population (in M)"
  ) +
  theme_void() +
  guides(colour = guide_legend()) +
  ggtitle("The 1000 biggest cities in the UK") +
  theme(
    legend.position = c(1, 0.6),
    text = element_text(color = "#22211d"),
    plot.margin = margin(r = 2, l = 2, unit = "cm"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#4e4d47"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )


#
# group points

#https://stackoverflow.com/questions/57843139/r-visualization-sensible-repel-points-on-map-beeswarm
library(sf)
#
places = sf::st_read("https://gist.githubusercontent.com/peeter-t2/9646a4169e993948fa97f6f503a0688b/raw/cb4e910bf153e51e3727dc9d1c73dd9ef86d2556/kih1897m.geojson", stringsAsFactors = FALSE)

schools <- read_tsv("https://gist.github.com/peeter-t2/34467636b3c1017e89f33284d7907b42/raw/6ea7dd6c005ef8577b36f5e84338afcb6c76b707/school_nums.tsv")
schools_geo <- merge(places,schools,by.x="KIHELKOND",by.y="Kihelkond") #94 matches

# p<- schools_geo %>% 
#   ggplot()+
#   geom_sf(data=schools_geo)+
#   geom_sf(data=st_centroid(schools_geo),aes(size=value))+
#   theme_bw()
# p

st_crs(schools_geo) <- 3301
# we can set type = 'hexagonal', 'regular' or 'random'
school_pts <- schools_geo %>% st_sample(size = .$value, type = 'hexagonal')


schools_geo %>% 
  ggplot()+
  geom_sf()+
  geom_sf(data=school_pts, size = .8)+
  theme_bw()


# put values on scale between 0 and 1
scale_fact <- (max(schools_geo$value) -  schools_geo$value) / (max(schools_geo$value) - min(schools_geo$value)) 
# re-scale between 0.2 and 0.9
scale_fact <- scale_fact * (0.9 - 0.2) + 0.2
# reverse the scale 
scale_fact <-  max(scale_fact) + min(scale_fact) - scale_fact 

# apply the scale factor
schools_centroid <- st_geometry(st_centroid(schools_geo))
schools_geo_rescaled <- (st_geometry(schools_geo) - schools_centroid) * scale_fact + schools_centroid

school_pts <- schools_geo_rescaled %>% 
  st_sf(crs = 3301) %>% 
  bind_cols(value = schools_geo$value) %>%
  st_sample(size = .$value, type = 'hexagonal')


# plot
schools_geo %>% 
  ggplot()+
  geom_sf()+
  geom_sf(data=school_pts, size = .8)+
  theme_bw()

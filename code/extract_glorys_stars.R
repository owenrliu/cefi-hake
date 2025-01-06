# trying to make a template raster with stars, for later plotting
# going off of https://r-spatial.github.io/stars/articles/stars1.html#cropping-a-rasters-extent
library(stars)
library(tidyverse)
library(here)
library(rnaturalearth)

# coastline, joined US and Canada
coast <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada','Idaho','Montana'))
coastcn <- ne_countries(country="Canada",scale=50,returnclass='sf')
coast <- st_union(coast,coastcn)

# Footprint of the FEAT survey data, as a polygon (created in 00_make_FEAT_grid.R)
feat_foot <- read_rds(here('data','grids','FEAT_footprint.rds'))
# in lat/lon coords
feat_foot_ll <- feat_foot %>% st_transform(4326)

# TESTING

glor_rast <- read_stars(here('data','glorys','phys','raw','glorys_physics_1995.nc'),sub="thetao") %>% 
  st_set_crs(4326)
slc <- slice(glor_rast,index=1,along='time') %>% slice("depth",1)
cropped <- slc[feat_foot_ll]
test <- cropped %>% st_as_sf()%>% 
  mutate(thetao=units::drop_units(thetao))
ggplot(test,aes(fill=thetao))+geom_sf(col=NA)+scale_fill_viridis(option='turbo')


sstmon <- glor_rast %>%
  slice("depth",1) %>% 
  as.tbl_cube.stars() %>% 
  group_by(x,y,time) %>% 
  summarise(mean_thetao=mean(thetao))


ggplot()+
  geom_stars(data=cropped)+
  scale_fill_viridis(option="turbo")

test_extract <- st_extract(slc,feat_xy %>% filter(year==1995)) %>% 
  mutate(thetao=units::drop_units(thetao))
test_extract2 <- st_extract(slc,feat_xy %>% filter(year==1995),bilinear = T)
test_extract %>% mutate(thetao_bl=test_extract2$thetao) %>% mutate(diff=thetao-thetao_bl) %>% ggplot(aes(diff))+geom_density()

# not much of a difference between bilinear interpolation and just using nearest neighbors
ggplot()+geom_sf(data=test_extract,aes(color=thetao))+geom_sf()+scale_color_viridis(option='turbo')
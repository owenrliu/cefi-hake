# Build matching keys to join GLORYS data to FEAT survey observations and projection grids
# There are two types of data to join here:
# physics and biogeochemistry, which are on different resolution grids
library(tidyverse)
library(here)
library(sf)
library(nngeo)
library(rnaturalearth)

# background map
coast <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada','Idaho','Montana'))
coastcn <- ne_countries(country="Canada",scale=50,returnclass='sf')
coast <- st_union(coast,coastcn)

# observational data and projection grid
# FEAT data with survey locations
# data were compiled from the original excel files in "explore FEAT.qmd"
feat <- read_rds(here('data','CONFIDENTIAL','un-kriged_aged_output_allyears.rds'))
glimpse(feat)
# template 5km projection grid (from 00_make_FEAT_grid.R)
feat_gr <- read_rds(here('data','grids','FEAT_5km_grid.rds'))

# get the xy locations of both the data and the grid
feat_xy <- feat %>% 
  distinct(year,Lat,Lon) %>% 
  st_as_sf(coords=c("Lon","Lat"),crs=4326)
feat_gr_xy <- feat_gr %>% 
  st_set_geometry(NULL) %>% 
  st_as_sf(coords=c("E_km","N_km"),crs=st_crs(feat_gr)) %>% 
  st_transform(4326)

# now, get the xy locations of the GLORYS grids from one of the filtered files from 02_extract_GLORYS
glor_phys_xy <- list.files(here('data','glorys','phys','filtered'),full.names = T)[1] %>% 
  read_rds() %>% 
  distinct(longitude,latitude) %>% 
  st_as_sf(coords=c("longitude","latitude"),crs=4326)
glor_bgc_xy <- list.files(here('data','glorys','bgc','filtered'),full.names = T)[1] %>% 
  read_rds() %>% 
  distinct(longitude,latitude) %>% 
  st_as_sf(coords=c("longitude","latitude"),crs=4326)

# find 3 nearest neighbors and save distances
# we probably need to do this 4 times, for the 4 combinations
# observations vs. grid, and GLORYS physics vs. bgc
feat_obs_glor_phys_join <- st_nn(feat_xy,glor_phys_xy,k=3,returnDist=T)

feat_gr_glor_phys_join <- st_nn(feat_gr_xy,glor_phys_xy,k=3,returnDist=T)

feat_obs_glor_bgc_join <- st_nn(feat_xy,glor_bgc_xy,k=3,returnDist=T)

feat_gr_glor_bgc_join <- st_nn(feat_gr_xy,glor_bgc_xy,k=3,returnDist=T)

# make the 3 nearest neighbor distances into a vector of weights (for inverse squared (power=2) distance weighting)
feat_obs_glor_phys_join[["invdist"]]<-feat_obs_glor_phys_join[["dist"]] %>% map(\(x) 1/x^2)
feat_obs_glor_phys_join[["weights"]]<-feat_obs_glor_phys_join[["invdist"]] %>% map(\(x)x/sum(x))

feat_gr_glor_phys_join[["invdist"]]<-feat_gr_glor_phys_join[["dist"]] %>% map(\(x) 1/x^2)
feat_gr_glor_phys_join[["weights"]]<-feat_gr_glor_phys_join[["invdist"]] %>% map(\(x)x/sum(x))

feat_obs_glor_bgc_join[["invdist"]]<-feat_obs_glor_bgc_join[["dist"]] %>% map(\(x) 1/x^2)
feat_obs_glor_bgc_join[["weights"]]<-feat_obs_glor_bgc_join[["invdist"]] %>% map(\(x)x/sum(x))

feat_gr_glor_bgc_join[["invdist"]]<-feat_gr_glor_bgc_join[["dist"]] %>% map(\(x) 1/x^2)
feat_gr_glor_bgc_join[["weights"]]<-feat_gr_glor_bgc_join[["invdist"]] %>% map(\(x)x/sum(x))

# Now we should have everything we need to match- 3 nearest-neighbor grid IDs, distances, and weights
# Let's collect this info into matching keys we can use later
feat_obs_glor_match <- feat %>% 
  # this is the same as the way we made 'feat_xy' above
  distinct(year,Lat,Lon) %>% 
  # matching info for GLORYS physics
  mutate(glorID_phys=feat_obs_glor_phys_join$nn,
         glor_phys_dists=feat_obs_glor_phys_join$dist,
         glor_phys_weights=feat_obs_glor_phys_join$weights) %>% 
  # matching info for GLORYS bgc
  mutate(glorID_bgc=feat_obs_glor_bgc_join$nn,
         glor_bgc_dists=feat_obs_glor_bgc_join$dist,
         glor_bgc_weights=feat_obs_glor_bgc_join$weights)

feat_gr_glor_match <- feat_gr %>% 
  # this is the same as the way we made 'feat_gr_xy' above
  st_set_geometry(NULL) %>% 
  # matching info for GLORYS physics
  mutate(glorID_bgc=feat_gr_glor_bgc_join$nn,
         glor_bgc_dists=feat_gr_glor_bgc_join$dist,
         glor_bgc_weights=feat_gr_glor_bgc_join$weights) %>% 
  # matching info for GLORYS bgc
  mutate(glorID_phys=feat_gr_glor_phys_join$nn,
         glor_phys_dists=feat_gr_glor_phys_join$dist,
         glor_phys_weights=feat_gr_glor_phys_join$weights)

glimpse(feat_obs_glor_match)
glimpse(feat_gr_glor_match)

write_rds(feat_obs_glor_match,here('data','feat_obs_glorys_matchkey.rds'))
write_rds(feat_gr_glor_match,here('data','feat_gr_glorys_matchkey.rds'))

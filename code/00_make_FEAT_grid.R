# Make a regular prediction grid for FEAT data
# Using the footprint of the FEAT survey data, make a regular grid
# on which to project our eventual distribution models

library(tidyverse)
library(sf)
library(here)
library(fmesher) # for nonconvex hull calculation
library(rnaturalearth)
library(units)
library(marmap)
library(viridis)

# background map
# load west coast land for mapping
coast <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada','Idaho','Montana')) %>%
  st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=km")
coastcn <- ne_countries(country="Canada",scale=50,returnclass='sf')%>%
  st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=km")
# join US and Canada
coast <- st_union(coast,coastcn)
ggplot(coast)+geom_sf()

# FEAT data with survey locations
# data were compiled from the original excel files in "explore FEAT.qmd"
feat <- read_rds(here('data','CONFIDENTIAL','un-kriged_aged_output_allyears.rds'))
glimpse(feat)

# Make into spatial object, and project to a flat projection (UTM Zone 10, with km as linear unit)
feat_sf <- feat %>% 
  distinct(Lon,Lat) %>% 
  st_as_sf(coords=c("Lon","Lat"),crs=4326) %>%
  st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=km")
bbox <- st_bbox(feat_sf)

# plot
ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=feat_sf,size=0.25)+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])

# make a buffered outline of the data points using fmesher package
dat_hull <- fm_nonconvex_hull(feat_sf,concave = -0.004,convex = -0.008)

# crop this buffer (i.e., crop at the coastline) so it doesn't overlap land areas
dat_buffer <- st_difference(dat_hull,summarise(coast))

# plot it on top of what we just mapped
ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=dat_buffer,fill=NA,color='red')+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])

# Now we can make a grid using this buffer
# this is the grid that covers the entire rectangle of the bounding box of the data
feat_grid_full <- st_make_grid(dat_buffer,
                               # square vs. hexagonal grid cells
                               square=TRUE,
                               # cellsize, presumably in units of the same crs as the data
                               cellsize = 5 ) # 5km
# we probably do not really need that full rectangle, so instead, use the same buffered outline
# from above to crop the prediction grid

dat_buffer_sf <- dat_buffer %>% st_as_sf()
# quicksave the polygon footprint of the FEAT survey; it will be useful later
write_rds(dat_buffer_sf,here('data','grids','FEAT_footprint.rds'))

# Crop our newly created grid using a spatial filter
# Only cells that overlap with the buffered FEAT footprint will be kept
feat_grid_cropped <- feat_grid_full %>% 
  st_as_sf() %>% 
  # filter such that we only retain grid cells that touch the buffer
  st_filter(dat_buffer_sf)

# plot again
ggplot()+
  # new cropped grid
  geom_sf(data=feat_grid_cropped,color='gray70')+
  # some sample data points
  geom_sf(data=feat_sf %>% sample_n(10000))+
  # coastline
  geom_sf(data=coast,fill='gray80')+
  # FEAT footprint
  geom_sf(data=dat_buffer,fill=NA,color='red')+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])

# this 5km cropped grid looks good
# include the coordinates of the grid centers; this might be helpful for other matching we need to do
xy <- st_centroid(feat_grid_cropped) %>% st_coordinates()
feat_grid_cropped <- feat_grid_cropped %>% 
  mutate(E_km=xy[,1],N_km=xy[,2])

# Could add bathymetry from marmap
### Go get bathymetric data from NOAA to overlay on the samples

# bounding box for downloading bathymetry
# transformed back to lat/lon 
limits.for.map <- st_bbox(feat_grid_cropped %>% st_transform(4326))

# grid cell centers, in lat/lon
xy_ll <- feat_grid_cropped %>% 
  st_set_geometry(NULL) %>% 
  st_as_sf(coords=c("E_km","N_km"),crs=st_crs(feat_grid_cropped)) %>% 
  st_transform(4326) %>% 
  st_coordinates()

# download bathymetry
b = getNOAA.bathy(lon1 = limits.for.map["xmin"],
                  lon2 = limits.for.map[ "xmax" ],
                  lat1 = limits.for.map["ymin"],
                  lat2 = limits.for.map["ymax"],
                  resolution = 1,keep = TRUE)

# We can pull depth from each grid point
xy.bathy <- get.depth(b,xy_ll,locator=F) %>% 
  rename(bathy=depth)

# add to the grid dataframe
feat_grid_final <- feat_grid_cropped %>% 
  mutate(bathy_m=xy.bathy$bathy)

# did this work?
feat_grid_bathy <- ggplot()+
  geom_sf(data=feat_grid_final,aes(fill=bathy_m),color=NA)+
  geom_sf(data=coast,fill='gray80')+
  scale_fill_viridis(option="A")+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])

feat_grid_bathy

# Distance to 200m isobath
# make the 200m contour line using the same bathymetry
bdf <- fortify(b) %>% 
  contoureR::getContourLines(levels=c(-200)) 

# make into a spatial object in the same CRS as the grid
bdf <- bdf %>% 
  st_as_sf(coords=c("x","y"),crs=4326) %>% 
  st_transform(st_crs(feat_grid_final)) %>% 
  mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2])

# save this for later!
write_rds(bdf,here('data','grids','200m_isobath.rds'))

# then, find nearest neighbors (i.e., minimum distance to the 200m isobath) for the prediction grid
library(nngeo)
nn_samp <- st_nn(st_centroid(feat_grid_cropped),bdf,k=1,returnDist = T)
nn_samp_dists <- nn_samp %>% pluck("dist") %>% unlist()

feat_grid_final <- feat_grid_final %>% 
  # distance calculated by st_nn() is in meters, so convert to km here
  mutate(dist_shelf_km=nn_samp_dists/1000)

# plot
feat_grid_shelf <- ggplot()+
  geom_sf(data=feat_grid_final,aes(fill=dist_shelf_km),color=NA)+
  geom_sf(data=coast,fill='gray80')+
  scale_fill_viridis(option="A",direction=-1)+
  xlim(bbox[1],bbox[3])+ylim(bbox[2],bbox[4])

feat_grid_shelf

# save the final grid
write_rds(feat_grid_final,here('data','grids','FEAT_5km_grid.rds'))

# Helper functions for the other scripts in this repo
library(tidyverse)
library(sf)
library(marmap)

# Take an sf object
# 1. find the bounding box and download bathymetry from marmap
# 2. extract depth using get.depth() function, and return

get_bathy_pts <- function(sf_dat){
  require(tidyverse)
  require(sf)
  require(marmap)
  # bounding box for downloading bathymetry
  # transformed back to lat/lon, with a small buffer
  limits.for.map <- sf_dat %>% st_transform(4326) %>% st_bbox()+c(-0.5,-0.5,0.5,0.5)
  
  # data points in lat/lon
  xy_ll <- sf_dat %>% 
    st_transform(4326) %>% 
    st_coordinates()
  
  # download bathymetry
  b = getNOAA.bathy(lon1 = limits.for.map["xmin"],
                    lon2 = limits.for.map[ "xmax" ],
                    lat1 = limits.for.map["ymin"],
                    lat2 = limits.for.map["ymax"],
                    resolution = 1,keep = TRUE)
  
  # We can pull depth for each data point
  xy.bathy <- get.depth(b,xy_ll,locator=F) %>% 
    rename(bathy=depth)
  
  # return the computed depths
  return(xy.bathy)
}

# Rotate map?? Just for visualizing facetted plots

# rotate function (see here: https://r-spatial.github.io/sf/articles/sf3.html#affine-transformations
rot <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)

rotate_sf <- function(sf_df,a){
  sf_df %>% 
    mutate(geom_rot = st_geometry(.)*rot(a)) %>%
    st_drop_geometry() %>%
    rename(geometry = geom_rot) %>%
    st_set_geometry("geometry")
}

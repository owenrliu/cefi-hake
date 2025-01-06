# Using the results of the previous scripts,
# build GLORYS-derived covariates for FEAT survey data and projection grid
library(tidyverse)
library(marmap)
library(sf)
library(nngeo)
library(here)
library(viridis)
library(stars)
options(dplyr.summarise.inform = FALSE)

pt <- theme_minimal()+
  theme(panel.border = element_rect(fill=NA,color='black'))
theme_set(pt)

library(rnaturalearth)
library(tictoc) # for timing code

# coastline, joined US and Canada
coast <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada','Idaho','Montana'))
coastcn <- ne_countries(country="Canada",scale=50,returnclass='sf')
coast <- st_union(coast,coastcn)

#### Datasets ####
# FEAT observations first
feat <- read_rds(here('data','CONFIDENTIAL','un-kriged_aged_output_allyears.rds'))
# template 5km projection grid (from 00_make_FEAT_grid.R)
feat_gr <- read_rds(here('data','grids','FEAT_5km_grid.rds'))

# helper scripts
source(here('code','helper_fxns.R'))

#### Bathymetry ####

# Bathymetry for the FEAT observations can be extracted from marmap
feat_xy <- feat %>% 
  distinct(year,Lat,Lon) %>% 
  st_as_sf(coords=c("Lon","Lat"),crs=4326,remove=F)
feat_obs_bathy <- get_bathy_pts(feat_xy) %>% 
  mutate(year=feat_xy$year)
feat <- feat %>% 
  left_join(feat_obs_bathy,by=join_by(year,Lon==lon,Lat==lat))
# plot to see how this worked
feat %>%
  filter(year==2019) %>% #arbitrary
  ggplot(aes(Lon,Lat,color=bathy))+
    geom_point(size=0.5)+
    coord_equal()
feat %>%
  ggplot(aes(-bathy))+ # swap scale from negative to positive depths
  geom_density(fill='gray50')+
  labs(x="Bottom Depth (m)",y="Observations")

feat %>%
  ggplot(aes(log(-bathy)))+
  geom_density(fill='gray50')+
  labs(x="Log(Bottom Depth)",y="Observations")

# Distance from the 200m isobath
# this is 200m isobath from the bathymetry downloaded previously for the FEAT grid
shelf_break <- read_rds(here('data','grids','200m_isobath.rds')) %>% 
  st_as_sf(coords=c("x","y"),crs=4326) %>%
  # projected coordinate system
  st_transform(st_crs(feat_gr)) %>% 
  mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2])

feat_xy_utm <- feat_xy %>% st_transform(st_crs(feat_gr))

# then, find nearest neighbors (i.e., minimum distance to the 200m isobath) for the FEAT obs
nn_samp <- st_nn(feat_xy_utm,shelf_break,k=1,returnDist = T)
nn_samp_dists <- nn_samp %>% pluck("dist") %>% unlist()

feat_obs_dist_shelf <- feat_xy %>% 
  # distance calculated by st_nn() is in meters, so convert to km here
  mutate(dist_shelf_km=nn_samp_dists/1000) %>% 
  st_set_geometry(NULL)
feat <- feat %>% 
  left_join(feat_obs_dist_shelf,by=join_by(year,Lon,Lat))
# we could use bathy to indicate a sign for the dist-shelf (inshore vs. offshore)
feat <- feat %>% 
  # if point is INSHORE from 200m (i.e., shallower), switch the shelf distance sign to negative
  mutate(dist_shelf_2=if_else(bathy< -200,dist_shelf_km,-dist_shelf_km))

# use these two new derived values to compare
# plot
feat_bathy_shelf <- feat %>% 
  sample_n(10000) %>% 
  ggplot(aes(bathy,dist_shelf_km))+
  geom_point()+
  coord_flip()+
  geom_vline(xintercept=-200,color='red',linetype=2)+
  labs(x="Depth",y="Distance to Shelf")
feat_bathy_shelf
# looks about right
# and with our sign switch:
feat_bathy_shelf2 <- feat %>% 
  sample_n(10000) %>% 
  ggplot(aes(bathy,dist_shelf_2))+
  geom_point()+
  coord_flip()+
  geom_vline(xintercept=-200,color='red',linetype=2)+
  labs(x="Depth",y="Distance to Shelf")
feat_bathy_shelf2

# and now just the overall distribution of distance
feat %>%
  ggplot(aes(dist_shelf_2))+ # swap scale from negative to positive depths
  geom_density(fill='gray50')+
  labs(x="Distance to Shelf (km)",y="Observations")
feat %>% 
  sample_n(10000) %>% 
  ggplot(aes(Lon,Lat,color=dist_shelf_2))+ # swap scale from negative to positive depths
  geom_point()+
  scale_color_gradient2()+
  coord_equal()

#### GLORYS Covariates ####
# Write a function to use the GLORYS summarized seasonal output (from 03_make_seasonal_glorys)
# The function will, for an input year, variable, and depth range, return the seasonal spatial fields
# For now, we extract the values at either
# 1. locations of the FEAT 5km projection grid (in UTM zone 10)
# 2. survey locations of the observed FEAT data (if there was a survey in that year)

# variables to choose from, and their file names
phys4dfls <- list.files(here('data','glorys','phys','seasonal'),full.names = T) %>% str_subset("4d")
phys4dvars <- phys4dfls[1] %>% read_rds() %>% 
  dplyr::select(-any_of(c("latitude","longitude","season","depth"))) %>% 
  names()
phys4dtbl <- crossing(varname=phys4dvars,year=1995:2023) %>% 
  mutate(file_name=rep(phys4dfls,length(phys4dvars)))%>% mutate(grid_type='phys')
phys3dfls <-list.files(here('data','glorys','phys','seasonal'),full.names = T) %>% str_subset("3d")
phys3dvars <- phys3dfls[1] %>% read_rds() %>% 
  dplyr::select(-any_of(c("latitude","longitude","season","depth"))) %>% 
  names() 
phys3dtbl <- crossing(varname=phys3dvars,year=1995:2023) %>% 
  mutate(file_name=rep(phys3dfls,length(phys3dvars)))%>% mutate(grid_type='phys')
bgcfls <- list.files(here('data','glorys','bgc','seasonal'),full.names = T)
bgcvars <- bgcfls[1] %>% 
  read_rds() %>% dplyr::select(-any_of(c("latitude","longitude","season","depth"))) %>% 
  names()
bgctbl <- crossing(varname=bgcvars,year=1995:2023) %>% 
  mutate(file_name=rep(bgcfls,length(bgcvars)))%>% mutate(grid_type='bgc')

# all year/file name combos
all_v_fls <- bind_rows(phys4dtbl,phys3dtbl,bgctbl)
all_v <- unique(all_v_fls$varname)
all_v

# all possible depth slices
glorys_depths_phys <- read_rds(here('data','glorys','phys','filtered','glorys_physics_2019_spatial_filtered.rds')) %>% 
  distinct(depth)
glorys_depths_bgc <- read_rds(here('data','glorys','bgc','filtered','glorys_bgc_2019_spatial_filtered.rds')) %>% 
  distinct(depth)

# template GLORYS grids
## physics first
# take an example year
glor_phys_stars <- read_stars(here('data','glorys','phys','raw','glorys_physics_1995.nc'),sub="thetao") %>% 
  st_set_crs(4326)
# take the first time/depth slice
slc <- slice(glor_phys_stars,index=1,along='time') %>% slice("depth",1)
# take a spatial subset/crop, using the FEAT footprint
cropped <- slc[feat_foot_ll]
glor_phys_templ <- cropped %>% 
  st_as_sf() %>% 
  mutate(gid=row_number()) %>% 
  select(gid)
# plot, with arbitrary value (this will get filled with what we extract from GLORYS)
ggplot(glor_phys_templ,aes(fill=gid))+geom_sf(col=NA)

## for BGC
# take an example year
glor_bgc_stars <- read_stars(here('data','glorys','bgc','raw','glorys_bgc_1995.nc'),sub="chl") %>% 
  st_set_crs(4326)
# take the first time/depth slice
slc <- slice(glor_bgc_stars,index=1,along='time') %>% slice("depth",1)
# take a spatial subset/crop, using the FEAT footprint
cropped <- slc[feat_foot_ll]
glor_bgc_templ <- cropped %>% 
  st_as_sf() %>% 
  mutate(gid=row_number()) %>% 
  select(gid)

# Matching function
# for any of the variables in all_v, we can use this function
# depths should be positive, in meters
make_glorys_covariate <- function(glorys_yr,variable='mean_thetao',
                                  # depth minimum and maximum (remember that the shallowest GLORYS depth is ~0.5m
                                  depth_min,depth_max=depth_min,
                                  # if a depth range is selected, take a mean or a sum (integral)?
                                  type="mean",
                                  # return FEAT survey data locations (survey) or FEAT projection grid (grid)?
                                  return_what='survey'
                                  ){
  
  # filter for the right GLORYS file using the table we made
  sub_glorys <- all_v_fls %>%
    # find the right year and the right variable
    filter(year==glorys_yr,varname==variable) %>% 
    pull(file_name) %>% 
    # load
    read_rds() %>% 
    # filter depth
    filter(depth >= depth_min, depth <= depth_max) %>% 
    dplyr::select(any_of(c(variable,"latitude","longitude","season","depth")))
  
  # make final calculation of the variable of interest
  if(type=="mean" & depth_max!=depth_min){
    sub_glorys <- sub_glorys %>% 
      group_by(latitude,longitude,season) %>% 
      summarise(across(all_of(variable),mean,.names="{.col}"))
  } else if(type=="sum" & depth_max!=depth_min){
    sub_glorys <- sub_glorys %>% 
      group_by(latitude,longitude,season) %>% 
      summarise(across(all_of(variable),sum,.names="{.col}"))
  }
  
  # cast to wide (one variable for each season)
  sub_glorys_wide <- sub_glorys %>% ungroup() %>% 
    pivot_wider(names_from="season",values_from=variable,names_prefix = paste0(variable,"_"))
  
  # match to the outputs
  # glorys points, spatial
  pts_glorys_sf <- sub_glorys_wide %>% 
    st_as_sf(coords=c('longitude','latitude'),crs=4326)
  
  # which template GLORYS grid to use?
  grid_type <- all_v_fls$grid_type[match(variable,all_v_fls$varname)]
  if(grid_type=="phys") glorys_sf=glor_phys_templ
  if(grid_type=="bgc") glorys_sf=glor_bgc_templ
  
  # fill the template with the GLORYS data
  glorys_sf_filled <- glorys_sf %>% 
    st_join(pts_glorys_sf)

  if(return_what=="survey"){
    # FEAT observation points
    sub_feat <- feat %>% filter(year==glorys_yr)
    feat_pts <- sub_feat %>%
      st_as_sf(coords=c("Lon","Lat"),crs=4326)
    # joined to GLORYS
    out <- st_join(feat_pts,glorys_sf_filled)
  }
  if(return_what=="grid"){
    # FEAT projection grid points
    feat_gr_pts <- feat_gr %>% 
      st_transform(4326) %>% 
      st_centroid()
    
    out <- feat_gr_pts %>%
      st_join(glorys_sf_filled) %>% 
      st_set_geometry(NULL) %>% 
      right_join(feat_gr) %>% 
      mutate(year=glorys_yr) %>% 
      st_as_sf()
  }
  
  return(out)
  }

# Try
# sst
test <- make_glorys_covariate(glorys_yr=2005,variable="mean_thetao",
                              depth_min=min(glorys_depths_phys$depth),type = "mean",
                              return_what="grid")
test %>%
  pivot_longer(contains("thetao")) %>% 
  ggplot()+
  geom_sf(aes(fill=value),col=NA)+
  facet_wrap(~name)+
  scale_fill_viridis(option="turbo")

# Try making a full time series?
tic("making sst")
sst_gridded <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_thetao",
                                                               depth_min=min(glorys_depths_phys$depth),type = "mean",
                                                               return_what="grid")) %>% 
  list_rbind()
toc()
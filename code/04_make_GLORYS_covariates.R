# Using the results of the previous scripts,
# build GLORYS-derived covariates for FEAT survey data and projection grid
library(tidyverse)
library(marmap)
library(sf)
library(nngeo)
library(here)
library(viridis)
pt <- theme_minimal()+
  theme(panel.border = element_rect(fill=NA,color='black'))
theme_set(pt)

#### DATASETS ####
# FEAT observations first
feat <- read_rds(here('data','CONFIDENTIAL','un-kriged_aged_output_allyears.rds'))
# template 5km projection grid (from 00_make_FEAT_grid.R)
feat_gr <- read_rds(here('data','grids','FEAT_5km_grid.rds'))

# matching keys for the two datasets above
feat_match <- read_rds(here('data','feat_obs_glorys_matchkey.rds'))
feat_gr_match <- read_rds(here('data','feat_gr_glorys_matchkey.rds'))
# correct GLORYS grids for the matching keys above
glor_phys_xy <- list.files(here('data','glorys','phys','filtered'),full.names = T)[1] %>% 
  read_rds() %>% 
  distinct(longitude,latitude) %>%
  mutate(glorID_phys=row_number())
glor_bgc_xy <- list.files(here('data','glorys','bgc','filtered'),full.names = T)[1] %>% 
  read_rds() %>% 
  distinct(longitude,latitude) %>% 
  mutate(glorID_bgc=row_number())

# helper scripts
source(here('code','helper_fxns.R'))

#### Covariate Construction ####
## Here is where we build covariates to then ##
## match to the FEAT data. Just a couple of  ##
## examples for now, we can probably expand  ##
## this later

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

## SEASONAL MEAN PHYSICS ###
# calculate seasonal means for covariates
# for now, at all depths
# Spring: MAM
# Summer: JJA
# Fall: SON
# Winter: DJF
# For a given year, winter will be from the PREVIOUS year (i.e. the year before plus January of the focal year)
glorys_depths <- read_rds(here('data','glorys','phys','filtered','glorys_physics_2019_spatial_filtered.rds')) %>% 
  distinct(depth)
gdf <- read_rds(here('data','glorys','phys','filtered','glorys_physics_2019_spatial_filtered.rds'))
summarise_seasonal_glorys_phys <- function(glorys_yr){
  # previous year
  gdf1 <- read_rds(here('data','glorys','phys','filtered',paste0('glorys_physics_',glorys_yr-1,'_spatial_filtered.rds')))
  # current year
  gdf2 <- read_rds(here('data','glorys','phys','filtered',paste0('glorys_physics_',glorys_yr,'_spatial_filtered.rds')))
  # summarise
  out <- bind_rows(gdf1,gdf2) %>% 
    mutate(year=year(date),mth=month(date)) %>% 
    # filter out December of the focal year, which will be included in next year's data
    # and filter out JF of the previous year, which was included in last year's data
    filter(!(year==glorys_yr&mth==12)) %>% 
    filter(!(year==(glorys_yr-1)& (mth%in%c(1,2)))) %>%
    #add seasons
    mutate(season=case_when(
      mth %in% c(3:5) ~ "spring",
      mth %in% c(6:8) ~ "summer",
      mth %in% c(9:11) ~ "autumn",
      mth %in% c(12,1,2) ~ "winter"
      )) %>% 
    left_join(glor_phys_xy,by=join_by(longitude,latitude)) %>% 
    # now, group by season and summarise
    group_by(glorID_phys,latitude,longitude,depth,season) %>% 
    summarise(across(all_of(c('so','thetao','uo','vo')),mean,.names="mean_{.col}")) %>% 
    ungroup()
  
  out
}
# test
seasonal_2019 <- summarise_seasonal_glorys_phys(glorys_yr=2019)
# test plot, sea surface temperature
seasonal_2019 %>% 
  st_as_sf(coords=c("longitude","latitude"),crs=4326) %>% 
  st_transform(st_crs(feat_gr)) %>% 
  filter(depth==unique(seasonal_2019$depth)[1]) %>% 
  ggplot(aes(color=mean_thetao))+
  geom_sf(size=0.25)+
  scale_color_viridis(option="turbo")+
  facet_wrap(~season,nrow=2)

# test plot, northward velocity between ~150 and 250m depths
seasonal_2019 %>% 
  filter(depth %in% glorys_depths$depth[25:28]) %>% 
  group_by(glorID_phys,season,latitude,longitude) %>% 
  summarise(northvel=mean(mean_vo)) %>% 
  st_as_sf(coords=c("longitude","latitude"),crs=4326) %>% 
  st_transform(st_crs(feat_gr)) %>% 
  ggplot(aes(color=northvel))+
  geom_sf(size=0.25)+
  scale_color_viridis(option="turbo")+
  facet_wrap(~season,nrow=2)

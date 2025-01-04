# Goal: extract data from already-downloaded GLORYS ncdf files (01_download_GLORYS.R)
library(ncdf4) 
library(tidync)
library(data.table)
library(tidyverse)
library(here)
library(geosphere)
library(sf)
library(rnaturalearth)
library(tictoc) # for timing code

# coastline, joined US and Canada
coast <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada','Idaho','Montana'))
coastcn <- ne_countries(country="Canada",scale=50,returnclass='sf')
coast <- st_union(coast,coastcn)

# Let's try extracting data from the large GLORYS netCDFs, but with coarse filtering to try to avoid excess data

# get the GLORYS grids from one of the nc files
# this is for physics, which is higher resolution than the biogeochemistry
gr <- tidync(here('data','glorys','phys','raw','glorys_physics_1995.nc'))
# this is for the bgc
gr_bgc <- tidync(here('data','glorys','bgc','raw','glorys_bgc_1995.nc'))

# Footprint of the FEAT survey data, as a polygon (created in 00_make_FEAT_grid.R)
feat_foot <- read_rds(here('data','grids','FEAT_footprint.rds'))
# in lat/lon coords
feat_foot_ll <- feat_foot %>% st_transform(4326)
# bounding box for the FEAT footprint
bb <- st_bbox(feat_foot_ll)

## GLORYS PHYSICS
# grab first time/depth slice, just to get grid
# note, this first part is only for the 4D variables (i.e. indexed by time, lat, lon, depth)
# other variables are stored in different grids because they are depth-invariant:
# i.e., mixed layer depth, ssh, and bottom temperature are 3D only (lat/lon/time)
slice1 <- gr %>% hyper_filter(time=index==2) %>% hyper_tibble()
# converted from dataframe to sf object
slice1_sf <- slice1 %>% st_as_sf(coords=c("longitude","latitude"),crs=4326)

# plot, relative to FEAT
ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=slice1_sf %>% sample_n(10000))+
  geom_sf(data=feat_foot,color='red',fill=NA)+
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])

# see, we only need a small subset of these values. Let's do a simple overlay and 
# keep track of which GLORYS points overlap with the FEAT footprint
glorys_in_feat <- slice1_sf %>% 
  # make a grid index. if we're careful, this should be consistent for every year of data
  mutate(glorID=row_number()) %>% 
  st_filter(feat_foot_ll)

# these are the grid indices for the rows of GLORYS data we will extract for each year/month
glorIDs <- glorys_in_feat$glorID

# lets try plotting again
ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=glorys_in_feat,size=0.25)+
  geom_sf(data=feat_foot,color='red',fill=NA)+
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])

# the spatial filter/overlay is slow for large data, so BECAUSE GLORYS is consistent in its grid, 
# we should be able to use non-spatial indexing to filter instead

# Do this for all year/months
# 1. open a year's ncdf
# 2. for each month, extract all data
# 3. filter by row number using the glorIDs from above
# then 
# 4. recombine across months

extract_glorys_feat_footprint <- function(nc,mth){
  # extract the data
  temp <- tidync(nc) %>% hyper_filter(time=index==mth) %>% hyper_tibble()
  # filter using the grid IDs above
  tempsub <- temp %>% slice(glorIDs)
  tempsub %>% 
    #convert to meaningful dates (instead of hours since 1950-01-01)
    mutate(date=as_datetime("1950-01-01 00:00:00")+hours(time))
}

# glorys physics files lists
fls <- list.files(here('data','glorys','phys','raw'),full.names = T)

tic("trying first year, first month")
test <- extract_glorys_feat_footprint(fls[1],mth=1)
toc()

tic("trying first year, all months")
test <- purrr::map(1:12,extract_glorys_feat_footprint,nc=fls[1]) %>% 
  list_rbind()
toc()
# ~10s for one year

# finally, apply to all years (except 2021 is a special case)
for(i in c(1:26,29,30)){
  tic(paste("processing",fls[i]))
  extr <- purrr::map(1:12,extract_glorys_feat_footprint,nc=fls[i]) %>% 
    list_rbind()
  fn_out <- fls[i] %>% str_replace(".nc","_spatial_filtered.rds") %>% str_replace("raw","filtered") 
  write_rds(extr,fn_out)
  toc()
}

# 2021 is weird because it splits two datasets
# the first has january to july, the second, august to december
tic(paste("processing",fls[27],"and",fls[28]))
extr2021.1 <- purrr::map(1:6,extract_glorys_feat_footprint,nc=fls[27]) %>% 
    list_rbind()
extr2021.2 <- purrr::map(1:6,extract_glorys_feat_footprint,nc=fls[28]) %>% 
  list_rbind()
extr <- bind_rows(extr2021.1,extr2021.2)
fn_out <- fls[27] %>% str_replace(".nc","_spatial_filtered.rds") %>% str_replace("raw","filtered")
write_rds(extr,fn_out)
toc()

### 3D VARIABLES #
# for the 3D variables, the process can be the same
# find gridIDs
slice3d <- gr %>% 
  activate("D2,D1,D3") %>% 
  hyper_filter(time=index==1) %>% 
  hyper_tibble()

# these are the row numbers (i.e., the indices) of the data we will extract for each year/month
slice3d_sf <- slice1 %>% st_as_sf(coords=c("longitude","latitude"),crs=4326)
glorys3d_in_feat <- slice3d_sf %>% 
  mutate(glor3dID=row_number()) %>% 
  st_filter(feat_foot_ll)
glor3dIDs <- glorys3d_in_feat$glor3dID

# same idea, extract filtered data for each year/month
extract_glorys3D_feat_footprint <- function(nc,mth){
  # extract the data
  temp <- tidync(nc) %>% activate("D2,D1,D3") %>%
    hyper_filter(time=index==mth) %>% hyper_tibble()
  # filter using the grid IDs above
  tempsub <- temp %>% slice(glor3dIDs)
  tempsub %>% 
    #real times (hours since 1950-01-01)
    mutate(date=as_datetime("1950-01-01 00:00:00")+hours(time))
}


for(i in c(1:26,29,30)){
  tic(paste("processing",fls[i]))
  extr <- purrr::map(1:12,extract_glorys3D_feat_footprint,nc=fls[i]) %>% 
    list_rbind()
  fn_out <- fls[i] %>% str_replace(".nc","_3d_filtered.rds") %>% str_replace("raw","filtered")
  write_rds(extr,fn_out)
  toc()
}

# for 2021
tic(paste("processing",fls[27],"and",fls[28]))
extr2021.1 <- purrr::map(1:6,extract_glorys3D_feat_footprint,nc=fls[27]) %>% 
  list_rbind()
extr2021.2 <- purrr::map(1:6,extract_glorys3D_feat_footprint,nc=fls[28]) %>% 
  list_rbind()
extr <- bind_rows(extr2021.1,extr2021.2)
fn_out <- fls[27] %>% str_replace(".nc","_3d_filtered.rds") %>% str_replace("raw","filtered")
write_rds(extr,fn_out)
toc()

## GLORYS BGC, with a different resolution than the physics
# grab first time/depth slice, just to get grid
slice_bgc <- gr_bgc %>% hyper_filter(time=index==1) %>% hyper_tibble()
# converted from dataframe to sf object
slice_bgc_sf <- slice_bgc %>% st_as_sf(coords=c("longitude","latitude"),crs=4326)

# plot, relative to FEAT
ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=slice_bgc_sf %>% sample_n(10000))+
  geom_sf(data=feat_foot,color='red',fill=NA)+
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])

glorys_bgc_in_feat <- slice_bgc_sf %>% 
  # make a grid index. if we're careful, this should be consistent for every year of data
  mutate(glorID_bgc=row_number()) %>% 
  st_filter(feat_foot_ll)

# these are the grid indices for the rows of GLORYS data we will extract for each year/month
glorIDs_bgc <- glorys_bgc_in_feat$glorID_bgc

# lets try plotting again
ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=glorys_bgc_in_feat,size=0.25)+
  geom_sf(data=feat_foot,color='red',fill=NA)+
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])

# apply to all bgc year/months
extract_glorys_bgc_feat_footprint <- function(nc,mth){
  # extract the data
  temp <- tidync(nc) %>% hyper_filter(time=index==mth) %>% hyper_tibble()
  # filter using the grid IDs above
  tempsub <- temp %>% slice(glorIDs_bgc)
  tempsub %>% 
    #convert to meaningful dates (instead of hours since 1950-01-01)
    mutate(date=as_datetime("1950-01-01 00:00:00")+hours(time))
}

# apply to all years
fls_bgc <- list.files(here('data','glorys','bgc','raw'),full.names = T)

for(i in 1:length(fls_bgc)){
  tic(paste("processing",fls_bgc[i]))
  extr <- purrr::map(1:12,extract_glorys_bgc_feat_footprint,nc=fls_bgc[i]) %>% 
    list_rbind()
  fn_out <- fls_bgc[i] %>% str_replace(".nc","_spatial_filtered.rds") %>% str_replace("raw","filtered") 
  write_rds(extr,fn_out)
  toc()
}

## At the end of this script, we should have extracted and coarsely spatially filtered both the physical and biogeochemical data
## In the next script, we can actually join these data to both our observations and to the projection grid
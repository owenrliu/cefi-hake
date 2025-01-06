## Summarize GLORYS data as seasonal spatial means

# Spring: MAM
# Summer: JJA
# Fall: SON
# Winter: DJF
# For a focal year, winter will include December of the previous year,
# and not December of the focal year

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
library(tictoc)
# Footprint of the FEAT survey data, as a polygon (created in 00_make_FEAT_grid.R)
feat_foot <- read_rds(here('data','grids','FEAT_footprint.rds'))
# in lat/lon coords
feat_foot_ll <- feat_foot %>% st_transform(4326)
bb_ll <- st_bbox(feat_foot_ll)


#### GLORYS Template Grids####
## Make template polygon grids for physics and bgc data ##

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
  mutate(thetao=NA) %>% 
  rename(val=thetao) %>% 
  mutate(gid=row_number())
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
  mutate(chl=NA) %>% 
  rename(val=chl) %>% 
  mutate(gid=row_number())
# this should look the same, with fewer cells (blockier/lower resolution)
ggplot(glor_bgc_templ,aes(fill=gid))+geom_sf(col=NA)
### seasonal mean GLORYS ####
# calculate seasonal means for covariates
glorys_depths_phys <- read_rds(here('data','glorys','phys','filtered','glorys_physics_2019_spatial_filtered.rds')) %>% 
  distinct(depth)
glorys_depths_bgc <- read_rds(here('data','glorys','bgc','filtered','glorys_bgc_2019_spatial_filtered.rds')) %>% 
  distinct(depth)

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
    mutate(season=factor(season,levels=c("winter","spring","summer","autumn"))) %>% 
    group_by(latitude,longitude,depth,season) %>% 
    summarise(across(all_of(c('so','thetao','uo','vo')),mean,.names="mean_{.col}")) %>% 
    ungroup()
  
  out
}

# for the only 3d physical variables
summarise_seasonal_glorys_3dphys <- function(glorys_yr){
  # previous year
  gdf1 <- read_rds(here('data','glorys','phys','filtered',paste0('glorys_physics_',glorys_yr-1,'_3d_filtered.rds')))
  # current year
  gdf2 <- read_rds(here('data','glorys','phys','filtered',paste0('glorys_physics_',glorys_yr,'_3d_filtered.rds')))
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
    mutate(season=factor(season,levels=c("winter","spring","summer","autumn"))) %>% 
    group_by(latitude,longitude,season) %>% 
    summarise(across(all_of(c('mlotst','zos','bottomT')),mean,.names="mean_{.col}")) %>% 
    ungroup()
  
  out
}

# similarly, for bgc
summarise_seasonal_glorys_bgc <- function(glorys_yr){
  # previous year
  gdf1 <- read_rds(here('data','glorys','bgc','filtered',paste0('glorys_bgc_',glorys_yr-1,'_spatial_filtered.rds')))
  # current year
  gdf2 <- read_rds(here('data','glorys','bgc','filtered',paste0('glorys_bgc_',glorys_yr,'_spatial_filtered.rds')))
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
    mutate(season=factor(season,levels=c("winter","spring","summer","autumn"))) %>% 
    group_by(latitude,longitude,depth,season) %>% 
    summarise(across(all_of(c("chl", "no3", "nppv", "o2", "ph", "phyc")),mean,.names="mean_{.col}")) %>% 
    ungroup()
  
  out
}
# the above function will produce point features, we can also always join them to our template polygons
match_glorys_pts_polys <- function(glorys_pt_df,grid="phys"){ # choose grid, phys or bgc
  pts_sf <- glorys_pt_df %>% 
    st_as_sf(coords=c('longitude','latitude'),crs=4326)
  if(grid=="phys") glorys_sf=glor_phys_templ
  else if(grid=="bgc") glorys_sf=glor_bgc_templ
  out <- st_join(glorys_sf,pts_sf,join=st_contains)
  out
}

# test, full depth profile
seasonal_2019_phys <- summarise_seasonal_glorys_phys(glorys_yr=2019)
seasonal_2019_3dphys <- summarise_seasonal_glorys_3dphys(glorys_yr=2019)
seasonal_2019_bgc <- summarise_seasonal_glorys_bgc(glorys_yr=2019)

# now, for SST
seasonal_2019_surface <- seasonal_2019_phys %>% 
  filter(depth==glorys_depths_phys$depth[1]) %>% 
  match_glorys_pts_polys()
seasonal_2019_surface

# test plot, sea surface temperature
seasonal_2019_surface %>% 
  ggplot()+
  geom_sf(aes(fill=mean_thetao),col=NA)+
  scale_fill_viridis(option="turbo")+
  facet_wrap(~season,nrow=2)

# test plot, sea surface salinity
seasonal_2019_surface %>% 
  ggplot()+
  geom_sf(aes(fill=mean_so),col=NA)+
  scale_fill_viridis(option="D")+
  facet_wrap(~season,nrow=2)

# test plot, northward velocity between ~150 and 250m depths
# i.e. CA undercurrent
seasonal_2019_CUC <- seasonal_2019_phys %>% 
  filter(depth %in% glorys_depths_phys$depth[25:28]) %>% 
  group_by(season,latitude,longitude) %>% 
  summarise(northvel=mean(mean_vo)) %>% 
  match_glorys_pts_polys() %>% 
  filter(!is.na(northvel))

ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=seasonal_2019_CUC,aes(fill=northvel),col=NA)+
  scale_fill_viridis(option="turbo")+
  geom_sf(data=feat_foot_ll,fill=NA,color='red')+
  facet_wrap(~season,nrow=2)+
  xlim(bb_ll[1],bb_ll[3])+ylim(bb_ll[2],bb_ll[4])

# test plot, chlorophyll, depth-integrated 50m (first 19 depth layers)
seasonal_2019_chl <- seasonal_2019_bgc %>%
  filter(depth %in% glorys_depths_bgc$depth[1:19]) %>%
  group_by(season,latitude,longitude) %>% 
  # integrate (sum) across depths
  summarise(chl_int=sum(mean_chl)) %>% 
  match_glorys_pts_polys(grid = 'bgc')

ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=seasonal_2019_chl,aes(fill=chl_int),col=NA)+
  scale_fill_viridis(option="D")+
  # geom_sf(data=feat_foot_ll,fill=NA,color='red')+
  facet_wrap(~season,nrow=2)+
  xlim(bb_ll[1],bb_ll[3])+ylim(bb_ll[2],bb_ll[4])

# test plot, phytoplankton, depth-integrated 50m (first 19 depth layers)
seasonal_2019_phyc <- seasonal_2019_bgc %>%
  filter(depth %in% glorys_depths_bgc$depth[1:19]) %>%
  group_by(season,latitude,longitude) %>% 
  # integrate (sum) across depths
  summarise(phyc_int=sum(mean_phyc)) %>% 
  match_glorys_pts_polys(grid = 'bgc')

ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=seasonal_2019_phyc,aes(fill=phyc_int),col=NA)+
  scale_fill_viridis(option="D")+
  # geom_sf(data=feat_foot_ll,fill=NA,color='red')+
  facet_wrap(~season,nrow=2)+
  xlim(bb_ll[1],bb_ll[3])+ylim(bb_ll[2],bb_ll[4])

# test plot, dissolved oxygen at 200m (depth slice 31 on BGC)
seasonal_2019_o2 <- seasonal_2019_bgc %>%
  filter(depth == glorys_depths_bgc$depth[31]) %>%
  match_glorys_pts_polys(grid = 'bgc') %>% 
  filter(!is.na(mean_o2))

ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=seasonal_2019_o2,aes(fill=mean_o2),col=NA)+
  scale_fill_viridis(option="D")+
  # geom_sf(data=feat_foot_ll,fill=NA,color='red')+
  facet_wrap(~season,nrow=2)+
  xlim(bb_ll[1],bb_ll[3])+ylim(bb_ll[2],bb_ll[4])+
  labs(title="Oxygen at 200m")

# test plot, MLD
seasonal_2019_mld <- seasonal_2019_3dphys %>%
  match_glorys_pts_polys(grid = 'phys')

ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=seasonal_2019_mld,aes(fill=mean_mlotst),col=NA)+
  scale_fill_viridis(option="C",direction=-1)+
  # geom_sf(data=feat_foot_ll,fill=NA,color='red')+
  facet_wrap(~season,nrow=2)+
  xlim(bb_ll[1],bb_ll[3])+ylim(bb_ll[2],bb_ll[4])+
  labs(title="Mixed Layer Depth")

# test plot, SSH

ggplot()+
  geom_sf(data=coast,fill='gray80')+
  geom_sf(data=seasonal_2019_mld,aes(fill=mean_zos),col=NA)+
  scale_fill_viridis(option="turbo")+
  # geom_sf(data=feat_foot_ll,fill=NA,color='red')+
  facet_wrap(~season,nrow=2)+
  xlim(bb_ll[1],bb_ll[3])+ylim(bb_ll[2],bb_ll[4])+
  labs(title="Sea Surface Height")

# Seems to work; let's generically make and save seasonal summaries for all years
# 4d physics
yrs <- 1995:2023
walk(yrs,function(x){
  tic(paste('seasonal summary',x))
  out <- summarise_seasonal_glorys_phys(x)
  fn_out <- here('data','glorys','phys','seasonal',paste0('seasonal_4d_phys_',x,'.rds'))
  write_rds(out,fn_out)
  toc()
})

# 3d physics
walk(yrs,function(x){
  tic(paste('seasonal summary 3d',x))
  out <- summarise_seasonal_glorys_3dphys(x)
  fn_out <- here('data','glorys','phys','seasonal',paste0('seasonal_3d_phys_',x,'.rds'))
  write_rds(out,fn_out)
  toc()
})

# bgc
walk(yrs,function(x){
  tic(paste('seasonal summary bgc',x))
  out <- summarise_seasonal_glorys_bgc(x)
  fn_out <- here('data','glorys','bgc','seasonal',paste0('seasonal_bgc_',x,'.rds'))
  write_rds(out,fn_out)
  toc()
})

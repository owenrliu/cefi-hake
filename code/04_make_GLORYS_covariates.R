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
# Footprint of the FEAT survey data, as a polygon (created in 00_make_FEAT_grid.R)
feat_foot <- read_rds(here('data','grids','FEAT_footprint.rds'))
# in lat/lon coords
feat_foot_ll <- feat_foot %>% st_transform(4326)

bb <- st_bbox(feat_gr)
coast_utm <- coast %>% st_transform(st_crs(feat_gr)) %>% summarise()

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

# by life stage
# cast to long form (so that each observation/age combo is a row)
# summarise by year
feat_long <- feat %>% 
  filter(wgt>0) %>% 
  mutate(sampnum=row_number()) %>% 
  pivot_longer(`1`:`20`,values_to = "n",names_to="age") %>% 
  mutate(age=as.integer(age)) %>% 
  mutate(binclass=ifelse(age%in%c(1,2,3),"juvenile","adult")) %>% 
  # weight depth by fish (sensu Agostini 2006)
  group_by(binclass,year) %>% 
  summarise(z_weighted=sum(n*bathy)/sum(n)) %>% 
  ungroup()

feat_long %>% 
  ggplot(aes(year,-z_weighted,color=binclass))+
  geom_line(size=1.5)+
  labs(color="Stage",x="Year",y="Bottom Depth Weighted by Hake Abundance")

#### Distance to Shelf ####

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
  mutate(dist_shelf_km=if_else(bathy< -200,dist_shelf_km,-dist_shelf_km))

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
# looks about right (i.e., obs with negative distance to shelf are shallower than 200m)

# and now just the overall distribution of distance
feat %>%
  ggplot(aes(dist_shelf_km))+ # swap scale from negative to positive depths
  geom_density(fill='gray50')+
  labs(x="Distance to Shelf (km)",y="Observations")

# map
p_feat_shelfdist <- feat %>% 
  sample_n(10000) %>% 
  st_as_sf(coords=c("Lon","Lat"),crs=4326) %>% st_transform(st_crs(shelf_break)) %>% 
  ggplot()+
  geom_sf(aes(color=dist_shelf_km),size=0.5)+
  geom_sf(data=shelf_break,size=0.25,color='gray50')+
  geom_sf(data=coast_utm,fill='gray80')+
  scale_color_gradient2()+
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
  labs(color="Distance to\n200m shelf break")
p_feat_gr_shelfdist <- feat_gr %>% 
  st_as_sf() %>% 
  ggplot()+
  geom_sf(aes(fill=dist_shelf_km,col=dist_shelf_km))+
  geom_sf(data=shelf_break,size=0.25,color='gray50')+
  geom_sf(data=coast_utm,fill='gray80')+
  scale_fill_gradient2()+
  scale_color_gradient2()+
  xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
  labs(fill="Distance to\n200m shelf break",col="Distance to\n200m shelf break")

ggsave(here('data','grids','feat_grid_dist_shelf.png'),p_feat_gr_shelfdist,h=8,w=6)

feat_long_shelf <- feat %>% 
  filter(wgt>0) %>% 
  mutate(sampnum=row_number()) %>% 
  pivot_longer(`1`:`20`,values_to = "n",names_to="age") %>% 
  mutate(age=as.integer(age)) %>% 
  mutate(binclass=ifelse(age%in%c(1,2,3),"juvenile","adult")) %>% 
  # weight depth by fish (sensu Agostini 2006)
  group_by(binclass,year) %>% 
  summarise(dist_weighted=sum(n*dist_shelf_km)/sum(n)) %>% 
  ungroup()

feat_long_shelf %>% 
  ggplot(aes(year,dist_weighted,color=binclass))+
  geom_line(size=1.5)+
  geom_hline(yintercept=0,linetype=2)+
  labs(color="Stage",x="Year",y="Shelf Distance Weighted by Hake Abundance")

# save the bathymetry and shelf distance covariates
bathy_shelf_feat <- feat %>% dplyr::select(Lat,Lon,sex,year,bathy,dist_shelf_km)
write_rds(bathy_shelf_feat,here('data','grids','feat_obs_bathy_shelfdist.rds'))

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
                                  # depth minimum and maximum (remember that the shallowest GLORYS depth is ~0.5m)
                                  # easiest to use the glorys_depths dataframes above
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
    read_rds() 
  
  # 3D physical variables don't have depth, so make depth=0 so it works in this flows
  if(variable %in% phys3dvars) sub_glorys=mutate(sub_glorys,depth=0)
  
  sub_glorys <- sub_glorys %>% 
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
      dplyr::select(all_of(c("Lon","Lat","sex","year"))) %>% 
      st_as_sf(coords=c("Lon","Lat"),crs=4326)
    # joined to GLORYS
    out <- st_join(feat_pts,glorys_sf_filled)
  }
  if(return_what=="grid"){
    # FEAT projection grid points
    feat_gr_pts <- feat_gr %>% 
      st_transform(4326) %>% 
      st_centroid()
    
    out <- st_join(feat_gr_pts,glorys_sf_filled) %>% 
      st_set_geometry(NULL) %>% 
      right_join(feat_gr,by = join_by(E_km, N_km, bathy_m, dist_shelf_km)) %>% 
      mutate(year=glorys_yr) %>% 
      st_as_sf()
  }
  
  return(out)
  }

#### Make SST ####

tic("making sst, FEAT grid")
sst_gridded <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_thetao",
                                                               depth_min=min(glorys_depths_phys$depth),type = "mean",
                                                               return_what="grid")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_thetao","sst"),contains("mean_thetao"))
toc() # this took about 1m for the full time series- pretty good!

# make and save plots for later
walk(1995:2023,function(y){
  sst_sub <- sst_gridded %>%
    filter(year==y,!is.na(sst_summer)) %>% 
    st_as_sf() %>% 
    pivot_longer(contains("sst"),values_to="SST") %>% 
    mutate(name=str_replace(name,"sst_","")) %>% 
    mutate(name=factor(name,levels=c("winter","spring","summer","autumn")))
  p <- ggplot()+
    geom_sf(data=coast_utm,fill='gray70')+
    geom_sf(data=sst_sub,aes(fill=SST,col=SST))+
    facet_wrap(~name)+
    xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
    labs(title=paste("SST:",y))+
    scale_fill_viridis(option="turbo",limits=c(3,22),breaks=seq(5,20,by=5))+
    scale_color_viridis(option="turbo",limits=c(3,22),breaks=seq(5,20,by=5))
  ggsave(here('data','glorys','maps',paste0("SST_seasons_",y,".png")),p,w=8,h=12)
})

# SST as a covariate

tic("making sst, FEAT obs")
sst_feat <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_thetao",
                                                               depth_min=min(glorys_depths_phys$depth),type = "mean",
                                                               return_what="survey")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_thetao","sst"),contains("mean_thetao"))
toc() # this took about 28s

write_rds(sst_feat,here('data','glorys','phys','derived covariates','sst_feat_obs.rds'))
write_rds(sst_gridded,here('data','glorys','phys','derived covariates','sst_gridded.rds'))



#### Make California Undercurrent ####

tic("making cuc, FEAT grid")
cuc_gridded <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_vo",
                                                               depth_min=glorys_depths_phys$depth[24],
                                                               depth_max=glorys_depths_phys$depth[29],type = "mean",
                                                               return_what="grid")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_vo","cuc"),contains("mean_vo"))
toc() # this took about 45s

# make and save plots for later
walk(1995:2023,function(y){
  cuc_sub <- cuc_gridded %>%
    filter(year==y,!is.na(cuc_summer)) %>% 
    st_as_sf() %>% 
    pivot_longer(contains("cuc"),values_to="cuc") %>% 
    mutate(name=str_replace(name,"cuc_","")) %>% 
    mutate(name=factor(name,levels=c("winter","spring","summer","autumn")))
  p <- ggplot()+
    geom_sf(data=coast_utm,fill='gray70')+
    geom_sf(data=cuc_sub,aes(fill=cuc,col=cuc))+
    facet_wrap(~name)+
    xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
    labs(title=paste("California Undercurrent:",y))+
    scale_fill_viridis(option="turbo",limits=c(-0.2,0.4),breaks=seq(-0.2,0.4,by=0.2))+
    scale_color_viridis(option="turbo",limits=c(-0.2,0.4),breaks=seq(-0.2,0.4,by=0.2))
  ggsave(here('data','glorys','maps',paste0("CUC_seasons_",y,".png")),p,w=8,h=12)
})

# CUC with latitude
feat_gr_lat <- cuc_gridded %>% filter(year==1995) %>% st_as_sf() %>% st_transform(4326) %>% st_centroid() %>% mutate(lat=st_coordinates(.)[,2]) %>% 
  st_set_geometry(NULL) %>% dplyr::select(N_km,lat)
p_cuc_lat <- cuc_gridded %>%  
  left_join(feat_gr_lat,by=join_by(N_km)) %>% 
  group_by(lat,year) %>% 
  summarise(zonal_cuc_spring=mean(cuc_spring,na.rm=T)) %>% 
  ungroup() %>% 
  ggplot(aes(lat,zonal_cuc_spring))+
  geom_point(size=0.5)+
  geom_vline(xintercept=44,linetype=2)+
  geom_hline(yintercept=0,linetype=2,color='red')+
  geom_smooth(se=F)+
  labs(x="Latitude",y="Zonal Mean Northward Velocity (m/s)")+
  facet_wrap(~year)
ggsave(here('data','glorys','phys','derived covariates','spring_undercurrent_annual_by_lat.png'),p_cuc_lat,w=11,h=8.5)

# location (inshore/offshore) of max flow
p_cuc_loc <- cuc_gridded %>%  
  select(-x) %>% 
  left_join(feat_gr_lat,by=join_by(N_km)) %>%
  filter(!is.na(cuc_spring)) %>% 
  group_by(lat,year) %>% 
  reframe(max_cuc_loc=dist_shelf_km[cuc_spring==max(cuc_spring,na.rm=T)],
          max_cuc_spring=max(cuc_spring,na.rm=T)) %>% 
  ungroup() %>% 
  ggplot(aes(lat,max_cuc_loc))+
  geom_point(size=0.5)+
  geom_vline(xintercept=44,linetype=2)+
  geom_hline(yintercept=0,linetype=2,color='red')+
  geom_smooth(se=F)+
  labs(x="Latitude",y="Distance from Shelf (km) of Maximum Northward Flow")+
  facet_wrap(~year)
ggsave(here('data','glorys','phys','derived covariates','max_spring_undercurrent_dist_shelf.png'),p_cuc_loc,w=11,h=8.5)

# for summer
p_cuc_loc2 <- cuc_gridded %>%  
  select(-x) %>% 
  left_join(feat_gr_lat,by=join_by(N_km)) %>%
  filter(!is.na(cuc_summer)) %>% 
  group_by(lat,year) %>% 
  reframe(max_cuc_loc=dist_shelf_km[cuc_summer==max(cuc_summer,na.rm=T)],
          max_cuc_spring=max(cuc_summer,na.rm=T)) %>% 
  ungroup() %>% 
  ggplot(aes(lat,max_cuc_loc))+
  geom_point(size=0.5)+
  geom_vline(xintercept=44,linetype=2)+
  geom_hline(yintercept=0,linetype=2,color='red')+
  geom_smooth(se=F)+
  labs(x="Latitude",y="Distance from Shelf (km) of Maximum Northward Flow")+
  facet_wrap(~year)

ggsave(here('data','glorys','phys','derived covariates','max_summer_undercurrent_dist_shelf.png'),p_cuc_loc2,w=11,h=8.5)

# CUC as a covariate

tic("making cuc, FEAT obs")
cuc_feat <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_vo",
                                                            depth_min=glorys_depths_phys$depth[24],
                                                            depth_max=glorys_depths_phys$depth[29],type = "mean",
                                                            return_what="survey")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_vo","cuc"),contains("mean_vo"))
toc() # this took about 28s

write_rds(cuc_feat,here('data','glorys','phys','derived covariates','cuc_feat_obs.rds'))
write_rds(cuc_gridded,here('data','glorys','phys','derived covariates','cuc_gridded.rds'))

#### Make MLD ####

# mixed layer depth
tic("making mld, FEAT grid")
mld_gridded <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_mlotst",
                                                               depth_min=0,type = "mean",
                                                               return_what="grid")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_mlotst","mld"),contains("mean_mlotst"))
toc()

# make and save plots for later
walk(1995:2023,function(y){
  mld_sub <- mld_gridded %>%
    filter(year==y,!is.na(mean_mlotst_summer)) %>% 
    st_as_sf() %>% 
    pivot_longer(contains("mld"),values_to="mld") %>% 
    mutate(name=str_replace(name,"mld_","")) %>% 
    mutate(name=factor(name,levels=c("winter","spring","summer","autumn")))
  p <- ggplot()+
    geom_sf(data=coast_utm,fill='gray70')+
    geom_sf(data=mld_sub,aes(fill=mld,col=mld))+
    facet_wrap(~name)+
    xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
    labs(title=paste("Mean Mixed Layer Depth:",y))+
    scale_fill_viridis(option="turbo",limits=c(0,60),breaks=seq(0,60,by=15))+
    scale_color_viridis(option="turbo",limits=c(0,60),breaks=seq(0,60,by=15))
  ggsave(here('data','glorys','maps',paste0("mld_seasons_",y,".png")),p,w=8,h=12)
})

# mld as a covariate
tic("making mld, FEAT obs")
mld_feat <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_mlotst",
                                                               depth_min=0,type = "mean",
                                                               return_what="survey")) %>% 
  list_rbind() %>% 
  # rename
  rename_with(~str_replace_all(.,"mean_mlotst","mld"),contains("mean_mlotst"))
toc()

write_rds(mld_feat,here('data','glorys','phys','derived covariates','mld_feat_obs.rds'))
write_rds(mld_gridded,here('data','glorys','phys','derived covariates','mld_gridded.rds'))


#### Make T_200 ####

# Temperature at 200m
tic("making t200, FEAT grid")
t200_gridded <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_thetao",
                                                               depth_min=glorys_depths_phys$depth[27],type = "mean",
                                                               return_what="grid")) %>% 
  list_rbind() %>% 
  # rename
  rename_with(~str_replace_all(.,"mean_thetao","t200"),contains("mean_thetao"))
toc()

# make and save plots for later
walk(1995:2023,function(y){
  t200_sub <- t200_gridded %>%
    filter(year==y,!is.na(mean_thetao_summer)) %>% 
    st_as_sf() %>% 
    pivot_longer(contains("thetao"),values_to="t200") %>% 
    mutate(name=str_replace(name,"mean_thetao_","")) %>% 
    mutate(name=factor(name,levels=c("winter","spring","summer","autumn")))
  p <- ggplot()+
    geom_sf(data=coast_utm,fill='gray70')+
    geom_sf(data=t200_sub,aes(fill=t200,col=t200))+
    facet_wrap(~name)+
    xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
    labs(title=paste("Temperature at 200m:",y))+
    scale_fill_viridis(option="turbo",limits=c(4,12),breaks=seq(4,12,by=2))+
    scale_color_viridis(option="turbo",limits=c(4,12),breaks=seq(4,12,by=2))
  ggsave(here('data','glorys','maps',paste0("t200_seasons_",y,".png")),p,w=8,h=12)
})

# mld as a covariate

tic("making t200, FEAT obs")
t200_feat <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_thetao",
                                                             depth_min=glorys_depths_phys$depth[27],type = "mean",
                                                            return_what="survey")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_thetao","t200"),contains("mean_thetao"))
toc()

write_rds(t200_feat,here('data','glorys','phys','derived covariates','t200_feat_obs.rds'))
write_rds(t200_gridded,here('data','glorys','phys','derived covariates','t200_gridded.rds'))

#### Make Chlorophyll ####

# Chlorophyll, integrated from 0-50m
tic("making chl, FEAT grid")
chl_gridded <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_chl",
                                                               depth_min=glorys_depths_bgc$depth[1],
                                                               depth_max= glorys_depths_bgc$depth[19],
                                                               type = "sum",
                                                               return_what="grid")) %>% 
  list_rbind() %>% 
  # rename
  rename_with(~str_replace_all(.,"mean_chl","chl"),contains("mean_chl"))
toc()

# make and save plots for later
walk(1995:2023,function(y){
  chl_sub <- chl_gridded %>%
    filter(year==y,!is.na(chl_summer)) %>% 
    st_as_sf() %>% 
    pivot_longer(contains("chl"),values_to="chl") %>% 
    mutate(name=str_replace(name,"chl_","")) %>% 
    mutate(name=factor(name,levels=c("winter","spring","summer","autumn")))
  p <- ggplot()+
    geom_sf(data=coast_utm,fill='gray70')+
    geom_sf(data=chl_sub,aes(fill=chl,col=chl))+
    facet_wrap(~name)+
    xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
    labs(title=paste("50m Integrated Chlorophyll:",y))+
    scale_fill_viridis(option="turbo",limits=c(0,80),breaks=seq(0,80,by=10))+
    scale_color_viridis(option="turbo",limits=c(0,80),breaks=seq(0,80,by=10))
  ggsave(here('data','glorys','maps',paste0("chl_seasons_",y,".png")),p,w=8,h=12)
})

# chl as a covariate

tic("making chl, FEAT obs")
chl_feat <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_chl",
                                                            depth_min=glorys_depths_bgc$depth[1],
                                                            depth_max= glorys_depths_bgc$depth[19],
                                                            type = "sum",
                                                            return_what="survey")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_chl","chl"),contains("mean_chl"))
toc()

write_rds(chl_feat,here('data','glorys','bgc','derived covariates','chl_feat_obs.rds'))
write_rds(chl_gridded,here('data','glorys','bgc','derived covariates','chl_gridded.rds'))

#### Make Cross Shelf Transport ####

# (maybe we won't use this but) find E-W transport
# for now, at the surface

tic("making zonal flow, FEAT grid")
zonal_gridded <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_uo",
                                                               depth_min=glorys_depths_phys$depth[1],
                                                               type = "mean",
                                                               return_what="grid")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_uo","zonal_u"),contains("mean_uo"))
toc() 

# make and save plots for later
walk(1995:2023,function(y){
  zonal_u_sub <- zonal_gridded %>%
    filter(year==y,!is.na(zonal_u_summer)) %>% 
    st_as_sf() %>% 
    pivot_longer(contains("zonal_u"),values_to="zonal_u") %>% 
    mutate(name=str_replace(name,"zonal_u_","")) %>% 
    mutate(name=factor(name,levels=c("winter","spring","summer","autumn"))) %>% 
    #squish ends of the color palette
    mutate(zonal_u=if_else(zonal_u> 0.3,0.3,zonal_u)) %>% 
    mutate(zonal_u=if_else(zonal_u< -0.3,0.3,zonal_u))
  p <- ggplot()+
    geom_sf(data=coast_utm,fill='gray70')+
    geom_sf(data=zonal_u_sub,aes(fill=zonal_u,col=zonal_u))+
    facet_wrap(~name)+
    xlim(bb[1],bb[3])+ylim(bb[2],bb[4])+
    labs(title=paste("Zonal Transport:",y))+
    scale_fill_viridis(option="turbo",limits=c(-0.3,0.3),breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3))+
    scale_color_viridis(option="turbo",limits=c(-0.3,0.3),breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3))
  ggsave(here('data','glorys','maps',paste0("zonal_u_seasons_",y,".png")),p,w=8,h=12)
})

# Zonal u as a covariate

tic("making zonal_u, FEAT obs")
zonal_u_feat <- map(1995:2023,function(x) make_glorys_covariate(glorys_yr=x,variable="mean_uo",
                                                                depth_min=glorys_depths_phys$depth[1],
                                                                type = "mean",
                                                                return_what="survey")) %>% 
  list_rbind()%>% 
  # rename
  rename_with(~str_replace_all(.,"mean_uo","zonal_u"),contains("mean_uo"))
toc() #

write_rds(zonal_u_feat,here('data','glorys','phys','derived covariates','zonal_u_feat_obs.rds'))
write_rds(zonal_gridded,here('data','glorys','phys','derived covariates','zonal_u_gridded.rds'))

# Go Get GLORYS data and download it, using python workaround
# from https://github.com/pepijn-devries/CopernicusMarine/issues/42#issuecomment-2080179599
library(tidyverse)
library(sf)
library(here)
library(reticulate)
library(tidync)
library(tictoc)

# install_python() 
use_python("C:/Users/Owen.Liu/AppData/Local/r-reticulate/r-reticulate/pyenv/pyenv-win/versions/3.10.11/python.exe")
virtualenv_create(envname = "thiscopernicus")
virtualenv_install("thiscopernicus", packages = c("copernicusmarine"))
use_virtualenv("thiscopernicus", required = TRUE)
py_install("copernicusmarine")

cm <- (import("copernicusmarine"))
# cm$login() # I was having trouble with this, so I did it in the R Terminal with cmd 'copernicusmarine login'
# then i could put in my credentials, which saved to my registry

# copernicusmarine subset -i cmems_mod_glo_phy_my_0.083deg_P1M-m -x -126 -X -115 -y 32 -Y 50 -z 120. -Z 150. -v thetao -t 1994-01-01 -T 2024-08-20 -o ./copernicus-data -f glorys_data_for_hake_temperature.nc
# copernicusmarine subset -i cmems_mod_glo_phy_myint_0.083deg_P1M-m -x -126 -X -115 -y 32 -Y 50 -z 120. -Z 150. -v thetao -t 1994-01-01 -T 2024-08-20 -o ./copernicus-data -f glorys_data_for_hake_temperature2.nc

# bounding box for the FEAT footprint (including Canada), in lat/lon
bb <- read_rds(here('data','grids','FEAT_5km_grid.rds')) %>%
  st_transform(4326) %>% 
  st_bbox()
# xmin       ymin       xmax       ymax 
#-139.04065   32.44635 -117.23297   58.37017 

# FEAT data with survey locations and dates
feat_times <- read_rds(here('data','CONFIDENTIAL','un-kriged_aged_output_allyears.rds')) %>% 
  pull(year) %>% 
  range()
startyr <- feat_times[1]
endyr <- feat_times[2]

## DOWNLOAD PHYSICS DATA FROM GLORYS HINDCAST
# 1993 to 06/2021: cmems_mod_glo_phy_my_0.083deg_P1M-M
# 07/2021-present: cmems_mod_glo_phy_myint_0.083deg_P1M-m
# variables: (not exhaustive, check back if something important is missing later)
# mixed layer depth mlotst
# salinity so
# temperature thetao
# velocities, uo and vo
# ssh zos
# bottom temperature, bottomT

for(i in startyr:2021){
  tic(paste("Downloading GLORYS: ",i))
  dat <- cm$subset(
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1M-m",
    variables=list("mlotst", "so", "thetao", "uo", "vo", "zos", "bottomT"),
    minimum_longitude=-140,
    maximum_longitude=-117,
    minimum_latitude=32,
    maximum_latitude=59,
    start_datetime=paste0(i,"-01-01T00:00:00"),
    end_datetime=paste0(i,"-12-31T23:59:00"),
    minimum_depth=0.49402499198913574,
    maximum_depth=2000,
    output_filename = paste0("data/glorys/phys/raw/glorys_physics_",i,".nc"),
    force_download = TRUE
  )
  toc()
}

# dataset switch since 07-01-2021
# cmems_mod_glo_phy_myint_0.083deg_P1M
for(i in 2021:endyr){
  tic(paste("Downloading GLORYS: ",i))
  dat <- cm$subset(
    dataset_id="cmems_mod_glo_phy_myint_0.083deg_P1M-m",
    variables=list("mlotst", "so", "thetao", "uo", "vo", "zos", "bottomT"),
    minimum_longitude=-140,
    maximum_longitude=-117,
    minimum_latitude=32,
    maximum_latitude=59,
    start_datetime=paste0(i,"-01-01T00:00:00"),
    end_datetime=paste0(i,"-12-31T23:59:00"),
    minimum_depth=0.49402499198913574,
    maximum_depth=2000,
    output_filename = paste0("data/glorys/phys/raw/glorys_physics_",i,".nc"),
    force_download = TRUE
  )
  toc()
}

## DOWNLOAD BIO-GEO-CHEMICAL DATA FROM GLORYS HINDCAST
# for BGC data, different resolution and different datasets
# 1993 to end of 2022: cmems_mod_glo_bgc_my_0.25deg_P1M-m
# 2023-present: cmems_mod_glo_bgc_myint_0.25deg_P1M-m
# variables: (not exhaustive, check back if something important is missing later)
# chlorophyll
# no3
# nppv
# o2
# ph
# phytoplankton

for(i in startyr:endyr){
  tic(paste("Downloading GLORYS BGC: ",i))
  dat <- cm$subset(
    dataset_id=ifelse(i %in% 1993:2022,"cmems_mod_glo_bgc_my_0.25deg_P1M-m","cmems_mod_glo_bgc_myint_0.25deg_P1M-m"),
    variables=list("chl", "no3", "nppv", "o2", "ph", "phyc"),
    minimum_longitude=-140,
    maximum_longitude=-117,
    minimum_latitude=32,
    maximum_latitude=59,
    start_datetime=paste0(i,"-01-01T00:00:00"),
    end_datetime=paste0(i,"-12-31T23:59:00"),
    minimum_depth=0.49402499198913574,
    maximum_depth=2000,
    output_filename = paste0("data/glorys/bgc/raw/glorys_bgc_",i,".nc"),
    force_download = TRUE
  )
  toc()
}
## END
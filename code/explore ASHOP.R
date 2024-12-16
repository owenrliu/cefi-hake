# Explore ASHOP data
library(tidyverse)
library(here)
library(viridis)
library(sf)
library(rnaturalearth)

# plot theme
pt <- theme_minimal()+theme(panel.border=element_rect(color="black",fill=NA))
theme_set(pt)

# backgroundf map
# spatial background map
# load west cost land for mapping
coast <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada')) %>%
  st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=km")


dat <- readxl::read_xlsx(here('data','CONFIDENTIAL','2002-2022_Haul_Catch_Hunsicker_062023.xlsx'))
glimpse(dat)

dat_age  <- readxl::read_xlsx(here('data','CONFIDENTIAL','Hake_ages_1991-2022_062323.xlsx'))

unique(dat$COMMON_NAME)

# obs by month and day of year
dat <- dat %>% 
  mutate(year=year(DEPLOYMENT_DATE),month=month(DEPLOYMENT_DATE),doy=yday(DEPLOYMENT_DATE))

obs_mon <- dat %>% 
  count(month) %>% 
  ggplot(aes(month,n))+
  geom_col()+
  scale_x_continuous(limits=c(1,12),breaks=seq(1,12,by=1))
obs_mon

obs_doy <- dat %>% 
  group_by(doy) %>% 
  summarise(n=n(),totwt=sum(EXTRAPOLATED_WEIGHT_KG,na.rm=T)) %>% 
  ggplot(aes(doy,n,size=totwt/1000))+
  geom_point()
obs_doy

# Map it!
dat_sf <- dat %>% 
  st_as_sf(coords=c("LONDD_END","LATDD_END"),crs=4326) %>% 
  st_transform(26910)
bbox <- st_bbox(dat_sf)

 
ggplot()+
  geom_sf(data=dat_sf)+
  geom_sf(data=coast,shape=4)+
  facet_wrap(~month)+
  scale_color_viridis()+
  coord_sf(xlim=c(bbox[1],bbox[3]),ylim=c(bbox[2],bbox[4]))

## Age/Length data
# Bulk length-freq
dat_age %>% 
  ggplot(aes(LENGTH))+
  geom_histogram(bins=15)
dat_age %>% 
  ggplot(aes(AGE))+
  geom_histogram(bins=10)
dat_age %>% 
  ggplot(aes(LENGTH,fill=SEX))+
  geom_histogram(bins=15,position='dodge')+
  facet_wrap(~YEAR)

# mean size year lat
dat_age %>% 
  filter(SEX !="U") %>% 
  mutate(doy=yday(HAUL_DATE)) %>% 
  mutate(doybin=cut(doy,breaks=c(50,100,150,200,250,300,365))) %>% 
  group_by(doybin,SEX) %>% 
  summarise(meanlen=mean(LENGTH,na.rm=T)) %>%
  ggplot(aes(doybin,meanlen,fill=SEX))+
  geom_col(position='dodge')

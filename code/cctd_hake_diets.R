#https://oceanview.pfeg.noaa.gov/cctd/
library(tidyverse)
library(here)
# annoying because the 2nd line of the csv screws with the parsing
all_content = readLines(here('data','cctd','cctd.csv'))
skip_second = all_content[-2]
dat = read.csv(textConnection(skip_second), header = TRUE, stringsAsFactors = FALSE)

hakedat <- dat %>% filter(predator_scientific_name=="Merluccius productus") %>% 
  # thin
  dplyr::select(collection_id,date,month,day,year,latitude,longitude,region,fishing_depth:surface_temp,predator_id:prey_contents,prey_id:prey_size_comments)
glimpse(hakedat)

hakedat %>% 
  count(year) %>% 
  ggplot(aes(year,n))+geom_line()+
  labs(x="Year",y="Observations")

# Unique prey items
length(unique(hakedat$prey_scientific_name)) #214 unique

prey_count <- hakedat %>% 
  count(prey_scientific_name) %>% 
  arrange(desc(n)) 

prey_count%>% 
  ggplot(aes(forcats::fct_reorder(prey_scientific_name,desc(n)),n))+
  geom_col()
# take the top 20
prey_count%>% 
  slice(1:20) %>% 
  ggplot(aes(forcats::fct_reorder(prey_scientific_name,desc(n)),n))+
  geom_col()+
  labs(x='Prey taxa',y='Observations')+
  coord_flip()+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=0))

# Arrowtooth??
atfdat <- dat %>% filter(predator_scientific_name=="Atheresthes stomias") %>% 
  # thin
  dplyr::select(collection_id,date,month,day,year,latitude,longitude,region,fishing_depth:surface_temp,predator_id:prey_contents,prey_id:prey_size_comments)
glimpse(atfdat)

atfdat %>% 
  count(year) %>% 
  ggplot(aes(year,n))+geom_col()+
  labs(x="Year",y="Observations")

# Unique prey items
length(unique(atfdat$prey_scientific_name)) #35 unique

prey_count <- atfdat %>% 
  count(prey_scientific_name) %>% 
  arrange(desc(n)) 

prey_count%>% 
  ggplot(aes(forcats::fct_reorder(prey_scientific_name,desc(n)),n))+
  geom_col()
# take the top 20
prey_count%>% 
  slice(1:20) %>% 
  ggplot(aes(forcats::fct_reorder(prey_scientific_name,desc(n)),n))+
  geom_col()+
  labs(x='Prey taxa',y='Observations')+
  coord_flip()+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=0))

# Generate world basemap

library(mapdata)
map = map_data("world")

# out antartica and big lakes
map = map[!map$region=="Antarctica",]
lakes1 = grep(unique(map$region),pattern = "Lake",value = T,ignore.case = T)
lakes2 = grep(unique(map$subregion),pattern = "lake",value = T,ignore.case = T)
for (i in lakes1){
  map = map[!map$region==i,]  
}

rm(lakes1)
rm(lakes2)
rm(i)

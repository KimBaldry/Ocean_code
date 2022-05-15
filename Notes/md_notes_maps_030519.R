ggplot() + geom_sf(data = sf::st_as_sf(wm)) + coord_sf(crs = "+proj=ortho +lat_0=-90") + xlim(c(-6e6, 6e6)) + ylim(c(-6e6, 6e6))


ggplot() + geom_sf(data = sf::st_as_sf(fronts), aes(colour = NAME)) +
  geom_sf(data = sf::st_as_sf(wm)) + 
  coord_sf(crs = "+proj=ortho +lat_0=-90")

ggplot() + geom_sf(data = sf::st_as_sf(orsifronts::orsifronts), aes(colour = front)) +
  +   geom_sf(data = sf::st_as_sf(wm)) + 
  +   coord_sf(crs = "+proj=stere +lat_0=-90")

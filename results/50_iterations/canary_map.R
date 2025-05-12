# # # 1. Install/load the necessary packages
# # install.packages(c("rnaturalearth","rnaturalearthhires",
# #                    "geogrid","rmapshaper"))
# # 
# # install.packages(
# #   "terra",
# #   repos = "https://cran.rstudio.com",
# #   type  = "binary"
# # )
# 
# # if you haven’t already installed it
# #install.packages("geogrid")  
# # then load it
# library(geogrid)
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthhires)    # for higher‐res admin units
# library(dplyr)
# library(ggplot2)
# library(geogrid)
# library(rmapshaper)
# 
# # 2. Fetch Spain’s first‐level admin boundaries and filter to the Canaries
# spain_states <- ne_states(country = "Spain", returnclass = "sf")
# canaries <- spain_states %>%
#   filter(name_en %in% 
#            c("Las Palmas", "Santa Cruz de Tenerife")) %>%
#   st_union(by_feature = FALSE) %>%    # if you want them as one poly,
#   st_cast("POLYGON")                  #   or drop these two lines to keep each island.
# 
# # 3. (Optional) Simplify the geometry for a cleaner schematic look
# canaries_simp <- ms_simplify(canaries, keep = 0.1)
# 
# class(canaries_simp)
# 
# # wrap your geometry into an sf data.frame
# can_sf <- st_as_sf(
#   data.frame(id = seq_along(canaries_simp),
#              geometry = canaries_simp),
#   crs = st_crs(canaries_simp)
# )
# 
# new_cells <- calculate_grid(
#   shape     = can_sf,
#   grid_type = "hexagonal",
#   seed      = 42
# )
# 
# # 4. If you want a HEX‐grid schematic instead of true shapes:
# #new_cells <- calculate_grid(canaries_simp, grid_type = "hexagonal", seed = 42)
# canaries_hex <- assign_polygons(can_sf, new_cells)
# 
# # 5. Prepare your island‐level data
# #    Suppose you have a data.frame like this:
# island_data <- tibble(
#   name = c("Tenerife","Gran Canaria","Lanzarote",
#            "Fuerteventura","La Palma","La Gomera","El Hierro"),
#   jaccard    = c(0.34, 0.28, 0.41, 0.22, 0.37, 0.19, 0.25),
#   net_size   = c(120,  95,   70,   50,   40,   35,   30)
# )
# 
# # 6. Join your metrics onto the geometries
# canaries_joined <- can_sf %>%
#   st_as_sf() %>%
#   mutate(name = c("Tenerife","Gran Canaria","Lanzarote",
#                   "Fuerteventura","La Palma","La Gomera","El Hierro")) %>%
#   left_join(island_data, by = "name")
# 
# # 7. Plot with ggplot2, coloring by Jaccard or network size
# ggplot(canaries_joined) +
#   geom_sf(aes(fill = jaccard), color = "grey30") +
#   scale_fill_viridis_c(name = "Jaccard\nindex") +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text  = element_blank(),
#     axis.ticks = element_blank()
#   )
# 
# #—— 0. (Re-)load packages ——#
# # install.packages(c("sf","rnaturalearth","rnaturalearthhires",
# #                    "dplyr","ggplot2","rmapshaper"))
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthhires)
# library(dplyr)
# library(ggplot2)
# library(rmapshaper)
# 
# #—— 1. Get & simplify all Canary Islands polygons ——#
# spain_states <- ne_states(country = "Spain", returnclass = "sf")
# can_all <- spain_states %>%
#   filter(name_en %in% c("Las Palmas", "Santa Cruz de Tenerife")) %>%
#   st_union() %>%
#   st_cast("MULTIPOLYGON") %>%
#   st_cast("POLYGON") %>%
#   st_as_sf()
# 
# # simplify for crisp plotting
# can_simp <- ms_simplify(can_all, keep = 0.1)
# 
# #—— 2. Keep only the 7 biggest polygons ——#
# # project to an equal‐area CRS so area is meaningful
# can_eq  <- st_transform(can_simp, 6933)
# areas   <- st_area(can_eq)
# idx7    <- order(areas, decreasing = TRUE)[1:7]
# can_main <- can_simp[idx7, ]
# 
# #—— 3. Auto‐name by matching centroids to known coords ——#
# # known island centroids (lon,lat)
# lookup <- tibble(
#   name = c("Tenerife","Gran Canaria","Lanzarote",
#            "Fuerteventura","La Palma","La Gomera","El Hierro"),
#   lon  = c(-16.6435,  -15.4314,  -13.5910,
#            -14.0168,  -17.8794,  -17.2586,  -18.0206),
#   lat  = c( 28.2916,   27.9959,   29.0469,
#             28.3587,   28.6833,   28.1084,   27.7264)
# )
# 
# # compute centroids of your 7 polygons
# cent <- st_centroid(can_main)
# coords <- st_coordinates(cent)
# 
# # for each centroid, find the nearest lookup point
# assigned <- sapply(1:nrow(coords), function(i) {
#   d <- sqrt((lookup$lon - coords[i,1])^2 +
#               (lookup$lat - coords[i,2])^2)
#   lookup$name[which.min(d)]
# })
# 
# can_main <- can_main %>% mutate(name = assigned)
# 
# #—— 4. Your island‐level data & join ——#
# island_data <- tibble(
#   name     = lookup$name,
#   jaccard  = c(0.34, 0.28, 0.41, 0.22, 0.37, 0.19, 0.25),
#   net_size = c(120,  95,   70,   50,   40,   35,   30)
# )
# 
# can_final <- left_join(can_main, island_data, by = "name")
# 
# # compute centroids again for labels
# cent_final <- st_centroid(can_final)
# 
# #—— 5. Get Western Sahara for context ——#
# ws <- ne_countries(scale = "large",
#                    country = "Western Sahara",
#                    returnclass = "sf")
# 
# #—— 6. Plot with pastel palette, no borders, and labels ——#
# ggplot() +
#   geom_sf(data = ws, fill = "grey90", color = NA) +
#   geom_sf(data = can_final,
#           aes(fill = jaccard),
#           color = NA) +
#   geom_sf_text(data = cent_final,
#                aes(label = name),
#                size = 3) +
#   scale_fill_gradient(
#     name = "Jaccard\nindex",
#     low  = "lightcyan",
#     high = "thistle"
#   ) +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text  = element_blank(),
#     axis.ticks = element_blank()
#   )
# 
# # 0. Install/load packages
# # install.packages(c("sf","rnaturalearth","rnaturalearthhires",
# #                    "dplyr","tibble","ggplot2","rmapshaper","ggrepel"))
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthhires)
# library(dplyr)
# library(tibble)
# library(ggplot2)
# library(rmapshaper)
# library(ggrepel)
# 
# # 1. Get & simplify Canary Islands
# spain_states <- ne_states(country = "Spain", returnclass = "sf")
# can_all <- spain_states %>%
#   filter(name_en %in% c("Las Palmas", "Santa Cruz de Tenerife")) %>%
#   st_union() %>%
#   st_cast("MULTIPOLYGON") %>%
#   st_cast("POLYGON") %>%
#   st_as_sf()
# can_simp <- ms_simplify(can_all, keep = 0.1)
# 
# # 2. Keep 7 largest polygons
# can_eq <- st_transform(can_simp, 6933)
# areas <- st_area(can_eq)
# idx7 <- order(areas, decreasing = TRUE)[1:7]
# can_main <- can_simp[idx7, ]
# 
# # 3. Auto-name by matching centroids to known coords
# lookup <- tibble(
#   name = c("Tenerife","Gran Canaria","Lanzarote",
#            "Fuerteventura","La Palma","La Gomera","El Hierro"),
#   lon  = c(-16.6435, -15.4314, -13.5910,
#            -14.0168, -17.8794, -17.2586, -18.0206),
#   lat  = c(28.2916, 27.9959, 29.0469,
#            28.3587, 28.6833, 28.1084, 27.7264)
# )
# cent <- st_centroid(can_main)
# coords <- st_coordinates(cent)
# assigned <- sapply(1:nrow(coords), function(i) {
#   d <- sqrt((lookup$lon - coords[i,1])^2 + (lookup$lat - coords[i,2])^2)
#   lookup$name[which.min(d)]
# })
# can_main <- can_main %>% mutate(name = assigned)
# 
# # 4. Join metrics
# tibble::tibble(
#   name     = lookup$name,
#   jaccard  = c(0.34, 0.28, 0.41, 0.22, 0.37, 0.19, 0.25),
#   net_size = c(633,   702,   570,   600,   650,   304,   350)
# ) -> island_data
# can_final <- left_join(can_main, island_data, by = "name")
# 
# # 5. Create label coords for repel
# cent_final <- st_centroid(can_final)
# cent_coords <- cent_final %>%
#   st_coordinates() %>%
#   as.data.frame() %>%
#   bind_cols(name = cent_final$name)
# 
# # 6. Plot without Western Sahara and with repelled labels
# p <- ggplot() +
#   geom_sf(data = can_final, aes(fill = net_size), color = NA) +
#   geom_text_repel(data = cent_coords, aes(X, Y, label = name),
#                   size = 3) +
#   scale_fill_gradient(
#     name = "Network\nsize",
#     low  = "lightsteelblue",
#     high = "thistle"
#   ) +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text  = element_blank(),
#     axis.ticks = element_blank()
#   )
# 
# print(p)
# 
# # 0. (Re-)load packages
# # install.packages(c("sf","rnaturalearth","rnaturalearthhires",
# #                    "dplyr","tibble","ggplot2","rmapshaper"))
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthhires)
# library(dplyr)
# library(tibble)
# library(ggplot2)
# library(rmapshaper)
# 
# # 0. (Re-)load packages
# # install.packages(c("sf","rnaturalearth","rnaturalearthhires",
# #                    "dplyr","tibble","ggplot2","rmapshaper"))
# library(sf)
# library(rnaturalearth); library(rnaturalearthhires)
# library(dplyr);    library(tibble)
# library(ggplot2);  library(rmapshaper)
# 
# # 1. Fetch & simplify Canary Islands (keep = 0.5 for more detail)
# spain_states <- ne_states(country = "Spain", returnclass = "sf")
# can_all <- spain_states %>%
#   filter(name_en %in% c("Las Palmas", "Santa Cruz de Tenerife")) %>%
#   st_union() %>%
#   st_cast("POLYGON") %>%
#   st_as_sf()
# can_simp <- ms_simplify(can_all, keep = 0.5)
# 
# # 2. Keep the 7 largest polygons
# can_eq   <- st_transform(can_simp, 6933)
# idx7     <- order(st_area(can_eq), decreasing = TRUE)[1:7]
# can_main <- can_simp[idx7, ]
# 
# # 3. Auto-name each polygon by nearest known centroid
# lookup <- tibble(
#   name = c("Tenerife","Gran Canaria","Lanzarote",
#            "Fuerteventura","La Palma","La Gomera","El Hierro"),
#   lon  = c(-16.6435, -15.4314, -13.5910,
#            -14.0168, -17.8794, -17.2586, -18.0206),
#   lat  = c( 28.2916,  27.9959,  29.0469,
#             28.3587,   28.6833,  28.1084,  27.7264)
# )
# cent_all <- st_centroid(can_main) %>% st_coordinates()
# assigned <- apply(cent_all, 1, function(pt) {
#   d <- sqrt((lookup$lon - pt[1])^2 + (lookup$lat - pt[2])^2)
#   lookup$name[which.min(d)]
# })
# can_main <- can_main %>% mutate(name = assigned)
# 
# # 4. Prepare non-Tenerife islands (La Palma & Lanzarote = no data)
# base_data <- tibble(
#   name     = c("Gran Canaria","Fuerteventura","La Gomera","El Hierro"),
#   net_size = c(702,600,304,350)
# )
# can_base <- can_main %>%
#   filter(name != "Tenerife") %>%
#   left_join(base_data, by = "name")
# 
# # 5. Split Tenerife with manual rectangles
# ten        <- can_main %>% filter(name == "Tenerife")
# lat_TB     <- 28.3531
# lat_F      <- 28.22228
# lat_split  <- mean(c(lat_TB, lat_F))
# bb         <- st_bbox(ten)
# 
# # build north rect polygon
# north_rect <- st_sfc(
#   st_polygon(list(matrix(
#     c(bb["xmin"], lat_split,
#       bb["xmax"], lat_split,
#       bb["xmax"], bb["ymax"],
#       bb["xmin"], bb["ymax"],
#       bb["xmin"], lat_split),
#     ncol = 2, byrow = TRUE
#   ))),
#   crs = st_crs(ten)
# )
# # build south rect polygon
# south_rect <- st_sfc(
#   st_polygon(list(matrix(
#     c(bb["xmin"], bb["ymin"],
#       bb["xmax"], bb["ymin"],
#       bb["xmax"], lat_split,
#       bb["xmin"], lat_split,
#       bb["xmin"], bb["ymin"]),
#     ncol = 2, byrow = TRUE
#   ))),
#   crs = st_crs(ten)
# )
# 
# # intersect to get the two halves
# net_TB   <- 420   # replace with your true Teno Bajo network size
# net_Fas  <- 375   # replace with your true Fasnia network size
# north_poly <- st_intersection(ten, north_rect) %>%
#   mutate(name = "Teno Bajo",  net_size = net_TB)
# south_poly <- st_intersection(ten, south_rect) %>%
#   mutate(name = "Fasnia",     net_size = net_Fas)
# 
# # 6. Combine all island pieces
# plot_isl <- bind_rows(can_base, north_poly, south_poly)
# 
# # 7. Add NW Sahara circle (area = Gran Canaria’s area)
# saha_net  <- base_data$net_size[base_data$name == "Gran Canaria"]
# ea_crs    <- 6933
# gr_ea     <- plot_isl %>% filter(name == "Gran Canaria") %>%
#   st_transform(ea_crs)
# r         <- sqrt(as.numeric(st_area(gr_ea)) / pi)
# center_pt <- st_sfc(st_point(c(-14.0, 26.5)), crs = 4326)
# circle_ea <- st_transform(center_pt, ea_crs) %>% st_buffer(r)
# circle_geo<- circle_ea %>% st_transform(4326) %>% st_as_sf() %>%
#   mutate(name = "NW Sahara", net_size = saha_net)
# 
# plot_all  <- bind_rows(plot_isl, circle_geo)
# 
# # 8. Compute external label positions
# centroids <- st_centroid(plot_all)
# cent_df   <- centroids %>% st_coordinates() %>% as.data.frame() %>%
#   bind_cols(name = plot_all$name)
# map_c     <- st_centroid(st_union(plot_all)) %>% st_coordinates()
# cent_df   <- cent_df %>%
#   mutate(dx = X - map_c[1], dy = Y - map_c[2],
#          dist = sqrt(dx^2 + dy^2),
#          label_x = X + 0.1 * dx / dist,
#          label_y = Y + 0.1 * dy / dist)
# 
# # 9. Final plot
# ggplot() +
#   # gray outlines for no‐data islands
#   geom_sf(data = filter(plot_all, is.na(net_size)),
#           fill = NA, color = "gray80", size = 0.4) +
#   # filled for islands + circle with data
#   geom_sf(data = filter(plot_all, !is.na(net_size)),
#           aes(fill = net_size), color = NA) +
#   # connector lines + external labels
#   geom_segment(data = cent_df,
#                aes(x = X, y = Y, xend = label_x, yend = label_y),
#                size = 0.3) +
#   geom_text(data = cent_df,
#             aes(x = label_x, y = label_y, label = name),
#             size = 3, hjust = 0) +
#   scale_fill_gradient(name = "Network\nsize",
#                       low  = "lightsteelblue",
#                       high = "thistle") +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text  = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank()    # no X/Y titles
#   )
# 
# # no label version
# 
# # 1. Fetch & simplify Canary Islands (keep = 0.5)
# spain_states <- ne_states(country = "Spain", returnclass = "sf")
# can_all <- spain_states %>%
#   filter(name_en %in% c("Las Palmas", "Santa Cruz de Tenerife")) %>%
#   st_union() %>% st_cast("POLYGON") %>% st_as_sf()
# can_simp <- ms_simplify(can_all, keep = 0.5)
# 
# # 2. Keep the 7 largest polygons & auto-name
# can_eq   <- st_transform(can_simp, 6933)
# idx7     <- order(st_area(can_eq), decreasing = TRUE)[1:7]
# can_main <- can_simp[idx7, ]
# lookup <- tibble(
#   name = c("Tenerife","Gran Canaria","Lanzarote",
#            "Fuerteventura","La Palma","La Gomera","El Hierro"),
#   lon  = c(-16.6435, -15.4314, -13.5910,
#            -14.0168, -17.8794, -17.2586, -18.0206),
#   lat  = c( 28.2916,  27.9959,  29.0469,
#             28.3587,   28.6833,  28.1084,  27.7264)
# )
# cent_all <- st_centroid(can_main) %>% st_coordinates()
# assigned <- apply(cent_all, 1, function(pt) {
#   d <- sqrt((lookup$lon - pt[1])^2 + (lookup$lat - pt[2])^2)
#   lookup$name[which.min(d)]
# })
# can_main <- can_main %>% mutate(name = assigned)
# 
# # 3. Base island data (La Palma & Lanzarote NA)
# base_data <- tibble(
#   name     = c("Gran Canaria","Fuerteventura","La Gomera","El Hierro"),
#   net_size = c(702,600,304,350)
# )
# can_base <- can_main %>%
#   filter(name != "Tenerife") %>%
#   left_join(base_data, by = "name")
# 
# # 4. Split Tenerife at midpoint latitude
# ten        <- can_main %>% filter(name == "Tenerife")
# lat_TB     <- 28.3531; lat_F <- 28.22228
# lat_split  <- mean(c(lat_TB, lat_F))
# bb         <- st_bbox(ten)
# north_rect <- st_sfc(st_polygon(list(matrix(
#   c(bb["xmin"], lat_split,
#     bb["xmax"], lat_split,
#     bb["xmax"], bb["ymax"],
#     bb["xmin"], bb["ymax"],
#     bb["xmin"], lat_split), ncol=2, byrow=TRUE
# ))), crs = st_crs(ten))
# south_rect <- st_sfc(st_polygon(list(matrix(
#   c(bb["xmin"], bb["ymin"],
#     bb["xmax"], bb["ymin"],
#     bb["xmax"], lat_split,
#     bb["xmin"], lat_split,
#     bb["xmin"], bb["ymin"]), ncol=2, byrow=TRUE
# ))), crs = st_crs(ten))
# net_TB  <- 420; net_Fas <- 375
# north_poly <- st_intersection(ten, north_rect) %>% mutate(name="Teno Bajo",  net_size=net_TB)
# south_poly <- st_intersection(ten, south_rect) %>% mutate(name="Fasnia",     net_size=net_Fas)
# 
# # 5. Combine islands
# plot_isl <- bind_rows(can_base, north_poly, south_poly)
# 
# # 6. Add NW Sahara circle
# saha_net  <- base_data$net_size[base_data$name=="Gran Canaria"]
# ea_crs    <- 6933
# gr_ea     <- plot_isl %>% filter(name=="Gran Canaria") %>% st_transform(ea_crs)
# r         <- sqrt(as.numeric(st_area(gr_ea))/pi)
# center_pt <- st_sfc(st_point(c(-14.0,26.5)), crs=4326)
# circle_geo<- st_transform(center_pt, ea_crs) %>% st_buffer(r) %>%
#   st_transform(4326) %>% st_as_sf() %>%
#   mutate(name="NW Sahara", net_size=saha_net)
# 
# plot_all <- bind_rows(plot_isl, circle_geo)
# 
# # 7. Fetch African coastline for context
# af_coast <- ne_coastline(scale="medium", returnclass="sf")
# 
# # 8. Final plot (no labels)
# ggplot() +
#   geom_sf(data = af_coast, color = "gray80", fill = NA) +
#   geom_sf(data = filter(plot_all, is.na(net_size)),
#           fill = NA, color = "gray80", size = 0.4) +
#   geom_sf(data = filter(plot_all, !is.na(net_size)),
#           aes(fill = net_size), color = NA) +
#   scale_fill_gradient(name="Network\nsize",
#                       low="lightsteelblue", high="thistle") +
#   coord_sf(xlim = c(-18, -12), ylim = c(26, 30)) +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text  = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank()
#   )

# 0. Load packages
library(sf)
library(rnaturalearth); library(rnaturalearthhires)
library(dplyr);    library(tibble)
library(ggplot2);  library(rmapshaper)

# 1. Fetch & simplify Canary Islands (keep = 0.5)
spain_states <- ne_states(country = "Spain", returnclass = "sf")
can_all <- spain_states %>%
  filter(name_en %in% c("Las Palmas", "Santa Cruz de Tenerife")) %>%
  st_union() %>%
  st_cast("POLYGON") %>%
  st_as_sf()
can_simp <- ms_simplify(can_all, keep = 0.5)

# 2. Keep the 7 largest polygons & auto‐name
can_eq   <- st_transform(can_simp, 6933)
idx7     <- order(st_area(can_eq), decreasing = TRUE)[1:7]
can_main <- can_simp[idx7, ]
lookup <- tibble(
  name = c("Tenerife","Gran Canaria","Lanzarote",
           "Fuerteventura","La Palma","La Gomera","El Hierro"),
  lon  = c(-16.6435, -15.4314, -13.5910,
           -14.0168, -17.8794, -17.2586, -18.0206),
  lat  = c( 28.2916,  27.9959,  29.0469,
            28.3587,   28.6833,  28.1084,  27.7264)
)
cent_all <- st_centroid(can_main) %>% st_coordinates()
assigned <- apply(cent_all, 1, function(pt) {
  d <- sqrt((lookup$lon - pt[1])^2 + (lookup$lat - pt[2])^2)
  lookup$name[which.min(d)]
})
can_main <- can_main %>% mutate(name = assigned)

# 3. Base data for non‐Tenerife islands (La Palma & Lanzarote = NA)
base_data <- tibble(
  name     = c("Gran Canaria","Fuerteventura","La Gomera","El Hierro"),
  net_size = c(702,570,304,350)
)
can_base <- can_main %>%
  filter(name != "Tenerife") %>%
  left_join(base_data, by = "name")

# 4. Split Tenerife at midpoint latitude
ten       <- filter(can_main, name == "Tenerife")
lat_TB    <- 28.3531; lat_F <- 28.22228
lat_split <- mean(c(lat_TB, lat_F))
bb        <- st_bbox(ten)
# north half
north_rect <- st_sfc(st_polygon(list(matrix(
  c(bb["xmin"], lat_split,
    bb["xmax"], lat_split,
    bb["xmax"], bb["ymax"],
    bb["xmin"], bb["ymax"],
    bb["xmin"], lat_split),
  ncol=2, byrow=TRUE
))), crs = st_crs(ten))
# south half
south_rect <- st_sfc(st_polygon(list(matrix(
  c(bb["xmin"], bb["ymin"],
    bb["xmax"], bb["ymin"],
    bb["xmax"], lat_split,
    bb["xmin"], lat_split,
    bb["xmin"], bb["ymin"]),
  ncol=2, byrow=TRUE
))), crs = st_crs(ten))
net_TB  <- 455; net_Fas <- 810
north_poly <- st_intersection(ten, north_rect) %>%
  mutate(name = "Teno Bajo",  net_size = net_TB)
south_poly <- st_intersection(ten, south_rect) %>%
  mutate(name = "Fasnia",     net_size = net_Fas)

# 5. Combine all island pieces + NW Sahara circle
plot_isl <- bind_rows(can_base, north_poly, south_poly)
#saha_net  <- base_data$net_size[base_data$name=="Gran Canaria"]
saha_net <- 350
# ea_crs    <- 6933
# gr_ea     <- filter(plot_isl, name=="Gran Canaria") %>% st_transform(ea_crs)
# r         <- sqrt(as.numeric(st_area(gr_ea))/pi)
# replace the old center_pt line with:
center_pt <- st_sfc(
  st_point(c(-14.42228, 26.1610)),
  crs = 4326
)
# compute radius from Gran Canaria’s area as before
ea_crs    <- 6933
gr_ea     <- plot_isl %>% filter(name=="Gran Canaria") %>% st_transform(ea_crs)
r         <- sqrt(as.numeric(st_area(gr_ea))/pi)

# make the correctly‐located circle
circle_geo <- center_pt %>%
  st_transform(ea_crs) %>%
  st_buffer(r) %>%
  st_transform(4326) %>%
  st_as_sf() %>%
  mutate(name = "NW Sahara", net_size = saha_net)

# circle_geo<- st_transform(center_pt, ea_crs) %>% st_buffer(r) %>%
#   st_transform(4326) %>% st_as_sf() %>%
#   mutate(name="NW Sahara", net_size=saha_net)
plot_all <- bind_rows(plot_isl, circle_geo)

# 6. Fetch Africa polygon for context
africa <- ne_countries(continent="Africa", scale="medium", returnclass="sf")

# 7. Plot: no labels, single outlines only for NA islands, no border on colored
ggplot() +
  # Africa outline for context
  geom_sf(data = africa, fill = NA, color = "gray80", linewidth = 0.7) +
  # islands with no data: gray outline only
  geom_sf(data = filter(plot_all, is.na(net_size)),
          fill = NA, color = "gray80", linewidth = 0.7) +
  # colored islands & circle: fill only, no border
  geom_sf(data = filter(plot_all, !is.na(net_size)),
          aes(fill = net_size), color = NA) +
  scale_fill_gradient(name="Network\nsize",
                      low="lightsteelblue", high="thistle") +
  coord_sf(xlim = c(-18, -12), ylim = c(26, 30)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

library(dplyr)
library(ggplot2)
library(sf)
library(grid)     # for arrow()

# --- 1. Compute island centroids (we’ll use these as arrow endpoints) ---
centroids <- st_centroid(plot_all)
cent_df   <- centroids %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  bind_cols(name = plot_all$name)

# --- 2. Define your edges with an f1 value for each link ---
#    (replace these with your real island pairs and f1s)
edge_df <- tibble(
  from = c("Tenerife","Gran Canaria","Fasnia"),
  to   = c("Gran Canaria","El Hierro","Teno Bajo"),
  f1   = c(0.12, 0.56, 0.33) # not real values
)

# --- 3. Join in the centroid coordinates for start/end points ---
edges_plot <- edge_df %>%
  left_join(cent_df, by = c("from" = "name"))  %>%
  rename(x   = X,  y   = Y) %>%
  left_join(cent_df, by = c("to"   = "name"))  %>%
  rename(xend = X, yend = Y)

# assume edges_plot has x,y,xend,yend already
short_frac <- 0.12

edges_short <- edges_plot %>%
  mutate(
    dx   = xend - x,
    dy   = yend - y,
    xs   = x   + short_frac * dx,      # new start
    ys   = y   + short_frac * dy,
    xe   = xend - short_frac * dx,     # new end
    ye   = yend - short_frac * dy
  )

# then in your ggplot replace geom_curve with:
geom_curve(
  data      = edges_short,
  aes(x = xs, y = ys, xend = xe, yend = ye, color = f1),
  curvature = 0.2,
  arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
  linewidth = 0.8
)

# --- 4. Add them to your existing ggplot as colored arrows ---
ggplot() +
  # (your background layers here: Africa coast, islands, circle…)  
  geom_sf(data = africa,         fill = NA, color = "gray80", linewidth = 0.8) +
  geom_sf(data = filter(plot_all, is.na(net_size)),
          fill = NA, color = "gray80", linewidth = 0.8) +
  geom_sf(data = filter(plot_all, !is.na(net_size)),
          aes(fill = net_size), color = NA) +
  
  # --- arrows, mapping color to f1 ---
  geom_curve(
    data      = edges_short,
    aes(x = xs, y = ys, xend = xe, yend = ye, color = f1),
    curvature = -0.5,
    arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.6
  ) +
  
  scale_color_gradient(
    name    = expression(f[1]),
    low     = "lightblue",
    high    = "salmon"
  ) +

  
  coord_sf(xlim = c(-18, -12), ylim = c(26, 30)) +
  scale_fill_gradient(
    name = "Network\nsize",
    low  = "lightsteelblue",
    high = "thistle"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

library(dplyr)

# 1. Recode the Tenerife names so they match the polygon names
metrics_edges <- result_summary_island %>% # that's from the analysis code
  mutate(
    from = recode(train_layer_name,
                  "Tenerife Teno"   = "Teno Bajo",
                  "Tenerife South"  = "Fasnia",
                  .default = train_layer_name),
    to   = recode(test_layer_name,
                  "Tenerife Teno"   = "Teno Bajo",
                  "Tenerife South"  = "Fasnia",
                  .default = test_layer_name)
  )

# 2. Filter out the diagonal comparisons
metrics_edges <- metrics_edges %>%
  filter(train_layer != test_layer)

# 3. Build edge_df
edge_df <- metrics_edges %>%
  select(from, to, f1 = f1_score)

# Inspect
print(edge_df)

edge_df <- edge_df %>%
  mutate(
    from = if_else(from %in% c("Western Sahara","Sahara"), "NW Sahara", from),
    to   = if_else(to   %in% c("Western Sahara","Sahara"), "NW Sahara", to)
  )

# real edges
edges_plot <- edge_df %>%
  left_join(cent_df, by = c("from" = "name")) %>% rename(x = X, y = Y) %>%
  left_join(cent_df, by = c("to"   = "name")) %>% rename(xend = X, yend = Y)

# --- 3. Join in the centroid coordinates for start/end points ---
edges_plot <- edge_df %>%
  left_join(cent_df, by = c("from" = "name"))  %>%
  rename(x   = X,  y   = Y) %>%
  left_join(cent_df, by = c("to"   = "name"))  %>%
  rename(xend = X, yend = Y)

# assume edges_plot has x,y,xend,yend already
short_frac <- 0.12

edges_short <- edges_plot %>%
  mutate(
    dx   = xend - x,
    dy   = yend - y,
    xs   = x   + short_frac * dx,      # new start
    ys   = y   + short_frac * dy,
    xe   = xend - short_frac * dx,     # new end
    ye   = yend - short_frac * dy
  )

# then in your ggplot replace geom_curve with:
geom_curve(
  data      = edges_short,
  aes(x = xs, y = ys, xend = xe, yend = ye, color = f1),
  curvature = 0.2,
  arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
  linewidth = 0.8
)

# --- 4. Add them to your existing ggplot as colored arrows ---
ggplot() +
  # (your background layers here: Africa coast, islands, circle…)  
  geom_sf(data = africa,         fill = NA, color = "gray80", linewidth = 0.8) +
  geom_sf(data = filter(plot_all, is.na(net_size)),
          fill = NA, color = "gray80", linewidth = 0.8) +
  geom_sf(data = filter(plot_all, !is.na(net_size)),
          aes(fill = net_size), color = NA) +
  
  # --- arrows, mapping color to f1 ---
  geom_curve(
    data      = edges_short,
    aes(x = xs, y = ys, xend = xe, yend = ye, color = f1),
    curvature = -0.5,
    arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.6
  ) +
  
  scale_color_gradient(
    name    = expression(f[1]),
    low     = "lightblue",
    high    = "salmon"
  ) +
  
  
  coord_sf(xlim = c(-18, -12), ylim = c(26, 30)) +
  scale_fill_gradient(
    name = "Network\nsize",
    low  = "lightsteelblue",
    high = "thistle"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

# selecting specific edges
library(dplyr)

# 1. pick only the three pairs you want
my_edges <- edge_df %>%
  filter(
    (from == "Fuerteventura" & to == "NW Sahara")    |
      (from == "Teno Bajo"   & to == "Fuerteventura") |
      (from == "Gran Canaria"    & to == "Fasnia")
  )

# 2. join on your centroid table (cent_df has columns name, X, Y)
my_edges_plot <- my_edges %>%
  left_join(cent_df, by = c("from" = "name")) %>%  rename(x   = X,  y   = Y) %>%
  left_join(cent_df, by = c("to"   = "name")) %>%  rename(xend= X,  yend= Y)

# 3. (optional) shorten them as before
short_frac <- 0.12
edges_short <- my_edges_plot %>%
  mutate(
    dx  = xend - x,
    dy  = yend - y,
    xs  = x   + short_frac * dx,
    ys  = y   + short_frac * dy,
    xe  = xend - short_frac * dx,
    ye  = yend - short_frac * dy
  )

edges_trimmed <- my_edges_plot %>%
  mutate(
    dx     = xend - x,
    dy     = yend - y,
    length = sqrt(dx^2 + dy^2),
    # custom trim: 0.3° for GC→Teno, 0.1° for the others
    trim   = case_when(
      from == "Gran Canaria" & to == "Tenerife Teno" ~ 0.4,
      TRUE                                         ~ 0.3
    ),
    # convert trim to fraction of each edge
    frac   = trim / length,
    # new start/end
    xs     = x   + frac * dx,
    ys     = y   + frac * dy,
    xe     = xend - frac * dx,
    ye     = yend - frac * dy
  )

# 4. plot just those arrows
p <- ggplot() +
  # (your background layers here: Africa coast, islands, circle…)  
  geom_sf(data = africa,         fill = NA, color = "gray80", linewidth = 0.8) +
  geom_sf(data = filter(plot_all, is.na(net_size)),
          fill = NA, color = "gray80", linewidth = 0.8) +
  geom_sf(data = filter(plot_all, !is.na(net_size)),
          aes(fill = net_size), color = NA) +
  
  # --- arrows, mapping color to f1 ---
  geom_curve(
    data      = edges_trimmed,
    aes(x = xs, y = ys, xend = xe, yend = ye, color = f1),
    curvature = -0.5,
    arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 1.2
  ) +
  
  scale_color_gradient(
    name    = expression(f[1]),
    low     = "lightblue",
    high    = "salmon"
  ) +
  
  
  coord_sf(xlim = c(-18, -12), ylim = c(26, 30)) +
  scale_fill_gradient(
    name = "Network\nsize",
    low  = "lightsteelblue",
    high = "thistle"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

# Base‐R PDF device
pdf(
  file   = "canary_network_map.pdf",
  width  = 10,    # inches
  height = 8,
  family = "Helvetica"   # or another installed font
)
print(p)      # draw your ggplot to the device
dev.off()     # close the file


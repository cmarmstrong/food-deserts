library(osmdata)
library(sf)
library(shiny)
library(units)

GETosm <- function(aabb, key, value='.') {
    aabb <- st_transform(aabb, 4326) # overpass requires 4326?
    query <- opq(st_bbox(aabb))
    query <- add_feature(query, key, value, value_exact=FALSE)
    osmdata_sf(query)
}

ui <- fluidPage(
    titlePanel('food-deserts')
    sidebarLayout(
        sidebarPanel(
            tableOutput(outputId='table')
        ),
        mainPanel(
            plotOutput(outputId='main',
                       click='click'))
    )
)

server <- function(input, output){}

shinyApp(ui=ui, server=server)

## GETosm(aabb, 'shop', c('bakery', 'butcher', 'cheese', 'deli', 'dairy', 'farm', 'greengrocer', 'froze_food', 'pasta', 'seafood', 'supermarket'))
query <- add_osm_feature(query, 'shop', 'bakery')
query <- add_osm_feature(query, 'shop', 'butcher')
query <- add_osm_feature(query, 'shop', 'cheese')
query <- add_osm_feature(query, 'shop', 'deli')
query <- add_osm_feature(query, 'shop', 'dairy')
query <- add_osm_feature(query, 'shop', 'farm')
query <- add_osm_feature(query, 'shop', 'greengrocer')
query <- add_osm_feature(query, 'shop', 'frozen_food')
query <- add_osm_feature(query, 'shop', 'pasta')
query <- add_osm_feature(query, 'shop', 'seafood')
query <- add_osm_feature(query, 'shop', 'supermarket')

foodOases <- osmdata_sf(query)
projPoints <- st_transform(foodOases $osm_points, 3857) # why 3857?
st_buffer(projPoints, ud_units $mi)
oases <- st_buffer(supermarkets, ud_units $mi)
## values(rUsa) <- as.numeric(log(values(rUsa)) > 3)
spAR <- as(AR, 'Spatial')
## rArGPW <- crop(rGPW, AR)

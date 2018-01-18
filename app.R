library(osmdata)
library(sf)
library(shiny)
library(units)

GETosm <- function(aabb, key, value='.') {
    aabb <- st_transform(aabb, 4326) # overpass requires 4326?
    query <- opq(st_bbox(aabb))
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
    # query <- add_feature(query, key, value, value_exact=FALSE)
    osmdata_sf(query)
}

rhumbDist2Degrees <- function(m, lat) {
    ## lat: latitude in decimal degrees
    ## m: distance in meters
    ## returns meters in 1 degree of longitude at latitude
    cos(lat) * 111321 # 111321 meters/degree at equator
}
## instead: rhumbDist2Degrees <- mymeters / distRhumb(pnt, pnt+c(1, 0))

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

server <- function(input, output) {
    rV <- reactiveValues(main=1:nstates, id=NA, condtions=NULL, click=NULL)
    ## click
    observe(rV $click <- input $click)
    observe({
        if(!is.null(rV $click)) {
            pnt <- st_point(as.numeric(c(rV $click $x, rV $click $y)))
            if(length(rV $main)==1) {
                pnt <- st_transform(st_sfc(pnt, crs=st_crs(usAdm1)), 4326) # OSM queries in 4326
                rV $id <- paste(st_coordinates(pnt)[2:1], collapse=',') # OSM queries in latlon
            } else {
                rV $main <- st_intersects(pnt, usAdm1)[[1]]
                rV $click <- NULL # cannot set input $click
            }
        }
    })
    ## double click
    observe({
        dblclick <- input $dblclick
        if(!is.null(dblclick)) { # reset
            rV $main <- 1:nstates
            rV $id <- NA
            rV $click <- NULL
        }
    })
    ## main
    output $main <- renderPlot({
        pnt <- rV $click
        aabb <- st_make_grid(pnt, what='center')
        ## bbox <- st_bbox(coords)
        ## aabb <- aabb(bbox)
        osm <- GETosm(aabb)
        layout(matrix(1:2, nrow=1), widths=c(5, 1))
        if(nrow(osm $osm_lines)>0) { # if query is not empty
            highways <- with(osm $osm_lines, osm $osm_lines[highway %in% highways, ])
            plot(st_geometry(st_transform(highways, espg)),
                 xlim=bbox[c(1, 3)], ylim=bbox[c(2, 4)], col='grey50', axes=TRUE)
        }
        ## plot urban
        ## plot roads?
        ## plot food sources
        ## plot food deserts
        ## each cell has a cost: rural 1 urban 10
    }, width=plotWidth, height=plotHeight)
}

shinyApp(ui=ui, server=server)

## GETosm(aabb, 'shop', c('bakery', 'butcher', 'cheese', 'deli', 'dairy', 'farm', 'greengrocer', 'froze_food', 'pasta', 'seafood', 'supermarket'))

## instead: rhumbDist2Degrees <- mymeters / distRhumb(pnt, pnt+c(1, 0))
pnt <- c(40.1165037, -88.2417297)
pnt <- st_sfc(st_point(pnt), crs=4326)
westPnt <- pnt + c(0, -1)
eastPnt <- pnt + c(0, 1)
northPnt <- pnt + c(1, 0)
southPnt <- pnt + c(-1, 0)

spPnt <- as(pnt, 'Spatial')
spPntN <- as(westPnt, 'Spatial')
spPntE <- as(eastPnt, 'Spatial')
spPntW <- as(northPnt, 'Spatial')
spPntS <- as(southPnt, 'Spatial')
rhumbDist2Degrees <- 10000 / distRhumb(spPnt, spPntE) # it works!: 10km to degrees at longlat

## aabb <- st_make_grid(pnt, what='centers')
foodOases <- GETosm(aabb)
projPoints <- st_transform(foodOases $osm_points, 3857) # why 3857?
st_buffer(projPoints, ud_units $mi)
oases <- st_buffer(supermarkets, ud_units $mi)
## values(rUsa) <- as.numeric(log(values(rUsa)) > 3)
spAR <- as(AR, 'Spatial')
## rArGPW <- crop(rGPW, AR)

library(geosphere)
library(osmdata)
library(RColorBrewer)
library(rnaturalearth)
library(sf)
library(shiny)
library(units)

## constants
nstates <- 51 # for selecting all states+DC from naturalearth
plotWidth <- 960
plotHeight <- 600
tol <- .Machine$double.eps^0.5
## proj4 string
albersEqualAreaConic <- '+proj=aea +lat_1=27.5 +lat_2=35 +lat_0=18 +lon_0=-100 +x_0=1500000 +y_0=6000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
webMercator <- '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

## geometries
statesUS <- ne_states(country='united states of america', returnclass='sf')
statesUS <- st_transform(statesUS, 3083)
statesUS <- with(statesUS, statesUS[order(name), ]) # 4-color assignment is in alpha order
statesUS $color <- factor(c(2, 1, 4, 3, 1, 2, 3, 3, 4, 1, 4, 4, 2, 4, 1, 3, 1, 2, 1, 1, 4, 2, 4, 4, 4, 2, 4, 4, 3, 4, 1, 3, 4, 2, 2, 3, 4, 4, 2, 1, 1, 1, 1, 2, 1, 1, 3, 1, 1, 1, 3))

urbanUS <- load('data/urbanUS.rda')
nullOSM <- load('data/nullOSM.rda')

## statesUS colors
col <- brewer.pal(4, 'Pastel1')[as.numeric(statesUS $color)]

## functions
dist2Degrees <- function(d, p1, p2) {
    if(abs(dist(rbind(p1 @coords, p2 @coords)) - 1) > tol) warning('euclidian dist != 1')
    d / distGeo(p1, p2)
}

GETosm <- function(aabb, key, value='.') {
    aabb <- st_transform(aabb, 4326) # overpass requires 4326?
    query <- opq(st_bbox(aabb))
    query <- add_osm_feature(query, key, value, value_exact=FALSE)
    osmdata_sf(query)
}

## main
ui <- fluidPage(
    titlePanel('food-deserts'),
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
    rV <- reactiveValues(click=NULL, pnt=NA)
    ## click
    observe(rV $click <- input $click)
    observe({
        if(!is.null(rV $click)) {
            pnt <- st_point(as.numeric(c(rV $click $x, rV $click $y)))
            rV $pnt <- st_transform(st_sfc(pnt, crs=st_crs(statesUS)), 4326) # WU queries in 4326
        } else rV $click <- NULL # us rV b/c cannot set input $click
    })
    ## double click
    observe({
        dblclick <- input $dblclick
        if(!is.null(dblclick)) { # reset
            rV $click <- NULL
            rV $pnt <- NA
        }
    })
    ## reactives
    GETosm <- reactive({
        if(!is.na(rV $pnt)) return nullOSM
        pnt <- rV $pnt
        pntE <- pnt + c(0, 1)
        pntN <- pnt + c(1, 0)
        spPntE <- as(pntE, 'Spatial')
        spPntN <- as(pntN, 'Spatial')
        dE <- dist2Degrees(7e3, spPnt, spPntE)
        dN <- dist2Degrees(7e3, spPnt, spPntN)
        pntNE <- pnt + c(dN, dE)
        pntNW <- pnt + c(dN, -dE)
        pntSE <- pnt + c(-dN, dE)
        pntSW <- pnt + c(-dN, -dE)
        ## make bbox for overpass
        aabb <- st_make_grid(st_sfc(c(pntNW, pntNE, pntSE, pntSW)), n=1)
        st_crs(aabb) <- 4326
        ## query
        GETosm(aabb, 'shop', 'supermarket')
    })
    ## table
    output $table <- renderTable({
        browser()
        osm <- GETosm()
        with(osm $osm_points,
             as.matrix(c(supermarkets=0,
                         bakery=0,
                         'long lat'="",
                         ))
             )
    }, rownames=TRUE, colnames=FALSE)
    ## main
    output $main <- renderPlot({
        osm <- GETosm()
        if(nrow(osm $osm_points)>0) { # if query is not empty
            osm <- st_transform(osm $osm_points, 3083) # buffer in projection
            food <- st_buffer(osm, ud_units $mi)
            food <- st_union(food)
            urbanFood <- crop(urbanUS, as(food, 'Spatial'))
            ## plot
            plot(urbanFood)
            plot(st_geometry(food), add=TRUE)
        } else {
            bbox <- st_bbox(statesUS)
            plot(statesUS[, 'color'], xlim=bbox[c(1, 3)], ylim=bbox[c(2, 4)],
                 col=col, main=NA, border=NA, graticule=st_crs(3083), axes=TRUE, key.pos=NULL, lwd.tick=0)
        }
    }, width=plotWidth, height=plotHeight)
}

shinyApp(ui=ui, server=server)

## query <- add_osm_feature(query, 'shop', 'bakery')
## query <- add_osm_feature(query, 'shop', 'butcher')
## query <- add_osm_feature(query, 'shop', 'cheese')
## query <- add_osm_feature(query, 'shop', 'deli')
## query <- add_osm_feature(query, 'shop', 'dairy')
## query <- add_osm_feature(query, 'shop', 'farm')
## query <- add_osm_feature(query, 'shop', 'greengrocer')
## query <- add_osm_feature(query, 'shop', 'frozen_food')
## query <- add_osm_feature(query, 'shop', 'pasta')
## query <- add_osm_feature(query, 'shop', 'seafood')
## query <- add_osm_feature(query, 'shop', 'supermarket')

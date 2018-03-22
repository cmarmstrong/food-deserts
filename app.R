library(buffy)
library(geosphere)
library(osmdata)
library(RColorBrewer)
library(raster)
library(rnaturalearth)
library(sf)
library(shiny)
library(units)

## constants
buffBbox1 <- 7e4
buffBbox2 <- 7e3
nstates <- 51 # for selecting all states+DC from naturalearth
plotWidth <- 960
plotHeight <- 600
tol <- .Machine$double.eps^0.5
## proj4 string
albersEqualAreaConic <- '+proj=aea +lat_1=27.5 +lat_2=35 +lat_0=18 +lon_0=-100 +x_0=1500000 +y_0=6000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
webMercator <- '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
highways <- c('motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'unclassified', 'residential')

## geometries
statesUS <- ne_states(country='united states of america', returnclass='sf')
statesUS <- st_transform(statesUS, 3083)
statesUS <- with(statesUS, statesUS[order(name), ]) # 4-color assignment is in alpha order
statesUS $color <- factor(c(2, 1, 4, 3, 1, 2, 3, 3, 4, 1, 4, 4, 2, 4, 1, 3, 1, 2, 1, 1, 4, 2, 4, 4, 4, 2, 4, 4, 3, 4, 1, 3, 4, 2, 2, 3, 4, 4, 2, 1, 1, 1, 1, 2, 1, 1, 3, 1, 1, 1, 3))

load('data/citiesUS.rda')
load('data/urbanUS.rda')
load('data/nullOSM.rda')

## statesUS colors
col <- brewer.pal(4, 'Pastel1')[as.numeric(statesUS $color)]

## functions
dist2Degrees <- function(d, p1, p2) { # d:= dist, p* := 2 points in long lat
    if(abs(dist(rbind(p1 @coords, p2 @coords)) - 1) > tol) warning('euclidian dist != 1')
    d / distGeo(p1, p2)
}

queryOSM <- function(aabb, key, value='.') {
    aabb <- st_transform(aabb, 4326) # overpass requires 4326?
    query <- opq(st_bbox(aabb))
    query <- add_osm_feature(query, key, value, value_exact=FALSE)
    osmdata_sf(query)
}

bufferSquare <- function(pnt, d) { # d:= half length of the square
    pntE <- pnt + c(0, 1)
    pntN <- pnt + c(1, 0)
    spPnt <- as(pnt, 'Spatial')
    spPntE <- as(pntE, 'Spatial')
    spPntN <- as(pntN, 'Spatial')
    dE <- dist2Degrees(d, spPnt, spPntE)
    dN <- dist2Degrees(d, spPnt, spPntN)
    pntNE <- pnt + c(dN, dE)
    pntNW <- pnt + c(dN, -dE)
    pntSE <- pnt + c(-dN, dE)
    pntSW <- pnt + c(-dN, -dE)
    st_make_grid(st_sfc(c(pntNW, pntNE, pntSE, pntSW)), n=1)
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
                       click='click',
                       dblclick='dblclick',
                       hover=hoverOpts(id='hover', delay=50, delayType='throttle'),
                       brush=brushOpts(id='brush', resetOnNew=TRUE))
        )
    )
)

server <- function(input, output) {
    rV <- reactiveValues(bbox1=NULL, click=NULL, pnt=NA)
    ## observers
    ## click
    observe(rV $click <- input $click)
    observe({
        if(!is.null(rV $click)) {
            pnt <- st_point(as.numeric(c(rV $click $x, rV $click $y)))
            rV $pnt <- st_transform(st_sfc(pnt, crs=st_crs(statesUS)), 4326) # WU queries in 4326
            if(is.null(rV $bbox1)) {
                aabb <- bufferSquare(rV $pnt, buffBbox1)
                st_crs(aabb) <- 4326
                rV $bbox1 <- st_transform(aabb, 3083)
            } else if(is.null(rV $bbox2)) {
                aabb <- bufferSquare(rV $pnt, buffBbox2)
                st_crs(aabb) <- 4326
                rV $bbox2 <- aabb
            }
            rV $click <- NULL # use rV b/c cannot set input $click
        }
    })
    ## double click
    observe({
        dblclick <- input $dblclick
        if(!is.null(dblclick)) { # reset
            rV $click <- NULL
            rV $pnt <- NA
            if(is.null(rV $bbox2)) rV $bbox1 <- NULL
            rV $bbox2 <- NULL
        }
    })
    ## reactives
    ## osm
    getFood <- reactive({
        if(is.na(rV $pnt)) osm <- nullOSM
        else {
            aabb <- bufferSquare(rV $pnt, buffBbox2)
            st_crs(aabb) <- 4326
            osm <- queryOSM(aabb, 'shop', 'supermarket')
        }
        osm
    })
    getStreets <- reactive({
        if(is.na(rV $pnt)) osm <- nullOSM
        else {
            aabb <- bufferSquare(rV $pnt, buffBbox1)
            st_crs(aabb) <- 4326
            osm <- queryOSM(aabb, 'highway', highways[1])
        }
        osm
    })
    ## outputs
    ## table
    output $table <- renderTable({
        osm <- getFood()
        if(nrow(osm $osm_points)>0) {
            with(osm,
                 as.matrix(c(supermarkets=nrow(osm_points),
                             'long lat'=st_as_text(rV $pnt)
                             )))
        } else as.matrix(c(supermarkets=NA))
    }, rownames=TRUE, colnames=FALSE)
    ## main
    output $main <- renderPlot({
        osmFood <- getFood()
        osmStreets <- getStreets()
        if(!is.null(rV $bbox2) && nrow(osmFood $osm_points)>0) { # if bbox2 & query not empty
            browser()
            osm <- osmFood $osm_points
            osm <- sf::st_transform(osm, 3083) # buffer in projection
            buffOsm <- sf::st_buffer(osm, 10*units::ud_units $mi) # max buffer for cropping
            urbanLocal <- crop(urbanUS, as(buffOsm, 'Spatial'))
            urbanPolys <- rasterToPolygons(urbanLocal, digits=6, dissolve=TRUE)
            urbanPolys <- st_as_sf(urbanPolys)
            st_crs(urbanPolys) <- 3083
            urbanPolys $s <- 10 # HERES THE PROBLEM, urbanPolys now an sfc object and assignment not working
            food <- buffy::surfBuff(st_geometry(osm), urbanPolys, 10*units::ud_units $mi)
            food <- st_union(food)
            ## plot
            ## bboxFood <- sf::st_bbox(food)
            ## plot(0, main='food deserts', xlab='', ylab='', # initialize plot window
            ##      xlim=bboxFood[c(1,3)], ylim=bboxFood[c(2,4)])
            plot(urbanLocal, main='food deserts', col='orange', legend=FALSE, add=TRUE)
            plot(st_geometry(st_as_sf(urbanPolys)), add=TRUE)
            plot(st_geometry(food), add=TRUE)
            plot(st_geometry(citiesUS), add=TRUE)
            plot(st_geometry(st_transform(osmStreets $osm_lines, 3083)), add=TRUE)
        }
        else if(!is.null(rV $bbox1)) { # need nrow(osmStreets $osm_lines)>0 ?
            urbanBox <- crop(urbanUS, as(rV $bbox1, 'Spatial'))
            plot(urbanBox, main='urban areas', col='orange', legend=FALSE)
            plot(st_geometry(statesUS), add=TRUE)
            plot(st_geometry(citiesUS[citiesUS $scale50, ]), add=TRUE)
            plot(st_geometry(st_transform(osmStreets $osm_lines, 3083)), add=TRUE)
        } else {
            bbox <- st_bbox(statesUS)
            plot(statesUS[, 'color'], xlim=bbox[c(1, 3)], ylim=bbox[c(2, 4)], col=col, main=NA,
                 border=NA, graticule=st_crs(3083), axes=TRUE, key.pos=NULL, lwd.tick=0)
            plot(st_geometry(citiesUS[citiesUS $scale110, ]), add=TRUE)
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

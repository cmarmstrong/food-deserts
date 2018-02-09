library(geosphere)
library(osmdata)
library(RColorBrewer)
library(raster)
library(rnaturalearth)
library(sf)
library(sp)
library(units)


buffBbox1 <- 1e5
buffBbox2 <- 7e3
tol <- .Machine$double.eps^0.5

load('data/citiesUS.rda')
load('data/nullOSM.rda')
load('data/urbanUS.rda')

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

getFood <- function(pnt) {
    if(is.na(pnt)) osm <- nullOSM
    else {
        aabb <- bufferSquare(pnt, buffBbox2)
        st_crs(aabb) <- 4326
        osm <- queryOSM(aabb, 'shop', 'supermarket')
    }
    osm
}

getStreets <- function(pnt) {
    if(is.na(pnt)) osm <- nullOSM
    else {
        aabb <- bufferSquare(pnt, buffBbox2)
        st_crs(aabb) <- 4326
        osm <- queryOSM(aabb, 'highway', 'motorway')
    }
    osm
}

## snap <- function(sfd, threshold) {
##     browser()
##     d <- st_distance(sfd)
##     hc <- hclust(as.dist(d>threshold), method='single')
##     groups <- cutree(hc, h=0.5)
##     sfdSnapped <- st_sf(geom=do.call(c, lapply(1:max(groups), function(g) {
##         st_union(sfd[groups==g, ])
##     })))
##     sfdSnapped $group <- 1:nrow(sfdSnapped)
##     sfdSnapped
## }

bufferAnalysis <- function() {
    browser()
    pnt <- st_sfc(st_point(as.numeric(c(-92.44744828628, 34.566107548536)))) # ~ little rock, AR
    osmFood <- getFood(pnt)
    osmStreets <- getStreets(pnt)
    osm <- st_transform(osmFood $osm_points, 3083) # buffer in projection
    food <- st_buffer(osm, 10*ud_units $mi)
    coordsFood <- st_coordinates(food)
    lCoords <- by(coordsFood, coordsFood[, 'L2'], apply, 1, function(M) {
        rbind(M[1:2], st_coordinates(osm[M[4], ]))
    })
    ## mls
    ## gMls <- lapply(lapply(lCoords, function(mls) {
    ##     lapply(split(t(mls), seq(NCOL(mls))), matrix, nrow=2)
    ## }), st_multilinestring)
    ## sfMls <- do.call(st_sfc, gMls)
    ## ls
    lM <- lapply(lCoords, function(coords) {
        lapply(split(t(coords), seq(NCOL(coords))), matrix, nrow=2)
    })
    lSf <- lapply(1:length(lM), function(i) {
        m <- lM[[i]]
        sfcLs <- do.call(st_sfc, lapply(m, st_linestring))
        st_sf(geometry=sfcLs, idLs=1:length(sfcLs), idFeature=i)
    })
    sfLs <- do.call(rbind, lSf)
    st_crs(sfLs) <- 3083
    
    food <- st_union(food) # union after buffers expanded for rural
    forCrop <- st_buffer(food, ud_units $mi) # plot slightly larger area
    urbanFood <- crop(urbanUS, as(forCrop, 'Spatial'))
    urbanPolys <- rasterToPolygons(urbanFood, digits=7, dissolve=TRUE) # truncate digits to snap polygons
    urbanPolys <- st_as_sf(urbanPolys)

    sfIntersections <- st_intersection(sfLs, urbanPolys)
    sfLs $urbanLengths <- st_length(st_cast(sfIntersections, 'MULTILINESTRING'))
    ## totalLength <- urbanLength + ruralLength/10
    ## xsLength <- totalLength - 16093.44
    ## urbanLength + ruralLength/10 = 1 mile = 16093.44 meteres
    ## cast to linestring and order by st_distance from center point
    ## order linestrings from furthest to closest from center
    ## for each linestring,
    ##   if length of linestring < xsLength
    ##     then xsLength <- xsLength - length of linestring
    ##     remove linestring from list
    ##   else
    ##     length of linestring <- length of linestring - xsLength
    ##     break

    ## plot
    plot(urbanFood, main='food deserts', col='orange', legend=FALSE)
    plot(st_geometry(urbanPolys), add=TRUE)
    plot(st_geometry(food), add=TRUE)
    plot(st_geometry(citiesUS), add=TRUE)
}

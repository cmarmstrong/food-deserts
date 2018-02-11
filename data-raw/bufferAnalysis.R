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
    
    ## food <- st_union(food) # union after buffers expanded for rural
    forCrop <- st_buffer(food, ud_units $mi) # plot slightly larger area
    urbanFood <- crop(urbanUS, as(forCrop, 'Spatial'))
    urbanPolys <- rasterToPolygons(urbanFood, digits=7, dissolve=TRUE) # digits to snap
    urbanPolys <- st_as_sf(urbanPolys)

    sfIn <- st_intersection(sfLs, urbanPolys)
    sfIn $urbanLength <- st_length(st_cast(sfIn, 'MULTILINESTRING'))
    sfIn $ruralLength <- set_units(10*ud_units $mi - sfIn $urbanLength, 'm')
    sfIn $totalLength <- with(sfIn, urbanLength + ruralLength/10)
    sfIn $xsLength <- set_units(ud_units $mi - sfIn $totalLength, 'm')
    ## urbanLength + ruralLength/10 = 1 mile = 16093.44 meteres
    newBuffer <- mapply(function(lstring, mls) {
        coordsLs <- st_coordinates(lstring)
        coordsLs[, 'L1'] <- c(0, 0)
        coordsLs <- cbind(coordsLs, L2=c(1, 1))
        coordsMls <- rbind(coordsLs[1, ], st_coordinates(mls), coordsLs[2, ])
        coordsMls <- coordsMls[, c('X', 'Y')]
        sfLs <- st_sf(urban=c(rep(c(0, 1), (nrow(coordsMls)-1)%/%2), 0),
                      geometry=st_sfc(lapply(2:nrow(coordsMls), function(i) {
                          st_linestring(rbind(coordsMls[i-1, ], coordsMls[i, ]))
                      })))
        st_crs(sfLs) <- 3083
        lenLs <- st_length(sfLs)
        csumLs <- cumsum(ifelse(sfLs $urban==1, lenLs, lenLs/10))
        xsLs <- with(ud_units, csumLs*m > mi)
        lsXs <- sfLs[xsLs, ][1, ]
        ## isUrban <- sfLs[xsLs, 'urban', drop=TRUE][1]
        ## isUrban <- lsXs $urban
        xsLen <- with(ud_units, csumLs[xsLs][1]*m - mi)
        startLen <- rev(csumLs[!xsLs])[1]
        okLen <- set_units(with(ud_units, mi - startLen*m), 'm')
        okLen <- ifelse(lsXs $urban==1, okLen, okLen*10)
        ## xsCoords <- st_coordinates(st_transform(sfLs[xsLs, ], 4326))[, c('X', 'Y')]
        xsCoords <- st_coordinates(st_transform(lsXs, 4326))[, c('X', 'Y')]
        b <- bearing(xsCoords)
        ## NOTE: newPnt may not be correct
        newPnt <- destPoint(xsCoords, b, okLen)[1, ]
        ## newLs <- st_linestring(rbind(xsCoords1[1, ], newPnt[1, ]))
        ## newSf <- st_sf(geometry=st_sfc(newLs), urban=isUrban)
        ## st_crs(newSf) <- 4326
        ## sfLs <- rbind(sfLs[!xsLs, ], st_transform(newSf, 3083))
        ## rev(st_coordinates(sfLs))[1]
        sfcPnt <- st_sfc(st_point(newPnt))
        st_crs(sfcPnt) <- 4326
        st_coordinates(st_transform(sfcPnt, 3083))
    }, split(sfLs, 1:nrow(sfLs)), split(sfIn, 1:nrow(sfIn)))
    newBuffer <- t(newBuffer)
    lNewBuffer <- lapply(split(newBuffer, rep(1:29, each=121)), matrix, ncol=2)
    sfcNewBuffers <- st_sfc(lapply(lNewBuffer, function(newBuffer) {
        st_polygon(list(rbind(newBuffer, newBuffer[1, ])))}))
    foodDeserts <- st_union(sfcNewBuffers)
    st_crs(foodDeserts) <- 3083

    ## plot
    plot(urbanFood, main='food deserts', col='orange', legend=FALSE)
    plot(st_geometry(urbanPolys), add=TRUE)
    plot(st_geometry(foodDeserts), add=TRUE)
    plot(st_geometry(food), add=TRUE)
    plot(st_geometry(citiesUS), add=TRUE)
}

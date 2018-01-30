buffBbox1 <- 1e5
buffBbox2 <- 7e3

load('data/nullOSM.rda')

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

bufferAnalysis <- function(osm, urbanUS) {
    pnt <- {}
    osm <- getFood(pnt)
    osm <- st_transform(osmFood $osm_points, 3083) # buffer in projection
    food <- st_buffer(osm, ud_units $mi)
    ## food <- st_union(food) # union after buffers expanded for rural
    forCrop <- st_buffer(food, ud_units $mi) # plot slightly larger area
    urbanFood <- crop(urbanUS, as(forCrop, 'Spatial'))
    urbanPolys <- rasterToPolygons(urbanFood, dissolve=TRUE)
    ## plot
    plot(urbanFood, main='food deserts', col='orange', legend=FALSE)
    plot(st_geometry(st_as_sf(urbanPolys)), add=TRUE)
    plot(st_geometry(food), add=TRUE)
    plot(st_geometry(citiesUS), add=TRUE)
    plot(st_geometry(st_transform(osmStreets $osm_lines, 3083)), add=TRUE)
}

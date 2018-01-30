library(rnaturalearth)
library(raster)
library(sf)

albersEqualAreaConic <- '+proj=aea +lat_1=27.5 +lat_2=35 +lat_0=18 +lon_0=-100 +x_0=1500000 +y_0=6000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

## bkGPW <- raster::brick('E:/data/GPW/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev10_2015_30_sec_tif/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev10_2015_30_sec.tif')
bkGPW <- raster::brick('E:/data/GPW/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev10_2015_30_sec_tif/gpw_popdensity_3083.tif')
names(bkGPW) <- 'density'

spUS <- ne_states(country='united states of america')
spUS <- spTransform(spUS, CRS(albersEqualAreaConic))
spAR <- spUS[spUS $fips=='US05', ]
bkUsGpw <- crop(bkGPW, spUS)
bkArGpw <- crop(bkGPW, spAR)

urbanAnalysis <- function(bk) {
    ## bk <- addLayer(bk, log(raster(bk, layer=1)) > 4)
    bk <- addLayer(bk, log(raster(bk, layer=1)) > 5)
    ## bk <- addLayer(bk, log(raster(bk, layer=1)) > 6)
    ## names(bk)[2:4] <- paste0('log_density_gt_', 4:6)
    names(bk)[2] <- 'log_density_gt_5'

    rClump <- clump(raster(bk, layer=2))
    bk <- addLayer(bk, rClump)
    ## omit clumps > 9 cells
    tooSmall <- freq(rClump)[, 'count'] < 9
    rClump[rClump %in% which(tooSmall)] <- NA
    rUrban <- !(is.na(rClump)) ## > 1 : 1, NA : 0
    ## values(rClump)[is.na(values(raster(bk, layer=1)))] <- NA # reset NA values to layer 1
    values(rUrban)[values(rUrban)==0] <- NA
    bk <- addLayer(bk, rUrban)
    names(bk)[3:4] <- c('clumpID', 'clump_n_gt_8')
    bk
}

bkUrban <- urbanAnalysis(bkUsGpw)
plot(bkUrban)
writeRaster(raster(bkUrban, layer=4), 'urbanUS.tif')

urbanUS <- readAll(raster('data-raw/urbanUS.tif'))
save(urbanUS, file='data/urbanUS.rda', compress='bzip2')

nullOSM <- st_transform(st_make_grid(
    ne_countries(scale=10, country='united states of america'), n=1), 3083)
save(nullOSM, file='data/nullOSM.rda')

citiesGLB <- ne_download(scale=10, type='populated_places', returnclass='sf')
citiesGLB50 <- ne_download(scale=50, type='populated_places', returnclass='sf')
citiesGLB110 <- ne_download(scale=110, type='populated_places', returnclass='sf')

citiesGLB $scale50 <- citiesGLB $GEONAMEID %in% citiesGLB50 $GEONAMEID
citiesGLB $scale110 <- citiesGLB $GEONAMEID %in% citiesGLB110 $GEONAMEID
citiesUS <- citiesGLB[citiesGLB $ADM0_A3=='USA', ]

citiesUS <- st_transform(citiesUS, 3083)
save(citiesUS, file='data/citiesUS.rda')

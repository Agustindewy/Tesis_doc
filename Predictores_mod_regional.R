
# Predictores scenopoéticos para los medelos regionales: descarga y preparación
# Predictores scenopoéticos para G. galeus y N. cepedianus
# Predictor biótico para N. cepedianus

library(raster)
library(rgdal)
library(sqldf)
library(maps)
library(ggplot2)
library(grec)
library(ncdf4)

setwd('DIRECTORIO DE TRABAJO')

# Proyección
crs <- CRS('+init=epsg:4326')

# Costa de la región (polígonos espaciales) - tomado de la base del GSHHG
coast <- readOGR(dsn = 'DESCARGAR Y LEER EL SHAPE', layer = 'GSHHS_f_L1_SouthAmerica')

# Estaciones
season <- c('summer', 'autumn', 'winter', 'spring')


#-------------------------------------- Extensión espacial ------------------------------------------

ext <- c(-71, -35, -60, -12) # (lon.min, lon.max, lat.min, lat.max)

# Raster de base
bb <- matrix(c(-70, -60, -35, -12), 2, 2) 
cs <- c(0.041, 0.041) 
cc <- bb[, 1] + (cs / 2) 
cd <- ceiling(diff(t(bb)) / cs) 
gri <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd) 
sp.gri <- SpatialGrid(gri, proj4string = crs) 
sp.ras <- raster(sp.gri) 


#-------------------------------------- Predictores estáticos ------------------------------------------------

# Distancia a la costa (km) - tomado de MARSPEC (30 arcseconds, ~1 km)
dc <- raster('DESCARGAR Y LEER EL RASTER/Dist_30s.tif')
dc <- crop(dc, ext) 
dc <- projectRaster(dc, sp.ras)
m <- matrix(1, ncol = 5, nrow = 5)
dc <- focal(dc, m, fun = mean, NAonly = T, pad = T, na.rm = T)
dc <- mask(dc, coast, inverse = T)

# Profundidad (m) - tomado de MARSPEC (30 arcseconds, ~1 km)
dep <- raster('DESCARGAR Y LEER EL RASTER/Bathy_30s.tif')
dep <- crop(dep, dc) 
dep <- projectRaster(dep, dc)
dep <- focal(dep, m, fun = mean, NAonly = T, pad = T, na.rm = T) 
dep <- mask(dep, dc)

static_pred <- stack(dc, dep)
writeRaster(static_pred, filename = 'static_pred.tif')


#-------------------------------------- Predictores dinámicos ------------------------------------------------

#-------------------------------------- Temperatura de superficie ----------------------------------------

# ºC, 01-2003 to 12-2020 mensual, ~4 km)

# Descargar en 2 periodos porque es muy pesado (1 = 2003-2011, 2 = 2012-2020)
sst <- read.csv('DESCARGAR Y LEER EL CSV/sst_1.csv', skip = 2, header = F)
sst <- sst[which(!is.na(sst$V4)), ] 
sst$Year <- substr(sst$V1, 1, 4)
sst$Month <- substr(sst$V1, 6, 7)
sst$Month <- as.numeric(sst$Month)
sst$Season2 <- ''
sst <- sst[, -1]

sst$Season2[which(sst$Month > 11 | sst$Month < 3)] <- 'summer'
sst$Season2[which(sst$Month > 2 & sst$Month < 6)] <- 'autumn'
sst$Season2[which(sst$Month > 5 & sst$Month < 9)] <- 'winter'
sst$Season2[which(sst$Month > 8 & sst$Month < 12)] <- 'spring'

# Annual
sst0 <- data.frame(lat = sst$V2, lon = sst$V3, sst = sst$V4)
sst.ras <- SpatialPoints(cbind(sst0$lon, sst0$lat), proj4string = crs)
sst.ras <- rasterize(sst.ras, sp.ras, sst0$sst, fun = mean) 
sst.ras <- projectRaster(sst.ras, sp.ras)
m <- matrix(1, ncol = 5, nrow = 5)
sst.ras <- focal(sst.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T) 
sst.ras <- mask(sst.ras, static_pred[[1]]) 
writeRaster(sst.ras, filename = 'sst_annual_1.tif') 

# Season 2 (dic-ene-feb)
sst0 <- data.frame(Season = sst$Season2, lat = sst$V2, lon = sst$V3, sst = sst$V4)
sst0 <- aggregate(sst ~ Season + lon + lat, data = sst0, mean)
sst.list <- vector('list', 4)
for(i in 1:4){
  sub <- subset(sst0, Season == season[i])
  sst.ras <- SpatialPoints(cbind(sub$lon, sub$lat), proj4string = crs)
  sst.ras <- rasterize(sst.ras, sp.ras, sub$sst, fun = mean) 
  sst.ras <- projectRaster(sst.ras, sp.ras)
  sst.ras <- focal(sst.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T)
  sst.ras <- mask(sst.ras, static_pred[[1]]) 
  sst.list[[i]] <- sst.ras
}
sst.season <- stack(sst.list)
writeRaster(sst.season, filename = 'sst_season2_1.tif')


#-------------------------------------- TSM agrupación temporal ----------------------------------------

# Annual
sst_1 <- raster('sst_annual_1.tif')
sst_2 <- raster('sst_annual_2.tif')
sst <- mean(sst_1, sst_2)
writeRaster(sst, filename = 'sst.tif')

# Season 2
sst_1 <- stack('sst_season2_1.tif')
sst_2 <- stack('sst_season2_2.tif')
sst_season <- stack(mean(sst_1[[1]], sst_2[[1]]), mean(sst_1[[2]], sst_2[[2]]), 
                    mean(sst_1[[3]], sst_2[[3]]), mean(sst_1[[4]], sst_2[[4]]))
writeRaster(sst_season, filename = 'sst_season2.tif') 


#-------------------------------------- Coefficient Kd490 ---------------------------------------

# 1/m, 01-2003 to 12-2020 mensual, ~4 km)

# Descargar en 2 periodos(1 = 2003-2011, 2 = 2012-2020)
tur <- read.csv('DESCARGAR Y LEER EL CSV/tur_1.csv', skip = 2, header = F)
tur <- tur[which(!is.na(tur$V4)), ] 
tur$Year <- substr(tur$V1, 1, 4)
tur$Month <- substr(tur$V1, 6, 7)
tur$Month <- as.numeric(tur$Month)
tur$Season2 <- ''
tur <- tur[, -1]

tur$Season2[which(tur$Month > 11 | tur$Month < 3)] <- 'summer'
tur$Season2[which(tur$Month > 2 & tur$Month < 6)] <- 'autumn'
tur$Season2[which(tur$Month > 5 & tur$Month < 9)] <- 'winter'
tur$Season2[which(tur$Month > 8 & tur$Month < 12)] <- 'spring'

# Annual
tur0 <- data.frame(lat = tur$V2, lon = tur$V3, tur = tur$V4)
tur.ras <- SpatialPoints(cbind(tur0$lon, tur0$lat), proj4string = crs) 
tur.ras <- rasterize(tur.ras, sp.ras, tur0$tur, fun = mean) 
tur.ras <- projectRaster(tur.ras, sp.ras)
m <- matrix(1, ncol = 5, nrow = 5)
tur.ras <- focal(tur.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T) 
tur.ras <- mask(tur.ras, static_pred[[1]]) 
writeRaster(tur.ras, filename = 'tur_annual_1.tif') 

# Season 2
tur0 <- data.frame(Season = tur$Season2, lat = tur$V2, lon = tur$V3, tur = tur$V4)
tur0 <- aggregate(tur ~ Season + lon + lat, data = tur0, mean)
tur.list <- vector('list', 4)
for(i in 1:4){
  sub <- subset(tur0, Season == season[i])
  tur.ras <- SpatialPoints(cbind(sub$lon, sub$lat), proj4string = crs)
  tur.ras <- rasterize(tur.ras, sp.ras, sub$tur, fun = mean) 
  tur.ras <- projectRaster(tur.ras, sp.ras)
  tur.ras <- focal(tur.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T)
  tur.ras <- mask(tur.ras, static_pred[[1]]) 
  tur.list[[i]] <- tur.ras
}
tur.season <- stack(tur.list)
writeRaster(tur.season, filename = 'tur_season2_1.tif') 


#-------------------------------------- Productividad primaria ------------------------------------

# mgC/m2/dia, 01-2003 to 12-2020 mensual, ~4 km)

# Descargar en 2 periodos (1 = 2003-2011, 2 = 2012-2020)
pro <- read.csv('DESCARGAR Y LEER EL CSV/pp_1.csv', skip = 2, header = F)
pro <- pro[which(!is.na(pro$V4)), ] 
pro$Year <- substr(pro$V1, 1, 4)
pro$Month <- substr(pro$V1, 6, 7)
pro$Month <- as.numeric(pro$Month)
pro$Season2 <- ''
pro <- pro[, -c(1:2)]

pro$Season2[which(pro$Month > 11 | pro$Month < 3)] <- 'summer'
pro$Season2[which(pro$Month > 2 & pro$Month < 6)] <- 'autumn'
pro$Season2[which(pro$Month > 5 & pro$Month < 9)] <- 'winter'
pro$Season2[which(pro$Month > 8 & pro$Month < 12)] <- 'spring'

# Annual
pro0 <- data.frame(lat = pro$V3, lon = pro$V4, pro = pro$V5)
pro.ras <- SpatialPoints(cbind(pro0$lon, pro0$lat), proj4string = crs)
pro.ras <- rasterize(pro.ras, sp.ras, pro0$pro, fun = mean) 
pro.ras <- projectRaster(pro.ras, sp.ras)
m <- matrix(1, ncol = 5, nrow = 5)
pro.ras <- focal(pro.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T) 
pro.ras <- mask(pro.ras, static_pred[[1]]) 
writeRaster(pro.ras/1000, filename = 'pro_annual_1.tif') 

# Season 2
pro0 <- data.frame(Season = pro$Season2, lat = pro$V3, lon = pro$V4, pro = pro$V5)
pro0 <- aggregate(pro ~ Season + lon + lat, data = pro0, mean)
pro.list <- vector('list', 4)
for(i in 1:4){
  sub <- subset(pro0, Season == season[i])
  pro.ras <- SpatialPoints(cbind(sub$lon, sub$lat), proj4string = crs)
  pro.ras <- rasterize(pro.ras, sp.ras, sub$pro, fun = mean) 
  pro.ras <- projectRaster(pro.ras, sp.ras)
  pro.ras <- focal(pro.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T)
  pro.ras <- mask(pro.ras, static_pred[[1]]) 
  pro.list[[i]] <- pro.ras
}
pro.season <- stack(pro.list)
writeRaster(pro.season/1000, filename = 'pro_season2_1.tif') 


#-------------------------------------- Kd490 y Pp agrupación temporal -------------------------------------

# Annual
tur_1 <- raster('tur_annual_1.tif')
tur_2 <- raster('tur_annual_2.tif')
tur <- mean(tur_1, tur_2)
writeRaster(tur, filename = 'tur.tif')
pro_1 <- raster('pro_annual_1.tif')
pro_2 <- raster('pro_annual_2.tif')
pro <- mean(pro_1, pro_2)
writeRaster(pro, filename = 'pro.tif')

# Season 2
tur_1 <- stack('tur_season2_1.tif')
tur_2 <- stack('tur_season2_2.tif')
tur_season <- stack(mean(tur_1[[1]], tur_2[[1]]), mean(tur_1[[2]], tur_2[[2]]), 
                    mean(tur_1[[3]], tur_2[[3]]), mean(tur_1[[4]], tur_2[[4]]))
writeRaster(tur_season, filename = 'tur_season2.tif') 
pro_1 <- stack('pro_season2_1.tif')
pro_2 <- stack('pro_season2_2.tif')
pro_season <- stack(mean(pro_1[[1]], pro_2[[1]]), mean(pro_1[[2]], pro_2[[2]]), 
                    mean(pro_1[[3]], pro_2[[3]]), mean(pro_1[[4]], pro_2[[4]]))
writeRaster(pro_season, filename = 'pro_season2.tif') 


#-------------------------------------- Frentes térmicos de TSM ------------------------------------------

# Derivado de TSM, en ºC

# Annual
sst <- raster('sst.tif')
Front <- detectFronts(sst, method = 'median_filter')
Front <- focal(Front, m, fun = mean, NAonly = T, pad = T, na.rm = T) 
Front <- mask(Front, sst)
writeRaster(Front, filename = 'sstf.tif') 

# Season 2
sst_season <- stack('sst_season2.tif')
Front0 <- stack()
for(i in 1:4){
  sst <- sst_season[[i]]
  Front <- detectFronts(sst)
  Front <- focal(Front, m, fun = mean, NAonly = T, pad = T, na.rm = T) 
  Front <- mask(Front, sst) 
  Front0 <- stack(Front0, Front)
}
writeRaster(Front0, filename = 'sstf_season2.tif')


#-------------------------------------- Rasters scenopoéticos finales ------------------------------------------

# Predictores scenopoéticos finales

Names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 
           'Primary_productivity')
# Annual
static_pred <- stack('static_pred.tif')
sst <- raster('sst.tif')
sstf <- raster('sstf.tif')
tur <- raster('tur.tif')
pro <- raster('pro.tif')
annual.stack <- stack(static_pred, sst, sstf, tur, pro)
names(annual.stack) <- Names
writeRaster(annual.stack, filename = 'env.tif') 

# Season 2
sst <- stack('sst_season2.tif')
sstf <- stack('sstf_season2.tif')
tur <- stack('tur_season2.tif')
pro <- stack('pro_season2.tif')
for(i in 1:4){
  sst0 <- sst[[i]]
  sstf0 <- sstf[[i]]
  tur0 <- tur[[i]]
  pro0 <- pro[[i]]
  seasonal.stack <- stack(static_pred, sst0, sstf0, tur0, pro0)
  names(seasonal.stack) <- Names
  writeRaster(seasonal.stack, filename = paste(season[i], '.tif', sep = ''))
}


#-------------------------------------- Rasters finales para modelos regionales ------------------------------

# Rasters finales para los modelos regionales de G. galeus y N. cepedianus

# Leer predictores finales para G. galeus
predictors <- stack('env.tif')

# Distancia a colonias de pinnípedos
mmcol <- read.csv('M. mamm. colonies.csv')
mmcol <- subset(mmcol, SPECIES != 'Pontoporia blainvillei')
coord <- cbind(mmcol$LONG, mmcol$LAT)
coord.p <- SpatialPoints(coord, proj4string = crs)
dist.x <- distanceFromPoints(predictors[[1]], coord.p)
mmcol.ras <- mask(dist.x/1000, predictors[[1]])

Names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 'Primary_productivity', 'Distance_to_colonies')

# Anual
annual.stack <- stack(predictors, mmcol.ras)
names(annual.stack) <- Names
writeRaster(annual.stack, filename = 'env.tif')

# Estacional
for(i in 1:4){
  Season <- stack(paste(season[[i]], '.tif', sep = ''))
  seasonal.stack <- stack(Season, mmcol.ras)
  names(seasonal.stack) <- Names
  writeRaster(seasonal.stack, filename = paste(season[i], '.tif', sep = ''))
}


#-------------------------------------- FIN ------------------------------------------

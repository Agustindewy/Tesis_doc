
# Tesis Doctoral De Wysiecki (2023) 

# Preparación de los datos y modelado de nicho global de Galeorhinus galeus

library(sp)
library(raster)
library(sf)
library(spThin)
library(rgdal)
library(sqldf)
library(maps)
library(ggplot2)
library(rgeos)
library(dismo)
library(blockCV)
library(kuenm)
library(ntbox)# disponible en GitHub (luismurao/ntbox)
library(rSDM) # disponible en GitHub (Pakillo/rSDM)
library(ENMeval) # disponible en GitHub (jamiemkass/ENMeval)
library(SDMtune) # disponible en GitHub (ConsBiol-unibern/SDMtune)

setwd('INDICAR DIRECTORIO DE TRABAJO')

# Proyección
crs <- CRS("+proj=longlat +datum=WGS84")

# Especie
species <- 'Galeorhinus galeus'

# Semilla
set.seed(111)

# Crear carpeta para la especie
dir.create(paste(species))

#-------------------------------- Espacio G --------------------------------------------

# Leer predictores y nombrarlos
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope', 'Thermal_fronts')
names(env) <- var_names

# Descartar batimetría y pendiente por discrepancias en el margen de la plataforma continental (De Wysiecki et al. 2022) 
var_set <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity', 
             'Distance_to_coast', 'Thermal_fronts')

# Espacio G
env.G <- env[[var_set]]


#-------------------------------- Preparación ----------------------------------

# Leer ocurrencias
dat <- read.csv('C:/Users/User1/OneDrive/Doc/Manuscrito/Ocurrencias_doc.csv') # no están todas ya que algunas no son de libre acceso
dat <- subset(dat, Species == species)
dat <- subset(dat, Source1 %in% c('GBIF', 'OBIS', 'Government data', 'Published literature', 'Grey literature', 'Social media', "Professional contact"))

# Eliminar duplicados
dups <- duplicated(dat[c('Longitude', 'Latitude')])
occ_cal <- dat[!dups, ]

# Filtrado espacial
train_thin <- thin(occ_cal, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                   thin.par = 50, reps = 1, locs.thinned.list.return = T, write.files = F, write.log.file = F)
occ_cal <- data.frame(Species = species, Longitude = round(train_thin[[1]]$Longitude, 2),
                      Latitude = round(train_thin[[1]]$Latitude, 2))

# Espacio M inicial para filtrado ambiental
# algunos datos caen en Nueva Zelanda y se necesitan buffers que superan el meridiano 180°
# por ello recentrar los puntos a 110° 
occ_cal$Longitude2 <- occ_cal$Longitude + 110
occ_cal$Longitude2 <- ifelse(occ_cal$Longitude2 > 180, (occ_cal$Longitude2 - 180) - 180, occ_cal$Longitude2)
coord <- data.frame(Longitude = occ_cal$Longitude2, Latitude = occ_cal$Latitude)
occ.tot = st_as_sf(coord, coords = c("Longitude", "Latitude"), crs = 4326)
occ.buff <- st_buffer(occ.tot, 1000000) # buffer de 1000 km
# recentrar el raster también a 110° 
env1 <- crop(env.G, extent(-180, 70, -90, 90))
env2 <- crop(env.G, extent(70, 180, -90, 90))   
extent(env1) <- c(-70, 180, -90, 90)
extent(env2) <- c(-180, -70, -90, 90)
env0 <- merge(env1, env2)
names(env0) <- var_set
env.M <- crop(env0, occ.buff) # cortar el raster al buffer
env.M <- mask(env.M, occ.buff) # enmascarar el raster al buffer
env.M <- stack(env.M)

# Rescatar registros que caigan fuera de alguna capa del raster y llevarlos al pixel con dato mas cercano 
for(i in 1:length(env.M@layers)){
  spp <- SpatialPoints(occ_cal[, c('Longitude2', 'Latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M, layer = i, move = T)
  occ_cal$Longitude2 <- round(spp_corrected@coords[, 1], 3)  # reemplazar las coordenadas, incluidas las nuevas, si las hubiera
  occ_cal$Latitude <- round(spp_corrected@coords[, 2], 3)
} # plots aparecen si se aplica alguna corrección

# Remover duplicados de nuevo en caso de que hayan surgido
dups <- duplicated(occ_cal[c('Longitude2', 'Latitude')])
occ_cal <- occ_cal[!dups, ]

# Filtrado ambiental
source('envSample.R')
coords <- SpatialPoints(occ_cal[, c('Longitude2', 'Latitude')], crs)
env.data <- extract(env.M, coords)
env.data <- as.data.frame(env.data)
coords <- envSample(coords, filters = list(env.data$Surface_temperature, env.data$Distance_to_coast),
                    res = list(0.5, 1), do.plot = T)
occ_cal <- merge(occ_cal, coords, by.x = c('Longitude2', 'Latitude'), by.y = c('lon', 'lat'))

# Outliers en espacio ambiental (Cobos et al., 2018)
variables_values <- na.omit(values(env.M))
occ_variables <- na.omit(cbind(as.numeric(row.names(occ_cal)), extract(env.M, SpatialPoints(occ_cal[, c('Longitude2', 'Latitude')], crs))))
occ_variables <- as.data.frame(occ_variables)
if(dim(variables_values)[1] > 10000) {
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ]
}

# Plotear todas las combinaciones e ir buscando outliers y escribir número abajo
for(i in 1:dim(env.M)[3]){
  for(j in 1:dim(env.M)[3]){
    par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
    plot(variables_values[, i], variables_values[, j], col = 'grey65', pch = 1,
         xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    points(occ_variables[, i + 1], occ_variables[, j + 1], col = 'blue', pch = 19, cex = 1.5)
    legend('bottomright', legend = c('Region of interest', 'Occurrences'),
           pch = c(1, 19), col = c('grey65', 'black'), bty = 'n')
    plot(variables_values[, i], variables_values[, j], col = 'grey65',
         pch = 1, xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    text(occ_variables[, i + 1], occ_variables[, j + 1], occ_variables[, 1], cex = 1, col = 'blue')
    legend('bottomright', legend = 'Occurrence ID', bty = 'n')
  }
}

# Remover outliers
occ_variables <- subset(occ_variables, !V1 %in% c(99999)) # poner los números aquí
occ_cal <- occ_cal[which(row.names(occ_cal) %in% as.character(occ_variables$V1)), ]

# Puntos de calibración finales
occ_cal <- rbind(occ_cal[, c(3, 4, 2)])
write.csv(occ_cal, paste(species, '/', 'Calibration_points.csv', sep = ''), row.names = F)


#-------------------------------- Área de calibración ----------------------------------------------

occ_cal <- read.csv(paste(species, '/', 'Calibration_points.csv', sep = ''))

# Espacio M final (área de calibración)
# Volvear a recentrar por datos en Nueva Zelanda
occ_cal$Longitude2 <- occ_cal$Longitude + 110
occ_cal$Longitude2 <- ifelse(occ_cal$Longitude2 > 180, (occ_cal$Longitude2 - 180) - 180, occ_cal$Longitude2)
coord <- data.frame(Longitude = occ_cal$Longitude2, Latitude = occ_cal$Latitude)
occ.tot = st_as_sf(coord, coords = c("Longitude", "Latitude"), crs = 4326)
occ.buff <- st_buffer(occ.tot, 1000000) # buffer de 1000 km
# Recentrar el raster
env1 <- crop(env.G, extent(-180, 70, -90, 90))
env2 <- crop(env.G, extent(70, 180, -90, 90))   
extent(env1) <- c(-70, 180, -90, 90)
extent(env2) <- c(-180, -70, -90, 90)
env0 <- merge(env1, env2)
names(env0) <- var_set
env.M <- crop(env0, occ.buff) # cortar el raster al buffer
env.M <- mask(env.M, occ.buff) # enmascarar el raster al buffer

# Una vez cortado volverlo a centrar en 0° 
env.M <- extend(env.M, extent(-180, 180, -62.75, 73.75)) # de vuelta a -180°, 180°
env.M1 <- crop(env.M, extent(-180, -70, -62.75, 73.75))
env.M2 <- crop(env.M, extent(-70, 180, -62.75, 73.75))   
extent(env.M1) <- c(70, 180, -62.75, 73.75)
extent(env.M2) <- c(-180, 70, -62.75, 73.75)
env.M <- merge(env.M2, env.M1)
names(env.M) <- var_set
env.M <- stack(env.M)

# quitar algunas zonas inchorentes que quedan dentro de los buffer al atravezar la tierra, por ej Mar Muerto
mar_muerto <- matrix(c(23.55, 31.32, 47.24, 37.60, 23.55, 43.28, 51.13, 43.70, 35.53, 43.28), 5, 2)
mar_muerto <- SpatialPolygons(list(Polygons(list(Polygon(mar_muerto)), ID = 1)), proj4string = crs)
env.M <- mask(env.M, mar_muerto, inverse = T)

# Chequear correlación en el espacio M
layerStats(env.M, 'pearson', na.rm = T) # sin correlación significativa

# Guardar al directorio
writeRaster(env.M, filename = paste(species, '/', 'Calibration_areas.tif', sep = ''))


#-------------------------------- Calibración, evaluación y selección de modelos ---------------------

# Puntos de calibración
occ_cal <- read.csv(paste(species, '/', 'Calibration_points.csv', sep = ''))
occ_cal <- occ_cal[, c('Longitude', 'Latitude')]
colnames(occ_cal) <- c('longitude', 'latitude')

# Áreas de calibración
env.M <- stack(paste(species, '/','Calibration_areas.tif', sep = ''))
names(env.M) <- var_set

# Puntos del background
bg_points <- as.data.frame(randomPoints(env.M[[1]], n = 10000))
colnames(bg_points) <- c('longitude', 'latitude')

# Puntos de calibración + background
all_pts <- rbind(occ_cal, bg_points)
all_pts <- SpatialPoints(all_pts, proj4string = crs)

# Chequear el rango efectivo de autocorrelación espacial
eff.range <- spatialAutoRange(rasterLayer = env.M, sampleNumber = 10000, doParallel = T, showPlots = T)
median_range <- 99999 # ingresar aquí el rango óptimo 

# Diseño espacial de bloques basado en 'blockCV' de Valavi et al. (2019)
sp_block <- spatialBlock(speciesData = all_pts, rasterLayer = env.M[[1]], theRange = median_range, k = 5, selection = 'random')

# Obtener la información de los bloques para el modelado
user.grp <- list(occs.grp = sp_block$foldID[1:nrow(occ_cal)],
                 bg.grp = sp_block$foldID[(nrow(occ_cal) + 1):length(sp_block$foldID)])

# Función para incluir estadistico del ROC parcial en la evaluación ENM
proc <- function(vars) {
  proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(proc_auc_ratio = proc$pROC_summary[1], 
                    proc_pval = proc$pROC_summary[2], row.names = NULL)
  return(out)
} 

# Correr y evaluar los modelos candidatos
mod_global <- ENMevaluate(occs = occ_cal, envs = env.G, bg = bg_points, user.grp = user.grp, 
                          algorithm = 'maxent.jar', partitions = 'user', user.eval = proc, doClamp = F, 
                          tune.args = list(fc = c('LQ', 'LP', 'LQP'), rm = c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10)))

# Selección de modelos según criterio de 3 pasos en Cobos et al. (2019)
res <- eval.results(mod_global)
opt.seq <- res %>% 
  filter(delta.AICc == min(delta.AICc)) %>%
  filter(proc_pval.avg == min(proc_pval.avg)) %>%
  filter(or.10p.avg == min(or.10p.avg))

# Modelo óptimo
mod.seq <- eval.models(mod_global)[[opt.seq$tune.args]]

# Contribución de los predictores al modelo
mod_global@variable.importance

# Predicción del modelo
pred.seq <- eval.predictions(mod_global)[[opt.seq$tune.args]]
writeRaster(pred.seq, filename = paste(species, '/','Predictions.tif', sep = ''))


#-------------------------------- Riesgo de extrapolación ---------------------------------

# Basado en el análisis de Mobility Oriented Parity (MOP) de Owens et al. (2013) 

# Espacio G
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope', 'Thermal_fronts')
names(env) <- var_names
var_set <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity', 
             'Distance_to_coast', 'Thermal_fronts')
env.G <- env[[var_set]]

# Línea de costa (polígonos espaciales) - tomado de la base de datos de costas GSHHG del 15 de junio de 2017
coastline <- 'C:/Users/User1/Desktop/D/ENV SDM/Coast shape'
coast <- st_read(dsn = coastline, layer = 'GSHHS_f_L1_World')
env.G <- crop(env.G, coast) # cortar regiones del Ártico y Antártida

# Espacio M
env.M <- stack(paste(species, '/', 'Calibration_areas.tif', sep = ''))
names(env.M) <- var_set
mop_analysis <- mop(M_stack = env.M, G_stack = env.G, percent = 10, comp_each = 2000, parallel = T)
writeRaster(mop_analysis, filename = paste(species, '/', 'Mop.tif', sep = ''))


#-------------------------------- FIN -------------------------------------

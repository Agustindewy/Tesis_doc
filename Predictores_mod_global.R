
# Predictores scenopoéticos para C. brachyurus y C. taurus: descarga y preparación
# Predictores a escala global y anual

library(sdmpredictors)
library(raster)
library(grec)

setwd('DIRECTORIO DE TRABAJO')

#----------------------------------------- Predictores ------------------------------------------------

# Bio-Oracle
list_layers(c('Bio-ORACLE'))[, 2:4] # explorar capas

# Chequear correlación a escala global
layers_correlation(c('BO2_tempmean_bdmean', 'BO_sstmean', 'BO2_ppmean_bdmean', 'BO_damean', 'BO_chlomean',
                     'BO2_salinitymean_bdmean', 'BO_bathymean'))
# BO_chlomean vs BO_damean corr > 0.8, descartamos BO_chlomean

bio <- load_layers(c('BO2_tempmean_bdmean', 'BO_sstmean', 'BO2_ppmean_bdmean', 'BO_damean',
                     'BO2_salinitymean_bdmean', 'BO_bathymean'))

# Distancia a la costa del Global Self-consistent, Hierarchical, High-resolution Geography Database
# Descargado desde https://www.soest.hawaii.edu/pwessel/gshhg/
dis <- raster('DESCARGAR Y LEER EL RASTER.tif') # leer el tif
dis <- aggregate(dis, fact = 2, fun = mean) # reducir la resolución para igualar los predictores del Bio-ORACLE

# Pendiente (creada del raster de batimetría)
slo <- terrain(bio[['BO_bathymean']], opt = 'slope', unit = 'degrees', neighbors = 4)
slo <- projectRaster(slo, bio)

# Frentes térmicos (creada del raster de temperatura superficial)
sst <- bio[['BO_sstmean']]
tfr <- detectFronts(sst, method = 'median_filter')
m <- matrix(1, ncol = 5, nrow = 5)
tfr <- focal(tfr, m, fun = mean, NAonly = T, pad = T, na.rm = T) 
tfr <- mask(tfr, sst)

# Stack final de predictores
env <- stack(bio, dis, slo, tfr)

# Chequear correlación otra vez con nuevos predictores
layerStats(env, 'pearson', na.rm = T) # sin correlación

# Chequear si los valores ambientales en cada caso tienen sentido
env # batimetría tiene algunos valores positivos, convertirlos a NA
env[[6]][env[[6]] > -1] <- NA

writeRaster(env, filename = 'predictors.tif') # guardar el tif

#-------------------------------------- FIN ------------------------------------------

####################################################
## PATTSI -  SPATIAL PATTERN SUITABILITY ANALYSIS ##
####################################################

##############################
## Load required R packages ##
##############################

require(raster)
require(Matrix)

####################################
## Distance-dependent suitability ##
####################################

# Goal:   focal pixel suitability index (SI) depends on the distance of a (numerical) habitat variable
# Method: generate a moving window where each position is weighted by its distance from the focal pixel
#         compute the focal pixel SI as the highest value of a habitat variable within this moving window
#         (here: circular window with five equal distance bands, radius/weights elicited from experts)


## Step 1: define a function that returns a circular matrix of given radius and resolution (Scroggie 2012)
# assigns value "1" if matrix position <= radius and value "NA" if matrix position > radius
make_circ_filter <- function(radius, res){
  circ_filter <- matrix(NA, nrow=1+(2*radius/res), ncol=1+(2*radius/res))
  dimnames(circ_filter)[[1]] <- seq(-radius, radius, by=res)
  dimnames(circ_filter)[[2]] <- seq(-radius, radius, by=res)
  sweeper <- function(mat){
    for(row in 1:nrow(mat)){
      for(col in 1:ncol(mat)){
        dist <- sqrt((as.numeric(dimnames(mat)[[1]])[row])^2 +
                     (as.numeric(dimnames(mat)[[1]])[col])^2)
        if(dist<=radius) {mat[row, col]<-1}
      }
    }
    return(mat)
  }
  out <- sweeper(circ_filter)
  return(out)
}

## Step 2: apply function to generate five matrices with different radii (= distance bands)
res <- your.pixel # resolution (= pixel size)
udt <- your.radius # upper distance threshold (= moving window radius, multiple of r)
m.vf <- make_circ_filter(udt, res) # distance band 'very far' (= udt)
m.f <- make_circ_filter((udt/5)*4, res) # distance band 'far' (= 4/5 of udt)
m.m <- make_circ_filter((udt/5)*3, res) # distance band 'medium' (= 3/5 of udt)
m.c <- make_circ_filter((udt/5)*2, res) # distance band 'close' (= 2/5 of udt)
m.vc <- make_circ_filter((udt/5), res) # distance band 'very close' (= 1/5 of udt)

# replace 'value==1' with unique temp value in ascending order from largest to smallest matrix
m.vf[m.vf == 1] <- 1 
m.f[m.f == 1] <- 2
m.m[m.m == 1] <- 3
m.c[m.c == 1] <- 4
m.vc[m.vc == 1] <- 5

## Step 3: combine the five matrices into one (two at a time starting with the smallest):
a.1 <- array(NA, dim(m.c), dimnames(m.c)) # create temp array of size = larger matrix
a.1[rownames(m.vc), colnames(m.vc)] <- m.vc # ... with values = smaller matrix
m.c <- pmax(m.c, a.1, na.rm = TRUE) # combine values: larger matrix + temp array
a.2 <- array(NA, dim(m.m), dimnames(m.m)) # repeat with: output + next-larger matrix
a.2[rownames(m.c), colnames(m.c)] <- m.c
m.m <- pmax(m.m, a.2, na.rm = TRUE)
a.3 <- array(NA, dim(m.f), dimnames(m.f))
a.3[rownames(m.m), colnames(m.m)] <- m.m
m.f <- pmax(m.f, a.3, na.rm = TRUE)
a.4 <- array(NA, dim(m.vf), dimnames(m.vf))
a.4[rownames(m.f), colnames(m.f)] <- m.f
m.band <- pmax(m.vf, a.4, na.rm = TRUE)
m.band

## Step 4: replace temp values with expert-elicited weight for each of five distance bands (in CSV file)
csv.D = read.csv("[...]/your-weights.csv") # must contain column header (numeric or text as below)
m.band[m.band == 1] <- mean(csv.D$very.far) # here, mean() takes average of multiple expert
m.band[m.band == 2] <- mean(csv.D$far)
m.band[m.band == 3] <- mean(csv.D$medium)
m.band[m.band == 4] <- mean(csv.D$close)
m.band[m.band == 5] <- mean(csv.D$very.close)

## Step 5: perform moving window analysis using function "focal {raster}" with parameters:
r = raster("[...]/your-habitat-variable.tif") # raster layer with numerical habitat variable (e.g. quality index)
w = m.band # banded weights matrix is used as moving window
fun = max # focal pixel takes weighted maximum value of habitat variable within moving window

# WARNING! The following process may take several hours depending on the size of r and w
r.D <- focal(r, w, fun, na.rm = TRUE, pad = FALSE, padValue = NA) # na.rm = TRUE ignores NoData
r.SI.D <- mask(r.D, r) # extract by r to remove padded edges
writeRaster(r.SI.D, filename = paste("[...]/your-distance-SI.tif", sep="")) # save ouput raster layer


#######################################
## Composition-dependent suitability ##
#######################################

# Goal:   focal pixel suitability index (SI) depends on the cover of a (numerical) habitat variable
# Method: for each focal pixel, compute the average value of a habitat variable within a moving window
#         reclassify the average value to a numerical focal pixel SI from a reclassification matrix
#         (here: circular window, five equal average value classes, radius/reclassifcation matrix elicited from experts)


## Step 1: perform moving window analysis using the "ESRI ArcGIS Focal Statistics" tool with parameters:
# in_raster = "your-habitat-variable.tif": numerical habitat variable (e.g. quality index)
# neighborhood = NbrCircle, radius = udt: circular moving window with radius = udt (see above)
# statistics_type = mean: focal pixel takes mean value of habitat variable within moving window
# ignore_nodata = TRUE: ignores NoData

######################################################
## # Currently not functioning: implementation using function "focal {raster}" with parameters:
## r = raster("your-habitat-variable.tif")
## w = m.vf # circular matrix with radius = udt (see above) is used as moving window
## fun = mean # focal pixel takes mean value of habitat variable within moving window
## r.D <- focal(r, w, fun, na.rm = TRUE, pad = FALSE, padValue = NA)
## # returns empty result due to "NA" values at matrix positions > radius
## # changing "NA" to "0" values distorts computations of mean value
######################################################

r.C = raster("[...]/your-mean-value.tif") # output raster layer with focal mean values 
r.C <- mask(r.C, r) # extract by r to remove padded edges

## Step 2: determine minimum and maximum focal mean values
fm.min <- minValue(r) # all habitat variables within moving window have lowest possible value
fm.max <- maxValue(r) # all habitat variables within moving window have highest possible value

## Step 3: replace focal mean values with expert-elicited SI for each of average value five classes (in CSV file)
csv.C = read.csv("[...]/your-reclass.csv") # must contain column header (text or numeric as below)
rcl.C <- c(minValue(r.C), fm.min + (fm.max - fm.min) * 0.20, mean(csv.C$X20),
  fm.min + (fm.max-fm.min) * 0.20, fm.min + (fm.max-fm.min) * 0.40, mean(csv.C$X40),
  fm.min + (fm.max-fm.min) * 0.40, fm.min + (fm.max-fm.min) * 0.60, mean(csv.C$X60),
  fm.min + (fm.max-fm.min) * 0.60, fm.min + (fm.max-fm.min) * 0.80, mean(csv.C$X80),
  fm.min + (fm.max-fm.min) * 0.80, fm.min + (fm.max-fm.min) * 1.00, mean(csv.C$X100))
m.rcl.C <- matrix(rcl.C, ncol=3, byrow=TRUE) # reclassification matrix (3 columns: "from" / "to" mean values, "SI" )
r.SI.C <- reclassify(r.C, m.rcl.C) # apply function "reclassify {raster}"
writeRaster(r.SI.C, filename = paste("[...]/your-composition-SI.tif", sep="")) # save ouput raster layer


#########################################################
## Combined composition/distance-dependent suitability ##
#########################################################

# Goal:   focal pixel suitability index (SI) depends on both distance and cover of a (numerical) habitat variable
# Method: generate a moving window where each position is weighted by its distance from the focal pixel
#         for each focal pixel, compute the average value of a habitat variable within this weighted moving window
#         reclassify the average value to a numerical focal pixel SI from a reclassification matrix
#         (here: circular window, five equal distance bands, five equal average value classes, radius/weights/reclassifcation matrix elicited from experts)


# Step 1: save banded weights matrix (as TXT, add header with number of rows/columns after export)
m.exp <- m.band
m.exp[is.na(m.exp)] <- 0 # change NAs to 0s
write.table(m.exp, file = "[...]/your-weights-matrix.txt", row.names = F, col.names = F) # remove row and column names

## Step 2: perform moving window analysis using the ESRI ArcGIS Focal Statistics tool with parameters:
# in_raster = "your-habitat-variable.tif": numerical habitat variable (e.g. quality index)
# neighborhood = NbrWeight: select exported banded weights matrix "your.txt"
# statistics_type = mean: focal pixel takes weighted mean value of habitat variable within moving window
# ignore_nodata = TRUE: ignores NoData

r.DC = raster("[...]/your-weightedmean-value.tif") # output raster layer with focal weighted mean values 
r.DC <- mask(r.DC, r) # extract by r to remove padded edges

## Step 3: determine minimum and maximum focal weightedmean values as
# number of values in each distance band * (minValue(r)/maxValue(r) * weight for distance band) / total number of values in moving window
m.t <- table(m.band) # number of pixels in each distance band
fwm.min <- ((m.t[5] * (minValue(r) * mean(csv.D$very.close))) + (m.t[4] * (minValue(r) * mean(csv.D$close))) 
            + (m.t[3] * (minValue(r) * mean(csv.D$medium))) + (m.t[2] * (minValue(r) * mean(csv.D$far))) 
            + (m.t[1] * (minValue(r) * mean(csv.D$very.far)))) / nnzero(m.band, na.counted = FALSE)
fwm.max <- ((m.t[5] * (maxValue(r) * mean(csv.D$very.close))) + (m.t[4] * (maxValue(r) * mean(csv.D$close))) 
            + (m.t[3] * (maxValue(r) * mean(csv.D$medium))) + (m.t[2] * (maxValue(r) * mean(csv.D$far))) 
            + (m.t[1] * (maxValue(r) * mean(csv.D$very.far)))) / nnzero(m.band, na.counted = FALSE)
as.vector(round(fwm.min, digits = 1)); as.vector(round(fwm.max, digits = 1)) # print as rounded vector

## Step 4: replace focal weightedmean values with expert-elicited SI for each of five average value classes (in CSV file)
rcl.DC <- c(minValue(r.DC), fwm.min + (fwm.max - fwm.min) * 0.20, mean(csv.C$X20),
            fwm.min + (fwm.max-fwm.min) * 0.20, fwm.min + (fwm.max-fwm.min) * 0.40, mean(csv.C$X40),
            fwm.min + (fwm.max-fwm.min) * 0.40, fwm.min + (fwm.max-fwm.min) * 0.60, mean(csv.C$X60),
            fwm.min + (fwm.max-fwm.min) * 0.60, fwm.min + (fwm.max-fwm.min) * 0.80, mean(csv.C$X80),
            fwm.min + (fwm.max-fwm.min) * 0.80, fwm.min + (fwm.max-fwm.min) * 1.00, mean(csv.C$X100))
m.rcl.DC <- matrix(rcl.DC, ncol=3, byrow=TRUE) # reclassification matrix (3 columns: "from" / "to" mean values, "SI" )
r.SI.DC <- reclassify(r.DC, m.rcl.DC) # apply function "reclassify {raster}"
writeRaster(r.SI.DC, filename = paste("[...]/your-distance-composition-SI.tif", sep="")) # save ouput raster layer

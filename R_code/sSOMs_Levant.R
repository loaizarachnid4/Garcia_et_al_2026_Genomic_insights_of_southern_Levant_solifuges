##############
#SuperSOMs (Pyron 2023)
#############
####The original code was derived from https://github.com/rpyron/delim-SOM, the only modifications 
###from the original code are to analyze less than 16 samples. Namely, those functions are "DNA.SOM.small", "match.labels.small"
###and "plotK.small". To use these functions refer to the "kohonen_code_modified.R" file to load these functions. 

###code to generate maps of study region from shape file is also included at the bottom of this script file. 

###This code was used in Garcia, E.L., Ganem, Z., Kulkarni, S.S., Perl, I., Sharma, P.P., Gavish-Regev, E. (2026). 
#Genomic insights unveil taxonomic incongruities and evolutionary origins of solifuges across the southern Levant
#Molecular Phylogenetics and Evolution###


setwd("/path/to/working/directory")
#dependencies must be installed
source("kohonen_code_modified.R")
library(vcfR)
library(RColorBrewer)
library(adegenet);library(maps);library(scales)
library(conStruct);library(poppr);library(kohonen)
library(lsr);library(combinat);library(viridis)
set.seed(1)

#read in structure file
sev <- "solp_adegenet_first.stru"

# Read the file
lines <- readLines(sev)

# Count individuals (each diploid individual has 2 rows)
n.ind <- length(lines) / 2

# Count loci (assumes first line is representative)
first.line <- strsplit(lines[1], "\\s+")[[1]]
n.loci <- length(first.line) - 1  # Assuming 1st column is ID or label

#read in structure file based on the parameters inferred above 
a <- read.structure(sev, n.ind=n.ind, n.loc=n.loci, onerowperind=FALSE, col.lab=1, col.pop=0, col.others=NULL, row.marknames=NULL, NA.char=-9, pop=NULL, ask=FALSE, quiet=FALSE)

##uncomment below to implement a missing data cut off if one wasn't applied before. If too much missing data, an error
##may occur

#a = missingno(a, type = "loci", cutoff = 0.20)

#Convert allele frequences to matrix
alleles <- makefreq(a)

#####################Parameters for runs - for samples more than 16##############
#Size of Grid
g <- round(sqrt(5*sqrt(length(rownames(alleles)))))#common rule of thumb

#Create an output grid of size sqrt(n)
som_grid <- somgrid(xdim = g,
                    ydim = g,
                    topo="hexagonal",
                    neighbourhood.fct = "gaussian")

#Number of Replicates - can increase if you like
n <- 100

#Number of steps - doesn't usually matter beyond ~100
m <- 100

##DEFAULT - if you get an error it is likely due to one of your samples having too much missing or too few samples
#run on single layer
res <- DNA.SOM()

###set up for structure-like plot
labels <- match.labels(alleles)#get DAPC labels
q_mat <- match.k(res,labels)#get admixture coefficients

##########################SMALL SAMPLE SIZES or if you get an error from above###################

#####MODIFIED CODE####
res <- DNA.SOM.small() 
###set up for small plots
labels <- match.labels.small(alleles)#get DAPC labels
q_mat <- match.k(res,labels)#get admixture coefficients


########## 
#change color according to family for maps
#daes
cols <- brewer.pal(n = 5, name = "Dark2")
#rhago
cols <- brewer.pal(n = 5, name = "Set2")
#solp
cols <- brewer.pal(n = 5, name = "Paired")
#karsch
cols <- brewer.pal(n = 5, name = "Set1")
#galeo
cols <- brewer.pal(n = 5, name = "Accent")
#gylip
cb_palette_17 <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442","#D55E00", "#CC79A7", "#999999", "#0072B2","#000000")
cols <- cb_palette_17
cols <- c(   "#4DAF4A" ,"#E41A1C","#377EB8","#984EA3", "#FF7F00")

tiff("daes75p_sSOM_k3_mr.tiff", units="in", width=12, height=8, res=300)
x <- q_mat[order(q_mat[,1]),]
z <- hclust(dist(x),"single")$order
make.structure.plot(admix.proportions = x[z,], 
                    sample.names = rownames(x[z,]), 
                    mar = c(8,4,2,2), 
                    layer.colors = cb_palette_17,
                    sort.by = 2)
dev.off()

###Optimize K plots####
tiff("daes75p_sSOM_k3_bicplots.tiff", units="in", width=12, height=8, res=300)
plotK(res)
dev.off()

###if supported clusters are >= 3 then delta BIC cannot be plotted
tiff("daes75p_sSOM_k3_bicplots.tiff", units="in", width=12, height=8, res=300)
plotK.small(res)
dev.off()

##############mapping admixture 
#this pulls coordinate data for matching sample names based on a seperate file with coordinate data with matching sample names
taxon <- read.csv("Daesiidae_localities_DNA.csv")
matches<- match(rownames(q_mat),taxon$Code...Record.Number)

#pulls only lat, long data
to.map <- taxon[matches,]
#xy<- cbind.data.frame(to.map$Y.Longitude,to.map$X.Latitude)
xy<- cbind(to.map$Longitude,to.map$Latitude)


#####run through this if NAs are produced (outgroups/samples out of the target region)
q_mat<- q_mat[!is.na(xy[,1]),]
xy <- xy[!is.na(xy[,1]),]

#remove individuals
q_mat <- q_mat[-1,]
xy <- xy[-1,]

##plot and export basic map with coodinates 
tiff("daes75p_sSOM_k3_admixplots.tiff", units="in", width=12, height=8, res=300)
par(mfrow=c(1,1),
    mar=c(0,0,0,0))
maps::map(database = 'world', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="white")
map.axes()
maps::map(database = 'world', xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), col="darkgreen",add=T)
maps::map(database = 'world', regions = c("Israel", "Gaza","Egypt", "Syria","Jordan"), xlim = range(xy[,1]) + c(-0.5,0.5), ylim = range(xy[,2]) + c(-0.5,0.5), add = T)
make.admix.pie.plot(q_mat,xy,layer.colors = cols,radii=2.5,add = T)
dev.off()

######updated mapping######
##this is to produce a seperate terrain map - must include the shape file for region of interest and
#figure out the bounding box 

library(ggplot2)
library(ggmap)
register_stadiamaps("657a5379-96bd-4ec7-a75f-73b6ff680fdf")
library(maps)
library(mapproj)
library(mapdata)
library(sp)
library(raster)
library(dismo)
library(sf)
library(readxl)

nc <- st_read("/World_Countries_6893716458424143584/World_Countries.shp")

# Transform nc to EPSG 3857 (Pseudo-Mercator, what Google uses)
nc_3857 <- st_transform(nc, 3857)

#map <- get_stadiamap(bbox =c(-14,6,75,43), zoom=7, maptype="stamen_terrain_background")
map <- get_stadiamap(bbox =c(32,29,37,34), zoom=7, maptype="stamen_terrain_background")

#function in modified kohonen_code.R
map <- ggmap_bbox(map)

map1 <- ggmap(map) + coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
  geom_sf(data = nc_3857 , inherit.aes = FALSE, color='black', alpha=0, size=0.5)

#export map to overlay
tiff("basicmap_Levant.tiff", units="in", width=8, height=8,res=300)
map1
dev.off()

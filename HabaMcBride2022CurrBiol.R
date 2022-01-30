### HabaMcBride2022CurrBiol.R ###

# R script to plot CQ11/bloodmeal data onto a map
# Haba and McBride 2022 Current Biology
# Origin and status of Culex pipiens mosquito ecotypes

#=====================================================================================
#  library
#=====================================================================================
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggimage)
library(reshape2)
library(rgdal)
library(ggthemes)
library("rnaturalearth");library("rnaturalearthdata") # for map

#=====================================================================================
#  working directory
#=====================================================================================
setwd("~/path/to/data/")

# ================================================================== #
# ##################  CQ11 analysis (Fig. 2a-d) #################### #
# ================================================================== #

#=====================================================================================
#  read CQ11 data + metadata
#=====================================================================================
data <- read.csv("Table S1. Published CQ11 Genotype Data Related to Figure 3.csv", 
                 header=T, stringsAsFactors = FALSE)

# check if import was successful
head(data)
data$Freq_mm # check if the imported data are numbers instead of strings
nrow(data)

#=====================================================================================
# Calculate heterozygosity deficit (F_IS) and minor allele frequency (MAF) and run HWE test
#=====================================================================================
# install.packages("HardyWeinberg")
library(HardyWeinberg)

# run test per row
for (i in 1:length(data$N_mm))
{
  x <- c(AA=data$N_mm[i],AB=data$N_pm[i],BB=data$N_pp[i])
  if (data$Freq_mm[i] == 1 | data$Freq_pp[i] == 1) # if there is only one allele, F_IS=0
  {
    # F_IS
    data$F_IS[i] <- 0
    # MAF
    data$MAF[i] <- 0
  } else
  {
    # F_IS
    HetCount_e <- HWChisq(x)$expected[2]
    HetCount_o <- x[2]
    
    F_IS= (1 - (HetCount_o/sum(x))/(HetCount_e/sum(x)))
    data$F_IS[i] <- F_IS
    # MAF
    data$MAF[i] <- min(data$Freq_m[i],data$Freq_p[i])
  }
  
  
  #HWE test
  HWstats = HWExact(x)
  data$HWpval[i] <- HWstats$pval
}
# head(data)

#=====================================================================================
#  store above vs below separately
#=====================================================================================
data.above <- subset(data, data$Microhabitat == "aboveground")
data.below <- subset(data, data$Microhabitat == "belowground")

#=====================================================================================
#  Fig. 2a. pip frequency pie charts on the map
#=====================================================================================
#===========================================
#  First, select a subset of individuals to plot
#===========================================
#==============
#  (For aboveground panel) selection for visualization purpose
#==============
pie.data.above.PieShown <- subset(data.above, data.above$Fig3A_PieShown == "X")
nrow(pie.data.above.PieShown)

#==============
#  (For belowground panel) plot all but a representative where several share the same GPS coordinates
#==============
pie.data.below <- data.below
pie.data.below.PieShown <- subset(pie.data.below, pie.data.below$Fig3A_PieShown == "X")
nrow(pie.data.below.PieShown)

#==============
#  get basemap
#==============
theme_set(theme_bw())

# had to download the land data from: https://www.naturalearthdata.com/downloads/
land <- 
  ne_load(scale = 50, type = "land", 
          category = "physical", 
          returnclass = "sf",
          destdir = "~/data/ne_50m_land/")

#==============
#  Eu-Afro basemap 
#==============
# generate baseline map
map <- 
  ggplot(data = land) +
  geom_sf(fill="grey", lwd = 0) + # lwd for borders
  coord_sf(xlim = c(-20, 57), ylim = c(22, 68), expand = FALSE) + # EU-North Africa 
  ggthemes::theme_map() # remove axis labels + grids

map

# add guiding dots to the place
# aboveground
map <- map + geom_point(data = pie.data.above.PieShown, aes(x = Longitude, y = Latitude), size = 1) 

# belowground
# map <- map + geom_point(data = pie.data.below.PieShown, aes(x = Longitude, y = Latitude), size = 1)

# add latitudes
map <- map + geom_hline(yintercept=c(30,40,50,60), linetype="dashed", alpha = 0.3)

# let's take a look
map

#==========================================
#  add pie charts to the basemap using subview
#==========================================
cols=c("#F7AC08","#27AAE1","#EF4136")

piesToBeVis <- pie.data.above.PieShown

for (i in 1:nrow(piesToBeVis)) {
  # subset the data for a pie
  pie_data<-piesToBeVis[i,c("Reference","Freq_pm","Freq_pp", "Freq_mm")]
  # gather data (previously melt)
  pie_data_melted <- melt(pie_data, id.vars = c("Reference"))
  
  # draw pie!
  pie <- 
    ggplot(pie_data_melted, aes(x="", value, fill = variable)) + 
    geom_bar(stat = "identity", color = "black", size = .1, width = .1) + # with border
    coord_polar(theta = "y") +
    scale_fill_manual(values=cols) +
    theme_void() +
    theme(legend.position="none") + 
    theme_transparent()
  
  # add the pie to map using geom_subview
  map <- map + geom_subview(subview=pie,
                            data = piesToBeVis[i,], 
                            aes(x = Longitude, y = Latitude), 
                            width=5, height=5)
}

# stop for a sec for R to process (sometimes the for loop takes a while)

### save pie plot!
ggsave(paste("~/your/fig/folder/",
             "Fig. 2a.pdf",
             sep = ""), 
       plot = map,
       width = 5,
       height = 4,
       units = "in"
)


#=====================================================================================
#  Fig. 2b Graph pipiens frequency vs latitudes
#=====================================================================================
### visualize all data points
# latitudinal gradient
latPipFreq <- 
  ggplot(data, aes(x=Freq_p, y=Latitude, color=Freq_p, shape=Microhabitat)) +
  geom_point(size=7.5) +
  scale_shape_manual(values=c(16, 17)) +
  ### outline North African samples ###
  geom_point(data=subset(data, data$Microhabitat == "aboveground" &
                           (data$Country == "Morocco" | data$Country == "Morocco" | data$Country == "Lybia" |
                              data$Country == "Algeria" | data$Country == "Tunisia")),
             size=7.5, fill=NA, color="black", shape=1) +
  geom_point(data=subset(data, data$Microhabitat == "belowground" &
                           (data$Country == "Morocco" | data$Country == "Morocco" | data$Country == "Lybia" |
                              data$Country == "Algeria" | data$Country == "Tunisia")),
             size=7.5, fill=NA, color="black", shape=2) +
  xlim(0,1)

### set a proper color scheme 
altcol=c("#EF4136","#F7AC08","#27AAE1")

latPipFreq <-
  latPipFreq + 
  scale_color_gradientn(colours = altcol,
                        values = scales::rescale(c(0, 0.5, 1))) +
  theme(
    # Hide panel borders and remove grid lines & background
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black")
  ) 
latPipFreq

### save latitude plot!
ggsave(paste("~/your/fig/folder/",
             "Fig. 2b.pdf",
             sep = ""), 
       plot = latPipFreq,
       width = 12,
       height = 7,
       units = "in"
)

#=====================================================================================
#  calculate correlation between latitudes and pipiens/molestus ancestry
#=====================================================================================
# just in case, subset data again
data.above <- subset(data, data$Microhabitat == "aboveground")
data.below <- subset(data, data$Microhabitat == "belowground")

# How many data points?
nrow(data)
nrow(data.above)
nrow(data.below)

# correlation 
library("ggpubr")

cor.test(data.above$Freq_p, data.above$Latitude, method=c("pearson"))
cor.test(data.below$Freq_p, data.below$Latitude, method=c("pearson"))

#=====================================================================================
#  4. visualize departure from expected allele frequencies
#=====================================================================================
#============================
# Fig. 2c. MAP with F_IS and HWE info
#============================
#===========================
#  Subset populations with N larger than a threshold and non-pooled
#===========================
F_IS.data.above.30 <- subset(data.above, 
                             data.above$N >= 30 &
                             data.above$MAF >= 0.05 &
                             data.above$PooledByAuthors == "no")
nrow(F_IS.data.above.30)

# =====
#  for the map, trim down populaitons by distance clustering
# =====
library(geosphere)
library(ggdendro)

pop.positions <- cbind(Longitude=F_IS.data.above.30$Longitude, Latitude=F_IS.data.above.30$Latitude)
dist <- distm(pop.positions,pop.positions, fun=distHaversine)/1000 # unit:km

rownames(dist) <- F_IS.data.above.30$Locality
colnames(dist) <- F_IS.data.above.30$Locality
head(dist)
dist

# Hierarchical Clustering with hclust
dist <- as.dist(dist)
hc <- hclust(dist)

## Plot the result
# ggdendrogram(hc, size = 5) +
#   geom_hline(yintercept=50, linetype="dashed", color = "red")

#### define clusters based on a distance (50km)
cluster50k <- cutree(hc, h = 50)
agg <- mutate(F_IS.data.above.30,
              cluster50k = as.factor(cluster50k))

### pick a population of largest N per distance group
agg.vis <-
  agg %>% group_by(cluster50k) %>% top_n(1, N)
# View(agg.vis)
nrow(agg.vis)

# =====
#  final F_IS map
# =====
theme_set(theme_bw())

# had to download the land data from: https://www.naturalearthdata.com/downloads/
land <- 
  ne_load(scale = 50, type = "land", 
          category = "physical", 
          returnclass = "sf",
          destdir = "~/data/ne_50m_land/")

# generate baseline map
basemap <- 
  ggplot(data = land) +
  geom_sf(fill="grey", lwd = 0) + # lwd for borders
  coord_sf(xlim = c(-20, 57), ylim = c(22, 68), expand = FALSE) + # EU-North Africa 
  ggthemes::theme_map() # remove axis labels + grids

# add dots with F_IS
map <- basemap + 
  geom_point(data = agg.vis, 
             aes(x = Longitude, y = Latitude, col=F_IS),
             # alpha=0.8,
             size = 10, shape=16) +
  # outline dots that violate HWE 
  geom_point(data=subset(agg.vis,
                         agg.vis$HWpval < 0.001),
             aes(x = Longitude, y = Latitude),
             size=10, fill=NA, color="red", shape=1) +
  scale_colour_gradient(high = "yellow", low = "black", na.value = NA)

# add latitudes
map <- map + geom_hline(yintercept=c(30,40,50,60), linetype="dashed", alpha = 0.3)

### save map with F_IS plot!
ggsave(paste("~/your/fig/folder/",
             "Fig. 2c.pdf",
             sep = ""), 
       plot = map,
       width = 8,
       height = 7,
       units = "in"
)

#============================
#  Fig. 2d. F_IS against latitudes
#============================
# latitudinal gradient
latF_IS <- 
  ggplot(F_IS.data.above.30,aes(x=F_IS, y=Latitude, color=F_IS)) + 
  geom_point(aes(x =F_IS, y = Latitude), size=10) +
  scale_colour_gradient(high = "yellow", low = "black", na.value = NA) +
  # outline dots that violate HWE 
  geom_point(data=subset(F_IS.data.above.30,F_IS.data.above.30$HWpval < 0.001),
             aes(x =F_IS, y = Latitude), size=10, 
             fill=NA, color="red", shape=1) +
  theme(
    # Hide panel borders and remove grid lines & background
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black")
  )
latF_IS

### save latitude plot!
ggsave(paste("~/your/fig/folder/",
             "Fig. 2d.pdf", 
             sep = ""), 
       plot = latF_IS,
       width = 8,
       height = 7,
       units = "in"
)

# ================================================================== #
# #################  bloodmeal analysis (Fig. 5a) ################## #
# ================================================================== #

#=====================================================================================
#  read bloodmeal data + metadata
#=====================================================================================
data <- read.csv("Table S2. Published Bloodmeal Data Related to Figure 5.csv", 
                 header=T,stringsAsFactors = FALSE)

#=====================================================================================
#  get basemap
#=====================================================================================
theme_set(theme_bw())

# had to download the land data from: https://www.naturalearthdata.com/downloads/
land <- ne_load(scale = 50, type = "land", category = "physical", returnclass = "sf",
                destdir = "~/data/ne_50m_land/")

# generate baseline map
map <- 
  ggplot(data = land) +
  geom_sf(fill="grey", lwd = 0) + # lwd for borders
  coord_sf(xlim = c(-20, 57), ylim = c(22, 68), expand = FALSE) + # EU-North Africa 
  ggthemes::theme_map() # remove axis labels + grids

# map

# add guiding dots to the place
map <- map + geom_point(data = pie.data, aes(x = Longitude, y = Latitude), size = 1)

# map

#=====================================================================================
#  add pie charts to the map using subview
#=====================================================================================
pie.data <- data

### add pie charts to the map (as subview)
### get pies for every row and add one pie at a time
#==============
#  pies for mammal vs birds
#==============
# set colors
color_mammal="#AA1F23"
color_bird="#65CBE4"
cols=c(color_mammal,color_bird)

for (i in 1:nrow(pie.data)) {
  # subset the data for a pie
  pie_data <- pie.data[i,c("Fig5A_PieID","Freq_MammalAll","Freq_Bird")]
  # gather data (previously melt)
  pie_data_melted <- melt(pie_data, id.vars = c("Fig5A_PieID"))
  
  # draw pie!
  pie <- 
    ggplot(pie_data_melted, aes(x="", value, fill = variable)) + 
    geom_bar(stat = "identity", color = "black", size = .3) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=cols) +
    theme_void() +
    theme(legend.position="none")
  
  # add the pie to map using geom_subview
  map <- map + geom_subview(subview=pie,
                            data = pie.data[i,], 
                            aes(x = Longitude, y = Latitude), 
                            width=5, height=5)
}
# stop for a sec for R to process (sometimes the for loop takes a while)
map



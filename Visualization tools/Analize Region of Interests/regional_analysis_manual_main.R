# Title:        regional_analysis_main.R
#
# Description:  Main script describing a step-wise analysis of a dataset
#                 with points including x and y-coordinates and a class 
#                 
# Usage:        Make sure you have installed all the required packages
#                 (try running to see which ones need installation)
#               Step through all the steps by running them (select 
#                 parts and ctrl+shift+enter to run snippets)
#
# Version:      2014-07-23 10:30:47 CEST
# Author:       Olle Nordesjö (olle.nordesjo[at]gmail.com)

#### Go to top #####
#
# Set the working directory (for loading and saving data)
setwd("~/Work/Projects/in-situ-sequencing-regional-analysis-tools/")
source("regional_analysis_source.R")
# Load the data
inputName<-"regional_analysis_sample_data_165"
data.raw<-read.csv(paste("input/",inputName,".csv",sep=""))

#### -------- subsetting (Selection & Preprocessing) ####
# If you want to remove some of the classes (transcripts), it's possible to do 
#   it here. Either use only the most abundant transcripts or select certain 
#   ones to keep or remove.

# Keep only the 20 most abundant transcripts and and remove NNNN.
#   %in% tests for membership in a group. ("1 %in% c(1,2)" gives TRUE)
includeTheseTranscripts<-rev(names(sort(table(data.raw$name))))[1:20]
data<-data.raw[(!data.raw$name %in% c("NNNN"))&
                 data.raw$name %in% includeTheseTranscripts,]

# Remove factors that are not used
data$name<-droplevels(data$name)

#### -------- polygon selection (Selection % Preprocessing) ####
# There are different alternatives for defining regions:
#
# 1. Select regions manually
# 2. Load predefined regions from a file
#
# Run what's in the section below after reading the how-to:
#
# How-to:
# 
# Click on the plot to select regions around some points. To close the 
#   polygon and continue with the next one, right-click and 
#   press "stop". Continue with any number of polygons until you are 
#   happy with the result. 
# Press escape when you've drawn all polygons. 
# Closing the window before pressing escape may make RStudio crash. 
# Now run the following lines of code:
# --- Select regions manually ####
select=F
restart=T
if(select==T){
  try(polygonsLast<-polygons,silent=T)
  if(restart==T){
    polygons<-c()}
  while(T){
    polygons<-addPolygon(data,polygons)
  }
  dev.off()
  polygons.df<-lapply(FUN=as.data.frame,polygons)  # convert to df for right format
  polygons.melted<-melt(polygons.df,id.vars=c('x','y'))
}

# A polygon is stored as a series of x and y points, and the polygons are 
#   stored in a list. Some actions for manipulating the polygons:
#     polygons[[3]]      - prints the points in polygon 3
#     polygons[[2]]<-c() - remove polygon 2 (polygon 3 will receive index 2)

### -------- 2. Load predefined regions from a file ####
# An alternative is to import polygons with these lines
import=F
if(import==T){
  polygons<-c()
  importedPolygons<-read.csv("input/regional_analysis_polygons_sample.csv")
  for(i in 1:max(importedPolygons$L1)){
    polygons[[i]]<-importedPolygons[which(importedPolygons$L1==i),1:2]
  }
  polygons.df<-lapply(FUN=as.data.frame,polygons)  # convert to df for right format
  polygons.melted<-melt(polygons.df,id.vars=c('x','y'))
}

### -------- save polygons to file ####
# if you want to save your polygons (for reproducibility),
# you can do it here. The polygons are reshaped in what is called a "melted"
# form. Type "polygons.melted" to see what it looks like.
polygons.df<-lapply(FUN=as.data.frame,polygons)  # convert to df for right format
polygons.melted<-melt(polygons.df,id.vars=c('x','y'))
write.table(file="input/regional_analysis_polygons_temp.csv",polygons.melted,sep=",",quote=F)

### -------- calculate blob-in-polygons (Transformation)####
# Calculate the information about which polygons the blobs
# are inside. If polygons are overlapping, the last selected polygon will
# be used. All points outside of any region will get region=0.

data$region<-0
for (i in 1:length(polygons.df)){
  ids<-pnt.in.poly(pnts=data[c("global_X_pos","global_Y_pos")],
                   poly.pnts=polygons.df[[i]])
  data[which(ids$pip==1),]$region<-i
}

# Remove the data that is not included in the polygons for accurate statistics
data.inregions<-data[data$region!=0,]
relativeTranscriptFrequencies<-hash(names(table(data.inregions$name)),
                      as.numeric(table(data.inregions$name))/sum(as.numeric(table(data.inregions$name))))
data.inregions.dt<-data.table(data.inregions) # Create a "data table" for quick calculations

#### --- statistical analysis (Data mining) ####
# From here, it is possible to use the values to do different things to
#   look for patterns in the expression profiles
#
# These are the implemented measures so far:
#   1. frequency         - the number of transcripts per region per class
#                   5 means 5 transcripts of some class in the region
#   2. relativeFrequency - frequency of one transcript divided by sum of 
#                             transcripts in the region. 
#                   0.01 means 1% of transcripts of some class in the region
#   3. deviation         - difference between relative frequency in a region 
#                             compared to all regions (presented as percentage points). 
#                  -0.01 means the frequency is 1 percentage point less 
#                   in the specified region
#                             

# Choose and run one of these below

#regionwiseExpression<-calculateRegionwiseExpression(data.inregions.dt,frequency)         #1
regionwiseExpression<-calculateRegionwiseExpression(data.inregions.dt,relativeFrequency)  #2
#regionwiseExpression<-calculateRegionwiseExpression(data.inregions.dt,deviation)         #3

#### --- kmeans (Data mining) ####
# kmeans can now be used to compare all the different regions to look 
# for common patterns. Choose a high number of clusters first, and make sure
# that "regionwiseExpression" doesn't contain information about regions or 
# clusters yet (regionwiseExpression$region shouldn't exist for example). 
# if it does, rerun from last step

if(length(regionwiseExpression[,1])>2){
  kmeansData<-calculateKmeans(regionwiseExpression,nClusters=3)
} else {
  cat("You have less than three regions, choosing two clusters")
  kmeansData<-calculateKmeans(regionwiseExpression,nClusters=2)
}

# After having run the kmeans, the clusters has to be associated with the
# regions to map them to each blob
regionwiseExpression$region<-seq(length(regionwiseExpression[,1]))
regionwiseExpression$cluster<-factor(c(kmeansData$cluster))
# Create a hash with regions as keys and clusters as values
hash.region.cluster<-hash(regionwiseExpression[,c("region")],
                          regionwiseExpression[,c("cluster")])
# ... and map the clusters to the blobs
data$cluster<-factor(unlist(lapply(data$region,FUN=getClusterID)))
data$region<-factor(data$region)



#### --- visualization (Data presentation) ####
# Some different examples for plotting the patterns that appear

# first a plot to show the regions qualitatively
plotRegions(data,polygons.melted,color="region")
ggsave(filename=paste("output/",inputName,"_regions_qualitatively.png",sep = ""),width=20,height=20)

# show the different clusters that have been generated.
plotRegions(data,polygons.melted,color="cluster")
ggsave(filename=paste("output/",inputName,"_clusters_qualitatively.png",sep = ""),width=20,height=20)

# parallell coordinates to display any obvious patterns from the clustering
# Each line represents a region, and the lines are colored by what cluster
# they are assigned to
ggparcoord(regionwiseExpression,columns=1:19,groupColumn="cluster",
           alphaLines=0.4,showPoints=F,scale='globalminmax')+
  guides(colour=guide_legend(override.aes=list(alpha = 1, size= 3)))+
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(filename=paste("output/",inputName,"_parallell_coordinates_relativeFrequency.png",sep = ""))

# Visualize all the different regions, and their individual expression profiles
regionwiseExpression.melted<-melt(regionwiseExpression,id.vars=c("cluster","region"))
ggplot(regionwiseExpression.melted,aes(variable,value,fill=cluster))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~region)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10))+
  scale_fill_manual(values=scale_fill_brewer(type="qual",palette="Set3")$palette(6)[2:5])
ggsave(filename=paste("output/",inputName,"_expression_profile_per_region.png",sep = ""),width=20,height=20)


# Visualize the different expression profiles for the different clusters
ggplot(regionwiseExpression.melted,aes(variable,value,color=cluster))+
  geom_point(stat="identity",position="dodge",size=3,alpha=0.6)+facet_wrap(~cluster)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10))+
  scale_color_manual(values=scale_fill_brewer(type="qual",palette="Set3")$palette(6)[2:5])+
  guides(colour=guide_legend(override.aes=list(alpha = 1, size= 5)))
ggsave(filename=paste("output/",inputName,"_expression_profile_per_cluster.png",sep = ""))

# Visualize a hierarchical clustering of the center vectors from kmeans clustering
#  This may contain information about if different expression profiles exists, but
#  should be interpreted only in association with all other data.
png(filename=paste("output/",inputName,"_dendrogram_heatmap_of_clustering_vectors.png",sep = ""))
heatmap.2(as.matrix(kmeansData$center),trace='none')
dev.off()

plotXY(data[!data$name%in%c(""),],facet=T,size=1,ncol = 2)
ggsave(filename=paste("output/",inputName,"_spatial_distribution_facet.png",sep = ""),width=14,height=14)

plotXY(data[!data$name%in%c(""),],facet=F,size=1,ncol = 2)
ggsave(filename=paste("output/",inputName,"_spatial_distribution.png",sep = ""),width=14,height=14)


# For visualization of two transcripts separately
#plotXY(data[data$name%in%c("ACTB","EZH2/MKI67"),],facet=F,size=3)
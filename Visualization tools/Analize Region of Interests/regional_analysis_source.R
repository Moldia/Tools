# Title:        regional_analysis_source.R
#
# Description:  Contains all source code for the functions 
#                 used in the regional analysis
#
# Usage:        run source("regional_analysis_source.R") once or 
#                 press ctrl+shift+enter to load all the functions
#                 After this, open regional_analysis_main.R for further 
#                 instructions
#
# Version:      2014-07-17 16:40:47 CEST
# Author:       Olle Nordesjö (olle.nordesjo[at]gmail.com)

#--- Set up the required packages  --------------
require(gplots)
require(SDMTools)
require(GGally)
require(ggplot2)
require(deldir)
require(data.table)
require(reshape)
require(gplots)
require(hash)

selectPoints<-function(data){
  plot.new()
  x11()
  plot(data$global_X_pos,data$global_Y_pos,pch='.')
  pos<-locator(type='p')
  return(pos)
}

reorderPolygonPoints<-function(uniquePoints){
  normalizedCoords<-NULL
  reorderedCoords<-NULL
  angles<-list(NULL)
  # Reorder a list containing polygon points so that they are plotted in positive direction
  nPolygons=length(uniquePoints)
  for (i in 1:nPolygons){  
    nPoints<-length(uniquePoints[[i]][,1])
    centerCoords<-colSums(uniquePoints[[i]])/nPoints
    normalizedCoords[[i]]<-uniquePoints[[i]]-matrix(centerCoords,byrow=T,nrow=nPoints,ncol=2)
    if(length(normalizedCoords[[i]])==0){
      next
    }
    angles[[i]]<-0
    for (j in 1:length(normalizedCoords[[i]][,1])){
      x=normalizedCoords[[i]][j,1]
      y=normalizedCoords[[i]][j,2]
      angles[[i]][j]=360*Arg(x+complex(imaginary=y)/(2*pi))+180
    }
    reorderedCoords[[i]]<-uniquePoints[[i]][order(angles[[i]]),]
  }
  return(reorderedCoords)
}

addPolygon<-function(data,polygons=polygons){
  i<-length(polygons)
  while(length(dev.list()>0)){dev.off()}
  plot.new()
  x11()
  plot(data$global_X_pos,data$global_Y_pos,pch='.')
  wID<-dev.cur()
  dev.set(wID)
  
  lapply(polygons,polygon)
  
  polygons[[i+1]]<-locator(type='o')
  
  return(polygons)
}  




deviation<-function(id){
  dt.subset<-data.inregions.dt[J(id)]
  a<-as.numeric(table(dt.subset$name))
  a=a/sum(a)
  a=a-as.numeric(values(transcriptFrequencies))
}

relativeFrequency<-function(id){
  dt.subset<-data.inregions.dt[J(id)]
  a<-as.numeric(table(dt.subset$name))
  a=a/sum(a)
}

frequency<-function(id){
  dt.subset<-data.inregions.dt[J(id)]
  a<-as.numeric(table(dt.subset$name))
}

getClusterID<-function(x){
  if(x!=0){
    hash.region.cluster[[as.character(x)]]
  } else {return(0
  )}
}

calculateRegionwiseExpression<-function(dt,method=absolutedensity){
  setkey(dt,region)
  arguments=seq(min(dt$region),max(dt$region))
  regionwise.raw<-lapply(FUN=method,X=arguments)
  regionwiseExpression<-data.frame(matrix(unlist(regionwise.raw),ncol=length(levels(data$name)),byrow=T))
  names(regionwiseExpression)<-levels(dt$name)
  return(regionwiseExpression)
}

calculateKmeans<-function(regionwiseExpression,nClusters){
  ns<-names(head(regionwiseExpression))
  km<-kmeans(regionwiseExpression,centers=nClusters,iter.max=400,nstart=10)
  return(km)  
}

plotRegions<-function(data,polygons,color="cluster"){
  #   Plots the regions that are defined by polygons
  #   
  #   Args: 
  #     data:        The dataset with the points
  #     polygons:    The data frame with polygons. 
  #     whatToColor: Either regions or clusters, whatever is interesting
  #
  #   Returns: 
  #     A plot with numbers of the regions overlayed by a polygon
  
  xcent=tapply(polygons$x,polygons$L1,mean) # finding the center of a polygon
  ycent=tapply(polygons$y,polygons$L1,mean)
  nRegions<-diff(range(as.numeric(data$region)))+1
  
  # has to be a data frame for ggplot
  centers<-as.data.frame(cbind(xcent,ycent,region=seq(1,nRegions-1)))
  
  j<-ggplot(data,aes(global_X_pos,global_Y_pos))+
    geom_polygon(data=polygons,aes(x,y,group=L1),size=1,alpha=0.4)+ # the polygons
    geom_point(data=data,aes_string(color=color),size=1)+ # the blobs
    geom_text(data=centers,aes(xcent,ycent,label=region),size=3)+# polygon text
    guides(colour = F)+ # remove legend
    coord_fixed() # fix aspect ratio
  if(color=="cluster"){
    j<-j+guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))+
      scale_color_brewer(type="qual",palette="Set3")
  } else {
    set.seed(1)
    j<-j+scale_color_manual(values=sample(terrain.colors(nRegions)))
  }
  return(j)
}

plotXY<-function(data,size=1,alpha=0.9,facet=T,ncol=1){
  #   Produces a faceted plot showing the locations of the transcripts
  #   
  #   Args: 
  #     data:   a dataframe
  #     size:   the size of the points
  #     alpha:  the opacity of the dots (0 to 1). Lower if data frame is dense
  #     facet:  whether to facet on transcript classes. Good for dense data
  #     ncol:   the number of columns to use in the legend. Good for many transcripts
  #
  #   Returns: 
  #     A transcript-faceted plot of x and y position
  
  plotObject<-ggplot(data,aes(global_X_pos,global_Y_pos))+
    geom_point(size=size,alpha=alpha,aes(color=factor(name)))+
    guides(color=guide_legend(ncol=ncol,override.aes = list(alpha=1)))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  if(facet==T){
    plotObject<-plotObject+facet_wrap(~name)
  }
  plotObject
  #ggsave(filename="facetPlot.png",plot=plotObject,width=12,height=12)  
}

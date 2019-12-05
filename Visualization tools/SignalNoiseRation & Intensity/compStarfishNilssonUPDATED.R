setwd("/Volumes/Christoffer")


library(data.table)
library(spatstat)
library(entropy)
library(RandomFields)
library(INLA)
library(fields)
library(sp)
library(methods)
library(iterators)
library(maptools)
library(rgeos)
library(tcltk)
library(FNN)
library(ggplot2)
library(LaplacesDemon)
library(MASS)
library(plotly)


nilsson_sorted = read.csv(file = "nilsson_export.csv", header = TRUE, sep = ",")
starfish_sorted = read.csv(file = "starfish_export.csv", header = TRUE, sep = ",")

nilsson_x = nilsson_sorted$nilsson_x_sorted_rounded;
nilsson_y = nilsson_sorted$nilsson_y_sorted_rounded;
nilsson_g = nilsson_sorted$Var1;

starfish_x = starfish_sorted$starfish_x_sorted_normalized_shift_rounded;
starfish_y = starfish_sorted$starfish_y_sorted_normalized_shift_rounded;
starfish_g = starfish_sorted$Var1;

nilssonDataTable = data.table(nilsson_sorted)
starfishDataTable = data.table(starfish_sorted)
nilssonDataFrame = data.frame(nilsson_sorted) 
starfishDataFrame = data.frame(starfish_sorted) 

nilssonOccurances = data.table(table(unlist(nilssonDataTable$Var1)) )
starfishOccurances = data.table(table(unlist(starfishDataTable$Var1)) )
mergedOccurances = merge(nilssonOccurances, starfishOccurances, by = "V1")
data.frame(mergedOccurances$V1)
correlation = cor(mergedOccurances$N.x, mergedOccurances$N.y)
correlation



n <- 32


ggplot(mergedOccurances, aes(x = N.x, y = N.y)) +
  geom_point(size = 10, aes(x = N.x, y = N.y, color = V1), show.legend = FALSE) +  
  geom_text(label = mergedOccurances$V1, size = 10, angle = 40, vjust = 0, nudge_y = 0.05) +
  ggtitle("Correlation") +
  xlab("Nilsson pipeline") + 
  ylab("Starfish pipeline") +
  theme(plot.title = element_text(hjust = 0.5)) +
  #theme(legend.position="bottom") + 
  theme(text = element_text(size=30)) + 
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10') +
  geom_smooth(method=lm, se=FALSE) 


range(nilsson_x) # determining the range to be specified in owin
range(nilsson_y)
nilsson_ppp = ppp(nilsson_x, nilsson_y, owin(xrange=c(800, 38500), yrange=c(800, 52700)))
#plot(nilsson_ppp_split, cex = 0.5, cex.main=0.5)

range(starfish_x) # determining the range to be specified in owin
range(starfish_y) 
starfish_ppp = ppp(starfish_x, starfish_y, owin(xrange=c(200, 39000), yrange=c(100, 53800)))

nilsson_ppp_split = split.ppp(nilsson_ppp, f=nilsson_sorted$Var1) 
starfish_ppp_split = split.ppp(starfish_ppp, f = starfish_sorted$Var1)
#plot(starfish_ppp_split, cex = 0.5, cex.main=0.5)

h_nilsson = hyperframe(pointPattern = nilsson_ppp_split) # the downside is that we loose the lables of the genes 
h_starfish = hyperframe(pointPattern = starfish_ppp_split) # the downside is that we loose the labels of the genes 

nilsson_DM = with(h_nilsson, data.matrix(data.frame(as.array(pointPattern))))
starfish_DM = with(h_starfish, data.matrix(data.frame(as.array(pointPattern))))

``
resultsmatrix = data.frame(mergedOccurances$V1, divergencevector)
resultsmatrix$quantdiff = (mergedOccurances$N.x-mergedOccurances$N.y)/(mean(mergedOccurances$N.x+mergedOccurances$N.y))

colnames(resultsmatrix)
names(resultsmatrix)[names(resultsmatrix) == "mergedOccurances.V1"] <- "Genes"
colnames(resultsmatrix)

dev.off()


OutputCompPlot = ggplot() + geom_point(data = resultsmatrix, size = 2, aes(x = divergencevector, y = quantdiff, color = Genes)) +
  ggtitle("KLD and difference in number of gene reads")+xlab("KL divergence")+ ylab("Difference in number of reads (Normalized)") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom") + theme(text = element_text(size=11))
OutputCompPlot


{plot(mergedOccurances$N.x, mergedOccurances$N.y, log = "xy", pch = c(1:8), col = color, xlab = "Nilsson pipeline", ylab = "Starfish pipeline", cex = 1)
  title(main = "Nilsson vs Starfish Pipeline")}
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("right", inset=c(-0.2,0), legend = mergedOccurances$V1, pch=c(1:8), col = color, title="Group", cex=0.7, horiz=FALSE)

# SLC17A7
plot(nilsson_ppp_split$SLC17A7$x, nilsson_ppp_split$SLC17A7$y, col = "green", xlab = "x values", ylab = "y values")
points(starfish_ppp_split$SLC17A7$x, starfish_ppp_split$SLC17A7$y, col = "red")
title(main = "SLC17A7", sub = "High KLD",
      cex.main = 3,   font.main= 4, col.main= "black",
      cex.sub = 1.2, font.sub = 4, col.sub = "black",
      col.lab ="darkblue"
)
legend(27500, 50000, legend=c("Nilsson", "Starfish"),
       col=c("green", "red"), pch=1:1, cex=1)


# ETNPPL
plot(nilsson_ppp_split$ETNPPL$x, nilsson_ppp_split$ETNPPL$y, col = "green", xlab = "x values", ylab = "y values")
points(starfish_ppp_split$ETNPPL$x, starfish_ppp_split$ETNPPL$y, col = "red")
title(main = "ETNPPL", sub = "High difference in number of gene reads",
      cex.main = 3,   font.main= 4, col.main= "black",
      cex.sub = 1.2, font.sub = 4, col.sub = "black",
      col.lab ="darkblue"
)
legend(32000, 50000, legend=c("Nilsson", "Starfish"),
       col=c("green", "red"), pch=1:1, cex=1)


# LY86
plot(nilsson_ppp_split$LY86$x, nilsson_ppp_split$LY86$y, col = "green", xlab = "x values", ylab = "y values")
points(starfish_ppp_split$LY86$x, starfish_ppp_split$LY86$y, col = "red")
title(main = "LY86", sub = "Neither high KDL nor high difference in number of gene reads",
      cex.main = 3,   font.main= 4, col.main= "black",
      cex.sub = 1.2, font.sub = 4, col.sub = "black",
      col.lab ="darkblue"
)
legend(27000, 49000, legend=c("Nilsson", "Starfish"),
       col=c("green", "red"), pch=1:1, cex=1)

# NOSTRIN
plot(nilsson_ppp_split$NOSTRIN$x, nilsson_ppp_split$NOSTRIN$y, col = "green", xlab = "x values", ylab = "y values")
points(starfish_ppp_split$NOSTRIN$x, starfish_ppp_split$NOSTRIN$y, col = "red")
title(main = "NOSTRIN", sub = "Neither high KDL nor high difference in number of gene reads",
      cex.main = 3,   font.main= 4, col.main= "black",
      cex.sub = 1.2, font.sub = 4, col.sub = "black",
      col.lab ="darkblue"
)
legend(30000, 49000, legend=c("Nilsson", "Starfish"),
       col=c("green", "red"), pch=1:1, cex=1)

# LAMP5
plot(nilsson_ppp_split$LAMP5$x, nilsson_ppp_split$LAMP5$y, col = "green", xlab = "x values", ylab = "y values")
points(starfish_ppp_split$LAMP5$x, starfish_ppp_split$LAMP5$y, col = "red")
title(main = "LAMP5", sub = "Neither high KDL nor high difference in number of gene reads",
      cex.main = 3,   font.main= 4, col.main= "black",
      cex.sub = 1.2, font.sub = 4, col.sub = "black",
      col.lab ="darkblue"
)
legend(27000, 49000, legend=c("Nilsson", "Starfish"),
       col=c("green", "red"), pch=1:1, cex=1)

# DCN
plot(nilsson_ppp_split$DCN$x, nilsson_ppp_split$DCN$y, col = "green", xlab = "x values", ylab = "y values")
points(starfish_ppp_split$DCN$x, starfish_ppp_split$DCN$y, col = "red")
title(main = "DCN", sub = "Neither high KDL nor high difference in number of gene reads",
      cex.main = 3,   font.main= 4, col.main= "black",
      cex.sub = 1.2, font.sub = 4, col.sub = "black",
      col.lab ="darkblue"
)
legend(27000, 49000, legend=c("Nilsson", "Starfish"),
       col=c("green", "red"), pch=1:1, cex=1)

# LINC00507
plot(nilsson_ppp_split$LINC00507$x, nilsson_ppp_split$LINC00507$y, col = "green", xlab = "x values", ylab = "y values")
points(starfish_ppp_split$LINC00507$x, starfish_ppp_split$LINC00507$y, col = "red")
title(main = "LINC00507", sub = "Relatively high KLD and small difference in number of gene read differences",
      cex.main = 3,   font.main= 4, col.main= "black",
      cex.sub = 1.2, font.sub = 4, col.sub = "black",
      col.lab ="darkblue"
)
legend(27000, 49000, legend=c("Nilsson", "Starfish"),
       col=c("green", "red"), pch=1:1, cex=1)





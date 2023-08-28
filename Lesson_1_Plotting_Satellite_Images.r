# LESSON 1: PLOTTING SATELLITE IMAGES

# Install the raster package
install.packages("raster")
library(raster)

# Setting the working directory
setwd("C:/Data_telerilevamento/lab")

# Import the data with the function brick()
img_2011 <- brick("p224r63_2011_masked.grd")

# Visualize information
img_2011

# Plot the data
plot(img_2011)

# Create a colour palette
cl_0 <- colorRampPalette(c("darkgoldenrod1","chartreuse3","darkcyan"))(100) 
# 100: n. of shades
cl_1 <- colorRampPalette(c("darkgoldenrod1","chartreuse3","darkcyan"))(3) 
# 3: n. of shades

# Plotting the data with cl_0
plot(img_2011, col=cl_0)
plot (img_2011, col=cl_1) # Less detail

# Plotting one element
plot(img_2011[[4]], col=cl_0)
# or
plot(img_2011$B4_sre, col=cl_0)
# or
nir <- img_2011[[4]]
plot(nir, col=cl_0)

# Change the colour palette
cl_2 <- colorRampPalette(c("aquamarine","cadetblue","darkblue"))(100)
plot(nir, col=cl_2)

# EXERCISE: plot the NIR band
# B1 = blue
# B2 = green
# B3 = red
# B4 = NIR

cl_3 <- colorRampPalette(c("pink","violet","darkorchid4"))(100)
plot(nir, col=cl_3)
dev.off # Clean plots

# Export graphs (PDF)
pdf("Graph_1.pdf")
plot(nir, col=cl_3)
dev.off()

# Export graphs (PNG)
png("Graph_1.png")
plot(nir, col=cl_3)
dev.off()

# Plot 2 or more graphs in a multiframe
par(mfrow=c(2,1)) # Multiframe with 2 rows and 1 column
plot(img_2011[[3]], col=cl_3)
plot(img_2011[[4]], col=cl_3)
devo.off()

# With a multiframe we can see the differences between color palette
par(mfrow=c(2,1))
plot(img_2011[[4]], col=cl_0) # n. of shadows = 100 (more detailed)
plot(img_2011[[4]], col=cl_1) # n. of shadows = 3 (less detailed)
dev.off

# Plot all bands/layers (PDF)
pdf("Grapf_2.pdf")

par(mfrow=c(2,2))

# Blue band
cl_blue <- colorRampPalette(c("blue4", "blue2", "lightblue"))(100)
plot(img_2011[[1]], col=cl_blue)

# Green band
cl_green <- colorRampPalette(c("chartreuse","chartreuse4","darkgreen"))(100)
plot(img_2011[[2]], col=cl_green)

# Red band
cl_red <- colorRampPalette(c("darkorange","brown3","darkred"))(100)
plot(img_2011[[3]], col=cl_red)

# NIR band
cl_NIR <- colorRampPalette(c("pink","violet","darkorchid4"))(100)
plot(img_2011[[4]], col=cl_NIR)

dev.off()

# Plot all bands/layers (PNG)
png("Grapf_2.png")

par(mfrow=c(2,2))

# Blue band
cl_blue <- colorRampPalette(c("blue4", "blue2", "lightblue"))(100)
plot(img_2011[[1]], col=cl_blue)

# Green band
cl_green <- colorRampPalette(c("chartreuse","chartreuse4","darkgreen"))(100)
plot(img_2011[[2]], col=cl_green)

# Red band
cl_red <- colorRampPalette(c("darkorange","brown3","darkred"))(100)
plot(img_2011[[3]], col=cl_red)

# NIR band
cl_NIR <- colorRampPalette(c("pink","violet","darkorchid4"))(100)
plot(img_2011[[4]], col=cl_NIR)

dev.off()

# RGB plotting, plot of multi-layered raster object
# 3 bands are combined such that they represent the red, green and blue channel
# This function can be used to make 'true (or false) color images
# from Landsat and other multi-band satellite images.

plotRGB(img_2011, r=3, g=2, b=1, stretch="Lin") # Natural colors
plotRGB(img_2011, r=4, g=2, b=1, stretch="Lin") # Vegetation = red
plotRGB(img_2011, r=3, g=4, b=1, stretch="Lin") # Vegetation = green
plotRGB(img_2011, r=3, g=2, b=4, stretch="Lin") # Vegetation = blue

# Multiframe with natural and false colours
par(mfrow=c(2,1))
plotRGB(img_2011, r=3, g=2, b=1, stretch="Lin") # Natural colors
plotRGB(img_2011, r=4, g=2, b=1, stretch="Lin") # Vegetation = red
dev.off()

# Histogram stretching
par(mfrow=c(2,1))
plotRGB(img_2011, r=3, g=2, b=1, stretch="Hist") # Natural colors
plotRGB(img_2011, r=4, g=2, b=1, stretch="Hist") # Vegetation = red
dev.off()

# Differences between the 2 type of stretch
par(mfrow=c(2,1))
plotRGB(img_2011, r=3, g=2, b=1, stretch="Lin") # Natural colors (Lin)
plotRGB(img_2011, r=3, g=2, b=1, stretch="Hist") # Natural colors (Hist)
dev.off()





# EXERCISE

# Plot the NIR band
cl_NIR <- colorRampPalette(c("pink","violet","darkorchid4"))(100)
plot(img_2011[[4]], col=cl_NIR)

# Import the 1988 image
img_1988 <- brick("p224r63_1988_masked.grd")
plot(img_1988)

# Plot in RGB space (natural colors)
plotRGB(img_1988, 4, 3, 2, stretch="Lin") # faster way (Lin)
plotRGB(img_1988, 4, 3, 2, stretch="Hist") # faster way (Hist)

# Multiframe to see the differences between 1988 and 2011
par(mfrow=c(2,1))
plotRGB(img_1988, 4, 3, 2, stretch="Lin")
plotRGB(img_2011, 4, 3, 2, stretch="Lin")
dev.off()

# Multiframe with 4 images
par(mfrow=c(2,2))
plotRGB(img_1988, 4, 3, 2, stretch="Lin")
plotRGB(img_2011, 4, 3, 2, stretch="Lin")
plotRGB(img_1988, 4, 3, 2, stretch="Hist")
plotRGB(img_2011, 4, 3, 2, stretch="Hist")
dev.off()

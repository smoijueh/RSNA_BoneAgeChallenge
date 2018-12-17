## Samuel Moijueh
##
## Preprocesses the raw image files of the RSNA Boneage dataset
## The normalized images are written to file
## For more details see https://smoijueh-boneage.netlify.com/
#############################################################################

# Load Packages
library(imager)
library(dplyr)
library(rsdepth)
library(EBImage)
library(dbscan)
library(fpc)

options("max.contour.segments"= 200000)


################ FUNCTIONS ################
# resize and pad the image to 500x500 resolution
squareImage<-function(img, desired_sq_size = 500){
  # get the intensity of the image background
  d <- as.data.frame(img)
  idx<-which(d$x == max(d$x) & d$y == min(d$y))
  bkg_intensity <- 0
  old_size <- dim(img)[1:2]
  resized_img <- NULL

  if (old_size[1] > old_size[2]){
    resized_img <- resize(img, w = desired_sq_size)
  } else {
    resized_img <- resize(img, h = desired_sq_size)
  }

  new_size <- dim(resized_img)[1:2]

  x_pad_adjust = (desired_sq_size - new_size[1])
  y_pad_adjust = (desired_sq_size - new_size[2])

  padded_x_img <- pad(resized_img, x_pad_adjust, "x", val=rep(bkg_intensity, spectrum(resized_img)))
  padded_final <- pad(padded_x_img, y_pad_adjust, "y",val=rep(bkg_intensity, spectrum(resized_img)))
  return(padded_final)
}

# Sobel Edge Detection using an Illumination-Gradient based Algorithm
detect.edges <- function(im,sigma=1){
  isoblur(im,sigma) %>% imgradient("xy") %>% enorm %>% imsplit("c") %>% add
}

# edge detection, extract hand contours, map onto Cartesian xy plane
segmentHand<-function(im){

  # apply a gradient transformation and detect edges on high contrast images
  if (dark_image(im) || RMS_contrast(im) > 0.10){
    grad.sq <- imgradient(im,"xy") %>% map_il(~ .^2) %>% add()
    edges <- imsplit(grad.sq,"c") %>% add
    edges <- detect.edges(im,2) %>% sqrt
    im <- edges
  }

  # fit a linear model to remove trend in illumination, threshold
  d <- as.data.frame(im)
  px <- d %>% lm(value ~ x*y,data=.) %>% predict(d) %>% { im - . } %>% threshold

  # Clean up
  px <- clean(px,3) %>% imager::fill(7)
  iso<-contours(px)

  # get the top two most densely populated contours
  max.c <- which.max(lapply(iso, function(x) sum(lengths(x))))
  max.c2 <- order(unlist(lapply(iso, function(x) sum(lengths(x)))), decreasing = T)[2]

  xval=c(iso[[max.c]]$x,iso[[max.c2]]$x)
  yval=c(iso[[max.c]]$y,iso[[max.c2]]$y)
  contour.df<-data.frame(x=xval,y=yval)

  crop_coordinates <- getHandContour(contour.df)
  return(crop_coordinates)
}

# segment the hand, apply DBSCAN algorithm
getHandContour<-function(contour.df){
  wrist_coordinate <- find_wrist(contour.df)
  wrist_coordinate <- ifelse(wrist_coordinate %>% is.null %>% rep(2), c(1,1), wrist_coordinate)

  if (max(contour.df$y)-wrist_coordinate[2] > 300){
    new_contour.df <- contour.df[contour.df[,2] < max(contour.df[,2]) - 200,]
  } else{
    new_contour.df<-contour.df[contour.df[,2] < wrist_coordinate[2] - 30,]
  }

  # Compute DBSCAN, Density-based spatial clustering
  db <- dbscan::dbscan(new_contour.df, eps = 3, minPts = 5)

  # extract the cluster information and join to the original dataframe
  new_contour.df$cluster <- db$cluster
  num_clusters = unique(new_contour.df$cluster) %>% length()

  if (num_clusters <= 2){
    if (diff(range(sort(table(new_contour.df$cluster),decreasing=T))) > 3000){
      cluster_of_interest <- names(sort(table(new_contour.df$cluster),decreasing=T)[1])
    } else {
      cluster_of_interest <- nearest_cluster(new_contour.df)
    }

    if(cluster_of_interest >= 1 &&
       diff(range(new_contour.df[new_contour.df$cluster == cluster_of_interest,]$y)) > 1000){
      new_contour.df <- new_contour.df[new_contour.df$cluster == cluster_of_interest,]
    }

    new_contour.df <- new_contour.df[,1:2]

    min.x<-min(new_contour.df$x); max.x<- max(new_contour.df$x);
    min.y <- min(new_contour.df$y); max.y <- max(new_contour.df$y)
    crop_coordinates <- c(min.x,max.x,min.y,max.y)
  } else {
    crop_coordinates <- analyze_DBSCAN_clusters(new_contour.df)
  }
  return(crop_coordinates)
}

# helper function
nearest_cluster<-function(df){
  centroid<-as.integer(centroid(as.matrix(df[1:2])))
  min_diff = c(Inf,Inf)
  for (clus in unique(df$cluster)){
    centroid_cluster_df <- df[df$cluster == clus,][1:2]
    xpts <- centroid_cluster_df[centroid_cluster_df$y == centroid[2],]$x
    if(length(xpts) == 0) return(-1)
    closeness <- c(diff(range(xpts)), clus)
    if (closeness[1] < min_diff[1]) min_diff = closeness
  }
  return(min_diff[2])
}

# obtain the coordinates of the hand cluster
analyze_DBSCAN_clusters<-function(df){

  # find the centroid of the hand cluster
  centroid<-as.integer(centroid(as.matrix(df[1:2])))

  xframe = range(df$x)
  yframe = range(df$y)

  # if the cluster is out of frame
  if (dplyr::between(centroid[1], xframe[1], xframe[2]) == 0 ||
      dplyr::between(centroid[2], yframe[1], yframe[2]) == 0){
    min.x <- min(df$x);  max.x <- max(df$x)
    min.y <- min(df$y);  max.y <- max(df$y)
    crop_coordinates <- c(min.x,max.x,min.y,max.y)
    return(crop_coordinates)
  }

  # filter for data points within 200 y coordinates
  a.df <- df[abs(centroid[2] - df$y) < 200,]

  cluster_pair_distances <- list()
  cluster_pairs <- list(); count = 1
  for (c1 in unique(a.df$cluster)){
    for (c2 in unique(a.df$cluster)){
      if(min(df[df$cluster == c1,]$x) < centroid[1] &
         max(df[df$cluster == c2,]$x) > centroid[1]){
        key <- paste0(c(c1,c2), collapse = "")
        dist = abs(max(a.df[a.df$cluster == c1,]$x) - min(a.df[a.df$cluster == c2,]$x))
        cluster_pair_distances[[key]] = list(distance = dist)
        cluster_pairs[[count]] <- key
        count <- count + 1
      }
    }
  }

  if (length(cluster_pair_distances) == 1){
    c1 <- get_cluster_coordinates(c1, df)
    return(c(c1[1],c1[2],c1[3],c1[4]))
  }

  frame_cluster = strsplit(names(which.max(sapply(cluster_pair_distances,'[[', "distance"))),"")[[1]]

  clusters_to_test <- which(lapply(cluster_pairs, grep, pattern=paste(frame_cluster,collapse="|"),
                                   value = F, invert = T) == 1)

  if (length(clusters_to_test) == 0 ){
    c1 <- get_cluster_coordinates(frame_cluster[1], df)
    c2 <- get_cluster_coordinates(frame_cluster[2], df)
    return(c(min(c1[1],c2[1]),max(c1[2],c2[2]),min(c1[3],c2[3]),max(c1[4],c2[4])))
  }

  max_distance = c(-Inf,0)
  for (c in clusters_to_test){
    pair = strsplit(cluster_pairs[[c]],"")[[1]]
    if(all(pair == pair[1]) & cluster_pair_distances[[c]]$distance < 600) next
    if(cluster_pair_distances[[c]]$distance > max_distance[1]){
      max_distance = c(cluster_pair_distances[[c]]$distance, c)
    }
  }

  if(max_distance[[2]] == 0){
    c1 <- get_cluster_coordinates(frame_cluster[1], df)
    c2 <- get_cluster_coordinates(frame_cluster[2], df)
    return(c(min(c1[1],c2[1]),max(c1[2],c2[2]),min(c1[3],c2[3]),max(c1[4],c2[4])))
  }

  optimal_clusters = strsplit(cluster_pairs[[max_distance[2]]], "")[[1]]

  # get the coordinates of the hand cluster
  if (diff(range(df[df$cluster == optimal_clusters[1],]$y)) < 1000 &
      diff(range(df[df$cluster == optimal_clusters[2],]$y)) < 1000){
    c1 <- get_cluster_coordinates(frame_cluster[1], df)
    c2 <- get_cluster_coordinates(frame_cluster[2], df)
    return(c(min(c1[1],c2[1]),max(c1[2],c2[2]),min(c1[3],c2[3]),max(c1[4],c2[4])))
  }

  c1 <- get_cluster_coordinates(optimal_clusters[1], df)
  c2 <- get_cluster_coordinates(optimal_clusters[2], df)

  crop_coordinates <- c(min(c1[1],c2[1]), max(c1[2],c2[2]), min(c1[3],c2[3]), max(c1[4],c2[4]))
  return(crop_coordinates)
}

# helper function
get_cluster_coordinates<-function(c,df){
  c.subset_df<-df[df$cluster == c,]
  min.x <- c.subset_df[which.min(c.subset_df$x),]$x
  max.x <- c.subset_df[which.max(c.subset_df$x),]$x
  min.y <- c.subset_df[which.min(c.subset_df$y),]$y
  max.y <- c.subset_df[which.max(c.subset_df$y),]$y
  crop_coordinates<-c(min.x,max.x,min.y,max.y)
  return(crop_coordinates)
}

# helper function
find_wrist<-function(df){
  for(i in 2:nrow(df)){
    p1 <- df[i-1,]
    p2 <- df[i+1,]
    slope = (p2$y - p1$y) / (p2$x - p1$x)
    if (is.infinite(slope) & slope < 0){
      return(unlist(df[i,]))
    }
  }
}

# computes the root-mean square contrast
RMS_contrast<-function(im){
  d <- as.data.frame(im)
  average_intensity <- mean(d$value)
  square_mean_difference <- (d$value - average_intensity)^2
  contrast <- sum(square_mean_difference)/ nrow(d)
  contrast <- sqrt(contrast)
  return(contrast)
}

# return true if the image is dark contrast
dark_image<-function(im){
  d <- as.data.frame(im)
  return (mean(d$value) < .3)
}

########## Image Preprocessing ##########
# # preprocess the raw image files
# folderpath = "/boneage-training-dataset/"
# piclist<- list.files(folderpath, pattern=".png", full.names=T)

# # read the images
# for(png in piclist){
#   im = load.image(png)
#   print(paste("filename:", basename(png)))
#   output_dir = "/normalized_images/"
#   destfile <- paste0(output_dir,basename(png))
#   if (!file.exists(destfile)){
#     crop_coor<-segmentHand(im)
#     cropped_image<-im[crop_coor[1]:crop_coor[2],crop_coor[3]:crop_coor[4]] %>% as.cimg

#     padded_im<-squareImage(cropped_image) # get pixel value of background, pad the image with that intensity
#     imager::save.image(padded_im,destfile)
#   }
# }

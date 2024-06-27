rm(list= ls())
pacman::p_load(tidyverse, glue, lubridate, dtwclust)
library(dplyr)
library(ggplot2)
library(fda)
library(caret)
library(sparrow)

# ----Set the working directory to the specified path
setwd("/Users/grace/Desktop/Projectone/RealData/DATA2")

# -----Read the raw dataset CSV file--
champ_word_data = read.csv("rowing_world_championships.csv")

# ------Display the variable names--
var_name = names(champ_word_data)

var_name
# If a race has a speed which is NA or less than 2, throw out the whole race
# There are GPS errors where the average speed recorded is less than the true speed

races_full <- champ_word_data %>%
  dplyr::select(race_id, gender, matches("speed"), matches("strokes")) %>%
  filter_at(vars(matches("speed")), all_vars(!is.na(.) & . > 2)) %>%
  filter_at(vars(matches("stroke")), all_vars(!is.na(.) & . > 0)) %>%
  mutate(row_number = row_number())

str(races_full)

# We have partitioned the dataset based on gender for model fitting purposes
# Filter for men
men_dataset <- races_full %>%
  filter(gender == "men")

# Filter for women
women_dataset <- races_full %>%
  filter(gender == "women")

#write.csv(women_dataset, file = "womendataset.csv")

#write.csv(men_dataset, file = "mendataset.csv")

# Create a vector of column names for speed and stroke
speed_columns <- paste0("speed_", seq(50, 2000, by = 50))
stroke_columns <- paste0("strokes_", seq(50, 2000, by = 50))

# Combine the two sets of column names in the desired order
desired_order <- c(stroke_columns, speed_columns)

# sort variables for Men
races_cleaned_M <-  men_dataset%>%
  dplyr::select(all_of(desired_order))

# sort variables for Women
races_cleaned_W <-  women_dataset%>%
  dplyr::select(all_of(desired_order))

# --------Separate into two datasets ---------
# (one for speed and one for stroke)
races_Y_M <- races_cleaned_M  %>%
  dplyr::select(starts_with("speed_"))

races_A_M <- races_cleaned_M  %>%
  dplyr::select(starts_with("strokes_"))

races_Y_W <- races_cleaned_W  %>%
  dplyr::select(starts_with("speed_"))

races_A_W <- races_cleaned_W  %>%
  dplyr::select(starts_with("strokes_"))

#Standardize the datasets
 speed_standardized_M <- scale_rows(races_Y_M, center = TRUE,scale = TRUE)  # 8137 x 40
 stroke_standardized_M <- scale_rows(races_A_M,center = TRUE,scale = TRUE) # 8137 x 40
# 
# 
# speed_standardized_W <- scale(races_Y_W,center = TRUE,scale = TRUE)  # 8137 x 40
# stroke_standardized_W <- scale(races_A_W,center = TRUE,scale = TRUE) # 8137 x 40
# 

# Standardize the datasets
# speed_standardized_M <- as.matrix( races_Y_M)  # 8137 x 40
# stroke_standardized_M <- as.matrix(races_A_M)  # 8137 x 40

# speed_standardized_W <- as.matrix(races_Y_W)   # 8137 x 40
# stroke_standardized_W <- as.matrix (races_A_W)# 8137 x 40

# ------Center Data-------
# speed_standardized_M <- speed_standardized_M  - apply(speed_standardized_M , 1, mean)
# stroke_standardized_M <- stroke_standardized_M  - apply(stroke_standardized_M ,1,mean)
# 
# speed_standardized_W <- speed_standardized_W  - apply(speed_standardized_W , 1, mean)
# stroke_standardized_W <- stroke_standardized_W  - apply(stroke_standardized_W ,1,mean)


# Add the "speed_0" column to the matrices
# races_speed_standardized <- cbind(speed_0 = zeros, races_speed_standardized)  # 8137 x 41
# races_stroke_standardized <- cbind(stroke_0 = zeros, races_stroke_standardized) # 8137 x 41

# Create a sequence of distances from 50 to 2000m with an interval of 50m
meter <- seq(50, 2000, by = 50)
meters <- seq(50, 2000, by = 50)

#meters = (meter-min(meter))/(max(meter)-min(meter))


# Create line plots for each row in the standardized speed Men dataset
for (i in 1:20) {
  plot(meters, stroke_standardized_M[i,], type = "l", xlab = "Distance (m)" ,
       ylab = "Stroke and Speed", col='red', ylim=c(0,50))
  lines(meters, speed_standardized_M[i,], type = "l", xlab = "Distance (m)", col='blue')
  title(paste("Stroke Speed vs Distance (Row", i, ")"))
}

# Create line plots for each row in the standardized speed Women dataset
for (i in 1000:1020) {
  plot(meters, stroke_standardized_W[i,], type = "l", xlab = "Distance" , col='red')
  lines(meters, speed_standardized_W[i,], type = "l", xlab = "Distance", col='blue')
}


# -------------------pre smoooth speed curve Men------------------------

# Define parameters for the basis
nbasis <- 20
norder <- 4
#rng <- c(0, 1)
rng <- c(50, 2000)
# Create the basis
bsplineBasis <- create.bspline.basis(rng, nbasis, norder)

lambda_values <- seq(1e-4, 1, length.out = 10)

# Initialize a vector to store errors for each lambda
errors <- numeric(length(lambda_values))

for (j in 1:length(lambda_values)) {
  smoothedCurves <- matrix(0, nrow = nrow(speed_standardized_M), ncol = ncol(speed_standardized_M))
  
  for(i in 1:nrow(speed_standardized_M)) {
    curveData <- speed_standardized_M[i,]
    fdParObj <- fdPar(bsplineBasis, Lfdobj=2, lambda=lambda_values[j])
    smoothObj <- smooth.basis(meters, curveData, fdParObj)
    smoothedCurves[i,] <- eval.fd(meters, smoothObj$fd)
  }
  
  # Calculate the error for this lambda value
  errors[j] <- mean((smoothedCurves - speed_standardized_M)^2)
}

# Find the lambda with the smallest error
best_lambda <- lambda_values[which.min(errors)]
cat("The best lambda value is:", best_lambda, "with an MSE of:", min(errors))


dense_points <- seq(50, 2000, length.out = 500)

speed_smoothed_M <- matrix(0, nrow = nrow(speed_standardized_M), ncol = ncol(speed_standardized_M))
dense_speed_smoothed_M <- matrix(0, nrow = nrow(speed_standardized_M), ncol = 500)
for(i in 1:nrow(speed_standardized_M)) {
  curveData <- speed_standardized_M[i,]
  fdParObj <-  fdPar(bsplineBasis, Lfdobj=2, lambda=1e-4) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(meters, curveData, fdParObj)
  speed_smoothed_M[i,] <- eval.fd(meters, smoothObj$fd)
  dense_speed_smoothed_M[i,] <- eval.fd(dense_points, smoothObj$fd)
}

for (i in 1:10) {
  plot(meters, speed_standardized_M[i,], col="red", 
       ylim=c(min(speed_standardized_M[1:80,]), max(speed_standardized_M[1:80,])),
       ylab="Value", 
       xlab="Time", 
       main="Sample Curve with Smoothing")
  lines(meters, speed_smoothed_M[i,], col="blue")
  legend("topright", legend=c("Original", "Smoothed"), 
         fill=c("red", "blue"))
}

for (i in 1:10) {
  plot(meters, speed_smoothed_M[i,], col="red", type = 'l',
       ylab="Speed",
       xlab="Distance(m)", 
       main="Sample Curve with Smoothing")
  lines(dense_points,dense_speed_smoothed_M[i,], col = 'blue')
  legend("topright", legend=c("Smoothed","Dense"), 
         fill=c("red", "blue"))
}


#png("MenSpeed.png", width = 1200, height = 800, res = 150)
set.seed(123)
random_rows <- sample(1:nrow(speed_smoothed_M), 20)

# Set up the plotting space using the first random curve
plot(dense_points, dense_speed_smoothed_M[random_rows[1],], type="n", 
     ylab="Men Speed", 
     xlab="Distance", 
     main="(a)")
for (i in random_rows) {
  lines(dense_points, dense_speed_smoothed_M[i,])
}
#dev.off()


# -------------------pre smoooth speed curve Women------------------------
# Define parameters for the basis
nbasis <- 20
norder <- 4
#rng <- c(0, 1)

# Create the basis
bsplineBasis <- create.bspline.basis(rng, nbasis, norder)

lambda_values <- seq(1e-4, 1, length.out = 10)

# Initialize a vector to store errors for each lambda
errors <- numeric(length(lambda_values))

for (j in 1:length(lambda_values)) {
  smoothedCurves <- matrix(0, nrow = nrow(speed_standardized_W), ncol = ncol(speed_standardized_W))
  
  for(i in 1:nrow(speed_standardized_W)) {
    curveData <- speed_standardized_W[i,]
    fdParObj <- fdPar(bsplineBasis, Lfdobj=2, lambda=lambda_values[j])
    smoothObj <- smooth.basis(meters, curveData, fdParObj)
    smoothedCurves[i,] <- eval.fd(meters, smoothObj$fd)
  }
  
  # Calculate the error for this lambda value
  errors[j] <- mean((smoothedCurves - speed_standardized_W)^2)
}

# Find the lambda with the smallest error
best_lambda <- lambda_values[which.min(errors)]
cat("The best lambda value is:", best_lambda, "with an MSE of:", min(errors))


speed_smoothed_W <- matrix(0, nrow = nrow(speed_standardized_W), ncol = ncol(speed_standardized_W))
dense_speed_smoothed_W <- matrix(0, nrow = nrow(speed_standardized_W), ncol = 500)
for(i in 1:nrow(speed_standardized_W)) {
  curveData <- speed_standardized_W[i,]
  fdParObj <-  fdPar(bsplineBasis, Lfdobj=2, lambda=1e-4) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(meters, curveData, fdParObj)
  speed_smoothed_W[i,] <- eval.fd(meters, smoothObj$fd)
  dense_speed_smoothed_W[i,] <- eval.fd(dense_points, smoothObj$fd)
}

for (i in 1:10) {
  plot(meters, speed_standardized_W[i,], col="red", 
       ylim=c(min(speed_standardized_W[1:80,]), max(speed_standardized_W[1:80,])),
       ylab="Value", 
       xlab="Time", 
       main="Sample Curve with Smoothing")
  lines(meters, speed_smoothed_W[i,], col="blue")
  lines(dense_points, dense_speed_smoothed_W[i,], col = 'green')
  legend("topright", legend=c("Original", "Smoothed"), 
         fill=c("red", "blue"))
}

for (i in 1:10) {
  plot(meters, speed_smoothed_W[i,], col="red", type = 'l',
       ylab="Speed",
       xlab="Distance(m)")
  lines(dense_points, dense_speed_smoothed_W[i,], col = 'blue')
  legend("topright", legend=c("Smoothed","Dense"), 
         fill=c("red", "blue"))
}

#dim(dense_speed_smoothed_W) # 2870 x 500


# -------------------pre smoooth stroke Men------------------------

# Define parameters for the basis
nbasis <- 20
norder <- 4
#rng <- c(0, 1)

# Create the basis
bsplineBasis <- create.bspline.basis(rng, nbasis, norder)

lambda_values <- seq(1e-5, 1e-3, length.out = 10)

# Initialize a vector to store errors for each lambda
errors <- numeric(length(lambda_values))

for (j in 1:length(lambda_values)) {
  smoothedCurves <- matrix(0, nrow = nrow(stroke_standardized_M), ncol = ncol(stroke_standardized_M))
  
  for(i in 1:nrow(stroke_standardized_M)) {
    curveData <- stroke_standardized_M[i,]
    fdParObj <- fdPar(bsplineBasis, Lfdobj=2, lambda=lambda_values[j])
    smoothObj <- smooth.basis(meters, curveData, fdParObj)
    smoothedCurves[i,] <- eval.fd(meters, smoothObj$fd)
  }
  
  # Calculate the error for this lambda value
  errors[j] <- mean((smoothedCurves - stroke_standardized_M)^2)
}

# Find the lambda with the smallest error
best_lambda <- lambda_values[which.min(errors)]
cat("The best lambda value is:", best_lambda, "with an MSE of:", min(errors))



stroke_smoothed_M <- matrix(0, nrow = nrow(stroke_standardized_M), ncol = ncol(stroke_standardized_M))
dense_stroke_smoothed_M <- matrix(0, nrow = nrow(stroke_standardized_M), ncol = 500)
for(i in 1:nrow(stroke_standardized_M)) {
  curveData <- stroke_standardized_M[i,]
  fdParObj <-  fdPar(bsplineBasis, Lfdobj=2, lambda=1e-4) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(meters, curveData, fdParObj)
  stroke_smoothed_M[i,] <- eval.fd(meters, smoothObj$fd)
  dense_stroke_smoothed_M[i,] <- eval.fd(dense_points, smoothObj$fd)
}

# for (i in 1:10) {
#   plot(meters, stroke_standardized_M[i,], col="red", 
#        ylim=c(min(stroke_standardized_M[1:80,]), max(stroke_standardized_M[1:80,])),
#        ylab="Stroke", 
#        xlab="Time", 
#        main="Sample Curve with Smoothing")
#   lines(meters, stroke_smoothed_M[i,], col="blue")
#   lines(dense_points, dense_stroke_smoothed_M[i,], col = 'green')
#   legend("topright", legend=c("Original", "Smoothed"), 
#          fill=c("red", "blue"))
# }

for (i in 1:10) {
  plot(meters,   stroke_smoothed_M[i,], col="red", type = 'l',
       ylab="Stroke",
       xlab="Distance(m)", 
       main="Sample Curve with Smoothing")
  lines(dense_points, dense_stroke_smoothed_M[i,], col = 'blue')
  legend("topright", legend=c("Smoothed","Dense"), 
         fill=c("red", "blue"))
}

#dim(dense_stroke_smoothed_M) # 5225 x 500

#png("MenStroke.png", width = 1200, height = 800, res = 150)
set.seed(123)

# Sample 10 random rows from speed_smoothed_M
random_rows <- sample(1:nrow(stroke_smoothed_M), 20)

# Set up the plotting space using the first random curve
plot(meters, stroke_smoothed_M[random_rows[1],], type="n", 
     ylab="Men Stroke", 
     xlab="Distance", 
     main="(b)")

# Loop to add the curves to the plot
for (i in random_rows) {
  lines(meters, stroke_smoothed_M[i,])  # Using rainbow colors for distinction
}
#dev.off

# -------------------pre smoooth stroke women-----------------------

stroke_smoothed_W <- matrix(0, nrow = nrow(stroke_standardized_W), ncol = ncol(stroke_standardized_W))
dense_stroke_smoothed_W <- matrix(0, nrow = nrow(stroke_standardized_W), ncol = 500)
for(i in 1:nrow(stroke_standardized_W)) {
  curveData <- stroke_standardized_W[i,]
  fdParObj <-  fdPar(bsplineBasis, Lfdobj=2, lambda=1e-4) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(meters, curveData, fdParObj)
  stroke_smoothed_W[i,] <- eval.fd(meters, smoothObj$fd)
  dense_stroke_smoothed_W[i,] <- eval.fd(dense_points, smoothObj$fd)
}

# for (i in 1:10) {
#   plot(meters, stroke_standardized_W[i,], col="red", 
#        ylab="stroke", 
#        xlab="Time", 
#        main="Sample Curve with Smoothing")
#   lines(meters, stroke_smoothed_W[i,], col="blue")
#   lines(dense_points, dense_stroke_smoothed_W[i,], col = 'green')
#   legend("topright", legend=c("Original", "Smoothed"), 
#          fill=c("red", "blue"))
# }

for (i in 1:10) {
  plot(meters,   stroke_smoothed_W[i,], col="red", type = 'l',
       ylab="Stroke",
       xlab="Distance(m)", 
       main="Sample Curve with Smoothing")
  lines(dense_points, dense_stroke_smoothed_W[i,], col = 'blue')
  legend("topright", legend=c("Smoothed","Dense"), 
         fill=c("red", "blue"))
}

#dim(dense_stroke_smoothed_W) # 2870 x 500

set.seed(123)

# Sample 10 random rows from stroke_smoothed_W
random_rows <- sample(1:nrow(stroke_smoothed_W), 10)

# Set up the plotting space using the first random curve
plot(meters, stroke_smoothed_W[random_rows[1],], type="n", 
     ylab="Women Stroke",
     xlab="Distance(m)", 
     main="(b)")

# Loop to add the curves to the plot
for (i in random_rows) {
  lines(meters, stroke_smoothed_W[i,], col=sample(rainbow(10), 1))  # Using rainbow colors for distinction
}


#save data
# write.csv(dense_speed_smoothed_M, file = "centerMenSpeed.csv")
# write.csv(dense_speed_smoothed_W, file = "centerWomenSpeed.csv")
# write.csv(dense_stroke_smoothed_M, file = "centerMenStroke.csv")
# write.csv(dense_stroke_smoothed_W, file = "centerWomenStroke.csv")


# write.csv(dense_speed_smoothed_M, file = "NormMenSpeed.csv")
# write.csv(dense_speed_smoothed_W, file = "NormWomenSpeed.csv")
# write.csv(dense_stroke_smoothed_M, file = "NormMenStroke.csv")
# write.csv(dense_stroke_smoothed_W, file = "NormWomenStroke.csv")


# ----------different types of races 
races_size <- champ_word_data %>%
  dplyr::select(gender, matches("speed"), matches("strokes"), size) %>%
  filter_at(vars(matches("speed")), all_vars(!is.na(.) & . > 2)) %>%
  filter_at(vars(matches("stroke")), all_vars(!is.na(.) & . > 0)) %>%
  mutate(row_number = row_number())

str(races_size)
names(races_size)
# Filter for men
men_class1 <- races_size %>%
  filter(gender == "men" & size %in% 1) %>%
  mutate(row_number = row_number())

# 2 people in boat for men
men_class2 <- races_size %>%
  filter(gender == "men" & size %in% 2) %>%
  mutate(row_number = row_number())

# 4 people in boat for men
men_class4 <- races_size %>%
  filter(gender == "men" & size %in% 4) %>%
  mutate(row_number = row_number())

# 8 people in boat for men
men_class8 <- races_size %>%
  filter(gender == "men" & size %in% 8) %>%
  mutate(row_number = row_number())


# sort variables for Men
Men_Class1 <-  men_class1%>%
  dplyr::select(all_of(desired_order))

Men_Class1_speed <- Men_Class1%>%
  dplyr::select(starts_with("speed_"))

Men_Class1_stroke<- Men_Class1%>%
  dplyr::select(starts_with("strokes_"))

Men_Class1_speed <- scale(Men_Class1_speed)
Men_Class1_stroke <- scale(Men_Class1_stroke)

# Create matrix plots in the standardized speed dataset
matplot(Men_Class1_speed, type ='l')
matplot(Men_Class1_stroke[1:20,], type ='l')

for (i in 1:30){
  # Assuming you have two data frames: races_speed_standardized and races_stroke_standardized
  data1 <- data.frame(Distance = meters, Speed = Men_Class1_speed[i,])
  data2 <- data.frame(Distance = meters, Stroke = Men_Class1_stroke[i,])
  
  # Create a scatterplot
  p <- ggplot() +
    geom_point(data = data1, aes(x = Distance, y = Speed), color = "blue", size = 2) +
    geom_point(data = data2, aes(x = Distance, y = Stroke), color = "red", size = 2) +
    labs(x = "Distance (m)", y = "Standardized Value") +
    ggtitle("Speed vs. Stroke")
  # Print the plot
  print(p)
}

# --------speed data for class 1 Men------------

nbasis <- 20
norder <- 4
rng <- c(50, 2000)

# Create the basis
bsplineBasis <- create.bspline.basis(rng, nbasis, norder)

dense_points <- seq(50, 2000, length.out = 500)

Men_Class1_speed_smooth <- matrix(0, nrow = nrow(Men_Class1_speed ), ncol = length(meters))

dense_Men_Class1_speed <- matrix(0, nrow = nrow(Men_Class1_speed ), ncol = 500)

for(i in 1:nrow(Men_Class1_speed )) {
  curveData <- Men_Class1_speed [i,]
  fdParObj <-  fdPar(bsplineBasis, Lfdobj=2, lambda = 1e-4) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(meters, curveData, fdParObj)
  Men_Class1_speed_smooth[i,] <- eval.fd(meters, smoothObj$fd)
  dense_Men_Class1_speed[i,] <- eval.fd(dense_points, smoothObj$fd)
}

for (i in 10:20) {
  plot(meters, Men_Class1_speed[i,], col="red",
       ylab="Speed", 
       xlab="distance(m)", 
       main="Smoothing Curve Class1")
  
  lines(meters, Men_Class1_speed_smooth[i,], col="blue")
  
  lines(dense_points, dense_Men_Class1_speed[i,], col = 'green')
  
  legend("topright", legend=c("Original", "Smoothed"), 
         fill=c("red", "blue"))
}

dim(dense_Men_Class1_speed) # 1574 x 500

matplot(dense_Men_Class1_speed[1:500,],type = 'l', 
        ylab="Speed",
        xlab="Distance(m)", 
        main="Men Smoothing Speed Curve" )

matplot(Men_Class1_speed_smooth[1:500,],type = 'l', 
        ylab="Speed",
        xlab="Distance(m)", 
        main="Men Smoothing Speed Curve" )


write.csv(dense_Men_Class1_speed, file = "MenSpeedClass1.csv")

# --------stroke data for class 1 Men------------

Men_Class1_stroke = as.matrix(Men_Class1_stroke)
Men_Class1_stroke_smooth <- matrix(0, nrow = nrow(Men_Class1_stroke), ncol = length(meters))

dense_Men_Class1_stroke <- matrix(0, nrow = nrow(Men_Class1_stroke ), ncol = 500)

for(i in 1:nrow(Men_Class1_stroke )) {
  curveData <- Men_Class1_stroke[i,]
  
  fdParObj <-  fdPar(bsplineBasis, Lfdobj=2, lambda = 1e-4) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(meters, curveData, fdParObj)
  Men_Class1_stroke_smooth[i,] <- eval.fd(meters, smoothObj$fd)
  dense_Men_Class1_stroke[i,] <- eval.fd(dense_points, smoothObj$fd)
}

for (i in 10:20) {
  plot(meters, Men_Class1_stroke[i,], col="red",
       ylab="stroke", 
       xlab="distance(m)", 
       main="Smoothing Curve Class1")
  
  lines(meters, Men_Class1_stroke_smooth[i,], col="blue")
  
  lines(dense_points, dense_Men_Class1_stroke[i,], col = 'green')
  
  legend("topright", legend=c("Original", "Smoothed"), 
         fill=c("red", "blue"))
}

dim(dense_Men_Class1_stroke) # 1574 x 500

matplot(dense_Men_Class1_stroke[1:500,],type = 'l', 
        ylab="Stroke",
        xlab="Distance(m)", 
        main="Men Smoothing Stroke Curve" )

matplot(Men_Class1_stroke_smooth[1:500,],type = 'l', 
        ylab="Stroke",
        xlab="Distance(m)", 
        main="Men Smoothing Stroke Curve" )


write.csv(dense_Men_Class1_stroke, file = "MenStrokeClass1.csv")

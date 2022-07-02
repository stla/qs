library(qs)

keyRotors <- rversor(5)
keyTimes <- as.double(1:5)
keyTimes <- numeric(0)
t <- c <- b <- 0
times <- c(1, 2, 3, 5)
nintertimes <- 3L
closed <- FALSE

qs:::KochanekBartels_cpp(
  as.matrix(keyRotors), keyTimes, t, c, b, times, nintertimes, closed
)

qsplines::KochanekBartels(
  keyRotors, keyTimes = NULL, tcb=c(t,c,b), times = times, 
  n_intertimes = nintertimes, endcondition = "natural")

#  bizarre que l'implémentation R réagit à nintertimes quad times est donné




library(qsplines)
# Using a Kochanek-Bartels quaternions spline to construct 
#   a spherical curve interpolating some key points on the 
#     sphere of radius 5

# helper function: spherical to Cartesian coordinates
sph2cart <- function(rho, theta, phi){
  return(c(
    rho * cos(theta) * sin(phi),
    rho * sin(theta) * sin(phi),
    rho * cos(phi)
  ))
}

# construction of the key points on the sphere
keyPoints <- matrix(nrow = 0L, ncol = 3L)
theta_ <- seq(0, 2*pi, length.out = 9L)[-1L]
phi <- 1.3
for(theta in theta_){
  keyPoints <- rbind(keyPoints, sph2cart(5, theta, phi))
  phi = pi - phi
}
n_keyPoints <- nrow(keyPoints)

# construction of the key rotors; the first key rotor 
#   is the identity quaternion and rotor i sends the 
#     key point i-1 to the key point i
keyRotors <- quaternion(length.out = n_keyPoints)
rotor <- keyRotors[1L] <- H1
for(i in seq_len(n_keyPoints - 1L)){
  keyRotors[i+1L] <- rotor <-
    quaternionFromTo(
      keyPoints[i, ]/5, keyPoints[i+1L, ]/5
    ) * rotor
}

# Kochanek-Bartels quaternions spline
rotors <- qs:::KochanekBartels_cpp(
  as.matrix(keyRotors), keyTimes = numeric(0), nintertimes = 25L, 
  closed = TRUE, t = -1, c = 5, b = 0, timesR = (2:14)/2
)
rotors <- as.quaternion(rotors)

# construction of the interpolating points on the sphere
points <- matrix(nrow = 0L, ncol = 3L)
keyPoint1 <- rbind(keyPoints[1L, ])
for(i in seq_along(rotors)){
  points <- rbind(points, rotate(keyPoint1, rotors[i]))
}

# visualize the result with the 'rgl' package
library(rgl)
spheres3d(0, 0, 0, radius = 5, color = "lightgreen")
spheres3d(points, radius = 0.2, color = "midnightblue")
spheres3d(keyPoints, radius = 0.25, color = "red")
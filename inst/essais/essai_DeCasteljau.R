library(qs)

segments <- list(
  rversor(3),
  rversor(4),
  rversor(5)
)
keyTimes <- as.double(1:4)
times <- c(1, 1.6, 1.7, 2.1, 2.3)

DeCasteljau(segments, keyTimes, times)
qsplines::DeCasteljau(segments, keyTimes, times)

ff <- function(q){
  x <- as.numeric(q)
  c(tail(x, -1), head(x, 1))
}
  

q1 <- rversor(1)
q2 <- rversor(1)
qsplines:::.slerp(q1,q2,0.3)
qs:::slerp_(ff(q1), ff(q2), 0.3)

pow <- function(q, t) qs:::slerp_(ff(H1), ff(q), t)

t = 0.3
v <- tail(as.numeric(q2), -1); w <- head(as.numeric(q2), 1)
c(v*sin(t*acos(w))/sqrt(1-w*w), cos(t*acos(w)))





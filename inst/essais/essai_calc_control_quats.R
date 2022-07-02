library(qs)

q1 <- as.quaternion(c(1/sqrt(2),1/sqrt(2),0,0), single = T)
q2 <- as.quaternion(c(0,0,0,1), single = T)
q3 <- as.quaternion(c(1,1,1,1)/2, single = T)
times <- c(0.2, 0.4, 0.6)
tcb <- c(1,1,1)/10
quats <- c(q1,q2,q3)
qsplines:::.calculate_control_quaternions(quats, times, tcb)
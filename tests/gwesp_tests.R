
library(ergm)
library(ergm.sp)
data("faux.dixon.high")

test.approx = function(a, b, message, tol=1e-6) {
  if (abs(a-b) > tol) stop(message)
}

# testing on a larger network (comparing with functions in the ergm package)
net <- faux.dixon.high
if (!all(summary(net ~ dgwesp(fixed=F, type='OTP'))[1:5] == c(516, 222, 65, 21, 3)))
  stop("ESP OTP count on large network incorrect")
test.approx(summary(net ~ dgwesp(fixed=T, type='OTP', alpha=0.1)), 857.4225, tol=1e-3, "GWESP summary")
test.approx(ergm(net ~ edges + dgwesp(fixed=T, alpha=0.1, type="OTP"), estimate = "MPLE")$coef[2], 1.807995, "MPLE estimate on larger nw")

if (summary(net ~ desp(type="OTP", d=2)) != 222) stop("desp OTP count error")
if (summary(net ~ desp(type="ITP", d=2)) != 185) stop("desp ITP count error")
if (summary(net ~ desp(type="OSP", d=2)) != 215) stop("desp OSP count error")
if (summary(net ~ desp(type="ISP", d=2)) != 228) stop("desp ISP count error")


# GWESP OTP test

net <- network.initialize(10, directed = T)
net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- net[1,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1
#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwesp(fixed=F, type="OTP")) )
  if (!all(espcounts[1:4]==c(1,1,0,0))) stop("ESP OTP mis-count")
}
espcounts

test.approx(summary(net~dgwesp(fixed=T, alpha=0.1, type="OTP")), 2.095163, "GWESP_OTP wrong stat")

test.approx(ergm(net ~ edges + dgwesp(fixed=T, alpha=0.1, type="OTP"), estimate = "MPLE")$coef[2], 0.9157069, "GWESP_OTP ergm MPLE wrong estimate")


# GWESP ITP test

net <- network.initialize(10, directed = T)
net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- net[4,1] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwesp(fixed=F, type="ITP")) )
  if (!all(espcounts[1:4]==c(7,1,0,0))) stop("ESP ITP mis-count")
}
espcounts

test.approx(summary(net~dgwesp(fixed=T, alpha=0.1, type="ITP")), 8.095163, "GWESP_ITP wrong stat")

test.approx(ergm(net ~ edges + dgwesp(fixed=T, alpha=0.1, type="ITP"), estimate = "MPLE")$coef[2], 1.97197, "GWESP_ITP ergm MPLE wrong estimate")



# GWESP OSP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- net[3,2] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[10,7] <-  net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwesp(fixed=F, type="OSP")) )
  if (!all(espcounts[1:4]==c(1,1,0,0))) stop("ESP_OSP mis-count")
}
espcounts

test.approx(summary(net~dgwesp(fixed=T, alpha=0.1, type="OSP")), 2.095163 , "GWESP_OSP wrong stat")

test.approx(ergm(net ~ edges + dgwesp(fixed=T, alpha=0.1, type="OSP"), estimate = "MPLE")$coef[2], 1.137139, "GWESP_OSP ergm MPLE wrong estimate")


# GWESP ISP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- net[1,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[7,10] <-  net[7,9] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwesp(fixed=F, type="ISP")) )
  if (!all(espcounts[1:4]==c(1,1,0,0))) stop("ESP ISP mis-count")
}
espcounts

test.approx(summary(net~dgwesp(fixed=T, alpha=0.1, type="ISP")), 2.095163, "GWESP_ISP wrong stat")

test.approx(ergm(net ~ edges + dgwesp(fixed=T, alpha=0.1, type="ISP"), estimate = "MPLE")$coef[2], 1.137139, "GWESP_ISP ergm MPLE wrong estimate")



library(ergm)
library(ergm.sp)
data("faux.dixon.high")

test.approx = function(a, b, message, tol=1e-5) {
  if (abs(a-b) > tol) stop(message)
}

# testing on a larger network (comparing with functions in the ergm package)
net <- faux.dixon.high
if (!all(summary(net ~ dgwnsp(fixed=F, type='OTP'))[1:5] == c(4355, 855,163,32,14)))
  stop("NSP OTP count on large network incorrect")
test.approx(summary(net ~ dgwnsp(fixed=T, type='OTP', alpha=0.1)), 5522.186 , "GWNSP summary", tol=1e-3)
test.approx(ergm(net ~ edges + dgwnsp(fixed=T, alpha=0.1, type="OTP"), estimate = "MPLE")$coef[2], -0.07421213 , "MPLE estimate on larger nw")

test.approx(summary(net ~ dgwdsp(fixed=T, type='OTP', alpha=0.1)), summary(net ~ dgwesp(fixed=T, type='OTP', alpha=0.1)) + summary(net ~ dgwnsp(fixed=T, type='OTP', alpha=0.1)),
            "GWDSP should equal GWESP + GWNSP")

if (summary(net ~ dnsp(type="OTP", d=2)) != 855) stop("dnsp OTP count error")
if (summary(net ~ dnsp(type="ITP", d=2)) != 892) stop("dnsp ITP count error")
if (summary(net ~ dnsp(type="OSP", d=2)) != 733) stop("dnsp OSP count error")
if (summary(net ~ dnsp(type="ISP", d=2)) != 904) stop("dnsp ISP count error")

# GWNSP OTP test

net <- network.initialize(10, directed = T)
net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1
#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwnsp(fixed=F, type="OTP")) )
  if (!all(espcounts[1:4]==c(3,1,0,0))) stop("NSP OTP mis-count")
}
espcounts

test.approx(summary(net~dgwnsp(fixed=T, alpha=0.1, type="OTP")), 4.095163, "GWNSP_OTP wrong stat")

test.approx(ergm(net ~ edges + dgwnsp(fixed=T, alpha=0.1, type="OTP"), estimate = "MPLE")$coef[2], -0.8563381 , "GWNSP_OTP ergm MPLE wrong estimate")



# GWNSP ITP test

net <- network.initialize(10, directed = T)

net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwnsp(fixed=F, type="ITP")) )
  if (!all(espcounts[1:4]==c(1,1,0,0))) stop("NSP ITP mis-count")
}
espcounts

test.approx(summary(net~dgwnsp(fixed=T, alpha=0.1, type="ITP")), 2.095163, "GWNSP_ITP wrong stat")

test.approx(ergm(net ~ edges + dgwnsp(fixed=T, alpha=0.1, type="ITP"), estimate = "MPLE")$coef[2], -1.760084 , "GWNSP_ITP ergm MPLE wrong estimate")



# GWNSP OSP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[10,7] <-  net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwnsp(fixed=F, type="OSP")) )
  if (!all(espcounts[1:4]==c(3,2,0,0))) stop("NSP_OSP mis-count")
}
espcounts

test.approx(summary(net~dgwnsp(fixed=T, alpha=0.1, type="OSP")), 5.190325 , "GWNSP_OSP wrong stat")

test.approx(ergm(net ~ edges + dgwnsp(fixed=T, alpha=0.1, type="OSP"), estimate = "MPLE")$coef[2], -0.2865438, "GWNSP_OSP ergm MPLE wrong estimate")


# GWNSP ISP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[7,10] <-  net[7,9] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwnsp(fixed=F, type="ISP")) )
  if (!all(espcounts[1:4]==c(3,2,0,0))) stop("NSP ISP mis-count")
}
espcounts

test.approx(summary(net~dgwnsp(fixed=T, alpha=0.1, type="ISP")), 5.190325, "GWNSP_ISP wrong stat")

test.approx(ergm(net ~ edges + dgwnsp(fixed=T, alpha=0.1, type="ISP"), estimate = "MPLE")$coef[2], -0.2865438, "GWNSP_ISP ergm MPLE wrong estimate")



require(ergm.terms.contrib)
data("faux.dixon.high")

test.approx = function(a, b, message, tol=1e-5) {
  if (abs(a-b) > tol) stop(message)
}

# testing on a larger network (comparing with functions in the ergm package)
net <- faux.dixon.high
if (!all(summary(net ~ dgwdsp(fixed=F, type='OTP'))[1:5] == c(4871, 1077, 228, 53, 17)))
  stop("DSP OTP count on large network incorrect")
test.approx(summary(net ~ dgwdsp(fixed=T, type='OTP', alpha=0.1)), 6379.609, tol=1e-3, "GWDSP summary")

test.approx(ergm(net ~ edges + dgwdsp(fixed=T, alpha=0.1, type="OTP"), estimate = "MPLE")$coef[2], 0.006467, "MPLE estimate on larger nw")

# OTP and ITP are the same for dyadwise SP
if (summary(net ~ ddsp(type="OTP", d=2)) != 1077) stop("ddsp OTP count error")
if (summary(net ~ ddsp(type="ITP", d=2)) != 1077) stop("ddsp ITP count error")
if (summary(net ~ ddsp(type="OSP", d=2)) != 948) stop("ddsp OSP count error")
if (summary(net ~ ddsp(type="ISP", d=2)) != 1132) stop("ddsp ISP count error")

# GWDSP OTP test

net <- network.initialize(10, directed = T)
net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1
#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwdsp(fixed=F, type="OTP")) )
  if (!all(espcounts[1:4]==c(4,1,0,0))) stop("DSP OTP mis-count")
}
espcounts

test.approx(summary(net~dgwdsp(fixed=T, alpha=0.1, type="OTP")), 5.095163, "GWDSP_OTP wrong stat")

test.approx(ergm(net ~ edges + dgwdsp(fixed=T, alpha=0.1, type="OTP"), estimate = "MPLE")$coef[2], -1.2014723, "GWDSP_OTP ergm MPLE wrong estimate")

# larger network


# GWDSP ITP test

net <- network.initialize(10, directed = T)

net[1,2] <- net[2,4] <- net[1,3] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[8,5] <- 1
net[10,7] <- net[10,9] <- net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwdsp(fixed=F, type="ITP")) )
  if (!all(espcounts[1:4]==c(4,1,0,0))) stop("DSP ITP mis-count")
}
espcounts

test.approx(summary(net~dgwdsp(fixed=T, alpha=0.1, type="ITP")), 5.095163, "GWDSP_ITP wrong stat")

test.approx(ergm(net ~ edges + dgwdsp(fixed=T, alpha=0.1, type="ITP"), estimate = "MPLE")$coef[2], -1.2014723, "GWDSP_ITP ergm MPLE wrong estimate")



# GWDSP OSP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[10,7] <-  net[9,7] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwdsp(fixed=F, type="OSP")) )
  if (!all(espcounts[1:4]==c(4,2,0,0))) stop("DSP_OSP mis-count")
}
espcounts

test.approx(summary(net~dgwdsp(fixed=T, alpha=0.1, type="OSP")), 6.190325, "GWDSP_OSP wrong stat")

test.approx(ergm(net ~ edges + dgwdsp(fixed=T, alpha=0.1, type="OSP"), estimate = "MPLE")$coef[2], -0.2104664, "GWDSP_OSP ergm MPLE wrong estimate")


# GWDSP ISP test

net <- network.initialize(10, directed = T)

net[2,1] <- net[2,4] <- net[3,1] <- net[3,4] <- 1
net[5,6] <- net[6,8] <- net[5,8] <- 1
net[7,10] <-  net[7,9] <- 1

#plot(net, edge.lwd=2, arrowhead.cex=2, label=1:network.size(net))

for (i in 1:100) {
  (espcounts <- summary(net ~ dgwdsp(fixed=F, type="ISP")) )
  if (!all(espcounts[1:4]==c(4,2,0,0))) stop("DSP ISP mis-count")
}
espcounts

test.approx(summary(net~dgwdsp(fixed=T, alpha=0.1, type="ISP")), 6.190325, "GWDSP_ISP wrong stat")

test.approx(ergm(net ~ edges + dgwdsp(fixed=T, alpha=0.1, type="ISP"), estimate = "MPLE")$coef[2], -0.2104664, "GWDSP_ISP ergm MPLE wrong estimate")


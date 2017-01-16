require(ergm.userterms)
require(sna)
require(testthat)

# par(mar=c(0,0,0,0))
net <- network.initialize(5,directed=F)
net%v%"pendants" <- c("zgreen","red","black", "zgreen", "red")
net[1,2] <- 1
net[1,3] <- 1
net[3,4] <- 1
net[4,2] <- 1
net[3,5] <- 1
net[2,5] <- 1

# plot(net, displaylabels = T,  vertex.col = as.factor(net%v%"pendants"),edge.lwd=.1, vertex.border=0, label.pos = 4, vertex.cex = 3)


test_that('gwnspattr term', {
  
# how many black-red dyads have ANY shared partners?
expect_equal(as.vector(summary(net ~ gwnspattr("pendants",1,2,0))), 1)

# how many shared partners are there for all black-red dyads?
expect_equal(as.vector(summary(net ~ gwnspattr("pendants",1,2,1000))), 3)

# non-extreme alphas weight each additional shared partner less than one before
expect_equal(as.vector(summary(net ~ gwnspattr("pendants",1,2,1))), 2.031697, tolerance=1e-5)

})



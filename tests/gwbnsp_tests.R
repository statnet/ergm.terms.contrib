
require(testthat)
require(ergm.terms.contrib)

bel <- data.frame(actor = c("a","a","a","b","b","b","b","c","c","c","c"), 
                  event = c("e1","e2","e3","e3","e1","e4","e5","e3","e1","e4","e6"))

bnet <- network(bel,bipartite=3,directed=F)


# plot for visualy checking the debugging
# bnet%v%"mode" <- c(rep("actors",3),rep("events",6))
# col <- as.numeric(as.factor(bnet%v%"mode")) + 1
# plot(bnet, vertex.cex = 5, vertex.col = col, displaylabels = T, label.pos = 5)

test_that('gwbnsp terms', {

  # how many actor dyads have ANY shared partners?
  expect_equal(as.vector( summary(bnet ~ gwb1nsp(0,T)) ), 3)
  
  # what is the TOTAL number of shared partners of all actor dyads?
  expect_equal(as.vector( summary(bnet ~ gwb1nsp(100,T)) ), 7) 
  
  # how many event dyads have ANY shared partners?
  expect_equal(as.vector( summary(bnet ~ gwb2nsp(0,T)) ), 11)  
  
  # what is the TOTAL number of shared partners of all event dyads?
  expect_equal(as.vector( summary(bnet ~ gwb2nsp(100,T)) ) , 15) 
  
  # non-extreme alphas weight each additional shared partner less than the one before
  expect_equal(as.vector( summary(bnet ~ gwb1nsp(.5,T)) ) , 4.335226, tolerance=1e-5) 
  
  # same as above but for 2nd mode
  expect_equal(as.vector( summary(bnet ~ gwb2nsp(.5,T)) ), 12.33523 , tolerance=1e-5)
})




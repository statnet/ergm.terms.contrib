require(testthat)
require(ergm)
require(ergm.terms.contrib)
data(florentine)

# check that it returns the expected results for example network
#degree(flobusiness,gmode = 'graph')
#[1] 0 0 4 3 3 2 2 4 5 1 4 0 0 1 0 1
expect_equal(as.vector(summary(flobusiness~maxdegree(0))) ,5)
expect_equal(as.vector(summary(flobusiness~maxdegree(4))) ,15)
expect_equal(as.vector(summary(flobusiness~maxdegree(30))) ,16)

# check that terms returned are named correctly
expect_equal(names(summary(flobusiness~maxdegree(20))), 'maxdegree20')

# check that the attribute matching works
# flobusiness%v%'priorates'
# [1] 53 65  0 12 22  0 21  0 53  0 42  0 38 35 74  0
expect_equal(as.vector(summary(flobusiness~maxdegree(2))) ,10)
# since there is only one edge across matching priorates, they will all below degree 2 if we count priorates
expect_equal(as.vector(summary(flobusiness~maxdegree(2,by='priorates'))) ,16)
# check naming
expect_equal(names(summary(flobusiness~maxdegree(2,by='priorates'))), 'maxdegree.priorates2')

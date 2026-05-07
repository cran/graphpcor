## ----setup--------------------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE, 
  warning = FALSE,
  fig.width = 9, 
  out.width = '0.49\\textwidth', 
  fig.align = 'center')
library(graphpcor)
chinclude = TRUE

## ----graphs, echo = FALSE, fig.width = 15, fig.height = 7,  out.width="99%", fig.cap = "Graphs representing different models. A model with one parent, (A), a model with two parents (B), and two different models with three parents (C) and (D)."----
tree1 <- treepcor(p1 ~ c1 + c2 - c3)
tree2 <- treepcor(p1 ~ p2 + c1 + c2, 
                  p2 ~ -c3 + c4)
tree3a <- treepcor(p1 ~ p2 + p3 + c1 + c2, 
                   p2 ~ -c3 + c4, 
                   p3 ~ c5 + c6)
tree3b <- treepcor(p1 ~ p2 + c1 + c2, 
                   p2 ~ p3 - c3 + c4, 
                   p3 ~ c5 + c6)
par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))
plot(tree1)
legend("topleft", "(A)", bty = "n")
plot(tree2)
legend("topleft", "(B)", bty = "n")
plot(tree3a)
legend("topleft", "(C)", bty = "n")
plot(tree3b)
legend("topleft", "(D)", bty = "n")

## ----graph1, include = chinclude----------------------------------------------
tree1 <- treepcor(p1 ~ c1 + c2 - c3)
tree1
summary(tree1)

## ----dim1---------------------------------------------------------------------
dim(tree1)

## ----visgraph1, include = chinclude-------------------------------------------
plot(tree1)

## ----qs1, include = chinclude-------------------------------------------------
tpcQ(tree1)

## ----q, include = chinclude---------------------------------------------------
q1 <- tpcQ(tree1, theta = 0)
q1

## ----v1, include = chinclude--------------------------------------------------
vcov(tree1) ## assume theta = 0 (\gamma_1 = 1)
vcov(tree1, theta = 0.5) # \gamma_1^2 = exp(2 * 0.5) = exp(1)
cov1a <- vcov(tree1, theta = 0) 
cov1a

## ----c1, include = chinclude--------------------------------------------------
c1 <- cov2cor(cov1a)
round(c1, 3)

## ----graph2-------------------------------------------------------------------
tree2 <- treepcor(
  p1 ~ p2 + c1 + c2, 
  p2 ~ -c3 + c4)
dim(tree2)
tree2
summary(tree2)

## ----visgraph2----------------------------------------------------------------
plot(tree2)

## ----drop1--------------------------------------------------------------------
drop1(tree2)

## ----q2-----------------------------------------------------------------------
q2 <- tpcQ(tree2, theta = c(0, 0))
q2

## ----c2-----------------------------------------------------------------------
cov2 <- vcov(tree2, theta = c(0, 0))
cov2
c2 <- cov2cor(cov2)
round(c2, 3)

## ----graph2b------------------------------------------------------------------
tree2b <- treepcor(
  p1 ~ -p2 + c1 + c2, 
  p2 ~ -c3 + c4)
tree2b
summary(tree2b)

## ----prec2--------------------------------------------------------------------
q2b <- tpcQ(tree2b, theta = c(0, 0))
q2b

## ----cov2b--------------------------------------------------------------------
all.equal(solve(q2)[1:4, 1:4], 
          solve(q2b)[1:4, 1:4])

## ----vfig1--------------------------------------------------------------------
vcov(tree1, raw = TRUE)
vcov(tree2, raw = TRUE)
tree3a <- treepcor(p1 ~ p2 + p3 + c1 + c2, 
                   p2 ~ -c3 + c4, 
                   p3 ~ c5 + c6)
vcov(tree3a, raw = TRUE)
tree3b <- treepcor(p1 ~ p2 + c1 + c2, 
                   p2 ~ p3 - c3 + c4, 
                   p3 ~ c5 + c6)
vcov(tree3b, raw = TRUE)

## ----tree2c-------------------------------------------------------------------
tree2
(np <- dim(tree2))

## ----sigmas-------------------------------------------------------------------
sigmas <- c(4, 2, 1, 0.5)

## ----pparams------------------------------------------------------------------
thetal <- c(0, 1)

## ----Vg-----------------------------------------------------------------------
theta1 <- c(log(sigmas), thetal)
Vg <- vcov(tree2, theta = theta1)
Vg

## ----data2--------------------------------------------------------------------
nrep <- 100
nd <- nrep * np[1]

xx <- matrix(rnorm(nd), nrep) %*% chol(Vg)
cov(xx)

theta.y <- log(10)
datar <- data.frame(
    r = rep(1:nrep, np[1]),
    i = rep(1:np[1], each = nrep),
    y = rnorm(nd, 1 + xx, exp(-2*theta.y))
)

## ----cg2----------------------------------------------------------------------
cmodel <- cgeneric(
  tree2, lambda = 5, 
  sigma.prior.reference = rep(1, np[1]), 
  sigma.prior.probability = rep(0.05, np[1]))

## ----haveINLA, echo = FALSE---------------------------------------------------
haveINLA <- FALSE

## ----ccheck, eval = haveINLA--------------------------------------------------
# cgeneric_cgeneric_graph(cmodel)
# cgeneric_initial(cmodel)
# prior(cmodel, theta = rep(0, sum(np)))
# prior(cmodel, theta = rep(1, sum(np)))
# 
# np
# tpcQ(cmodel, theta = rep(0, sum(np)))
# (Qc <- tpcQ(cmodel, theta = theta1))
# all.equal(Vg, as.matrix(solve(Qc)))

## ----mfit, eval = haveINLA----------------------------------------------------
# m1 <- y ~ f(i, model = cmodel, replicate = r)
# pfix <- list(prec = list(initial = 10, fixed = TRUE))
# fit <- inla(
#     formula = m1,
#     data = datar,
#     control.family = list(hyper = pfix)
# )
# fit$cpu.used

## ----Vcompare, eval = haveINLA------------------------------------------------
# round(Vg, 2)
# round(cov(xx), 2)
# round(Vfit <- vcov(tree2, theta = fit$mode$theta), 2)

## ----Cfit, eval = haveINLA----------------------------------------------------
# round(Cfit <- vcov(tree2, theta = fit$mode$theta[np[1]+1:np[2]]), 2)
# round(cov2cor(Vfit), 2)

## ----t4, echo = FALSE,  fig.width = 9, fig.height = 4, out.width="69%", fig.align = 'center', fig.cap = "Two different models for the same number of children and parents."----
t4a <- treepcor(
  p1 ~ p2 + p3 + c1 + c2, 
  p2 ~ p4 + c3 + c4, 
  p3 ~ p5 + c5 + c6, 
  p4 ~ c7 + c8, p5 ~ c9 + c10)
t4b <- treepcor(p1 ~ p2 + p3 + c1 + c2, 
  p2 ~ c3 + c4, p3 ~ p4 + p5 + c5 + c6, 
  p4 ~ c7 + c8, p5 ~ c9 + c10)
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
plot(t4a)
legend("topleft", "(E)", bty = "n")
plot(t4b)
legend("topleft", "(F)", bty = "n")


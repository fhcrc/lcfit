#!/usr/bin/env Rscript

library(ape)

N_TREES <- 10
N_LEAVES <- 10

dir.create('trees')

set.seed(1)
for(i in 1:N_TREES) {
  write.tree(rtree(N_LEAVES), sprintf('trees/tree%04d.tre', i))
}

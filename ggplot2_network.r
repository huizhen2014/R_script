##ggplot2 for network

library(GGally)
library(ggplot2)
library(sna)
library(intergraph)
library(RColorBrewer)

##example 1: random graph / undirected Bernoulli random graph

net = rgraph(10,mode="graph",tprob=0.5)

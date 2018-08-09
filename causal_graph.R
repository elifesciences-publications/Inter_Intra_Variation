library(pcalg)
library(igraph)
df              <- data.sc %>% dplyr::select( vm:fi, dvlocmm) %>% na.omit
N               <- nrow(df)
Rho_neurons     <- cor(df)
tol             <- 1.96/sqrt(N)           # rough approximation for

suffStat <- list( C = cor(df), n=nrow(df) )
pc.fit   <- pc( suffStat, indepTest = gaussCItest, alpha = 0.05, p=ncol(df) )
g        <- igraph.from.graphNEL( pc.fit@graph, name = TRUE, weight = TRUE, unlist.attrs = TRUE )
df.g     <- as_data_frame( g, "edges" ) %>% apply( 2, as.numeric ) %>% as.data.frame
tmp      <- vector( 'list', nrow(df.g) )
for(i in 1:dim(df.g)[1]) tmp[[i]] <- ida(df.g$from[i], df.g$to[i], cor(df), pc.fit@graph)

edgewidth    <- 5*( unlist(tmp))                        #width of edges
edgecolour   <- ifelse(edgewidth > 0, "black", "red")   #colours



## 
## Plot the graph 
##

## pdf("causal_graph.pdf")

nms <- names(df)
V(g)$name    <- names(df)
V(g)$size <- 25
E(g)$weight <- edgewidth
plot(g, order=V(g), layout=layout.circle, edge.width = abs(edgewidth), edge.color=edgecolour, edge.arrow.size=1.5*(abs(E(g)$weight)))

## dev.off()

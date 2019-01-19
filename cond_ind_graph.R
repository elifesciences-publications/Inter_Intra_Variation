library(ggraph)
library(igraph)
data.sc_neurons <- data.sc %>% dplyr::select(vm:fi, dvlocmm) %>%
    na.omit
id <- (data.sc %>% dplyr::select(vm:fi, dvlocmm, id) %>%
    na.omit %>% dplyr::select(id))$id
N               <- nrow(data.sc)
Rho_neurons     <- cor(data.sc_neurons)
tol <- 0
Q_neurons       <- corpcor::cor2pcor(Rho_neurons)
Q_neurons       <- Q_neurons %>% as.data.frame

## bootstrap, resampling rows of data matrix
M <- 1000
index <- 1:nrow(data.sc_neurons)
## index <- 1:length(unique(id))
## uniqueid <- unique(id)
Q_neurons_star <- array(dim=c(ncol(data.sc_neurons), ncol(data.sc_neurons), M))
tmp <- NULL
for(i in 1:M){
    index_star <- base::sample(x=index, length(index), replace=TRUE)
    tmp <- rbind(tmp, data.sc[index_star,])
    ## for(j in 1:length(uniqueid[index_star])){
    ##     tmp <- rbind(tmp, data.sc[which(id == uniqueid[index_star][j]),])
    ## }
    Rho_neurons_star <- cor(data.sc_neurons[index_star,])
    Q_neurons_star[,,i] <- as.matrix(corpcor::cor2pcor(Rho_neurons_star))
    tmp <- NULL
}

## bootstrap, resampling animals
if(FALSE){
    M <- 1000
    index <- 1:nrow(data.sc_neurons)
    index <- 1:length(unique(id))
    uniqueid <- unique(id)
    Q_neurons_star <- array(dim=c(ncol(data.sc_neurons), ncol(data.sc_neurons), M))
    tmp <- NULL
    for(i in 1:M){
        index_star <- base::sample(x=index, length(index), replace=TRUE)
        for(j in 1:length(uniqueid[index_star])){
            tmp <- rbind(tmp, data.sc[which(id == uniqueid[index_star][j]),])
        }
        Rho_neurons_star <- cor(data.sc_neurons[index_star,])
        Q_neurons_star[,,i] <- as.matrix(corpcor::cor2pcor(Rho_neurons_star))
        tmp <- NULL
    }
}


Q_neurons_low <- apply(Q_neurons_star, c(1,2), quantile, 0.025)
Q_neurons_upp <- apply(Q_neurons_star, c(1,2), quantile, 0.975)    
CI <- Q_neurons_low * Q_neurons_upp
CI[CI<0]  <- 0
CI[CI!=0] <- 1
CI <- as.data.frame(CI)
Q_neurons <- CI*Q_neurons
g              <- igraph::graph.adjacency(abs(Q_neurons)>0,
                                          mode="undirected", diag=FALSE)
## g              <- igraph::graph.adjacency(abs(CI)>0,
##                                   mode="undirected", diag=FALSE)
igraph::V(g)$class     <- c(data.sc_r$property, "dvlocmm")
igraph::V(g)$degree    <- igraph::degree(g)
edge_width    <- NULL
edge_colour      <- NULL
k <- 0
for(i in 1:(ncol(Q_neurons)-1))
{
  for(j in (i+1):ncol(Q_neurons))
    if(abs(Q_neurons[i,j])>0)
    {
      k<-k+1
      edge_width[k]    <- abs(Q_neurons[i,j])/2
      edge_colour[k]   <- ifelse(Q_neurons[i,j] > 0 , "black", "red")
    }
}
igraph::E(g)$width     <- edge_width
igraph::E(g)$colour    <- edge_colour

## graph          <- tidygraph::as_tbl_graph(g)
## ##
## p <- ggraph(graph, 'igraph', algorithm = 'circle') +
##   ## ggraph(graph, layout= "kk") +
##   geom_node_circle(aes( r=0.1), fill="orange", colour="black")+
##   geom_edge_fan(aes(width = edge_width, color= colour), 
##                 show.legend = FALSE,
##                 start_cap = circle(1.4, 'cm'),
##                 end_cap = circle(1.4,  'cm'),
##                 lineend = "butt") +
##   ## scale_edge_color_hue(colour=c("black","gray"))+
##   geom_node_text(aes(label=class),
##                  size=6, nudge_x=.0, nudge_y=0.0, colour="blue") +
##   coord_fixed()+ 
##   theme_graph(foreground = 'steelblue', fg_text_colour = 'white')
## ## 

Q <- as.matrix(Q_neurons)
colnames(Q) <- colnames(data.sc_neurons)
rownames(Q) <- colnames(data.sc_neurons)




## pdf("CIG_neurons.pdf", height=13, width=13)
par(mfrow=c(1,2))
set.seed(111020)
corrplot(Q)
g$layout <- layout_with_fr
plot(g, vertex.label=colnames(data.sc_neurons), vertex.size=20,
     edge.width=20*edge_width, edge.color=edge_colour)
## dev.off()


## pdf("CIG_neurons_circle.pdf", height=13, width=13)
par(mfrow=c(1,2))
set.seed(111020)
corrplot(Q)
g$layout <- layout_in_circle
plot(g, vertex.label=colnames(data.sc_neurons), vertex.size=20,
     edge.width=20*edge_width, edge.color=edge_colour)
## dev.off()















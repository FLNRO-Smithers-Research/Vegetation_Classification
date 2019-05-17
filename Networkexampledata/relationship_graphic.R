#####An example of a relationship diagram
rm(list=ls())
wd <- tk_choose.dir(); setwd(wd)

library(igraph)
#library(gcookbook) # For the data set

#nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
#links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

nodes <- read.csv("testingunit_nodes.csv", header=T, as.is=T)
nodes <- nodes[!duplicated(nodes$id),]
links <- read.csv("testingunit_links.csv", header=T, as.is=T)

# Examine the data:
head(nodes)
head(links)
nrow(nodes); length(unique(nodes$id))
nrow(links); nrow(unique(links[,c("from", "to")]))

# Collapse multiple links of the same type between the same two nodes
# by summing their weights, using aggregate() by "from", "to", & "type":
# (we don't use "simplify()" here so as not to collapse different link types)
links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
class(net)
net 

# We can look at the nodes, edges, and their attributes:
E(net)
V(net)
E(net)$type
V(net)$media

plot(net, edge.arrow.size=.4,vertex.label=NA)

# Removing loops from the graph:
net <- simplify(net, remove.multiple = F, remove.loops = T) 

# If you need them, you can extract an edge list or a matrix from igraph networks.
as_edgelist(net, names=T)
as_adjacency_matrix(net, attr="weight")

# Or data frames describing nodes and edges:
as_data_frame(net, what="edges")
as_data_frame(net, what="vertices")

# Replace the vertex label with the node names stored in "media"
plot(net, edge.arrow.size=.2, edge.curved=0,
     vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(net)$id, vertex.label.color="black",
     vertex.label.cex=.3) 

# We can also add a legend explaining the meaning of the colors we used:
plot(net) 
legend(x=-1.1, y=-1.1, c("Newspaper","Television", "Online News"), pch=5,
       col="#777777", pt.bg=colrs, pt.cex=0.5, bty="n", ncol=1)

#net.bg <- sample_pa(80, 1.2) 
V(net)$size <- 6
V(net)$frame.color <- "white"
V(net)$color <- "orange"
#V(net)$label <- "" 
E(net)$arrow.mode <- 0
#plot(net)
plot(net, edge.arrow.size=.2, edge.curved=0,
     vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(net)$id, vertex.label.color="black",
     vertex.label.cex=.3) 

l <- layout_with_fr(net)# Layout with fruchterman.reingold. Get the layout coordinates:
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)# Normalize them so that they are in the -1, 1 interval:
plot(net, layout=l*1.2,
     rescale=F,
     edge.arrow.size=.2,
     edge.curved=0,
     vertex.size        = 4,          # Smaller nodes
     vertex.label       = V(net)$id,  # Set the labels
     vertex.label.cex   = 0.4,        # Slightly smaller font
     vertex.label.dist  = 0.4,        # Offset the labels
     vertex.color="orange",
     vertex.frame.color="#555555",
     vertex.label.color = "black")

#################################
# Check out all available layouts in igraph:
?igraph::layout_

layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]

par(mfrow=c(3,3), mar=c(1,1,1,1))

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(net)) 
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }

dev.off()


# Copy madmen and drop every other row
m <- madmen[1:nrow(madmen) %% 2 == 1, ]
g <- graph.data.frame(m, directed=FALSE)

plot(g, layout=layout.fruchterman.reingold,
     vertex.size        = 4,          # Smaller nodes
     vertex.label       = V(g)$name,  # Set the labels
     vertex.label.cex   = 0.8,        # Slightly smaller font
     vertex.label.dist  = 0.4,        # Offset the labels
     vertex.label.color = "black")

g$layout <- layout.fruchterman.reingold

plot(g)

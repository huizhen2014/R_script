##特殊字符
##NA 缺失值或为定义数据
##NULL 空对象(null/empty lists)
##Inf/-Inf 正或负无限值
##NaN 不能定义结果
##is.na()/is.finite()/is.nan()
##(v1>2) | (v2>3)  '|' boolean OR, 返回向量
##(v1>2) & (v2>3)  '&' boolean AND,返回向量
##(v1>2)||(v2>3)  '||' boolean OR，返回单个值
##(v1>2)&&(v2>3)  '&&' boolean AND,ditto

###Netwroks in igraph
##10各节点,默认有向,1,2,3 vectrices are directed
g <- graph(edges=c(1,2,2,3,3,1),n=10)
plot(g)
g2 <- graph(c("John","Jim","Jim","Jill","Jill","John"))
plot(g2)
g3 <- graph(c("John","Jim","Jim","Jack","Jim","Jack","John","John"),
            isolates=c("Jesse","Janis","Jennifer","Justin"))
plot(g3,edge.arrow.size=0.5,vertex.color="gold",vertex.size=15,
     vertex.frame.color="gray",vertex.label.color="black",
     vertex.label.cex=0.8,vertex.label.dist=2,edge.curved=0.2)
##使用-描述无向联系，+—和-+用于有向连接，分别朝左和右，++为对称连接，
##: 用于一套向量
plot(graph_from_literal(a---b,b---c))
plot(graph_from_literal(a-+b,b+-c))
plot(graph_from_literal(a+-+b,b+-+c))
#plot(graph_from_literal(a:b:c-c:d:e))
#plot(graph_from_literal(a:b:c-a:b:c))

##edge/vertex/network属性
E(g3) #返回对象edges,
V(g3) #返回对象的vertices,V(g3)$name，返回名称
g3[] #返回对象的矩阵(dgCMatrix),从左到右一次mate
V(g3)$gender <- c("male","male","male","male","female","female","male")
E(g3)$type <- "email"
E(g3)$weight <- 10
edge_attr(g3) #返回edge属性
vertex_attr(g3) #返回vertex属性
graph_attr(g3) #返回图像属性
g3 <- set_graph_attr(g3,"name","Email Netwrok") #设置图形属性
g3 <- set_graph_attr(g3,"something","A thing") #设置图形属性
#或使用set_edge_attr()/set_vertex_attr()设置对应属性
g3 <- delete_graph_attr(g3,"something") #删除对应图形属性
plot(g3,edge.arrow.size=0.5,vertex.label.color="black",vertex.label.dist=1.5,
     vertex.color=c("pink","skyblue")[1+(V(g3)$gender=="male")])
##去除自身loop和multiple edges,使用edge.attr.comb来指明edge
##attributes合并的可能选项:sum,mean,prod(product),min,max,
##first/last(选择第一/最后一个edge属性)，ignore表示忽略或舍弃对应属性
g3s <- simplify(g3,remove.multiple = T,remove.loops = F,
                edge.attr.comb = c(weight="sum",type="ignore"))
plot(g3s,vertex.label.dist=1.5)
##igraph对象描述
##1. D/U，表示有向或无向的图形
##2. N表示命名了的图形(nodes拥有name属性)
##3. W表示weighted 图形(edges拥有weight属性)
##4. B表示bipartite(two-mode)图形(ndoes拥有type属性)
##随后的连个数字(7 3)指得是图形中node和edge数目
##同时还有对应nodes/edges属性
##(g/c) graph-level charcter attibute
##(v/c) vertex-level character attribute
##(e/n) edge-level numeric attribute

##特殊图形和图形模式
##Empty graph
eg <- make_empty_graph(40)
plot(eg,vertex.size=10,vertex.label=NA)
##Full graph
fg <- make_full_graph(40)
plot(fg,vertex.size=10,vertex.label=NA)
##Simple star graph
st <- make_star(40)
plot(st,vertex.size=10,vertex.label=NA)
##Tree graph
tr <- make_tree(40,children=3,mode="undirected")
plot(tr,vertex.size=10,vertex.label=NA)
##Ring graph
rn <- make_ring(40)
plot(rn, vertex.size=10,vertex.label=NA)
##Erdos-Renyi random graph mdoel
##n为nodes数目，m为edges数目
er <- sample_gnm(n=100,m=40)
plot(er,vertex.size=6,vertex.label=NA)
##Watts-Strongatz small-wrold model
##构建格子图,dim为维度，size为贯穿维度的nodes数目,指定概率p来rewire edges,
##相邻edges连接nei，允许loops和multiple edges
sw <- sample_smallworld(dim=2,size=10,nei=1,p=0.1)
plot(sw,vertex.size=6,vertex.label=NA,layout=layout_in_circle)
##...

##Reading network data from files
nodes <- read.csv("Dataset1-Media-Example-NODES.csv",header = T,as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv",header = T,as.is = T)
##存在两相同nodes间多重links，使用aggregate函数，collapse来将相同type合并
##起来，by summing their weight
links <- aggregate(links[,3],links[,-3],sum)
links <- links[order(links$from,links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL
##Reading network matrix
nodes2 <- read.csv("Dataset2-Media-User-Example-NODES.csv",
                   header=T,as.is=T)
links2 <- read.csv("Dataset2-Media-User-Example-EDGES.csv",
                   header=T,as.is=T)
##将links2看作two-mode的相邻矩阵
links2 <- as.matrix(links2)

##truning networks into igraph objects
##graph.data.frame/d/vertices
##d 指明network的edges，前两列为每个edge的：source IDs和
##target mode IDs;随后的列为edge属性(weight,type,label...)
##vertices以一个node IDs为起点，随后列为node属性
library(igraph)
net <- graph_from_data_frame(d=links,vertices = nodes,directed=T)

plot(net)
##去除loops
net <- simplify(net,remove.multiple = F,remove.loops = T)
plot(net)
edge_attr(net)
vertex_attr(net)
graph_attr(net)
##也可使用simplify(net,edge.attr.comb=list(weight="sum","Ignore"))
##该过程也会将多重edge types合并(例如，hyperlinks/mentions)
##获取edge list
as_edgelist(net,names=T)
as_adjacency_matrix(net,attr="weight")
##或以数据框形式展示nodes和edges
as_data_frame(net,what="edges")
as_data_frame(net,what="vertices")

###
net2 <- graph_from_incidence_matrix(links2)
net2[]
##在igraph中，bipartite双向网络拥有node属性type
##(FALSE/0为one mode vertices TRUE/1 为其他mode)
table(V(net2)$type)

##将one-mode network matrix转换成igraph对象, 使用graph_from_adjacency_matrix()
##从two-mode network获得二分图(co-memberships are easy to 
##calculate by nultiplying the network matrix by its transposed
##matrix, or using igrpah's vipartite.projection()function)
net2.bp <- bipartite.projection(net2)
##calculate the projections manually
as_incidence_matrix(net2) %*% t(as_incidence_matrix(net2))
t(as_incidence_matrix(net2)) %*% as_incidence_matrix(net2)

plot(net2.bp$proj1,vertex.label.color="black",vertex.label.dist=1,
     vertex.size=7,vertex.label=nodes$media[!is.na(nodes2$media.type)])
plot(net2.bp$proj2,vertext.label.color="black",vertext.label.dist=1,
     vertex.size=7,vertex.label=nodes2$media[is.na(nodes2$media.type)])

##plotting networkds with igraph 
##曲线edge
plot(net,edge.arrow.size=0.4,edge.curved=0.1)
##edge颜色，node颜色，使用node names 'media'取代vertex label
plot(net,edge.arrow.size=0.2,edge.curved=0,
     vertex.color="orange",vertex.frame.color="#555555",
     vertex.label=V(net)$media,vertex.label.color="black",
     vertex.label.cex=0.7)
##设置igraph对象属性。例如，根据media类型设置颜色，根据audience size
##设置size，根据weight设置edges的width
colors <- c("gray50","tomato","gold")
V(net)$color <- colors[V(net)$media.type]
V(net)$size <- V(net)$audience.size*0.7
V(net)$label.color <- "black"
V(net)$label <- NA
E(net)$width <- E(net)$weight/6
E(net)$arrow.size <- 0.2
E(net)$edge.color <- "gray80"
E(net)$width <- 1+E(net)$weight/12

##在plot参数中设置属性，覆盖之前设置
plot(net,edge.color="orange",vertex.color="gray50")

##使用legen设置图例
plot(net)
legend(x=-1.5,y=-1.1,c("Newspaper","Television","Online News"),
       pch=21,col="#777777",pt.bg=colors,pt.cex=2,cex=0.8,bty="n",
       ncol=1)
##或仅绘制nodes的标签图
plot(net,vertex.shape="none",vertex.label=V(net)$media,
     vertex.label.font=2,vertex.label.color="gray40",
     vertex.label.cex=0.7,edge.color="gray85")
##根据node来源对edge着色;graph,es(sequence of edges to query),names=T/F
edge.start <- ends(net,es=E(net),names=F)[,1]
edge.col <- V(net)$color[edge.start]
plot(net,edge.color=edge.col,edge.curved=0.1)

##network layouts
##network layouts are simply algorithms that return coordinates 
##for each node in a network
net.bg <- sample_pa(80)
V(net.bg)$size <- 8
V(net.bg)$frame.color <- "white"
V(net.bg)$color <- "orange"
V(net.bg)$label <- ""
E(net.bg)$arrow.mode <- 0
plot(net.bg)
##设置plot fuction layout
plot(net.bg,layout=layout_randomly)
##计算vertex坐标，绘制环妆图
l <- layout_in_circle(net.bg)
##l 为图中N各nodes的x，y简单矩阵坐标
plot(net.bg,layout=l)
##igraph内置layout
l <- layout_randomly(net.bg)
plot(net.bg,layout=l)
l <- layout_in_circle(net.bg)
plot(net.bg,layout=l)
l <- layout_on_sphere(net.bg)
plot(net.bg,layout=l)

##Fruchterman-Reingold 为最有用的force-directed layout algorithms
##Force-directed layouts，nice-looking graph/edges similar in length
##and cross each other as little as possible. They simulate the graph
##as a physical system. Nodes are electrically charged particles that
##repulse each other when they get too close. The edges act as springs 
##that attract connnected nodes closer together.
l <- layout_with_fr(net.bg)
plot(net.bg,layout=l)

##每次绘制，layout是不同的,同时显示其不同配置
par(mfrow=c(2,2),mar=c(0,0,0,0))
plot(net.bg,layout=layout_with_fr)
plot(net.bg,layout=layout_with_fr)
plot(net.bg,layout=l)
plot(net.bg,layout=l)

##默认，plots的坐标均rescale为[-1,1]范围。
##rescale=FALSE，取消该设置;使用norm_coords重置
l <- layout_with_fr(net.bg)
l <- norm_coords(l,ymin=-1,ymax=1,xmin=-1,xmax=1)
par(mfrow=c(2,2),mar=c(0,0,0,0))
plot(net.bg,rescale=F,layout=l*0.4)
plot(net.bg,rescale=F,layout=l*0.6)
plot(net.bg,rescale=F,layout=l*0.8)
plot(net.bg,rescale=F,layout=l*1.0)

##Kamada Kawai，forced-directed algorithm
##attempts to minimize the energy in a spring system
l <- layout_with_kk(net.bg)
plot(net.bg,layout=l)

##LGL algorithem意味着，large,connected graphs
plot(net.bg,layout=layout_with_lgl)
layouts <- grep("^layout_",ls("package:igraph"),value=T)[-1]
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree",layouts)]

par(mfrow=c(3,3),mar=c(1,1,1,1)) #设置图形配置和边缘距离
for(layout in layouts){
  print(layout)
  l <- do.call(layout,list(net))
  plot(net,edge.arrow.mode=0,layout=l,main=layout)
}
##根据type分别绘制不同edge颜色
E(net)$width <- 1.5
plot(net,edge.color=c("dark red","slategrey")[(E(net)$type=="hyperlink")+1],
     vertex.color="gray40",layout=layout.circle)

##根据link类型分别绘制
net.m <- net - E(net)[E(net)$type == "hyperlink"]
net.h <- net - E(net)[E(net)$type == "mention"]
par(mfrow=c(1,2))
plot(net.h, vertex.color="orange",main="Tie: Hyperlink")
plot(net.m, vertex.color="lightsteelblue2",main="Tie: Mention")

##绘制二分图,根据type(0/FALSE,1/TRUE)改变node形状
V(net2)$color <- c("steel blue","orange")[V(net2)$type + 1]
V(net2)$shape <- c("square","circle")[V(net2)$type + 1]
V(net2)$label <- ""
V(net2)$label[V(net2)$type == F] <- nodes2$media[V(net2)$type == F]
V(net2)$label.cex=0.4
V(net2)$label.font=2

##针对二分图，采用biparitite network
plot(net2,vertex.label=NA,vertex.size=7,layout=layout_as_bipartite)

##使用文本作为nodes
plot(net2,vertex.shape="none",vertex.label=nodes$media,
     vertex.label.color=V(net2)$color,vertex.label.font=2.5,
     vertex.label.cex=0.6,edge.color="gray70",edge.width=2)

##Hubs and authorities
##Hubs were expected to contain catalogs with a large number of 
##outgoning links ; while authorities would get many incoming 
##links from hubs, presumably because of their high-quality
##relevant information
hs <- hub_score(net,weight=NA)$verctor
hs <- hub_score(net_mf,weight=NA)$vector
as <- authority_score(net_mf,weight=NA)$vector
plot(net_mf,vertex.size=hs*50,main="hubs")
plot(net_mf,vertex.size=as*30,main="Authorities")

##Community detection
##Community detecton based on edge betweenness(Newman-Girvan)
ceb <- cluster_edge_betweenness(net_mf)
dendPlot(ceb,mode="hclust")
plot(ceb,net_mf)
##community detecton based on based on propagating labels
clp <- cluster_label_prop(net_mf)
plot(clp,net_mf)

##grep("^layout_",ls("package:igraph"),value=T)
##输出graph至文件
##write.graph(g,file="my_graph.dl",format="pajek")
##write.graph(g,file='my_graph.txt',fromat="edgelist")

##legned
##image.plot(legend.only = TRUE,zlim=range(rescale(-log10(vertices$qvalue))),
library(fields)
image.plot(legend.only = TRUE,zlim=range(rescale(-log10(vertices$qvalue))), ##将颜色范围均一到0，1
           col=rbPal(100)[cut(seq(0,1,by=0.001),breaks=100)], ##设定颜色范围
           horizontal = TRUE,legend.shrink=0.2,legend.width = 1,
           legend.args=list(text="-log10(Qvalue)",cex=0.8,line=0.1), ##legend标签，大小，距离
           axis.args=list(at=c(0,0.5,1),labels=c(0,0.5,1),cex=0.5,line=0.05), ##axis 位置，labels，距离
           smallplot = c(0.8,0.9,0.85,0.9)) ## 未知，擦





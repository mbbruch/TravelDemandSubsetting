rm(list=ls());
base.path <- "/ocean/projects/eng200002p/mbruchon/"
data.path <- paste0(base.path,"trips/")
.libPaths( c(paste0(base.path,"R_Libs"), .libPaths())) 
#install.packages(c("data.table","lubridate","bit64")
library(data.table)
library(lubridate)
library(bit64)
library(reshape2)
library(dtwclust)
library(BBmisc)
library(cluster)
library(FNN)
library(igraph)
library(ggplot2)
library(tidyr)
library(sf)
library(ggplot2)

alreadyBuilt <- TRUE
#Getting relevant columns
if(!alreadyBuilt){
	trips <- fread(paste0(data.path,"Transportation_Network_Providers_-_Trips.csv"))
	trips[,':='(
		`Trip ID`=NULL,
		`Pickup Centroid Location`=NULL,
		`Dropoff Centroid Location`=NULL,
		`Fare`=NULL,
		`Tip`=NULL,
		`Additional Charges`=NULL,
		`Trip Total`=NULL
		)]
	setnames(trips, "Pickup Community Area", "start.area")
	setnames(trips, "Dropoff Community Area", "end.area")
	setnames(trips, "Pickup Census Tract", "start.tract")
	setnames(trips, "Dropoff Census Tract", "end.tract")
	setnames(trips, "Pickup Centroid Latitude", "start.lat")
	setnames(trips, "Pickup Centroid Longitude", "start.long")
	setnames(trips, "Dropoff Centroid Latitude", "end.lat")
	setnames(trips, "Dropoff Centroid Longitude", "end.long")
	setnames(trips,"Shared Trip Authorized","pool.allowed")
	setnames(trips,"Trips Pooled","pool.size")
	setnames(trips,"Trip Miles","miles")
	setnames(trips,"Trip Seconds","seconds")
	column.subset <- trips[,.(start.date,start.datetime,start.area,end.area,start.tract,end.tract,start.long,start.lat,end.long,end.lat)]
	column.subset <- column.subset[year(start.date) == 2019]
	column.subset[,':='(hour=hour(start.datetime))]
	column.subset[,':='(shifted.date=fifelse(hour>4,start.date,start.date-1))]
	unique.dates <- data.table(shifted.date=column.subset[,unique(shifted.date)])
	unique.dates[,':='(week=as.Date(cut(as.Date(shifted.date), "week")))]
	column.subset[unique.dates,week.date:=week,on=.(shifted.date=shifted.date)]
	column.subset[,':='(volume.daily.od=.N),by=.(start.area,end.area,shifted.date)]
	saveRDS(column.subset, file = paste0(data.path,"chicago_pooling_subset.rds"))
} else{
	column.subset <- readRDS(paste0(data.path,"chicago_pooling_subset.rds"))
}

#Determining the most representative day(s)
if(!alreadyBuilt){
	column.subset[,':='(time.rank=frank(start.datetime,ties.method="dense")),by=week.date]
	column.subset[,':='(hour.rank=time.rank-(time.rank-1)%%4)]
	column.subset[,':='(hour.rank=frank(hour.rank,ties.method="dense")),by=week.date]
	start.rollup <- column.subset[,.(start.volume=.N),by=.(week.date,start.area,hour.rank)]
	end.rollup <- column.subset[,.(end.volume=.N),by=.(week.date,end.area,hour.rank)]
	all.combos <- data.table(crossing(week.date=unique(column.subset$week.date),area=unique(c(column.subset$end.area,column.subset$start.area)),hour.rank=unique(column.subset$hour.rank)))
	all.combos[start.rollup,start.volume:=start.volume,on=.(week.date,hour.rank,area=start.area)]
	all.combos[end.rollup,end.volume:=end.volume,on=.(week.date,hour.rank,area=end.area)]
	all.combos[is.na(end.volume),end.volume:=0]
	all.combos[is.na(start.volume),start.volume:=0]
	all.combos[,volume:=end.volume+start.volume]
	all.combos[,weekly.start.volume:=sum(start.volume),by=.(area,week.date)]
	all.combos[,weekly.end.volume:=sum(end.volume),by=.(area,week.date)]
	all.combos[,start.frac:=start.volume/max(start.volume),by=.(area,week.date)]
	all.combos[,end.frac:=end.volume/max(end.volume),by=.(area,week.date)]
	all.combos[,week.rank:=frank(week.date,ties.method="dense")]
	all.combos <- all.combos[week.rank != min(week.rank) & week.rank != max(week.rank)]
	for (i in names(all.combos)) all.combos[is.na(get(i)), (i):=0]
	start.curves.wide <- data.table(dcast(all.combos,area+week.date ~ hour.rank, value.var="start.frac"))
	end.curves.wide <- data.table(dcast(all.combos,area+week.date ~ hour.rank, value.var="end.frac"))
	vbh <- copy(start.curves.wide)
	vbh[,':='(area=NULL,week.date=NULL)]
	clust.pam <- tsclust(vbh, type="partitional", k=4, distance="dtw_basic", centroid="pam")
	start.curves.wide$cluster <- clust.pam@cluster
	write.csv(clust.pam@centroids,paste0(data.path,"start_curve_pam_clusters.csv"))
	vbh <- copy(end.curves.wide)
	vbh[,':='(area=NULL,week.date=NULL)]
	clust.pam <- tsclust(vbh, type="partitional", k=4, distance="dtw_basic", centroid="pam")
	end.curves.wide$cluster <- clust.pam@cluster
	write.csv(clust.pam@centroids,paste0(data.path,"end_curve_pam_clusters.csv"))
	
	all.combos[start.curves.wide,start.cluster:=cluster,on=.(area,week.date)]
	all.combos[end.curves.wide,end.cluster:=cluster,on=.(area,week.date)]
	start.cluster.volumes <- all.combos[,.(weekly.volume=sum(volume,na.rm=TRUE)),by=.(week.date,start.cluster)]
	end.cluster.volumes <- all.combos[,.(weekly.volume=sum(volume,na.rm=TRUE)),by=.(week.date,end.cluster)]
	start.wide <- data.table(dcast(start.cluster.volumes, week.date ~ start.cluster, value.var="weekly.volume"))
	end.wide <- data.table(dcast(end.cluster.volumes, week.date ~ end.cluster, value.var="weekly.volume"))
	for (i in names(start.wide)) end.wide[is.na(get(i)), (i):=0]
	for (i in names(end.wide)) end.wide[is.na(get(i)), (i):=0]
	total.wide <- start.wide[end.wide,':='(end.1=i.1,end.2=i.2,end.3=i.3,end.4=i.4),on=.(week.date)]
	write.csv(total.wide,paste0(data.path,"week_volumes_by_cluster.csv"))
	vbcw.copy <- copy(total.wide)
	vbcw.copy[,week.date:=NULL]
	set.seed(123)
	kmedioid.cluster <- pam(vbcw.copy, k=1, metric = "euclidean", stand = FALSE)
	week.date.to.use <- total.wide[kmedioid.cluster$id.med,]$week.date
	#pdf(file = "/ocean/projects/eng200002p/mbruchon/trips/plot.pdf", width = 1500, height = 800)
	#plot(clust.pam[[1]],type="sc")
	#dev.off()
	trips <- column.subset[week.date==week.date.to.use]
	saveRDS(trips, file = paste0(data.path,"representative_week.rds"))
} else{
	trips <- readRDS( paste0(data.path,"representative_week.rds"))
}
trips <- trips[(!is.na(start.tract) | !is.na(start.area)) & (!is.na(end.tract) | !is.na(end.area))]
trips[,':='(start.node=NULL,end.node=NULL)]
trips[,':='(start.node=as.integer64(NA),end.node=as.integer64(NA))]

#Process edge lists
swr = function(string, nwrap=10) {
     paste(strwrap(string, width=nwrap), collapse="\n")
}
nodes.path <- paste0(base.path,"trips/edgelists")
swr = Vectorize(swr)
files = list.files(path=nodes.path, pattern="node_list.csv",full.names=TRUE, recursive=TRUE) 
node.tract.combos=lapply(files, fread)
tracts <- as.integer64(substr(files, nchar(files)-24, nchar(files)-14))
for (i in 1:length(node.tract.combos)){node.tract.combos[[i]]$tract <- tracts[i]}
node.tract.combos <- rbindlist(node.tract.combos)
node.tract.combos <- unique(node.tract.combos[,.(osmid,tract,x,y)])
edges.path <- paste0(base.path,"trips/edgelists")
swr = Vectorize(swr)
files = list.files(path=edges.path, pattern="edge_list.csv",full.names=TRUE, recursive=TRUE) 
edges.tract.combos=lapply(files, fread)
tracts <- as.integer64(substr(files, nchar(files)-24, nchar(files)-14))
for (i in 1:length(edges.tract.combos)){edges.tract.combos[[i]]$tract <- tracts[i]}
edges.tract.combos <- rbindlist(edges.tract.combos)
edge.tract.combos <- unique(edges.tract.combos[,.(osmid,tract,u,v,highway,length)])
nodes.ua <- fread("/ocean/projects/eng200002p/mbruchon/trips/16264_Chicago_IL--IN_Urbanized_Area/node_list.csv")
edges.ua <- fread("/ocean/projects/eng200002p/mbruchon/trips/16264_Chicago_IL--IN_Urbanized_Area/edge_list.csv")
edges <- unique(rbindlist(list(
	unique(edge.tract.combos[,.(osmid,u,v)]),
	unique(edges.ua[,.(osmid,u,v)]))))
nodes <- unique(rbindlist(list(
	unique(node.tract.combos[,.(osmid)]),
	unique(nodes.ua[,.(osmid)]))))
node.outdegree <- edges[,.(degree=.N),by=u]
node.indegree <- edges[,.(degree=.N),by=v]
node.degrees <- rbindlist(list(
	edges[,.(degree=.N),by=u][,osmid:=u][,u:=NULL],
	edges[,.(degree=.N),by=v][,osmid:=v][,v:=NULL]
))
node.degrees <- node.degrees[,degree:=sum(degree,na.rm=TRUE),by=osmid]
nodes[node.degrees,degree:=degree,on=.(osmid=osmid)][is.na(degree),degree:=0]
node.tract.combos[node.degrees,degree:=degree,on=.(osmid=osmid)][is.na(degree),degree:=0]
nodes.ua[node.degrees,degree:=degree,on=.(osmid=osmid)][is.na(degree),degree:=0]
node.locations.ua <- unique(nodes.ua[,.(osmid,x,y)])
node.locations.tract <- unique(node.tract.combos[,.(osmid,x,y)])
nodes[node.locations.ua,':='(x=x,y=y),on=.(osmid=osmid)]
nodes[node.locations.tract,':='(x=fifelse(is.na(x),i.x,x),y=fifelse(is.na(y),i.y,y)),on=.(osmid=osmid)]
rm(node.locations.ua); rm(node.locations.tract); rm(node.outdegree); rm(node.indegree); rm(node.degrees); 

edges.more.info <- unique(rbindlist(list(
	edges.tract.combos[,.(u=pmin(u,v),v=pmax(u,v),highway,length)],
	edges.ua[,.(u=pmin(u,v),v=pmax(u,v),highway,length)]	
)))
edges.more.info <- edges.more.info[,.(
	length=min(length,na.rm=TRUE),
	residential=fifelse(max(highway!='residential',na.rm=TRUE)==1,"Other Road Type(s)","Residential Only")
	),by=.(u,v)]


if(!alreadyBuilt){
	combos <- unique(rbindlist(list(
		unique(column.subset[,.(tract=start.tract,area=start.area)]),
		unique(column.subset[,.(tract=end.tract,area=end.area)])
		)))
	combos <- combos[!is.na(area)]
	saveRDS(combos, file = paste0(data.path,"tract_area_combos.rds"))
 } else {
	combos <- readRDS( paste0(data.path,"tract_area_combos.rds"))
 }

node.tract.combos[combos,':='(area=area),on=.(tract)]

chicago.tracts <- fread(paste0(data.path,"all_chicago_census_tracts.csv"))
chicago.tracts[,tract:=as.integer64(GEOID10)]
tract.trips <- rbindlist(list(
	trips[,.N,by=start.tract][,tract:=start.tract][,start.tract:=NULL],
	trips[,.N,by=end.tract][,tract:=end.tract][,end.tract:=NULL]
))
tract.trips <- tract.trips[,.(trips=sum(N,na.rm=TRUE)),by=tract]
tract.areas <- rbindlist(list(
	trips[,.N,by=start.area][,area:=start.area][,start.area:=NULL],
	trips[,.N,by=end.area][,area:=end.area][,end.area:=NULL]
))
tract.areas <- tract.areas[,.(trips=sum(N,na.rm=TRUE)),by=area]
#tract.trips <- tract.trips[trips >= 5][,tract.character:=as.character(tract)]
node.tract.combos[tract.trips,in.trip.tracts:=fifelse(is.na(i.tract),FALSE,TRUE),on=.(tract=tract)][is.na(in.trip.tracts),in.trip.tracts:=FALSE]
node.tract.combos[tract.areas,in.trip.areas:=fifelse(is.na(i.area),FALSE,TRUE),on=.(area=area)][is.na(in.trip.areas),in.trip.areas:=FALSE]
node.tract.combos[chicago.tracts,in.chicago:=fifelse(is.na(i.tract),FALSE,TRUE),on=.(tract=tract)][is.na(in.chicago),in.chicago:=FALSE]
node.tract.combos[,':='(in.trips= in.trip.tracts | in.trip.areas)][,':='(to.keep = in.trips)] #& in.chicago)]

library(tripack)
tr<-tri.mesh(unique(node.tract.combos[to.keep==TRUE,.(x,y)]))
nodes.ua[,in.hull:=in.convex.hull(tr, nodes.ua$x, nodes.ua$y)]
node.tract.combos[,in.hull:=in.convex.hull(tr, node.tract.combos$x, node.tract.combos$y)]
nodes[,in.hull:=in.convex.hull(tr, nodes$x, nodes$y)]
png('/ocean/projects/eng200002p/mbruchon/trips/rplot1.png')
pdf(file = "/ocean/projects/eng200002p/mbruchon/trips/plot1.pdf")
ggplot(nodes[in.hull==TRUE], aes(x=x, y=y,color=in.hull)) + geom_point()
dev.off()

combos.trip <- unique(rbindlist(list(
	unique(trips[,.(tract=start.tract,area=start.area)]),
	unique(trips[,.(tract=end.tract,area=end.area)])
	)))
for(i in 1:nrow(combos.trip)){
	this.tract <- combos.trip[i,tract]
	this.area <- combos.trip[i,area]
	if(is.na(this.tract)){
		starting <- which(is.na(trips$start.tract) & trips$start.area==this.area)
		trips[starting, start.node := node.tract.combos[degree !=2 & area==this.area][sample(.N,length(starting),replace=TRUE),osmid]]
		ending <- which(is.na(trips$end.tract) & trips$end.area==this.area)
		trips[ending, end.node := node.tract.combos[degree !=2 & area==this.area][sample(.N,length(ending),replace=TRUE),osmid]]			
	} else{
		starting <- which(trips$start.tract==this.tract)
		start.to.use <- node.tract.combos[degree !=2 & tract==this.tract][sample(.N,length(starting),replace=TRUE),osmid]
		if(length(start.to.use)==0) start.to.use <- node.tract.combos[tract==this.tract][sample(.N,length(starting),replace=TRUE),osmid]
		if(length(start.to.use)==0  & length(starting)>0){
			index <- get.knnx(node.tract.combos[,.(x,y)],
				unique(trips[start.tract==this.tract,.(start.lat,start.long)]),10)$nn.index[1:5]
			start.to.use <- node.tract.combos[index,][sample(.N,length(starting),replace=TRUE),osmid]
		}
		trips[starting, start.node:= start.to.use]
		
		ending <- which(trips$end.tract==this.tract)
		end.to.use <- node.tract.combos[degree !=2 & tract==this.tract][sample(.N,length(ending),replace=TRUE),osmid]	
		if(length(end.to.use)==0) end.to.use <- node.tract.combos[tract==this.tract][sample(.N,length(ending),replace=TRUE),osmid]	
		if(length(end.to.use)==0 & length(ending)>0){
			index <- get.knnx(node.tract.combos[,.(x,y)],
				unique(trips[end.tract==this.tract,.(end.lat,end.long)]),10)$nn.index[1:5]
			end.to.use <- node.tract.combos[index,][sample(.N,length(ending),replace=TRUE),osmid]
		}
		trips[ending, end.node:=end.to.use]
	}
}
trips[nodes,':='(start.lat=y, start.long=x),on=.(start.node=osmid)]
trips[nodes,':='(end.lat=y, end.long=x),on=.(end.node=osmid)]
trips[,':='(start.node=NULL,end.node=NULL)]
nodes.to.keep <- nodes[in.hull==TRUE,osmid] 
edges.more.info[nodes,':='(x.start=x,y.start=y),on=.(u=osmid)]
edges.more.info[nodes,':='(x.end=x,y.end=y),on=.(v=osmid)]
pdf('/ocean/projects/eng200002p/mbruchon/trips/segment_types.pdf')
ggplot() +theme_light() +geom_segment(data=edges.more.info[u %in% nodes.to.keep | v %in% nodes.to.keep],mapping=
                         aes(x=x.start,y=y.start,xend=x.end,yend=y.end,color=residential),size=0.15) #+ scale_color_brewer(palette='Dark2')
dev.off()

edges <- edges.more.info[(u %in% nodes.to.keep | v %in% nodes.to.keep) & u!=v & residential=="Other Road Type(s)"] # highway !='residential']
nodes.to.keep <- unique(c(edges$u,edges$v))
areas <- st_read(paste0(data.path,"Boundaries - Community Areas (current)/geo_export_6c2c6aa8-3ec3-445f-83dd-85a3ae13bfd3.shp"))
nodes <- unique(nodes[osmid %in% nodes.to.keep,.(osmid,x,y)])
nodes <- nodes[,newid:=frank(osmid,ties.method="dense")][order(newid)]
makestring_one <- function(input){
  this.point <- st_point(rbind(c(input[1],input[2])))
  return(areas$area_num_1[st_within(this.point,areas)[[1]]])
}
#temp <- apply(nodes[,.(x,y)],1,makestring_one)
#saveRDS(temp,paste0(data.path,"makestring_one.RDS")
temp <- readRDS(paste0(data.path,"makestring_one.RDS"))
nodes$area <- 0
for(i in 1:nrow(nodes)){
	nodes[i,]$area <- max(c(as.integer(temp[[i]]),0),na.rm=TRUE)
}
train <- nodes[area>0]
test <- nodes[area==0]
index <- get.knnx(train[,.(x,y)],test[,.(x,y)],1)$nn.index
test$new.area <- train[index,area]
nodes[test,area:=new.area,on=.(newid)]
pdf(file = "/ocean/projects/eng200002p/mbruchon/trips/node_areas.pdf")
ggplot(nodes, aes(x=x, y=y,color=as.factor(area))) + geom_point(size=0.1,show.legend=FALSE) + theme_light()
dev.off()


edges[nodes,node1:=newid,on=.(u=osmid)]
edges[nodes,node2:=newid,on=.(v=osmid)]
node.tract.combos[nodes,newid:=newid,on=.(osmid=osmid)]
node_areas <- node.tract.combos[osmid %in% nodes.to.keep,.(area=max(area)),by=newid][is.na(area),area:=0]
nodes[node_areas,area:=area,on=.(newid=newid)]
nodes[is.na(area),area:=0]

g <- graph_from_data_frame(d=data.frame(source=edges$node1,target=edges$node2),
                           vertices=data.frame(name=nodes$newid),directed=FALSE)
clusters <- clusters(g)
connected.nodes <- as.numeric(V(g)$name[which(clusters$membership==which.max(clusters$csize)[1])])
#png('/ocean/projects/eng200002p/mbruchon/trips/rplot.png')
#pdf(file = "/ocean/projects/eng200002p/mbruchon/trips/plot.pdf")
#ggplot(nodes[connected.nodes,], aes(x=x, y=y)) + geom_point()
#dev.off()
nodes <- nodes[newid %in% connected.nodes]
edges <- edges[node1 %in% connected.nodes | node2 %in% connected.nodes]
nodes[,newid:=frank(newid,ties.method="dense")]
edges[nodes,':='(node1=newid,x.start=x,y.start=y),on=.(u=osmid)]
edges[nodes,':='(node2=newid,x.end=x,y.end=y),on=.(v=osmid)]
unique <- unique(edges[,.(node1,node2)])
unique[,':='(min=pmin(node1,node2),max=pmax(node1,node2))]
unique <- unique(unique[,.(min,max)])
edges[,':='(min=pmin(node1,node2),max=pmax(node1,node2))]
avg <- edges[,.(dist=mean(length,na.rm=TRUE)),by=.(min,max)]
avg[nodes,':='(x.start=x,y.start=y),on=.(min=newid)]
avg[nodes,':='(x.end=x,y.end=y),on=.(max=newid)]

if(!alreadyBuilt){
	these.tracts <- st_read(paste0("/ocean/projects/eng200002p/mbruchon/trips/tl_2015_17_tract/","tl_2015_17_tract.shp"))
	makestring <- function(input){
		this.arc <- st_linestring(rbind(c(input[1],input[2]),
									  c(input[3],input[4])))
		return(these.tracts$GEOID[st_intersects(this.arc, these.tracts)[[1]]])
	}

	temp <- apply(avg[,.(x.start,y.start,x.end,y.end)],1,makestring)
	mins <- as.integer(c()); maxs <- as.integer(c()); tracts <- as.integer64(c())
	for(i in 1:nrow(avg)){
	  for(j in 1:length(temp[[i]])){
		mins <- as.integer(c(mins,as.integer(avg[i,as.integer(min)])))
		maxs <- as.integer(c(maxs,as.integer(avg[i,as.integer(max)])))
		tracts <- as.integer64(c(tracts,as.integer64(temp[[i]][j])))
	  }
	}
	overlaps <- data.table(min=as.integer(mins),max=as.integer(maxs),tracts=as.integer64(tracts))
	saveRDS(overlaps,paste0(data.path,"chicago_nores_tract_arc_overlaps.RDS")
} else{
	overlaps <- readRDS(paste0(data.path,"chicago_nores_tract_arc_overlaps.RDS"))
}

edges.out <- avg[,.(node1=min,node2=max,dist)][order(node1,node2)]
nodes.out <- nodes[,.(id=newid,area,x,y)][order(id)]
saveRDS(edges.out,paste0(data.path,"chicago_edges_nores.RDS"))
saveRDS(nodes.out,paste0(data.path,"chicago_nodes_nores.RDS"))
write.table(data.frame(n=nrow(nodes.out),m=nrow(edges.out)*2),
            file=paste0(data.path,"chicago_nores.edge"), 
            col.names=FALSE, row.names=FALSE, append=FALSE, sep=" ")
write.table(edges.out[,.(node1,node2,dm=round(dist*10))],
            file=paste0(data.path,"chicago_nores.edge"), 
            col.names=FALSE, row.names=FALSE, append=TRUE, sep=" ") #length in decimeters
write.table(nodes.out,
            file=paste0(data.path,"chicago_withores.co"), 
            col.names=FALSE, row.names=FALSE, append=FALSE, sep=" ")

train <- nodes[,.(long=x,lat=y)]
trips <- trips[!is.na(end.long) & !is.na(end.lat) & !is.na(start.long) & !is.na(start.lat)]
test <- unique(rbindlist(list(
  trips[,.(long=end.long,lat=end.lat)],
  trips[,.(long=start.long,lat=start.lat)])))
  
test$node <- as.numeric(FNN::knn(train, test, nodes$newid, k = 1))

trips[test,newnode1:=node, on=.(start.long=long,start.lat=lat)]
trips[test,newnode2:=node, on=.(end.long=long,end.lat=lat)]
trips[,second:=start.datetime-min(start.datetime)]
requests.out <- trips[,.(node1=newnode1,node2=newnode2,second)][order(second,node1,node2)]
requests.out[,second:=pmax(0,second+sample(-450:450, .N, replace=T))]
write.table(requests.out[,.(second,node1,node2)][order(second,node1,node2)],
            file=paste0(data.path,"requests_chicago_nores.csv"), col.names=FALSE, row.names=FALSE, append=FALSE, sep=","
)


if(!alreadyBuilt){
	week.to.use <- '2019-09-09'
	column.subset[,week.rank:=frank(week.date,ties.method="dense")]
	week.rank.to.use <- unique(column.subset[week.date==week.to.use,week.rank])
	forecast.subset <- column.subset[week.rank<week.rank.to.use & week.rank >= week.rank.to.use-4 & (!is.na(start.lat))]
	#fix missing start areas
	train <- nodes
	test <- forecast.subset[(is.na(start.area) | start.area==0)]
	test$new.area <- as.numeric(FNN::knn(train[,.(x,y)], test[,.(start.long,start.lat)], train[,as.numeric(area)], k = 1))
	forecast.subset[(is.na(start.area) | start.area==0),]$start.area <- test$new.area
	#get medoids, step 1: "snap" to a nearby tract centroid
	with.tract <- forecast.subset[!is.na(start.tract)]
	with.tract.unique <- unique(with.tract[,.(start.area,start.tract,start.long,start.lat)])
	index <- get.knnx(nodes[,.(x,y)], with.tract.unique[,.(start.long,start.lat)],1)$nn.index
	with.tract.unique$node <-  nodes[index,newid]
	with.tract.unique[nodes,':='(node.long=x,node.lat=y),on=.(node=newid)]
	with.tract[with.tract.unique,node:=node,on=.(start.tract)]
	with.tract[nodes,':='(node.long=x,node.lat=y),on=.(node=newid)]
	areas.to.center <- data.table(area=unique(with.tract$start.area))[order(area)]
	areas.to.center$medoid.node <- as.numeric(NA)
	for(i in 1:nrow(areas.to.center)){
		this.area <- with.tract[start.area==areas.to.center[i,area]]
		if(nrow(this.area) > 1000){
			this.area <- this.area[sample(.N,1000)]
		}
		idx.of.medoid <- pam(this.area[,.(start.long,start.lat)], k=1, metric = "manhattan", stand = FALSE)$id.med
		areas.to.center[i,]$medoid.node <- this.area[idx.of.medoid,node]
	}
	write.table(areas.to.center,
            file=paste0(data.path,"chicago_nores_area_medoids.txt"), 
            col.names=FALSE, row.names=FALSE, append=FALSE, sep=" ")
	saveRDS(areas.to.center,paste0(data.path,"chicago_nores_area_medoids.RDS")

	forecast.subset <- forecast.subset[,.N, by=.(start.datetime,start.area)]
	forecast.subset[,time.rank:=frank(start.datetime,ties.method="dense")]
	forecast.subset[,time.rank:=1+(time.rank-1)%%(24*7*4)]
	forecast.subset[,':='(wday=wday(start.datetime),time=as.ITime(start.datetime))]
	forecast.subset <- forecast.subset[,.(N=sum(N)),by=.(time.rank,wday,time,start.area)]
	forecasts <- data.table(crossing(
		unique(forecast.subset[,.(time.rank,wday,time)]),
		unique(forecast.subset[,.(start.area)])))
	forecasts[forecast.subset,N:=N,on=.(time.rank,wday,time,start.area)]
	forecasts[is.na(N),N:=0]
	forecasts[is.na(start.area),start.area:=0]
	forecasts[,time.rank:=(time.rank*3)-2]
	forecasts.stacked <- rbindlist(list(
		forecasts,
		forecasts[,.(time.rank=time.rank+1,wday=wday,time=time+300,start.area=start.area,N=N)],
		forecasts[,.(time.rank=time.rank+2,wday=wday,time=time+600,start.area=start.area,N=N)]		
	))
	setorder(forecasts.stacked,time.rank,wday,time,start.area)
	forecasts.stacked[,N:=N/3]
	forecasts.stacked[,':='(
		demand.5=N,
		demand.10=N +
		shift(N,n=1,fill=0,type="lead"),
		demand.15=N +
		shift(N,n=1,fill=0,type="lead") +
		shift(N,n=2,fill=0,type="lead"),
		demand.30=N +
		shift(N,n=1,fill=0,type="lead") +
		shift(N,n=2,fill=0,type="lead") +
		shift(N,n=3,fill=0,type="lead") +
		shift(N,n=4,fill=0,type="lead") +
		shift(N,n=5,fill=0,type="lead"))
		,by=.(start.area)]	
	
	saveRDS(forecasts.stacked, file = paste0(data.path,"chicago_forecasts_stacked.rds"))
	write.table(forecasts.stacked[,.(300*(time.rank-1),start.area,round(demand.5),round(demand.10),round(demand.15),round(demand.30))],
            file=paste0(data.path,"chicago_forecasts.txt"), 
            col.names=FALSE, row.names=FALSE, append=FALSE, sep=" ")
}

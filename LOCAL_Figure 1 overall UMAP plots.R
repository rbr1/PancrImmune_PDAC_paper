Draw_box_plot<-function(box,x,width,c,lwd,line_col){
	segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
	segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
	segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
	rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
	segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}
	
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha)) }

concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}

######################
batch = "PDAC150Ka"
type = "all"
analysis = "Overall_UMAP_plots"
names = c("all","B_cell","T_cell","Myeloid_cells")
output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"

library(yarrr)
library(RColorBrewer)
colx = c(piratepal(palette = "info2"), piratepal(palette = "basel"), piratepal(palette = "google"),piratepal(palette = "espresso"),brewer.pal(8, "Dark2"))
colx = c(piratepal(palette = "google"),"grey", piratepal(palette = "espresso")[c(2,4,5)],piratepal(palette = "info2"),  brewer.pal(8, "Dark2"), piratepal(palette = "basel"))
colx = unique(colx)
cols =  add.alpha (colx, alpha = 0.95)
cols1 = add.alpha (cols,alpha = 0.5)
cols2 = add.alpha (cols, alpha = 0.5)

for(c in c(1:length(names))){
  #colx = c(piratepal(palette = "info2"), piratepal(palette = "basel"), piratepal(palette = "google"),piratepal(palette = "espresso"),brewer.pal(8, "Dark2"))
  colx = c(piratepal(palette = "google"),"grey", piratepal(palette = "espresso")[c(2,4,5)],piratepal(palette = "info2"),  brewer.pal(8, "Dark2"), piratepal(palette = "basel"))
  colx = unique(colx)
  cols =  add.alpha (colx, alpha = 0.95)
  cols1 = add.alpha (cols,alpha = 0.5)
  cols2 = add.alpha (cols, alpha = 0.5)
  
	file = concat(c(input_directory,"Overall_UMAP_annotations_PDAC150Ka_", names[c],".txt"))
	p <- as.matrix(read.csv(file, head=T, sep="\t"))
	p=p[which(p[,3]!="Acinar/Ductal"),]
	p=p[which(p[,3]!="-"),]
	xy = cbind(as.numeric(p[,1]),as.numeric(p[,2]))
	cluster = p[,3]
	print(names[c])
	print(length(cluster))
	
	if(c == 1){
		cluster[grep("B cell", cluster)]= "B cell"
		cluster[grep("T cell", cluster)]= "T cell"
		cluster[grep("NK", cluster)]= "NK cell"
		cluster[grep("Myeloid", cluster)]= "Myeloid"
		cluster[grep("ILC", cluster)]= "ILC"
	}
	print(concat(c(names[c]," ", length(xy[,1]))))
	clusters = sort(unique(cluster))
	if(c == 4){
		w = grep(" b", cluster, invert = T)
		xy = xy[w,]
		cluster = cluster[w]
		clusters = sort(unique(cluster))
		groups = NULL		
		broad_cell_type= clusters
		broad_cell_type[grep("DC", broad_cell_type)] = "DC"
		broad_cell_type[grep("momac", broad_cell_type)] = "momac"
		broad_cell_type[grep("mono", broad_cell_type)] = "monocyte"
		broad_cell_type[grep("ILC", broad_cell_type)] = "ILC"
		broad_cell_type[grep("Mast", broad_cell_type)] = "Mast"
		broad_cell_type[grep("granul", broad_cell_type)] = "granul"
		broad_cell_type[which(broad_cell_type=="megakaryocyte b")] = "megakaryocyte"
		broad_cell_types = unique(broad_cell_type)
		for(i in c(1:length(broad_cell_types))){
			w = which(broad_cell_type==broad_cell_types[i])
			groups = c(groups, list(clusters[w]))
		}
		clusters %in% unlist(groups)
		col1 = c(piratepal(palette = "google"),"grey", piratepal(palette = "espresso")[c(2,4,5)],piratepal(palette = "info2"),  brewer.pal(8, "Dark2"), piratepal(palette = "basel"))
		col1 = col1[c(1,3,2,4,5,7,6,8,9,10)]
		plot(c(1:length(col1)), c(1:length(col1)), bg= col1,col = col1, pch = 21)
		
		clusters = unlist(groups)
		cols_plot = rep("", length(clusters))
		names(cols_plot) = clusters
		for(i in c(1:length(groups))){
			w = match(groups[[i]],clusters)
			co = col1[i]
			seq = ceiling(seq(from = 30, to = 100, length = length(groups[[i]])))
			cols=colorRampPalette(c('white', co))(100)
			cols_plot [w] = cols[seq]
		}
	}else{
		if(c == 3){
			w = grep(" b", cluster, invert = T)
			xy = xy[w,]
			cluster = cluster[w]
			clusters = sort(unique(cluster))
			groups = NULL		
			broad_cell_type= clusters
			broad_cell_type[grep("CD4", clusters)] = "CD4"
			broad_cell_type[grep("Treg", clusters)] = "Treg"
			broad_cell_type[grep("Th", clusters)] = "Th"
			broad_cell_type[grep("Tfh", clusters)] = "Th"
			broad_cell_type[grep("CD8", clusters)] = "CD8"
			broad_cell_type[grep("ILC", clusters)] = "ILC"
			broad_cell_type[grep("NK", clusters)] = "NK"
			broad_cell_type[grep("iNKT", clusters)] = "Non-conventional"
			broad_cell_type[grep("gdT", clusters)] = "Non-conventional"
			broad_cell_type[grep("MAIT", clusters)] = "Non-conventional"
			broad_cell_types = unique(broad_cell_type)
			for(i in c(1:length(broad_cell_types))){
				w = which(broad_cell_type==broad_cell_types[i])
				groups = c(groups, list(clusters[w]))
			}
			clusters %in% unlist(groups)
			col1 = c(piratepal(palette = "info2")[c(1:3)], piratepal(palette = "google")[c(1:4)],piratepal(palette = "espresso"),brewer.pal(8, "Dark2"))
			col1 = c(piratepal(palette = "google"),add.alpha("black", alpha = 0.75), piratepal(palette = "espresso")[c(2,4,5)],piratepal(palette = "info2"),  brewer.pal(8, "Dark2"), piratepal(palette = "basel"))
			col1 = col1[c(1, 2, 3, 4, 5, 8, 6, 7)]
			col1[6] = "purple"
			plot(c(1:length(col1)), c(1:length(col1)), bg= col1,col = col1, pch = 21)
			  
			clusters = unlist(groups)
			cols_plot = rep("", length(clusters))
			names(cols_plot) = clusters
			for(i in c(1:length(groups))){
				w = match(groups[[i]],clusters)
				co = col1[i]
				seq = ceiling(seq(from = 30, to = 100, length = length(groups[[i]])))
				cols=colorRampPalette(c('white', co))(100)
				cols_plot [w] = cols[seq]
			}
		}else{
			if(c == 2){
				clusters = sort(unique(cluster))
				groups = NULL		
				broad_cell_type= clusters
				broad_cell_type[grep("plasma", clusters)] = "ASC"
				broad_cell_type[grep("activated", clusters)] = "activated"
				broad_cell_type[grep("memory", clusters)] = "memory"
				broad_cell_type[grep("GC", clusters)] = "GC"
				broad_cell_type[grep("activated pre-memory", clusters)] = "activated"
				broad_cell_types = unique(broad_cell_type)
				for(i in c(1:length(broad_cell_types))){
					w = which(broad_cell_type==broad_cell_types[i])
					groups = c(groups, list(clusters[w]))
				}
				clusters %in% unlist(groups)
				col1 = c(piratepal(palette = "info2")[c(1:3)], piratepal(palette = "google")[c(1:4)],piratepal(palette = "espresso"),brewer.pal(8, "Dark2"))
				clusters = unlist(groups)
				cols_plot = rep("", length(clusters))
				names(cols_plot) = clusters
				for(i in c(1:length(groups))){
					w = match(groups[[i]],clusters)
					co = col1[i]
					seq = ceiling(seq(from = 30, to = 100, length = length(groups[[i]])))
					if(length(groups[[i]])==2){seq = ceiling(seq(from = 60, to = 100, length = length(groups[[i]])))}
					if(length(groups[[i]])==1){seq = ceiling(seq(from = 90, to = 100, length = length(groups[[i]])))}
					cols=colorRampPalette(c('white', co))(100)
					cols_plot [w] = cols[seq]
				}
			}else{
				cols_plot = cols
				clusters = sort(unique(cluster))}}}
	
	fileout1=concat(c(output_directory, analysis,"_", batch,"_", names[c],".pdf"))
	w=3.45
	pdf(file=fileout1, height=w*1.1*2, width=w*1.1)
	par(mfrow= c(2,1), mar = c(5,5,3,3))

	cluster_match = match(cluster, clusters)
	pches = rep(c(21:25), 100)
	plot(xy[,1], xy[,2],pch = pches[cluster_match], col = cols_plot[cluster_match],bg = cols_plot[cluster_match] ,main = "original UMAP clustering", cex = 0.1,xlab = "umap dim1",ylab = "umap dim2")
	plot(c(0,1),c(0,1), col = "white" ,main = "", cex = 0.25,axes = F, col.axis = "white",xlab = '',ylab = '')
	clusters1 = clusters
	if(length(clusters1)<25){
		clusters1 = gsub("_"," ", clusters1)
		legend("top", clusters1, pch = pches,cex= 0.4, bty="n", pt.bg = cols_plot, col = cols_plot)
	}else{
		seq1 = c(1:16)
		seq2 = c(17:40)
		clusters1 = gsub("_"," ", clusters1)
		legend("topleft", clusters1[seq1], pch = pches,cex= 0.4, bty="n", pt.bg = cols_plot[seq1], col = cols_plot[seq1])
		legend("topright", clusters1[seq2], pch = pches,cex= 0.4, bty="n", pt.bg = cols_plot[seq2], col = cols_plot[seq2])
	}
	dev.off()
	print(c)
}



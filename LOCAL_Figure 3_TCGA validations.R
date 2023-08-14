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
batch = "TCGA_PAAD"
analysis = "TCGA"

######################
## get patient groups for plotting: 
input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"

input_directory_groups  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"

library(yarrr)
library(RColorBrewer)
library(BayesPrism)
load("~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/TCGA/PDAC_filtered_cell_fractions_Level_1.rda")


## group proportions
file= "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/Cell_annotations_PDAC150Ka_levels.txt"
p <- as.matrix(read.csv(file, head=T, sep="\t"))
## check
which(colnames(subset_theta) %in% p[,6] ==F)

p=p[which(is.na(p[,6])==F),]
cell_type_level0 =trimws(p[,6], "r")
cell_type_level1 = apply(cbind(trimws(p[,4], "r"), ""), 1, paste, collapse = "")
cell_type_level2 = apply(cbind(trimws(p[,5], "r"), ""), 1, paste, collapse = "")
names(cell_type_level1) = cell_type_level0
names(cell_type_level2) = cell_type_level0

cells_type_level0 = sort(unique(cell_type_level0))
cells_type_level1 = sort(unique(cell_type_level1))
cells_type_level2 = sort(unique(cell_type_level2))

list_proportions = NULL
names = NULL
for(c in c(1:length(cells_type_level1))){
  w = which(cell_type_level1==cells_type_level1[c])
  cell_type_use = sort(unique(cell_type_level0[w]))
  if(length(cell_type_use)>=2){
    list_proportions = c(list_proportions, list(subset_theta[,cell_type_use]))
    names = c(names, cells_type_level1[c])}}
for(c in c(1:length(cells_type_level2))){
  w = which(cell_type_level2==cells_type_level2[c])
  cell_type_use = sort(unique(cell_type_level0[w]))
  if(length(cell_type_use)>=2){
    list_proportions = c(list_proportions, list(subset_theta[,cell_type_use]))
    names = c(names, cells_type_level2[c])}}


names(list_proportions) = names

a = cbind(rowSums(list_proportions[["B cell" ]]), rowSums(list_proportions[["Myeloid" ]]), 
          rowSums(list_proportions[["NK" ]]), rowSums(list_proportions[["T cell" ]]))
colnames(a) = c("B cell", "Myeloid", "NK", "T cell")

list_proportions = c(list_proportions, list(a))
names(list_proportions) = c(names, "all")

## normalise per cell group
list_proportions_norm = NULL
for(c in c(1:length(list_proportions))){
  m = list_proportions[[c]]
  for(i in c(1:length(m[,1]))){m[i,] = m[i,]*100/sum(m[i,])}
  list_proportions_norm = c(list_proportions_norm, list(m))
}
names(list_proportions_norm) = names(list_proportions)

## get corr plots all against all
Plot_all_correlations<-function(list_proportions_norm, output_directory,batch){
  # plot correlations
  library(RColorBrewer)
  library(psych)
  library(corrplot)
  mat_all = NULL
  for(c in c(1:length(list_proportions_norm))){
    m = list_proportions_norm[[c]]
    colnames(m) = apply(cbind(colnames(m), "_",names(list_proportions_norm)[c]), 1, paste, collapse = "")
    if(length(mat_all)==0){mat_all = m
    }else{
      mat_all = cbind(mat_all , m)
    }}
  
  Treg_total_T_cell = mat_all[,"T cell CD4 Treg_T cell"]+mat_all[,"T cell CD4 Treg Activated_T cell"]
  mat_all = cbind(mat_all,Treg_total_T_cell)
  
  cortest = corr.test(mat_all,adjust="holm")
  pval = cortest $ p
  rval = cortest $ r
  fileout1=concat(c(output_directory,"Correlation_cell_types_pairs_plot_biopsy_overall_", batch,".pdf"))
  w=30
  pdf(file=fileout1, height=w*1.1, width=w*1.1)
  par(mfrow= c(1,1), mar = c(5,5,3,3))
  corrplot(rval, type="upper", p.mat=pval, insig="label_sig", tl.pos="td", sig.level=0.05, title ="", method = "ellipse", diag = F, order = "hclust",, cl.pos = 'n')
  dev.off()
  
  #### specific plots
 
  from = c("Myeloid_all", "Myeloid_all", "Myeloid_all", "T cell_all", "T cell_all")
  to = c("B cell  plasma cell_B cell","T cell CD4 Treg_T cell", "Treg_total_T_cell", "T cell CD8 EM_T cell","Myeloid_all")
  fileout1=concat(c(output_directory,"Correlation_cell_types_pairs_plot_biopsy_overall_", batch,"_2.pdf"))
  w=4
  pdf(file=fileout1, height=w*1*3, width=w*1*2)
  par(mfrow= c(3,2), mar = c(5,5,3,3))
  for(c in c(1:length(from))){
    plot(mat_all[,from[c]], mat_all[,to[c]], ylab = to[c], xlab =from[c],main = batch)
  }
  dev.off()
  
  ## split by high mid low myeloid
  x =  mat_all[,"Myeloid_all"]
  groups_id = c(list(which(x<quantile(x, 0.33))), list(intersect(which(x>=quantile(x, 0.33)),which(x<=quantile(x, 0.66)))), 
                     list(which(x>quantile(x, 0.66))))
  names(groups_id) = c("low % myeloid", "mid % myeloid", "high % myeloid")
  groups_id = rev(groups_id)
  
  plot_cell_types = sort(unique(c(from, to)))
  fileout1=concat(c(output_directory,"Correlation_cell_types_pairs_plot_biopsy_overall_", batch,"_3.pdf"))
  w=3
  pdf(file=fileout1, height=w*1*3, width=w*0.6*3)
  par(mfrow= c(3,3), mar = c(8,5,3,3))
  log = c("Treg_total_T_cell"  )
  for(c in c(1:length(plot_cell_types))){
    groups = NULL
    for(i in c(1:length(groups_id))){
      if(plot_cell_types[c] %in% log){groups = c(groups, list(log10(mat_all[groups_id[[i]],plot_cell_types[c]])))
      }else{groups = c(groups, list(mat_all[groups_id[[i]],plot_cell_types[c]]))}
    }
    names(groups) = names(groups_id)
    factors = names(groups) 
    main = plot_cell_types[c]
    max = max(c(unlist(groups), unlist(groups)))
    min = 0
    if(plot_cell_types[c] %in% log){min = min(c(unlist(groups), unlist(groups))*1.05)}
    b = (max-min)*0.05
    max = max+15*b
    ylab = "%"
    if(plot_cell_types[c] %in% log){ylab = "log10(%)"}
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 3
    max_scale = min(c(100, max))
    range = max-min
    if(range>50){scale = c(-100:100)*20}
    if(range<=50){scale = c(-100:100)*5}
    if(range <10){scale = c(-100:100)*2.5}
    if(range <5){scale = c(-100:100)*1}
    if(range <4){scale = c(-100:100)*0.5}
    if(range <1.5){scale = c(-100:1000)*0.2}
    if(range <0.5){scale = c(-100:100)*0.1}
    if(range <0.1){scale = c(-100:100)*0.01}
    if(range <0.01){scale = c(-100:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.36
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.2/max(shift)
    
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,4,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(i in c(1:l)){
        points1=as.numeric(groups[[i]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i,width,cols[i],1, cols1[i])
        points(rep(i, length(points1)),points1, pch =21, bg = cols[i], col=cols1[i], cex = 0.7)
      }
    
    y = max(unlist(groups))+b
    for(i1 in c(1:l)){	
      for(i2 in c(i1:l)){	
        if(i1<i2){
          p = wilcox.test(groups[[i1]],y=groups[[i2]])$p.value
          pval1 = signif(p, digits = 3)
          segments(i1, y, i2, y, lwd = 2, col = "grey")
          text(mean(c(i1,i2)), y+1.5*b, labels = pval1, cex = 1)
          y = y+3*b
        }}}
    }
  dev.off()
}
Plot_all_correlations(list_proportions_norm, output_directory,batch)

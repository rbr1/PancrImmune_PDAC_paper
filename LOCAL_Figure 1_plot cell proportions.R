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
analysis = "PCA_cell_counts"
analysis_all = "PCA_cell_counts"
######################

output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
input_directory  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"


Get_cell_counts_overall<-function(input_directory, output_directory, batch){
	files = concat(c(input_directory,"Overall_all_cell_counts_PDAC150Ka.txt"))
	file_type = c("CD45")
	list_counts = NULL
	for(f in c(1:length(files))){
		p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
		cell_type = p[,1]
		p=p[,c(2:length(p[1,]))]
		rownames(p) = cell_type
		p=p[which(rownames(p)!='-'),]
		sample = colnames(p)
		pat_sample = gsub("_CD45p1","", sample)
		pat_sample = gsub("_CD45p2","", pat_sample)
		pat_sample = gsub("X","", pat_sample)
		pat_sample = gsub("01_","", pat_sample)
		pat_sample = gsub("02_","", pat_sample)
		pat_sample = gsub("03_","", pat_sample)
		pat_sample = gsub("04_","", pat_sample)
		pat_sample = gsub("05_","", pat_sample)
		pat_sample = gsub("06_","", pat_sample)
		pat_sample = gsub("07_","", pat_sample)
		cell_type = rownames(p)
		pat_samples = sort(unique(pat_sample))
		p1 = matrix(data = 0,nrow = length(cell_type), ncol = length(pat_samples), dimnames = c(list(cell_type), list(pat_samples)))
		for(i in c(1:length(pat_samples))){
			w = which(pat_sample== pat_samples[i])
			for(j in w){p1[, pat_samples[i]] = p1[, pat_samples[i]]+as.numeric(p[,j])}
		}
		p1 = p1[order(rownames(p1)),]
		list_counts = c(list_counts, list(p1))
	}
	list_names = file_type
	for(f in c(1:length(files))){
		p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
		p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
		cell_type = p[,1]
		p=p[,c(2:length(p[1,]))]
		rownames(p) = cell_type
		p=p[which(rownames(p)!='-'),]
		sample = colnames(p)
		pat_sample = gsub("_CD45p1","", sample)
		pat_sample = gsub("_CD45p2","", pat_sample)
		pat_sample = gsub("X","", pat_sample)
		pat_sample = gsub("01_","", pat_sample)
		pat_sample = gsub("02_","", pat_sample)
		pat_sample = gsub("03_","", pat_sample)
		pat_sample = gsub("04_","", pat_sample)
		pat_sample = gsub("05_","", pat_sample)
		pat_sample = gsub("06_","", pat_sample)
		pat_sample = gsub("07_","", pat_sample)
		cell_type = rownames(p)
		pat_samples = sort(unique(pat_sample))
		p1 = matrix(data = 0,nrow = length(cell_type), ncol = length(pat_samples), dimnames = c(list(cell_type), list(pat_samples)))
		for(i in c(1:length(pat_samples))){
			w = which(pat_sample== pat_samples[i])
			for(j in w){p1[, pat_samples[i]] = p1[, pat_samples[i]]+as.numeric(p[,j])}
		}
		p1 = p1[order(rownames(p1)),]
		cell_type = rownames(p)
		group = cell_type
		group [grep("^NK", cell_type)] = "NK"
		group [grep("^Myeloid", cell_type)] = "Myeloid"
		group [grep("^B cell", cell_type)] = "B cell"
		group [grep("^T cell", cell_type)] = "T cell"
		group [grep("^ILC", cell_type)] = "ILC"
		groups = sort(unique(group))
		list_group = NULL
		nam = NULL
		for(i in c(1:length(groups))){
			w = cell_type [which(group==groups[i])]
			if(length(w)>1){
				list_group = c(list_group, list(w))
				list_counts = c(list_counts, list(p1[w,]))
				list_names = c(list_names, concat(c(file_type[f],":",groups[i] )))
			}
		}
		names(list_group) = groups
		group_count = NULL
		for(i in c(1:length(groups))){
			w = cell_type [which(group==groups[i])]
			x = p[w,]
			if(length(w)>1){
				x = colSums(p1[w,])
			}
			if(length(group_count)==0){group_count = x
			}else{group_count = rbind(group_count, x)}
		}
		rownames(group_count) = groups
		list_counts = c(list_counts, list(group_count))
		list_names = c(list_names, concat(c(file_type[f],":","broad" )))
	}
	names(list_counts) = list_names
	##### print out
	out = NULL
	for(i in c(1:length(list_counts))){
		x = cbind(rownames(list_counts[[i]]), list_names[i], list_counts[[i]])
		if(length(out)==0){out = x
		}else{out = rbind(out, x)}
	}
	colnames(out)[1] = "Cell_type"
	colnames(out)[2] = "Group"	
	out_file_table = concat(c(output_directory,"Seurat_cell_counts_",batch,"_all.txt"))
	write.table(out, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}
Get_cell_counts_overall(input_directory, output_directory, batch)

file = concat(c(output_directory,"Seurat_cell_counts_",batch,"_all.txt"))
p <- as.matrix(read.csv(file, head=T, sep="\t"))
p_blood = unique(colnames(p)[c(grep("blood",colnames(p)), grep("HC",colnames(p)), grep("pbmc",colnames(p)))])
p_blood = p[,p_blood]

sample = colnames(p)[-1]
sample = sample[-1]
sites = sample
disease = sample
sites[grep("biopsy", sites)] = "PDAC tissue"
sites[grep("blood", sites)] = "PDAC PBMC"
names(sites) = sample

cell_types = p[,1]
cell_group = p[,2]
cell_groups = sort(unique(cell_group))
cell_groups = cell_groups[which(cell_groups!="Overall")]
list_normalised = NULL

include = NULL
for(ind in c(1:length(cell_groups))){
	w = which(cell_group== cell_groups[ind])
	if(cell_groups[ind] =="CD45:broad"){
		w = intersect(w, which(p[,1] %in% c( "B cell"  ,"ILC", "Myeloid","NK" , "T cell"   )))
	}
	cell_type_sub = p[w,1]
	p1 = p[w,]
	m_count_norm = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(sample), dimnames = c(list(cell_type_sub), list(sample)))
	
	for(s in c(1:length(sample))){
		m_count_norm[, sample[s]]= as.numeric(p1[, sample[s]])*100/sum(as.numeric(p1[, sample[s]]))
		if(sum(as.numeric(p1[, sample[s]]))<7){m_count_norm[, sample[s]] = -1}
	}
	if(length(which(apply(m_count_norm, 2, max)!=-1))>3){include = c(include, ind)}
	m_count_norm= m_count_norm[order(rownames(m_count_norm)),]
	list_normalised= c(list_normalised, list(t(m_count_norm)))
}
names(list_normalised) = cell_groups
list_normalised = list_normalised[include]


for(i in c(1:length(list_normalised))){
	print (names(list_normalised)[i])
	print (colnames(list_normalised[[i]]))
}

x = list_normalised[["CD45:broad"]]
B_T_cell = rowSums(x[,c("B cell", "T cell")])
x = cbind(x[,c("Myeloid", "NK"  )], B_T_cell)
colnames(x)[3] = "B/T cells"

y = list_normalised[["CD45:B cell" ]]
broad = c("plasma", " memory","activated pre-memory", "GC", "MZ", "naive")
y1 =NULL
for(i in c(1:length(broad))){
	w = colnames(y)[grep(broad[i], colnames(y))]
	if(length(w)>1){yi = rowSums(y[,w])
	}else{yi =(y[,w])}
	if(length(y1)==0){y1 = yi
	}else{y1=cbind(y1, yi)}
}
colnames(y1) = gsub(" ","", gsub(" pre-memory","", gsub("plasma","ASC", broad, fixed = T), fixed = T),fixed = T)
y = y1

z = list_normalised[["CD45:T cell" ]]
broad = c("CD4", "CD8","MAIT", "gd")
z1 =NULL
for(i in c(1:length(broad))){
	w = colnames(z)[grep(broad[i], colnames(z))]
	if(length(w)>1){yi = rowSums(z[,w])
	}else{yi =(z[,w])}
	if(length(z1)==0){z1 = yi
	}else{z1=cbind(z1, yi)}
}
colnames(z1) = gsub("gd","gamma/delta", broad, fixed = T)
z = z1

s = list_normalised[["CD45:Myeloid" ]]
broad = c("DC", "granulo","Mast", "megakaryocyte","momac","monocyte")
s1 =NULL
for(i in c(1:length(broad))){
	w = colnames(s)[grep(broad[i], colnames(s))]
	print(length(w))
	if(length(w)>1){yi = rowSums(s[,w])
	}else{yi =(s[,w])}
	if(length(s1)==0){s1 = yi
	}else{s1=cbind(s1, yi)}
}
colnames(s1) = gsub("granulo","granulocyte", broad, fixed = T)
s = s1

l = length(list_normalised)
list_normalised1 = c(list_normalised, list(x), list(y), list(z), list(s))
names(list_normalised1)[c((l+1):length(list_normalised1))] = c("CD45:broad condensed", "CD45:broad B cell","CD45:broad T cell", "CD45:broad Myeloid")

list_normalised = list_normalised1

counts = list_normalised[["CD45:T cell"]]
CD8 = c(grep("CD8",colnames(counts)),grep("MAIT",colnames(counts)),grep("gdT",colnames(counts)),grep("NK",colnames(counts)))
CD4 = c(grep("CD4",colnames(counts)))
CD8 = counts[,CD8]
CD4 = counts[,CD4]
for(i in c(1:length(counts[,1]))){
  CD8[i,] = CD8[i,]*100/sum(CD8[i,])
  CD4[i,] = CD4[i,]*100/sum(CD4[i,])
}

list_normalised = c(list_normalised, list(CD4),list(CD8))
names(list_normalised) = c(names(list_normalised1), "CD4","CD8")

for(i in c(1:length(list_normalised))){
	print (names(list_normalised)[i])
	print (colnames(list_normalised[[i]]))
}

T_cells = list_normalised[["CD45:T cell"]]
Treg_CD8_ratio = rowSums(T_cells [,grep("Treg",colnames(T_cells))])/rowSums(T_cells [,grep("CD8",colnames(T_cells))])

T_cells = list_normalised[["CD45"]]
CD4_CD8_ratio = rowSums(T_cells [,grep("CD4",colnames(T_cells))])/rowSums(T_cells [,grep("CD8",colnames(T_cells))])

cd4cd8 = c(list(Treg_CD8_ratio), list(CD4_CD8_ratio))
names(cd4cd8) = c("Treg-CD8 T cell ratio", "CD4-CD8 ratio")

T_cells = list_normalised[["CD45:T cell"]]
Treg_act_ratio = T_cells [,"T cell CD4 Treg Activated"]/rowSums(T_cells [,grep("Treg",colnames(T_cells))])

T_ratios = cbind(Treg_CD8_ratio, CD4_CD8_ratio, Treg_act_ratio)

saveRDS(file = concat(c(output_directory, "Cell_counts_normalised_PDAC150Ka.RDS")), list_normalised)
############################################ split patients into groups
all_proportions = list_normalised[["CD45:broad"]]
ws = rownames(all_proportions)[grep("biopsy", rownames(all_proportions))]
group1 = names(which(all_proportions[ws,"Myeloid"] > quantile(all_proportions[ws,"Myeloid"], 0.6)))
group2 = names(which(all_proportions[ws,"Myeloid"] <= quantile(all_proportions[ws,"Myeloid"], 0.6)))

groups_PCA = c(list(group1), list(group2))
saveRDS(file = concat(c(output_directory, "Cell_counts_normalised_PDAC150Ka.RDS")), groups_PCA)

############################################
groups_PCA = readRDS(file = concat(c(output_directory, "Cell_counts_normalised_PDAC150Ka.RDS")))
group1 = groups_PCA[[1]]
group2 = groups_PCA[[2]]

############################################ plot PCA of CD45 %s
Plot_PCA<-function(list_normalised, output_directory){
  all_proportions = (list_normalised[["CD45:broad" ]])
  type_to_plot = "PDAC tissue"
  
  ws = names(which(sites == type_to_plot))
  x <- princomp(all_proportions[ws,],dim=2)
  PCA_coordinates_sub = x$scores
  
  x1 = x$scores[,1]
  x2 = x$scores[,2]
  
  y1 = x$loadings[,1]
  y2 = x$loadings[,2]
  
  factor1 = rep(0,length(x1))
  names(factor1) = names(x1)
  factor1[group1] = 1
  factor1[group2] = 2
  factor2 = gsub("-biopsy","",names(x1))
  
  factors1 = sort(unique(factor1))
  factors2 = sort(unique(factor2))
  matching1 = match(factor1, factors1)
  matching2 = match(factor2, factors2)
  
  library(RColorBrewer)
  cols =  add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95)
  cols1 = add.alpha (c( cols[1],"white"),alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)
  
  pches = c(21:25)
  fileout1=concat(c(output_directory,"", analysis_all,"_biopsy_PCA_scores_", batch,".pdf"))
  w=2.2
  pdf(file=fileout1, height=w*1, width=w*3)
  par(mfrow= c(1,3), mar = c(5,5,3,3))
  
  plot(x1,x2,pch = pches[matching1], col = cols[matching1],bg = cols2[matching1], main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  
  xy = cbind(x1,x2)
  library(car)
  for(pca in c(1:length(groups_PCA))){
    w = which(factor1 ==pca)
    if(length(w)>1){
      x = xy[w,1]
      y = xy[w,2]
      dataEllipse(x, y, levels=c(0.65), add = TRUE, col = cols2[pca], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
    }
  }
  
  plot(x1,x2,pch = pches[matching1], col = "white",bg = "white", main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  legend("topleft", paste("group",factors1), pch = pches,cex= 0.8, bty="n", pt.bg = cols2, col = cols, pt.lwd = 2, text.font = 2)
  
  plot(x1,x2,pch = pches[matching1], col = cols[matching1],bg = cols2[matching1], main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  text(x1, y = x2, labels = factor2,cex = 0.4,font = 1)
  
  dev.off()
}
Plot_PCA(list_normalised, output_directory)
############################################ plot boxplots of cell %s in tumour between groups
Plot_tumour_cell_composition_between_groups<-function(list_normalised, output_directory){
  fileout1=concat(c(output_directory,"", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_", batch,".pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*6, width=w*3*2.9)
  par(mfrow= c(6,2), mar = c(16,5,2,3))
  
  summary_tables = NULL
  analysis1 = "cell_proportions"
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  for(c in c(1:length(list_normalised))){
    mat = list_normalised[[c]]
    name = names(list_normalised)[c]
    analysis = concat(c(analysis1,"_", name))	
    mat_stat = mat
    for(i in c(1:length(mat[1,]))){
      w = which(mat[,i]==-1)
      mat_stat[w,i] =NA
    }
    
    mat_biopsy = mat_stat[unlist(groups_PCA),]
    mat_blood = mat_stat[gsub("biopsy","blood",unlist(groups_PCA)),]
    
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    mat_biopsy = mat_biopsy[,which(colSums(mat_biopsy)>0)]
    ##### biopsy
    mat1 = mat_biopsy[group_PCA_list,]
    factor = factor(factor_PCA)
    fit = manova(formula = mat1^0.5 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
    Means_factor = function(factor, x){
      m = NULL
      sd = NULL
      for(i1 in c(1:length(levels(factor)))){
        x1 = x[which(factor==levels(factor)[i1])]
        x1 = x1[which(x1!=-1)]
        sd = c(sd, sd(x1))
        m = c(m, mean(x1))}
      return(c(m, sd))}
    i1 = 0
    for(i in p1){
      i1 = i1+1
      p_value = c(p_value, i$'Pr(>F)'[1]) 
      if(length(mean)==0){means = Means_factor(factor, mat1[,i1])
      }else{means = rbind(means, Means_factor(factor, mat1[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    colnames(means) = c(paste("mean.group.", c(1:length(levels(factor)))), paste("sd.group.", c(1:length(levels(factor)))))
    combined_p_value = cbind(p_value ,means)
    rownames(combined_p_value) = nam
    p.site = rep("biopsy", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    if(length(summary_tables)==0){summary_tables = x
    }else{summary_tables = rbind(summary_tables ,x)}
    ##### blood
    mat1 = mat_blood
    factor = factor(factor_PCA)
    fit = manova(formula = mat1 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
    i1 = 0
    for(i in p1){
      i1 = i1+1
      p_value = c(p_value, i$'Pr(>F)'[1]) 
      if(length(mean)==0){means = Means_factor(factor, mat1[,i1])
      }else{means = rbind(means, Means_factor(factor, mat1[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    colnames(means) = c(paste("mean.group.", c(1:length(levels(factor)))), paste("sd.group.", c(1:length(levels(factor)))))
    combined_p_value_blood = cbind(p_value ,means)
    rownames(combined_p_value_blood) = nam
    p.site = rep("blood", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value_blood)
    if(length(summary_tables)==0){summary_tables = x
    }else{summary_tables = rbind(summary_tables ,x)}
    
    groups = NULL
    mat = mat_biopsy[,which(colSums(mat_biopsy)>0)]
    p_values = combined_p_value[,1]
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA))){
        data = c(data, mat [groups_PCA[[g]], i])
        factor = c(factor, rep(g, length(mat [groups_PCA[[g]], i])))
        g1 = c(g1, list(mat [groups_PCA[[g]], i]))}
      groups = c(groups, list(g1))
      # group = factor(factor)
      # fit = aov(formula = data ~ group)
      # pvalues = summary.lm(fit)$ coefficients[,4]
      # pVal <- anova(fit)$'Pr(>F)'[1]
      # p_values = c(p_values, pVal)
    }
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = gsub("_"," ",colnames(mat))
    main = concat(c("% of ",name))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 70
    max_scale = min(c(max, 100))
    range = max-min
    if(range>55){scale = c(0:100)*20}
    if(range<=55){scale = c(0:100)*10}
    if(range<=40){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(2.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.18
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.25/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg = cols[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        # if(p_values[i] <signif_threshold/10){pval1 = "**"}
        # if(p_values[i] <signif_threshold/100){pval1 = "***"}
        y = max(unlist(groups[[i]]))
        y = y+1*b
        # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
        text(i, y+2*b, labels = pval1, cex = 1.7)
      }
    }
  }
  
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
  legend("topright", factors1, pch = 22,cex= 0.9, bty="n", pt.bg = add.alpha("white",alpha =0), col = add.alpha("white",alpha =0), text.col = cols1, text.font = 2)
  dev.off()
  return(summary_tables)
}
summary_tables = Plot_tumour_cell_composition_between_groups(list_normalised, output_directory)

outfile = concat(c(output_directory,"Stats_PCA_cell_types_by_groups_PDAC150Ka.txt"))
write.table(summary_tables, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")

############################################ plot FC versus p-values of cell %s in tumour between groups
Plot_FC_p_value<-function(summary_tables, output_directory){
  summary_tables[which(as.numeric(summary_tables[,"p_value"])<0.05),]
  FC = log10(as.numeric(summary_tables[,"mean.group. 1"])/as.numeric(summary_tables[,"mean.group. 2"]))
  pval = -log10(as.numeric(summary_tables[,"p_value"]))
  plot(FC, pval)
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  types_include = c("cell_proportions_CD45:T cell", "cell_proportions_CD45:NK", "cell_proportions_CD45:ILC", "cell_proportions_CD45:broad T cell","cell_proportions_CD45:B cell")
  w = intersect(which(summary_tables[,"p.site"]=="biopsy"), which(summary_tables[,"p.analysis"] %in% types_include))
  pval = as.numeric(summary_tables[w,"p_value"])
  
  names =rownames(summary_tables[w])
  names = cbind(names, summary_tables[w,"p.analysis"])
  names[,1] = gsub("cell_proportions_CD45:", "% of ",names[,1])
  names[,1] = gsub("broad T cell", "T/NK/ILC cell",names[,1])
  rownames(names) = gsub("  "," ",rownames(names))
  rownames(names) = gsub("T cell ","",rownames(names))
  rownames(names) = gsub("NK cell ","",rownames(names))
  rownames(names) = gsub("B cell ","",rownames(names))
  rownames(names) = gsub("Myeloid ","",rownames(names))
  names = apply(cbind(rownames(names), " ",names,"s"), 1, paste, collapse = "")
  names = gsub("gamma/delta","gdT",names)
  cbind(names)
  
  FC = as.numeric(summary_tables[w,"mean.group. 2"]) - as.numeric(summary_tables[w,"mean.group. 1"])
  x = FC
  y = -log10(pval)
  names(x) = names
  names(y) = names
  b = max(y)*0.1
  range = max(abs(FC))*5
  cols = rep(NA, length(pval))
  cols[which(FC <0)] = cols1[1]
  cols[which(FC >0)] = cols1[2]
  w = which(pval>0.05)
  cols[w] = add.alpha(cols[w], alpha = 0.4)
  names(cols)= names
  
  fileout1=concat(c(output_directory,"", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_FC_plot_", batch,".pdf"))
  w=2.5
  pdf(file=fileout1, height=w*1.3, width=w*2.8)
  par(mfrow= c(1,1), mar = c(5,5,2,2))
  
  plot(x,y, xlim = c(-1*range, range),ylim = c(0,5), pch =21, bg = cols, col = add.alpha("white", alpha = 0), cex = 1.2, xlab = "fold change", ylab = "-log10(p-value)", main = "cell proportions combined")
  segments(-100,-log10(0.05), 100, -log10(0.05),col = "grey", lty = 3, lwd= 2)
  segments(0,-100,0,1000,col = "grey", lty = 1, lwd= 2)
  
  w = intersect(which(pval<0.05), which (FC>0))
  if(length(w)>0){
    order = names[w][order(y[w])]
    cex =0.9
    range1 = range(y)
    seq =(seq(from = -log10(0.035),to = 4.9, length = length(order)))
    text( 0.65*(range/1.9), seq, labels = gsub("DC CD207","CD207",gsub("memory","mem.",gsub("_"," ",order))), cex= cex-0.2,las = 1, font = 1, col = cols[order], pos = 4)
    # mtext(side = 4, text = gsub("_"," ",order), line = 0.2,cex= cex-0.2, at = seq, las = 1, font = 1, col = cols[order])
    for(i in c(1:length(order))){
      segments(x[order[i]], y[order[i]], 0.7*(range/1.9), seq[i], lwd = 1.5, lty = 1, col = add.alpha ("grey", alpha = 0.5))
    }
  }
  w = intersect(which(pval<0.05), which (FC<0))
  if(length(w)>0){
    order = names[w][order(y[w])]
    cex =0.9
    range1 = range(y)
    seq =(seq(from = -log10(0.035),to = 4.9, length = length(order)))
    text( -0.75*(range/1.9), seq, labels = gsub("DC CD207","CD207",gsub("memory","mem.",gsub("_"," ",order))), cex= cex-0.2,las = 1, font = 1, col = cols[order], pos = 2)
    # mtext(side = 4, text = gsub("_"," ",order), line = 0.2,cex= cex-0.2, at = seq, las = 1, font = 1, col = cols[order])
    for(i in c(1:length(order))){
      segments(x[order[i]], y[order[i]], -0.8*(range/1.9), seq[i], lwd = 1.5, lty = 1, col = add.alpha ("grey", alpha = 0.5))
    }
  }
  dev.off()
}
Plot_FC_p_value(summary_tables,output_directory)  

############################################ plot boxplots of cell %s in blood between groups
Plot_blood_cell_composition_between_groups<-function(list_normalised, output_directory){
  fileout1=concat(c(output_directory, analysis_all,"_blood_PCA_groups_boxplots_cell_proportions_", batch,".pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*6, width=w*3*2.9)
  par(mfrow= c(6,2), mar = c(16,5,2,3))
  
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  summary_tables = NULL
  for(c in c(1:length(list_normalised))){
    mat = list_normalised[[c]]
    name = names(list_normalised)[c]
    mat_stat = mat
    for(i in c(1:length(mat[1,]))){
      w = which(mat[,i]==-1)
      mat_stat[w,i] =NA
    }
    analysis = concat(c(analysis_all,"_", name))
    
    mat_biopsy = mat_stat[unlist(groups_PCA),]
    mat_blood = mat_stat[gsub("biopsy","blood",unlist(groups_PCA)),]
    
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    group_PCA_list = gsub("biopsy","blood", group_PCA_list)
    mat_blood = mat_blood[,which(colSums(mat_blood[which(is.na(mat_blood[,1])==F),])>0)]
    ##### biopsy
    mat1 = mat_blood[group_PCA_list,]
    
    factor = factor(factor_PCA)
    fit = manova(formula = mat1^0.5 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
    Means_factor = function(factor, x){
      m = NULL
      for(i1 in c(1:length(levels(factor)))){
        x1 = x[which(factor==levels(factor)[i1])]
        x1 = x1[which(x1!=-1)]
        m = c(m, mean(x1))}
      return(m)}
    i1 = 0
    for(i in p1){
      i1 = i1+1
      p_value = c(p_value, i$'Pr(>F)'[1]) 
      if(length(mean)==0){means = Means_factor(factor, mat1[,i1])
      }else{means = rbind(means, Means_factor(factor, mat1[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    colnames(means) = paste("mean.group.", c(1:length(means[1,])))
    combined_p_value = cbind(p_value ,means)
    rownames(combined_p_value) = nam
    p.site = rep("biopsy", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    summary_tables = rbind(summary_tables, x)
    ##### blood
    mat1 = mat_blood
    factor = factor(factor_PCA)
    fit = manova(formula = mat1 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
    Means_factor = function(factor, x){
      m = NULL
      for(i1 in c(1:length(levels(factor)))){
        x1 = x[which(factor==levels(factor)[i1])]
        x1 = x1[which(x1!=-1)]
        m = c(m, mean(x1))}
      return(m)}
    i1 = 0
    for(i in p1){
      i1 = i1+1
      p_value = c(p_value, i$'Pr(>F)'[1]) 
      if(length(mean)==0){means = Means_factor(factor, mat1[,i1])
      }else{means = rbind(means, Means_factor(factor, mat1[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    colnames(means) = paste("mean.group.", c(1:length(means[1,])))
    combined_p_value_blood = cbind(p_value ,means)
    rownames(combined_p_value_blood) = nam
    p.site = rep("blood", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    summary_tables = rbind(summary_tables, x)
    
    groups = NULL
    p_values = combined_p_value_blood[,1]
    mat = mat_blood
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA))){
        x = mat [gsub("biopsy","blood", groups_PCA[[g]]) , i]
        x = x[which(is.na(x)==F)]
        data = c(data, x)
        factor = c(factor, rep(g, length(x)))
        g1 = c(g1, list(x))}
      groups = c(groups, list(g1))
      # group = factor(factor)
      # fit = aov(formula = data ~ group)
      # pvalues = summary.lm(fit)$ coefficients[,4]
      # pVal <- anova(fit)$'Pr(>F)'[1]
      # p_values = c(p_values, pVal)
    }
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = gsub("_"," ",colnames(mat))
    main = concat(c("% of ",name))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 70
    max_scale = min(c(max, 100))
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(2.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.18
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.25/max(shift)
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        if(p_values[i] <signif_threshold/10){pval1 = "**"}
        # if(p_values[i] <signif_threshold/100){pval1 = "***"}
        y = max(unlist(groups[[i]]))
        y = y+1*b
        # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
        text(i, y+2*b, labels = pval1, cex = 1.7)
      }
    }
  }
  
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
  legend("topright", factors1, pch = 22,cex= 0.9, bty="n", pt.bg = add.alpha("white",alpha =0), col = add.alpha("white",alpha =0), text.col = cols1, text.font = 2)
  dev.off()
  return(summary_tables)
}
summary_tables = Plot_blood_cell_composition_between_groups(list_normalised,output_directory)

outfile = concat(c(output_directory,"Stats_PCA_cell_types_by_groups_PDAC150Ka.txt"))
write.table(summary_tables, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")

############################################ plot boxplots of cell %s in tumour and blood
Plot_cell_composition_between_blood_and_tumour<-function(list_normalised,output_directory){
  fileout1=concat(c(output_directory,analysis_all,"_blood_versus_biopsy_PCA_groups_boxplots_cell_proportions_", batch,".pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*6, width=w*3*3.9)
  par(mfrow= c(6,2), mar = c(16,5,2,3))
  
   summary_tables = NULL
  analysis1 = "Blood_versus_biopsy_cell_proportions"
  
  library(RColorBrewer)
  cols = add.alpha (brewer.pal(8, "Set1")[c(2:1)], alpha = 0.95)
  cols1 = add.alpha (c(cols,alpha = 0.5))
  cols2 =  brewer.pal(8, "Set1")[c(2:1)]
  
  p_values_summary = NULL
  analysis_summary = NULL
  analysis_group = NULL
  analysis_population = NULL
  mean_biopsy = NULL
  mean_blood = NULL
  list_normalised1 = list_normalised
  #list_normalised1[[1]] = list_normalised_overall
  for(c in c(1:length(list_normalised1))){
    if(c==3){	plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ="", axes = FALSE, ylim = c(min, max))}
    mat = list_normalised1[[c]]
    name = names(list_normalised1)[c]
    analysis = concat(c(analysis1,"_", name))	
    mat_stat = mat
    for(i in c(1:length(mat[1,]))){
      w = which(mat[,i]==-1)
      mat_stat[w,i] =NA
    }
    mat_biopsy = mat_stat[unlist(groups_PCA),]
    mat_blood = mat_stat[gsub("biopsy","blood",unlist(groups_PCA)),]
    
    w = intersect(which(is.na(mat_biopsy[,1])==F), which(is.na(mat_blood[,1])==F))
    mat_biopsy = mat_biopsy[w,]
    mat_blood = mat_blood[w,]
    
    groups = NULL
    p_values = NULL
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      x1 = mat_biopsy[,i]
      x2 = mat_blood[,i]
      x1 = x1[which(is.na(x1)==F)]
      x2 = x2[which(is.na(x2)==F)]
      g1 = c(list(x1),list(x2))
      pv = wilcox.test(x1, y = x2,alternative = c("two.sided"),paired = T)$p.value
      p_values  = c(p_values , pv)
      groups = c(groups,list(g1))
      p_values_summary = c(p_values_summary, pv)
      analysis_summary = c(analysis_summary, analysis1)
      analysis_population = c(analysis_population, colnames(mat)[i])
      analysis_group = c(analysis_group, name)
      mean_biopsy = c(mean_biopsy, mean(x1))
      mean_blood = c(mean_blood, mean(x2))
    }
   
    factors1 = paste("group",c(1:length(g1)))
    factors = gsub("_"," ",colnames(mat))
    main = concat(c("% of ",name))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 90
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(3.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.15
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.2/max(shift)
    
    library(RColorBrewer)
    cols = add.alpha (brewer.pal(8, "Set1")[c(2:1)], alpha = 0.95)
    cols1 = add.alpha (cols,alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols1[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21,bg = cols1[i1], col=cols[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        if(p_values[i] <signif_threshold/10){pval1 = "**"}
        # if(p_values[i] <signif_threshold/100){pval1 = "***"}
        y = max(unlist(groups[[i]]))
        y = y+1*b
        # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
        text(i, y+2*b, labels = pval1, cex = 1.3)
      }
    }
  }
  
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
  legend("topright", c("PDAC tumour", "blood"), pch = 22,cex= 0.9, bty="n", pt.bg = add.alpha("white",alpha =0), col = add.alpha("white",alpha =0), text.col = cols, text.font = 2)
  dev.off()
  
  return(summary_tables)
}
Plot_cell_composition_between_blood_and_tumour(list_normalised,output_directory)


############################################ plot boxplots of tumour CD4-CD8 %
Plot_cell_CD4_CD8_ratio_tumour<-function(list_normalised,output_directory){
  analysis_all = "T_cell_CD4_CD8_proportions"
  fileout1=concat(c(output_directory, analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_", batch,".pdf"))
  w=4
  pdf(file=fileout1, height=w*0.8*1, width=w*1*0.8)
  par(mfrow= c(1,2), mar = c(6,4,2,1))
  
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  list_CD8_prop= NULL
  for(c in c(1:length(cd4cd8))){
    mat = cd4cd8[[c]]
    groups = NULL
    for(pca in c(1:length(groups_PCA))){
      groups = c(groups, list(mat[gsub("-","_",groups_PCA[[pca]], fixed = T)]))
    }
    names(groups) = c(1:2)
    # list_CD8_prop = c(list_CD8_prop, list(unlist(groups)))
    p_val = wilcox.test(groups[[1]], y = groups[[2]])$p.value
    
    factors =c("AL group", "AE group")
    main = concat(c(names(cd4cd8)[c]))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 2
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max),cex.main = 0.7)
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.25
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.25/max(shift)
    
    for(i in c(1:l)){
      points1=as.numeric(groups[[i]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i,width,cols[i],1, cols1[i])
      points(rep(i, length(points1)),points1, pch =21, col=cols[i],bg = cols[i], cex = 0.7)
    }
    b = max*0.045
    signif_threshold = 0.05
    # if(p_val <signif_threshold/10){pval1 = "**"}
    # if(p_val <signif_threshold/100){pval1 = "***"}
    y = max(unlist(groups))
    y = y+1*b
    segments(1,y+b, 2,y+b,lwd = 3, col = "darkgrey")
    pval1 = "NS"
    if(p_val <signif_threshold){pval1 = "*"
    text(1.5, y+2.5*b, labels = pval1, cex = 1.5)
    }else{
      text(1.5, y+2.5*b, labels = pval1, cex = 0.7)}
  }
  
  dev.off()
}
Plot_cell_CD4_CD8_ratio_tumour(list_normalised,output_directory)

############################################ plot boxplots of tumour activated Treg of all Tregs %
Plot_cell_act_Treg_ratio_tumour<-function(T_ratios,output_directory){
  analysis_all = "Treg_act_of_Tregs_proportions"
  fileout1=concat(c(output_directory, analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_", batch,".pdf"))
  w=4
  pdf(file=fileout1, height=w*0.65*1, width=w*1*0.8)
  par(mfrow= c(1,3), mar = c(6,4,2,1))
  
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  for(c in c(1:length(T_ratios[1,]))){
    mat = T_ratios[,c]
    name = colnames(T_ratios)[c]
    groups = NULL
    for(pca in c(1:length(groups_PCA))){
      groups = c(groups, list(mat[gsub("-","_",groups_PCA[[pca]], fixed = T)]))
    }
    names(groups) = c("ME","AE")
    # list_CD8_prop = c(list_CD8_prop, list(unlist(groups)))
    p_val = wilcox.test(groups[[1]], y = groups[[2]])$p.value
    
    factors =c("ME group", "AE group")
    main = name
    max = max(c(unlist(groups), unlist(groups))*1.25)
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 2
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max),cex.main = 0.7)
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.45
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.25/max(shift)
    
    for(i in c(1:l)){
      points1=as.numeric(groups[[i]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i,width,cols[i],1, cols1[i])
      points(rep(i, length(points1)),points1, pch =21, col=cols[i],bg = cols[i], cex = 0.7)
    }
    b = max*0.045
    signif_threshold = 0.05
    # if(p_val <signif_threshold/10){pval1 = "**"}
    # if(p_val <signif_threshold/100){pval1 = "***"}
    y = max(unlist(groups))
    y = y+1*b
    segments(1,y+b, 2,y+b,lwd = 3, col = "darkgrey")
    pval1 = "NS"
    if(p_val <signif_threshold){pval1 = "*"
    text(1.5, y+2.5*b, labels = pval1, cex = 1.5)
    }else{
      text(1.5, y+2.5*b, labels = pval1, cex = 0.7)}
  }
  
  dev.off()
  
}
Plot_cell_act_Treg_ratio_tumour(T_ratios,output_directory)

############################################ plot correlation between blood and tumour proportions
Plot_correlation_tumour_blood_CD4_CD8<-function(list_normalised, output_directory){
  blood = names(cd4cd8[[1]])[intersect(grep("blood", names(cd4cd8[[1]])),grep("PDAC", names(cd4cd8[[1]]))) ]
  biopsy = names(cd4cd8[[1]])[grep("biopsy", names(cd4cd8[[1]]))]
  cbind(blood, biopsy)
  
  c = 2
  plot(cd4cd8[[c]][blood],cd4cd8[[c]][biopsy], main = names(cd4cd8)[c], xlab = "blood",ylab = "tumour")
  
  analysis_all = "T_cell_CD4_CD8_proportions"
  fileout1=concat(c(output_directory, analysis_all,"_correlations_blood_biopsy_", batch,".pdf"))
  w=3.
  pdf(file=fileout1, height=w*4, width=w*1*2)
  par(mfrow= c(4,2), mar = c(4,4,3,3))
  cex = 0.9
  for(c in c(1:length(cd4cd8))){	
    x = cd4cd8[[c]][blood]
    y = cd4cd8[[c]][biopsy]
    factors = c(rep(1,length(groups_PCA[[1]])), rep(2,length(groups_PCA[[2]])))
    
    max = max(c(x,y))*1.1
    min = min(c(x,y))
    min = 0
    range = max-min
    scale = c(0:100)*10
    if(range>50){scale = c(0:100)*20}
    if(range<=20){scale = c(0:100)*5}
    if(range <15){scale = c(0:100)*2.5}
    if(range <10){scale = c(0:100)*2}
    if(range <5){scale = c(0:100)*1}
    if(range <2){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    scale = scale[intersect(which(scale<= max), which(scale>=min))]
    range = range(c(x,y,min, scale))
    # range1= range(c(max(x)*1.4,x, min))
    # range2= range(c(max(y)*1.2,y, min))
    x11 = x
    y11 = y
    fit1 <- lm( y11~poly(x11,2))
    fit1 <- lm( y11~x11)
    p_value = summary.lm(fit1)$ coefficients["x11",4]
    
    plot(range, range, pch=20, col="white",main = concat(c(names(cd4cd8)[c],"\np-value:",signif(p_value,digits = 3))), xlab = "blood",ylab = "tumour",cex=cex, cex.lab=cex+0.1,cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), axes = FALSE)
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    pches = c(21:25)
    
    Fun<-function(x){x}
    segments(min,Fun(scale),max(max)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    segments(Fun(scale),min,Fun(scale),max(max)+0.5,col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.35,cex= cex,  at =Fun(scale), las = 2, font = 1)
    scale1 = scale#[which(scale!=0)]
    # scale1 = scale1[which(scale1<=max(x)*1.2)]
    mtext(side = 1, text = scale1, line = 0.35,cex= cex,  at =Fun(scale1), las = 2, font = 1)
    points(x,y, pch = pches[factors], bg=cols[factors], cex = 1, col= add.alpha("white", alpha = 0))	
    
    new = seq(from= min(x11), to = max(x11), length = 100)
    predicted.intervals <- predict(fit1, data.frame(x11= new),interval='prediction',level = 0.85)
    cols2 =  add.alpha ("darkblue", alpha = 0.05)
    points(new,predicted.intervals[,1],col= cols2,lwd=4, type = "l")
    polygon(x = c(new, rev(new)),y=c(predicted.intervals[,2],rev(predicted.intervals[,3]) ),  col= add.alpha(cols2,alpha = 0.15),border = NA)
  }
  
  mat = list_normalised [["CD45:broad" ]]
  for(c in c(1:length(mat[1,]))){	
    x = mat[blood,c]
    y = mat[biopsy,c]
    factors = c(rep(1,length(groups_PCA[[1]])), rep(2,length(groups_PCA[[2]])))
    feature = colnames(mat)[c]
    
    max = max(c(x,y))*1.1
    min = min(c(x,y))
    min = 0
    range = max-min
    scale = c(0:100)*10
    if(range>50){scale = c(0:100)*20}
    if(range<=20){scale = c(0:100)*5}
    if(range <15){scale = c(0:100)*2.5}
    if(range <10){scale = c(0:100)*2}
    if(range <5){scale = c(0:100)*1}
    if(range <2){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    scale = scale[intersect(which(scale<= max), which(scale>=min))]
    range = range(c(x,y,min, scale))
    # range1= range(c(max(x)*1.4,x, min))
    # range2= range(c(max(y)*1.2,y, min))
    x11 = x
    y11 = y
    fit1 <- lm( y11~poly(x11,2))
    fit1 <- lm( y11~x11)
    p_value = summary.lm(fit1)$ coefficients["x11",4]
    
    plot(range, range, pch=20, col="white",main = concat(c("% ", feature,"of CD45+\np-value:",signif(p_value,digits = 3))), xlab = "blood",ylab = "tumour",cex=cex, cex.lab=cex+0.1,cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), axes = FALSE)
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    pches = c(21:25)
    
    Fun<-function(x){x}
    segments(min,Fun(scale),max(max)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    segments(Fun(scale),min,Fun(scale),max(max)+0.5,col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.35,cex= cex,  at =Fun(scale), las = 2, font = 1)
    scale1 = scale#[which(scale!=0)]
    # scale1 = scale1[which(scale1<=max(x)*1.2)]
    mtext(side = 1, text = scale1, line = 0.35,cex= cex,  at =Fun(scale1), las = 2, font = 1)
    points(x,y, pch = pches[factors], bg=cols[factors], cex = 1, col= add.alpha("white", alpha = 0))	
    
    new = seq(from= min(x11), to = max(x11), length = 100)
    predicted.intervals <- predict(fit1, data.frame(x11= new),interval='prediction',level = 0.85)
    cols2 =  add.alpha ("darkblue", alpha = 0.05)
    points(new,predicted.intervals[,1],col= cols2,lwd=4, type = "l")
    polygon(x = c(new, rev(new)),y=c(predicted.intervals[,2],rev(predicted.intervals[,3]) ),  col= add.alpha(cols2,alpha = 0.15),border = NA)
    
  }
  
  dev.off()
  
}
Plot_correlation_tumour_blood_CD4_CD8(list_normalised,output_directory)

############################################ plot correlation between blood and tumour proportions
Plot_correlation_tumour_blood<-function(list_normalised,output_directory){
  biopsy = rownames(list_normalised[[1]])[grep("biopsy", rownames(list_normalised[[1]]))]
  blood = gsub("biopsy","blood", biopsy)
  library(RColorBrewer)
  cols1 =  c(add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95),add.alpha (brewer.pal(8, "Paired"), alpha = 0.95))
  cols2 =  add.alpha (cols1, alpha = 0.75)
  cols =  add.alpha (cols1, alpha = 0.25)
  cols_grey =  add.alpha ("grey", alpha = 0.5)
  
  analysis_all = "Cell_proportions_blood_vv_biopsy"
  fileout1=concat(c(output_directory, analysis_all,"_correlations_blood_biopsy_", batch,".pdf"))
  w=3.
  pdf(file=fileout1, height=w*4, width=w*1*4)
  par(mfrow= c(4,4), mar = c(4,4,3,3))

  library(psych)
  for(c in c(1:length(list_normalised))){	
    x = list_normalised[[c]][blood,]
    y = list_normalised[[c]][biopsy,]
    
    cell_type_sub = colnames(x)
    corr_mat = matrix(data = NA,nrow = length(cell_type_sub), ncol = length(cell_type_sub), dimnames = c(list(cell_type_sub), list(cell_type_sub)))
    for(i in c(1:length(cell_type_sub))){
      for(j in c(1:length(cell_type_sub))){
        if(i<j){
          m = cbind(x[,cell_type_sub[i]], y[,cell_type_sub[j]])
          corr_mat[i,j] = corr.test(m,adjust="holm")$p[1,2]
        }
      }
    }
    corr_mat1 = corr_mat
    table(corr_mat1[which(is.na(corr_mat1)==F)])
    
    from = NULL
    to = NULL
    for(i in c(1:length(cell_type_sub))){
      for(j in c(1:length(cell_type_sub))){
        if(is.na(corr_mat1[i,j])==F){
          corr_mat1[i,j] = p.adjust(corr_mat[i,j], method = "holm", n = length(corr_mat[which(is.na(corr_mat)==F)]))
          if(corr_mat1[i,j]<0.05){
            from = c(from, cell_type_sub[i])
            to = c(to, cell_type_sub[j])
          }}}}
    if(length(from)!=0){
      for(ind in c(1:length(from))){
        main = concat(c(names(list_normalised)[c], "\n",from[ind],"-",to[ind] ))
        x1 = x[,from[ind]]
        y1 = y[,to[ind]]
        models = NULL
        
        cortest = corr.test(cbind(x1,y1))
        pval = corr_mat1[from[ind],to[ind]]
        rval = cortest $ r[1,2]
        fit =  lm( y1~x1)
        
        main = concat(c(main,"\np:", signif(pval, digits = 3)))
        plot(x1,y1,pch = 21, bg = cols[c],col = cols1[c], main = main, xlab =from[ind], ylab = to[ind],cex.main = 0.7)
        
        new = seq(from= min(x1)*1.5, to = max(x1)*1.5, length = 100)
        predicted.intervals <- predict(fit, data.frame(x1= new),interval='prediction',level = 0.85)
        points(new,predicted.intervals[,1],col= cols2[c],lwd=4, type = "l")
        polygon(x = c(new, rev(new)),y=c(predicted.intervals[,2],rev(predicted.intervals[,3]) ), col= cols[c],border = NA)
      }
    }
  }
  
  dev.off()
  
}
Plot_correlation_tumour_blood(list_normalised,output_directory)

############################################ plot healthy versus PDAC blood
Get_list_bloods<-function(p_blood, output_directory){
  file = concat(c(output_directory,"Seurat_cell_counts_",batch,"_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  cell_types = p[,1]
  cell_group = p[,2]
  cell_groups = sort(unique(cell_group))
  cell_groups = cell_groups[which(cell_groups!="Overall")]
  
  p_blood = unique(colnames(p)[c(grep("blood",colnames(p)), grep("HC",colnames(p)), grep("pbmc",colnames(p)))])
  p_blood = p[,p_blood]
  
  p = p_blood
  p = cbind(cell_types, cell_group, p)
  rownames(p) = cell_types
  sample = colnames(p)[-1]
  sample = sample[-1]
  sites = sample
  disease = sample
  sites[grep("PDAC", sites)] = "PDAC PBMC"
  sites[grep("PDAC", sites, invert = T)] = "HC PBMC"
  names(sites) = sample
  
  list_normalised = NULL
  
  include = NULL
  for(ind in c(1:length(cell_groups))){
    w = which(cell_group== cell_groups[ind])
    if(cell_groups[ind] =="CD45:broad"){
      w = intersect(w, which(p[,1] %in% c( "B cell"  ,"ILC", "Myeloid","NK" , "T cell"   )))
    }
    cell_type_sub = p[w,1]
    p1 = p[w,]
    m_count_norm = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(sample), dimnames = c(list(cell_type_sub), list(sample)))
    
    for(s in c(1:length(sample))){
      m_count_norm[, sample[s]]= as.numeric(p1[, sample[s]])*100/sum(as.numeric(p1[, sample[s]]))
      if(sum(as.numeric(p1[, sample[s]]))<7){m_count_norm[, sample[s]] = -1}
    }
    if(length(which(apply(m_count_norm, 2, max)!=-1))>3){include = c(include, ind)}
    m_count_norm= m_count_norm[order(rownames(m_count_norm)),]
    list_normalised= c(list_normalised, list(t(m_count_norm)))
  }
  names(list_normalised) = cell_groups
  list_normalised = list_normalised[include]
  
  
  for(i in c(1:length(list_normalised))){
    print (names(list_normalised)[i])
    print (colnames(list_normalised[[i]]))
  }
  
  x = list_normalised[["CD45:broad"]]
  B_T_cell = rowSums(x[,c("B cell", "T cell")])
  x = cbind(x[,c("Myeloid", "NK"  )], B_T_cell)
  colnames(x)[3] = "B/T cells"
  
  y = list_normalised[["CD45:B cell" ]]
  broad = c("plasma", " memory","activated pre-memory", "GC", "MZ", "naive")
  y1 =NULL
  for(i in c(1:length(broad))){
    w = colnames(y)[grep(broad[i], colnames(y))]
    if(length(w)>1){yi = rowSums(y[,w])
    }else{yi =(y[,w])}
    if(length(y1)==0){y1 = yi
    }else{y1=cbind(y1, yi)}
  }
  colnames(y1) = gsub(" ","", gsub(" pre-memory","", gsub("plasma","ASC", broad, fixed = T), fixed = T),fixed = T)
  y = y1
  
  z = list_normalised[["CD45:T cell" ]]
  broad = c("CD4", "CD8","MAIT", "gd")
  z1 =NULL
  for(i in c(1:length(broad))){
    w = colnames(z)[grep(broad[i], colnames(z))]
    if(length(w)>1){yi = rowSums(z[,w])
    }else{yi =(z[,w])}
    if(length(z1)==0){z1 = yi
    }else{z1=cbind(z1, yi)}
  }
  colnames(z1) = gsub("gd","gamma/delta", broad, fixed = T)
  z = z1
  
  s = list_normalised[["CD45:Myeloid" ]]
  broad = c("DC", "granulo","Mast", "megakaryocyte","momac","monocyte")
  s1 =NULL
  for(i in c(1:length(broad))){
    w = colnames(s)[grep(broad[i], colnames(s))]
    print(length(w))
    if(length(w)>1){yi = rowSums(s[,w])
    }else{yi =(s[,w])}
    if(length(s1)==0){s1 = yi
    }else{s1=cbind(s1, yi)}
  }
  colnames(s1) = gsub("granulo","granulocyte", broad, fixed = T)
  s = s1
  
  l = length(list_normalised)
  list_normalised1 = c(list_normalised, list(x), list(y), list(z), list(s))
  names(list_normalised1)[c((l+1):length(list_normalised1))] = c("CD45:broad condensed", "CD45:broad B cell","CD45:broad T cell", "CD45:broad Myeloid")
  
  list_normalised = list_normalised1
  
  for(i in c(1:length(list_normalised))){
    print (names(list_normalised)[i])
    print (colnames(list_normalised[[i]]))
  }
  
  T_cells = list_normalised[["CD45:T cell"]]
  Treg_CD8_ratio = rowSums(T_cells [,grep("Treg",colnames(T_cells))])/rowSums(T_cells [,grep("CD8",colnames(T_cells))])
  
  T_cells = list_normalised[["CD45"]]
  CD4_CD8_ratio = rowSums(T_cells [,grep("CD4",colnames(T_cells))])/rowSums(T_cells [,grep("CD8",colnames(T_cells))])
  
  cd4cd8 = c(list(Treg_CD8_ratio), list(CD4_CD8_ratio))
  names(cd4cd8) = c("Treg-CD8 T cell ratio", "CD4-CD8 ratio")
  return(list_normalised)
}
  
Plot_PCA<-function(list_normalised,sites, output_directory){
  all_proportions = (list_normalised[["CD45" ]])
  type_to_plot = "PDAC versus Health blood"
  
  x <- princomp(all_proportions,dim=2)
  PCA_coordinates_sub = x$scores
  
  x1 = x$scores[,1]
  x2 = x$scores[,2]
  
  y1 = x$loadings[,1]
  y2 = x$loadings[,2]
  
  factor1 = rep(0,length(x1))
  names(factor1) = names(x1)
  factor1[which(sites=="HC PBMC")] = 3
  factor1[which(names(x1) %in% gsub("_biopsy","_blood",group1))] = 1
  factor1[which(names(x1) %in% gsub("_biopsy","_blood",group2))] = 2
  
  factor2 = names(x1)
  
  factors1 = sort(unique(factor1))
  factors2 = sort(unique(factor2))
  matching1 = match(factor1, factors1)
  matching2 = match(factor2, factors2)
  
  library(RColorBrewer)
  cols =  add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95)
  cols1 = add.alpha (c( cols[1],"white"),alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8, 3)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)
  
  pches = c(21:25)
  fileout1=concat(c(output_directory, analysis_all,"_blood_PCA_scores_", batch,".pdf"))
  w=2.2
  pdf(file=fileout1, height=w*1, width=w*3)
  par(mfrow= c(1,3), mar = c(5,5,3,3))
  
  plot(x1,x2,pch = pches[matching1], col = cols[matching1],bg = cols2[matching1], main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  
  xy = cbind(x1,x2)
  library(car)
  for(pca in c(1:length(factors1))){
    w = which(factor1 ==pca)
    if(length(w)>1){
      x = xy[w,1]
      y = xy[w,2]
      dataEllipse(x, y, levels=c(0.65), add = TRUE, col = cols2[pca], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
    }
  }
  
  plot(x1,x2,pch = pches[matching1], col = "white",bg = "white", main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  legend("topleft", paste("group",factors1), pch = pches,cex= 0.8, bty="n", pt.bg = cols2, col = cols, pt.lwd = 2, text.font = 2)
  
  plot(x1,x2,pch = pches[matching1], col = cols[matching1],bg = cols2[matching1], main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  text(x1, y = x2, labels = factor2,cex = 0.4,font = 1)
  
  dev.off()
}
Plot_PCA(list_normalised)
############################################ plot boxplots of cell %s in tumour between groups
Plot_blood_cell_composition_between_PDAC_health<-function(list_normalised, output_directory){
  group1 = names(sites)[which(sites=="HC PBMC")]
  group2 = names(sites)[which(sites!="HC PBMC")]
  groups_PCA = c(list(group1), list(group2))
  
  fileout1=concat(c(output_directory, analysis_all,"_blood_PDAC_health_boxplots_cell_proportions_", batch,".pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*6, width=w*3*2.9)
  par(mfrow= c(6,2), mar = c(16,5,2,3))
  
  summary_tables = NULL
  analysis1 = "cell_proportions"
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(3,6)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  for(c in c(1:length(list_normalised))){
    mat = list_normalised[[c]]
    name = names(list_normalised)[c]
    analysis = concat(c(analysis1,"_", name))	
    mat_stat = mat
    for(i in c(1:length(mat[1,]))){
      w = which(mat[,i]==-1)
      mat_stat[w,i] =NA
    }
    
    mat_biopsy = mat_stat[unlist(groups_PCA),]

    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    #mat_biopsy = mat_biopsy[,which(colSums(mat_biopsy)>0)]
    ##### biopsy
    mat1 = mat_biopsy[group_PCA_list,]
    factor = factor(factor_PCA)
    fit = manova(formula = mat1^0.5 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
    Means_factor = function(factor, x){
      m = NULL
      sd = NULL
      for(i1 in c(1:length(levels(factor)))){
        x1 = x[which(factor==levels(factor)[i1])]
        x1 = x1[which(x1!=-1)]
        sd = c(sd, sd(x1))
        m = c(m, mean(x1))}
      return(c(m, sd))}
    i1 = 0
    for(i in p1){
      i1 = i1+1
      p_value = c(p_value, i$'Pr(>F)'[1]) 
      if(length(mean)==0){means = Means_factor(factor, mat1[,i1])
      }else{means = rbind(means, Means_factor(factor, mat1[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    colnames(means) = c(paste("mean.group.", c(1:length(levels(factor)))), paste("sd.group.", c(1:length(levels(factor)))))
    combined_p_value = cbind(p_value ,means)
    rownames(combined_p_value) = nam
    p.site = rep("biopsy", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    if(length(summary_tables)==0){summary_tables = x
    }else{summary_tables = rbind(summary_tables ,x)}
    
    groups = NULL
    mat = mat_biopsy#[,which(colSums(mat_biopsy)>0)]
    p_values = combined_p_value[,1]
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA))){
        x = mat [groups_PCA[[g]], i]
        x = x[which(is.na(x)==F)]
        data = c(data, x)
        factor = c(factor, rep(g, length(x)))
        g1 = c(g1, list(x))}
      groups = c(groups, list(g1))
      # group = factor(factor)
      # fit = aov(formula = data ~ group)
      # pvalues = summary.lm(fit)$ coefficients[,4]
      # pVal <- anova(fit)$'Pr(>F)'[1]
      # p_values = c(p_values, pVal)
    }
    factors1 = c("HC blood", "PDAC blood")
    factors = gsub("_"," ",colnames(mat))
    main = concat(c("% of ",name))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 70
    max_scale = min(c(max, 100))
    range = max-min
    if(range>55){scale = c(0:100)*20}
    if(range<=55){scale = c(0:100)*10}
    if(range<=40){scale = c(0:100)*5}
    if(range <10){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(2.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.18
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.25/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg = cols[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        # if(p_values[i] <signif_threshold/10){pval1 = "**"}
        # if(p_values[i] <signif_threshold/100){pval1 = "***"}
        y = max(unlist(groups[[i]]))
        y = y+1*b
        # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
        text(i, y+2*b, labels = pval1, cex = 1.7)
      }
    }
  }
  
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
  legend("topright", factors1, pch = 22,cex= 0.9, bty="n", pt.bg = add.alpha("white",alpha =0), col = add.alpha("white",alpha =0), text.col = cols1, text.font = 2)
  dev.off()
  
  return(summary_tables)
}
summary_tables = Plot_tumour_cell_composition_between_groups(list_normalised, output_directory)

outfile = concat(c(output_directory,"Stats_PCA_cell_types_by_groups_PDAC150Ka.txt"))
write.table(summary_tables, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")

  
  
Means_factor = function(factor, x){
	m = NULL
	for(i1 in c(1:length(levels(factor)))){
		x1 = x[which(factor==levels(factor)[i1])]
		x1 = x1[which(x1!=-1)]
		m = c(m, mean(x1))}
	return(m)}

Medians_factor = function(factor, x){
  m = NULL
  for(i1 in c(1:length(levels(factor)))){
    x1 = x[which(factor==levels(factor)[i1])]
    x1 = x1[which(x1!=-1)]
    m = c(m, median(x1))}
  return(m)}

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
batch = "PDAC150K"
analysis_all = "PSEUDO_BULK_per_cell_type"
PLOTS = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Psuedo bulk/Psuedobulk per cell type/"
input_directory_groups  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"
output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/RECEPTOR_LIGAND_DGE/"
groups_PCA = readRDS(file = concat(c(input_directory_groups, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))
######################

Boxplot_B_T_cell_biopsy_apoptosis_per_subset_clonal<-function(input_directory, output_directory,groups_PCA){
  ###################### get cell counts
  clone_sizes=readRDS(file = concat(c(input_directory,"Clone_sizes_all.rds")))
  clone_sizes_all = c(clone_sizes[[1]],clone_sizes[[2]])
  
  clone_sizes_all_catagorical = clone_sizes_all
  clone_sizes_all_catagorical[which(clone_sizes_all>=2)] = "expanded"
  clone_sizes_all_catagorical[which(clone_sizes_all<2)] = "unexpanded"
  
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  analysis = "apoptosis"
  raw_pathway_scores = readRDS(file = concat(c(input_directory, "apoptosis_UMAP_reclustering_PDAC150Ka_raw_pathway_level.rds")))
  str(raw_pathway_scores)
  
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA[[i]], fixed = T))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  
  ######
  cell_type= as.character(raw_pathway_scores[,"pbmc@meta.data$cell_refined_annotation"])
  pat_sample  = raw_pathway_scores[,"pat_sample"]
  pat_samples = sort(unique(pat_sample))
  
  inter = intersect(names(clone_sizes_all_catagorical), rownames(raw_pathway_scores))
  clone_size_match = rep("UNKNOWN",length(cell_type))
  names(clone_size_match) = rownames(raw_pathway_scores)
  clone_size_match[inter] = clone_sizes_all_catagorical[inter]
  
  cell_type_clone_size = apply(cbind(cell_type, clone_size_match), 1 , paste, collapse = "-")
  cell_types = sort(unique(cell_type_clone_size))
  #### pathways
  
  pws = colnames(raw_pathway_scores)
  pathways2 = pws[c(grep("KEGG",pws), grep("GO_",pws), grep("REACTOME_",pws), grep("HALLMARK_", pws))]
  list_mean_per_cell_type_pathway= NULL
  for(g in c(1:length(pathways2))){
    mean_per_cell_type_pathway = matrix(data = NA, nrow = length(pat_samples), ncol = length(cell_types), dimnames = c(list(pat_samples), list(cell_types)))
    score = raw_pathway_scores[,pathways2[g]]
    for(s in c(1:length(pat_samples))){
      for(c in c(1:length(cell_types))){
        w = intersect(which (pat_sample==pat_samples[s]), which(cell_type_clone_size==cell_types[c]))
        if(length(w)>=5){
          mean_per_cell_type_pathway[pat_samples[s], cell_types[c]]= mean(score[w])
        }
      }
    }
    list_mean_per_cell_type_pathway = c(list_mean_per_cell_type_pathway, list(mean_per_cell_type_pathway))
    print(g)
  }
  names(list_mean_per_cell_type_pathway) = pathways2
  
  #
  
  
  # ki67 expression = MKI67
  MKI67 = raw_pathway_scores[,"MKI67[rownames(pbmc@meta.data)]"]
  w = which(MKI67!=0)
  t0 = table(cell_type)
  t1 = t0*0
  ki67_per_cell_type = matrix(data = 0, nrow = length(pat_samples), ncol = length(cell_types), dimnames = c(list(pat_samples), list(cell_types)))
  for(s in c(1:length(pat_samples))){
    w1 = which (pat_sample==pat_samples[s])
    if(length(w1)>=10){
      w2 = intersect(w,w1)
      t1a = table(cell_type_clone_size[w2])
      t2a = table(cell_type_clone_size[w1])
      w3 = which(t2a<5)
      ki67_per_cell_type[pat_samples[s],names(t1a)] = t1a
      ki67_per_cell_type[pat_samples[s],names(t2a)] = ki67_per_cell_type[pat_samples[s],names(t2a)]*100/t2a
      ki67_per_cell_type[pat_samples[s],names(w3)] = NA
    }}
  
  list_mean_per_cell_type_pathway = c(list_mean_per_cell_type_pathway, list(ki67_per_cell_type))
  names(list_mean_per_cell_type_pathway)[length(names(list_mean_per_cell_type_pathway))] = "% Ki67+"
  
  # cell cycle prediction
  
  cell_cycle = raw_pathway_scores[,"pbmc@meta.data$Phase"]
  cell_cycles = sort(unique(cell_cycle))
  
  t0 = table(cell_cycle)
  t1 = t0*0
  S_per_cell_type = matrix(data = 0, nrow = length(pat_samples), ncol = length(cell_types), dimnames = c(list(pat_samples), list(cell_types)))
  w = which(cell_cycle=="S")
  for(s in c(1:length(pat_samples))){
    w1 = which (pat_sample==pat_samples[s])
    if(length(w1)>=10){
      w2 = intersect(w,w1)
      t1a = table(cell_type_clone_size[w2])
      t2a = table(cell_type_clone_size[w1])
      w3 = which(t2a<5)
      S_per_cell_type[pat_samples[s],names(t1a)] = t1a
      S_per_cell_type[pat_samples[s],names(t2a)] = S_per_cell_type[pat_samples[s],names(t2a)]*100/t2a
      S_per_cell_type[pat_samples[s],names(w3)] = NA
    }}
  
  list_mean_per_cell_type_pathway = c(list_mean_per_cell_type_pathway, list(S_per_cell_type))
  names(list_mean_per_cell_type_pathway)[length(names(list_mean_per_cell_type_pathway))] = "% cells in S phase"
  
  saveRDS(file = concat(c(output_directory, "apoptosis_UMAP_reclustering_PDAC150Ka_per_cell_type_expanded_pathway_level.rds")), list_mean_per_cell_type_pathway)

  ####
  list_counts=readRDS(file = concat(c(output_directory, "apoptosis_UMAP_reclustering_PDAC150Ka_per_cell_type_expanded_pathway_level.rds")))
  
  
  fileout1 = concat(c(output_directory,"B_T_cell_apoptosis_boxplot_counts_per_patient_grouped_per_subest_expanded.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1.6*2, width=w*2.1*3)
  par(mfrow= c(2,3), mar = c(25,5,2,3))
  
  for(subset in c(1:length(list_counts))){
    types_sub = c("B cell", "T cell", "NK cell")
    for (tt in c(1:length(types_sub))){
      name = names(list_counts)[subset]
      proportions = list_counts[[subset]]
      w = grep(types_sub[tt], colnames(proportions))
      mat_stat = proportions[group_PCA_list,w]
      w = which(apply(mat_stat, 2, function(x){length(which(is.na(x)==F))})>=10)
      mat_stat = mat_stat[,w]
      factor = factor(factor_PCA)
      mat_stat1 = mat_stat-min(mat_stat[which(is.na(mat_stat)==F)])
      fit = manova(formula = mat_stat1^1 ~ factor)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      
      i1 = 0
      for(i in p1){
        i1 = i1+1
        p_value = c(p_value, i$'Pr(>F)'[1]) 
        if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      print(min(p_value))
      
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.site = rep(concat(c("biopsy ", types_sub[tt])), length(p_value))
      p.analysis = rep(analysis, length(nam))
      x = cbind(p.site, p.analysis, combined_p_value)
      
      #out_file_table = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Bregs_Tregs/Breg_boxplot_counts_per_patient_grouped_tumour_.txt"))
      #write.table(x, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
      
      groups = NULL
      for(i in c(1:length(mat_stat[1,]))){
        x1 = mat_stat[gsub("_","-", groups_PCA[[1]],fixed=T),i]
        x2 = mat_stat[gsub("_","-", groups_PCA[[2]],fixed=T),i]
        x1 = x1[which(is.na(x1)==F)]
        x2 = x2[which(is.na(x2)==F)]
        g1 = c(list(x1),
               list(x2))
        names(g1) = c("ME","AE")
        groups = c(groups, list(g1))
      }
      names(groups)=colnames(mat_stat)
      
      library(RColorBrewer)
      cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
      cols =  add.alpha (cols1, alpha = 0.5)
      
      factors1 = c("ME","AE")
      factors = gsub("_"," ",colnames(mat_stat))
      main = concat(c(name," biopsy ", types_sub[tt]))
      max = max(c(unlist(groups), unlist(groups))*1.15)
      min = min(c(unlist(groups), unlist(groups))*1.15)
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 40
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(-100:100)*20}
      if(range<=50){scale = c(-100:100)*10}
      if(range<=30){scale = c(-100:100)*5}
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
      mtext(side = 1, text = gsub("."," ",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
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
        pval1= "NS"
        if(p_value[i]<signif_threshold){pval1 = "*"
        y = max(unlist(groups[[i]]))
        y = y+3*b
        text(i, y+2*b, labels = pval1, cex = 1.3)
        }}
    }
  }
  dev.off()	
  
}
Boxplot_B_T_cell_biopsy_apoptosis_per_subset_clonal(input_directory, output_directory,groups_PCA)

## not used
Boxplot_B_T_cell_biopsy_apoptosis_per_subset_immunosurveilling<-function(input_directory, output_directory,groups_PCA){
  
  #### get immunosurveillance status
  info = readRDS(file = concat(c(input_directory, "Cell_info_immunosurveilling.rds")))
  clone_sizes_all_catagorical = info[,"type"]
  
  
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  analysis = "apoptosis"
  raw_pathway_scores = readRDS(file = concat(c(input_directory, "apoptosis_UMAP_reclustering_PDAC150Ka_raw_pathway_level.rds")))
  str(raw_pathway_scores)
  
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA[[i]], fixed = T))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  
  ######
  cell_type= as.character(raw_pathway_scores[,"pbmc@meta.data$cell_refined_annotation"])
  pat_sample  = raw_pathway_scores[,"pat_sample"]
  pat_samples = sort(unique(pat_sample))
  
  inter = intersect(names(clone_sizes_all_catagorical), rownames(raw_pathway_scores))
  clone_size_match = rep("UNKNOWN",length(cell_type))
  names(clone_size_match) = rownames(raw_pathway_scores)
  clone_size_match[inter] = clone_sizes_all_catagorical[inter]
  
  cell_type_clone_size = apply(cbind(cell_type, clone_size_match), 1 , paste, collapse = "-")
  cell_type_clone_size = gsub(" blood","", cell_type_clone_size)
  cell_type_clone_size = gsub(" TIL","", cell_type_clone_size)
  cell_types = sort(unique(cell_type_clone_size))
  cell_types = cell_types[grep("-UNKNOWN", cell_types, invert = T)]
  
  t = table(cell_type_clone_size[which(cell_type_clone_size %in% cell_types)])
  exclude = names(which(t<25))
  exclude = c(exclude, gsub("-immunosurveilling","-private", exclude))
  exclude = c(exclude, gsub("-private","-immunosurveilling", exclude))
  exclude = sort(unique(exclude))
  t[exclude]
  cell_types = setdiff(cell_types, exclude)
  #### pathways
  
  pws = colnames(raw_pathway_scores)
  pathways2 = pws[c(grep("KEGG",pws), grep("GO_",pws), grep("REACTOME_",pws), grep("HALLMARK_", pws))]
  list_mean_per_cell_type_pathway= NULL
  for(g in c(1:length(pathways2))){
    mean_per_cell_type_pathway = matrix(data = NA, nrow = length(pat_samples), ncol = length(cell_types), dimnames = c(list(pat_samples), list(cell_types)))
    score = raw_pathway_scores[,pathways2[g]]
    for(s in c(1:length(pat_samples))){
      print(s)
      for(c in c(1:length(cell_types))){
        w = intersect(which (pat_sample==pat_samples[s]), which(cell_type_clone_size==cell_types[c]))
        if(length(w)>=5){
          mean_per_cell_type_pathway[pat_samples[s], cell_types[c]]= mean(score[w])
        }
      }
    }
    list_mean_per_cell_type_pathway = c(list_mean_per_cell_type_pathway, list(mean_per_cell_type_pathway))
    print(g)
  }
  names(list_mean_per_cell_type_pathway) = pathways2
  
  #
  
  
  # ki67 expression = MKI67
  MKI67 = raw_pathway_scores[,"MKI67[rownames(pbmc@meta.data)]"]
  w = which(MKI67!=0)
  t0 = table(cell_type)
  t1 = t0*0
  ki67_per_cell_type = matrix(data = 0, nrow = length(pat_samples), ncol = length(cell_types), dimnames = c(list(pat_samples), list(cell_types)))
  for(s in c(1:length(pat_samples))){
    w1 = which (pat_sample==pat_samples[s])
    if(length(w1)>=10){
      w2 = intersect(w,w1)
      t1a = table(cell_type_clone_size[w2])
      t2a = table(cell_type_clone_size[w1])
      t1a = t1a[which(names(t1a) %in% cell_types)]
      t2a = t2a[which(names(t2a) %in% cell_types)]
      w3 = which(t2a<5)
      ki67_per_cell_type[pat_samples[s],names(t1a)] = t1a
      ki67_per_cell_type[pat_samples[s],names(t2a)] = ki67_per_cell_type[pat_samples[s],names(t2a)]*100/t2a
      ki67_per_cell_type[pat_samples[s],names(w3)] = NA
    }}
  
  list_mean_per_cell_type_pathway = c(list_mean_per_cell_type_pathway, list(ki67_per_cell_type))
  names(list_mean_per_cell_type_pathway)[length(names(list_mean_per_cell_type_pathway))] = "% Ki67+"
  
  # cell cycle prediction
  
  cell_cycle = raw_pathway_scores[,"pbmc@meta.data$Phase"]
  cell_cycles = sort(unique(cell_cycle))
  
  t0 = table(cell_cycle)
  t1 = t0*0
  S_per_cell_type = matrix(data = 0, nrow = length(pat_samples), ncol = length(cell_types), dimnames = c(list(pat_samples), list(cell_types)))
  w = which(cell_cycle=="S")
  for(s in c(1:length(pat_samples))){
    w1 = which (pat_sample==pat_samples[s])
    if(length(w1)>=10){
      w2 = intersect(w,w1)
      t1a = table(cell_type_clone_size[w2])
      t2a = table(cell_type_clone_size[w1])
      t1a = t1a[which(names(t1a) %in% cell_types)]
      t2a = t2a[which(names(t2a) %in% cell_types)]
      w3 = which(t2a<5)
      S_per_cell_type[pat_samples[s],names(t1a)] = t1a
      S_per_cell_type[pat_samples[s],names(t2a)] = S_per_cell_type[pat_samples[s],names(t2a)]*100/t2a
      S_per_cell_type[pat_samples[s],names(w3)] = NA
    }}
  
  list_mean_per_cell_type_pathway = c(list_mean_per_cell_type_pathway, list(S_per_cell_type))
  names(list_mean_per_cell_type_pathway)[length(names(list_mean_per_cell_type_pathway))] = "% cells in S phase"
  
  saveRDS(file = concat(c(output_directory, "apoptosis_UMAP_reclustering_PDAC150Ka_per_cell_type_immunosurveilling_pathway_level.rds")), list_mean_per_cell_type_pathway)
  
  
  ####
  list_counts=readRDS(file = concat(c(output_directory, "apoptosis_UMAP_reclustering_PDAC150Ka_per_cell_type_immunosurveilling_pathway_level.rds")))
  
  
  fileout1 = concat(c(output_directory,"B_T_cell_apoptosis_boxplot_counts_per_patient_grouped_per_subest_immunosurveilling.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1.6*2, width=w*2.1*2)
  par(mfrow= c(2,2), mar = c(25,5,2,3))
  
  for(subset in c(1:length(list_counts))){
    types_sub = c("B cell", "T cell")
    for (tt in c(1:length(types_sub))){
      name = names(list_counts)[subset]
      proportions = list_counts[[subset]]
      w = grep(types_sub[tt], colnames(proportions))
      mat_stat = proportions[group_PCA_list,w]
      w = which(apply(mat_stat, 2, function(x){length(which(is.na(x)==F))})>=10)
      mat_stat = mat_stat[,w]
      factor = factor(factor_PCA)
      mat_stat1 = mat_stat-min(mat_stat[which(is.na(mat_stat)==F)])
      fit = manova(formula = mat_stat1^1 ~ factor)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      
      i1 = 0
      for(i in p1){
        i1 = i1+1
        p_value = c(p_value, i$'Pr(>F)'[1]) 
        if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      print(min(p_value))
      
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.site = rep(concat(c("biopsy ", types_sub[tt])), length(p_value))
      p.analysis = rep(analysis, length(nam))
      x = cbind(p.site, p.analysis, combined_p_value)
      
      #out_file_table = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Bregs_Tregs/Breg_boxplot_counts_per_patient_grouped_tumour_.txt"))
      #write.table(x, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
      
      groups = NULL
      for(i in c(1:length(mat_stat[1,]))){
        x1 = mat_stat[gsub("_","-", groups_PCA[[1]],fixed=T),i]
        x2 = mat_stat[gsub("_","-", groups_PCA[[2]],fixed=T),i]
        x1 = x1[which(is.na(x1)==F)]
        x2 = x2[which(is.na(x2)==F)]
        g1 = c(list(x1),
               list(x2))
        names(g1) = c("ME","AE")
        groups = c(groups, list(g1))
      }
      names(groups)=colnames(mat_stat)
      
      library(RColorBrewer)
      cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
      cols =  add.alpha (cols1, alpha = 0.5)
      
      factors1 = c("ME","AE")
      factors = gsub("_"," ",colnames(mat_stat))
      main = concat(c(name," biopsy ", types_sub[tt]))
      max = max(c(unlist(groups), unlist(groups))*1.15)
      min = min(c(unlist(groups), unlist(groups))*1.15)
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 40
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(-100:100)*20}
      if(range<=50){scale = c(-100:100)*10}
      if(range<=30){scale = c(-100:100)*5}
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
      
      
      library(RColorBrewer)
      colx =  add.alpha (brewer.pal(8, "Set1")[c(2,1)], alpha = 0.95)
      
      col.names = rep("black", length(factors))
      col.names[grep("-immunosurveilling", factors)] = colx[1]
      col.names[grep("-private", factors)] = colx[2]
      mtext(side = 1, text = gsub("immunosurveilling","circ.",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1,col = col.names)
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
        pval1= "NS"
        if(p_value[i]<signif_threshold){pval1 = "*"
        y = max(unlist(groups[[i]]))
        y = y+3*b
        text(i, y+2*b, labels = pval1, cex = 1.3)
        }}
    }
  }
  dev.off()	
  
  #### blood
  
  fileout1 = concat(c(output_directory,"B_T_cell_apoptosis_boxplot_counts_per_patient_grouped_per_subest_immunosurveilling_blood.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1.6*2, width=w*2.1*2)
  par(mfrow= c(2,2), mar = c(25,5,2,3))
  
  for(subset in c(1:length(list_counts))){
    types_sub = c("B cell", "T cell")
    for (tt in c(1:length(types_sub))){
      name = names(list_counts)[subset]
      proportions = list_counts[[subset]]
      w = grep(types_sub[tt], colnames(proportions))
      mat_stat = proportions[gsub("biopsy", "blood",group_PCA_list),w]
      
      exclude = names(which(apply(mat_stat, 2, function(x){length(which(is.na(x)==F))})<9))
      exclude = c(exclude, gsub("-immunosurveilling","-private", exclude))
      exclude = c(exclude, gsub("-private","-immunosurveilling", exclude))
      exclude = sort(unique(exclude))
      
      mat_stat = mat_stat[,setdiff(colnames(mat_stat),exclude)]
      factor = factor(factor_PCA)
      mat_stat1 = mat_stat-min(mat_stat[which(is.na(mat_stat)==F)])
      fit = manova(formula = mat_stat1^1 ~ factor)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      
      i1 = 0
      for(i in p1){
        i1 = i1+1
        p_value = c(p_value, i$'Pr(>F)'[1]) 
        if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      print(min(p_value))
      
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.site = rep(concat(c("biopsy ", types_sub[tt])), length(p_value))
      p.analysis = rep(analysis, length(nam))
      x = cbind(p.site, p.analysis, combined_p_value)
      
      #out_file_table = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Bregs_Tregs/Breg_boxplot_counts_per_patient_grouped_tumour_.txt"))
      #write.table(x, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
      
      groups = NULL
      for(i in c(1:length(mat_stat[1,]))){
        x1 = mat_stat[gsub("_biopsy", "-blood", groups_PCA[[1]],fixed=T),i]
        x2 = mat_stat[gsub("_biopsy", "-blood", groups_PCA[[2]],fixed=T),i]
        x1 = x1[which(is.na(x1)==F)]
        x2 = x2[which(is.na(x2)==F)]
        g1 = c(list(x1),
               list(x2))
        names(g1) = c("ME","AE")
        groups = c(groups, list(g1))
      }
      names(groups)=colnames(mat_stat)
      
      library(RColorBrewer)
      cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
      cols =  add.alpha (cols1, alpha = 0.5)
      
      factors1 = c("ME","AE")
      factors = gsub("_"," ",colnames(mat_stat))
      main = concat(c(name," biopsy ", types_sub[tt]))
      max = max(c(unlist(groups), unlist(groups))*1.15)
      min = min(c(unlist(groups), unlist(groups))*1.15)
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 40
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(-100:100)*20}
      if(range<=50){scale = c(-100:100)*10}
      if(range<=30){scale = c(-100:100)*5}
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
      
      
      library(RColorBrewer)
      colx =  add.alpha (brewer.pal(8, "Set1")[c(2,1)], alpha = 0.95)
      
      col.names = rep("black", length(factors))
      col.names[grep("-immunosurveilling", factors)] = colx[1]
      col.names[grep("-private", factors)] = colx[2]
      mtext(side = 1, text = gsub("immunosurveilling","circ.",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1,col = col.names)
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
        pval1= "NS"
        if(p_value[i]<signif_threshold){pval1 = "*"
        y = max(unlist(groups[[i]]))
        y = y+3*b
        text(i, y+2*b, labels = pval1, cex = 1.3)
        }}
    }
  }
  dev.off()	
  
}
Boxplot_B_T_cell_biopsy_apoptosis_per_subset_immunosurveilling(input_directory, output_directory,groups_PCA)




  
   
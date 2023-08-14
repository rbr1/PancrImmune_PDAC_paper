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


Get_cell_counts_overall<-function(){
  files = c("Overall_UMAP_blood_COMBINED_PENG_STEELE_Non_immune_cells.txt",
            "Overall_UMAP_blood_COMBINED_PENG_STEELE_CAF.txt",
            "Overall_UMAP_blood_COMBINED_PENG_STEELE_B_cells.txt",
            "Overall_UMAP_blood_COMBINED_PDAC150K_Myeloid_cells.txt",
            "Overall_UMAP_blood_COMBINED_PDAC150K_T_cells.txt"
            )
  file_type = c("Non-immune","CAF","B cells","Myeloid","T cells")
  list_counts = NULL
  all_patients = NULL
  for(f in c(1:5)){
  #for(f in c(1:2)){
    p <- as.matrix(read.csv(concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/",files[f])), head=T, sep="\t"))
    p=p[which(rownames(p)!='-'),]
    t = table(p[,"cell.type"], p[,"sample"])
    colnames(t)
    colnames(t) = gsub("GSM47107PDAC_TISSUE_12","GSM4710701_PDAC_TISSUE_12",colnames(t))
    colnames(t) = gsub("GSM47107PDAC_TISSUE_13","GSM4710702_PDAC_TISSUE_13",colnames(t))
    colnames(t) = gsub("GSM47107PDAC_TISSUE_14","GSM4710703_PDAC_TISSUE_14",colnames(t))
    colnames(t) = gsub("GSM47107PDAC_TISSUE_15","GSM4710704_PDAC_TISSUE_15",colnames(t))
    colnames(t) = gsub("GSM47107PDAC_TISSUE_16","GSM4710705_PDAC_TISSUE_16",colnames(t))
    list_counts = c(list_counts, list(t))
    all_patients = c(all_patients, colnames(t))
    print(concat(c(f," ",length(colnames(t)))))
  }
  sort(table(all_patients))
  names(list_counts) = file_type[c(1:length(list_counts))]
  all_patients = sort(unique(all_patients))
  
  list_counts_full = NULL
  for(f in c(1:length(list_counts))){
    m = list_counts[[f]]
    all_cell_types = rownames(m)
    mat_counts = matrix(data = 0,nrow = length(all_cell_types),ncol = length(all_patients), dimnames = c(list(all_cell_types),list(all_patients)))
    mat_counts[rownames(m), colnames(m)] = m
    list_counts_full = c(list_counts_full, list(mat_counts))
    print (dim(mat_counts))
  }    
  names(list_counts_full) = names(list_counts)
  T_NK_cell = list_counts_full[["T cells"]]
  T_NK_cell[grep("T cell", rownames(T_NK_cell)),]
  
  broad_immune = rbind(colSums(list_counts_full[["B cells"]]), colSums(list_counts_full[["Myeloid"]]),colSums(T_NK_cell[grep("T cell", rownames(T_NK_cell)),]), colSums(T_NK_cell[grep("NK cell", rownames(T_NK_cell)),]))
  rownames(broad_immune) = c('B cells',"Myeloid cells","T cells","NK cells")
  
  broad_all = rbind(list_counts_full[["Non-immune"]], broad_immune)
  
  list_counts_full = c(list_counts_full, list(broad_immune), list(broad_all))
  names(list_counts_full) = c(names(list_counts) , "broad immune", "broad all")
  
 ##### print out
  out = NULL
  for(i in c(1:length(list_counts_full))){
    x = cbind(rownames(list_counts_full[[i]]), names(list_counts_full)[i], list_counts_full[[i]])
    if(length(out)==0){out = x
    }else{out = rbind(out, x)}
  }
  colnames(out)[1] = "Cell_type"
  colnames(out)[2] = "Group"	
  out_file_table = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/Seurat_cell_counts_COMBINED_PENG_STEELE_PDAC150K_all.txt"
  write.table(out, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}

output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
input_directory  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"
batch = "COMBINED_PENG_STEELE_PDAC150K"
library(psych)
library(corrplot)

PENG_data_plot<-function(){
  ######################
  batch = "2019_PUBLISHED_PDAC"
  batch1 = "PENG"
  type = "all"
  analysis = "PCA_cell_counts"
  analysis_all = "PCA_cell_counts"
  ######################
  clinical.info= readRDS(file = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC Chinese reanalysis/PLOTS_reannotation/PCA groups/Clinical_inforamtion.rds")
  
  ### Plot Peng
  w_include= clinical.info$samples_use
  
  sample = colnames(p)[-1]
  sample = sample[-1]
  sites = sample
  disease = sample
  sites[grep("PDAC_TISSUE", sites)] = "PDAC tissue"
  sites[grep("AdjNorm_TISSUE", sites)] = "AdjNorm tissue"
  sites[grep("PDAC_PBMC", sites)] = "PDAC PBMC"
  sites[grep("Healthy_PBMC", sites)] = "Healthy PBMC"
  names(sites) = sample
  
  cell_types = p[,1]
  cell_group = p[,2]
  cell_groups = sort(unique(cell_group))
  cell_groups = cell_groups[which(cell_groups!="Overall")]
  list_normalised = NULL
  
  include = NULL
  totals = NULL
  for(ind in c(1:length(cell_groups))){
    w = which(cell_group== cell_groups[ind])
    if(cell_groups[ind] =="CD45:broad"){
      w = intersect(w, which(p[,1] %in% c( "B cell"  , "Myeloid","NK" , "T cell"   )))
    }
    cell_type_sub = p[w,1]
    p1 = p[w, w_include]
    m_count_norm = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(w_include), dimnames = c(list(cell_type_sub), list(w_include)))
    m_total= matrix(data = 0,nrow = length(cell_type_sub), ncol = length(w_include), dimnames = c(list(cell_type_sub), list(w_include)))
    for(s in c(1:length(w_include))){
      m_count_norm[, w_include[s]]= as.numeric(p1[, w_include[s]])*100/sum(as.numeric(p1[, w_include[s]]))
      m_total[, w_include[s]]=as.numeric(p1[, w_include[s]])
      if(sum(as.numeric(p1[, w_include[s]]))<7){m_count_norm[, w_include[s]] = -1}
    }
    if(length(which(apply(m_count_norm, 2, max)!=-1))>3){include = c(include, ind)}
    m_count_norm= m_count_norm[order(rownames(m_count_norm)),]
    #m_count_norm  = m_count_norm [,rownames(group_mat)]
    list_normalised= c(list_normalised, list(m_count_norm))
    totals = c(totals, list(m_total))
  }
  names(list_normalised) = cell_groups
  list_normalised = list_normalised[include]
  names(totals) = cell_groups
  
  
  for(i in c(1:length(list_normalised))){
    print (names(list_normalised)[i])
    print (rownames(list_normalised[[i]]))
  }
  
  T_cells = t(list_normalised[[ "T cells"]])
  Treg_CD8_ratio = rowSums(T_cells [,grep("Treg",colnames(T_cells))])/rowSums(T_cells [,grep("CD8",colnames(T_cells))])
  
  T_cells = t(list_normalised[["T cells"]])
  CD4_CD8_ratio = rowSums(T_cells [,grep("CD4",colnames(T_cells))])/rowSums(T_cells [,grep("CD8",colnames(T_cells))])
  
  cd4cd8 = c(list(Treg_CD8_ratio), list(CD4_CD8_ratio))
  names(cd4cd8) = c("Treg-CD8 T cell ratio", "CD4-CD8 ratio")
  
  T_cells = t(list_normalised[["T cells"]])
  Treg_act_ratio = T_cells [,"T cell CD4 Treg Activated"]/rowSums(T_cells [,grep("Treg",colnames(T_cells))])
  
  T_ratios = cbind(Treg_CD8_ratio, CD4_CD8_ratio, Treg_act_ratio)
  
  samples_run= colnames(list_normalised[[1]])
  
  
  
  ############################################
  all_proportions = t(list_normalised[["broad immune"]])
  
  ########### split by stage
  plot_ids1 = clinical.info$samples_use[which(clinical.info$Staging %in% "IB")]
  plot_ids2 = clinical.info$samples_use[which(clinical.info$Staging %in% c( "IIA","IIB","III"))]
  plot_ids_all = c(list(plot_ids1), list(plot_ids2))
  names(plot_ids_all) = c("Stage IB-IIA", "Stage IIB-III")
  
  groups_PCA = readRDS(file = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC Chinese reanalysis/PDAC150K validation/Plots for PDAC150K validation/Cell Annotations/Immunogroups_2019_PUBLISHED_PDAC.rds")
  group1 = groups_PCA[[1]]
  group2 = groups_PCA[[2]]
  
  all_proportions = t(list_normalised[["broad all"]])
  x = all_proportions[,"T cells"]+all_proportions[,"B cells"]
  
  #group1 = names(which(x > quantile(x, 0.5)))
  #group2 = names(which(x <= quantile(x, 0.5)))
  #groups_PCA = c(list(group1), list(group2))
  #groups_PCA1 = groups_PCA

  
  g1 = c(list(intersect(plot_ids1, groups_PCA[[1]])), list(intersect(plot_ids1, groups_PCA[[2]])))
  g2 = c(list(intersect(plot_ids2, groups_PCA[[1]])), list(intersect(plot_ids2, groups_PCA[[2]])))
  
  ga = c(list(g1[[1]]), list(g2[[1]]))
  gb = c(list(g1[[2]]), list(g2[[2]]))
  names(ga) = names(plot_ids_all)
  names(gb) = names(plot_ids_all)
  groups_PCR_by_stage = c(list(ga), list(gb))
  names(groups_PCR_by_stage) = c("ME", "AE")
  
  ############################################ plot PCA of CD45 %s
  Plot_PCA<-function(list_normalised){
    all_proportions = (list_normalised[["broad immune" ]][,c(group1,group2)])
    x <- princomp(t(all_proportions),dim=2)
    PCA_coordinates_sub = x$scores
    
    x1 = x$scores[,1]
    x2 = x$scores[,2]
    
    #x1 = x$loadings[,1]
    #x2 = x$loadings[,2]
    
    factor1 = rep(0,length(x1))
    names(factor1) = names(x1)
    factor1[group1] = 1
    factor1[group2] = 2
    factor2 = gsub("-biopsy","",names(x1))
    
    factors1 = sort(unique(factor1))
    factors2 = sort(unique(factor2))
    matching1 = match(factor1, factors1)
    matching2 = match(factor2, factors2)
    names(matching1) = names(factor1)
    names(matching2) = names(factor2)
    
    library(RColorBrewer)
    cols =  add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95)
    cols1 = add.alpha (c( cols[1],"white"),alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)
    
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)
    
    pches = c(21:25)
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_scores_", batch1,".pdf"))
    w=2.2
    pdf(file=fileout1, height=w*2, width=w*3)
    par(mfrow= c(2,3), mar = c(5,5,3,3))
    
    
    # stage IB
    wx = c(1:length(x1))
    main= "All"
    plot(x1[wx],x2[wx],pch = pches[matching1[wx]], col = cols[matching1[wx]],bg = cols2[matching1[wx]], main = main, xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
    
    xy = cbind(x1,x2)
    library(car)
    for(pca in c(1:length(groups_PCA))){
      w = which(factor1[wx] ==pca)
      if(length(w)>1){
        x = xy[wx[w],1]
        y = xy[wx[w],2]
        dataEllipse(x, y, levels=c(0.65), add = TRUE, col = cols2[pca], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
      }
    }
    
    
    # stage IB
    wx = unlist(g1)
    main= "Stage IB"
    plot(x1[wx],x2[wx],pch = pches[matching1[wx]], col = cols[matching1[wx]],bg = cols2[matching1[wx]], main = main, xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
    
    xy = cbind(x1,x2)
    library(car)
    for(pca in c(1:length(groups_PCA))){
      w = which(factor1[wx] ==pca)
      if(length(w)>1){
        x = xy[wx[w],1]
        y = xy[wx[w],2]
        dataEllipse(x, y, levels=c(0.65), add = TRUE, col = cols2[pca], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
      }
    }
    
    # stage II-III
    wx = unlist(g2)
    main= "Stage II-III"
    plot(x1[wx],x2[wx],pch = pches[matching1[wx]], col = cols[matching1[wx]],bg = cols2[matching1[wx]], main = main, xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
    
    xy = cbind(x1,x2)
    library(car)
    for(pca in c(1:length(groups_PCA))){
      w = which(factor1[wx] ==pca)
      if(length(w)>1){
        x = xy[wx[w],1]
        y = xy[wx[w],2]
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
  Plot_tumour_cell_composition_between_groups<-function(list_normalised){
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_", batch1,".pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*6, width=w*3*2.9)
    par(mfrow= c(6,2), mar = c(16,5,2,3))
    
    summary_tables = NULL
    analysis1 = "cell_proportions"
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(c in c(1:length(list_normalised))){
      mat = t(list_normalised[[c]])
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
      
      ##### biopsy
      w = which(is.na(mat_biopsy[,1])==F)
      mat1 = mat_biopsy[group_PCA_list[w],]
      factor = factor(factor_PCA[w])
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
      mat = mat_biopsy
      p_values = combined_p_value[,1]
      for(i in c(1:length(mat[1,]))){
        g1 = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(groups_PCA))){
          x = mat [groups_PCA[[g]], i]
          x = x[which(is.na(x)!=T)]
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
  summary_tables = Plot_tumour_cell_composition_between_groups(list_normalised)
  
  outfile = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/Stats_PCA_cell_types_by_groups_",batch1,".txt"))
  write.table(summary_tables, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ############################################ plot FC versus p-values of cell %s in tumour between groups
  Plot_FC_p_value<-function(summary_tables){
    summary_tables[which(as.numeric(summary_tables[,"p_value"])<0.05),]
    FC = log10(as.numeric(summary_tables[,"mean.group. 1"])/as.numeric(summary_tables[,"mean.group. 2"]))
    pval = -log10(as.numeric(summary_tables[,"p_value"]))
    plot(FC, pval)
    
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    types_include = c("cell_proportions_B cells","cell_proportions_Myeloid","cell_proportions_T cells")
    w = intersect(which(summary_tables[,"p.site"]=="biopsy"), which(summary_tables[,"p.analysis"] %in% types_include))
    pval = as.numeric(summary_tables[w,"p_value"])
    
    names =rownames(summary_tables[w])
    names = cbind(names, summary_tables[w,"p.analysis"])
    names[,1] = gsub("cell_proportions_", "% of ",names[,1])
    names[,1] = gsub("_"," ",names[,1])
    rownames(names) = gsub("  "," ",rownames(names))
    rownames(names) = gsub("_"," ",rownames(names), fixed = T)
    names = apply(cbind(rownames(names), " ",names,""), 1, paste, collapse = "")
    cbind(names)
    
    FC = as.numeric(summary_tables[w,"mean.group. 2"]) - as.numeric(summary_tables[w,"mean.group. 1"])
    x = FC
    y = -log10(pval)
    names(x) = names
    names(y) = names
    b = max(y)*0.1
    range = max(abs(FC))*4
    cols = rep(NA, length(pval))
    cols[which(FC <0)] = cols1[1]
    cols[which(FC >0)] = cols1[2]
    w = which(pval>0.05)
    cols[w] = add.alpha(cols[w], alpha = 0.4)
    names(cols)= names
    
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_FC_plot_", batch1,".pdf"))
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
      text( 0.35*(range/1.9), seq, labels = gsub("DC CD207","CD207",gsub("memory","mem.",gsub("_"," ",gsub("Myeloid ","",order) ))), cex= cex-0.2,las = 1, font = 1, col = cols[order], pos = 4)
      # mtext(side = 4, text = gsub("_"," ",order), line = 0.2,cex= cex-0.2, at = seq, las = 1, font = 1, col = cols[order])
      for(i in c(1:length(order))){
        segments(x[order[i]], y[order[i]], 0.4*(range/1.9), seq[i], lwd = 1.5, lty = 1, col = add.alpha ("grey", alpha = 0.5))
      }
    }
    w = intersect(which(pval<0.05), which (FC<0))
    if(length(w)>0){
      order = names[w][order(y[w])]
      cex =0.9
      range1 = range(y)
      seq =(seq(from = -log10(0.035),to = 4.9, length = length(order)))
      text( -0.35*(range/1.9), seq, labels = gsub("DC CD207","CD207",gsub("memory","mem.",gsub("_"," ",gsub("Myeloid ","",order) ))), cex= cex-0.2,las = 1, font = 1, col = cols[order], pos = 2)
      # mtext(side = 4, text = gsub("_"," ",order), line = 0.2,cex= cex-0.2, at = seq, las = 1, font = 1, col = cols[order])
      for(i in c(1:length(order))){
        segments(x[order[i]], y[order[i]], -0.4*(range/1.9), seq[i], lwd = 1.5, lty = 1, col = add.alpha ("grey", alpha = 0.5))
      }
    }
    dev.off()
    
  }
  Plot_FC_p_value(summary_tables)  
  
  ############################################ plot boxplots of cell %s in tumour between groups
  ### stage IB
  Plot_tumour_cell_composition_between_groups<-function(list_normalised){
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_", batch1,"_StageIB.pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*6, width=w*3*2.9)
    par(mfrow= c(6,2), mar = c(16,5,2,3))
    
    summary_tables = NULL
    analysis1 = "cell_proportions"
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(c in c(1:length(list_normalised))){
      mat = t(list_normalised[[c]])
      g1 = c(list(intersect(plot_ids1, groups_PCA1[[1]])), list(intersect(plot_ids1, groups_PCA1[[2]])))
      g2 = c(list(intersect(plot_ids2, groups_PCA1[[1]])), list(intersect(plot_ids2, groups_PCA1[[2]])))
      groups_PCA = g1
      
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
      
      ##### biopsy
      w = which(is.na(mat_biopsy[,1])==F)
      mat1 = mat_biopsy[group_PCA_list[w],]
      factor = factor(factor_PCA[w])
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
      mat = mat_biopsy
      p_values = combined_p_value[,1]
      for(i in c(1:length(mat[1,]))){
        g1 = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(groups_PCA))){
          x = mat [groups_PCA[[g]], i]
          x = x[which(is.na(x)!=T)]
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
    
    
    
    #####
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_", batch1,"_StageII_III.pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*6, width=w*3*2.9)
    par(mfrow= c(6,2), mar = c(16,5,2,3))
    
    summary_tables = NULL
    analysis1 = "cell_proportions"
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(c in c(1:length(list_normalised))){
      mat = t(list_normalised[[c]])
      g1 = c(list(intersect(plot_ids1, groups_PCA1[[1]])), list(intersect(plot_ids1, groups_PCA1[[2]])))
      g2 = c(list(intersect(plot_ids2, groups_PCA1[[1]])), list(intersect(plot_ids2, groups_PCA1[[2]])))
      groups_PCA = g2
      
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
      
      ##### biopsy
      w = which(is.na(mat_biopsy[,1])==F)
      mat1 = mat_biopsy[group_PCA_list[w],]
      factor = factor(factor_PCA[w])
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
      mat = mat_biopsy
      p_values = combined_p_value[,1]
      for(i in c(1:length(mat[1,]))){
        g1 = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(groups_PCA))){
          x = mat [groups_PCA[[g]], i]
          x = x[which(is.na(x)!=T)]
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
  summary_tables = Plot_tumour_cell_composition_between_groups(list_normalised)
  
  outfile = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/Stats_PCA_cell_types_by_groups_PDAC150Ka.txt"))
  write.table(summary_tables, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ############################################ plot FC versus p-values of cell %s in tumour between groups
  Plot_FC_p_value<-function(summary_tables){
    summary_tables[which(as.numeric(summary_tables[,"p_value"])<0.05),]
    FC = log10(as.numeric(summary_tables[,"mean.group. 1"])/as.numeric(summary_tables[,"mean.group. 2"]))
    pval = -log10(as.numeric(summary_tables[,"p_value"]))
    plot(FC, pval)
    
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    types_include = c("cell_proportions_B cells","cell_proportions_Myeloid","cell_proportions_T cells")
    w = intersect(which(summary_tables[,"p.site"]=="biopsy"), which(summary_tables[,"p.analysis"] %in% types_include))
    pval = as.numeric(summary_tables[w,"p_value"])
    
    names =rownames(summary_tables[w])
    names = cbind(names, summary_tables[w,"p.analysis"])
    names[,1] = gsub("cell_proportions_", "% of ",names[,1])
    names[,1] = gsub("_"," ",names[,1])
    rownames(names) = gsub("  "," ",rownames(names))
    rownames(names) = gsub("_"," ",rownames(names), fixed = T)
    names = apply(cbind(rownames(names), " ",names,"s"), 1, paste, collapse = "")
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
    
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_FC_plot_", batch1,".pdf"))
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
  Plot_FC_p_value(summary_tables)  
  
  
  ############################################ plot boxplots of tumour activated Treg of all Tregs %
  Plot_cell_act_Treg_ratio_tumour<-function(T_ratios){
    stage1 = clinical.info$samples_use[which(clinical.info$Staging %in% c("IB"))]
    stage2 = clinical.info$samples_use[which(clinical.info$Staging %in% c("IIB"))]
    stages = c(list(stage1),list(stage2))
    names(stages) = c("Stage IB","Stage IIB")
    analysis_all = "Treg_act_of_Tregs_proportions"
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_", batch,".pdf"))
    w=4
    pdf(file=fileout1, height=w*0.65*1, width=w*1*0.8)
    par(mfrow= c(1,3), mar = c(6,4,2,1))
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(s in c(1:length(stages))){
      ids_use = stages[[s]]
      for(c in c(1:length(T_ratios[1,]))){
        mat = T_ratios[,c]
        name = concat(c(colnames(T_ratios)[c], " ",names(stages)[s]))
        groups = NULL
        for(pca in c(1:length(groups_PCA))){
          groups = c(groups, list(mat[intersect(gsub("-","_",groups_PCA[[pca]], fixed = T),ids_use)]))
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
    }
   
    dev.off()
    
  }
  Plot_cell_act_Treg_ratio_tumour(T_ratios)
  
  
}

STEELE_data_plot<-function(){
  
  file = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/Seurat_cell_counts_COMBINED_PENG_STEELE_PDAC150K_all.txt"))
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  ####################################################### subset to FNAs
  ########## STEELE et al
  batch = "PDAC_GSE155698"
  batch1 = "STEELE"
  type = "all"
  analysis = "PCA_cell_counts"
  analysis_all = "PCA_cell_counts"
  ####################################################### subset to FNAs
  file = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC_GSE155698/Clinical data/6a_SampleAttributes_DS_pasca_v4.txt"
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[which(p[,"Tissue.GEO.ID"] %in% colnames(p1)),]
  p=p[which(p[,"tissue"]=="Tumor"),]
  p=p[which(p[,"Tissue.GEO.ID"]!=''),]
  sample = p[,"Tissue.GEO.ID"]
  sample1 = p[,"PBMC.GEO.ID"]
  sex =p[,"sex"]
  sample_type =p[,"sample_type"]
  sample_type[which(sample_type=="Surgical\xca" )] = "Surgical"
  stage =p[,"stage"]
  stage[grep("sectable", stage)] = "Resectable"
  recist_category =p[,"recist_category"]
  recist_category[which(recist_category %in% c("Complete response","Non-progressive disease"))]="CR/PR"
  group_mat = cbind(sex, sample_type, stage, recist_category, sample1)
  rownames(group_mat) = sample
  w_include = rownames(group_mat)[which(group_mat[,"sample_type"]!="Surgical" )]
  group_mat = group_mat[w_include,]
  
  p1 = p[which(p[,"Tissue.GEO.ID"] %in% w_include),]
  p1 = p1[match(p1[,"Tissue.GEO.ID"], w_include),]
  w_include==p1[,"Tissue.GEO.ID"]
  w_include_tumour = w_include
  w_include_blood = p1[,"PBMC.GEO.ID"]
  w_include = c(w_include, w_include_blood)
  ####################
  
  sample = colnames(p)[-1]
  sample = sample[-1]
  sites = sample
  disease = sample
  sites[grep("PDAC_TISSUE", sites)] = "PDAC tissue"
  sites[grep("AdjNorm_TISSUE", sites)] = "AdjNorm tissue"
  sites[grep("PDAC_PBMC", sites)] = "PDAC PBMC"
  sites[grep("Healthy_PBMC", sites)] = "Healthy PBMC"
  names(sites) = sample
  
  list_normalised = NULL
  file = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/Seurat_cell_counts_COMBINED_PENG_STEELE_PDAC150K_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  cell_types = p[,1]
  cell_group = p[,2]
  cell_groups = sort(unique(cell_group))
  cell_groups = cell_groups[which(cell_groups!="Overall")]

    
 
  w_include=w_include_tumour
  include = NULL
  for(ind in c(1:length(cell_groups))){
    w = which(cell_group== cell_groups[ind])
    if(cell_groups[ind] =="CD45:broad"){
      w = intersect(w, which(p[,1] %in% c( "B cell"  , "Myeloid","NK" , "T cell"   )))
    }
    cell_type_sub = p[w,1]
    p1 = p[w, w_include]
    m_count_norm = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(w_include), dimnames = c(list(cell_type_sub), list(w_include)))
    
    for(s in c(1:length(w_include))){
      m_count_norm[, w_include[s]]= as.numeric(p1[, w_include[s]])*100/sum(as.numeric(p1[, w_include[s]]))
      if(sum(as.numeric(p1[, w_include[s]]))<7){m_count_norm[, w_include[s]] = -1}
    }
    if(length(which(apply(m_count_norm, 2, max)!=-1))>3){include = c(include, ind)}
    m_count_norm= m_count_norm[order(rownames(m_count_norm)),]
    #m_count_norm  = m_count_norm [,rownames(group_mat)]
    list_normalised= c(list_normalised, list(m_count_norm))
  }
  names(list_normalised) = cell_groups
  list_normalised = list_normalised[include]
  
  
  for(i in c(1:length(list_normalised))){
    print (names(list_normalised)[i])
    print (rownames(list_normalised[[i]]))
  }
  
  '''
x = list_normalised[["CD45:broad"]]
B_T_cell = colSums(x[c("B cell", "T cell"),])
x = rbind(x[c("Myeloid", "NK"  ),], B_T_cell)
rownames(x)[3] = "B/T cells"
list_normalised = c(list_normalised, list(x))
names(list_normalised)[length(list_normalised)] = "CD45:broad condensed"

x = list_normalised[["All"]]
B_T_cell = colSums(x[c("Acinar", "Ductal 1", "Ductal 2" ),])
x = rbind(x[c("Acinar", "Ductal 1", "Ductal 2"  ),])
list_normalised = c(list_normalised, list(x))
names(list_normalised)[length(list_normalised)] = "All:tumour"
'''
  
  T_cells = t(list_normalised[[ "T cells"]])
  Treg_CD8_ratio = rowSums(T_cells [,grep("Treg",colnames(T_cells))])/rowSums(T_cells [,grep("CD8",colnames(T_cells))])
  
  T_cells = t(list_normalised[["T cells"]])
  CD4_CD8_ratio = rowSums(T_cells [,grep("CD4",colnames(T_cells))])/rowSums(T_cells [,grep("CD8",colnames(T_cells))])
  
  cd4cd8 = c(list(Treg_CD8_ratio), list(CD4_CD8_ratio))
  names(cd4cd8) = c("Treg-CD8 T cell ratio", "CD4-CD8 ratio")
  
  
  samples_run= colnames(list_normalised[[1]])
  
  
  ######### ratios of T cells
  groups_PCA = readRDS(file = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC_GSE155698/Plots/Immunogroups_PDAC_GSE155698.rds")
  
  '''
  all_proportions = t(list_normalised[["broad immune" ]])
  x = all_proportions[,"B cell"]+all_proportions[,"T cell"]
  group2 = names(which(x[w_include_tumour] > quantile(x[w_include_tumour], 0.6)))
  group1 = names(which(x[w_include_tumour] <= quantile(x[w_include_tumour], 0.6)))
  #group1 = names(which(all_proportions[w_include_tumour,"Myeloid"] > quantile(all_proportions[w_include_tumour,"Myeloid"], 0.5)))
  #group2 = names(which(all_proportions[w_include_tumour,"Myeloid"] <= quantile(all_proportions[w_include_tumour,"Myeloid"], 0.5)))
  groups_PCA = c(list(group1), list(group2))
  groups_PCA1 = groups_PCA
  '''

  groups_PCA_blood = groups_PCA
  groups_PCA_blood[[1]] = w_include_blood [match(groups_PCA[[1]], w_include_tumour)]
  groups_PCA_blood[[2]] = w_include_blood [match(groups_PCA[[2]], w_include_tumour)]
  
  g1 = groups_PCA
  g2 = groups_PCA_blood
  group1 = groups_PCA[[1]]#, groups_PCA_blood[[1]])
  group2 = groups_PCA[[2]]#, groups_PCA_blood[[2]])
  
  
  ############################################ plot PCA of CD45 %s
  Plot_PCA<-function(list_normalised){
    all_proportions = (list_normalised[["broad immune" ]][,c(group1,group2)])
    x <- princomp(t(all_proportions),dim=2)
    PCA_coordinates_sub = x$scores
    
    x1 = x$scores[,1]
    x2 = x$scores[,2]
    
    #x1 = x$loadings[,1]
    #x2 = x$loadings[,2]
    
    factor1 = rep(0,length(x1))
    names(factor1) = names(x1)
    factor1[group1] = 1
    factor1[group2] = 2
    factor2 = gsub("-biopsy","",names(x1))
    
    factors1 = sort(unique(factor1))
    factors2 = sort(unique(factor2))
    matching1 = match(factor1, factors1)
    matching2 = match(factor2, factors2)
    names(matching1) = names(factor1)
    names(matching2) = names(factor2)
    
    library(RColorBrewer)
    cols =  add.alpha (brewer.pal(6, "Dark2"), alpha = 0.95)
    cols1 = add.alpha (c( cols[1],"white"),alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)
    
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)
    
    pches = c(21:25)
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_scores_", batch1,".pdf"))
    w=2.2
    pdf(file=fileout1, height=w*2, width=w*3)
    par(mfrow= c(2,3), mar = c(5,5,3,3))
    
    # stage IB
    wx = unlist(g1)
    main= "Stage IB"
    plot(x1[wx],x2[wx],pch = pches[matching1[wx]], col = cols[matching1[wx]],bg = cols2[matching1[wx]], main = main, xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
    
    xy = cbind(x1,x2)
    library(car)
    for(pca in c(1:length(groups_PCA))){
      w = which(factor1[wx] ==pca)
      if(length(w)>1){
        x = xy[wx[w],1]
        y = xy[wx[w],2]
        dataEllipse(x, y, levels=c(0.65), add = TRUE, col = cols2[pca], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
      }
    }
    
    # stage II-III
    wx = unlist(g2)
    main= "Stage II-III"
    plot(x1[wx],x2[wx],pch = pches[matching1[wx]], col = cols[matching1[wx]],bg = cols2[matching1[wx]], main = main, xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
    
    xy = cbind(x1,x2)
    library(car)
    for(pca in c(1:length(groups_PCA))){
      w = which(factor1[wx] ==pca)
      if(length(w)>1){
        x = xy[wx[w],1]
        y = xy[wx[w],2]
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
  Plot_tumour_cell_composition_between_groups<-function(list_normalised){
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_", batch1,".pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*6, width=w*3*2.9)
    par(mfrow= c(6,2), mar = c(16,5,2,3))
    
    summary_tables = NULL
    analysis1 = "cell_proportions"
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(c in c(1:length(list_normalised))){
      mat = t(list_normalised[[c]])
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
      
      ##### biopsy
      w = which(is.na(mat_biopsy[,1])==F)
      mat1 = mat_biopsy[group_PCA_list[w],]
      factor = factor(factor_PCA[w])
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
      mat = mat_biopsy
      p_values = combined_p_value[,1]
      for(i in c(1:length(mat[1,]))){
        g1 = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(groups_PCA))){
          x = mat [groups_PCA[[g]], i]
          x = x[which(is.na(x)!=T)]
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
  summary_tables = Plot_tumour_cell_composition_between_groups(list_normalised)
  
  outfile = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/Stats_PCA_cell_types_by_groups_",batch1,".txt"))
  write.table(summary_tables, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ############################################ plot FC versus p-values of cell %s in tumour between groups
  Plot_FC_p_value<-function(summary_tables){
    summary_tables[which(as.numeric(summary_tables[,"p_value"])<0.05),]
    FC = log10(as.numeric(summary_tables[,"mean.group. 1"])/as.numeric(summary_tables[,"mean.group. 2"]))
    pval = -log10(as.numeric(summary_tables[,"p_value"]))
    plot(FC, pval)
    
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    types_include = c("cell_proportions_B cells","cell_proportions_Myeloid","cell_proportions_T cells")
    w = intersect(which(summary_tables[,"p.site"]=="biopsy"), which(summary_tables[,"p.analysis"] %in% types_include))
    pval = as.numeric(summary_tables[w,"p_value"])
    
    names =rownames(summary_tables[w])
    names = cbind(names, summary_tables[w,"p.analysis"])
    names[,1] = gsub("cell_proportions_", "% of ",names[,1])
    names[,1] = gsub("_"," ",names[,1])
    rownames(names) = gsub("  "," ",rownames(names))
    rownames(names) = gsub("_"," ",rownames(names), fixed = T)
    names = apply(cbind(rownames(names), " ",names,""), 1, paste, collapse = "")
    cbind(names)
    
    FC = as.numeric(summary_tables[w,"mean.group. 2"]) - as.numeric(summary_tables[w,"mean.group. 1"])
    x = FC
    y = -log10(pval)
    names(x) = names
    names(y) = names
    b = max(y)*0.1
    range = max(abs(FC))*4
    cols = rep(NA, length(pval))
    cols[which(FC <0)] = cols1[1]
    cols[which(FC >0)] = cols1[2]
    w = which(pval>0.05)
    cols[w] = add.alpha(cols[w], alpha = 0.4)
    names(cols)= names
    
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Cross-dataset analyses/", analysis_all,"_biopsy_PCA_groups_boxplots_cell_proportions_FC_plot_", batch1,".pdf"))
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
      text( 0.35*(range/1.9), seq, labels = gsub("DC CD207","CD207",gsub("memory","mem.",gsub("_"," ",gsub("Myeloid ","",order) ))), cex= cex-0.2,las = 1, font = 1, col = cols[order], pos = 4)
      # mtext(side = 4, text = gsub("_"," ",order), line = 0.2,cex= cex-0.2, at = seq, las = 1, font = 1, col = cols[order])
      for(i in c(1:length(order))){
        segments(x[order[i]], y[order[i]], 0.4*(range/1.9), seq[i], lwd = 1.5, lty = 1, col = add.alpha ("grey", alpha = 0.5))
      }
    }
    w = intersect(which(pval<0.05), which (FC<0))
    if(length(w)>0){
      order = names[w][order(y[w])]
      cex =0.9
      range1 = range(y)
      seq =(seq(from = -log10(0.035),to = 4.9, length = length(order)))
      text( -0.35*(range/1.9), seq, labels = gsub("DC CD207","CD207",gsub("memory","mem.",gsub("_"," ",gsub("Myeloid ","",order) ))), cex= cex-0.2,las = 1, font = 1, col = cols[order], pos = 2)
      # mtext(side = 4, text = gsub("_"," ",order), line = 0.2,cex= cex-0.2, at = seq, las = 1, font = 1, col = cols[order])
      for(i in c(1:length(order))){
        segments(x[order[i]], y[order[i]], -0.4*(range/1.9), seq[i], lwd = 1.5, lty = 1, col = add.alpha ("grey", alpha = 0.5))
      }
    }
    dev.off()
    
  }
  Plot_FC_p_value(summary_tables)  
  
  
}

Plot_correlations_between_all<-function(input_directory, output_directory){
  batch = "COMBINED_PENG_STEELE_PDAC150K"
  library(psych)
  library(corrplot)
  file = concat(c(input_directory,"Seurat_cell_counts_COMBINED_PENG_STEELE_PDAC150K_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  groups_PCA1 = readRDS(file = concat(c(input_directory, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))
  groups_PCA2 = readRDS(file = concat(c(input_directory, "Immunogroups_2019_PUBLISHED_PDAC.rds")))
  groups_PCA3 = readRDS(file = concat(c(input_directory,"Immunogroups_PDAC_GSE155698.rds")))
  
  w_include = unlist(c(groups_PCA1, groups_PCA2,groups_PCA3))
  
  cell_types = p[,1]
  cell_group = p[,2]
  cell_groups = sort(unique(cell_group))
  cell_groups = cell_groups[which(cell_groups!="Overall")]
  list_normalised = NULL
  
  include = NULL
  for(ind in c(1:length(cell_groups))){
    w = which(cell_group== cell_groups[ind])
    if(cell_groups[ind] =="CD45:broad"){
      w = intersect(w, which(p[,1] %in% c( "B cell"  , "Myeloid","NK" , "T cell"   )))
    }
    cell_type_sub = p[w,1]
    p1 = p[w, w_include]
    m_count_norm = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(w_include), dimnames = c(list(cell_type_sub), list(w_include)))
    
    for(s in c(1:length(w_include))){
      m_count_norm[, w_include[s]]= as.numeric(p1[, w_include[s]])*100/sum(as.numeric(p1[, w_include[s]]))
      if(sum(as.numeric(p1[, w_include[s]]))<7){m_count_norm[, w_include[s]] = -1}
    }
    if(length(which(apply(m_count_norm, 2, max)!=-1))>3){include = c(include, ind)}
    m_count_norm= m_count_norm[order(rownames(m_count_norm)),]
    #m_count_norm  = m_count_norm [,rownames(group_mat)]
    list_normalised= c(list_normalised, list(m_count_norm))
  }
  names(list_normalised) = cell_groups
  list_normalised = list_normalised[include]
  
  
  for(i in c(1:length(list_normalised))){
    print (names(list_normalised)[i])
    print (rownames(list_normalised[[i]]))
  }  
  
  
  # plot correlations
  library(RColorBrewer)
  library(psych)
  library(corrplot)
  for(i1 in c(1:length(list_normalised))){
    mat = list_normalised[[i1]]
    name = names(list_normalised)[i1]
    if(name %in% c( "broad all" , "CAF",  "Non-immune" )){
      use = colnames(mat)[grep("_biopsy", colnames(mat), invert = T)]
      mat = mat[,use]
    }
    a = apply(mat, 1, function(x){length(which(x!=0))})
    mat = mat[which(a>=1), ]
    print(length(mat[,1]))
    cortest = corr.test(t(mat),adjust="holm")
    pval = cortest $ p
    rval = cortest $ r
    fileout1=concat(c(output_directory,"Correlation_cell_types_pairs_plot_biopsy_overall_", batch,"_",name,".pdf"))
    w=3.5
    if(length(mat[,1])>15){w=8}
    if(length(mat[,1])<7){w=1.7}
    
    pdf(file=fileout1, height=w*1.1, width=w*1.1)
    par(mfrow= c(1,1), mar = c(5,5,3,3))
    corrplot(rval, type="upper", p.mat=pval, insig="label_sig", tl.pos="td", sig.level=0.05, title ="", method = "ellipse", diag = F, order = "hclust",, cl.pos = 'n')
    dev.off()
    
  }
  
  ### specific cell types
  cell_types_of_interest = c("Total T cell CD4 Treg Activated","Total T cell CD4 Treg","Total CD4 Treg","Total T cell CD8 EM","Total Myeloid DC mregDC","Total Myeloid momac","Total B cell plasma cell","Total ASC","Total B cell GC")
  cell_groups_of_interest = c("T cells","T cells","T cells","T cells","Myeloid cells","Myeloid cells","B cells","B cells","B cells")
  w1 = grep("momac", rownames(list_normalised[["Myeloid"]]))
  w2 = grep("plasma", rownames(list_normalised[["B cells"]]))
  w3 = grep("GC", rownames(list_normalised[["B cells"]]))
  specific_cell_frequencies = cbind(list_normalised[["T cells"]]["T cell CD4 Treg Activated",], list_normalised[["T cells"]]["T cell CD4 Treg",],
        colSums(list_normalised[["T cells"]][c("T cell CD4 Treg Activated","T cell CD4 Treg"),]),
        list_normalised[["T cells"]]["T cell CD8 EM",],
        list_normalised[["Myeloid"]]["Myeloid DC mregDC",],
        colSums(list_normalised[["Myeloid"]][ w1,]),
        list_normalised[["B cells"]]["B cell  plasma cell",],
        colSums(list_normalised[["B cells"]][ w2,]),
        colSums(list_normalised[["B cells"]][ w3,])
    )
  colnames(specific_cell_frequencies) = cell_types_of_interest
  ### multiply out for overall frequencies
  totals = list_normalised[["broad all"]][cell_groups_of_interest,]
  specific_cell_frequencies_norm = specific_cell_frequencies
  for(i in c(1:length(cell_types_of_interest))){
    specific_cell_frequencies_norm[,i] = specific_cell_frequencies[,i]*totals[i,]/100
  }
  plot(specific_cell_frequencies_norm, t(totals))
  a = specific_cell_frequencies_norm/t(totals)
  apply(a, 2, max)
  
  
  ### plot correlations
  for(i1 in c(1:length(list_normalised))){
    mat = list_normalised[[i1]]
    mat1 = rbind(mat, t(specific_cell_frequencies_norm))
    name = names(list_normalised)[i1]
    if(name %in% c( "broad all" , "CAF",  "Non-immune" )){
      use = colnames(mat1)[grep("_biopsy", colnames(mat1), invert = T)]
      mat1 = mat1[,use]
    }
    a = apply(mat1, 1, function(x){length(which(x!=0))})
    mat1 = mat1[which(a>=1), ]
    print(length(mat1[,1]))
    cortest = corr.test(t(mat1),adjust="holm")
    pval = cortest $ p
    rval = cortest $ r
    pval = pval[colnames(specific_cell_frequencies_norm), rownames(mat)]
    rval = rval[colnames(specific_cell_frequencies_norm), rownames(mat)]
    
    
    fileout1=concat(c(output_directory,"Correlation_cell_types_pairs_plot_biopsy_overall_", batch,"_",name,"_2.pdf"))
    w=3.5
    if(length(mat[,1])>15){w=8}
    if(length(mat[,1])<7){w=1.7}
    
    pdf(file=fileout1, height=w*1.1, width=w*1.1)
    par(mfrow= c(1,1), mar = c(5,5,3,3))
    corrplot(rval, type="full", p.mat=pval, insig="label_sig",  sig.level=0.05, title ="", method = "ellipse", diag = T, order = "original", cl.pos = 'n')
    dev.off()
    
  }
  
  
  
    
  ### broad immune
    
  mat = rbind(list_normalised[[ "broad immune"]], (list_normalised[["B cells"]]["B cell  plasma cell",]+list_normalised[["B cells"]]["B cell  plasmablast",]))
  rownames(mat)[length(rownames(mat))] = "PB/PC of B cells"
  cell_types = rownames(mat)
  
  fileout1=concat(c(output_directory,"Correlation_broad_immune_frequencies", batch,".pdf"))
  w=2.3
  pdf(file=fileout1, height=w*1, width=w*1*7)
  par(mfrow= c(1,7), mar = c(5,5,5,5))
  
  for(i1 in c(1:length(cell_types))){
    for(i2 in c(i1:length(cell_types))){
      if(i1<i2){
        x = mat[cell_types[i1],]
        y = mat[cell_types[i2],]
        factors = rep(1,length(x))
        names(factors) =names(x)
        factors[which(names(factors) %in% unlist(groups_PCA2))] = 2
        factors[which(names(factors) %in% unlist(groups_PCA3))] = 3
        
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
        
        cortest = corr.test(cbind(x,y),adjust="holm")
        pval = cortest $ p[1,2]
        rval = cortest $ r[1,2]
        
        
        name = concat(c(cell_types[i1]," vs ",cell_types[i2]))
        range1 = range
        if(max(x)<0.75*max(y)){
          range1 = range(x)
        }
        cex = 1
        plot(range1, range, pch=20, col="white",main = concat(c(name,"\n\np-value:",signif(p_value,digits = 3),"\nr-value:",signif(rval,digits = 3))), xlab = concat(c("% ",cell_types[i1])),ylab = concat(c("% ",cell_types[i2])),cex=cex, cex.lab=cex+0.1,cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), axes = FALSE)
        
        library(yarrr)
        library(RColorBrewer)
        cols1 =  add.alpha (piratepal(palette = "google"), alpha = 0.95)
        cols =  add.alpha (cols1, alpha = 0.5)
        pches = c(21:25)
        Fun<-function(x){x}
        segments(min,Fun(scale),max(max)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
        segments(Fun(scale),min,Fun(scale),max(max)+0.5,col = "grey",lwd = 1,lty = 3 )
        mtext(side = 2, text = scale, line = 0.45,cex= cex-0.2,  at =Fun(scale), las = 2, font = 1)
        scale1 = scale#[which(scale!=0)]
        # scale1 = scale1[which(scale1<=max(x)*1.2)]
        mtext(side = 1, text = scale1, line = 0.45,cex= cex-0.2,  at =Fun(scale1), las = 2, font = 1)
        
        new = seq(from= min(x11), to = max(x11), length = 100)
        predicted.intervals <- predict(fit1, data.frame(x11= new),interval='prediction',level = 0.85)
        cols2 =  add.alpha ("grey", alpha = 0.85)
        points(new,predicted.intervals[,1],col= cols2,lwd=4, type = "l")
        polygon(x = c(new, rev(new)),y=c(predicted.intervals[,2],rev(predicted.intervals[,3]) ),  col= add.alpha(cols2,alpha = 0.25),border = NA)
        points(x,y, pch = pches[factors], bg=cols[factors], cex = 1.5, col= add.alpha("white", alpha = 0))	
        points(x,y, pch = pches[factors], bg=cols[factors], cex = 1.5, col= add.alpha("white", alpha = 0))	
        
      }
    }
  }
  plot(range1, range, col = "white",bg = "white", main = "", xlab = "", ylab = "", cex = 1,lwd = 2, axes = F)
  leg = c("PDAC150K", "PENG","STEELE")
  legend("topleft", leg, pch = pches,cex= 0.8, bty="n", pt.bg = cols1, col = cols1, pt.lwd = 2, text.font = 2)
  
  
  dev.off()
  
  ##### get p-values per cohort
  
  cell_type1 = NULL
  cell_type2 = NULL
  p_val1 = NULL
  p_val2 = NULL
  p_val3 = NULL
  
  for(i1 in c(1:length(cell_types))){
    for(i2 in c(i1:length(cell_types))){
      if(i1<i2){
        x = mat[cell_types[i1],]
        y = mat[cell_types[i2],]
        factors = rep(1,length(x))
        names(factors) =names(x)
        factors[which(names(factors) %in% unlist(groups_PCA2))] = 2
        factors[which(names(factors) %in% unlist(groups_PCA3))] = 3
        w1 = which(factors==1)
        w2 = which(factors==2)
        w3 = which(factors==3)
        
        x11 = x[w1]
        y11 = y[w1]
        fit1 <- lm( y11~poly(x11,2))
        fit1 <- lm( y11~x11)
        p_value = summary.lm(fit1)$ coefficients["x11",4]
        p_val1 = c(p_val1, p_value)
        
        x11 = x[w2]
        y11 = y[w2]
        fit1 <- lm( y11~poly(x11,2))
        fit1 <- lm( y11~x11)
        p_value = summary.lm(fit1)$ coefficients["x11",4]
        p_val2 = c(p_val2, p_value)
        
        x11 = x[w3]
        y11 = y[w3]
        fit1 <- lm( y11~poly(x11,2))
        fit1 <- lm( y11~x11)
        p_value = summary.lm(fit1)$ coefficients["x11",4]
        p_val3 = c(p_val3, p_value)
        cell_type1 = c(cell_type1, cell_types[i1])
        cell_type2 = c(cell_type2, cell_types[i2])
      }
    }}
  x = cbind(cell_type1,cell_type2,p_val1 , p_val2 ,p_val3)
  colnames(x) = c("cell type 1","cell type 2","p-val (PancrImmune)", "p-val (PENG)", "p-val (STEELE)")
     
  out_file_table = concat(c(output_directory, "Correlations_per_dataset_COMBINED_PENG_STEELE_PDAC150K_all.txt"))
  write.table(x, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  
  
}
Plot_correlations_between_all(input_directory, output_directory)
############################################ plot PCA of CD45 %s
Plot_PCA<-function(input_directory, output_directory){
  
  groups_PCA1 = readRDS(file = concat(c(input_directory, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))
  groups_PCA2 = readRDS(file = concat(c(input_directory, "Immunogroups_2019_PUBLISHED_PDAC.rds")))
  groups_PCA3 = readRDS(file = concat(c(input_directory,"Immunogroups_PDAC_GSE155698.rds")))
  
  all_samples_plot = sort(unique(c(unlist(groups_PCA1), unlist(groups_PCA2), unlist(groups_PCA3))))
  w_include = all_samples_plot
  file = concat(c(input_directory,"Seurat_cell_counts_COMBINED_PENG_STEELE_PDAC150K_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  sample = colnames(p)[-1]
  sample = sample[-1]
  sites = sample
  disease = sample
  sites[grep("PDAC_TISSUE", sites)] = "PDAC tissue"
  sites[grep("AdjNorm_TISSUE", sites)] = "AdjNorm tissue"
  sites[grep("PDAC_PBMC", sites)] = "PDAC PBMC"
  sites[grep("Healthy_PBMC", sites)] = "Healthy PBMC"
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
      w = intersect(w, which(p[,1] %in% c( "B cell"  , "Myeloid","NK" , "T cell"   )))
    }
    cell_type_sub = p[w,1]
    p1 = p[w, w_include]
    m_count_norm = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(w_include), dimnames = c(list(cell_type_sub), list(w_include)))
    
    for(s in c(1:length(w_include))){
      m_count_norm[, w_include[s]]= as.numeric(p1[, w_include[s]])*100/sum(as.numeric(p1[, w_include[s]]))
      if(sum(as.numeric(p1[, w_include[s]]))<7){m_count_norm[, w_include[s]] = -1}
    }
    if(length(which(apply(m_count_norm, 2, max)!=-1))>3){include = c(include, ind)}
    m_count_norm= m_count_norm[order(rownames(m_count_norm)),]
    #m_count_norm  = m_count_norm [,rownames(group_mat)]
    list_normalised= c(list_normalised, list(m_count_norm))
  }
  names(list_normalised) = cell_groups
  list_normalised = list_normalised[include]
  
  
  for(i in c(1:length(list_normalised))){
    print (names(list_normalised)[i])
    print (rownames(list_normalised[[i]]))
  }
  
  
  all_proportions = t(list_normalised[[ "broad immune" ]])

  x <- princomp(all_proportions,dim=2)
  PCA_coordinates_sub = x$scores
  
  x1 = x$scores[,1]
  x2 = x$scores[,2]
  
  y1 = x$loadings[,1]
  y2 = x$loadings[,2]
  
  factor1 = rep(0,length(x1))
  names(factor1) = names(x1)
  factor2 = factor1
  factor1[unlist(groups_PCA1)] = "PDAC150K"
  factor1[unlist(groups_PCA2)] = "PENG"
  factor1[unlist(groups_PCA3)] = "STEELE"
  factor2[c(groups_PCA1[[1]],groups_PCA2[[1]],groups_PCA3[[1]] )] = "ME"
  factor2[c(groups_PCA1[[2]],groups_PCA2[[2]],groups_PCA3[[2]] )] = "AE"
  
  
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
  analysis_all = "PDAC150K_PENG_STEELE"
  fileout1=concat(c(output_directory,"", analysis_all,"_biopsy_PCA_scores_", batch,".pdf"))
  w=2.2
  pdf(file=fileout1, height=w*1, width=w*3)
  par(mfrow= c(1,3), mar = c(5,5,3,3))
  
  pches = c(21:25)
  plot(x1,x2,pch = pches[matching1], col = NA,bg = cols2[matching2], main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1.3,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  
  xy = cbind(x1,x2)
  library(car)
  for(pca in c(1:length(factors2))){
    w = which(factor2 ==factors2[pca])
    if(length(w)>1){
      x = xy[w,1]
      y = xy[w,2]
      dataEllipse(x, y, levels=c(0.65), add = TRUE, col = cols2[pca], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
    }
  }
  
  plot(x1,x2,pch = pches[matching1], col = "white",bg = "white", main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  legend("topleft", paste("",factors1), pch = pches,cex= 0.8, bty="n", pt.bg = "black", col = "black", pt.lwd = 2, text.font = 2)
  
  plot(x1,x2,pch = pches[matching1], col = cols[matching1],bg = cols2[matching1], main = "all cells", xlab = "PCA1", ylab = "PCA2", cex = 1,lwd = 2, xlim = range(x1*1.1), ylim = range(x2*1.1))
  text(x1, y = x2, labels = factor2,cex = 0.4,font = 1)
  
  dev.off()
  
}
Plot_PCA(input_directory, output_directory)

############################################ Make summary plot for Figure 1
Summary_plot<-function(input_directory, output_directory){
  files = c("Stats_PCA_cell_types_by_groups_PDAC150Ka.txt","Stats_PCA_cell_types_by_groups_PENG.txt","Stats_PCA_cell_types_by_groups_STEELE.txt")
  names(files) = c("PancrImmune","PENG","STEELE")
  all_stats = NULL
  all_cell_types = NULL
  for(f in c(1:length(files))){
    p <- as.matrix(read.csv(concat(c(input_directory,files[f])), head=F, sep="\t"))
    headers = c("cell_type",p[1,c(1:(length(p[1,])-1))])
    p=p[-1,]
    colnames(p) = headers
    w = which(p[,"p.analysis"] %in% c("cell_proportions_CD45:B cell","cell_proportions_CD45:T cell","cell_proportions_CD45:Myeloid",
      "cell_proportions_B cells","cell_proportions_Myeloid","cell_proportions_T cells"))
    p=p[w,]
    p=p[which(p[,"p.site"]=="biopsy"),]
    all_stats= c(all_stats, list(p))
    all_cell_types = sort(unique(c(all_cell_types, p[,1])))
  }
  
  mat_pval = matrix(data =0, nrow = length(all_cell_types), ncol = length(files), dimnames = c(list(all_cell_types), list(names(files))))
  mat_dir = matrix(data =0, nrow = length(all_cell_types), ncol = length(files), dimnames = c(list(all_cell_types), list(names(files))))
  for(f in c(1:length(files))){
    mat = all_stats[[f]]
    table(table(mat[,1]))
    dir = (as.numeric(mat[,"mean.group. 1"])-as.numeric(mat[,"mean.group. 2"]))/(as.numeric(mat[,"mean.group. 1"])+as.numeric(mat[,"mean.group. 2"]))
    mat_dir[mat[,1],f] = dir
    table(mat_dir[mat[,1],f] == dir)
    mat_pval[as.character(mat[,1]),names(files)[f]] = as.numeric(as.character(mat[,"p_value"]))
    print(table(mat_pval[mat[,1],names(files)[f]] == as.numeric(as.character(mat[,"p_value"]))))
  }
  
  a = apply(mat_pval, 1, function(x){length(which(x<0.05))})
  b = apply(mat_dir, 1, function(x){length(which(x!=0))})
  c = apply(mat_dir, 1, function(x){length(which(x<0))})
  d = apply(mat_dir, 1, function(x){length(which(x>0))})
  
  pairs(mat_dir[which(b==3),])
  
  mat_dir_plot = mat_dir[which(b==3),]
  
  w1 = intersect(which(a>=1), which(c>=3))
  w2 = intersect(which(a>=1), which(d>=3))
  w = sort(c(w1,w2))
  
  
  #### plot heatmap (col = dir, * = signif)
  mat1 = mat_dir[w,]
  mat2 = mat_pval[w,]
  ylab1 = rownames (mat1)
  ylab = gsub("  pl"," pl", ylab1, fixed = T)
  ylab = gsub("_"," ", ylab, fixed = T)
  xlab1 = colnames (mat1)
  xlab = xlab1
  
  posx = c(1:length(xlab))
  posy = c(1:length(ylab))
  #posy[grep("T cell", ylab)] = posy[grep("T cell", ylab)]+0.5
  
  max_pos1 = max(posx)+0.5
  max_pos2 = max(posy)-0.
  
  fileout1=concat(c(output_directory,"Comparison_between_datasets_cell_counts_significance_heatmap.pdf"))
  w=4.7
  pdf(file=fileout1, height=w*1., width=w*1.3)
  cex = 0.9
  par(mfrow= c(1,2), mar = c(1.5,12,9,1))
  
  plot(c(1, max_pos1), c(-1, max_pos2), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "", 
       ylim = c(-1, max(max_pos2)),xlim = c(0.5, max(max_pos1)), axes = F)
  mtext(side = 2, text = ylab, line = 0.35,cex=0.7,  at = posy, las = 1, font = 1)
  mtext(side = 3, text = xlab, line = 0.35,cex=0.7,  at = posx, las = 3, font = 1)
  
  library(RColorBrewer)
  cols1=rev(add.alpha(colorRampPalette(c("darkorange",'orange','yellow',"white", "lightgrey","grey","black"))(201), alpha = 0.5))
  cols2 =  add.alpha (cols1, alpha = c(0.5))
  mat = ((mat1+1)*100)+1
  
  x = 0
  for(i in c(1:length(posx))){
    for(j in c(1:length(posy))){
      rect(posx[i]-0.5+x, posy[j]-0.5, posx[i]+0.5+x, posy[j]+0.5,col =cols1[mat[j,i]], border = "grey")
      if(mat2[j,i]<0.05){
        text(posx[i], posy[j], "*")}
    }}
  
  plot(c(1, max_pos1), c(-1, max_pos2), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "", 
       ylim = c(-1, (max(max_pos2)-1)),xlim = c(0.5, max(max_pos1)), axes = F)
  seq = seq(from = -1, to = 1, length = 5)
  seq1= seq(from = -1, to = 1, length = 200)
  mat = ((seq1+1)*100)+1
  posx1 = seq(from = 1.1, to = max_pos2*0.5, length = 200)
  posx2 = seq(from = 1.1, to = max_pos2*0.5, length = 5)
  wid = posx1[2]-posx1[1]
  for(j in c(1:length(mat))){
    rect(0.5,posx1[j]-wid, 1.5, posx1[j]+wid, col =cols1[mat[j]], border = NA)
  }
  for(j in c(1:length(seq))){
    text(2.6, posx2[j], seq[j], cex = 0.85)
  }
  
  dev.off()
  
  
  rangex = range(mat_dir_plot)
    
  fileout1=concat(c(output_directory,"Comparison_between_datasets_cell_counts_significance_correlation.pdf"))
  w=2.5
  pdf(file=fileout1, height=w*1, width=w*1*3.5)
  cex = 0.9
  par(mfrow= c(1,4), mar = c(4,4,4,4))
  
  for(f1 in c(1:length(files))){
    for(f2 in c(1:length(files))){
      if(f1<f2){
        x = mat_dir_plot[,f1]
        y = mat_dir_plot[,f2]
        fit1 <- lm( y~x)
        print(summary.lm(fit1)$coefficients)
        
        main= concat(c("p-value:",signif(summary.lm(fit1)$coefficients[1,4], digits = 3),
                       "\nr-value:",signif(summary.lm(fit1)$r.squared, digits = 3)))
        plot(rangex, rangex, pch=20, col="white",xlab=concat(c("cell enrichment (",names(files)[f1],")")),ylab =concat(c("cell enrichment (",names(files)[f2],")")),cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "black",tck=0, mgp = c(2,0,0), main = main, axes = T)
        segments(0,-1,0,1,col = "grey",lwd = 1.5, lty = 3)
        segments(1,0,-1,0,col = "grey",lwd = 1.5, lty = 3)
        typ = rep(1, length(x))
        typ[grep("T cell", rownames(mat_dir_plot))] = 2
        typ[grep("Myeloid", rownames(mat_dir_plot))] = 3
        pches = c(21:25)
        col = add.alpha("black",alpha = 0.5)
        
        fit1 <- lm( y~x)
        print(summary.lm(fit1)$coefficients)
        new = seq(from= min(x), to = max(x), length = 100)
        predicted.intervals <- predict(fit1, data.frame(x= new),interval='prediction',level = 0.85)
        cols2 =  add.alpha ("grey", alpha = 0.85)
        points(new,predicted.intervals[,1],col= cols2,lwd=4, type = "l")
        polygon(x = c(new, rev(new)),y=c(predicted.intervals[,2],rev(predicted.intervals[,3]) ),  col= add.alpha(cols2,alpha = 0.25),border = NA)
        points (x,y,col = col,bg = col, pch = pches[typ], cex = 0.8)
      }
    }
  }
  
  leg = c("B cell","T cell","Myeloid cell")
  plot(c(1, max_pos1), c(-1, max_pos2), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "", 
       ylim = c(-1, (max(max_pos2)-1)),xlim = c(0.5, max(max_pos1)), axes = F)
  legend("topleft", leg, pch = pches,cex= 0.8, bty="n", pt.bg = col, col = col, pt.lwd = 2, text.font = 2)
  
  dev.off()
  
  
  
  
}
Summary_plot(input_directory, output_directory)
  
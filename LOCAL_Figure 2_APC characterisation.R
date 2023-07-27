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
output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/APCs/"
input_directory_groups  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"
groups_PCA = readRDS(file = concat(c(input_directory_groups, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))

##### get meta data (cell type) with VDJ information
Get_BCR_information<-function(input_directory_groups){
  type = "BCR"
  type1 = "B_cells"
  file = concat(c(input_directory_groups,"VDJ_information_",type,"_PDAC150Ka.txt"))
  VDJ <- as.matrix(read.csv(file, head=T, sep="\t"))
  clone = VDJ[,"clone1"]
  
  file = concat(c(input_directory_groups,"Cell_annotation_ALL_PDAC150Ka.txt"))
  ### this file can be replaced with the metadata information from your Seurat object
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  p1 <- p1[grep("B cell", p1[,"cell_type"]), ]
  id1 = p1[,"barcode"]
  sample = p1[,"sample"]
  id = apply(cbind(id1, sample), 1, function(x){
    x1 = gsub(concat(c(x[2], "_")), "", x[1])
    return(concat(c(x1,"||", x[2])))})
  id = gsub("-1","", id, fixed = T)
  length(intersect(id, names(clone)))
  cell_type = p1[,"cell_type"]
  Patient = p1[,"sample1"]
  Sample.Type= p1[,"sample"]
  Patient = gsub("_CD45p1","", Patient)
  Patient = gsub("_CD45p2","", Patient)
  Sample.Type[grep("biopsy", Sample.Type)] = "biopsy"
  Sample.Type[grep("blood", Sample.Type)] = "blood"
  names(cell_type) = id
  sample = apply(cbind(Patient, Sample.Type), 1, paste, collapse = "-")
  table(sample)
  samples = sort(unique(sample))
  cell_type1 = cell_type
  cell_type[grep("memory", cell_type)] = "B cell memory"
  cell_type[grep("activated", cell_type)] = "B cell activated"
  cell_type[grep("plasma", cell_type)] = "PB/PC"
  cell_type[grep("GC", cell_type)] = "B cell GC"
  
  names(cell_type) = id
  names(cell_type1) = id
  names(sample) = id
  
  list_return = c(list(VDJ), list(cell_type), list(cell_type1), list(sample), list(Patient), list(Sample.Type))
  names(list_return) = c("VDJ_object", "cell_type","cell_type_broad","sample","Patient","Sample.Type")
  return(list_return)
}

VDJ_list_BCR = Get_BCR_information()


APC_analysis_PDAC150K<-function(output_directory,input_directory, batch, groups_PCA, VDJ_list_BCR){
  type = "All"
  type1 = "APCs"
  file = concat(c(input_directory,"Pathway_levels_per_cell.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  cell_type_level0 = p[,"cell_type"]
  cell_type_level1 = p[,"cells_type_level1"]
  cell_type_level2 = p[,"cells_type_level2"]
  sample = p[,"sample"]
  samples = sort(unique(sample))
  pathways = c("ANTIGEN.PROCESSING.AND.PRESENTATION.MHC.CLASS.I1",
                    "ANTIGEN.PROCESSING.AND.PRESENTATION..MCHII1" ,    
                    "APC.MCHII_inhibitory1"                        ,   
                    "Activation1"                                   ,  
                    "Apoptosis.1"                                    , 
                    "COX.IS1"                                         ,
                    "Cytokine1"                                       ,
                    "Inhibitory1"                                     ,
                    "NFKB.PATHWAY1")
  
  mat_pws = NULL
  for(i in c(1:length(pathways))){
    if(length(mat_pws)==0){
      mat_pws = as.numeric(p[,pathways[i]])
    }else{mat_pws=cbind(mat_pws , as.numeric(p[,pathways[i]]))}
  }
  colnames(mat_pws) = pathways
  rownames(mat_pws) = rownames(p)
  
  genes = colnames(p)[grep("expr.", colnames(p))]
  mat_genes = NULL
  for(i in c(1:length(genes))){
    if(length(mat_genes)==0){
      mat_genes = as.numeric(p[,genes[i]])
    }else{mat_genes=cbind(mat_genes , as.numeric(p[,genes[i]]))}
  }
  colnames(mat_genes) = genes
  rownames(mat_genes) = rownames(p)
  
  cell_type = cell_type_level0
  cell_types = sort(unique(cell_type))
  
  ## difference between groups
  groups_PCA1 = groups_PCA
  groups_PCA1[[1]] = gsub("_biopsy","-blood",groups_PCA1[[1]],fixed = T)
  groups_PCA1[[2]] = gsub("_biopsy","-blood",groups_PCA1[[2]],fixed = T)
  groups_PCA[[1]] = gsub("_","-",groups_PCA[[1]],fixed = T)
  groups_PCA[[2]] = gsub("_","-",groups_PCA[[2]],fixed = T)
  
  groups_PCA_all = c(list(groups_PCA),list(groups_PCA1))
  names(groups_PCA_all) = c("tumour", "blood")
  
  Get_pAPCs<-function(output_directory,input_directory, batch, groups_PCA, mat_pws, pathways){
    g = 2
    x = mat_pws[,pathways[g]]
    
    cells_sub = c(grep("T cell CD8", cell_type_level2), grep("DC", cell_type_level2))
    y = x[cells_sub]
    x1 = data.frame(y) ### subset of data
    class = cell_type_level2[cells_sub]
    class[which(class != "T cell CD8")] = "APC"
    class[which(class == "T cell CD8")] = "non-APC"
    x1$class = factor(class)
    
    model <- glm(class ~ ., x1,family=binomial(link='logit'))
    p <- predict(model, newdata=x1, type="response")
    glm.pred=rep("APC", length(p))
    glm.pred[p>0.5]="non-APC"
    length(which(glm.pred==class))*100/length(p)
    table(glm.pred,class)
    ### get decision point (it should happen between the two means of the distributions)
    a1 = mean(x[grep("T cell CD8", cell_type_level2)])
    a2 = mean(c(x[grep("B cell", cell_type_level2)], x[grep("DC", cell_type_level2)]))
    
    x_trial = data.frame(seq(from = a1, to = a2, length = 10000))
    x_trial = data.frame(x_trial)
    colnames(x_trial) = "y"
    p <- predict(model, newdata=x_trial, type="response")
    threshold_pathway = x_trial[,1][max(which(p>0.5))]
    
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",pathways[g],".pdf"))
    w=2.5
    pdf(file=fileout1, height=w*1.1, width=w*1.15)
    par(mfrow= c(1,1), mar = c(5,5,2,3))
    hist(x, breaks = 100, border = NA, col = "grey", freq = F, ylim = c(0,1.8), xlab = "score", xlim = c(-1,3))
    h1 = hist(x[grep("T cell CD8", cell_type_level2)], breaks = 50, plot = F)
    points(h1$mids,h1$density , col = add.alpha("red", 0.6),  type = "l", lwd = 4)
    
    h1 = hist(x[grep("DC", cell_type_level2)], breaks = 50, plot = F)
    points(h1$mids,h1$density , col = add.alpha("lightskyblue", 0.6),  type = "l", lwd = 4)
    
    segments(threshold_pathway, -10, threshold_pathway, 1000, col = "black", lwd= 2, lty = 2)
    dev.off()
    
    library("mixtools")
    library(tibble)
    library(ggplot2)
    
    plot_mix_comps <- function(x, mu, sigma, lam) {
      lam * dnorm(x, mu, sigma)
    }
    
    mixmdl <- normalmixEM(x, k = 2)
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",pathways[g],"1.pdf"))
    w=3.2
    pdf(file=fileout1, height=w*1, width=w*1.9)
    par(mfrow= c(1,1), mar = c(5,5,2,3))
    
    data.frame(x = mixmdl$x) %>%
      ggplot() +
      geom_histogram(aes(x, ..density..), binwidth = 0.05, colour = "black", 
                     fill = "white") +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                    colour = "red", lwd = 1.5) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                    colour = "blue", lwd = 1.5) +
      ylab("Density") + 
      xlab(concat(c(pathways[g], " score")))+
      geom_vline(xintercept= c(a <- c(threshold_pathway)), 
                 colour = c("red")[c(1)])
    dev.off()
    
    ### proportions of APCs per cell type per sample in the tumour
    y = rep(0, length(x))
    y[which(x>threshold_pathway)] = 1
    t = table(cell_type_level2,y)
    pAPC_status = y
    percentage_pos = t[,2]*100/rowSums(t)
    type_names = names(percentage_pos)
    col = rep(1, length(type_names))
    col[grep("DC", type_names)] =2
    col[grep("MoMac", type_names)] =3
    col[grep("Mono", type_names)] =4
    col[grep("granulocyte", type_names)] =5
    col[grep("NK", type_names)] =6
    col[grep("T cell", type_names)] =7
    col[grep("Mast", type_names)] =8
    col[grep("ILC", type_names)] =8
    col[grep("iNKT", type_names)] =8
    col[grep("Megakaryocyte", type_names)] =8
    names(col)=type_names
    names(percentage_pos) = gsub("B cell antigen-experienced", "B cell Ag.-exp.",names(percentage_pos))
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",pathways[g],"_proportions_overall.pdf"))
    w=3.2
    pdf(file=fileout1, height=w*1.1, width=w*1.6)
    par(mfrow= c(1,1), mar = c(12,5,2,3))
    barplot(percentage_pos, las = 2, col = col, ylab = "% positive", border = NA)
    dev.off()
    
    ##### get broad proportion of APCs per sample
    cell_type = cell_type_level2
    cell_types= sort(unique(cell_type))
    mat_prop_of_APC = matrix(data = 0, nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
    for(s in c(1:length(samples))){
      w = intersect(which(sample==samples[s]), which(y==1))
      if(length(w)>=10){
        t = table(cell_type[w])
        t = t*100/sum(t)
        mat_prop_of_APC[samples[s],names(t)] = t
      }
    }
    
   
    
    a = cbind(apply(mat_prop_of_APC[groups_PCA[[1]],], 2, mean), apply(mat_prop_of_APC[groups_PCA[[2]],], 2, mean))
    b = colSums(a[grep("DC", rownames(a)),])
    c = colSums(a[grep("B cell", rownames(a)),])
    d =rbind(b, c, a["MoMac",])
    rownames(d) = c("DC", "B cell", "MoMac")
    colnames(d) = c("Mean cellular contribution to pAPCs (ME)","Mean cellular contribution to pAPCs (AE)")
    
    e = rbind(d, a)
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",batch,"_",pathways[g],"_proportions_between_groups.txt"))
    write.table(e, file = fileout1, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",pathways[g],"_proportions_between_groups.pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*2.5, width=w*1*1.7)
    par(mfrow= c(2,1), mar = c(12,5,3,3))
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(ind in c(1:length(groups_PCA_all))){
      groups_PCA_use = groups_PCA_all[[ind]]
      
      group_PCA_list = NULL
      factor_PCA = NULL
      for(i in c(1:length(groups_PCA_use))){
        group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA_use[[i]]))
        factor_PCA = c(factor_PCA, rep(i, length(groups_PCA_use[[i]])))
      }
      name = names(groups_PCA_all)[ind]
      mat = mat_prop_of_APC[group_PCA_list,]
      factor = factor(factor_PCA)
      fit = manova(formula = mat ~ factor)
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
        if(length(mean)==0){means = Means_factor(factor, mat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      names(p_value) = nam
      sort(p_value)
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.site = rep(name, length(nam))
      p.analysis = rep(pathways[g], length(nam))
      x = cbind(p.site, p.analysis, combined_p_value)
      summary_tables = x
      
      groups = NULL
      p_values = combined_p_value[,1]
      for(i in c(1:length(mat[1,]))){
        g1 = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(groups_PCA_use))){
          x = mat [gsub("_","-",groups_PCA_use[[g]]), i]
          x = x[which(x!=-1)]
          data = c(data,x )
          factor = c(factor, rep(g, length(x)))
          g1 = c(g1, list(x))}
        groups = c(groups, list(g1))
        # group = factor(factor)
        # fit = aov(formula = data ~ group)
        # pvalues = summary.lm(fit)$ coefficients[,4]
        # pVal <- anova(fit)$'Pr(>F)'[1]
        # p_values = c(p_values, pVal)
      }
      factors1 = paste("group",c(1:length(groups_PCA_use)))
      factors = gsub("_"," ",colnames(mat))
      main = concat(c(pathways[g], " ",name))
      max = max(c(unlist(groups), unlist(groups))*1.35)
      min = 0
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 20
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
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = gsub("B cell antigen-experienced", "B cell Ag.-exp.",factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
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
          points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg=cols1[i1], cex = 0.7)
        }}
      
      for(i in c(1:l)){	
        b = max*0.035
        signif_threshold = 0.05
        if(p_values[i]<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))
          y = y+1*b
          # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
          text(i, y+2*b, labels = pval1, cex = 1.7)
        }
      }
      
    }
    dev.off()
    names(pAPC_status) = names(cell_type)
    
    return(pAPC_status)
  }
  pAPC_status = Get_pAPCs(output_directory,input_directory, batch, groups_PCA, mat_pws, pathways)
  
  Plot_isotype_usage_APC<-function(pAPC_status, output_directory,input_directory, batch, groups_PCA, VDJ_list_BCR){
    VDJ_list = VDJ_list_BCR
    type = "BCR"
    type1 = "B_cells"
    
    VDJ = VDJ_list[[1]]
    cell_type_broad = VDJ_list[[2]]
    cell_type = VDJ_list[[3]]
    sample = VDJ_list[[4]]
    Patient = VDJ_list[[5]]
    Sample.Type = VDJ_list[[5]]
    
    samples = sort(unique(sample))
    Patients = sort(unique(Patient))
    Sample.Types = sort(unique(Sample.Type))
    cell_types = sort(unique(cell_type))
    cell_types_broad = sort(unique(cell_type_broad))
    clone = VDJ[,"clone1"]
    w_clone = names(which(clone!='-'))
    
    isotype = VDJ[,"constant_region1"]
    isotypes = sort(unique(isotype))
    
    
    cell_type_use = cell_types[grep("B cell",cell_types)]
    p_values = NULL
    for(c in c(1:length(cell_type_use))){
      w = names(which(cell_type==cell_type_use[c]))
      w = intersect(w, names(isotype))
      t = table(isotype[w],pAPC_status[w])
      t = t[which(rownames(t)!="-"),]
      for(i in c(1:length(t[1,]))){
        t[,i] = t[,i]*100/sum(t[,i])
      }
      p_val = signif(fisher.test(t, alternative = "greater",simulate.p.value = T,B = 1000000)$p.value, digits = 4)
      p_values = c(p_values, p_val)
    }
    names(p_values) = cell_type_use
    print("P-values of isotype usage by APC status: ")
    print(cbind(p_values))
  }
  Plot_isotype_usage_APC(pAPC_status, output_directory,input_directory, batch, groups_PCA, VDJ_list_BCR)
  ## no significant differences in isotype usages
  
  Plot_pathways_per_APC_population<-function(pAPC_status, output_directory,input_directory, batch, groups_PCA, mat_pws, pathways){
    ### only include cell types that are observed enough times in the pAPC pool
    w =which(pAPC_status==1)
    t = table(cell_type[w], sample[w])
    t = t[,grep("biopsy", colnames(t))]
    a1 = apply(t[,groups_PCA[[1]]], 1, function(x){length(which(x>=10))})
    a2 = apply(t[,groups_PCA[[2]]], 1, function(x){length(which(x>=10))})
    cell_types_use = names(which(apply(cbind(a1,a2), 1, min)>=3))
    
    samples_use = samples[grep("biopsy", samples)]
    
    cell_type = cell_type_level2
    cell_types = sort(unique(cell_type))
    
    ### mean pathway per patient per cell type
    list_pw = NULL
    for(g in c(1:length(pathways))){
      mat_pw_mean = matrix(data = 0, nrow = length(samples_use),ncol = length(cell_types_use), dimnames = c(list(samples_use), list(cell_types_use)))
      for(c in c(1:length(cell_types_use))){
        print(c)
        w = names(which(cell_type==cell_types_use[c]))
        for(s in c(1:length(samples_use))){
          w2 = intersect(w, names(which(sample==samples_use[s])))
          mat_pw_mean[samples_use[s],cell_types_use[c]] = median(mat_pws[w2,pathways[g]])
         }
      }
      list_pw = c(list_pw, list(mat_pw_mean))
    }
    names(list_pw) = pathways
    
    fileout1= concat(c(output_directory,"Pathway_levels_pAPC_pathway_levels_between_groups.pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*4, width=w*1*3)
    par(mfrow= c(4,3), mar = c(12,5,3,3))
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(ind in c(1:length(list_pw))){
      groups_PCA_use = groups_PCA_all[[1]]
      
      group_PCA_list = NULL
      factor_PCA = NULL
      for(i in c(1:length(groups_PCA_use))){
        group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA_use[[i]]))
        factor_PCA = c(factor_PCA, rep(i, length(groups_PCA_use[[i]])))
      }
      name = names(list_pw)[ind]
      mat = list_pw[[ind]]
      mat=mat[group_PCA_list,]
      factor = factor(factor_PCA)
      fit = manova(formula = mat ~ factor)
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
        if(length(mean)==0){means = Means_factor(factor, mat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      names(p_value) = nam
      sort(p_value)
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.site = rep(name, length(nam))
      p.analysis = rep(pathways[g], length(nam))
      x = cbind(p.site, p.analysis, combined_p_value)
      summary_tables = x
      
      groups = NULL
      p_values = combined_p_value[,1]
      for(i in c(1:length(mat[1,]))){
        g1 = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(groups_PCA_use))){
          x = mat [gsub("_","-",groups_PCA_use[[g]]), i]
          x = x[which(x!=-1)]
          data = c(data,x )
          factor = c(factor, rep(g, length(x)))
          g1 = c(g1, list(x))}
        groups = c(groups, list(g1))
        # group = factor(factor)
        # fit = aov(formula = data ~ group)
        # pvalues = summary.lm(fit)$ coefficients[,4]
        # pVal <- anova(fit)$'Pr(>F)'[1]
        # p_values = c(p_values, pVal)
      }
      factors1 = paste("group",c(1:length(groups_PCA_use)))
      factors = gsub("_"," ",colnames(mat))
      main = concat(c(name, "\n","pAPC only"))
      max = max(c(unlist(groups), unlist(groups))*1.35)
      min = min(c(unlist(groups), unlist(groups)))
      b = (max-min)*0.034
      min = min-b
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 10
      max_scale = max
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
      if(range <0.001){scale = c(-100:100)*0.0001}
      if(range <0.0001){scale = c(-100:100)*0.00001}
      cex = 0.9
      Fun<-function(x){x}
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = gsub("B cell antigen-experienced", "B cell Ag.-exp.",factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
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
          points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg=cols1[i1], cex = 0.7)
        }}
      
      for(i in c(1:l)){	
        b = max*0.035
        signif_threshold = 0.05
        if(p_values[i]<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))
          y = y+1*b
          # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
          text(i, y+2*b, labels = pval1, cex = 1.7)
        }
      }
      
    }
    dev.off()
    
    
    
    
  }
  Plot_pathways_per_APC_population(pAPC_status, output_directory,input_directory, batch, groups_PCA,mat_pws, pathways)
    
  Plot_genes_per_APC_population<-function(pAPC_status, output_directory,input_directory, batch, groups_PCA){
    ### only include cell types that are observed enough times in the pAPC pool
    w =which(pAPC_status==1)
    t = table(cell_type[w], sample[w])
    t = t[,grep("biopsy", colnames(t))]
    a1 = apply(t[,groups_PCA[[1]]], 1, function(x){length(which(x>=10))})
    a2 = apply(t[,groups_PCA[[2]]], 1, function(x){length(which(x>=10))})
    cell_types_use = names(which(apply(cbind(a1,a2), 1, min)>=3))
    
    samples_use = samples[grep("biopsy", samples)]
    ### mean pathway per patient per cell type
    list_pw = NULL
    for(g in c(1:length(genes))){
      mat_pw_mean = matrix(data = 0, nrow = length(samples_use),ncol = length(cell_types_use), dimnames = c(list(samples_use), list(cell_types_use)))
      for(c in c(1:length(cell_types_use))){
        print(c)
        w = names(which(cell_type==cell_types_use[c]))
        for(s in c(1:length(samples_use))){
          w2 = intersect(w, names(which(sample==samples_use[s])))
          mat_pw_mean[samples_use[s],cell_types_use[c]] = median(mat_genes[w2,genes[g]])
        }
      }
      list_pw = c(list_pw, list(mat_pw_mean))
    }
    names(list_pw) = genes
    
    fileout1= concat(c(output_directory,"Pathway_levels_pAPC_gene_levels_between_groups.pdf"))
    w=3.35
    pdf(file=fileout1, height=w*0.8*4, width=w*1*3)
    par(mfrow= c(4,3), mar = c(12,5,3,3))
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    summary_tables = NULL
    for(ind in c(1:length(list_pw))){
      groups_PCA_use = groups_PCA_all[[1]]
      
      group_PCA_list = NULL
      factor_PCA = NULL
      for(i in c(1:length(groups_PCA_use))){
        group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA_use[[i]]))
        factor_PCA = c(factor_PCA, rep(i, length(groups_PCA_use[[i]])))
      }
      name = names(list_pw)[ind]
      mat = list_pw[[ind]]
      mat=mat[group_PCA_list,]
      a = apply(mat, 2, sd)
      w = intersect(which(is.na(a)==F), which(a>0))
      if(length(w)>=2){
        mat = mat[,w]
        factor = factor(factor_PCA)
        fit = manova(formula = mat ~ factor)
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
          if(length(mean)==0){means = Means_factor(factor, mat[,i1])
          }else{means = rbind(means, Means_factor(factor, mat[,i1]))}
        }
        p_value[which(is.na(p_value))] = 2
        names(p_value) = nam
        sort(p_value)
        colnames(means) = paste("mean.group.", c(1:length(means[1,])))
        combined_p_value = cbind(p_value ,means)
        rownames(combined_p_value) = nam
        p.site = rep(name, length(nam))
        p.analysis = rep(pathways[g], length(nam))
        x = cbind(p.site, p.analysis, combined_p_value)
        if(length(summary_tables)==0){
          summary_tables = x
        }else{summary_tables= rbind(summary_tables, x)}
        
        
        groups = NULL
        p_values = combined_p_value[,1]
        for(i in c(1:length(mat[1,]))){
          g1 = NULL
          data = NULL
          factor = NULL
          for(g in c(1:length(groups_PCA_use))){
            x = mat [gsub("_","-",groups_PCA_use[[g]]), i]
            x = x[which(x!=-1)]
            data = c(data,x )
            factor = c(factor, rep(g, length(x)))
            g1 = c(g1, list(x))}
          groups = c(groups, list(g1))
          # group = factor(factor)
          # fit = aov(formula = data ~ group)
          # pvalues = summary.lm(fit)$ coefficients[,4]
          # pVal <- anova(fit)$'Pr(>F)'[1]
          # p_values = c(p_values, pVal)
        }
        factors1 = paste("group",c(1:length(groups_PCA_use)))
        factors = gsub("_"," ",colnames(mat))
        n = gsub("ANTIGEN.PROCESSING.AND.PRESENTATION", "Ag. pres.", name)
        n = gsub("Negative.regulatory.of.ealry.B.cell.activation", "Neg.reg.early.Bcell.act.", n)
        n = gsub("Negative.regulatory.of.ealry.B.cell.activation", "Neg.reg.early.Bcell.act.", n)
        main = concat(c(n, "\n","pAPC only"))
        max = max(c(unlist(groups), unlist(groups))*1.35)
        min = min(c(unlist(groups), unlist(groups)))
        b = (max-min)*0.034
        min = min-b
        ylab = ""
        draw_signif_lines = TRUE
        y = max(c(unlist(groups), unlist(groups))*1)+b
        max_width = 10
        max_scale = max
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
        if(range <0.001){scale = c(-100:100)*0.0001}
        if(range <0.0001){scale = c(-100:100)*0.00001}
        cex = 0.9
        Fun<-function(x){x}
        scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
        plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
        mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
        mtext(side = 1, text = gsub("B cell antigen-experienced", "B cell Ag.-exp.",factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
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
            points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg=cols1[i1], cex = 0.7)
          }}
        
        for(i in c(1:l)){	
          b = max*0.035
          signif_threshold = 0.05
          if(p_values[i]<signif_threshold){
            pval1 = "*"
            y = max(unlist(groups[[i]]))
            y = y+1*b
            # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
            text(i, y+2*b, labels = pval1, cex = 1.7)
          }
        }
        
      }
    }
    dev.off()
    file=concat(c(output_directory,"Pathway_levels_pAPC_gene_levels_between_groups.txt"))
    write.table(summary_tables, file = file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
  }
  Plot_genes_per_APC_population(pAPC_status, output_directory,input_directory, batch, groups_PCA)
  
}

APC_analysis_PDAC150K(output_directory,input_directory, batch, groups_PCA, VDJ_list_BCR)

######################
batch = "PENG"
APC_analysis_PENG<-function(output_directory,input_directory, batch, groups_PCA){
  type = "All"
  type1 = "APCs"
  file = concat(c(input_directory, "Pathway_levels_per_cell_PENG.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[which(p[,"cell_type"] %in% c("-","doublet/unclear")==F),]
  cell_type_level0 = p[,"cell_type"]
  cell_type_level1 = p[,"cells_type_level1"]
  cell_type_level2 = p[,"cells_type_level2"]
  sample = p[,"sample"]
  samples = sort(unique(sample))
  pathways = c("ANTIGEN.PROCESSING.AND.PRESENTATION.MHC.CLASS.I1",
               "ANTIGEN.PROCESSING.AND.PRESENTATION..MCHII1",     
               "NFKB.PATHWAY1" )
  
  
  mat_pws = NULL
  for(i in c(1:length(pathways))){
    if(length(mat_pws)==0){
      mat_pws = as.numeric(p[,pathways[i]])
    }else{mat_pws=cbind(mat_pws , as.numeric(p[,pathways[i]]))}
  }
  colnames(mat_pws) = pathways
  
  genes = colnames(p)[grep("expr.", colnames(p))]
  mat_genes = NULL
  for(i in c(1:length(genes))){
    if(length(mat_genes)==0){
      mat_genes = as.numeric(p[,genes[i]])
    }else{mat_genes=cbind(mat_genes , as.numeric(p[,genes[i]]))}
  }
  colnames(mat_genes) = genes
  rownames(mat_pws) = rownames(p)
  rownames(mat_genes) = rownames(p)
  
  cell_type = cell_type_level0
  cell_types = sort(unique(cell_type))
  
  for(g in c(2)){
    g = 2
    x = mat_pws[,pathways[g]]
    
    cells_sub = c(grep("T cell CD8", cell_type_level2), grep("DC", cell_type_level2))
    y = x[cells_sub]
    x1 = data.frame(y) ### subset of data
    class = cell_type_level2[cells_sub]
    class[which(class != "T cell CD8")] = "APC"
    class[which(class == "T cell CD8")] = "non-APC"
    x1$class = factor(class)
    
    model <- glm(class ~ ., x1,family=binomial(link='logit'))
    p <- predict(model, newdata=x1, type="response")
    glm.pred=rep("APC", length(p))
    glm.pred[p>0.5]="non-APC"
    length(which(glm.pred==class))*100/length(p)
    table(glm.pred,class)
    ### get decision point (it should happen between the two means of the distributions)
    a1 = mean(x[grep("T cell CD8", cell_type_level2)])
    a2 = mean(c(x[grep("B cell", cell_type_level2)], x[grep("DC", cell_type_level2)]))
    
    x_trial = data.frame(seq(from = a1, to = a2, length = 10000))
    x_trial = data.frame(x_trial)
    colnames(x_trial) = "y"
    p <- predict(model, newdata=x_trial, type="response")
    threshold_pathway = x_trial[,1][max(which(p>0.5))]
    
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",batch,"_",pathways[g],".pdf"))
    w=2.5
    pdf(file=fileout1, height=w*1.1, width=w*1.15)
    par(mfrow= c(1,1), mar = c(5,5,2,3))
    hist(x, breaks = 100, border = NA, col = "grey", freq = F, ylim = c(0,1.8), xlab = "score", xlim = c(-1,3))
    h1 = hist(x[grep("T cell CD8", cell_type_level2)], breaks = 50, plot = F)
    points(h1$mids,h1$density , col = add.alpha("red", 0.6),  type = "l", lwd = 4)
    
    h1 = hist(x[grep("DC", cell_type_level2)], breaks = 50, plot = F)
    points(h1$mids,h1$density , col = add.alpha("lightskyblue", 0.6),  type = "l", lwd = 4)
    
    segments(threshold_pathway, -10, threshold_pathway, 1000, col = "black", lwd= 2, lty = 2)
    dev.off()
    
    library("mixtools")
    library(tibble)
    library(ggplot2)
    
    plot_mix_comps <- function(x, mu, sigma, lam) {
      lam * dnorm(x, mu, sigma)
    }
    
    mixmdl <- normalmixEM(x, k = 2)
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",batch,"_",pathways[g],"1.pdf"))
    w=3.2
    pdf(file=fileout1, height=w*1, width=w*1.9)
    par(mfrow= c(1,1), mar = c(5,5,2,3))
    
    data.frame(x = mixmdl$x) %>%
      ggplot() +
      geom_histogram(aes(x, ..density..), binwidth = 0.05, colour = "black", 
                     fill = "white") +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                    colour = "red", lwd = 1.5) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                    colour = "blue", lwd = 1.5) +
      ylab("Density") + 
      xlab(concat(c(pathways[g], " score")))+
      geom_vline(xintercept= c(a <- c(threshold_pathway)), 
                 colour = c("red")[c(1)])
    dev.off()
    
    ### proportions of APCs per cell type per sample in the tumour
    y = rep(0, length(x))
    y[which(x>threshold_pathway)] = 1
    t = table(cell_type_level2,y)
    percentage_pos = t[,2]*100/rowSums(t)
    type_names = names(percentage_pos)
    col = rep(1, length(type_names))
    col[grep("DC", type_names)] =2
    col[grep("MoMac", type_names)] =3
    col[grep("Mono", type_names)] =4
    col[grep("granulocyte", type_names)] =5
    col[grep("NK", type_names)] =6
    col[grep("T cell", type_names)] =7
    col[grep("Mast", type_names)] =8
    col[grep("ILC", type_names)] =8
    col[grep("iNKT", type_names)] =8
    col[grep("Megakaryocyte", type_names)] =8
    names(col)=type_names
    names(percentage_pos) = gsub("B cell antigen-experienced", "B cell Ag.-exp.",names(percentage_pos))
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",batch,"_",pathways[g],"_proportions_overall.pdf"))
    w=3.2
    pdf(file=fileout1, height=w*1.1, width=w*1.9)
    par(mfrow= c(1,1), mar = c(12,5,2,3))
    barplot(percentage_pos, las = 2, col = col, ylab = "% positive", border = NA)
    dev.off()
    
    ##### get broad proportion of APCs per sample
    cell_type = cell_type_level2
    cell_types= sort(unique(cell_type))
    mat_prop_of_APC = matrix(data = 0, nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
    for(s in c(1:length(samples))){
      w = intersect(which(sample==samples[s]), which(y==1))
      if(length(w)>=10){
        t = table(cell_type[w])
        t = t*100/sum(t)
        mat_prop_of_APC[samples[s],names(t)] = t
      }
    }
    
    ## difference between groups
    groups_PCA = readRDS(file = concat(c(input_directory_groups, "Immunogroups_2019_PUBLISHED_PDAC.rds")))
   
    groups_PCA_all = c(list(groups_PCA))
    names(groups_PCA_all) = c("tumour-PENG")
    
    
    a = cbind(apply(mat_prop_of_APC[groups_PCA[[1]],], 2, mean), apply(mat_prop_of_APC[groups_PCA[[2]],], 2, mean))
    a1 = a[which(rownames(a) %in% c("Acinar","Centroacinar","Endocrine","Endothelial cell","Epithelial","Fibroblast","Stellate cell")==F),]
    a1[,1] = a1[,1]*100/sum(a1[,1])
    a1[,2] = a1[,2]*100/sum(a1[,2])
    b = colSums(a[grep("DC", rownames(a1)),])
    c = colSums(a[grep("B cell", rownames(a1)),])
    d =rbind(b, c, a1["MoMac",])
    rownames(d) = c("DC", "B cell", "MoMac")
    colnames(d) = c("Mean cellular contribution to pAPCs (ME)","Mean cellular contribution to pAPCs (AE)")
    
    e = rbind(d, a)
    
    fileout1= concat(c(output_directory,"Pathway_levels_distribution_modelling_",batch,"_",pathways[g],"_proportions_between_groups.txt"))
    write.table(e, file = fileout1, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    fileout1= concat(c(output_directory, "Pathway_levels_distribution_modelling_",batch,"_",pathways[g],"_proportions_between_groups_",batch,".pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*2.5, width=w*1*2.3)
    par(mfrow= c(2,1), mar = c(12,5,3,3))
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(ind in c(1:length(groups_PCA_all))){
      groups_PCA_use = groups_PCA_all[[ind]]
      
      group_PCA_list = NULL
      factor_PCA = NULL
      for(i in c(1:length(groups_PCA_use))){
        group_PCA_list = c(group_PCA_list, gsub("_","_",groups_PCA_use[[i]]))
        factor_PCA = c(factor_PCA, rep(i, length(groups_PCA_use[[i]])))
      }
      name = names(groups_PCA_all)[ind]
      mat = mat_prop_of_APC[group_PCA_list,]
      factor = factor(factor_PCA)
      fit = manova(formula = mat ~ factor)
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
        if(length(mean)==0){means = Means_factor(factor, mat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      names(p_value) = nam
      sort(p_value)
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.site = rep(name, length(nam))
      p.analysis = rep(pathways[g], length(nam))
      x = cbind(p.site, p.analysis, combined_p_value)
      summary_tables = x
      
      groups = NULL
      p_values = combined_p_value[,1]
      for(i in c(1:length(mat[1,]))){
        g1 = NULL
        data = NULL
        factor = NULL
        for(g in c(1:length(groups_PCA_use))){
          x = mat [gsub("_","_",groups_PCA_use[[g]]), i]
          x = x[which(x!=-1)]
          data = c(data,x )
          factor = c(factor, rep(g, length(x)))
          g1 = c(g1, list(x))}
        groups = c(groups, list(g1))
        # group = factor(factor)
        # fit = aov(formula = data ~ group)
        # pvalues = summary.lm(fit)$ coefficients[,4]
        # pVal <- anova(fit)$'Pr(>F)'[1]
        # p_values = c(p_values, pVal)
      }
      factors1 = paste("group",c(1:length(groups_PCA_use)))
      factors = gsub("_"," ",colnames(mat))
      main = concat(c(pathways[g], " ",name))
      max = max(c(unlist(groups), unlist(groups))*1.35)
      min = 0
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 24
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
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = gsub("B cell antigen-experienced", "B cell Ag.-exp.",factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
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
          points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg=cols1[i1], cex = 0.7)
        }}
      
      for(i in c(1:l)){	
        b = max*0.035
        signif_threshold = 0.05
        if(p_values[i]<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))
          y = y+1*b
          # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
          text(i, y+2*b, labels = pval1, cex = 1.7)
        }
      }
      
    }
    dev.off()
    
  }
  
  ### only include cell types that are observed enough times in the pAPC pool
  w =which(pAPC_status==1)
  t = table(cell_type[w], sample[w])
  t = t[,grep("PDAC_T", colnames(t))]
  a1 = apply(t[,groups_PCA[[1]]], 1, function(x){length(which(x>=10))})
  a2 = apply(t[,groups_PCA[[2]]], 1, function(x){length(which(x>=10))})
  cell_types_use = names(which(apply(cbind(a1,a2), 1, min)>=3))
  
  samples_use = samples[grep("PDAC_T", samples)]
  ### mean pathway per patient per cell type
  list_pw = NULL
  for(g in c(1:length(genes))){
    mat_pw_mean = matrix(data = 0, nrow = length(samples_use),ncol = length(cell_types_use), dimnames = c(list(samples_use), list(cell_types_use)))
    for(c in c(1:length(cell_types_use))){
      print(c)
      w = names(which(cell_type==cell_types_use[c]))
      for(s in c(1:length(samples_use))){
        w2 = intersect(w, names(which(sample==samples_use[s])))
        mat_pw_mean[samples_use[s],cell_types_use[c]] = median(mat_genes[w2,genes[g]])
      }
    }
    list_pw = c(list_pw, list(mat_pw_mean))
  }
  names(list_pw) = genes
  
  fileout1= concat(c(output_directory,"Pathway_levels_pAPC_pathway_levels_between_groups_",batch,".pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*4, width=w*1*3)
  par(mfrow= c(4,3), mar = c(12,5,3,3))
  
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  for(ind in c(1:length(list_pw))){
    groups_PCA_use = groups_PCA_all[[1]]
    
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA_use))){
      group_PCA_list = c(group_PCA_list, gsub("-","_",groups_PCA_use[[i]]))
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA_use[[i]])))
    }
    name = names(list_pw)[ind]
    mat = list_pw[[ind]]
    mat=mat[group_PCA_list,]
    factor = factor(factor_PCA)
    fit = manova(formula = mat ~ factor)
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
      if(length(mean)==0){means = Means_factor(factor, mat[,i1])
      }else{means = rbind(means, Means_factor(factor, mat[,i1]))}
    }
    p_value[which(is.na(p_value))] = 2
    names(p_value) = nam
    sort(p_value)
    colnames(means) = paste("mean.group.", c(1:length(means[1,])))
    combined_p_value = cbind(p_value ,means)
    rownames(combined_p_value) = nam
    p.site = rep(name, length(nam))
    p.analysis = rep(pathways[g], length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    summary_tables = x
    
    groups = NULL
    p_values = combined_p_value[,1]
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA_use))){
        x = mat [gsub("-","_",groups_PCA_use[[g]]), i]
        x = x[which(x!=-1)]
        data = c(data,x )
        factor = c(factor, rep(g, length(x)))
        g1 = c(g1, list(x))}
      groups = c(groups, list(g1))
      # group = factor(factor)
      # fit = aov(formula = data ~ group)
      # pvalues = summary.lm(fit)$ coefficients[,4]
      # pVal <- anova(fit)$'Pr(>F)'[1]
      # p_values = c(p_values, pVal)
    }
    factors1 = paste("group",c(1:length(groups_PCA_use)))
    factors = gsub("_"," ",colnames(mat))
    main = concat(c(name, "\n","pAPC only"))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = min(c(unlist(groups), unlist(groups)))
    b = (max-min)*0.034
    min = min-b
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 10
    max_scale = max
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
    if(range <0.001){scale = c(-100:100)*0.0001}
    if(range <0.0001){scale = c(-100:100)*0.00001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = gsub("B cell antigen-experienced", "B cell Ag.-exp.",factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
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
        points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg=cols1[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        y = max(unlist(groups[[i]]))
        y = y+1*b
        # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
        text(i, y+2*b, labels = pval1, cex = 1.7)
      }
    }
    
  }
  dev.off()

  
}
APC_analysis_PENG(output_directory,input_directory, batch, groups_PCA)

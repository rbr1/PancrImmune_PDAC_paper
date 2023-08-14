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
analysis = "VDJ_clonality"

######################
## get patient groups for plotting: 
input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"

input_directory_groups  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"
groups_PCA = readRDS(file = concat(c(input_directory_groups, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))

library(yarrr)
library(RColorBrewer)

Plot_txt_files_BCR<-function(input_directory,output_directory,groups_PCA, batch,chain){
  chain = "BCR"
  types = c("nduplicate_clones_intra", "nduplicate_clones_inter","gini", "renyi","d5")
  files = apply(cbind(input_directory, "VDJ_Clonality_information_",chain,"_",batch,"_", types,".txt"), 1, paste, collapse = "")
  names(files) = types
  
   
  
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  
  Means_factor = function(factor, x){
    m = NULL
    for(i1 in c(1:length(levels(factor)))){
      x1 = x[which(factor==levels(factor)[i1])]
      x1 = x1[which(x1!=-1)]
      m = c(m, mean(x1))}
    return(m)}
  
  fileout1=concat(c(output_directory, "VDJ_Clonality_information_",chain,"_boxplots_between_groups_", batch,"_tumour.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*5, width=w*1*1.8)
  par(mfrow= c(5,1), mar = c(15,5,3,3))
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  summary_tables = NULL
  for(f in c(1:length(files))){
    analysis = names(files)[f]
    p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
    rownames(p) = gsub("-","_",rownames(p))
    mat = p[group_PCA_list,]
    mat[which(mat==-1)]= NA
    a = apply(mat, 2, function(x){length(which(is.na(x)))})
    mat = mat[,which(a<=4)]
    
    ## stats here
    factor = factor(factor_PCA)
    fit = manova(formula = mat ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
   
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
    p.site = rep("biopsy", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    if(length(summary_tables)==0){summary_tables = x
    }else{summary_tables = rbind(summary_tables, x)}
    
    ## plot here
    groups = NULL
    p_values = combined_p_value[,1]
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA))){
        x = mat[groups_PCA[[g]], i]
        x = x[which(x!=-1)]
        data = c(data,x )
        factor = c(factor, rep(g, length(x)))
        g1 = c(g1, list(x))}
      groups = c(groups, list(g1))
    }
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = gsub("."," ",colnames(mat), fixed = T)
    main = concat(c(analysis, " ", chain))
    max = 90
    if(max(unlist(groups))<35){max = 40}
    #max = max(c(unlist(groups), unlist(groups))*1.35)
    
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 26
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*10}
    if(range<=20){scale = c(0:100)*5}
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
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.16
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.2/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg = cols[i1], col=cols1[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        if(p_values[i] <signif_threshold/10){pval1 = "*"}
        y = max(unlist(groups[[i]]))
        y = y+1*b
        text(i, y+2*b, labels = pval1, cex = 1.7)
      }
    }
    
  }
  dev.off()
  
  
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("biopsy","blood",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  
  fileout1=concat(c(output_directory, "VDJ_Clonality_information_",chain,"_boxplots_between_groups_", batch,"_blood.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*5, width=w*1*1.8)
  par(mfrow= c(5,1), mar = c(15,5,3,3))
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  summary_tables = NULL
  for(f in c(1:length(files))){
    analysis = names(files)[f]
    p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
    rownames(p) = gsub("-","_",rownames(p))
    mat = p[gsub("biopsy","blood",group_PCA_list),]
    mat[which(mat==-1)]= NA
    a = apply(mat, 2, function(x){length(which(is.na(x)))})
    mat = mat[,which(a<=4)]
    
    ## stats here
    factor = factor(factor_PCA)
    fit = manova(formula = mat ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
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
    p.site = rep("biopsy", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    if(length(summary_tables)==0){summary_tables = x
    }else{summary_tables = rbind(summary_tables, x)}
    
    ## plot here
    groups = NULL
    p_values = combined_p_value[,1]
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA))){
        x = mat[gsub("biopsy","blood",groups_PCA[[g]]), i]
        x = x[which(x!=-1)]
        data = c(data,x )
        factor = c(factor, rep(g, length(x)))
        g1 = c(g1, list(x))}
      groups = c(groups, list(g1))
    }
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = gsub("."," ",colnames(mat), fixed = T)
    main = concat(c(analysis, " ", chain))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    max = 90
    if(max(unlist(groups))<35){max = 40}
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 26
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*10}
    if(range<=20){scale = c(0:100)*5}
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
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.16
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.2/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg = cols[i1], col=cols1[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        if(p_values[i] <signif_threshold/10){pval1 = "*"}
        y = max(unlist(groups[[i]]))
        y = y+1*b
        text(i, y+2*b, labels = pval1, cex = 1.7)
      }
    }
    
  }
  dev.off()

}
Plot_txt_files_BCR(input_directory,output_directory,groups_PCA, batch,chain)
  
Plot_rds_files_BCR<-function(input_directory,output_directory,groups_PCA, batch,chain){
  chain = "BCR"
  file = concat(c(input_directory, "VDJ_repertoire_feature_information_",chain,"_",batch,".rds"))
  VDJ = readRDS(file= file)

  
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  
  Means_factor = function(factor, x){
    m = NULL
    for(i1 in c(1:length(levels(factor)))){
      x1 = x[which(factor==levels(factor)[i1])]
      x1 = x1[which(x1!=-1)]
      m = c(m, mean(x1))}
    return(m)}
  
  fileout1=concat(c(output_directory, "VDJ_Clonality_information_",chain,"_boxplots_between_groups_", batch,"_tumour2.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*5, width=w*1*20.8*2)
  par(mfrow= c(5,2), mar = c(15,5,3,3))
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  summary_tables = NULL
  for(f in c(1:length(VDJ))){
    p =VDJ[[f]]
    name = names(VDJ)[f]
    if(length(grep("usage", name))!=0){ ## need to normalise
      for(c in c(1:length(p))){
        m = p[[c]]
        m=m[,which(colnames(m)!='-')]
        for(i in c(1:length(m[,1]))){
          if(sum(m[i,])>5){m[i,]=m[i,]*100/sum(m[i,])
          }else{m[i,] = NA}
        }
        p[[c]] = m
      }
    }else{
      for(c in c(1:length(p))){
        m = p[[c]]
        m=m[,which(colnames(m)!='-')]
        p[[c]] = m
      }
    }
    for(c in c(1:length(p))){
      mat = p[[c]]
      rownames(mat) = gsub("-","_",rownames(mat))
      mat = mat[group_PCA_list,]
      mat[which(mat==-1)]= NA
      a = apply(mat, 2, function(x){length(which(is.na(x)))})
      if(length(which(a<=4))>=2){
        mat = mat[,which(a<=4)]
        analysis = concat(c(name," ",names(p)[c]))
        ## stats here
        factor = factor(factor_PCA)
        fit = manova(formula = mat ~ factor)
        p1 = summary.aov(fit)
        nam = gsub(" Response ","",names(p1))
        p_value = NULL
        means = NULL
        
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
        p.site = rep("biopsy", length(nam))
        p.analysis = rep(analysis, length(nam))
        x = cbind(p.site, p.analysis, combined_p_value)
        if(length(summary_tables)==0){summary_tables = x
        }else{summary_tables = rbind(summary_tables, x)}
        
        ## plot here
        groups = NULL
        p_values = combined_p_value[,1]
        for(i in c(1:length(mat[1,]))){
          g1 = NULL
          data = NULL
          factor = NULL
          for(g in c(1:length(groups_PCA))){
            x = mat[groups_PCA[[g]], i]
            x = x[which(x!=-1)]
            data = c(data,x )
            factor = c(factor, rep(g, length(x)))
            g1 = c(g1, list(x))}
          groups = c(groups, list(g1))
        }
        factors1 = paste("group",c(1:length(groups_PCA)))
        factors = gsub("."," ",colnames(mat), fixed = T)
        main = concat(c(analysis, " ", chain))
        max = max(c(unlist(groups), unlist(groups))*1.35)
        min = 0
        b = (max-min)*0.034
        ylab = ""
        draw_signif_lines = TRUE
        y = max(c(unlist(groups), unlist(groups))*1)+b
        max_width = 500
        max_scale = max
        range = max-min
        if(range>50){scale = c(0:100)*20}
        if(range<=50){scale = c(0:100)*10}
        if(range<=20){scale = c(0:100)*5}
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
        plot(c(19, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
        mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
        mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
        segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
        mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
        width = 0.16
        index = 1
        l = length(groups)
        l1 = length(groups[[1]])
        shift = c(1:l1)
        shift = (mean(shift)-shift)
        shift = shift*0.2/max(shift)
        
        for(i in c(1:l)){
          for(i1 in c(1:l1)){
            points1=as.numeric(groups[[i]][[i1]])
            box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
            Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
            points(rep(i-shift[i1], length(points1)),points1, pch =21, bg = cols[i1], col=cols1[i1], cex = 0.7)
          }}
        
        for(i in c(1:l)){	
          b = max*0.035
          signif_threshold = 0.05
          if(p_values[i]<signif_threshold){
            pval1 = "*"
            if(p_values[i] <signif_threshold/10){pval1 = "*"}
            y = max(unlist(groups[[i]]))
            y = y+1*b
            text(i, y+2*b, labels = pval1, cex = 1.7)
          }
        }
      }
    }
  }
  dev.off()
  
  
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("biopsy","blood",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  

  fileout1=concat(c(output_directory, "VDJ_Clonality_information_",chain,"_boxplots_between_groups_", batch,"_blood2.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*5, width=w*1*20.8*2)
  par(mfrow= c(5,2), mar = c(15,5,3,3))
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  summary_tables = NULL
  for(f in c(1:length(VDJ))){
    p =VDJ[[f]]
    name = names(VDJ)[f]
    if(length(grep("usage", name))!=0){ ## need to normalise
      for(c in c(1:length(p))){
        m = p[[c]]
        m=m[,which(colnames(m)!='-')]
        for(i in c(1:length(m[,1]))){
          if(sum(m[i,])>5){m[i,]=m[i,]*100/sum(m[i,])
          }else{m[i,] = NA}
        }
        p[[c]] = m
      }
    }else{
      for(c in c(1:length(p))){
        m = p[[c]]
        m=m[,which(colnames(m)!='-')]
        p[[c]] = m
      }
    }
    for(c in c(1:length(p))){
      mat = p[[c]]
      rownames(mat) = gsub("-","_",rownames(mat))
      mat = mat[group_PCA_list,]
      mat[which(mat==-1)]= NA
      a = apply(mat, 2, function(x){length(which(is.na(x)))})
      if(length(which(a<=4))>=2){
        mat = mat[,which(a<=4)]
        analysis = concat(c(name," ",names(p)[c]))
        ## stats here
        factor = factor(factor_PCA)
        fit = manova(formula = mat ~ factor)
        p1 = summary.aov(fit)
        nam = gsub(" Response ","",names(p1))
        p_value = NULL
        means = NULL
        
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
        p.site = rep("biopsy", length(nam))
        p.analysis = rep(analysis, length(nam))
        x = cbind(p.site, p.analysis, combined_p_value)
        if(length(summary_tables)==0){summary_tables = x
        }else{summary_tables = rbind(summary_tables, x)}
        
        ## plot here
        groups = NULL
        p_values = combined_p_value[,1]
        for(i in c(1:length(mat[1,]))){
          g1 = NULL
          data = NULL
          factor = NULL
          for(g in c(1:length(groups_PCA))){
            x = mat[gsub("biopsy","blood",groups_PCA[[g]]), i]
            x = x[which(x!=-1)]
            data = c(data,x )
            factor = c(factor, rep(g, length(x)))
            g1 = c(g1, list(x))}
          groups = c(groups, list(g1))
        }
        factors1 = paste("group",c(1:length(groups_PCA)))
        factors = gsub("."," ",colnames(mat), fixed = T)
        main = concat(c(analysis, " ", chain))
        max = max(c(unlist(groups), unlist(groups))*1.35)
        min = 0
        b = (max-min)*0.034
        ylab = ""
        draw_signif_lines = TRUE
        y = max(c(unlist(groups), unlist(groups))*1)+b
        max_width = 500
        max_scale = max
        range = max-min
        if(range>50){scale = c(0:100)*20}
        if(range<=50){scale = c(0:100)*10}
        if(range<=20){scale = c(0:100)*5}
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
        plot(c(19, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
        mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
        mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
        segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
        mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
        width = 0.16
        index = 1
        l = length(groups)
        l1 = length(groups[[1]])
        shift = c(1:l1)
        shift = (mean(shift)-shift)
        shift = shift*0.2/max(shift)
        
        for(i in c(1:l)){
          for(i1 in c(1:l1)){
            points1=as.numeric(groups[[i]][[i1]])
            box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
            Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
            points(rep(i-shift[i1], length(points1)),points1, pch =21, bg = cols[i1], col=cols1[i1], cex = 0.7)
          }}
        
        for(i in c(1:l)){	
          b = max*0.035
          signif_threshold = 0.05
          if(p_values[i]<signif_threshold){
            pval1 = "*"
            if(p_values[i] <signif_threshold/10){pval1 = "*"}
            y = max(unlist(groups[[i]]))
            y = y+1*b
            text(i, y+2*b, labels = pval1, cex = 1.7)
          }
        }
      }
    }
  }
  dev.off()
  
  
  
}
Plot_rds_files_BCR(input_directory,output_directory,groups_PCA, batch,chain)

Plot_txt_files_TCR<-function(input_directory,output_directory,groups_PCA, batch,chain){
  overall_chain = "TCR"
  chain = c("TCR_CD4", "TCR_CD8")
  types = c("nduplicate_clones_intra", "nduplicate_clones_inter")
  files = c(apply(cbind(input_directory, "VDJ_Clonality_information_",chain,"_",batch,"_", types[1],".txt"), 1, paste, collapse = ""), 
            apply(cbind(input_directory, "VDJ_Clonality_information_",chain,"_",batch,"_", types[2],".txt"), 1, paste, collapse = ""))
  names(files) = c(apply(cbind(chain, types[1]), 1, paste, collapse = "_"), apply(cbind(chain, types[2]), 1, paste, collapse = "_"))
  chains_plot = c(chain, chain)
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  
  Means_factor = function(factor, x){
    m = NULL
    for(i1 in c(1:length(levels(factor)))){
      x1 = x[which(factor==levels(factor)[i1])]
      x1 = x1[which(x1!=-1)]
      m = c(m, mean(x1))}
    return(m)}
  
  fileout1=concat(c(output_directory, "VDJ_Clonality_information_",overall_chain,"_boxplots_between_groups_", batch,"_tumour.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*5, width=w*1*1.8)
  par(mfrow= c(5,1), mar = c(15,5,3,3))
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  summary_tables = NULL
  for(f in c(1:length(files))){
    analysis = names(files)[f]
    p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
    rownames(p) = gsub("-","_",rownames(p))
    mat = p[group_PCA_list,]
    mat[which(mat==-1)]= NA
    
    if(chains_plot[f]=="TCR_CD4"){ mat = mat[,grep("CD4", colnames(mat))]}
    if(chains_plot[f]=="TCR_CD8"){ mat = mat[,c(grep("CD8", colnames(mat)), grep("MAIT", colnames(mat)))]}
    
    ## stats here
    a = apply(mat, 2, function(x){length(which(is.na(x)))})
    mat1 = mat[,which(a<=4)]
    factor = factor(factor_PCA)
    fit = manova(formula = mat1^0.5 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
    
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
    p.site = rep("biopsy", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    if(length(summary_tables)==0){summary_tables = x
    }else{summary_tables = rbind(summary_tables, x)}
    
    ## plot here
    groups = NULL
    p_values = combined_p_value[,1]
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA))){
        x = mat[groups_PCA[[g]], i]
        x[which(x==-1)] = 0
        x[which(is.na(x))] = 0
        data = c(data,x )
        factor = c(factor, rep(g, length(x)))
        g1 = c(g1, list(x))}
      groups = c(groups, list(g1))
    }
    names(groups ) = colnames(mat)
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = gsub("."," ",colnames(mat), fixed = T)
    main = concat(c(analysis))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    max = 90
    if(max(unlist(groups))<35){max = 40}
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 26
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*10}
    if(range<=20){scale = c(0:100)*5}
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
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.16
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.2/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg = cols[i1], col=cols1[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(colnames(mat)[i] %in% names(p_value)){
        p = p_value[colnames(mat)[i]]
        if(p<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))
          y = y+1*b
          text(i, y+2*b, labels = pval1, cex = 1.7)
        }}
    }
    
  }
  dev.off()
  
  
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("biopsy","blood",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  
  fileout1=concat(c(output_directory, "VDJ_Clonality_information_",overall_chain,"_boxplots_between_groups_", batch,"_blood.pdf"))
  w=3.35
  pdf(file=fileout1, height=w*1*5, width=w*1*1.8)
  par(mfrow= c(5,1), mar = c(15,5,3,3))
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  summary_tables = NULL
  for(f in c(1:length(files))){
    analysis = names(files)[f]
    p <- as.matrix(read.csv(files[f], head=T, sep="\t"))
    rownames(p) = gsub("-","_",rownames(p))
    mat = p[gsub("biopsy","blood",group_PCA_list),]
    mat[which(mat==-1)]= NA
    
    if(chains_plot[f]=="TCR_CD4"){ mat = mat[,grep("CD4", colnames(mat))]}
    if(chains_plot[f]=="TCR_CD8"){ mat = mat[,c(grep("CD8", colnames(mat)), grep("MAIT", colnames(mat)))]}
    
    ## stats here
    a = apply(mat, 2, function(x){length(which(is.na(x)))})
    mat1 = mat[,which(a<=4)]
    factor = factor(factor_PCA)
    fit = manova(formula = mat1^0.5 ~ factor)
    p1 = summary.aov(fit)
    nam = gsub(" Response ","",names(p1))
    p_value = NULL
    means = NULL
   
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
    p.site = rep("biopsy", length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    if(length(summary_tables)==0){summary_tables = x
    }else{summary_tables = rbind(summary_tables, x)}
    
    ## plot here
    groups = NULL
    p_values = combined_p_value[,1]
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      data = NULL
      factor = NULL
      for(g in c(1:length(groups_PCA))){
        x = mat[gsub("biopsy","blood",groups_PCA[[g]]), i]
        x[which(x==-1)] = 0
        x[which(is.na(x))] = 0
        data = c(data,x )
        factor = c(factor, rep(g, length(x)))
        g1 = c(g1, list(x))}
      groups = c(groups, list(g1))
    }
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = gsub("."," ",colnames(mat), fixed = T)
    main = concat(c(analysis, " ", chain))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    max = 90
    if(max(unlist(groups))<35){max = 40}
    min = 0
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 26
    max_scale = max
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range<=50){scale = c(0:100)*10}
    if(range<=20){scale = c(0:100)*5}
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
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.16
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.2/max(shift)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg = cols[i1], col=cols1[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      if(colnames(mat)[i] %in% names(p_value)){
        p = p_value[colnames(mat)[i]]
        if(p<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))
          y = y+1*b
          text(i, y+2*b, labels = pval1, cex = 1.7)
        }}
    }
    
  }
  dev.off()
  
}
Plot_txt_files_TCR(input_directory,output_directory,groups_PCA, batch,chain)

## not written
Plot_rds_files_TCR<-function(input_directory,output_directory,groups_PCA, batch,chain){
  overall_chain = "TCR"
  
  chains = c("TCR_CD4", "TCR_CD8")
  for(ch in c(1:length(chains))){
    chain = chains[ch]
    file = concat(c(input_directory, "VDJ_repertoire_feature_information_",chain,"_",batch,".rds"))
    VDJ = readRDS(file= file)
    
    
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    
    Means_factor = function(factor, x){
      m = NULL
      for(i1 in c(1:length(levels(factor)))){
        x1 = x[which(factor==levels(factor)[i1])]
        x1 = x1[which(x1!=-1)]
        m = c(m, mean(x1))}
      return(m)}
    
    fileout1=concat(c(output_directory, "VDJ_Clonality_information_",chain,"_boxplots_between_groups_", batch,"_tumour2.pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*5, width=w*1*20.8*2)
    par(mfrow= c(5,2), mar = c(15,5,3,3))
    
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    summary_tables = NULL
    for(f in c(1:length(VDJ))){
      p =VDJ[[f]]
      name = names(VDJ)[f]
      if(length(grep("usage", name))!=0){ ## need to normalise
        for(c in c(1:length(p))){
          m = p[[c]]
          m=m[,which(colnames(m)!='-')]
          for(i in c(1:length(m[,1]))){
            if(sum(m[i,])>5){m[i,]=m[i,]*100/sum(m[i,])
            }else{m[i,] = NA}
          }
          p[[c]] = m
        }
      }else{
        for(c in c(1:length(p))){
          m = p[[c]]
          m=m[,which(colnames(m)!='-')]
          p[[c]] = m
        }
      }
      for(c in c(1:length(p))){
        mat = p[[c]]
        rownames(mat) = gsub("-","_",rownames(mat))
        mat = mat[group_PCA_list,]
        mat[which(mat==-1)]= NA
        a = apply(mat, 2, function(x){length(which(is.na(x)))})
        if(length(which(a<=4))>=2){
          mat = mat[,which(a<=4)]
          analysis = concat(c(name," ",names(p)[c]))
          ## stats here
          factor = factor(factor_PCA)
          fit = manova(formula = mat ~ factor)
          p1 = summary.aov(fit)
          nam = gsub(" Response ","",names(p1))
          p_value = NULL
          means = NULL
          
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
          p.site = rep("biopsy", length(nam))
          p.analysis = rep(analysis, length(nam))
          x = cbind(p.site, p.analysis, combined_p_value)
          if(length(summary_tables)==0){summary_tables = x
          }else{summary_tables = rbind(summary_tables, x)}
          
          ## plot here
          groups = NULL
          p_values = combined_p_value[,1]
          for(i in c(1:length(mat[1,]))){
            g1 = NULL
            data = NULL
            factor = NULL
            for(g in c(1:length(groups_PCA))){
              x = mat[groups_PCA[[g]], i]
              x = x[which(x!=-1)]
              data = c(data,x )
              factor = c(factor, rep(g, length(x)))
              g1 = c(g1, list(x))}
            groups = c(groups, list(g1))
          }
          factors1 = paste("group",c(1:length(groups_PCA)))
          factors = gsub("."," ",colnames(mat), fixed = T)
          main = concat(c(analysis, " ", chain))
          max = max(c(unlist(groups), unlist(groups))*1.35)
          min = 0
          b = (max-min)*0.034
          ylab = ""
          draw_signif_lines = TRUE
          y = max(c(unlist(groups), unlist(groups))*1)+b
          max_width = 500
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
          plot(c(19, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
          mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
          mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
          segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
          mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
          width = 0.16
          index = 1
          l = length(groups)
          l1 = length(groups[[1]])
          shift = c(1:l1)
          shift = (mean(shift)-shift)
          shift = shift*0.2/max(shift)
          
          for(i in c(1:l)){
            for(i1 in c(1:l1)){
              points1=as.numeric(groups[[i]][[i1]])
              box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
              Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
              points(rep(i-shift[i1], length(points1)),points1, pch =21, bg = cols[i1], col=cols1[i1], cex = 0.7)
            }}
          
          for(i in c(1:l)){	
            b = max*0.035
            signif_threshold = 0.05
            if(p_values[i]<signif_threshold){
              pval1 = "*"
              if(p_values[i] <signif_threshold/10){pval1 = "*"}
              y = max(unlist(groups[[i]]))
              y = y+1*b
              text(i, y+2*b, labels = pval1, cex = 1.7)
            }
          }
        }
      }
    }
    dev.off()
    
    
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, gsub("biopsy","blood",groups_PCA[[i]]))
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    
    
    fileout1=concat(c(output_directory, "VDJ_Clonality_information_",chain,"_boxplots_between_groups_", batch,"_blood2.pdf"))
    w=3.35
    pdf(file=fileout1, height=w*1*5, width=w*1*20.8*2)
    par(mfrow= c(5,2), mar = c(15,5,3,3))
    
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    summary_tables = NULL
    for(f in c(1:length(VDJ))){
      p =VDJ[[f]]
      name = names(VDJ)[f]
      if(length(grep("usage", name))!=0){ ## need to normalise
        for(c in c(1:length(p))){
          m = p[[c]]
          m=m[,which(colnames(m)!='-')]
          for(i in c(1:length(m[,1]))){
            if(sum(m[i,])>5){m[i,]=m[i,]*100/sum(m[i,])
            }else{m[i,] = NA}
          }
          p[[c]] = m
        }
      }else{
        for(c in c(1:length(p))){
          m = p[[c]]
          m=m[,which(colnames(m)!='-')]
          p[[c]] = m
        }
      }
      for(c in c(1:length(p))){
        mat = p[[c]]
        rownames(mat) = gsub("-","_",rownames(mat))
        mat = mat[group_PCA_list,]
        mat[which(mat==-1)]= NA
        a = apply(mat, 2, function(x){length(which(is.na(x)))})
        if(length(which(a<=4))>=2){
          mat = mat[,which(a<=4)]
          analysis = concat(c(name," ",names(p)[c]))
          ## stats here
          factor = factor(factor_PCA)
          fit = manova(formula = mat ~ factor)
          p1 = summary.aov(fit)
          nam = gsub(" Response ","",names(p1))
          p_value = NULL
          means = NULL
          
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
          p.site = rep("biopsy", length(nam))
          p.analysis = rep(analysis, length(nam))
          x = cbind(p.site, p.analysis, combined_p_value)
          if(length(summary_tables)==0){summary_tables = x
          }else{summary_tables = rbind(summary_tables, x)}
          
          ## plot here
          groups = NULL
          p_values = combined_p_value[,1]
          for(i in c(1:length(mat[1,]))){
            g1 = NULL
            data = NULL
            factor = NULL
            for(g in c(1:length(groups_PCA))){
              x = mat[gsub("biopsy","blood",groups_PCA[[g]]), i]
              x = x[which(x!=-1)]
              data = c(data,x )
              factor = c(factor, rep(g, length(x)))
              g1 = c(g1, list(x))}
            groups = c(groups, list(g1))
          }
          factors1 = paste("group",c(1:length(groups_PCA)))
          factors = gsub("."," ",colnames(mat), fixed = T)
          main = concat(c(analysis, " ", chain))
          max = max(c(unlist(groups), unlist(groups))*1.35)
          min = 0
          b = (max-min)*0.034
          ylab = ""
          draw_signif_lines = TRUE
          y = max(c(unlist(groups), unlist(groups))*1)+b
          max_width = 500
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
          plot(c(19, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
          mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
          mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
          segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
          mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
          width = 0.16
          index = 1
          l = length(groups)
          l1 = length(groups[[1]])
          shift = c(1:l1)
          shift = (mean(shift)-shift)
          shift = shift*0.2/max(shift)
          
          for(i in c(1:l)){
            for(i1 in c(1:l1)){
              points1=as.numeric(groups[[i]][[i1]])
              box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
              Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
              points(rep(i-shift[i1], length(points1)),points1, pch =21, bg = cols[i1], col=cols1[i1], cex = 0.7)
            }}
          
          for(i in c(1:l)){	
            b = max*0.035
            signif_threshold = 0.05
            if(p_values[i]<signif_threshold){
              pval1 = "*"
              if(p_values[i] <signif_threshold/10){pval1 = "*"}
              y = max(unlist(groups[[i]]))
              y = y+1*b
              text(i, y+2*b, labels = pval1, cex = 1.7)
            }
          }
        }
      }
    }
    dev.off()
    
  }
  
}
Plot_rds_files_TCR(input_directory,output_directory,groups_PCA, batch,chain)


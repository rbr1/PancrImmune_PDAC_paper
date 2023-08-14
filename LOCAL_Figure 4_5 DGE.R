Draw_box_plot<-function(box,x,width,c,lwd,line_col){
  segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
  segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
  segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
  rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
  segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}

concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste0(res,v[i])}
	res
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha)) }


batch = "PDAC150K"
analysis_all = "PSEUDO_BULK_per_cell_type"
PLOTS = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Psuedo bulk/Psuedobulk per cell type/"
input_directory_groups  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"
output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/IMMUNOSURVEILLANCE/"
groups_PCA = readRDS(file = concat(c(input_directory_groups, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))
groups_PCA_Peng = readRDS(file = concat(c(input_directory_groups, "Immunogroups_2019_PUBLISHED_PDAC.rds")))
######################################
library(edgeR)
library(limma)
library(Glimma) ####
library(gplots)
library(RColorBrewer)

Perform_DGE_blood_versus_tumour_overall_PDAC150K<-function(batch, analysis_all, output_directory, input_directory, groups_PCA){
  list_matrices_gene_counts_sub= readRDS(file = concat(c(input_directory,"Pseudobulk_PDAC150K_immunosurveilling.rds")))
  all_cell_types = names(list_matrices_gene_counts_sub)
  ############ get list of cell type matrices which have enough samples per immunogroup
  list_matrices_gene_counts = NULL
  list_include_names = NULL
  all_cell_types_names = gsub(" ","_",gsub("+","p",gsub("-","_",gsub("/","_",all_cell_types, fixed = T), fixed = T), fixed = T), fixed = T)
  for(c in c(1:length(all_cell_types))){
    file = concat(c(input_directory,"Psuedobulk per cell type_PDAC150K/Raw_gene_counts_all_",all_cell_types_names[c],".txt"))
    if(file.exists(file)==T){
      p <- as.matrix(read.csv(file, head=T, sep="\t"))
      colnames(p) = gsub(".", "_",colnames(p),fixed = T)
      nz = which(colSums(p)!=0)
      n_group = c(length(which(gsub(".", "_",names(nz),fixed = T) %in% groups_PCA[[1]])), length(which(gsub(".", "_",names(nz),fixed = T) %in% groups_PCA[[2]])))
      if(min(n_group)>=3){
        list_matrices_gene_counts = c(list_matrices_gene_counts, list(p))
        list_include_names = c(list_include_names, all_cell_types[c])
        print (c)
      }
    }else{print(concat(c("FAIL ", all_cell_types_names[c])))}
  }
  names(list_matrices_gene_counts)=list_include_names
  
  ############# filter samples
  countdata = list_matrices_gene_counts[[1]]
  samples = colnames(countdata)
  genes = rownames(countdata)
  samples1 = gsub(".","_", samples,fixed = T)
  names(samples)= samples1
  
  ############# data table
  samples = colnames(countdata)
  genes = rownames(countdata)
  samples1 = gsub(".","_", samples,fixed = T)
  names(samples)= samples1
  
  group = rep(0,length(samples1))
  names(group) = samples1
  group[which(samples1 %in% c(groups_PCA[[1]], gsub("_biopsy","_blood",groups_PCA[[1]])))] = "ME"
  group[which(samples1 %in% c(groups_PCA[[2]], gsub("_biopsy","_blood",groups_PCA[[2]])))] = "AE"
  type = samples
  type[grep("blood",samples)] = "blood"
  type[grep("biopsy",samples)] = "tumour"
  
  sampleinfo = data.frame(samples)
  sampleinfo$group = group
  sampleinfo$type = type
  groups = rev(sort(unique(group)))
  
  ##### get mean expression per cell type
  cell_types = names(list_matrices_gene_counts)
  mean_expression_per_cell_type = matrix(data = 0, nrow = length(genes),ncol = length(cell_types), dimnames = c(list(genes), list(cell_types)))
  for(c in c(1:length(list_matrices_gene_counts))){
    countdata= list_matrices_gene_counts[[c]]
    countdata = countdata[,intersect(colnames(countdata), rownames(sampleinfo))]
    w = which(colSums(countdata)!=0)
    if(length(w)>=2){
      countdata = countdata[,w]
      mean_expression_per_cell_type[rownames(countdata),names(list_matrices_gene_counts)[c]] = apply(countdata, 1, sum)
    }}
  mean_expression_per_cell_type1 = mean_expression_per_cell_type[,which(colSums(mean_expression_per_cell_type)>100)]
  mean_expression_per_cell_type1 <- edgeR::cpm(mean_expression_per_cell_type1,log=F)
  
  ############# plot PCAs per cell type per group
  for(g in c(1:length(groups))){
    sampleinfo1 = sampleinfo[which(sampleinfo$group==groups[g]),]
    top_degs_list = NULL
    include = NULL
    for(c in c(1:length(list_matrices_gene_counts))){
      countdata= list_matrices_gene_counts[[c]][,rownames(sampleinfo1)]
      w = which(colSums(countdata)!=0)
      countdata = countdata[,w]
      genes = rownames(countdata)
      exclude= genes[c(grep("^TRAV",genes), grep("^TRBV",genes), grep("^IGHV",genes), grep("^IGKV",genes), grep("^IGLV",genes))]
      countdata = countdata[setdiff(genes, exclude),]
      sample_info2= sampleinfo1$type[w]
      if(length(unique(sampleinfo1$type[w]))>1){
        include = c(include, c)
        dge <- DGEList(counts = countdata, group = factor(sample_info2))
        keep <- filterByExpr(y = dge)
        dge <- dge[keep, , keep.lib.sizes=FALSE]
        dge <- calcNormFactors(object = dge)
        dge <- estimateDisp(y = dge)
        et <- exactTest(object = dge)
        top_degs = topTags(object = et, n = "Inf")
        summary(decideTests(object = et, lfc = 1))
        top_degs_list = c(top_degs_list, list(top_degs))
        print(c)
      }}
    names(top_degs_list) = names(list_matrices_gene_counts)[include]
    saveRDS(file = concat(c(output_directory, "DGE_raw_results_blood_vv_biopsy_PCAgroup",g,".rds")), top_degs_list)
  }
  
  genes = rownames(countdata)
  exclude= genes[c(grep("^TRAV",genes), grep("^TRBV",genes), grep("^IGHV",genes), grep("^IGKV",genes), grep("^IGLV",genes))]
  
  list_DGE_genes_p_value = NULL
  list_DGE_genes_FC = NULL
  for(g in c(1:length(groups))){
    top_degs_list = readRDS(file = concat(c(output_directory,"DGE_raw_results_blood_vv_biopsy_PCAgroup",g,".rds")))
    genes_use = setdiff(genes, exclude)
    m_DGE_genes_p_value = matrix(data = 2, nrow = length(genes_use),ncol = length(names(top_degs_list)), dimnames = c(list(genes_use), list(names(top_degs_list))))
    m_DGE_genes_FC = matrix(data = 0, nrow = length(genes_use),ncol = length(names(top_degs_list)), dimnames = c(list(genes_use), list(names(top_degs_list))))
    sig_threshold= 0.05#/length(top_degs_list)
    for(c in c(1:length(top_degs_list))){
      toptable = top_degs_list[[c]]
      toptable = toptable@.Data[[1]]
      w = which(toptable$FDR<sig_threshold)
      m_DGE_genes_p_value[rownames(toptable)[w],c] = toptable$FDR[w]
      m_DGE_genes_FC[rownames(toptable),c] = toptable$logFC
      print(c)
    }
    list_DGE_genes_p_value = c(list_DGE_genes_p_value, list(m_DGE_genes_p_value))
    list_DGE_genes_FC = c(list_DGE_genes_FC, list(m_DGE_genes_FC))
  }
  names(list_DGE_genes_p_value) = groups
  names(list_DGE_genes_FC) = groups
  
  #### plot heatmap of differences
  
  genes_of_interest = c(genes[c(grep("CCR",genes), grep("CXCR",genes),grep("ACKR",genes))])
  extra = c("CD69","PDCD1")#"IL2RA","IL1R2","SIGIRR","IL1RN", "IL36RN","CCRL2", "PITPNM3")
  # get genes that are expressed reasonably highly
  genes_of_interest = names(which(apply(mean_expression_per_cell_type1[genes_of_interest,], 1, max)>quantile(apply(mean_expression_per_cell_type1[genes_of_interest,], 1, max), 0.5)))
  genes_of_interest = c(sort(unique(genes_of_interest)), extra)
  #heatmap(mean_expression_per_cell_type1[genes_of_interest,])
  
  w = colnames(list_DGE_genes_FC[[1]])[grep("grouped",colnames(list_DGE_genes_FC[[1]]), invert = T)]
  w = setdiff(w,  w[c(grep("NK",w), grep("Myeloid",w))])
  comparison_FC = list_DGE_genes_FC[[1]][genes_of_interest,w]*0
  comparison_p = list_DGE_genes_FC[[1]][genes_of_interest,w]*0
  for(g in c(1:length(genes_of_interest))){
    x1 = list_DGE_genes_p_value[[1]][genes_of_interest[g],w]
    x2 = list_DGE_genes_p_value[[2]][genes_of_interest[g],w]
    y1 = list_DGE_genes_FC[[1]][genes_of_interest[g],w]
    y2 = list_DGE_genes_FC[[2]][genes_of_interest[g],w]
    
    w1 = intersect(which(x1<0.05), which(x2>0.05))
    w2 = intersect(which(x1>0.05), which(x2<0.05))
    w3 = intersect(which(x1<0.05), which(x2<0.05))
    comparison_p[genes_of_interest[g], w1] = 1
    comparison_p[genes_of_interest[g], w2] = 2
    comparison_p[genes_of_interest[g], w3] = 3
    comparison_FC[genes_of_interest[g], w1] = y1[w1]
    comparison_FC[genes_of_interest[g], w2] = y2[w2]
    comparison_FC[genes_of_interest[g], w3] = y1[w3]+y2[w3]
  }
  ### plot heatmap  
  matrix = comparison_p[which(rowSums(comparison_p)>0),]
  
  #matrix = (matrix[(sort(rownames(matrix))),])
  matrix = (matrix[,rev(sort(colnames(matrix)))])
  w = rev(c(grep("B cell", colnames(matrix)), grep("T cell", colnames(matrix))))
  matrix = matrix[,w]
  matrix_mean = comparison_FC[rownames(matrix),colnames(matrix)]
  
  matrix  = t(matrix)
  matrix_mean  = t(matrix_mean)
  mat_size = t(mean_expression_per_cell_type1[colnames(matrix),rownames(matrix)])
  main = "DGE per cell type - cytokines"
  m_col = matrix+2
  library(RColorBrewer)
  cols1 =  add.alpha (c(brewer.pal(8, "Dark2")[c(2,8)], brewer.pal(8, "Paired")[c(2)]), alpha = 0.5)
  cols1 = c("green", cols1,"red")
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)	
  
  ynames = rownames(matrix)
  xnames = colnames(matrix)
  ### plot heatmap
  x_range = c(1,length(xnames))
  y_range = c(1,length(ynames))
  cex = 0.9
  
  analysis = "Cytokine_receptors_total"
  fileout1=concat(c(output_directory,analysis,"_DGE_per_cell_type_", batch,"_cytokines.pdf"))
  w=6
  pdf(file=fileout1, height=w*1.0, width=w*1.65)
  par(mfrow= c(1,2), mar = c(5,12,6,5))
  
  plot(c(0.75,max(x_range)), c(0.75,(length(y_range))+0.5), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, ylim = c(1.5, max(y_range)-0),xlim = c(1,max(x_range)+0.5), axes = F)
  mtext(side = 3, text = gsub("-grouped1","",xnames), line = 0.01,cex=0.7,  at =c(1:length(xnames)), las = 2, font = 1)
  mtext(side = 2, text = gsub("-grouped1","",ynames), line = 0.01,cex=0.7,  at =c(1:length(ynames)), las = 1, font = 1)
  
  m_col = round(matrix, digits = 0)+1
  
  segments(0.5, c(1:length(ynames)) , max(x_range), c(1:length(ynames)), col = "grey", lty = 3, lwd = 0.5)
  segments(c(1:length(xnames)),1, c(1:length(xnames)) , max(y_range)+0.5,  col = "grey", lty = 3, lwd = 0.5)
  size = mat_size^0.25
  size = 1.8*size/max(size)
  pches = c(24, 25)
  for(i in c(1:length(m_col[,1]))){
    for(j in c(1:length(m_col[1,]))){
      if(matrix[i,j]!=0){
        pch = 21
        if(matrix_mean[i,j]>0){pch = pches[1]}
        if(matrix_mean[i,j]<0){pch = pches[2]}
        points(j,i, pch = pch, bg = cols[m_col[i,j]], col = cols[m_col[i,j]], cex = size[i,j])
      }}}
  
  legs = c("significant in ME only", "significant in AE only", "significant in ME and AE")
  plot(c(0,1), c(0,1), pch=20, col="white",xlab="",ylab ="",col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE)
  legend("topright", legs, pch = 21,cex= 0.8, bty="n", pt.bg = cols1[c(2:length(cols1))], col = NA, pt.cex = 1.2)
  legs = c("significantly up in tumour", "significantly down in tumour")
  legend("bottomright", legs, pch = pches,cex= 0.8, bty="n", pt.bg = "black", col = NA, pt.cex = 1.2)
  
  dev.off()
  
  ####### plot raw boxplots
  #### genes to plot: genes_of_interest
  #### cell types to plot: rownames(matrix)
  #### cell types to plot: rownames(matrix)
  cell_types_plot = rownames(matrix)
  groups = rev(sort(unique(sampleinfo$group)))
  types = (sort(unique(sampleinfo$type)))
  groups_exp = NULL
  genes_of_interest= colnames(matrix)
  for(c in c(1:length(cell_types_plot))){
    g2 = NULL
    nam = NULL
    for(g in c(1:length(groups))){
      sampleinfo1 = sampleinfo[which(sampleinfo$group==groups[g]),]
      top_degs_list = NULL
      include = NULL
      countdata= list_matrices_gene_counts[[cell_types_plot[c]]][,rownames(sampleinfo1)]
      w = which(colSums(countdata)!=0)
      countdata = countdata[,w]
      genes = rownames(countdata)
      sample_info2= sampleinfo1$type[w]
      lcpm <- edgeR::cpm(countdata,log=F)
      for(t in c(1:length(types))){
        w = which(sample_info2==types[t])
        g2 = c(g2, list(lcpm[genes_of_interest,w]))
        nam = c(nam, concat(c(groups[g]," ",types[t])))
      }
    }
    names(g2) = nam
    groups_exp = c(groups_exp, list(g2))
  }
  names(groups_exp) = cell_types_plot
  
  fileout1=concat(c(output_directory,analysis,"_DGE_per_cell_type_", batch,"_cytokines_boxplots_full.pdf"))
  w=3.2
  pdf(file=fileout1, height=w*5, width=w*9)
  par(mfrow= c(4,3), mar = c(13,5.5,3.5,0.5))
  summary_stats = NULL
  
  cell_types_plot = rownames(matrix)
  genes_of_interest = intersect(genes_of_interest, rownames(m_DGE_genes_p_value))
  p_values_sub = m_DGE_genes_p_value[genes_of_interest,cell_types_plot]
  
  for(f in c(1:length(genes_of_interest))){
    name = genes_of_interest[f]
    groups = NULL
    use = NULL
    for(c in c(1:length(groups_exp))){
      g1 = NULL
      exp = groups_exp[[c]]
      fail = 0
      for(g in c(1:length(exp))){
        if(length(exp[[g]])>length(genes_of_interest)){
          g1 = c(g1, list(exp[[g]][genes_of_interest[f],]))
        }else{fail = 1}}
      if(fail==0){
        names(g1) = names(exp)
        groups = c(groups, list(g1))
        use = c(use, c)
      }}
    names(groups) = names(groups_exp)[use]
    
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = names(groups)
    main = name
    max = max(c(unlist(groups), unlist(groups))*1.2)
    min = 0
    b = (max-min)*0.035
    ylab = concat(c(name, " expression (CPM)"))
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 26#length(groups)
    max_scale = min(c(max,max))
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range>200){scale = c(0:100)*50}
    if(range>500){scale = c(0:100)*100}
    if(range>1000){scale = c(0:100)*200}
    if(range>2000){scale = c(0:100)*500}
    if(range>5000){scale = c(0:100)*1000}
    if(range<=50){scale = c(0:100)*10}
    if(range<=30){scale = c(0:100)*5}
    if(range <15){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = name, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 3.5,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = gsub("activated","act.",factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.07
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = rev(sort(mean(shift)-shift))
    shift = shift*0.35/max(shift)
    
    cols11 =  brewer.pal(8, "Dark2")[c(2,8)]
    cols1 = c("gold",cols11[1], "darkblue",cols11[2])
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(i in c(1:l)){
      for(i1 in c(1:length(groups[[i]]))){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
      }
    }
    y = max(unlist(groups))
    x1 = NULL
    x2 = NULL
    if(name %in% rownames(list_DGE_genes_p_value[[1]])){x1 = list_DGE_genes_p_value[[1]][name,]}
    if(name %in% rownames(list_DGE_genes_p_value[[2]])){x2 = list_DGE_genes_p_value[[2]][name,]}

    
    for(i in c(1:length(groups))){	
      b = max*0.035
      signif_threshold = 0.05
      if(length(x1)!=0){
        pval = x1[names(groups)[i]]
        if(pval<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))+b
          segments(i-shift[1],y+0.5*b, i-shift[2],y+0.5*b,lwd = 3, col = cols11[1])
          text(mean(c(i-shift[c(1,2)])), y+2*b, labels = pval1, cex = 1.4,col = cols11[1])
        }
      }
      if(length(x1)!=0){
        pval = x2[names(groups)[i]]
        if(pval<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))+b
          segments(i-shift[3],y+0.5*b, i-shift[4],y+0.5*b,lwd = 3, col = cols11[2])
          text(mean(c(i-shift[c(3,4)])), y+2*b, labels = pval1, cex = 1.4,col = cols11[2])
        }
      }
    }
    
  }
  
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
  legend("topright", names(groups[[1]]), pch = 21,cex= 0.9, bty="n", pt.bg = cols1, col = cols1, text.col = cols1, text.font = 2)
  dev.off()
}
Perform_DGE_blood_versus_tumour_overall_PDAC150K(batch, analysis_all, output_directory, input_directory, groups_PCA_Peng)

batch = "PDAC150K"
Perform_DGE_blood_versus_tumour_immunosurveilling_PDAC150K<-function(batch, analysis_all, output_directory, input_directory, groups_PCA){
  list_matrices_gene_counts= readRDS(file = concat(c(input_directory,"Pseudobulk_PDAC150K_immunosurveilling.rds")))
  names(list_matrices_gene_counts)[which(names(list_matrices_gene_counts)=="Ductal")] = "Epithelial"
  
  
  ############# filter samples
  countdata = list_matrices_gene_counts[[1]]
  samples = colnames(countdata)
  genes = rownames(countdata)
  samples1 = gsub(".","_", samples,fixed = T)
  names(samples)= samples1
  
  group = rep(0,length(samples1))
  names(group) = samples1
  group[which(samples %in% groups_PCA[[1]])] = "ME"
  group[which(samples %in% groups_PCA[[2]])] = "AE"
  group[which(samples %in% gsub("biopsy","blood",groups_PCA[[1]]))] = "ME"
  group[which(samples %in% gsub("biopsy","blood",groups_PCA[[2]]))] = "AE"
  
  sampleinfo = data.frame(samples)
  sampleinfo$group = group
  sampleinfo = sampleinfo[which(sampleinfo$group!=0),]
  source = sampleinfo$samples
  source[grep("biopsy",sampleinfo$samples)] = "biopsy"
  source[grep("blood",sampleinfo$samples)] = "blood"
  sampleinfo$source = source
  sampleinfo$source1 = source
  sampleinfo$source1[which(sampleinfo$source1=="biopsy")] = "2.biopsy"
  sampleinfo$source1[which(sampleinfo$source1=="blood")] = "1.blood"
  ############# mean exp
  genes_all = NULL
  for(c in c(1:length(list_matrices_gene_counts))){
    genes_all= sort(unique(c(genes_all, rownames(list_matrices_gene_counts[[c]]))))}
  
  cell_types = names(list_matrices_gene_counts)
  mean_expression_per_cell_type = matrix(data = 0, nrow = length(genes_all),ncol = length(cell_types), dimnames = c(list(genes_all), list(cell_types)))
  for(c in c(1:length(list_matrices_gene_counts))){
    countdata= list_matrices_gene_counts[[c]]
    countdata = countdata[,intersect(colnames(countdata), rownames(sampleinfo))]
    w = which(colSums(countdata)!=0)
    if(length(w)>=2){
      countdata = countdata[,w]
      mean_expression_per_cell_type[rownames(countdata),names(list_matrices_gene_counts)[c]] = apply(countdata, 1, sum)
    }}
  mean_expression_per_cell_type1 = mean_expression_per_cell_type[,which(colSums(mean_expression_per_cell_type)>100)]
  mean_expression_per_cell_type1 <- edgeR::cpm(mean_expression_per_cell_type1,log=F)
  
  ############# plot PCAs per cell type
  groups_names = c("ME", "AE")
  for(gr in c(1:length(groups_names))){
    group_use = groups_names[gr]
    top_degs_list1 = NULL
    include = NULL
    for(c in c(1:length(names(list_matrices_gene_counts)))){
      countdata= list_matrices_gene_counts[[colnames(mean_expression_per_cell_type1)[c]]]
      if(length(countdata)>0){
        countdata = countdata[,intersect(colnames(countdata), rownames(sampleinfo))]
        w = intersect(names(which(colSums(countdata)!=0)), sampleinfo$samples[which(sampleinfo$group==group_use)])
        min = c(length(which(w %in% sampleinfo$samples[which(sampleinfo$source =="biopsy")])), 
                length(which(w %in% sampleinfo$samples[which(sampleinfo$source =="blood")])))
        if(min(min)>=3){
          countdata = countdata[,w]
          genes = rownames(countdata)
          sample_info1= sampleinfo[w,"source1"]
          if(length(unique(sample_info1))>1){
            print(colnames(mean_expression_per_cell_type1)[c])
            print(c)
            include = c(include, c)
            dge <- DGEList(counts = countdata, group = factor(sample_info1))
            keep <- filterByExpr(y = dge,min.count = 2, min.total.count = 5, large.n = 5, min.prop = 0.4)
            dge <- dge[keep, , keep.lib.sizes=FALSE]
            dge <- calcNormFactors(object = dge)
            dge <- estimateDisp(y = dge)
            et <- exactTest(object = dge)
            top_degs = topTags(object = et, n = "Inf")
            summary(decideTests(object = et, lfc = 1))
            top_degs_list1 = c(top_degs_list1, list(top_degs))
            print(c)
          }
        }
      }}
    names(top_degs_list1) = colnames(mean_expression_per_cell_type1)[include]
    saveRDS(file = concat(c(output_directory,"DGE_raw_results_PDAC150K_immunosurveilling_",group_use,"_level0.rds")), top_degs_list1)
  }
  
  ##### get meta p_values
  genes = NULL
  for(gr in c(1:length(groups_names))){
    group_use = groups_names[gr]
    top_degs_list1 = readRDS(file = concat(c(output_directory, "DGE_raw_results_PDAC150K_immunosurveilling_",group_use,"_level0.rds")))
    cell_types = sort(unique(c(names(top_degs_list1))))
    top_degs_lists = c(list(top_degs_list1))
    names(top_degs_lists) = c("PENG")
    for(s in c(1:length(top_degs_lists))){
      top_degs_list = top_degs_lists[[s]]
      for(c in c(1:length(top_degs_list))){
        toptable = top_degs_list[[c]][[1]]
        genes =  sort(unique(c(genes, rownames(toptable))))
        print(length(genes))
      }}}
  
  groups_names = c("ME", "AE")
  list_DGE_genes_p_value = NULL
  list_DGE_genes_FC = NULL
  for(gr in c(1:length(groups_names))){
    group_use = groups_names[gr]
    top_degs_list1 = readRDS(file = concat(c(output_directory, "DGE_raw_results_PDAC150K_immunosurveilling_",group_use,"_level0.rds")))
    cell_types = sort(unique(c(names(top_degs_list1))))
    top_degs_lists = c(list(top_degs_list1))
    names(top_degs_lists) = c("PENG")
    list_p_values = NULL
    list_FCs = NULL
    for(s in c(1:length(top_degs_lists))){
      top_degs_list = top_degs_lists[[s]]
      m_DGE_genes_p_value = matrix(data = 2, nrow = length(genes),ncol = length(cell_types), dimnames = c(list(genes), list(cell_types)))
      m_DGE_genes_FC = matrix(data = 0, nrow = length(genes),ncol = length(cell_types), dimnames = c(list(genes), list(cell_types)))
      sig_threshold= 0.05#/length(top_degs_list)
      for(c in c(1:length(top_degs_list))){
        toptable = top_degs_list[[c]][[1]]
        ct = names(top_degs_list)[c]
        w = which(toptable[,"PValue"]<sig_threshold)
        if(length(w)>0){
          m_DGE_genes_p_value[rownames(toptable)[w],ct] = toptable[w,"PValue"]
          m_DGE_genes_FC[rownames(toptable)[w],ct] = toptable[w,"logFC"]
        }
        print(c)
      }
      list_p_values = c(list_p_values, list(m_DGE_genes_p_value))
      list_FCs = c(list_FCs, list(m_DGE_genes_FC))
    }
    names(list_p_values) =   names(top_degs_lists)
    names(list_FCs) =   names(top_degs_lists)
    list_DGE_genes_p_value = c(list_DGE_genes_p_value, list(list_p_values))
    list_DGE_genes_FC = c(list_DGE_genes_FC, list(list_FCs))
  }
  names(list_DGE_genes_p_value) = groups_names
  names(list_DGE_genes_FC) = groups_names
  ### get common p_values
  
  library(scran)
  
  #### plot heatmap of differences
  
  genes_of_interest = c(genes[c(grep("CCR",genes), grep("CXCR",genes),grep("ACKR",genes))],"IL1R2","SIGIRR","IL1RN", "IL36RN","CCRL2", "PITPNM3")
  # get genes that are expressed reasonably highly
  genes_of_interest = names(which(rowSums(mean_expression_per_cell_type1[genes_of_interest,])>quantile(rowSums(mean_expression_per_cell_type1[genes_of_interest,]), 0.25)))
  genes_of_interest = sort(unique(genes_of_interest))
  #heatmap(mean_expression_per_cell_type1[genes_of_interest,])
  genes_of_interest = intersect(genes_of_interest, rownames(list_DGE_genes_p_value[[1]][[1]]))
  genes_of_interest = intersect(genes_of_interest, rownames(list_DGE_genes_p_value[[2]][[1]]))
  
  w1 = colnames(list_DGE_genes_FC[[1]][[1]])[grep("grouped",colnames(list_DGE_genes_FC[[1]][[1]]), invert = T)]
  w2 = colnames(list_DGE_genes_FC[[2]][[1]])[grep("grouped",colnames(list_DGE_genes_FC[[2]][[1]]), invert = T)]
  w = intersect(w1,w2)
  w = setdiff(w,  w[c(grep("NK",w), grep("Myeloid",w))])
  comparison_FC = list_DGE_genes_FC[[1]][[1]][genes_of_interest,w]*0
  comparison_p = list_DGE_genes_FC[[1]][[1]][genes_of_interest,w]*0
  for(g in c(1:length(genes_of_interest))){
    x1 = list_DGE_genes_p_value[[1]][[1]][genes_of_interest[g],w]
    x2 = list_DGE_genes_p_value[[2]][[1]][genes_of_interest[g],w]
    y1 = list_DGE_genes_FC[[1]][[1]][genes_of_interest[g],w]
    y2 = list_DGE_genes_FC[[2]][[1]][genes_of_interest[g],w]
    
    w1 = intersect(which(x1<0.05), which(x2>0.05))
    w2 = intersect(which(x1>0.05), which(x2<0.05))
    w3 = intersect(which(x1<0.05), which(x2<0.05))
    comparison_p[genes_of_interest[g], w1] = 1
    comparison_p[genes_of_interest[g], w2] = 2
    comparison_p[genes_of_interest[g], w3] = 3
    comparison_FC[genes_of_interest[g], w1] = y1[w1]
    comparison_FC[genes_of_interest[g], w2] = y2[w2]
    comparison_FC[genes_of_interest[g], w3] = y1[w3]+y2[w3]
  }
  ### plot heatmap  
  matrix = comparison_p[which(rowSums(comparison_p)>0),]
  
  matrix = (matrix[(sort(rownames(matrix))),])
  matrix = (matrix[,rev(sort(colnames(matrix)))])
  w = rev(c(grep("B cell", colnames(matrix)), grep("T cell", colnames(matrix))))
  matrix = matrix[,w]
  matrix_mean = comparison_FC[rownames(matrix),colnames(matrix)]
  
  matrix  = t(matrix)
  matrix_mean  = t(matrix_mean)
  mat_size = t(mean_expression_per_cell_type1[colnames(matrix),rownames(matrix)])
  main = "DGE per cell type - cytokines\nimmunosurveilling"
  m_col = matrix+2
  library(RColorBrewer)
  cols1 =  add.alpha (c(brewer.pal(8, "Dark2")[c(2,8)], brewer.pal(8, "Paired")[c(2)]), alpha = 0.5)
  cols1 = c("green", cols1,"red")
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)	
  
  ynames = rownames(matrix)
  xnames = colnames(matrix)
  ### plot heatmap
  x_range = c(1,length(xnames))
  y_range = c(1,length(ynames))
  cex = 0.9
  analysis = "Cytokine_receptors_immunosurveilling_only"
  fileout1=concat(c(output_directory,analysis,"_DGE_per_cell_type_", batch,"_cytokines.pdf"))
  w=6
  pdf(file=fileout1, height=w*0.8, width=w*1.6)
  par(mfrow= c(1,2), mar = c(5,12,6,5))
  
  plot(c(0.75,max(x_range)), c(0,(length(y_range))+0.5), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, ylim = c(0.5, max(y_range)-0),xlim = c(1,max(x_range)+0.5), axes = F)
  mtext(side = 3, text = gsub("-grouped1","",xnames), line = 0.01,cex=0.7,  at =c(1:length(xnames)), las = 2, font = 1)
  mtext(side = 2, text = gsub("-grouped1","",ynames), line = 0.01,cex=0.7,  at =c(1:length(ynames)), las = 1, font = 1)
  
  m_col = round(matrix, digits = 0)+1
  
  segments(0.5, c(1:length(ynames)) , max(x_range), c(1:length(ynames)), col = "grey", lty = 3, lwd = 0.5)
  segments(c(1:length(xnames)),1, c(1:length(xnames)) , max(y_range)+0.5,  col = "grey", lty = 3, lwd = 0.5)
  size = mat_size^0.25
  size = 1.8*size/max(size)
  pches = c(24, 25)
  for(i in c(1:length(m_col[,1]))){
    for(j in c(1:length(m_col[1,]))){
      if(matrix[i,j]!=0){
        pch = 21
        if(matrix_mean[i,j]>0){pch = pches[1]}
        if(matrix_mean[i,j]<0){pch = pches[2]}
        points(j,i, pch = pch, bg = cols[m_col[i,j]], col = cols[m_col[i,j]], cex = size[i,j])
      }}}
  
  legs = c("significant in ME only", "significant in AE only", "significant in ME and AE")
  plot(c(0,1), c(0,1), pch=20, col="white",xlab="",ylab ="",col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE)
  legend("topright", legs, pch = 21,cex= 0.8, bty="n", pt.bg = cols1[c(2:length(cols1))], col = NA, pt.cex = 1.2)
  legs = c("significantly up in tumour", "significantly down in tumour")
  legend("bottomright", legs, pch = pches,cex= 0.8, bty="n", pt.bg = "black", col = NA, pt.cex = 1.2)
  
  dev.off()
  
  ####### plot raw boxplots
  #### genes to plot: genes_of_interest
  
  #### cell types to plot: rownames(matrix)
  cell_types_plot = sort(rownames(matrix))
  groups = rev(sort(unique(sampleinfo$group)))
  types = rev(sort(unique(sampleinfo$source)))
  groups_exp = NULL
  for(c in c(1:length(cell_types_plot))){
    g2 = NULL
    nam = NULL
    for(g in c(1:length(groups))){
      sampleinfo1 = sampleinfo[which(sampleinfo$group==groups[g]),]
      top_degs_list = NULL
      include = NULL
      countdata= list_matrices_gene_counts[[cell_types_plot[c]]][,rownames(sampleinfo1)]
      w = which(colSums(countdata)!=0)
      countdata = countdata[,w]
      genes = rownames(countdata)
      sample_info2= sampleinfo1$source[w]
      lcpm <- edgeR::cpm(countdata,log=F)
      for(t in c(1:length(types))){
        w = which(sample_info2==types[t])
        g2 = c(g2, list(lcpm[genes_of_interest,w]))
        nam = c(nam, concat(c(groups[g]," ",types[t])))
      }
    }
    names(g2) = nam
    groups_exp = c(groups_exp, list(g2))
  }
  names(groups_exp) = cell_types_plot
  
  fileout1=concat(c(output_directory,analysis,"_DGE_per_cell_type_", batch,"_cytokines_boxplots_full.pdf"))
  w=3.2
  pdf(file=fileout1, height=w*5, width=w*8.5)
  par(mfrow= c(4,3), mar = c(13,5.5,3.5,0.5))
  summary_stats = NULL
  
  cell_types_plot = rownames(matrix)
  genes_of_interest = intersect(genes_of_interest, rownames(m_DGE_genes_p_value))
  p_values_sub = m_DGE_genes_p_value[genes_of_interest,cell_types_plot]
  
  for(f in c(1:length(genes_of_interest))){
    name = genes_of_interest[f]
    groups = NULL
    use = NULL
    for(c in c(1:length(groups_exp))){
      g1 = NULL
      exp = groups_exp[[c]]
      fail = 0
      for(g in c(1:length(exp))){
        if(length(exp[[g]])>length(genes_of_interest)){
          g1 = c(g1, list(exp[[g]][genes_of_interest[f],]))
        }else{fail = 1}}
      if(fail==0){
        names(g1) = names(exp)
        groups = c(groups, list(g1))
        use = c(use, c)
      }}
    names(groups) = names(groups_exp)[use]
    
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = names(groups)
    main = name
    max = max(c(unlist(groups), unlist(groups))*1.2)
    min = 0
    b = (max-min)*0.035
    ylab = concat(c(name, " expression (CPM)"))
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 16#length(groups)
    max_scale = min(c(max,max))
    range = max-min
    if(range>50){scale = c(0:100)*20}
    if(range>200){scale = c(0:100)*50}
    if(range>500){scale = c(0:100)*100}
    if(range>1000){scale = c(0:100)*200}
    if(range>2000){scale = c(0:100)*500}
    if(range>5000){scale = c(0:100)*1000}
    if(range<=50){scale = c(0:100)*10}
    if(range<=30){scale = c(0:100)*5}
    if(range <15){scale = c(0:100)*2.5}
    if(range <5){scale = c(0:100)*1}
    if(range <4){scale = c(0:100)*0.5}
    if(range <1.5){scale = c(0:1000)*0.2}
    if(range <0.5){scale = c(0:100)*0.1}
    if(range <0.1){scale = c(0:100)*0.01}
    if(range <0.01){scale = c(0:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = name, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 3.5,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.08
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = rev(sort(mean(shift)-shift))
    shift = shift*0.38/max(shift)
    
    cols11 =  brewer.pal(8, "Dark2")[c(2,8)]
    cols1 = c("gold",cols11[1], "darkblue",cols11[2])
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(i in c(1:l)){
      for(i1 in c(1:length(groups[[i]]))){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
      }
    }
    y = max(unlist(groups))
    x1 = NULL
    x2 = NULL
    if(name %in% rownames(list_DGE_genes_p_value[[1]][[1]])){x1 = list_DGE_genes_p_value[[1]][[1]][name,]}
    if(name %in% rownames(list_DGE_genes_p_value[[2]][[1]])){x2 = list_DGE_genes_p_value[[2]][[1]][name,]}

    for(i in c(1:length(groups))){	
      b = max*0.035
      signif_threshold = 0.05
      if(length(x1)!=0){
        pval = x1[names(groups)[i]]
        if(pval<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))+b
          segments(i-shift[1],y+0.5*b, i-shift[2],y+0.5*b,lwd = 3, col = cols11[1])
          text(mean(c(i-shift[c(1,2)])), y+2*b, labels = pval1, cex = 1.4,col = cols11[1])
        }
      }
      if(length(x1)!=0){
        pval = x2[names(groups)[i]]
        if(pval<signif_threshold){
          pval1 = "*"
          y = max(unlist(groups[[i]]))+b
          segments(i-shift[3],y+0.5*b, i-shift[4],y+0.5*b,lwd = 3, col = cols11[2])
          text(mean(c(i-shift[c(3,4)])), y+2*b, labels = pval1, cex = 1.4,col = cols11[2])
        }
      }
    }
    
    }

  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
  legend("topright", names(groups[[1]]), pch = 21,cex= 0.9, bty="n", pt.bg = cols1, col = cols1, text.col = cols1, text.font = 2)
  dev.off()
  
}
Perform_DGE_blood_versus_tumour_immunosurveilling_PDAC150K()

###
batch = "PENG"
Perform_DGE_Peng_only_PDAC_PCA_groups_Level1<-function(batch, analysis_all, output_directory, input_directory, groups_PCA_Peng){
  batch = "PENG"
  groups_PCA = groups_PCA_Peng
  
  names = c("all")
  names = c("all")
  list_matrices_gene_counts = NULL
  list_include_names = NULL
  for(f in c(1:length(names))){
    list_matrices_gene_counts_sub= readRDS(file = concat(c(input_directory,"Seurat_reannotated_PENG_",names[f],"_gene_counts_level2.rds")))
    inc = NULL
    all_cell_types = names(list_matrices_gene_counts_sub)
    for(c in c(1:length(list_matrices_gene_counts_sub))){
      m = list_matrices_gene_counts_sub[[c]]
      nz = which(colSums(m)!=0)
      n_group = c(length(which(gsub(".", "_",names(nz),fixed = T) %in% groups_PCA[[1]])), length(which(gsub(".", "_",names(nz),fixed = T) %in% groups_PCA[[2]])))
      if(min(n_group)>=3){
        list_matrices_gene_counts = c(list_matrices_gene_counts, list(m))
        list_include_names = c(list_include_names, all_cell_types[c])
        print (c)
      }else{print(concat(c("FAIL ", all_cell_types[c])))}
    }
  }
  names(list_matrices_gene_counts) =list_include_names
  names(list_matrices_gene_counts)[which(names(list_matrices_gene_counts)=="Ductal")] = "Epithelial"
  
  
  ############# filter samples
  countdata = list_matrices_gene_counts[[1]]
  samples = colnames(countdata)
  genes = rownames(countdata)
  samples1 = gsub(".","_", samples,fixed = T)
  names(samples)= samples1
  
  
  group = rep(0,length(samples1))
  names(group) = samples1
  group[grep("PDAC_TISSUE",samples1)] = "PDAC"
  group[grep("AdjNorm",samples1)] = "AdjNorm"
  group[grep("PDAC_T",samples1)] = "PDAC"
  group[grep("PDAC_N",samples1)] = "AdjNorm"
  
  sampleinfo = data.frame(samples)
  sampleinfo$group = group
  sampleinfo = sampleinfo[which(sampleinfo$group!=0),]
  source = sampleinfo$samples
  source[grep("PDAC_N",sampleinfo$samples)] = "PENG"
  source[grep("PDAC_T",sampleinfo$samples)] = "PENG"
  source[grep("_TISSUE_",sampleinfo$samples)] = "STEELE"
  sampleinfo$source = source
  
  ############# mean exp
  genes_all = NULL
  for(c in c(1:length(list_matrices_gene_counts))){
    genes_all= sort(unique(c(genes_all, rownames(list_matrices_gene_counts[[c]]))))}
  
  group_use = "PENG"
  cell_types = names(list_matrices_gene_counts)
  mean_expression_per_cell_type = matrix(data = 0, nrow = length(genes_all),ncol = length(cell_types), dimnames = c(list(genes_all), list(cell_types)))
  selected_genes = NULL
  selected_cell_types = NULL
  for(c in c(1:length(list_matrices_gene_counts))){
    countdata= list_matrices_gene_counts[[c]]
    countdata = countdata[,intersect(colnames(countdata), rownames(sampleinfo))]
    w = intersect(which(colSums(countdata)!=0), which(sampleinfo$source==group_use))
    if(length(w)>=2){
      countdata = countdata[,w]
      mean_expression_per_cell_type[rownames(countdata),names(list_matrices_gene_counts)[c]] = apply(countdata, 1, sum)
      countdata1 = edgeR::cpm(countdata)
      selected_genes = c(selected_genes, (countdata1["CCL1",]))
      selected_cell_types = c(selected_cell_types, rep(names(list_matrices_gene_counts)[c], length(countdata1["CCL1",])))
    }}
  mean_expression_per_cell_type1 = mean_expression_per_cell_type[,which(colSums(mean_expression_per_cell_type)>100)]
  mean_expression_per_cell_type1 <- edgeR::cpm(mean_expression_per_cell_type1,log=F)
  
  #### plot CCL1 expression
  library(tidyverse)
  theme_set(theme_bw(16))
  
  selected_index = c(1:length(selected_genes))
  df = data.frame(selected_index )
  df$CCL1_expression = selected_genes
  df$cell_type = as.factor(selected_cell_types)
  
  library(ggplot2)
  p <- ggplot(df, aes(x=cell_type, y=CCL1_expression)) + 
    geom_boxplot()+
    geom_jitter(width=0.15)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  
  fileout1=concat(c(output_directory,"Expression_levels", batch,"_CCL1_Peng_level1.pdf"))
  w=4
  pdf(file=fileout1, height=w*1.2, width=w*1.3)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  p
  dev.off()
  
  
  ############# plot PCAs per cell type
  
  group_use = "PENG"
  top_degs_list1 = NULL
  include = NULL
  for(c in c(1:length(colnames(mean_expression_per_cell_type1)))){
    countdata= list_matrices_gene_counts[[colnames(mean_expression_per_cell_type1)[c]]]
    countdata = countdata[,intersect(colnames(countdata), rownames(sampleinfo))]
    w = intersect(which(colSums(countdata)!=0), which(sampleinfo$source==group_use))
    if(length(w)>=6){
      countdata = countdata[,w]
      genes = rownames(countdata)
      exclude= genes[c(grep("^TRAV",genes), grep("^TRBV",genes), grep("^IGHV",genes), grep("^IGKV",genes), grep("^IGLV",genes))]
      countdata = countdata[setdiff(genes, exclude),]
      sample_info1= sampleinfo$group[w]
      
      if(length(unique(sample_info1))>1){
        include = c(include, c)
        dge <- DGEList(counts = countdata, group = factor(sample_info1))
        keep <- filterByExpr(y = dge,min.count = 2, min.total.count = 5, large.n = 5, min.prop = 0.4)
        dge <- dge[keep, , keep.lib.sizes=FALSE]
        dge <- calcNormFactors(object = dge)
        dge <- estimateDisp(y = dge)
        et <- exactTest(object = dge)
        top_degs = topTags(object = et, n = "Inf")
        summary(decideTests(object = et, lfc = 1))
        top_degs_list1 = c(top_degs_list1, list(top_degs))
        print(c)
      }
    }
  }
  names(top_degs_list1) = names(list_matrices_gene_counts)[include]
  saveRDS(file = concat(c(output_directory,"DGE_raw_results_PDAC_PCA_",group_use,"_PENG_level1.rds")), top_degs_list1)
  "NECTIN3" %in% rownames(top_degs)
  
  
  ##### get meta p_values
  top_degs_list1 = readRDS(file = concat(c(output_directory,"DGE_raw_results_PDAC_PCA_",group_use,"_PENG_level1.rds")))
  "NECTIN2" %in% rownames(top_degs_list1[[3]])
  cell_types = sort(unique(c(names(top_degs_list1))))
  
  top_degs_lists = c(list(top_degs_list1))
  names(top_degs_lists) = c("PENG")
  genes = NULL
  for(s in c(1:length(top_degs_lists))){
    top_degs_list = top_degs_lists[[s]]
    for(c in c(1:length(top_degs_list))){
      toptable = top_degs_list[[c]][[1]]
      genes = sort(unique(c(genes, rownames(toptable))))}}
  
  list_p_values = NULL
  list_FCs = NULL
  for(s in c(1:length(top_degs_lists))){
    top_degs_list = top_degs_lists[[s]]
    m_DGE_genes_p_value = matrix(data = 2, nrow = length(genes),ncol = length(cell_types), dimnames = c(list(genes), list(cell_types)))
    m_DGE_genes_FC = matrix(data = 0, nrow = length(genes),ncol = length(cell_types), dimnames = c(list(genes), list(cell_types)))
    sig_threshold= 0.05#/length(top_degs_list)
    for(c in c(1:length(top_degs_list))){
      toptable = top_degs_list[[c]][[1]]
      ct = names(top_degs_list)[c]
      w = which(toptable[,"PValue"]<sig_threshold)
      if(length(w)>0){
        m_DGE_genes_p_value[rownames(toptable)[w],ct] = toptable[w,"PValue"]
        m_DGE_genes_FC[rownames(toptable)[w],ct] = toptable[w,"logFC"]
      }
      print(c)
    }
    list_p_values = c(list_p_values, list(m_DGE_genes_p_value))
    list_FCs = c(list_FCs, list(m_DGE_genes_FC))
  }
  names(list_p_values) =   names(top_degs_lists)
  names(list_FCs) =   names(top_degs_lists)
  
  ### get common p_values
  
  library(scran)
  
  ### Plot_checkpoint_genes_only_boxplot
  w1 = which(sampleinfo$source=="PENG")
  w2 = which(sampleinfo$source=="STEELE")
  w3 = which(sampleinfo$group=="PDAC")
  w4 = which(sampleinfo$group=="AdjNorm")
  
  groups_PCA = c(list(sampleinfo$samples[intersect(w1,w3)]), list(sampleinfo$samples[intersect(w1,w4)]))
  library(RColorBrewer)
  cols1 =  rep(add.alpha (brewer.pal(8, "Dark2")[c(1,5)], alpha = 0.95), 3)
  cols =  add.alpha (cols1, alpha = 0.5)
  library(yarrr)
  cols1 =  rep(add.alpha (piratepal(palette = "google")[c(3:4)], alpha = 0.95), 3)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  matrix = list_p_values[[1]]
  genes_sub = rownames(matrix)
  genes_of_interest = (c("TIGIT","NECTIN2","NECTIN3","PVR","ICOS","ICOSLG","CTLA4","CD80", "CD86","CD28", "CCL1","CCL8","CCL16"))
  genes_of_interest = intersect(genes_of_interest, genes_sub)
  matrix_p_values_corrected = matrix[genes_of_interest,]
  
  #CTLA-4 is a highly endocytic molecule that binds to two distinct ligands, CD80 and CD86 
  # CD28/CTLA-4 receptor interactions stems from the fact that there are two natural ligands CD80 (B7-1) and CD86 (B7-2)
  # CCR8 ligands
  analysis   = "Peng"
  fileout1=concat(c(output_directory,"DGE_boxplots_PDAC_PCA_", batch,"_",analysis,"_Peng_level1.pdf"))
  w=3.2
  pdf(file=fileout1, height=w*0.48*6, width=w*3*1.6)
  par(mfrow= c(4,6), mar = c(7,3,2,0.2))
  p_val_grouped= NULL
  FC_grouped= NULL
  genes_grouped= NULL
  cell_type_grouped = NULL
  
  genes_of_interest = (c("TIGIT","NECTIN2","NECTIN3","PVR","ICOS","ICOSLG","CTLA4","CD80", "CD86","CD28", "CCL1","CCL8","CCL16"))
  genes_of_interest = intersect(genes_of_interest, genes_sub)
  
  pvals = list_p_values[[1]][genes_of_interest,cell_types]
  FC =  list_FCs[[1]][genes_of_interest,cell_types[c]]
  for(i in c(1:length(pvals[,1]))){
    pvals[i,] = p.adjust(pvals[i,], method = "fdr", n = length(pvals[i,]))}
  
  for(c in c(1:length(cell_types))){
    pval_merged = pvals[,cell_types[c]]
    FC =  list_FCs[[1]][,cell_types[c]]
    countdata= list_matrices_gene_counts[[cell_types[c]]]
    w = which(colSums(countdata)!=0)
    countdata = countdata[,w]
    sample_info1= factor(sampleinfo$group[w])
    lcpm <- edgeR::cpm(countdata,log=T)
    lcpm1 = lcpm[genes_of_interest,]
    
    #p_val_grouped= c(p_val_grouped, pval_merged[w])
    #genes_grouped= c(genes_grouped, genes_of_interest)
    #FC_grouped= c(FC_grouped, apply(FC, 1, mean)[w])
    #cell_type_grouped = c(cell_type_grouped, rep(cell_types[c], length(w)))
    groups = NULL
    p_value = NULL
    for(g in c(1:length(genes_of_interest))){
      g1 = NULL
      for(i in c(1:length(groups_PCA))){
        x = lcpm1[genes_of_interest[g],intersect(colnames(lcpm1), groups_PCA[[i]])]
        x = x[which(is.na(x)==F)]
        g1 = c(g1, list(x))
      }
      groups = c(groups, list(g1))
    }
    names(groups) = genes_of_interest
    pval_merged = pval_merged[genes_of_interest]
    FC = FC[genes_of_interest]
    
    if(length(unique(unlist(groups)))>=3){
      factors1 = c("tumour","adj. normal")
      factors = names(groups)
      main = concat(c(cell_types[c]))
      max = max(c(unlist(groups), unlist(groups))*1.45)
      min = 0
      min = min(c(unlist(groups), unlist(groups)))
      b = (max-min)*0.034
      min = min-b
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 12
      max_scale = min(c(max, 100))
      range = max-min
      if(range>55){scale = c(-100:100)*20}
      if(range<=55){scale = c(-100:100)*10}
      if(range<=40){scale = c(-100:100)*5}
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
          points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.7)
        }}
      
      for(i1 in c(1:length(pval_merged))){	
        b = max*0.035
        signif_threshold = 0.05
        if(pval_merged[i1]<signif_threshold){
          pval1 = "*"
          #if(pval_merged[i] <signif_threshold/10){pval1 = "**"}
          #if(pval_merged[i] <signif_threshold/100){pval1 = "***"}
          y = max(unlist(groups[[i1]]))
          y = y+1*b
          # segments(i-shift[1],y+b, i-shift[3],y+b,lwd = 3, col = "darkgrey")
          text(i1, y+2*b, labels = pval1, cex = 1.3)
        }
      }
    }
  }
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE, ylim = c(min, max))
  legend("topright", factors1, pch = 22,cex= 0.9, bty="n", pt.bg = add.alpha("white",alpha =0), col = add.alpha("white",alpha =0), text.col = cols1, text.font = 2)
  dev.off()
    
   
  ### Plot_checkpoint_genes_only_heatmap
    
  matrix = list_p_values[[1]]
  for(i in c(1:length(matrix[,1]))){
    matrix[i,] = p.adjust(matrix[i,], method = "fdr", n = length(matrix[i,]))}
  w1 = colnames(matrix)[which(colnames(matrix) %in% c("Acinar" ,"Centroacinar" ,"Endocrine","Endothelial.cell" , "Epithelial","Fibroblast" ,"Stellate.cell"))]
  w2 = setdiff(colnames(matrix), w1)
  matrix = matrix[,c(w1,w2)]
  
  genes_sub = rownames(matrix)
  genes_of_interest = (c("TIGIT","NECTIN2","NECTIN3","PVR","ICOS","ICOSLG","CTLA4","CD80", "CD86","CD28", "CCL1","CCL8","CCL16","CCL18"))
  genes_of_interest = intersect(genes_of_interest, genes_sub)
  matrix = matrix[genes_of_interest,]
  matrix_mean = list_FCs[[1]][genes_of_interest,colnames(matrix)]
  matrix_mean_exp = mean_expression_per_cell_type1[genes_of_interest,colnames(matrix)]
  
  matrix = t(matrix)
  matrix_mean = t(matrix_mean)
  matrix_mean_exp = t(matrix_mean_exp)
  
  main = "DGE per cell type - checkpoint"
  m_col = matrix+2
  library(RColorBrewer)
  cols1 =  add.alpha (c(brewer.pal(8, "Dark2")[c(2,8)]), alpha = 0.95)
  cols1 = c("red","darkgreen","lightblue", cols1,"red")
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)	
  
  ynames = rownames(matrix)
  xnames = colnames(matrix)
  ### plot heatmap
  x_range = c(1,length(xnames))
  y_range = c(1,length(ynames))
  cex = 0.9
  
  library(yarrr)
  cols1 =  rep(add.alpha (c(piratepal(palette = "google")[c(2,4)],"black"), alpha = 0.95), 3)
  cols =  add.alpha (cols1, alpha = 0.75)
  
  
  fileout1=concat(c(output_directory,"DGE_heatmap_PDAC_PCA_", batch,"_",analysis,"_Peng_level1.pdf"))
  w=5.8
  pdf(file=fileout1, height=w*1.2, width=w*1.5)
  par(mfrow= c(1,2), mar = c(5,7,19,5))
  
  plot(c(0.5,max(x_range)), c(0.5,(length(y_range))+0.5), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, ylim = c(0.5, max(y_range)+0.5), xlim = c(0.5,max(x_range)+5), axes = F)
  mtext(side = 3, text = xnames, line = 0.01,cex=0.7,  at =c(1:length(xnames)), las = 2, font = 1)
  mtext(side = 2, text = gsub("antigen-experienced","Ag-exp.",gsub("."," ",gsub("-grouped1","",ynames ), fixed = T)), line = 0.01,cex=0.7,  at =c(1:length(ynames)), las = 1, font = 1)
  
  m_col = round(matrix, digits = 0)+1
  
  segments(0.5, c(1:length(ynames)) , max(x_range), c(1:length(ynames)), col = "grey", lty = 3, lwd = 0.5)
  segments(c(1:length(xnames)),1, c(1:length(xnames)) , max(y_range)+0.5,  col = "grey", lty = 3, lwd = 0.5)
  size = matrix_mean_exp^0.25
  size = size*1.6/max(size)
  
  pches = c(21,21,21)
  for(i in c(1:length(m_col[,1]))){
    for(j in c(1:length(m_col[1,]))){
      ind = 3
      if(matrix[i,j]<0.05){
        if(matrix_mean[i,j]>0){ind = 1}
        if(matrix_mean[i,j]<0){ind = 2}}
      points(j,i, pch = pches[ind], bg = cols[ind], col = cols[ind], cex = size[i,j])
    }}
  
  legs = c("significantly higher in tumour than adj. normal tissue", "significantly higher in adj. normal than tumour","NS")
  order = c(1,3,2)
  plot(x_range, y_range, pch=20, col="white",xlab="",ylab ="",col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE)
  legend("topright", legs[order], pch = pches,cex= 0.4, bty="n", pt.bg = cols1[order], col = cols[order], pt.cex = 0.6)
  
  dev.off()
  
  }
Perform_DGE_Peng_only_PDAC_PCA_groups_Level1(batch, analysis_all, output_directory, input_directory, groups_PCA_Peng)


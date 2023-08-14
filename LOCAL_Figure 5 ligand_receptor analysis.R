Means_factor = function(factor, x){
	m = NULL
	for(i1 in c(1:length(levels(factor)))){
		x1 = x[which(factor==levels(factor)[i1])]
		x1 = x1[which(x1!=-1)]
		m = c(m, mean(x1))}
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

###############
batch = "PDAC150Ka"
analysis_all = "PSEUDO_BULK_per_cell_type"
PLOTS = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Psuedo bulk/Psuedobulk per cell type/"
input_directory_groups  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"
output_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/OUTPUTS/"
input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/RECEPTOR_LIGAND_DGE/"
groups_PCA = readRDS(file = concat(c(input_directory_groups, "PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))

############################## STEP 1 #################################################

Step_1_read_input_and_generate_cell_linkage_list<-function(input_directory, batch, output_directory){
  
  ###################### get cell counts
  analysis = "Chemokine_analysis"

  file = concat(c(input_directory, "Seurat_cell_counts_",batch,"_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[which(p[,"Group"]=="CD45"),]
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
  list_count = NULL
  
  include = NULL
  for(ind in c(1:length(cell_groups))){
    w = which(cell_group== cell_groups[ind])
    if(cell_groups[ind] =="CD45:broad"){
      w = intersect(w, which(p[,1] %in% c( "B cell"  , "Myeloid","NK" , "T cell"   )))
    }
    cell_type_sub = p[w,1]
    p1 = p[w,]
    m_count_norm = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(sample), dimnames = c(list(cell_type_sub), list(sample)))
    m_count = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(sample), dimnames = c(list(cell_type_sub), list(sample)))
    
    for(s in c(1:length(sample))){
      m_count_norm[, sample[s]]= as.numeric(p1[, sample[s]])*100/sum(as.numeric(p1[, sample[s]]))
      m_count[, sample[s]]= as.numeric(p1[, sample[s]])*100/sum(as.numeric(p1[, sample[s]]))
      if(sum(as.numeric(p1[, sample[s]]))<7){m_count_norm[, sample[s]] = -1}
    }
    if(length(which(apply(m_count_norm, 2, max)!=-1))>3){include = c(include, ind)}
    list_count = c(list_count, list(t(m_count[order(rownames(m_count)),])))
    m_count_norm= m_count_norm[order(rownames(m_count_norm)),]
    list_normalised= c(list_normalised, list(t(m_count_norm)))
  }
  names(list_normalised) = cell_groups
  names(list_count) = cell_groups
  list_normalised = list_normalised[include]
  
  m_count = list_normalised[[1]]
  cell_types_biopsy = names(which(colSums(m_count[grep("biopsy", rownames(m_count)),])>0))
  #################################
  type = "all"
  analysis = "Receptor_ligand_analysis"
  cell_types = c("all")
  list_counts = NULL
  list_expression = NULL
  
  file = concat(c(input_directory,"Cytokines_counts_ligand_receptor_genes_all_between_groups_",batch,".txt"))
  p <- as.matrix(read.csv(file, head=F, sep="\t"))
  headers = p[1,]
  #check
  headers[grep("T cell CD8 Senescent", headers)]
  
  headers = as.character(headers[c(3:(length(headers)-1))])
  headers1 = headers[grep(" .perc_expressing_per_cell_type", headers)]
  headers2 = headers[grep(" .total_RNA_per_cell_type", headers)]
  headers3 = headers[grep(" .prop_per_CD45_expressing", headers)]
  headers4 = headers[grep(" .n_expressing_per_cell_type", headers)]
  headers5 = headers[grep(" .raw_per_cell_type", headers)]
  
  headers_trim = gsub(" .perc_expressing_per_cell_type","", headers1)
  #check
  headers_trim2 = gsub(" .total_RNA_per_cell_type","", headers2)
  headers_trim3 = gsub(" .prop_per_CD45_expressing","", headers3)
  headers_trim4 = gsub(" .n_expressing_per_cell_type","", headers4)
  headers_trim5 = gsub(" .raw_per_cell_type","", headers5)
  headers_trim==headers_trim2
  headers_trim==headers_trim3
  headers_trim==headers_trim4
  headers_trim==headers_trim5
  #
  headers_trim_intersect = intersect(headers_trim,cell_types_biopsy)
  p = p[which(p[,2]=="pancreatic"),]
  colnames(p) = c("gene","site","pat", headers)
  gene = p[,1]
  site = p[,2]
  patient =p[,3]
  genes = sort(unique(gene))
  patients = sort(unique(patient))
  genes = sort(unique(gene))
  perc_expressing_per_cell_type = NULL
  total_RNA_per_cell_type = NULL
  prop_per_CD45_expressing = NULL
  n_expressing_per_cell_type = NULL
  raw_per_cell_type = NULL
  for(g in c(1:length(genes))){
    m1 = matrix(data = 0, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    m2 = matrix(data = -1, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    m3 = matrix(data = 0, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    m4 = matrix(data = 0, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    m5 = matrix(data = 0, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    w = which(gene== genes[g])
    for(i in c(1:length(headers_trim))){
      m1[patient[w], gsub(" .perc_expressing_per_cell_type","",headers1[i])] = as.numeric(p[w, headers1[i]])
      m2[patient[w],  gsub(" .total_RNA_per_cell_type","",headers2[i])] = as.numeric(p[w, headers2[i]])
      m3[patient[w],  gsub(" .prop_per_CD45_expressing","",headers3[i])] = as.numeric(p[w, headers3[i]])
      m4[patient[w],  gsub(" .n_expressing_per_cell_type","",headers4[i])] = as.numeric(p[w, headers4[i]])
      m5[patient[w],  gsub(" .raw_per_cell_type","",headers5[i])] = as.numeric(p[w, headers5[i]])
    }
    m1 = m1 [,which(headers_trim!="-")]
    m2 = m2 [,which(headers_trim!="-")]
    m3 = m3 [,which(headers_trim!="-")]
    m4 = m4 [,which(headers_trim!="-")]
    m5 = m5 [,which(headers_trim!="-")]
    cn= colnames(m1)
    #w = setdiff(c(1:length(cn)), c(grep("B cell", cn), grep("T cell", cn)))
    #cn [w] = paste("myeloid", cn [w])
    colnames(m1) = cn
    colnames(m2) = cn
    colnames(m3) = cn
    colnames(m4) = cn
    colnames(m5) = cn
    perc_expressing_per_cell_type = c(perc_expressing_per_cell_type, list(m1[,headers_trim_intersect]))
    total_RNA_per_cell_type = c(total_RNA_per_cell_type, list(m2[,headers_trim_intersect]))
    prop_per_CD45_expressing =c(prop_per_CD45_expressing,list(m3[,headers_trim_intersect]) )
    n_expressing_per_cell_type =c(n_expressing_per_cell_type,list(m4[,headers_trim_intersect]) )
    raw_per_cell_type =c(raw_per_cell_type,list(m5[,headers_trim_intersect]) )
  }
  names(perc_expressing_per_cell_type)= genes
  names(total_RNA_per_cell_type)= genes
  names(prop_per_CD45_expressing)= genes
  names(n_expressing_per_cell_type)= genes
  names(raw_per_cell_type)= genes
  
  saveRDS(file = concat(c(output_directory, "perc_expressing_per_cell_type_", batch,".rds")), perc_expressing_per_cell_type)
  
  ### get receptors and ligands
  file =  concat(c(input_directory, "Receptor_ligand_fantom.gsc.riken.jp.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  ligand = p[,"Ligand.ApprovedSymbol"]
  receptor = p[,"Receptor.ApprovedSymbol"]
  
  w = which(ligand %in% genes)
  cytokine_receptor = cbind(ligand, receptor)[w,]
  
  w = which(receptor %in% genes)
  cytokine_receptor = rbind(cytokine_receptor, cbind(ligand, receptor)[w,])
  cytokine_receptor = unique(cytokine_receptor)
  
  # make unique
  a1 = apply(cbind(cytokine_receptor), 1, paste, collapse = "|")
  a2 = apply(cbind(cytokine_receptor[,2], cytokine_receptor[,1]), 1, paste, collapse = "|")
  w = which(a1 %in% a2)
  length(w) ### all unique!
  
  # exclude pairs missing any gene
  w = intersect(which(cytokine_receptor[,1] %in% genes), which(cytokine_receptor[,2] %in% genes))
  ligand_receptor= cytokine_receptor[w,]
  ligands = sort(unique(ligand_receptor[,1]))
  receptors = sort(unique(ligand_receptor[,2]))
  
  # Print out receptor ligand list used
  out_file_table = concat(c(output_directory, "Receptor_ligand_list_used_PDAC150K.txt"))
  write.table(ligand_receptor, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")

  
  
  
  ###################### biopsy expression per gene
  ## generate matrix per rec-ligand pair for nxn
  ind1 = 1
  cell_types = colnames(perc_expressing_per_cell_type[[ligand_receptor[ind1,1]]])
  all_cell_types_combined = NULL
  all_cell_types_from = NULL
  all_cell_types_to = NULL
  for(c1 in c(1:length(cell_types))){
    for(c2 in c(1:length(cell_types))){
      all_cell_types_combined = c(all_cell_types_combined, concat(c(cell_types[c1]," - ", cell_types[c2])))
      all_cell_types_from = c(all_cell_types_from, cell_types[c1])
      all_cell_types_to = c(all_cell_types_to, cell_types[c2])
    }
  }
  
  rec_ligand_list = NULL
  list_names = NULL
  for(ind1 in c(1:length(ligand_receptor[,1]))){
    lig = perc_expressing_per_cell_type[[ligand_receptor[ind1,1]]]
    rec = perc_expressing_per_cell_type[[ligand_receptor[ind1,2]]]
    number_lig = raw_per_cell_type[[ligand_receptor[ind1,1]]]
    number_rec = raw_per_cell_type[[ligand_receptor[ind1,2]]]
    m_links_type = matrix(data = 0, nrow = length(all_cell_types_combined), ncol = length(patients), dimnames = c(list(all_cell_types_combined), list(patients)))
    w1 = intersect(which(apply(rec, 1, function(x){sum(x[which(x!=-1)])})>0), which(apply(lig, 1, function(x){sum(x[which(x!=-1)])})>0))
    if(length(w1)>0){
      for(i1 in c(1:length(w1))){
        a = lig[patients[w1[i1]],]
        b = rec[patients[w1[i1]],]
        n = number_lig[w1[i1],]
        wn = names(n)[intersect(which(n>=3), intersect(which(a>=0), which(b>=0)))]
        if(length(wn)>0){
          for(j1 in c(1:length(wn))){
            for(j2 in c(1:length(wn))){
              nam = concat(c(wn[j1], " - ", wn[j2]))
              m_links_type[nam,patients[w1[i1]]] = a[wn[j1]]*b[wn[j2]]
            }}}}}
    rec_ligand_list = c(rec_ligand_list, list(m_links_type))
    list_names = c(list_names, concat(c(ligand_receptor[ind1,1], "_", ligand_receptor[ind1,2] )))
    print (ind1)
  }
  names(rec_ligand_list) = list_names
  
  saveRDS(file = concat(c(output_directory,"Receptor_ligand_analysis_biopsy_per_pair_", batch,".rds")), rec_ligand_list)
}

Step_1_read_input_and_generate_cell_linkage_list(input_directory, batch, output_directory)

############################## STEP 2 #################################################
rec_ligand_list = readRDS(file = concat(c(output_directory,"Receptor_ligand_analysis_biopsy_per_pair_", batch,".rds")))
list_names = names(rec_ligand_list)
library(RColorBrewer)

perc_expressing_per_cell_type = readRDS(file = concat(c(output_directory, "perc_expressing_per_cell_type_", batch,".rds")))

cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
cols =  add.alpha (cols1, alpha = 0.5)
cols2 =  add.alpha (cols, alpha = 0.5)	
group_PCA_list = NULL
factor_PCA = NULL
for(i in c(1:length(groups_PCA))){
	group_PCA_list = c(group_PCA_list, gsub("_biopsy","",groups_PCA[[i]]))
	factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
}

Plot_diff_network_level1_annotations<-function(input_directory, batch, output_directory, rec_ligand_list, list_names, cols1, cols, cols2, group_PCA_list, factor_PCA){
  file = concat(c(input_directory, "Seurat_annotation_matching_PDAC150Ka.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  cell_type_level0 =trimws(p[,1], "r")
  cell_type_level1 = apply(cbind(trimws(p[,"Level.1"], "r"), ""), 1, paste, collapse = "")
  cell_type_level2 = apply(cbind(trimws(p[,"Level.2"], "r"), ""), 1, paste, collapse = "")
  names(cell_type_level1) = cell_type_level0
  names(cell_type_level2) = cell_type_level0
  
  cells_type_level0 = sort(unique(cell_type_level0))
  cells_type_level1 = sort(unique(cell_type_level1))
  cells_type_level2 = sort(unique(cell_type_level2))
  
  cells_type_level1_unique = setdiff(cells_type_level1, cells_type_level0)
  cells_type_level2_unique = setdiff(setdiff(cells_type_level2, cells_type_level0), cells_type_level1)
  
  all_cell_types = (unique(c(cells_type_level0, cells_type_level1_unique, cells_type_level2_unique)))
  
  names(cell_type_level1) = gsub("B cell plasm","B cell  plasm",names(cell_type_level1))
  cell_types[which(cell_types %in% names(cell_type_level1)==F)]
  
  Very.broad.annotation = cell_type_level1
  ######
  sites = names(rec_ligand_list)
  
  cell_types = sort(unique(colnames(perc_expressing_per_cell_type[[ligand_receptor[1,1]]])))
  
  cell_combined_all = rownames(rec_ligand_list[[1]])# sort(unique(p[,"cell-combination"]))
  str = strsplit(cell_combined_all," - ", fixed = T)
  new_cell_combined = NULL
  exc = NULL
  for(i in c(1:length(str))){
    if(str[[i]][1] %in% names(Very.broad.annotation)==F){
      exc = c(exc, i)
    }
    new_cell_combined = c(new_cell_combined , concat(c(Very.broad.annotation[str[[i]][1]], " - ",Very.broad.annotation[str[[i]][2]])))}
  
  
  names(new_cell_combined) = cell_combined_all
  cell_combined_all = unique(new_cell_combined)
  n_cell_combined_all = rep(0,length(cell_combined_all))
  names(n_cell_combined_all) = cell_combined_all
  list_linkages = matrix(data = 0, nrow = length(cell_combined_all), ncol = length(patients), dimnames = c(list(cell_combined_all), list(patients)))
  for(s in c(1:length(patients))){
    for(ind1 in c(1:length(sites))){
      edge = rec_ligand_list[[ind1]][,patients[s]]
      edge = edge[which(edge!=0)]
      if(length(edge)>0){
        cell_combined = names(edge)
        cell_combined1= new_cell_combined[cell_combined]
        t = table(cell_combined1)
        list_linkages[names(t),patients[s]] = list_linkages[names(t),patients[s]]+1
      }}
    print(s)
  }

  list_linkages1 = list_linkages[which(rowSums(list_linkages)!=0),]
  mat_stat= t(list_linkages1[,group_PCA_list])
  factor = factor(factor_PCA)
  fit = manova(formula = mat_stat ~ factor)
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
  nam = colnames(mat_stat)
  p_value[which(is.na(p_value))] = 2
  names(p_value) = nam
  signif_group1 = nam[intersect(which(p_value<0.05), which(means[,1]>means[,2]))]
  signif_group2 = nam[intersect(which(p_value<0.05), which(means[,1]<means[,2]))]
  x = cbind("receptor-ligand numbers between AE and ME patients",nam, p_value, means)
  colnames(x) = c("Analysis","cell group", "p_value", "mean_ME","mean_AE")
  
  fileout = concat(c(output_directory,"Stats_R_L_numbers_between_cell_type_pairs_PCA_groups.txt"))
  write.table(x, file = fileout, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  

  file = concat(c(output_directory,"Stats_R_L_numbers_between_cell_type_pairs_PCA_groups.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
 
  library(igraph)
  
  n_cell_combined_all = rep(0,length(new_cell_combined))
  names(n_cell_combined_all) = new_cell_combined
  mean_groups = cbind(as.numeric(p[,"mean_ME"]), as.numeric(p[,"mean_AE"]))
  mean_groups_all = apply(mean_groups, 1, mean)
  names(mean_groups_all) = p[,"cell.group"]

  cell_combined = names(mean_groups_all)[which(mean_groups_all>0)]
  cell_from = strsplit(cell_combined, " - ", fixed = T)
  cell_to = NULL
  for(i in c(1:length(cell_from))){
    cell_to = c(cell_to, cell_from[[i]][2])
    cell_from[i]= cell_from[[i]][1]
  }
  cell_from = unlist(cell_from)
  all_cell_types = sort(unique(c(cell_from, cell_to)))
  edge_strength = mean_groups_all
  from <- cell_from
  to <- cell_to
  g <- graph.empty( n=0, directed=T)
  g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
  names <- V(g)$name
  ids <- 1:length(names)
  names(ids) <- names
  edges <- matrix(c(ids[from], ids[to]), nc=2)
  g <- add.edges(g, t(edges), weight= edge_strength)		
  
  d1 = degree(g, v = V(g), mode = c( "out"),loops = TRUE, normalized = FALSE)^2
  library(RColorBrewer)
  cols=add.alpha(colorRampPalette(c('white','yellow','orange', "red"))(100), alpha = 0.8)
  m1 = d1*100/max(d1)
  m1 = round(m1, digits = 0)
  range(m1)
  sizes = 1
  sizes_scaled = sizes^0.3
  sizes_scaled = sizes_scaled*5
  
  V(g)$size<-sizes_scaled
  V(g)$label.cex<-0.5
  V(g)$color = add.alpha(cols[m1],alpha = 0.25)
  layout1 =layout_in_circle(g)
  layout1 = layout_with_graphopt(g,niter = 500)
  edge_strength_plot = edge_strength ^1
  edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05

  Ecolor = rep(add.alpha("grey",alpha = 0.5), length(edge_strength_plot))
  names(Ecolor) = names(edge_strength_plot)
  Ecolor[signif_group1] = add.alpha("red",alpha = 0.5)
  Ecolor[signif_group2] = add.alpha("blue",alpha = 0.5)
  Ecolor=Ecolor[names(edge_strength_plot)]
  E(g)$color = Ecolor

  layout1 =layout_in_circle(g)
  split_point = 4
  l = length(layout1[,1])
  layout1 = layout1[c(split_point:l, 1:(split_point-1)),]
  layout1 = layout_with_graphopt(g,start = layout1 , niter = 300,charge = 0.00010,mass = 30,spring.length = 200,spring.constant = 1.2,max.sa.movement = 4)
  
  fileout1=concat(c(output_directory, "Receptor_ligand_analysis_network_biopsy_pca_1vs2_", batch,"_simplified.pdf"))
  w=15
  pdf(file=fileout1, height=w*1, width=w*1)
  par(mfrow= c(1,1), mar = c(1,1,3,1))
  V(g)$name = gsub("myeloid ","",V(g)$name)
  V(g)$name = gsub("T cell ","", V(g)$name)
  plot(g, layout=layout1, edge.color= E(g)$color, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 1.9,vertex.label.font = 2, edge.arrow.size = 1, edge.lty = 1, main = concat(c("tumour")),xlim = c(-1, 1)*1.25,ylim = c(-1, 1),cex.main =3)
  
  g1= g
  E(g1)$color[which(E(g1)$color ==add.alpha("red",alpha = 0.5))] = add.alpha("grey",alpha = 0.5)
  plot(g1, layout=layout1, edge.color= E(g1)$color, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 1.9,vertex.label.font = 2, edge.arrow.size = 1, edge.lty = 1, main = concat(c("tumour group 1")),xlim = c(-1, 1)*1.25,ylim = c(-1, 1),cex.main =3)

  g1= g
  E(g1)$color[which(E(g1)$color ==add.alpha("blue",alpha = 0.5))] = add.alpha("grey",alpha = 0.5)
  plot(g1, layout=layout1, edge.color= E(g1)$color, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 1.9,vertex.label.font = 2, edge.arrow.size = 1, edge.lty = 1, main = concat(c("tumour group 2")),xlim = c(-1, 1)*1.25,ylim = c(-1, 1),cex.main =3)
  
  dev.off()	
  
  
}
Plot_diff_network_level1_annotations(input_directory, batch, output_directory, rec_ligand_list, list_names, cols1, cols, cols2, group_PCA_list, factor_PCA)

### Not used in paper
Plot_get_in_out_degree<-function(input_directory, batch, output_directory,rec_ligand_list, list_names, cols1, cols, cols2, group_PCA_list, factor_PCA){
  sites = names(rec_ligand_list)
  cell_types = sort(unique(colnames(perc_expressing_per_cell_type[[ligand_receptor[1,1]]])))
  
  mean_pca1_in = matrix(data = -1,nrow = length(cell_types), ncol = length(sites), dimnames = c(list(cell_types), list(sites)))
  mean_pca2_in = matrix(data = -1,nrow = length(cell_types), ncol = length(sites), dimnames = c(list(cell_types), list(sites)))
  mean_p_value_in = matrix(data = 2,nrow = length(cell_types), ncol = length(sites), dimnames = c(list(cell_types), list(sites)))
  
  mean_pca1_out = matrix(data = -1,nrow = length(cell_types), ncol = length(sites), dimnames = c(list(cell_types), list(sites)))
  mean_pca2_out = matrix(data = -1,nrow = length(cell_types), ncol = length(sites), dimnames = c(list(cell_types), list(sites)))
  mean_p_value_out = matrix(data = 2,nrow = length(cell_types), ncol = length(sites), dimnames = c(list(cell_types), list(sites)))
  
  for(ind1 in c(1:length(sites))){
    m_out = matrix(data = 0,nrow = length(cell_types), ncol = length(patients), dimnames = c(list(cell_types), list(patients)))
    m_in = matrix(data = 0,nrow = length(cell_types), ncol = length(patients), dimnames = c(list(cell_types), list(patients)))
    for(s in c(1:length(patients))){
      edge = rec_ligand_list[[ind1]][,patients[s]]
      edge = edge[which(edge!=0)]
      if(length(edge)>0){
        cell_combined = names(edge)
        cell_from = strsplit(cell_combined, " - ", fixed = T)
        cell_to = NULL
        for(i in c(1:length(cell_from))){
          cell_to = c(cell_to, cell_from[[i]][2])
          cell_from[i]= cell_from[[i]][1]
        }
        cell_from = unlist(cell_from)
        all_cell_types = sort(unique(c(cell_from, cell_to)))
        library(igraph)
        g <- graph.empty( n=0, directed=T)	
        g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
        names <- V(g)$name
        ids <- 1:length(names)
        names(ids) <- names
        from <- cell_from
        to <- cell_to
        edges <- matrix(c(ids[from], ids[to]), nc=2)
        g <- add.edges(g, t(edges), weight= edge)
        d1 = degree(g, v = V(g), mode = c( "out"),loops = TRUE, normalized = FALSE)
        d2 = degree(g, v = V(g), mode = c( "in"),loops = TRUE, normalized = FALSE)
        m_out[names(d1), patients[s] ] = d1*100/length(edge)
        m_in[names(d2), patients[s] ] = d2*100/length(edge)
      }}
    m_out1 = m_out[which(apply(m_out, 1, function(x){length(which(x>0))})>=3), ]
    m_in1 = m_in[which(apply(m_in, 1, function(x){length(which(x>0))})>=3), ]
    if(length(m_out1)>length(patients)*1.8){
      mat_stat= t(m_out1[,group_PCA_list])
      factor = factor(factor_PCA)
      fit = manova(formula = mat_stat ~ factor)
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
      nam = colnames(mat_stat)
      p_value[which(is.na(p_value))] = 2
      names(p_value) = nam
      # head(sort(p_value),5)
      mean_pca1_out[nam, sites[ind1]] = means[,1]
      mean_pca2_out[nam, sites[ind1]] = means[,2]
      mean_p_value_out[nam, sites[ind1]] = p_value
    }
    if(length(m_in1)>length(patients)*1.8){
      mat_stat= t(m_in1[,group_PCA_list])
      factor = factor(factor_PCA)
      fit = manova(formula = mat_stat ~ factor)
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
      nam = colnames(mat_stat)
      p_value[which(is.na(p_value))] = 2
      # names(p_value) = nam
      # head(sort(p_value),5)
      mean_pca1_in[nam, sites[ind1]] = means[,1]
      mean_pca2_in[nam, sites[ind1]] = means[,2]
      mean_p_value_in[nam, sites[ind1]] = p_value
    }
    print (ind1)
  }
  list_in_out = c(list(c(list(mean_pca1_in), list(mean_pca2_in), list(mean_p_value_in)) ), list(c(list(mean_pca1_out), list(mean_pca2_out), list(mean_p_value_out)) ))
  names(list_in_out) = c("in","out")
  
  #### plot differential frequency plot
  summary_stats = NULL
  fileout1=concat(c(output_directory,"Receptor_ligand_analysis_centrality_correlation_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*1.6, width=w*1.95*2)
  par(mfrow= c(1,2), mar = c(5,5,2.5,15))
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(3,8, 2)], alpha = 0.75)
  for(s in c(1:length(list_in_out))){
    mean_pca1 = list_in_out[[s]][[1]]
    mean_pca2 = list_in_out[[s]][[2]]
    pval_mat = list_in_out[[s]][[3]]
    x = NULL
    y = NULL
    names = NULL
    pval = NULL
    for(i in c(1:length(mean_pca1[1,]))){
      x = c(x, mean_pca1[,i])
      y = c(y, mean_pca2[,i])
      pval = c(pval, pval_mat[,i])
      names = c(names, apply(cbind(rownames(mean_pca1), colnames(mean_pca1)[i]), 1, paste, collapse = " "))
    }
    names(x) = names
    names(y) = names
    names(pval) = names
    b = max(y)*0.1
    range = max(abs(c(x,y)))*1.1
    cols = rep(1, length(pval))
    w = which(pval <0.05)
    cols[intersect(w, which(x <y))] = 2
    cols[intersect(w, which(x >y))] = 3
    plot(c(0, range), c(0, range), xlim = c(0, range),ylim = c(0, range), pch =21, bg = add.alpha("white", alpha = 0), col = add.alpha("white", alpha = 0), cex = 1.2, , xlab = "mean (group 1)", ylab = "mean (group 2)", main = names(list_in_out)[s])
    w = order(cols)
    points(x[w],y[w], pch =21, bg = cols1[cols[w]], col = add.alpha("white", alpha = 0), cex = 0.8)
    w = names(which(pval<0.05))
    w = names(head(rev(sort(abs(x[w]-y[w]))), 60))	
    lab = gsub("_"," ",w)
    order = w[order(y[w])]
    cex =0.9
    seq =(seq(from = 0,to = max(range), length = length(order)))
    mtext(side = 4, text = order, line = 0.2,cex= cex-0.5, at = seq, las = 1, font = 1)
    for(i in c(1:length(order))){
      segments(x[order[i]], y[order[i]], max(range)+b*0.2, seq[i], lwd = 1, lty = 1, col = add.alpha ("grey", alpha = 0.5))
    }
  }
  dev.off()
  
}
Plot_get_in_out_degree(input_directory, batch, output_directory,rec_ligand_list, list_names, cols1, cols, cols2, group_PCA_list, factor_PCA)


Get_network_centrality_between_patient_groups<-function(input_directory, batch, output_directory,rec_ligand_list, list_names, cols1, cols, cols2, group_PCA_list, factor_PCA){
  sites = names(rec_ligand_list)
  cell_types = sort(unique(colnames(perc_expressing_per_cell_type[[ligand_receptor[1,1]]])))
  m_out = matrix(data = -1,nrow = length(cell_types), ncol = length(patients), dimnames = c(list(cell_types), list(patients)))
  m_in = matrix(data = -1,nrow = length(cell_types), ncol = length(patients), dimnames = c(list(cell_types), list(patients)))
  ### self-feedback loops
  m_self = matrix(data = -1,nrow = length(cell_types), ncol = length(patients), dimnames = c(list(cell_types), list(patients)))
  
  for(s in c(1:length(patients))){
    print (s)
    library(igraph)
    g <- graph.empty( n=0, directed=T)	
    g <- igraph::add.vertices(g, length(cell_types), name= cell_types) 
    names <- V(g)$name
    ids <- 1:length(names)
    names(ids) <- names
    for(ind1 in c(1:length(sites))){
      edge = rec_ligand_list[[ind1]][,patients[s]]
      edge = edge[which(edge!=0)]
      if(length(edge)>0){
        cell_combined = names(edge)
        cell_from = strsplit(cell_combined, " - ", fixed = T)
        cell_to = NULL
        for(i in c(1:length(cell_from))){
          cell_to = c(cell_to, cell_from[[i]][2])
          cell_from[i]= cell_from[[i]][1]
        }
        cell_from = unlist(cell_from)
        all_cell_types = sort(unique(c(cell_from, cell_to)))
        from <- cell_from
        to <- cell_to
        edges <- matrix(c(ids[from], ids[to]), nc=2)
        g <- add.edges(g, t(edges), weight= edge)
      }}
    d1 = degree(g, v = V(g), mode = c( "out"),loops = TRUE, normalized = FALSE)
    d2 = degree(g, v = V(g), mode = c( "in"),loops = TRUE, normalized = FALSE)
    m_out[names(d1), patients[s] ] = d1
    m_in[names(d2), patients[s] ] = d2
    el <- as_edgelist(g)
    el = el[which(el[,1]==el[,2]), ]
    self_loops = table(el)
    m_self[names(self_loops), patients[s]] = self_loops
  }
  
  lists_degree = c(list(m_in), list(m_out))
  names(lists_degree) = c("in","out")
  
  saveRDS(file = concat(c(output_directory,"Receptor_ligand_analysis_m_self_pca_per_cell_type_", batch,".rds")),m_self)
  
  headers = c("in","out")
  mean_pca1 = matrix(data = -1,nrow = length(cell_types), ncol = length(headers), dimnames = c(list(cell_types), list(headers)))
  mean_pca2 = matrix(data = -1,nrow = length(cell_types), ncol = length(headers), dimnames = c(list(cell_types), list(headers)))
  mean_p_value = matrix(data = 2,nrow = length(cell_types), ncol = length(headers), dimnames = c(list(cell_types), list(headers)))
  
  for(ind in c(1:length(lists_degree))){
    mat_stat= t(lists_degree[[ind]][,group_PCA_list])
    w = which(factor_PCA ==1)
    mean_pca1[colnames(mat_stat), ind] = colSums(mat_stat [w, ])/length(w)
    w = which(factor_PCA ==2)
    mean_pca2[colnames(mat_stat), ind] = colSums(mat_stat [w, ])/length(w)
    w = which(apply(mat_stat, 2, function(x){length(unique(x))})>=4)
    mat_stat =(mat_stat[,w])
    # for(i in c(1:length(mat_stat[,1]))){
    # x = mat_stat[i,]
    # x[which(x==0)] = NA
    # mat_stat[i,] = x
    # }
    # w1 = which(apply(mat_stat[which(factor_PCA ==1),], 2, function(x){length(which(is.na(x)==F))})>=3)
    # w2 = which(apply(mat_stat[which(factor_PCA ==2),], 2, function(x){length(which(is.na(x)==F))})>=3)
    # mat_stat  = mat_stat [,intersect(w1,w2)]
    
    factor = factor(factor_PCA)
    # fit = manova(formula = log10(mat_stat+1) ~ factor)
    fit = manova(formula = mat_stat^1 ~ factor)
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
    nam = colnames(mat_stat)
    p_value[which(is.na(p_value))] = 2
    names(p_value) = nam
    mean_p_value[nam, ind] = p_value
    #head(sort(p_value),20)
    #p_value["T cell CD4 Treg Activated"]
    length(which(p_value<0.05))
    
  }
  length(which(mean_p_value[,1]<0.05))
  
  list_in_out = c(list(mean_p_value), list(mean_pca1), list(mean_pca2))
  names(list_in_out) = c("p-value","group1","group2")
  saveRDS(file = concat(c(output_directory,"Receptor_ligand_analysis_in_out_degree_pca_per_cell_type_", batch,".rds")),list_in_out)
  
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)	
  use = names(which(rowSums(cbind(mean_pca1, mean_pca2))!=0))
  mean_pca1 = mean_pca1[use,]
  mean_pca2 = mean_pca2[use,]
  mean_p_value = mean_p_value[use,]
  #o = rownames(mean_pca1 )[order(rowSums(mean_pca1[,1]))]
  o = rownames(mean_pca1 )[order(mean_pca1[,1])]
  o = o[c(grep("B cell",o),grep("T cell",o),grep("NK cell ", o), grep("ILC", o),grep("Myeloid",o))]
  o = unique(o)
  range = max(c(mean_pca1,mean_pca2))
  range_x = length(rownames(mean_pca1))
  mean_pca1 = mean_pca1[o,]
  mean_pca2 = mean_pca2[o,]
  
  fileout1=concat(c(output_directory, "Receptor_ligand_analysis_in_out_degree_pca_per_cell_type_", batch,".pdf"))
  w=3.7
  pdf(file=fileout1, height=w*1.6*1, width=w*3)
  par(mfrow= c(1,1), mar = c(15,4,2.5,1))
  
  #names(groups) = colnames(mat_stat)
  factors = o
  main = concat(c("in/out degree"))
  b = (range)*0.035
  max_width = length(o)
  min = -1*range
  max = range
  if(range>50){scale = c(-100:100)*20}
  if(range>1000){scale = c(-100:100)*250}
  if(range>2000){scale = c(-100:100)*500}
  if(range>5000){scale = c(-100:100)*1000}
  if(range>20000){scale = c(-100:100)*5000}
  if(range<=50){scale = c(-100:100)*10}
  if(range<=30){scale = c(-100:100)*5}
  if(range <15){scale = c(-100:100)*2.5}
  if(range <5){scale = c(-100:100)*1}
  if(range <4){scale = c(-100:100)*0.5}
  if(range <1.5){scale = c(-100:1000)*0.2}
  if(range <0.5){scale = c(-100:100)*0.1}
  if(range <0.1){scale = c(-100:100)*0.01}
  if(range <0.01){scale = c(-100:100)*0.001}
  cex = 0.9
  Fun<-function(x){x}
  scale = scale[intersect(which(scale<= max), which(scale>=min))]
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max)*1.2)
  mtext(side = 2, text = "out degree", line = 2.8,cex= cex-0.1, at = range/2, las = 3, font = 1)
  mtext(side = 2, text = "in degree", line = 2.8,cex= cex-0.1, at = range/-2, las = 3, font = 1)
  segments(Fun(scale),length(groups)+0.5,Fun(scale),0.5,col = "grey",lwd = 1,lty = 3 )
  segments(0.5,Fun(scale),max_width +0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  segments(0.5,0,max_width +0.5,0,col = "grey",lwd = 3,lty = 1 )
  scale1 = scale
  w = which(scale<0) 
  scale1[w] = scale1[w]*-1
  mtext(side = 2, text = scale1, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  
  width = 0.18
  index = 1
  l1 = 2
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.2/max(shift)
  
  for(i in c(1:length(mean_pca1[,1]))){
    points1=as.numeric(mean_pca1[o[i],1])
    points2=as.numeric(mean_pca2[o[i],1])
    rect(i-shift[1]+ width,0,i-shift[1]- width, points1, col = cols[1],border = NA)
    rect(i-shift[2]+ width,0,i-shift[2]- width, points2, col = cols[2],border = NA)
  }
  
  for(i in c(1:length(mean_pca1[,1]))){
    points1=as.numeric(mean_pca1[o[i],2])
    points2=as.numeric(mean_pca2[o[i],2])
    rect(i-shift[1]+ width,0,i-shift[1]- width, -1*points1, col = cols[1],border = NA)
    rect(i-shift[2]+ width,0,i-shift[2]- width, -1*points2, col = cols[2],border = NA)
  }
  o1 = gsub("myeloid ","",o)
  mtext(side = 1, text = o1, line = 0.05,cex= cex-0.1,  at =c(1:length(o)), las = 2, font = 1)
  
  for(i in c(1:length(mean_pca1[,1]))){
    pval = mean_p_value[o[i],1]
    if(pval<0.05){
      points1=as.numeric(mean_pca1[o[i],1])
      points2=as.numeric(mean_pca2[o[i],1])
      b = max*0.045
      y = max(c(points1, points2))+b
      pval1 = "*"
      text(i, y, labels = pval1, cex = 1.4)
    }
  }
  for(i in c(1:length(mean_pca1[,1]))){
    pval = mean_p_value[o[i],2]
    if(pval<0.05){
      points1=as.numeric(mean_pca1[o[i],2])
      points2=as.numeric(mean_pca2[o[i],2])
      b = max*0.045
      y = max(c(points1, points2))+b
      pval1 = "*"
      text(i, -1*y, labels = pval1, cex = 1.4)
    }
  }
  
  dev.off()
}
Get_network_centrality_between_patient_groups(input_directory, batch, output_directory,rec_ligand_list, list_names, cols1, cols, cols2, group_PCA_list, factor_PCA)
  
  
Treg_analyses<-function(input_directory, batch, output_directory,rec_ligand_list, list_names, cols1, cols, cols2, group_PCA_list, factor_PCA){
  analysis = "Chemokine_analysis"

  file = concat(c(input_directory, "Seurat_cell_counts_",batch,"_all.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[which(p[,"Group"]=="CD45"),]
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
  list_count = NULL
  
  include = NULL
  for(ind in c(1:length(cell_groups))){
    w = which(cell_group== cell_groups[ind])
    if(cell_groups[ind] =="CD45:broad"){
      w = intersect(w, which(p[,1] %in% c( "B cell"  , "Myeloid","NK" , "T cell"   )))
    }
    cell_type_sub = p[w,1]
    p1 = p[w,]
    m_count_norm = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(sample), dimnames = c(list(cell_type_sub), list(sample)))
    m_count = matrix(data = 0,nrow = length(cell_type_sub), ncol = length(sample), dimnames = c(list(cell_type_sub), list(sample)))
    
    for(s in c(1:length(sample))){
      m_count_norm[, sample[s]]= as.numeric(p1[, sample[s]])*100/sum(as.numeric(p1[, sample[s]]))
      m_count[, sample[s]]= as.numeric(p1[, sample[s]])*100/sum(as.numeric(p1[, sample[s]]))
      if(sum(as.numeric(p1[, sample[s]]))<7){m_count_norm[, sample[s]] = -1}
    }
    if(length(which(apply(m_count_norm, 2, max)!=-1))>3){include = c(include, ind)}
    list_count = c(list_count, list(t(m_count[order(rownames(m_count)),])))
    m_count_norm= m_count_norm[order(rownames(m_count_norm)),]
    list_normalised= c(list_normalised, list(t(m_count_norm)))
  }
  names(list_normalised) = cell_groups
  names(list_count) = cell_groups
  list_normalised = list_normalised[include]
  
  m_count = list_normalised[[1]]
  cell_types_biopsy = names(which(colSums(m_count[grep("biopsy", rownames(m_count)),])>0))
  
  cell_types = c("all")
  list_counts = NULL
  list_expression = NULL
  
  file = concat(c(input_directory,"Cytokines_counts_ligand_receptor_genes_all_between_groups_",batch,".txt"))
  p <- as.matrix(read.csv(file, head=F, sep="\t"))
  headers = p[1,]
  headers[grep("T cell CD8 Senescent", headers)]
  
  headers = as.character(headers[c(3:(length(headers)-1))])
  headers1 = headers[grep(" .perc_expressing_per_cell_type", headers)]
  headers2 = headers[grep(" .total_RNA_per_cell_type", headers)]
  headers3 = headers[grep(" .prop_per_CD45_expressing", headers)]
  headers4 = headers[grep(" .n_expressing_per_cell_type", headers)]
  headers5 = headers[grep(" .raw_per_cell_type", headers)]
  
  headers_trim = gsub(" .perc_expressing_per_cell_type","", headers1)
  #check
  headers_trim2 = gsub(" .total_RNA_per_cell_type","", headers2)
  headers_trim3 = gsub(" .prop_per_CD45_expressing","", headers3)
  headers_trim4 = gsub(" .n_expressing_per_cell_type","", headers4)
  headers_trim5 = gsub(" .raw_per_cell_type","", headers5)
  headers_trim==headers_trim2
  headers_trim==headers_trim3
  headers_trim==headers_trim4
  headers_trim==headers_trim5
  #
  headers_trim_intersect = intersect(headers_trim,cell_types_biopsy)
  p = p[which(p[,2]=="pancreatic"),]
  colnames(p) = c("gene","site","pat", headers)
  gene = p[,1]
  site = p[,2]
  patient =p[,3]
  genes = sort(unique(gene))
  patients = sort(unique(patient))
  genes = sort(unique(gene))
  perc_expressing_per_cell_type = NULL
  total_RNA_per_cell_type = NULL
  prop_per_CD45_expressing = NULL
  n_expressing_per_cell_type = NULL
  raw_per_cell_type = NULL
  for(g in c(1:length(genes))){
    m1 = matrix(data = 0, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    m2 = matrix(data = -1, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    m3 = matrix(data = 0, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    m4 = matrix(data = 0, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    m5 = matrix(data = 0, nrow = length(patients), ncol = length(headers_trim), dimnames = c(list(patients), list(headers_trim)))
    w = which(gene== genes[g])
    for(i in c(1:length(headers_trim))){
      m1[patient[w], gsub(" .perc_expressing_per_cell_type","",headers1[i])] = as.numeric(p[w, headers1[i]])
      m2[patient[w],  gsub(" .total_RNA_per_cell_type","",headers2[i])] = as.numeric(p[w, headers2[i]])
      m3[patient[w],  gsub(" .prop_per_CD45_expressing","",headers3[i])] = as.numeric(p[w, headers3[i]])
      m4[patient[w],  gsub(" .n_expressing_per_cell_type","",headers4[i])] = as.numeric(p[w, headers4[i]])
      m5[patient[w],  gsub(" .raw_per_cell_type","",headers5[i])] = as.numeric(p[w, headers5[i]])
    }
    m1 = m1 [,which(headers_trim!="-")]
    m2 = m2 [,which(headers_trim!="-")]
    m3 = m3 [,which(headers_trim!="-")]
    m4 = m4 [,which(headers_trim!="-")]
    m5 = m5 [,which(headers_trim!="-")]
    cn= colnames(m1)
    #w = setdiff(c(1:length(cn)), c(grep("B cell", cn), grep("T cell", cn)))
    #cn [w] = paste("myeloid", cn [w])
    colnames(m1) = cn
    colnames(m2) = cn
    colnames(m3) = cn
    colnames(m4) = cn
    colnames(m5) = cn
    perc_expressing_per_cell_type = c(perc_expressing_per_cell_type, list(m1[,headers_trim_intersect]))
    total_RNA_per_cell_type = c(total_RNA_per_cell_type, list(m2[,headers_trim_intersect]))
    prop_per_CD45_expressing =c(prop_per_CD45_expressing,list(m3[,headers_trim_intersect]) )
    n_expressing_per_cell_type =c(n_expressing_per_cell_type,list(m4[,headers_trim_intersect]) )
    raw_per_cell_type =c(raw_per_cell_type,list(m5[,headers_trim_intersect]) )
  }
  names(perc_expressing_per_cell_type)= genes
  names(total_RNA_per_cell_type)= genes
  names(prop_per_CD45_expressing)= genes
  names(n_expressing_per_cell_type)= genes
  names(raw_per_cell_type)= genes
  
  a = perc_expressing_per_cell_type[[1]]
  for(j in c(1:length(a[,1]))){a[j,which(a[j,]<0)] = 0}
  x = a
  for(i in c(2:length(perc_expressing_per_cell_type))){
    a = perc_expressing_per_cell_type[[i]]
    for(j in c(1:length(a[,1]))){a[j,which(a[j,]<0)] = 0}
    x = x+a
  }
  
  
  #### genes of interest
  ligand = c("CCL18","CCL1", "CCL8", "CCL16","CXCL16","CCL17","CCL22","CCL19","CCL21","CCL20", "CXCL12")
  receptor = c("CCR8","CCR8","CCR8","CCR8","CXCR6","CCR4","CCR4", "CCR7", "CCR7", "CCR6", "CXCR4")
  genes_sub = names(total_RNA_per_cell_type)
  w = intersect(which(ligand %in% genes_sub),which(receptor %in% genes_sub))
  ligand = ligand[w]
  receptor = receptor[w]
  
  dir = c(1,-1)
  
  which(ligand %in% names(perc_expressing_per_cell_type))
  which(receptor %in% names(perc_expressing_per_cell_type))
  
  Means_factor = function(factor, x){
    m = NULL
    for(i1 in c(1:length(levels(factor)))){
      x1 = x[which(factor==levels(factor)[i1])]
      x1 = x1[which(x1!=-1)]
      m = c(m, mean(x1))}
    return(m)}
  
  library(RColorBrewer)	
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)	
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("_biopsy","",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  summary_stats = NULL
  fileout1=concat(c(output_directory,"Receptor_ligand_analysis_biopsy_per_pair_Treg_subset", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*2.5*15, width=w*6.5)
  par(mfrow= c(25,2), mar = c(20,4,3.5,1))
  
  genes = NULL
  for(i in c(1:length(ligand))){
    genes =c(genes, ligand[i], receptor[i])}
  headers_trim = colnames(total_RNA_per_cell_type[[1]])
  
  for(ind in c(1:length(genes))){
    print(ind)
    mean_expressing1 = total_RNA_per_cell_type[[genes[ind]]]/raw_per_cell_type[[genes[ind]]]
    mean_expressing1 = mean_expressing1 [,colnames(mean_expressing1)[grep(" b",colnames(mean_expressing1), invert = T)]]
    mat_stat= mean_expressing1[group_PCA_list,]
    w1 = which(apply(mat_stat[which(factor_PCA ==1),], 2, function(x){length(intersect(which(x!=0), which(is.na(x)==F)))})>=3)
    w2 = which(apply(mat_stat[which(factor_PCA ==2),], 2, function(x){length(intersect(which(x!=0), which(is.na(x)==F)))})>=3)
    w = intersect(w1,w2)
    w1 = which(apply(mat_stat[which(factor_PCA ==2),], 2, function(x){length(which(is.na(x)==T))})<=1)
    w = intersect(w,w1)
    if(length(w)>=2){
      mat_stat1 = mat_stat[,w]
      factor = factor(factor_PCA)
      mat_stat1 = mat_stat1^0.5
      fit = manova(formula = mat_stat1 ~ factor)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      inc = NULL
      i1 = 0
      for(i in p1){
        i1 = i1+1
        if(length(i$'Pr(>F)'[1])!=0){
          inc = c(inc, i1)
          p_value = c(p_value, i$'Pr(>F)'[1]) 
          if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
          }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
        }
      }
      if(length(p_value)>2){
        nam = colnames(mat_stat1)[inc]
        p_value[which(is.na(p_value))] = 2
        names(p_value) = nam
        head(p_value,5)
        colnames(means) = paste("mean.group.", c(1:length(means[1,])))
        combined_p_value = cbind(p_value ,means)
        rownames(combined_p_value) = nam
        p.site = rep(concat(c("biopsy ", genes[ind], length(nam))))
        p.analysis = rep(analysis, length(nam))
        x = cbind(p.site, p.analysis, combined_p_value)
        if(length(summary_stats)==0){summary_stats = x
        }else{summary_stats = rbind(summary_stats,x)}
      }
      groups = NULL
      for(g1 in c(1:length(mat_stat[1,]))){
        groups_sub = NULL
        for(g2 in c(1:length(groups_PCA))) {
          x = mat_stat[gsub("_biopsy","",groups_PCA[[g2]]),g1]
          x=x[which(is.na(x)==F)]
          groups_sub  = c(groups_sub , list(x))
        }
        groups = c(groups, list(groups_sub))
      }
      names(groups) = colnames(mat_stat)
      p_values = p_value
      factors1 = paste("group",c(1:length(groups_PCA)))
      factors = names(groups)
      main = concat(c("biopsy\n",genes[ind]))
      max = max(c(unlist(groups), unlist(groups))*1.2)
      min = 0
      b = (max-min)*0.035
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 70
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(0:100)*20}
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
      plot(c(2, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = gsub("myeloid ","",factors), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
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
        for(i1 in c(1:length(groups[[i]]))){
          points1=as.numeric(groups[[i]][[i1]])
          box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
          Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
          points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
        }
      }
      y = max(unlist(groups))
      if(length(p_values)>0){
        for(i1 in c(1:length(p_values))){
          i = which(factors==names(p_values)[i1])
          b = max*0.035
          signif_threshold = 0.05
          if(p_values[i1]<signif_threshold){
            pval1 = "*"
            y = max(unlist(groups[[i]]))+b
            text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.4)
          }
        }
      }
    }
  }
  
  dev.off()
  
  summary_stats = NULL
  fileout1=concat(c(output_directory,"Receptor_ligand_analysis_biopsy_per_pair_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*2.5*15, width=w*6)
  par(mfrow= c(25,2), mar = c(20,4,3.5,1))
  
  genes = c(ligand, receptor)
  genes = names(total_RNA_per_cell_type)
  headers_trim = colnames(total_RNA_per_cell_type[[1]])
  means_group_1 = matrix(data = 0, nrow = length(genes), ncol = length(headers_trim), dimnames = c(list(genes), list(headers_trim)))
  means_group_2 = matrix(data = 0, nrow = length(genes), ncol = length(headers_trim), dimnames = c(list(genes), list(headers_trim)))
  
  for(ind in c(1:length(genes))){
    print(ind)
    mean_expressing1 = total_RNA_per_cell_type[[genes[ind]]]/raw_per_cell_type[[genes[ind]]]
    mean_expressing1 = mean_expressing1 [,colnames(mean_expressing1)[grep(" b",colnames(mean_expressing1), invert = T)]]
    mat_stat= mean_expressing1[group_PCA_list,]
    w1 = which(apply(mat_stat[which(factor_PCA ==1),], 2, function(x){length(intersect(which(x!=0), which(is.na(x)==F)))})>=3)
    w2 = which(apply(mat_stat[which(factor_PCA ==2),], 2, function(x){length(intersect(which(x!=0), which(is.na(x)==F)))})>=3)
    w = intersect(w1,w2)
    w1 = which(apply(mat_stat[which(factor_PCA ==2),], 2, function(x){length(which(is.na(x)==T))})<=1)
    w = intersect(w,w1)
    means_group_1[genes[ind], colnames(mat_stat)] = apply(mat_stat[which(factor_PCA ==1),], 2, function(x){mean( x[which(is.na(x)==F)])})
    means_group_2[genes[ind], colnames(mat_stat)] = apply(mat_stat[which(factor_PCA ==2),], 2, function(x){mean( x[which(is.na(x)==F)])})
    if(length(w)>=3){
      mat_stat1 = mat_stat[,w]
      factor = factor(factor_PCA)
      mat_stat1 = mat_stat1^0.5
      fit = manova(formula = mat_stat1 ~ factor)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      inc = NULL
      i1 = 0
      for(i in p1){
        i1 = i1+1
        if(length(i$'Pr(>F)'[1])!=0){
          inc = c(inc, i1)
          p_value = c(p_value, i$'Pr(>F)'[1]) 
          if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
          }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
        }
      }
      if(length(p_value)>2){
        nam = colnames(mat_stat1)[inc]
        p_value[which(is.na(p_value))] = 2
        names(p_value) = nam
        head(p_value,5)
        colnames(means) = paste("mean.group.", c(1:length(means[1,])))
        combined_p_value = cbind(p_value ,means)
        rownames(combined_p_value) = nam
        p.site = rep(concat(c("biopsy ", genes[ind], length(nam))))
        p.analysis = rep(analysis, length(nam))
        x = cbind(p.site, p.analysis, combined_p_value)
        if(length(summary_stats)==0){summary_stats = x
        }else{summary_stats = rbind(summary_stats,x)}
      }
      groups = NULL
      for(g1 in c(1:length(mat_stat[1,]))){
        groups_sub = NULL
        for(g2 in c(1:length(groups_PCA))) {
          x = mat_stat[gsub("_biopsy","",groups_PCA[[g2]]),g1]
          x=x[which(is.na(x)==F)]
          groups_sub  = c(groups_sub , list(x))
        }
        groups = c(groups, list(groups_sub))
      }
      names(groups) = colnames(mat_stat)
      p_values = p_value
      factors1 = paste("group",c(1:length(groups_PCA)))
      factors = names(groups)
      main = concat(c("biopsy\n",genes[ind]))
      max = max(c(unlist(groups), unlist(groups))*1.2)
      min = 0
      b = (max-min)*0.035
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 63
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(0:100)*20}
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
      plot(c(2, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
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
        for(i1 in c(1:length(groups[[i]]))){
          points1=as.numeric(groups[[i]][[i1]])
          box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
          Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
          points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
        }
      }
      y = max(unlist(groups))
      if(length(p_values)>0){
        for(i1 in c(1:length(p_values))){
          i = which(factors==names(p_values)[i1])
          b = max*0.035
          signif_threshold = 0.05
          if(p_values[i1]<signif_threshold){
            pval1 = "*"
            y = max(unlist(groups[[i]]))+b
            text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.4)
          }
        }
      }
    }
  }
  
  dev.off()
  
  
  outfile = concat(c(output_directory,"Receptor_ligand_analysis_biopsy_per_pair_", batch,".txt"))
  write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  
  
  for(i in c(1:length(genes))){
    w = which(is.na(means_group_1[i,]))
    means_group_1[i,w] = 0
    w = which(is.na(means_group_2[i,]))
    means_group_2[i,w] = 0
  }
  
  means_group = (means_group_1+ means_group_2)/2
  # ligand = c("CCL18","CXCL16","CCL17","CCL22")
  # receptor = c("CCR8","CXCR6","CCR4","CCR4")
  pair_names = apply(cbind(ligand, receptor), 1, paste, collapse = "-")
  
  library(igraph)
  
  from = NULL
  to = NULL
  strength = NULL
  type = NULL
  
  ligand_receptor = apply(cbind(ligand, receptor),1,paste, collapse = "-")
  
  cell_type = names(means_group[1,])
  cell_type = cell_type [grep(" b", cell_type, invert = T)]
  rank_mat = matrix(data = 0, nrow = length(ligand_receptor), ncol = length(cell_type), dimnames = c(list(ligand_receptor), list(cell_type)))
  cell_type_of_interset= "T cell CD4 Treg Activated" 
  
  for(ind in c(1:length(ligand))){
    rec = means_group[receptor[ind],]
    lig = means_group[ligand[ind],]
    cell_type = names(lig)
    for(i1 in c(1:length(rec))){
      for(i2 in c(1:length(rec))){
        weight = rec[i1]*lig[i2]
        if(weight!=0){
          from = c(from, cell_type[i2])
          to = c(to, cell_type[i1])
          strength = c(strength, weight)
          type = c(type, ind)
          if(cell_type[i1] ==cell_type_of_interset){
            rank_mat[ligand_receptor[ind],cell_type[i2]] = weight
          }
        }
      }
    }
  }
  apply(rank_mat, 1, max)
  
  max(strength[which(type==1)])
  
  # w = which(strength>=quantile(strength, 0.5))
  w = which(to == cell_type_of_interset)
  from = from[w]
  to = to[w]
  strength = strength[w]
  type = type[w]
  
  w = grep(" b", from, invert = T)
  from = from[w]
  to = to[w]
  strength = strength[w]
  type = type[w]
  
  
  g <- graph.empty( n=0, directed=T)
  cell_type1 = cell_type [grep(" b", cell_type, invert = T)]
  cell_type1 = c(cell_type1, "N")
  g <- igraph::add.vertices(g, length(cell_type1), name= cell_type1) 
  names <- V(g)$name
  ids <- 1:length(names)
  names(ids) <- names
  edges <- matrix(c(ids[from], ids[to]), nc=2)
  g <- add.edges(g, t(edges), weight= strength)		
  
  d1 = degree(g, v = V(g), mode = c( "out"),loops = TRUE, normalized = FALSE)^2
  library(RColorBrewer)
  cols=add.alpha(colorRampPalette(c('white','yellow','orange', "red"))(100), alpha = 0.8)
  m1 = d1*100/max(d1)
  m1 = round(m1, digits = 0)
  range(m1)
  sizes = 1
  sizes_scaled = sizes^0.3
  sizes_scaled = sizes_scaled*5
  
  V(g)$size<-0
  V(g)$label.cex<-0.5
  V(g)$color = "grey"
    layout1 =layout_in_circle(g)
    # layout1 = layout_with_graphopt(g,niter = 500)
    edge_strength_plot = strength^0.5
    edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
    
    cols=add.alpha(colorRampPalette(c("gold",'orange', "red","darkred","black"))(100), alpha = 0.6)
    
    library(RColorBrewer)	
    cols1 =  add.alpha (brewer.pal(8, "Set2"), alpha = 0.95)
    cols1[6] = add.alpha ("black", alpha = 0.45)
    E(g)$color = cols1[type]
    
    layout1 = cbind(0,c(1:length(cell_type1)))
    rownames(layout1) = cell_type1
    w = which(cell_type1 %in% c(cell_type_of_interset,"N")==F)
    layout1[w,1] = -0.5
    l= length(which(cell_type1 %in% c(cell_type_of_interset,"N")==F))
    layout1[w,2] = seq(from = -1, t = 1, length = l)
    w = which(cell_type1 == cell_type_of_interset)
    # layout1[w,2] = 0
    layout1[w,2] = mean(layout1[which(cell_type1 %in% c(cell_type_of_interset,"N")==F),2])
    w = which(cell_type1 == "N")
    layout1[w,] = c(1,0)
    
    
    fileout1=concat(c(output_directory,"Receptor_ligand_analysis_network_biopsy_Treg_", batch,".pdf"))
    w=20
    pdf(file=fileout1, height=w*0.5, width=w*1)
    par(mfrow= c(1,1), mar = c(1,1,3,1))
    V(g)$name = gsub("myeloid ","",V(g)$name)
    plot(g, layout=layout1, edge.color= E(g)$color, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.000000001, edge.arrow.size = 1, edge.lty = 1, main = "Treg connections",cex.main =3,xlim = range(layout1[,1])+c(-1.3,2.5),ylim = range(layout1[,2])+c(0,0))
    
    w = which(cell_type1 != cell_type_of_interset)
    text(layout1[w,1]-0.5, layout1[w,2], labels = V(g)$name[w], cex = 0.7, pos = 2, font = 1)
    w = which(cell_type1 == cell_type_of_interset)
    text(layout1[w,1]-0.2, layout1[w,2]+0.05, labels = gsub("T cell","", V(g)$name[w]), cex = 0.7, pos = 3, font = 2,offset = -2.5)
    
    legend("topright", pair_names, pch = 21,cex= 1, bty="n", pt.bg = cols1, col = cols1, pt.lwd = 2, text.font = 1)
    dev.off()
    
    ##### only show top X% of cell cell interactions
    cell_type2 = c(unique(rev(from[order(strength)])[c(1:30)]), unique(to))
    cell_type2 = c(cell_type2, "N")
    g <- graph.empty( n=0, directed=T)
    g <- igraph::add.vertices(g, length(cell_type2), name= cell_type2) 
    names <- V(g)$name
    ids <- 1:length(names)
    names(ids) <- names
    w1 = which(from %in% cell_type2)
    edges <- matrix(c(ids[from[w1]], ids[to[w1]]), nc=2)
    g <- add.edges(g, t(edges), weight= strength[w1])		
    
    d1 = degree(g, v = V(g), mode = c( "out"),loops = TRUE, normalized = FALSE)^2
    library(RColorBrewer)
    cols=add.alpha(colorRampPalette(c('white','yellow','orange', "red"))(100), alpha = 0.8)
    m1 = d1*100/max(d1)
    m1 = round(m1, digits = 0)
    range(m1)
    sizes = 1
    sizes_scaled = sizes^0.3
    sizes_scaled = sizes_scaled*5
    
    V(g)$size<-0
    V(g)$label.cex<-0.5
    V(g)$color = "grey"
      layout1 =layout_in_circle(g)
      # layout1 = layout_with_graphopt(g,niter = 500)
      edge_strength_plot = strength[w1]^0.25
      edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
      
      cols=add.alpha(colorRampPalette(c("gold",'orange', "red","darkred","black"))(100), alpha = 0.6)
      
      library(RColorBrewer)	
      cols1 =  add.alpha (brewer.pal(8, "Set2"), alpha = 0.95)
      cols1[6] = add.alpha ("black", alpha = 0.45)
      E(g)$color = cols1[type[w1]]
      
      layout1 = cbind(0,c(1:length(cell_type2)))
      rownames(layout1) = cell_type2
      w = which(cell_type2 %in% c(cell_type_of_interset,"N")==F)
      layout1[w,1] = -0.5
      l= length(which(cell_type2 %in% c(cell_type_of_interset,"N")==F))
      layout1[w,2] = seq(from = -1, t = 1, length = l)
      w = which(cell_type2 == cell_type_of_interset)
      # layout1[w,2] = 0
      layout1[w,2] = mean(layout1[which(cell_type2 %in% c(cell_type_of_interset,"N")==F),2])
      w = which(cell_type2 == "N")
      layout1[w,] = c(1,0)
      
      fileout1=concat(c(output_directory,"Receptor_ligand_analysis_network_biopsy_Treg_", batch,"_top.pdf"))
      w=8
      pdf(file=fileout1, height=w*0.5, width=w*1)
      par(mfrow= c(1,1), mar = c(1,1,3,1))
      plot(g, layout=layout1, edge.color= E(g)$color, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.000000001, edge.arrow.size = 0, edge.lty = 1, main = "Treg connections",cex.main =3,xlim = range(layout1[,1])+c(-1.8,2.5),ylim = range(layout1[,2])+c(0,0))
      
      w = which(cell_type2 %in% c("N",cell_type_of_interset)==F)
      text(layout1[w,1]-0.5, layout1[w,2], labels = V(g)$name[w], cex = 0.7, pos = 2, font = 1)
      w = which(cell_type2 == cell_type_of_interset)
      text(layout1[w,1]+0.1, layout1[w,2]+0.1, labels = gsub("reg ","reg\n", gsub("T cell","", V(g)$name[w])), cex = 0.7, pos = 3, font = 2,offset = -2.5)
      
      legend("topright", pair_names, pch = 21,cex= 1, bty="n", pt.bg = cols1, col = cols1, pt.lwd = 2, text.font = 1)
      dev.off()
      
      
      ### contribution to Treg cytokine signalling
      
      nam_sub = names(rev(sort(apply(rank_mat,2,  max))))[c(1:20)]
      mat_sub = rank_mat[, nam_sub]^0.5
      
      factors = gsub("myeloid ","",colnames(mat_sub))
      main = "interaction strength"
      max = max(mat_sub*1.01)
      min = 0
      b = (max-min)*0.035
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 40
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(0:100)*20}
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
      
      fileout1=concat(c(output_directory,"Receptor_ligand_analysis_network_biopsy_Treg_", batch,"_top_connected.pdf"))
      w=3
      pdf(file=fileout1, height=w*1.25*2, width=w*2.3*2)
      par(mfrow= c(2,2), mar = c(14,4,4,1))
      
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(factors)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      
      library(RColorBrewer)	
      cols1 =  add.alpha (brewer.pal(8, "Set2"), alpha = 0.95)
      cols1[6] = add.alpha ("black", alpha = 0.45)
      
      for(i1 in c(1:length(mat_sub[,1]))){
        points1=mat_sub[i1,]
        points(c(1:length(points1)), points1, type = "l", col= cols1[i1], lwd= 1.3, lty = 1)
        points(c(1:length(points1)), points1, pch =21, bg= cols1[i1],col= cols1[i1], cex = 0.9)
      }
      
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ="", axes = FALSE, ylim = c(min, max))
      legend("topleft", rownames(mat_sub), pch = 21,cex= 0.5, bty="n", pt.bg = cols1, col = cols1, pt.lwd = 2, text.font = 1)
      
      
      ### show only CCR8 links
      w = grep("CCR8",rownames(rank_mat))
      nam_sub = names(rev(sort(apply(rank_mat[w,],2,  max))))[c(1:20)]
      mat_sub = rank_mat[w, nam_sub]^0.5
      
      factors = gsub("myeloid ","",colnames(mat_sub))
      main = "interaction strength"
      max = max(mat_sub*1.01)
      min = 0
      b = (max-min)*0.035
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 40
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(0:100)*20}
      if(range<=50){scale = c(0:100)*10}
      if(range<=30){scale = c(0:100)*5}
      if(range <15){scale = c(0:100)*2.5}
      if(range <5){scale = c(0:100)*1}
      if(range <4){scale = c(0:100)*0.5}
      if(range <1.5){scale = c(0:1000)*0.2}
      if(range <1){scale = c(0:1000)*0.1}
      if(range <0.5){scale = c(0:100)*0.1}
      if(range <0.1){scale = c(0:100)*0.01}
      if(range <0.01){scale = c(0:100)*0.001}
      cex = 0.9
      Fun<-function(x){x}
      
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "CCR8-Tregs", axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(factors)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      
      library(RColorBrewer)	
      cols1 =  add.alpha (brewer.pal(8, "Set2"), alpha = 0.95)
      cols1[6] = add.alpha ("black", alpha = 0.45)
      
      for(i1 in c(1:length(mat_sub[,1]))){
        points1=mat_sub[i1,]
        points(c(1:length(points1)), points1, type = "l", col= cols1[i1], lwd= 1.3, lty = 1)
        points(c(1:length(points1)), points1, pch =21, bg= cols1[i1],col= cols1[i1], cex = 0.9)
      }
      
      plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main ="", axes = FALSE, ylim = c(min, max))
      legend("topleft", rownames(mat_sub), pch = 21,cex= 0.5, bty="n", pt.bg = cols1, col = cols1, pt.lwd = 2, text.font = 1)
      
      dev.off()
      
      
      
  
}
Treg_analyses(input_directory, batch, output_directory,rec_ligand_list, list_names, cols1, cols, cols2, group_PCA_list, factor_PCA)


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
batch = "PDAC150Ka"
analysis = "VDJ_clonality"

## get patient groups for plotting: 
groups_PCA =readRDS(file = concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/Overall UMAPs and annotations/PCA_cell_counts_blood_PCA_groups_PDAC150Ka.RDS")))
input_directory_groups  = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/PROPORTIONS/"

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

Get_TCR_information<-function(input_directory_groups){
  type = "TCR"
  type1 = "T_cells"
  file = concat(c(input_directory_groups,"VDJ_information_",type,"_PDAC150Ka.txt"))
  VDJ <- as.matrix(read.csv(file, head=T, sep="\t"))
  clone = VDJ[,"clone1"]
  
  file = concat(c(input_directory_groups,"Cell_annotation_ALL_PDAC150Ka.txt"))
  p1 <- as.matrix(read.csv(file, head=T, sep="\t"))
  p1 <- p1[grep("T cell", p1[,"cell_type"]), ]
  id1 = p1[,"barcode"]
  sample = p1[,"sample"]
  id = apply(cbind(id1, sample), 1, function(x){
    x1 = gsub(concat(c(x[2], "_")), "", x[1])
    return(concat(c(x1,"||", x[2])))})
  id = gsub("-1","", id, fixed = T)
  names(id) = id1
  length(intersect(id, names(clone)))
  cell_type = p1[,"cell_type"]
  Patient = p1[,"sample1"]
  Sample.Type= p1[,"sample"]
  Patient = gsub("_CD45p1","", Patient)
  Patient = gsub("_CD45p2","", Patient)
  Sample.Type[grep("biopsy", Sample.Type)] = "biopsy"
  Sample.Type[grep("blood", Sample.Type)] = "blood"
  sample = apply(cbind(Patient, Sample.Type), 1, paste, collapse = "-")
  
  file = concat(c(input_directory_groups,"Cell_names_broad.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p=p[which(p[,"Full.annotation"]%in% cell_type),]
  p=p[which(p[,"Very.broad.annotation"]%in% c("ILC","T cell gdT","T cell MAIT")==F),]
  Very.broad.annotation = p[,"Very.broad.annotation"]
  Broad.annotation = p[,"Broad.annotation"]
  names(Very.broad.annotation) = p[,"Full.annotation"]
  names(Broad.annotation) = p[,"Full.annotation"]
  
  cell_type1 = Broad.annotation[cell_type]
  
  names(cell_type) = id
  names(cell_type1) = id
  names(sample) = id
  names(Patient) = id
  
  list_return = c(list(VDJ), list(cell_type), list(cell_type1), list(sample), list(Patient), list(Sample.Type))
  names(list_return) = c("VDJ_object", "cell_type","cell_type_broad","sample","Patient","Sample.Type")
  return(list_return)
}

VDJ_list_BCR = Get_BCR_information(input_directory_groups)
VDJ_list_TCR = Get_TCR_information(input_directory_groups)

output_directory = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Data and references/"

### BCR
BCR_duplicate_counts_subsampled_intra<-function(VDJ_list, output_directory, batch){
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
  
  ### get total cell numbers for total and cell subsets
  count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
  totals = rowSums(count_per_sample)
  count_per_sample = cbind(count_per_sample, totals)
  
  print("Calculating duplicate rates")
  m_n_duplicate_clones = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  table(sample[w_clone], cell_type_broad[w_clone])
  repeats = 1000
  thresholds1 = 5
  for(c in c(1:length(cell_types_broad))){
    print (concat(c(c, " of ",length(cell_types_broad)," done" )))
    for(s in c(1:length(samples))){
      w2 = intersect(which(sample == samples[s]), which(cell_type_broad== cell_types_broad[c]))
      id_sub = intersect((names(sample)[w2]), w_clone)
      clon = clone[id_sub]
      clon = clon[which(clon!='-')]
      if(length(clon)>= thresholds1){
        n_dupl = NULL
        if(max(table(clon))==1){n_dupl=0
        }else{
          for(r in c(1:repeats)){
            rand = sample(clon, thresholds1)
            t = table(rand)
            n_dupl = c(n_dupl, length(which(t>1)))
          }}
        m_n_duplicate_clones[samples[s], cell_types_broad[c]] = mean(n_dupl)
      }}
  }
  all = c(cell_types )
  non_naive = setdiff(cell_types,"naive B cells")
  
  groups_cell_type = c(list(non_naive), list(all))
  groups_cell_type_names = c("non-naive","all")
  
  m_n_duplicate_clones_sub = matrix(data = -1,nrow = length(samples), ncol = length(groups_cell_type_names), dimnames = c(list(samples), list(groups_cell_type_names)))
  
  repeats = 1000
  thresholds1 = 5
  for(c in c(1:length(groups_cell_type_names))){
    print (concat(c(c, " of ",length(groups_cell_type_names)," done" )))
    for(s in c(1:length(samples))){
      w2 = intersect(which(sample == samples[s]), which(cell_type %in% groups_cell_type[[c]]))
      id_sub = intersect((names(sample)[w2]), w_clone)
      clon = clone[id_sub]
      clon = clon[which(clon!='-')]
      if(length(clon)>= thresholds1){
        n_dupl = NULL
        if(max(table(clon))==1){n_dupl=0
        }else{
          for(r in c(1:repeats)){
            rand = sample(clon, thresholds1)
            t = table(rand)
            n_dupl = c(n_dupl, length(which(t>1)))
          }}
        m_n_duplicate_clones_sub[samples[s], groups_cell_type_names[c]] = mean(n_dupl)
      }}
  }
  m_n_duplicate_clones1 = cbind(m_n_duplicate_clones, m_n_duplicate_clones_sub)
  
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_nduplicate_clones_intra.txt"))
  write.table(m_n_duplicate_clones1, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}
BCR_duplicate_counts_subsampled_intra(VDJ_list, output_directory,batch)

BCR_duplicate_counts_subsampled_inter<-function(VDJ_list, output_directory, batch){
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
  
  ### get total cell numbers for total and cell subsets
  count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
  totals = rowSums(count_per_sample)
  count_per_sample = cbind(count_per_sample, totals)
  
  print("Calculating inter duplicate rates")
  m_n_duplicate_clones = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  repeats = 500
  thresholds1 = min(c(50,quantile(totals, 0.05)))
  for(s in c(1:length(samples))){
    print (concat(c(s, " of ",length(samples)," done" )))
    w1 = intersect(names(sample)[which(sample == samples[s])], w_clone)
    repeated_measures = matrix(data = 0,nrow = length(c(1:repeats)), ncol = length(cell_types_broad), dimnames = c(list(c(1:repeats)), list(cell_types_broad)))
    if(length(w1)>= thresholds1){
      for(r in c(1:repeats)){
        rand = sample(w1, thresholds1)
        t = table(clone[rand], cell_type_broad[rand])
        expanded = which(rowSums(t)>1)
        if(length(expanded)==1){
          repeated_measures[r,] = repeated_measures[r,]+t[expanded,]
        }else{
          if(length(expanded)>1){
            repeated_measures[r,] = repeated_measures[r,]+colSums(t[expanded,])
          }}}
      if(sum(repeated_measures)>=5){
        m_n_duplicate_clones[samples[s],]= colSums(repeated_measures)*100/sum(repeated_measures)
      }}
    }
      
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_nduplicate_clones_inter.txt"))
  write.table(m_n_duplicate_clones, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}
BCR_duplicate_counts_subsampled_inter(VDJ_list, output_directory,batch)

BCR_renyi_gini_d5_subsampled<-function(VDJ_list, output_directory){
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
  
  count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
  totals = rowSums(count_per_sample)
  count_per_sample = cbind(count_per_sample, totals)
  
  #### clonal expansion from subsampled cells. 
  m_gini = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  m_renyi = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  m_d5 = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
  
  print("Calculating Gini, Renyi and D5")
  thresholds = 15
  repeats = 1000
  library(ineq)
  library(edgeR)
  for(c in c(1:length(cell_types_broad))){
    print (concat(c(c, " of ",length(cell_types_broad)," done" )))
    g1_gini = NULL
    for(s in c(1:length(samples))){
      w2 = intersect(which(sample == samples[s]), which(cell_type_broad== cell_types_broad[c]))
      id_sub = intersect((names(sample)[w2]), w_clone)
      clon = clone[id_sub]
      clon = clon[which(clon!='-')]
      if(length(clon)>= thresholds){
        ginis = NULL
        reynis = NULL
        d5s = NULL
        for(r in c(1:repeats)){
          rand = sample(clon, thresholds)
          dist = table(rand)
          g = ineq::Gini(dist)
          ginis = c(ginis, g)
          dist = dist/sum(dist)
          r = sum(-(dist*log(dist)))
          renyis = c(reynis, r)
          d5 = rev(sort(table(rand)))
          if(length(d5)>=3){d5 = sum(d5[c(1:3)])*100/sum(d5)
          }else{d5 = 100}
          d5s = c(d5s, d5)
        }
        g1_gini = c(g1_gini,list(ginis))
        m_gini[samples[s], cell_types_broad[c]] = mean(ginis)
        m_renyi[samples[s], cell_types_broad[c]] = mean(renyis)
        m_d5[samples[s], cell_types_broad[c]] = mean(d5s)
      }}
  }
  
  gini_all = rep(-1,length(samples))
  names(gini_all)= samples
  renyi_all = gini_all
  d5_all = gini_all
  thresh = floor(min(totals)*0.95)
  for(s in c(1:length(samples))){
    w2 =which(sample == samples[s])
    id_sub = names(sample)[w2]
    clon = clone[id_sub]
    clon = clon[which(clon!='-')]
    if(length(clon)>= thresh){
      ginis = NULL
      reynis = NULL
      for(r in c(1:repeats)){
        rand = sample(clon, thresh)
        dist = table(rand)
        g = ineq::Gini(dist)
        ginis = c(ginis, g)
        dist = dist/sum(dist)
        r = sum(-(dist*log(dist)))
        renyis = c(reynis, r)
        d5 = rev(sort(table(rand)))
        if(length(d5)>5){d5 = sum(d5[c(1:5)])*100/sum(d5)
        }else{d5 = 100}
        d5s = c(d5s, d5)
      }
      gini_all[samples[s]] = mean(ginis)
      renyi_all[samples[s]] = mean(renyis)
      d5_all[samples[s]] = mean(d5s)
    }}	
  
  m_gini = cbind(m_gini, gini_all)
  m_renyi = cbind(m_renyi, renyi_all)
  m_d5 = cbind(m_d5, d5_all)
  colnames(m_d5)[length(colnames(m_d5))] = "all"
  colnames(m_gini)[length(colnames(m_gini))] = "all"
  colnames(m_renyi)[length(colnames(m_renyi))] = "all"
  
  
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_gini.txt"))
  write.table(m_gini, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_renyi.txt"))
  write.table(m_renyi, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_d5.txt"))
  write.table(m_d5, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
}
BCR_renyi_gini_d5_subsampled(VDJ_list, output_directory)

BCR_isotype_SHM_VJ<-function(VDJ_list, output_directory){
  VDJ_list = VDJ_list_BCR
  type = "BCR"
  type1 = "B_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[2]]
  cell_type = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  isotype = VDJ[,"constant_region1"]
  IGHKL = VDJ[,"chain2"]
  SHM1 = VDJ[,"V_mm1"]
  SHM2 = VDJ[,"V_mm2"]
  V_gene1 = VDJ[,"V_gene1"]
  V_gene2 = VDJ[,"V_gene2"]
  VJ_gene1 = apply(cbind(VDJ[,"V_gene1"], VDJ[,"J_gene1"]), 1, paste, collapse = ":")
  VJ_gene2 = apply(cbind(VDJ[,"V_gene2"], VDJ[,"J_gene2"]), 1, paste, collapse = ":")
  SHM1[which(SHM1=="-")] = 0
  SHM2[which(SHM2=="-")] = 0
  SHM1 = as.numeric(SHM1)
  SHM2 = as.numeric(SHM2)
  names(SHM1) = rownames(VDJ)
  names(SHM2) = rownames(VDJ)
  
  isotypes = sort(unique(isotype))
  IGHKLs = sort(unique(IGHKL))
  V_gene1s = sort(unique(V_gene1))
  V_gene2s = sort(unique(V_gene2))
  VJ_gene1s = sort(unique(VJ_gene1))
  VJ_gene2s = sort(unique(VJ_gene2))
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  
  group_ids = NULL
  for(c in c(1:length(cell_types_broad))){
    w = intersect(w_clone, names(cell_type_broad)[which(cell_type_broad == cell_types_broad[c])])
    group_ids = c(group_ids, list(w))
  }
  for(c in c(1:length(cell_types))){
    w = intersect(w_clone, names(cell_type)[which(cell_type == cell_types[c])])
    group_ids = c(group_ids, list(w))
  }
  group_ids = c(group_ids, list(w_clone))
  names(group_ids) = c(apply(cbind(cell_types_broad,"broad"),1, paste, collapse = ":"), 
                       apply(cbind(cell_types,"detailed"),1, paste, collapse = ":") , "All" )
  
  list_isotype_percentage = NULL
  list_IGHKL_percentage = NULL
  list_meanSHM1 = NULL ## by isotype
  list_meanSHM2 = NULL ## by isotype
  list_meanSHM12 = NULL ## by isotype
  list_IGHV_gene_percentage = NULL
  list_IGKLV_gene_percentage = NULL
  list_IGHVJ_gene_percentage = NULL
  list_IGKLVJ_gene_percentage = NULL
  
  for(c in c(1:length(group_ids))){
    m_isotype_percentage = matrix(data = 0,nrow = length(samples), ncol = length(isotypes), dimnames = c(list(samples), list(isotypes)))
    m_IGHKL_percentage = matrix(data = 0,nrow = length(samples), ncol = length(IGHKLs), dimnames = c(list(samples), list(IGHKLs)))
    m_meanSHM1 = matrix(data = NA,nrow = length(samples), ncol = length(isotypes), dimnames = c(list(samples), list(isotypes)))
    m_meanSHM2 = matrix(data = NA,nrow = length(samples), ncol = length(isotypes), dimnames = c(list(samples), list(isotypes)))
    m_meanSHM12 = matrix(data = NA,nrow = length(samples), ncol = length(isotypes), dimnames = c(list(samples), list(isotypes)))
    m_IGHV_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(V_gene1s), dimnames = c(list(samples), list(V_gene1s)))
    m_IGKLV_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(V_gene2s), dimnames = c(list(samples), list(V_gene2s)))
    m_IGHVJ_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(VJ_gene1s), dimnames = c(list(samples), list(VJ_gene1s)))
    m_IGKLVJ_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(VJ_gene2s), dimnames = c(list(samples), list(VJ_gene2s)))
    
    w  = group_ids[[c]]
    t = table(sample[w],isotype[w])
    m_isotype_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],IGHKL[w])
    m_IGHKL_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],V_gene1[w])
    m_IGHV_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],V_gene2[w])
    m_IGKLV_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],VJ_gene1[w])
    m_IGHVJ_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],VJ_gene2[w])
    m_IGKLVJ_gene_percentage[rownames(t), colnames(t)] = t
    
    df = data.frame(w)
    df$SHM1 =SHM1[w]
    df$SHM2 =SHM2[w]
    df$SHM12 =SHM1[w]+SHM2[w]
    df$sample =sample[w]
    df$isotype =isotype[w]
    
    dd<-with(df, aggregate(SHM1, by=list(sample=sample, isotype=isotype), mean))
    dn<-with(df, aggregate(SHM1, by=list(sample=sample, isotype=isotype), length))
    dd = dd[which(dn$x>=5),]
    ddd<-reshape(dd, idvar='sample', timevar='isotype', direction='wide')
    colnames(ddd) = gsub("x.", "", colnames(ddd))
    rownames(ddd) = ddd$sample
    ddd = ddd[,intersect(colnames(ddd), isotypes)]
    m_meanSHM1[rownames(ddd), colnames(ddd)] = as.matrix(ddd)
    
    dd<-with(df, aggregate(SHM2, by=list(sample=sample, isotype=isotype), mean))
    dn<-with(df, aggregate(SHM2, by=list(sample=sample, isotype=isotype), length))
    dd = dd[which(dn$x>=5),]
    ddd<-reshape(dd, idvar='sample', timevar='isotype', direction='wide')
    colnames(ddd) = gsub("x.", "", colnames(ddd))
    rownames(ddd) = ddd$sample
    ddd = ddd[,intersect(colnames(ddd), isotypes)]
    m_meanSHM2[rownames(ddd), colnames(ddd)] = as.matrix(ddd)
    
    dd<-with(df, aggregate(SHM12, by=list(sample=sample, isotype=isotype), mean))
    dn<-with(df, aggregate(SHM12, by=list(sample=sample, isotype=isotype), length))
    dd = dd[which(dn$x>=5),]
    ddd<-reshape(dd, idvar='sample', timevar='isotype', direction='wide')
    colnames(ddd) = gsub("x.", "", colnames(ddd))
    rownames(ddd) = ddd$sample
    ddd = ddd[,intersect(colnames(ddd), isotypes)]
    m_meanSHM12[rownames(ddd), colnames(ddd)] = as.matrix(ddd)
    
    list_isotype_percentage = c(list_isotype_percentage, list( m_isotype_percentage))
    list_IGHKL_percentage = c(list_IGHKL_percentage, list( m_IGHKL_percentage))
    list_meanSHM1 = c(list_meanSHM1, list( m_meanSHM1))
    list_meanSHM2 = c(list_meanSHM2, list( m_meanSHM2))
    list_meanSHM12 = c(list_meanSHM12, list( m_meanSHM12))
    list_IGHV_gene_percentage = c(list_IGHV_gene_percentage, list( m_IGHV_gene_percentage))
    list_IGKLV_gene_percentage = c(list_IGKLV_gene_percentage, list( m_IGKLV_gene_percentage))
    list_IGHVJ_gene_percentage = c(list_IGHVJ_gene_percentage, list( m_IGHVJ_gene_percentage))
    list_IGKLVJ_gene_percentage = c(list_IGKLVJ_gene_percentage, list( m_IGKLVJ_gene_percentage))
  }
  
  names(list_isotype_percentage) = names(group_ids) 
  names(list_IGHKL_percentage) = names(group_ids) 
  names(list_meanSHM1) = names(group_ids) 
  names(list_meanSHM2) = names(group_ids) 
  names(list_meanSHM12) = names(group_ids) 
  names(list_IGHV_gene_percentage) = names(group_ids) 
  names(list_IGKLV_gene_percentage) = names(group_ids) 
  names(list_IGHVJ_gene_percentage ) = names(group_ids) 
  names(list_IGKLVJ_gene_percentage) = names(group_ids) 
  
  list_VDJ_freatures_aggregated  = c(list(list_isotype_percentage),
                                     list(list_IGHKL_percentage),
                                     list(list_meanSHM1),
                                     list(list_meanSHM2),
                                     list(list_meanSHM12),
                                     list(list_IGHV_gene_percentage),
                                     list(list_IGKLV_gene_percentage),
                                     list(list_IGHVJ_gene_percentage),
                                     list(list_IGKLVJ_gene_percentage)
  )
  names(list_VDJ_freatures_aggregated)  = c("isotype usage","IGK/L usage","mean SHM IGH","mean SHM IGK/L","mean SHM IGH+K/L",
                                            "IGHV gene usage", "IGK/L V gene usage","IGHVJ gene usage","IGK/K VJ gene usage")
  
  saveRDS(file = concat(c(output_directory,"VDJ_repertoire_feature_information_",type,"_",batch,".rds")), list_VDJ_freatures_aggregated)
  writeLines(concat(c("\nOutput as RDS object with the following features calculated per cell type:\n",
                      paste(c(" ", names(list_VDJ_freatures_aggregated)), collapse = "\n\t"),
                      "\n\nOutput location:", concat(c(output_directory,"VDJ_repertoire_feature_information_",type,"_",batch,".rds")))))
}
BCR_isotype_SHM_VJ(VDJ_list, output_directory)

### TCR

## identify CD4 and CD8 T cell subtypes: the user may need to edit this depending on annotations
VDJ_list = VDJ_list_TCR
cell_type_broad = VDJ_list[[2]]
cell_type = VDJ_list[[3]]
cell_type_broads = sort(unique(cell_type_broad))
cell_types = sort(unique(cell_type))
CD4 = sort(unique(c(cell_type_broads[grep("CD4",cell_type_broads)], cell_types[grep("CD4",cell_types)])))
CD8 = sort(unique(c(cell_type_broads[grep("CD8",cell_type_broads)], cell_types[grep("CD8",cell_types)], 
                    cell_type_broads[grep("MAIT",cell_type_broads)], cell_types[grep("MAIT",cell_types)])))
CD48 = c(list(CD4),list(CD8))
names(CD48)= c("CD4","CD8")

TCR_duplicate_counts_subsampled_intra<-function(VDJ_list, output_directory, batch, CD48){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  
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
  
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    wx = names(sample)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[wx])
    w_clone = intersect(names(which(clone!='-')),wx)
    
    ### get total cell numbers for total and cell subsets
    count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
    totals = rowSums(count_per_sample)
    count_per_sample = cbind(count_per_sample, totals)
    
    print("Calculating duplicate rates")
    m_n_duplicate_clones = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
    table(sample[w_clone], cell_type_broad[w_clone])
    repeats = 500
    thresholds1 = 5
    for(c in c(1:length(cell_types_broad))){
      print (concat(c(c, " of ",length(cell_types_broad)," done" )))
      for(s in c(1:length(samples))){
        w2 = intersect(which(sample == samples[s]), which(cell_type_broad== cell_types_broad[c]))
        id_sub = intersect((names(sample)[w2]), w_clone)
        clon = clone[id_sub]
        clon = clon[which(clon!='-')]
        if(length(clon)>= thresholds1){
          n_dupl = NULL
          if(max(table(clon))==1){n_dupl=0
          }else{
            for(r in c(1:repeats)){
              rand = sample(clon, thresholds1)
              t = table(rand)
              n_dupl = c(n_dupl, length(which(t>1)))
            }}
          m_n_duplicate_clones[samples[s], cell_types_broad[c]] = mean(n_dupl)
        }}
    }
    all = c(cell_types )
    non_naive = setdiff(cell_types,"naive B cells")
    
    groups_cell_type = c(list(non_naive), list(all))
    groups_cell_type_names = c("non-naive","all")
    
    m_n_duplicate_clones_sub = matrix(data = -1,nrow = length(samples), ncol = length(groups_cell_type_names), dimnames = c(list(samples), list(groups_cell_type_names)))
    
    repeats = 1000
    thresholds1 = 5
    for(c in c(1:length(groups_cell_type_names))){
      print (concat(c(c, " of ",length(groups_cell_type_names)," done" )))
      for(s in c(1:length(samples))){
        w2 = intersect(which(sample == samples[s]), which(cell_type %in% groups_cell_type[[c]]))
        id_sub = intersect((names(sample)[w2]), w_clone)
        clon = clone[id_sub]
        clon = clon[which(clon!='-')]
        if(length(clon)>= thresholds1){
          n_dupl = NULL
          if(max(table(clon))==1){n_dupl=0
          }else{
            for(r in c(1:repeats)){
              rand = sample(clon, thresholds1)
              t = table(rand)
              n_dupl = c(n_dupl, length(which(t>1)))
            }}
          m_n_duplicate_clones_sub[samples[s], groups_cell_type_names[c]] = mean(n_dupl)
        }}
    }
    m_n_duplicate_clones1 = cbind(m_n_duplicate_clones, m_n_duplicate_clones_sub)
    
    out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_nduplicate_clones_intra.txt"))
    write.table(m_n_duplicate_clones1, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  }
}
TCR_duplicate_counts_subsampled_intra(VDJ_list, output_directory,batch,CD48)

TCR_duplicate_counts_subsampled_inter<-function(VDJ_list, output_directory, batch,CD48){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  
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
  
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    wx = names(sample)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[wx])
    w_clone = intersect(names(which(clone!='-')),wx)
    
    table(cell_type[w_clone])
    table(cell_type[w1])
    
    ### get total cell numbers for total and cell subsets
    count_per_sample = table(sample [names(clone)], cell_type [names(clone)])
    totals = rowSums(count_per_sample)
    count_per_sample = cbind(count_per_sample, totals)
    
    print("Calculating inter duplicate rates")
    m_n_duplicate_clones = matrix(data = -1,nrow = length(samples), ncol = length(cell_types_broad), dimnames = c(list(samples), list(cell_types_broad)))
    repeats = 500
    thresholds1 = min(c(50,quantile(totals, 0.05)))
    for(s in c(1:length(samples))){
      print (concat(c(s, " of ",length(samples)," done" )))
      w1 = intersect(names(sample)[which(sample == samples[s])], w_clone)
      repeated_measures = matrix(data = 0,nrow = length(c(1:repeats)), ncol = length(cell_types_broad), dimnames = c(list(c(1:repeats)), list(cell_types_broad)))
      if(length(w1)>= thresholds1){
        for(r in c(1:repeats)){
          rand = sample(w1, thresholds1)
          t = table(clone[rand], cell_type_broad[rand])
          expanded = which(rowSums(t)>1)
          if(length(expanded)==1){
            repeated_measures[r,colnames(t)] = repeated_measures[r,colnames(t)]+t[expanded,]
          }else{
            if(length(expanded)>1){
              repeated_measures[r,colnames(t)] = repeated_measures[r,colnames(t)]+colSums(t[expanded,])
            }}}
        if(sum(repeated_measures)>=5){
          m_n_duplicate_clones[samples[s],]= colSums(repeated_measures)*100/sum(repeated_measures)
        }}
    }
    
    out_file_table = concat(c(output_directory,"VDJ_Clonality_information_",type,"_",batch,"_nduplicate_clones_inter.txt"))
    write.table(m_n_duplicate_clones, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  }
}
TCR_duplicate_counts_subsampled_inter(VDJ_list, output_directory,batch,CD48)

### TCR v gene usages
TCR_VJ<-function(VDJ_list, output_directory,CD48){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type_broad = VDJ_list[[3]]
  cell_type = VDJ_list[[2]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  V_gene1 = VDJ[,"V_gene1"]
  V_gene2 = VDJ[,"V_gene2"]
  VJ_gene1 = apply(cbind(VDJ[,"V_gene1"], VDJ[,"J_gene1"]), 1, paste, collapse = ":")
  VJ_gene2 = apply(cbind(VDJ[,"V_gene2"], VDJ[,"J_gene2"]), 1, paste, collapse = ":")
  
  V_gene1s = sort(unique(V_gene1))
  V_gene2s = sort(unique(V_gene2))
  VJ_gene1s = sort(unique(VJ_gene1))
  VJ_gene2s = sort(unique(VJ_gene2))
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone = VDJ[,"clone1"]
  w_clone = names(which(clone!='-'))
  
  group_ids = NULL
  for(c in c(1:length(cell_types_broad))){
    w = intersect(w_clone, names(cell_type_broad)[which(cell_type_broad == cell_types_broad[c])])
    group_ids = c(group_ids, list(w))
  }
  for(c in c(1:length(cell_types))){
    w = intersect(w_clone, names(cell_type)[which(cell_type == cell_types[c])])
    group_ids = c(group_ids, list(w))
  }
  group_ids = c(group_ids, list(w_clone),CD48)
  
  names(group_ids) = c(apply(cbind(cell_types_broad,"broad"),1, paste, collapse = ":"), 
                       apply(cbind(cell_types,"detailed"),1, paste, collapse = ":") ,names(CD48) )
  
  list_IGHV_gene_percentage = NULL
  list_IGKLV_gene_percentage = NULL
  list_IGHVJ_gene_percentage = NULL
  list_IGKLVJ_gene_percentage = NULL
  for(c in c(1:length(group_ids))){
    m_IGHV_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(V_gene1s), dimnames = c(list(samples), list(V_gene1s)))
    m_IGKLV_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(V_gene2s), dimnames = c(list(samples), list(V_gene2s)))
    m_IGHVJ_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(VJ_gene1s), dimnames = c(list(samples), list(VJ_gene1s)))
    m_IGKLVJ_gene_percentage = matrix(data = 0,nrow = length(samples), ncol = length(VJ_gene2s), dimnames = c(list(samples), list(VJ_gene2s)))
    
    t = table(sample[w],V_gene1[w])
    m_IGHV_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],V_gene2[w])
    m_IGKLV_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],VJ_gene1[w])
    m_IGHVJ_gene_percentage[rownames(t), colnames(t)] = t
    
    t = table(sample[w],VJ_gene2[w])
    m_IGKLVJ_gene_percentage[rownames(t), colnames(t)] = t
    
    list_IGHV_gene_percentage = c(list_IGHV_gene_percentage, list( m_IGHV_gene_percentage))
    list_IGKLV_gene_percentage = c(list_IGKLV_gene_percentage, list( m_IGKLV_gene_percentage))
    list_IGHVJ_gene_percentage = c(list_IGHVJ_gene_percentage, list( m_IGHVJ_gene_percentage))
    list_IGKLVJ_gene_percentage = c(list_IGKLVJ_gene_percentage, list( m_IGKLVJ_gene_percentage))
  }
  
  names(list_IGHV_gene_percentage) = names(group_ids) 
  names(list_IGKLV_gene_percentage) = names(group_ids) 
  names(list_IGHVJ_gene_percentage ) = names(group_ids) 
  names(list_IGKLVJ_gene_percentage) = names(group_ids) 
  
  list_VDJ_freatures_aggregated  = c(list(list_IGHV_gene_percentage),
                                     list(list_IGKLV_gene_percentage),
                                     list(list_IGHVJ_gene_percentage),
                                     list(list_IGKLVJ_gene_percentage)
  )
  names(list_VDJ_freatures_aggregated)  = c("TRAV gene usage", "TRBV gene usage","TRAVJ gene usage","TRBVJ gene usage")
  
  saveRDS(file = concat(c(output_directory,"VDJ_repertoire_feature_information_",type,"_",batch,".rds")), list_VDJ_freatures_aggregated)
  writeLines(concat(c("\nOutput as RDS object with the following features calculated per cell type:\n",
                      paste(c(" ", names(list_VDJ_freatures_aggregated)), collapse = "\n\t"),
                      "\n\nOutput location:", concat(c(output_directory,"VDJ_repertoire_feature_information_",type,"_",batch,".rds")))))
  
  
  }
TCR_VJ(VDJ_list, output_directory)

print("You can now plot between patient groups using these outputs")
  
################# clonal overlap between phenotypes

### BCR overlap
BCR_clonal_overlap_per_site<-function(VDJ_list, output_directory,groups_PCA){
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
  
  m_cell_types = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
  m_cell_types_proportions = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
  for(s in c(1:length(samples))){
    w = intersect(w_clone , names(which(sample == samples[s])))
    t = table(cell_type[w])
    m_cell_types[samples[s], names(t)]= t
    m_cell_types_proportions[samples[s], names(t)]= t*100/sum(t)
  }
  totals = rowSums(m_cell_types)
  
  ############# get overlap between compartments
  c1 = NULL
  c2 = NULL
  for(i1 in c(1:length(cell_types))){
    for(i2 in c(i1:length(cell_types))){
      if(i1<i2){
        c1 = c(c1, cell_types[i1])
        c2 = c(c2, cell_types[i2])}}}
  sharing_names = apply(cbind(c1,c2), 1, paste,collapse = " - ")
  
  print("Running without subsampling")
  m_overlap_raw = matrix(data = 0,nrow = length(samples), ncol = length(sharing_names), dimnames = c(list(samples), list(sharing_names)))
  
  for(s in c(1:length(samples))){
    print (s)
    w_sample = intersect(w_clone , names(which(sample == samples[s])))
    for(i in c(1:length(c1))){
      w1 = intersect(w_sample, names(which(cell_type==c1[i])))
      w2 = intersect(w_sample, names(which(cell_type==c2[i])))
      shared_clones = intersect(clone[w1],clone[w2])
      m_overlap_raw[samples[s], sharing_names[i]] = length(shared_clones)
    }}
  
  list_overlapping_per_sample = NULL
  for(s in c(1:length(samples))){
    m_overlapping = matrix(data = 0,nrow = length(cell_types), ncol = length(cell_types), dimnames = c(list(cell_types), list(cell_types)))
    m_overlapping_total = matrix(data = 0,nrow = length(cell_types), ncol = length(cell_types), dimnames = c(list(cell_types), list(cell_types)))
    wa = intersect(names(sample)[which(sample== samples[s])], w_clone)
    t = table(clone[wa], cell_type[wa])
    w = which(apply(t, 1, function(x){length(which(x!=0))})>1)
    if(length(w)>0){
      # print(r)
      t=rbind(t[w,])
      for(i in c(1:length(w))){
        x = names(t[i,][which(t[i,]!=0)])
        for(i1 in c(1:(length(x)-1))){
          for(i2 in c((i1+1):length(x))){
            m_overlapping_total[x[i1],x[i2]] = m_overlapping_total[x[i1],x[i2]]+1
          }}
      }
    }
    
    scale = NULL
    if(sum(m_overlapping_total)==0){scale = 1
    m_overlapping_total = m_overlapping_total#/sum(m_overlapping_total)
    }else{m_overlapping_total = m_overlapping_total/sum(m_overlapping_total)}
   
    scale = 1
    list_overlapping_per_sample = c(list_overlapping_per_sample, list(m_overlapping_total*mean(scale)))
    print (s)
  }
  names(list_overlapping_per_sample) = samples
    
  groups = NULL
  links = sharing_names
  
  overlapping_per_link = matrix(data = 0,nrow = length(samples), ncol = length(links), dimnames = c(list(samples), list(links)))
  for(s in c(1:length(samples))){
    mat = list_overlapping_per_sample[[s]]
    for(c1 in c(1:length(cell_types))){
      for(c2 in c(c1:length(cell_types))){
        if(c1<c2){
          if(is.na(mat[cell_types[c1],cell_types[c2]])==F){
            link = concat(c(cell_types[c1], " - ", cell_types[c2]))
            overlapping_per_link[samples[s],link] = mat[cell_types[c1],cell_types[c2]]
          }}}}}
  
  m_overlap_sampled1 = overlapping_per_link
  saveRDS(file = concat(c(output_directory,"Clonal_overlap_cell_types_", batch,"_", type1,".RDS")), m_overlap_sampled1)
  
  
  #################### plot by PCA group
  m_overlap_sampled1  = readRDS(file = concat(c(output_directory,"Clonal_overlap_cell_types_", batch,"_", type1,".RDS")))
  
  Means_factor = function(factor, x){
    m = NULL
    for(i1 in c(1:length(levels(factor)))){
      x1 = x[which(factor==levels(factor)[i1])]
      x1 = x1[which(x1!=-1)]
      m = c(m, mean(x1))}
    return(m)}
  
  mat = m_overlap_sampled1
  nz = apply(mat,2,function(x){length(which(x!=0))})
  mat = mat[,which(nz>=1)]
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  w = which(group_PCA_list %in% rownames(mat))
  mat_stat = mat[group_PCA_list[w],]
  
  factor = factor(factor_PCA[w])
  mat_stat1 = mat_stat^0.25
  fit = manova(formula = mat_stat1 ~ factor)
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
  length(which(p_value<0.05))
  
  colnames(means) = paste("mean.group.", c(1:length(means[1,])))
  combined_p_value = cbind(p_value ,means)
  name = colnames(mat_stat)
  rownames(combined_p_value) = nam
  p.site = rep(concat(c("biopsy")), length(nam))
  analysis = "overlap"
  p.analysis = rep(analysis, length(nam))
  x = cbind(p.site, p.analysis, combined_p_value)
  summary_stats = x
  
  
  groups = NULL
  for(i in c(1:length(mat_stat[1,]))){
    g1 = NULL
    for(s in c(1:length(groups_PCA))){
      g1 = c(g1, list(mat_stat[ gsub("_","-",groups_PCA[[s]]),i]))}
    groups = c(groups, list(g1))
  }
  names(groups) = colnames(mat_stat)
  
  ######################################
  analysis = "Clonal_overlap_B_cell_subsets"
  
  outfile = concat(c(output_directory,"Stats_",analysis,"_biopsy_", batch,".txt"))
  
  write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  summary_stats_biopsy = summary_stats
  
  fileout1=concat(c(output_directory,"",analysis,"_biopsy_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*1.8, width=w*3.7)
  par(mfrow= c(1,1), mar = c(20,4,2.5,1))
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  p_values = p_value
  factors1 = paste("group",c(1:length(groups_PCA)))
  factors = names(groups)
  main = concat(c("clonal overlap:biopsy"))
  max = max(c(unlist(groups), unlist(groups))*1.3)
  min = 0
  b = (max-min)*0.035
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
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
  plot(c(1.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
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
  for(i in c(1:length(groups))){	
    b = max*0.035
    signif_threshold = 0.05
    if(p_values[i]<signif_threshold){
      pval1 = "*"
      y = max(unlist(groups[[i]]))+b
      # if(p_values[i] <signif_threshold/10){pval1 = "**"}
      # if(p_values[i] <signif_threshold/100){pval1 = "***"}
      # segments(1,y+b, i+1,y+b,lwd = 3, col = "darkgrey")
      text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.4)
    }
  }
  
  
  dev.off()
  
  #################### plot by PCA group blood
  
  mat = m_overlap_sampled1
  group_PCA_list = NULL
  factor_PCA = NULL
  for(i in c(1:length(groups_PCA))){
    group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA[[i]]))
    factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
  }
  group_PCA_list = gsub("biopsy","blood",group_PCA_list)
  w = which(group_PCA_list %in% rownames(mat))
  mat_stat = mat[group_PCA_list[w],]
  nz = apply(mat_stat ,2,function(x){length(which(x!=0))})
  mat_stat  = mat_stat [,which(nz>=1)]
  
  factor = factor(factor_PCA[w])
  mat_stat1 = mat_stat^0.25
  fit = manova(formula = mat_stat1 ~ factor)
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
  name = colnames(mat_stat)
  rownames(combined_p_value) = nam
  p.site = rep(concat(c("biopsy")), length(nam))
  p.analysis = rep(analysis, length(nam))
  x = cbind(p.site, p.analysis, combined_p_value)
  summary_stats = x
  
  analysis = "Clonal_overlap_B_cell_subsets"
  
  outfile = concat(c(output_directory,"Stats_",analysis,"_blood_", batch,".txt"))
  
  write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  groups = NULL
  for(i in c(1:length(mat_stat[1,]))){
    g1 = NULL
    for(s in c(1:length(groups_PCA))){
      g1 = c(g1, list(mat_stat[gsub("_biopsy","-blood",groups_PCA[[s]]),i]))}
    groups = c(groups, list(g1))
  }
  names(groups) = colnames(mat_stat)
  
  ######################################
  analysis = "Clonal_overlap_B_cell_subsets"
  fileout1=concat(c(output_directory,"",analysis,"_blood_", batch,".pdf"))
  w=3.2
  pdf(file=fileout1, height=w*1.8, width=w*3.7)
  par(mfrow= c(1,1), mar = c(20,4,2.5,1))
  p_values = p_value
  factors1 = paste("group",c(1:length(groups_PCA)))
  factors = names(groups)
  main = concat(c("clonal overlap:blood"))
  max = max(c(unlist(groups), unlist(groups))*1.2)
  min = 0
  b = (max-min)*0.035
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
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
    for(i1 in c(1:length(groups[[i]]))){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =21, bg=cols[i1],col=cols[i1], cex = 0.4)
    }
  }
  
  y = max(unlist(groups))
  for(i in c(1:length(groups))){	
    b = max*0.035
    signif_threshold = 0.05
    if(p_values[i]<signif_threshold){
      pval1 = "*"
      y = max(unlist(groups[[i]]))+b
      # if(p_values[i] <signif_threshold/10){pval1 = "**"}
      # if(p_values[i] <signif_threshold/100){pval1 = "***"}
      # segments(1,y+b, i+1,y+b,lwd = 3, col = "darkgrey")
      text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.4)
    }
  }
  
  
  dev.off()
  
  ###################### plot network plot for clonal sharing
  pca_groups = paste("Group", c(1:length(groups_PCA)))
  mean_cell_type_proportions = matrix(data = 0,nrow = length(pca_groups), ncol = length(cell_types), dimnames = c(list(pca_groups), list(cell_types)))
  
  mean_cell_overlap = matrix(data = 0,nrow = length(pca_groups), ncol = length(colnames(m_overlap_sampled1)), dimnames = c(list(pca_groups), list(colnames(m_overlap_sampled1))))
  
  for(i in c(1:length(pca_groups))){
    mean_cell_type_proportions [i,]= apply(m_cell_types_proportions[gsub("_","-",groups_PCA[[i]]), ], 2, median)
    mean_cell_overlap [i,] = apply(m_overlap_sampled1[gsub("_","-",groups_PCA[[i]]), ], 2, mean)
  }
  
  
  ### 
  analysis = "Clonal_overlap_B_cell_subsets"
  file = concat(c(output_directory,"Stats_",analysis,"_biopsy_", batch,".txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  p_value = as.numeric(p[,"p_value"])
  transition = rownames(p)
  names(p_value) = transition
  mean1 = as.numeric(p[,"mean.group..1"])
  mean2 = as.numeric(p[,"mean.group..2"])
  names(mean1) = transition
  names(mean2) = transition
  pca_group_means = c(list(mean1), list(mean2))
  transition_matrix = cbind(mean1, mean2)
  transition_split=strsplit(transition, " - ", fixed = T)
  transition1 = NULL
  transition2 = NULL
  for(i in c(1:length(transition_split))){
    transition1 = c(transition1, transition_split[[i]][[1]])
    transition2 = c(transition2, transition_split[[i]][[2]])
  }
  all_cell_types = cell_types
  
  ##################
  library(igraph)
  
  signif = which(p_value<0.05)
  
  fileout1=concat(c(output_directory,"",analysis,"_network_biopsy_", batch,".pdf"))
  w=2.7
  pdf(file=fileout1, height=w*1, width=w*3)
  par(mfrow= c(1,3), mar = c(1,1,3,1))
  
  c = 1
  mean_edge_strength = (pca_group_means[[c]]+pca_group_means[[c+1]])/2
  main = concat(c("immuno-group ",c))
  g <- graph.empty( n=0, directed=FALSE)
  names_plot = all_cell_types
  names_plot = gsub("B cell  ","",names_plot)
  names_plot = gsub("B cell ","",names_plot)
  names_plot = gsub("activated pre-","activated\npre-",names_plot)
  names_plot = gsub("memory activated","memory\nactivated",names_plot)
  
   g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
  names <- V(g)$name
  ids <- 1:length(names)
  names(ids) <- names
  
  w = which(mean_edge_strength!=0)
  w = c(1:length(mean_edge_strength))
  edge_strength = mean_edge_strength[w]
  p_value_sub = p_value[w]
  from <- transition1[w]
  to <- transition2[w]
  w = intersect(which(from %in% names),which(to %in% names))
  p_value_sub = p_value_sub[w]
  edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
  g <- add.edges(g, t(edges), weight= edge_strength[w])
  signif_edges = names(which(p_value_sub<0.05))
  t = transition_matrix[signif_edges,]
  
  sizes = mean_cell_type_proportions[c, all_cell_types]
  sizes_scaled = sizes^0.3
  sizes_scaled = sizes_scaled*5
  
  V(g)$size<-sizes_scaled
  V(g)$label.cex<-0.5
  V(g)$name = names_plot
  V(g)$color = "grey"
    layout1 =layout_in_circle(g)
    
    edge_strength_plot = edge_strength^0.5
    edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
    
    col = rgb(0.5,0.5,0.5,alpha = 0.5)
    cols = rep(col, length(edge_strength_plot))
    names(cols) = names(p_value_sub)
    cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, min))))] = add.alpha("blue",alpha = 0.65)
    cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, max))))] = add.alpha("red",alpha = 0.65)
    
    # col = "black"
    plot(g, layout=layout1, edge.color= cols, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.00000001, edge.arrow.size = 1, edge.lty = 1,xlim = range(layout1[,1]*1.5),ylim = range(layout1[,2]*1.5))
    
    ## Apply labels manually
    #Specify x and y coordinates of labels, adjust outward as desired
    x = layout1[,1]*1.5
    y = layout1[,2]*1.5
    
    #create vector of angles for text based on number of nodes 
    # (flipping the orientation of the words half way around so none appear 
    # upside down)
    angle = ifelse(atan(-(layout1[,1]/layout1[,2]))*(180/pi) < 0,  
                   90 + atan(- (layout1[,1]/layout1[,2]))*(180/pi), 270 + atan(-layout1[,1]/layout1[,2])*(180/pi))
    
    #Apply the text labels with a loop with angle as srt
    for (i in 1:length(x)) {
      text(x=x[i], y=y[i], labels=V(g)$name[i], adj=NULL, 
           pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
    }
    dev.off()
    
    
    ####### for blood
    ### 
    analysis = "Clonal_overlap_B_cell_subsets"
    file = concat(c(output_directory,"Stats_",analysis,"_blood_", batch,".txt"))
    p <- as.matrix(read.csv(file, head=T, sep="\t"))
    p_value = as.numeric(p[,"p_value"])
    transition = rownames(p)
    names(p_value) = transition
    mean1 = as.numeric(p[,"mean.group..1"])
    mean2 = as.numeric(p[,"mean.group..2"])
    names(mean1) = transition
    names(mean2) = transition
    pca_group_means = c(list(mean1), list(mean2))
    transition_matrix = cbind(mean1, mean2)
    transition_split=strsplit(transition, " - ", fixed = T)
    transition1 = NULL
    transition2 = NULL
    for(i in c(1:length(transition_split))){
      transition1 = c(transition1, transition_split[[i]][[1]])
      transition2 = c(transition2, transition_split[[i]][[2]])
    }
    all_cell_types = cell_types
    
    ##################
    library(igraph)
    
    signif = which(p_value<0.05)
    
    fileout1=concat(c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_PDAC150K/Plots new annotations/BCR_TCR/",analysis,"_network_blood_", batch,".pdf"))
    w=2.7
    pdf(file=fileout1, height=w*1, width=w*3)
    par(mfrow= c(1,3), mar = c(1,1,3,1))
    
    c = 1
    mean_edge_strength = (pca_group_means[[c]]+pca_group_means[[c+1]])/2
    main = concat(c("immuno-group ",c))
    g <- graph.empty( n=0, directed=FALSE)
    names_plot = all_cell_types
    names_plot = gsub("B cell  ","",names_plot)
    names_plot = gsub("B cell ","",names_plot)
    names_plot = gsub("activated pre-","activated\npre-",names_plot)
    names_plot = gsub("memory activated","memory\nactivated",names_plot)
    g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
    names <- V(g)$name
    ids <- 1:length(names)
    names(ids) <- names
    
    w = which(mean_edge_strength!=0)
    w = c(1:length(mean_edge_strength))
    edge_strength = mean_edge_strength[w]
    p_value_sub = p_value[w]
    from <- transition1[w]
    to <- transition2[w]
    w = intersect(which(from %in% names),which(to %in% names))
    p_value_sub = p_value_sub[w]
    edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
    g <- add.edges(g, t(edges), weight= edge_strength[w])
    signif_edges = names(which(p_value_sub<0.05))
    t = transition_matrix[signif_edges,]
    
    sizes = mean_cell_type_proportions[c, all_cell_types]
    sizes_scaled = sizes^0.3
    sizes_scaled = sizes_scaled*5
    
    V(g)$size<-sizes_scaled
    V(g)$label.cex<-0.5
    V(g)$name = names_plot
    V(g)$color = "grey"
      layout1 =layout_in_circle(g)
      edge_strength_plot = edge_strength^0.5
      edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
      
      col = rgb(0.5,0.5,0.5,alpha = 0.5)
      cols = rep(col, length(edge_strength_plot))
      names(cols) = names(p_value_sub)
      if(length(signif_edges)<1){
        cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, min))))] = add.alpha("blue",alpha = 0.65)
        cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, max))))] = add.alpha("red",alpha = 0.65)
      }
      
      plot(g, layout=layout1, edge.color= cols, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.00000001, edge.arrow.size = 1, edge.lty = 1,xlim = range(layout1[,1]*1.5),ylim = range(layout1[,2]*1.5))
      
      ## Apply labels manually
      #Specify x and y coordinates of labels, adjust outward as desired
      x = layout1[,1]*1.5
      y = layout1[,2]*1.5
      
      #create vector of angles for text based on number of nodes 
      # (flipping the orientation of the words half way around so none appear 
      # upside down)
      angle = ifelse(atan(-(layout1[,1]/layout1[,2]))*(180/pi) < 0,  
                     90 + atan(- (layout1[,1]/layout1[,2]))*(180/pi), 270 + atan(-layout1[,1]/layout1[,2])*(180/pi))
      
      #Apply the text labels with a loop with angle as srt
      for (i in 1:length(x)) {
        text(x=x[i], y=y[i], labels=V(g)$name[i], adj=NULL, 
             pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
      }
      dev.off()
}
BCR_clonal_overlap_per_site(VDJ_list, output_directory)

### TCR overlap
TCR_clonal_overlap_per_site<-function(VDJ_list, output_directory,CD48, groups_PCA){
  VDJ_list = VDJ_list_TCR
  type = "TCR"
  type1 = "T_cells"
  
  VDJ = VDJ_list[[1]]
  cell_type = VDJ_list[[2]]
  cell_type_broad = VDJ_list[[3]]
  sample = VDJ_list[[4]]
  Patient = VDJ_list[[5]]
  Sample.Type = VDJ_list[[5]]
  
  samples = sort(unique(sample))
  Patients = sort(unique(Patient))
  Sample.Types = sort(unique(Sample.Type))
  cell_types = sort(unique(cell_type))
  cell_types_broad = sort(unique(cell_type_broad))
  clone1 = VDJ[,"clone1"]
  clone2 = VDJ[,"clone2"]
  clone = clone1
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type1 = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    w_cd48 = names(cell_type)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[w_cd48])
    w_clone = intersect(names(which(clone!='-')),w_cd48)
    chain = names(CD48)[cd48]
    cell_types1 = cell_types_use

    ## cell counts per cell type
    w_clone = intersect(names(which(clone1!='-')),names(which(clone2!='-')))
    w_clone = intersect(w_clone,w_cd48)
    cell_types = sort(unique(cell_type[w_clone]))
    m_cell_types = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
    m_cell_types_proportions = matrix(data = 0,nrow = length(samples), ncol = length(cell_types), dimnames = c(list(samples), list(cell_types)))
    for(s in c(1:length(samples))){
      w = intersect(w_clone , names(which(sample == samples[s])))
      t = table(cell_type[w])
      m_cell_types[samples[s], names(t)]= t
      m_cell_types_proportions[samples[s], names(t)]= t*100/sum(t)
    }
    totals = rowSums(m_cell_types)
    
    ############# get overlap between compartments
    c1 = NULL
    c2 = NULL
    for(i1 in c(1:length(cell_types1))){
      for(i2 in c(i1:length(cell_types1))){
        if(i1<i2){
          c1 = c(c1, cell_types1[i1])
          c2 = c(c2, cell_types1[i2])}}}
    sharing_names = apply(cbind(c1,c2), 1, paste,collapse = " - ")
    
    print("Running without subsampling")
    m_overlap_raw = matrix(data = 0,nrow = length(samples), ncol = length(sharing_names), dimnames = c(list(samples), list(sharing_names)))
    
    for(s in c(1:length(samples))){
      print(s)
      w_sample = intersect(w_clone , names(which(sample == samples[s])))
      t = table(clone[w_sample], cell_type[w_sample])
      a = apply(t, 1, function(x){length(which(x!=0))})
      t[which(a>1),]
      for(i in c(1:length(c1))){
        w1 = intersect(w_sample, names(which(cell_type==c1[i])))
        w2 = intersect(w_sample, names(which(cell_type==c2[i])))
        shared_clones = intersect(clone[w1],clone[w2])
        m_overlap_raw[samples[s], sharing_names[i]] = length(shared_clones)
      }
    }
    
    list_overlapping_per_sample = NULL
    for(s in c(1:length(samples))){
      m_overlapping = matrix(data = 0,nrow = length(cell_types), ncol = length(cell_types), dimnames = c(list(cell_types), list(cell_types)))
      m_overlapping_total = matrix(data = 0,nrow = length(cell_types), ncol = length(cell_types), dimnames = c(list(cell_types), list(cell_types)))
      wa = intersect(names(sample)[which(sample== samples[s])], w_clone)
      t = table(clone[wa], cell_type[wa])
      w = which(apply(t, 1, function(x){length(which(x!=0))})>1)
      if(length(w)>0){
        # print(r)
        t=rbind(t[w,])
        for(i in c(1:length(w))){
          x = names(t[i,][which(t[i,]!=0)])
          for(i1 in c(1:(length(x)-1))){
            for(i2 in c((i1+1):length(x))){
              m_overlapping_total[x[i1],x[i2]] = m_overlapping_total[x[i1],x[i2]]+1
            }}
        }
      }
      
      scale = NULL
      if(sum(m_overlapping_total)==0){scale = 1
      m_overlapping_total = m_overlapping_total#/sum(m_overlapping_total)
      }else{m_overlapping_total = m_overlapping_total/sum(m_overlapping_total)}
      
      scale = 1
      list_overlapping_per_sample = c(list_overlapping_per_sample, list(m_overlapping_total*mean(scale)))
      print (s)
    }
    names(list_overlapping_per_sample) = samples
    
    groups = NULL
    links = sharing_names
    
    overlapping_per_link = matrix(data = 0,nrow = length(samples), ncol = length(links), dimnames = c(list(samples), list(links)))
    for(s in c(1:length(samples))){
      mat = list_overlapping_per_sample[[s]]
      for(c1 in c(1:length(cell_types))){
        for(c2 in c(c1:length(cell_types))){
          if(c1<c2){
            if(is.na(mat[cell_types[c1],cell_types[c2]])==F){
              link = concat(c(cell_types[c1], " - ", cell_types[c2]))
              overlapping_per_link[samples[s],link] = mat[cell_types[c1],cell_types[c2]]
            }}}}}
    
    m_overlap_sampled1 = overlapping_per_link
    saveRDS(file = concat(c(output_directory,"Clonal_overlap_cell_types_", batch,"_", type1,".RDS")), m_overlap_sampled1)
    print(concat(c(type1, " done")))
  }
  
  print("Plots and statistics running")
  
  #################### plot by PCA group
  for(cd48 in c(1:length(CD48))){
    cell_types_use = CD48[[cd48]]
    type1 = concat(c("TCR_",names(CD48)[cd48]))
    print(concat(c("Starting ", type)))
    w_cd48 = names(sample)[sort(unique(c(which(cell_type %in% cell_types_use), which(cell_type_broad %in% cell_types_use))))]
    ## check table(cell_type[wx])
    w_clone = intersect(names(which(clone!='-')),w_cd48)
    chain = names(CD48)[cd48]
    cell_types1 = cell_types_use
    
    
    m_overlap_sampled1  =readRDS(file = concat(c(output_directory,"Clonal_overlap_cell_types_", batch,"_", type1,".RDS")))
    
    Means_factor = function(factor, x){
      m = NULL
      for(i1 in c(1:length(levels(factor)))){
        x1 = x[which(factor==levels(factor)[i1])]
        x1 = x1[which(x1!=-1)]
        m = c(m, mean(x1))}
      return(m)}
    
    mat = m_overlap_sampled1
    nz = apply(mat,2,function(x){length(which(x!=0))})
    mat = mat[,which(nz>=1)]
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, gsub("_","-",groups_PCA[[i]]))
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    w = which(group_PCA_list %in% rownames(mat))
    mat_stat = mat[group_PCA_list[w],]
    
    factor = factor(factor_PCA[w])
    mat_stat1 = mat_stat^1
    fit = manova(formula = mat_stat1 ~ factor)
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
    names(p_value) =nam
    print(head(sort(p_value), 10))
    analysis= concat(c("Clonal_overlap_T_cell_subsets_",type1))
    colnames(means) = paste("mean.group.", c(1:length(means[1,])))
    combined_p_value = cbind(p_value ,means)
    name = colnames(mat_stat)
    rownames(combined_p_value) = nam
    p.site = rep(concat(c("biopsy")), length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    summary_stats = x
    
    
    groups = NULL
    for(i in c(1:length(mat_stat[1,]))){
      g1 = NULL
      for(s in c(1:length(groups_PCA))){
        g1 = c(g1, list(mat_stat[gsub("_","-",groups_PCA[[s]]),i]))}
      groups = c(groups, list(g1))
    }
    names(groups) = colnames(mat_stat)
    
    ######################################
    analysis = "Clonal_overlap_T_cell_subsets"
    
    outfile = concat(c(output_directory,"Stats_",analysis,"_biopsy_", type1,"_", batch,".txt"))
    
    write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    summary_stats_biopsy = summary_stats
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(2,8)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    fileout1=concat(c(output_directory,"",analysis,"_biopsy_", type1,"_", batch,".pdf"))
    w=3.2
    pdf(file=fileout1, height=w*1.8, width=w*6.7)
    par(mfrow= c(1,1), mar = c(20,4,2.5,1))
    p_values = p_value
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = names(groups)
    main = concat(c("clonal overlap:biopsy"))
    max = max(c(unlist(groups), unlist(groups))*1.2)
    min = 0
    b = (max-min)*0.035
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 80#length(groups)
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
    plot(c(3, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
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
    for(i in c(1:length(groups))){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        y = max(unlist(groups[[i]]))+b
        # if(p_values[i] <signif_threshold/10){pval1 = "**"}
        # if(p_values[i] <signif_threshold/100){pval1 = "***"}
        # segments(1,y+b, i+1,y+b,lwd = 3, col = "darkgrey")
        text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.4)
      }
    }
    
    
    dev.off()
    
    #################### plot by PCA group blood
    
    mat = m_overlap_sampled1
    group_PCA_list = NULL
    factor_PCA = NULL
    for(i in c(1:length(groups_PCA))){
      group_PCA_list = c(group_PCA_list, groups_PCA[[i]])
      factor_PCA = c(factor_PCA, rep(i, length(groups_PCA[[i]])))
    }
    group_PCA_list = gsub("_biopsy","-blood",group_PCA_list)
    w = which(group_PCA_list %in% rownames(mat))
    mat_stat = mat[group_PCA_list[w],]
    nz = apply(mat_stat ,2,function(x){length(which(x!=0))})
    mat_stat  = mat_stat [,which(nz>=2)]
    
    factor = factor(factor_PCA[w])
    mat_stat1 = mat_stat^1
    fit = manova(formula = mat_stat1 ~ factor)
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
    name = colnames(mat_stat)
    rownames(combined_p_value) = nam
    p.site = rep(concat(c("biopsy")), length(nam))
    p.analysis = rep(analysis, length(nam))
    x = cbind(p.site, p.analysis, combined_p_value)
    summary_stats = x
    
    analysis = "Clonal_overlap_T_cell_subsets"
    
    outfile = concat(c(output_directory,"Stats_",analysis,"_blood_", type1,"_", batch,".txt"))
    
    write.table(summary_stats, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    
    
    groups = NULL
    for(i in c(1:length(mat_stat[1,]))){
      g1 = NULL
      for(s in c(1:length(groups_PCA))){
        g1 = c(g1, list(mat_stat[gsub("_biopsy","-blood",groups_PCA[[s]]),i]))}
      groups = c(groups, list(g1))
    }
    names(groups) = colnames(mat_stat)
    
    ######################################
    analysis = "Clonal_overlap_T_cell_subsets"
    fileout1=concat(c(output_directory,"",analysis,"_blood_", type1,"_", batch,".pdf"))
    w=3.2
    pdf(file=fileout1, height=w*1.8, width=w*8.7)
    par(mfrow= c(1,1), mar = c(20,4,2.5,1))
    p_values = p_value
    factors1 = paste("group",c(1:length(groups_PCA)))
    factors = names(groups)
    main = concat(c("clonal overlap:blood"))
    max = max(c(unlist(groups), unlist(groups))*1.2)
    min = 0
    b = (max-min)*0.035
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 60#length(groups)
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
    plot(c(3, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
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
    for(i in c(1:length(groups))){	
      b = max*0.035
      signif_threshold = 0.05
      if(p_values[i]<signif_threshold){
        pval1 = "*"
        y = max(unlist(groups[[i]]))+b
        # if(p_values[i] <signif_threshold/10){pval1 = "**"}
        # if(p_values[i] <signif_threshold/100){pval1 = "***"}
        # segments(1,y+b, i+1,y+b,lwd = 3, col = "darkgrey")
        text(mean(c(i)), y+1.5*b, labels = pval1, cex = 1.7)
      }
    }
    
    
    dev.off()
    
    ###################### plot network plot for clonal sharing
    pca_groups = paste("Group", c(1:length(groups_PCA)))
    mean_cell_type_proportions = matrix(data = 0,nrow = length(pca_groups), ncol = length(cell_types1), dimnames = c(list(pca_groups), list(cell_types1)))
    for(i in c(1:length(pca_groups))){
      mean_cell_type_proportions [i,colnames(m_cell_types_proportions)]= apply(m_cell_types_proportions[gsub("_","-",groups_PCA[[i]]), ], 2, median)
    }
    
    
    ### 
    analysis = "Clonal_overlap_T_cell_subsets"
    file = concat(c(output_directory,"Stats_",analysis,"_biopsy_", type1,"_", batch,".txt"))
    p <- as.matrix(read.csv(file, head=T, sep="\t"))
    p_value = as.numeric(p[,"p_value"])
    transition = rownames(p)
    names(p_value) = transition
    mean1 = as.numeric(p[,"mean.group..1"])
    mean2 = as.numeric(p[,"mean.group..2"])
    names(mean1) = transition
    names(mean2) = transition
    pca_group_means = c(list(mean1), list(mean2))
    transition_matrix = cbind(mean1, mean2)
    transition_split=strsplit(transition, " - ", fixed = T)
    transition1 = NULL
    transition2 = NULL
    for(i in c(1:length(transition_split))){
      transition1 = c(transition1, transition_split[[i]][[1]])
      transition2 = c(transition2, transition_split[[i]][[2]])
    }
    all_cell_types = unique(c(transition1,transition2))
    
    ##################
    library(igraph)
    
    signif = which(p_value<0.05)
    
    c = 1
    mean_edge_strength = (pca_group_means[[c]]+pca_group_means[[c+1]])/2
    main = concat(c("immuno-group ",c))
    g <- graph.empty( n=0, directed=FALSE)
    if(cd48==1){
      all_cell_types_broad = gsub("Activated ", "", all_cell_types)
      order = order(all_cell_types_broad)
      all_cell_types = all_cell_types[order]}
    # order= c(grep("CD4", all_cell_types),grep("Treg", all_cell_types),grep("Tfh", all_cell_types)
    # order = intersect(cell_classifications [which(cell_classifications[,1]=="CD4"), 2], all_cell_types)
    # all_cell_types = c(order , setdiff(all_cell_types, order))
    # all_cell_types = setdiff(all_cell_types, "MAIT")
    names_plot = all_cell_types
    names_plot = gsub("_"," ",names_plot)
    names_plot = gsub(" EffectorMem","\nEffectorMem",names_plot)
    names_plot = gsub("Exhausted ","Exhausted\n",names_plot)
    names_plot = gsub("chemokine ","chemokine\n",names_plot)
    
    g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
    names <- V(g)$name
    ids <- 1:length(names)
    names(ids) <- names
    
    w = which(mean_edge_strength!=0)
    w = c(1:length(mean_edge_strength))
    edge_strength = mean_edge_strength[w]
    p_value_sub = p_value[w]
    from <- transition1[w]
    to <- transition2[w]
    w = intersect(which(from %in% names),which(to %in% names))
    p_value_sub = p_value_sub[w]
    edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
    g <- add.edges(g, t(edges), weight= edge_strength[w])
    signif_edges = names(which(p_value_sub<0.05))
    t = transition_matrix[signif_edges,]
    
    sizes = mean_cell_type_proportions[c, all_cell_types]
    sizes_scaled = sizes^0.3
    sizes_scaled = sizes_scaled*5
    
    V(g)$size<-sizes_scaled
    V(g)$label.cex<-0.5
    V(g)$name = names_plot
    V(g)$color = "grey"
      layout1 =layout_in_circle(g)
      edge_strength_plot = edge_strength[w]^0.5
      edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
      
      col = rgb(0.5,0.5,0.5,alpha = 0.5)
      cols = rep(col, length(edge_strength_plot))
      names(cols) = names(p_value_sub)
      c = 1
      cols [intersect(signif_edges, names(which(transition_matrix[,c]==apply(transition_matrix, 1, min))))] = add.alpha("blue",alpha = 0.65)
      cols [intersect(signif_edges, names(which(transition_matrix[,c]==apply(transition_matrix, 1, max))))] = add.alpha("red",alpha = 0.65)
      #cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, min))))] = add.alpha("purple",alpha = 0.65)
      
      fileout1=concat(c(output_directory,"",analysis,"_network_biopsy_", type1,"_", batch,".pdf"))
      w=2.7
      pdf(file=fileout1, height=w*1, width=w*3)
      par(mfrow= c(1,3), mar = c(2,1,2,1))
      # col = "black"
      plot(g, layout=layout1, edge.color= cols, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.00000001, edge.arrow.size = 1, edge.lty = 1,xlim = range(layout1[,1]*1.7),ylim = range(layout1[,2]*1.7))
      
      ## Apply labels manually
      #Specify x and y coordinates of labels, adjust outward as desired
      x = layout1[,1]*1.65
      y = layout1[,2]*1.65
      
      #create vector of angles for text based on number of nodes 
      # (flipping the orientation of the words half way around so none appear 
      # upside down)
      angle = ifelse(atan(-(layout1[,1]/layout1[,2]))*(180/pi) < 0,  
                     90 + atan(- (layout1[,1]/layout1[,2]))*(180/pi), 270 + atan(-layout1[,1]/layout1[,2])*(180/pi))
      
      #Apply the text labels with a loop with angle as srt
      lab = gsub("T cell ","",V(g)$name)
      lab = gsub("Activated ","Activated\n", lab)
      lab = gsub("NK cell ","", lab)
      lab = gsub("CD4","CD4\n", lab)
      lab = gsub("CD8","CD8\n", lab)
      
      for (i in 1:length(x)) {
        text(x=x[i], y=y[i], labels=lab[i], adj=NULL, 
             pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
      }
      dev.off()
      
      
      ####### for blood
      ### 
      analysis = "Clonal_overlap_T_cell_subsets"
      file = concat(c(output_directory,"Stats_",analysis,"_blood_", type1,"_", batch,".txt"))
      p <- as.matrix(read.csv(file, head=T, sep="\t"))
      p_value = as.numeric(p[,"p_value"])
      transition = rownames(p)
      names(p_value) = transition
      mean1 = as.numeric(p[,"mean.group..1"])
      mean2 = as.numeric(p[,"mean.group..2"])
      names(mean1) = transition
      names(mean2) = transition
      pca_group_means = c(list(mean1), list(mean2))
      transition_matrix = cbind(mean1, mean2)
      transition_split=strsplit(transition, " - ", fixed = T)
      transition1 = NULL
      transition2 = NULL
      for(i in c(1:length(transition_split))){
        transition1 = c(transition1, transition_split[[i]][[1]])
        transition2 = c(transition2, transition_split[[i]][[2]])
      }
      all_cell_types = unique(c(transition1,transition2))
      
      ##################
      library(igraph)
      
      signif = which(p_value<0.05)
      if(cd48==1){
        all_cell_types_broad = gsub("Activated ", "", all_cell_types)
        order = order(all_cell_types_broad)
        all_cell_types = all_cell_types[order]}
      
      
      c = 1
      mean_edge_strength = (pca_group_means[[c]]+pca_group_means[[c+1]])/2
      main = concat(c("immuno-group ",c))
      g <- graph.empty( n=0, directed=FALSE)
      names_plot = all_cell_types
      names_plot = gsub("_"," ",names_plot)
      names_plot = gsub(" EffectorMem","\nEffectorMem",names_plot)
      names_plot = gsub("Exhausted ","Exhausted\n",names_plot)
      names_plot = gsub("chemokine ","chemokine\n",names_plot)
      g <- igraph::add.vertices(g, length(all_cell_types), name= all_cell_types) 
      names <- V(g)$name
      ids <- 1:length(names)
      names(ids) <- names
      
      w = which(mean_edge_strength!=0)
      w = c(1:length(mean_edge_strength))
      edge_strength = mean_edge_strength[w]
      p_value_sub = p_value[w]
      from <- transition1[w]
      to <- transition2[w]
      w = intersect(which(from %in% names),which(to %in% names))
      p_value_sub = p_value_sub[w]
      edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
      g <- add.edges(g, t(edges), weight= edge_strength[w])
      signif_edges = names(which(p_value_sub<0.05))
      t = transition_matrix[signif_edges,]
      
      sizes = mean_cell_type_proportions[c, all_cell_types]
      sizes_scaled = sizes^0.3
      sizes_scaled = sizes_scaled*5
      
      V(g)$size<-sizes_scaled
      V(g)$label.cex<-0.5
      V(g)$name = names_plot
      V(g)$color = "grey"
        layout1 =layout_in_circle(g)
        edge_strength_plot = edge_strength[w]^0.5
        edge_strength_plot = (edge_strength_plot*7/max(edge_strength_plot))+0.05
        
        col = rgb(0.5,0.5,0.5,alpha = 0.5)
        cols = rep(col, length(edge_strength_plot))
        names(cols) = names(p_value_sub)
        c = 1
        cols [intersect(signif_edges, names(which(transition_matrix[,c]==apply(transition_matrix, 1, min))))] = add.alpha("blue",alpha = 0.65)
        cols [intersect(signif_edges, names(which(transition_matrix[,c]==apply(transition_matrix, 1, max))))] = add.alpha("red",alpha = 0.65)
        #cols [intersect(signif_edges, names(which(t[,c]==apply(t, 1, min))))] = add.alpha("purple",alpha = 0.65)
        
        fileout1=concat(c(output_directory,"",analysis,"_network_blood_", type1,"_", batch,".pdf"))
        w=2.7
        pdf(file=fileout1, height=w*1, width=w*3)
        par(mfrow= c(1,3), mar = c(2,1,2,1))
        # col = "black"
        plot(g, layout=layout1, edge.color= cols, main="", edge.width= edge_strength_plot,vertex.label.family="sans", vertex.label.cex = 0.00000001, edge.arrow.size = 1, edge.lty = 1,xlim = range(layout1[,1]*1.7),ylim = range(layout1[,2]*1.7))
        
        ## Apply labels manually
        #Specify x and y coordinates of labels, adjust outward as desired
        x = layout1[,1]*1.65
        y = layout1[,2]*1.65
        
        #create vector of angles for text based on number of nodes 
        # (flipping the orientation of the words half way around so none appear 
        # upside down)
        angle = ifelse(atan(-(layout1[,1]/layout1[,2]))*(180/pi) < 0,  
                       90 + atan(- (layout1[,1]/layout1[,2]))*(180/pi), 270 + atan(-layout1[,1]/layout1[,2])*(180/pi))
        
        #Apply the text labels with a loop with angle as srt
        lab = gsub("T cell ","",V(g)$name)
        lab = gsub("Activated ","Activated\n", lab)
        lab = gsub("NK cell ","", lab)
        
        for (i in 1:length(x)) {
          text(x=x[i], y=y[i], labels=lab[i], adj=NULL, 
               pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
        }
        dev.off()
  }
  
}
TCR_clonal_overlap_per_site(VDJ_list, output_directory,CD48, groups_PCA)





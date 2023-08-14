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

file = concat(c(input_directory,"Clinical information PDAC150K.txt"))
p <- as.matrix(read.csv(file, head=T, sep="\t"))

sample = p[,1]
df = data.frame(sample)
df$source = as.factor(p[,"Source"])
df$Sex = as.factor(p[,"Sex"])
df$Age = as.numeric(p[,"Age"])
df$BMI = as.numeric(p[,"BMI"])
df$Past.medical.history = as.factor(p[,"Past.medical.history"])
df$Group = as.factor(p[,"Group"])

### check continuous variables
mat = cbind(df[,c("Age")],df[,"BMI"])
fit = manova(formula = mat ~ df$Group)
p1 = summary.aov(fit)
p1
# not significant

### check catagorical variables
vars = c( "source","Sex"  ,"Past.medical.history")
p_values = NULL
for (v in c(1:length(vars))){
  t = table(df[,vars[v]], df$Group)
  p_val = signif(fisher.test(t, alternative = "greater",simulate.p.value = T,B = 1000000)$p.value, digits = 4)
  p_values = c(p_values, p_val)
}
names(p_values) = vars

p_values
# not significant










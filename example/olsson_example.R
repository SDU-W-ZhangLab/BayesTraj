
library(BayesTraj)

# data <- read.csv('data_Olsson_wishbone.csv',row.names = 1)
data <- read.csv('data/multifurcation/data_Olsson.tsv.gz',sep = '\t',row.names = 1)

gene_names <- rownames(data)
data_counts <- t(as.matrix(data)) 

ery_genes <- c('Gata1','Klf1','Epor','Itga2b','Pf4','Tal1')
gran_genes <- c('S100a8','S100a9','Gfi1','Camp','Cebpe')
mono_genes <- c('Ly86','Irf8','Ccr2','Cx3cr1')

ery_genes %in% gene_names
gran_genes %in% gene_names
mono_genes %in% gene_names


gene_expression_matrix=data_counts[,c(ery_genes,gran_genes,mono_genes)]
gene_expression_matrix <- t(gene_expression_matrix)


cell_labels <- read.csv('cell_label.tsv.gz',sep = '\t',header = F)$V1
desired_order <- c('HSCP-1','HSCP-2','Multi-Lin','Meg','Eryth','Gran','Myelocyte','MDP','Mono')
annotation <- data.frame(Label = factor(cell_labels, levels = desired_order))
rownames(annotation) <- colnames(t(data_counts))


sorted_order <- order(annotation$Label)
gene_expression_matrix <- gene_expression_matrix[, sorted_order]
annotation <- annotation[sorted_order, , drop = FALSE]


unique_labels <- levels(annotation$Label)
# colors <- rainbow(length(unique_labels)) 
colors <- c('#FEE551','#E4BE00','#E84037','#DB1135','#9B278A','#6081C0','#55B1E4','#00973D','#71BD58')
names(colors) <- unique_labels
ann_colors <- list(Label = colors)

library(RColorBrewer)
bluecolors <- brewer.pal(9,"Blues")
redcolors <- brewer.pal(9,"Reds")
color_palette <- c(colorRampPalette(colors = c(bluecolors[9], bluecolors[7]))(20),
                   colorRampPalette(colors = c(bluecolors[7], bluecolors[5]))(20),
                   colorRampPalette(colors = c(bluecolors[5], "white"))(10),
                   colorRampPalette(colors = c("white", redcolors[5]))(10),
                   colorRampPalette(colors = c(redcolors[5], redcolors[6]))(20),
                   colorRampPalette(colors = c(redcolors[6], redcolors[7]))(20))

library(pheatmap)
pheatmap(gene_expression_matrix,
         annotation_col = annotation,      
         annotation_colors = ann_colors,   
         color = color_palette,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F
)


b1_genes <- ery_genes
b2_genes <- gran_genes
b3_genes <- mono_genes

b1_genes %in% gene_names
b2_genes %in% gene_names
b3_genes %in% gene_names

sub_data_counts <- normalized_expression_data(data_counts)
sub_data_counts <- sub_data_counts[,c(b1_genes,b2_genes,b3_genes)]
Y1 <- sub_data_counts[,c(b1_genes),drop = F]
Y2 <- sub_data_counts[,c(b2_genes),drop = F]
Y3 <- sub_data_counts[,c(b3_genes),drop = F]

ery_df <- sub_data_counts[cell_labels %in% c('Eryth','Meg'),]
mono_df <- sub_data_counts[cell_labels %in% c('Mono','MDP'),]
gran_df <- sub_data_counts[cell_labels %in% c('Gran','Myelocyte'),]
par(mfrow = c(3,1))
boxplot(ery_df)
boxplot(gran_df)
boxplot(mono_df)


quantile_ery <- apply(ery_df,2,function(x){
  quantile(x,probs = c(0.25,0.5,0.75))
})
quantile_gran <- apply(gran_df,2,function(x){
  quantile(x,probs = c(0.25,0.5,0.75))
})
quantile_mono <- apply(mono_df,2,function(x){
  quantile(x,probs = c(0.25,0.5,0.75))
})

quantile_all <- rbind(quantile_ery[2,],quantile_gran[2,],quantile_mono[2,])
quantile_all <- quantile_all[,unique(c(ery_genes,gran_genes,mono_genes))]

transient_mu <- apply(quantile_all, 2, function(x){
  sort(x, decreasing = TRUE)[2]
})

quantile_all <- rbind(quantile_ery[3,],quantile_gran[3,],quantile_mono[3,])
quantile_all <- quantile_all[,unique(c(ery_genes,gran_genes,mono_genes))]

switch_mu <- apply(quantile_all, 2, function(x){
  sort(x, decreasing = TRUE)[1]
})


b1_mu0_s_mu = switch_mu[ery_genes]/2
b1_mu0_s_sd = rep(0.05,length(ery_genes))
b2_mu0_s_mu = switch_mu[c(gran_genes)]/2
b2_mu0_s_sd = rep(0.05,length(c(gran_genes)))
b3_mu0_s_mu = switch_mu[c(mono_genes)]/2
b3_mu0_s_sd = rep(0.05,length(c(mono_genes)))

num_genes <- length(unique(c(ery_genes,gran_genes,mono_genes)))

b1_mu0_t_mu <- transient_mu[!(names(transient_mu) %in% ery_genes)]/2
b1_mu0_t_sd <- rep(0.02,num_genes-length(ery_genes))

b2_mu0_t_mu <- transient_mu[!(names(transient_mu) %in% c(gran_genes))]/2
b2_mu0_t_sd <- rep(0.02,num_genes-length(c(gran_genes)))

b3_mu0_t_mu <- transient_mu[!(names(transient_mu) %in% c(mono_genes))]/2
b3_mu0_t_sd <- rep(0.02,num_genes-length(c(mono_genes)))


bt <- BayesTraj(iter = 100,n_branch = 3,
                Y1 = Y1,
                Y2 = Y2,
                Y3 = Y3,
                alpha = 0.1,eplison = 0.04,
                lower_k_s = 10,lower_k_t = 4,lower_t0 = 0.3,upper_b1_t0 = 0.5,lower_b2_t0 = 0.5,
                lower_mu0_s = 0.2,upper_mu0_s = 0.5,lower_mu0_t = 0,upper_mu0_t=0.25,
                
                b1_mu0_s_mu = b1_mu0_s_mu,b1_mu0_s_sd = b1_mu0_s_sd,
                b2_mu0_s_mu = b2_mu0_s_mu,b2_mu0_s_sd = b2_mu0_s_sd ,
                b3_mu0_s_mu = b3_mu0_s_mu,b3_mu0_s_sd = b3_mu0_s_sd,
                b1_mu0_t_mu = b1_mu0_t_mu,b1_mu0_t_sd = b1_mu0_t_sd,
                b2_mu0_t_mu = b2_mu0_t_mu,b2_mu0_t_sd = b2_mu0_t_sd ,
                b3_mu0_t_mu = b3_mu0_t_mu,b3_mu0_t_sd = b3_mu0_t_sd,
                
                chains = 1)

bt <- infer_pseudotime(bt)
bt <- infer_mu0(bt)
bt <- infer_k(bt)
bt <- infer_t0(bt)
bt <- branch_probability(bt)
bt <- calculate_branches(bt,limit_prob = 0.9)

unname(bt$t)
bt$para_mu0
bt$para_k
bt$para_t0

label <- read.csv('cell_label.tsv.gz',sep = '\t',header = F)
rownames(label) <- rownames(data_counts)

all(rownames(Y1) == rownames(label))
label <- label[match(rownames(Y1),rownames(label)),]


plot_branches(bt,branch = 1,ncol = 4,label = label,label_color = colors)
plot_branches(bt,branch = 2,ncol = 4,label = label,label_color = colors)
plot_branches(bt,branch = 3,ncol = 4,label = label,label_color = colors)

t1 <- bt$all_branch_data[[1]]$t
expr1 <- data[,bt$all_branch_data[[1]]$index]
# plot(t1,expr1['Gfi1',])

filtered_genes1 <- filter_de_genes(t1,as.matrix(expr1),limit_cor = 0.3,limit_bf = 1e5,limit_min_t0 = 0.15,limit_max_t0 = 1,
                                   limit_min_k = 0,limit_min_mu0 = 0)

genes1 <- filtered_genes1[['selected_genes']]

para1 <- filtered_genes1$bf
ordered_para1 <- para1[order(para1[,1],decreasing = T),]
head(ordered_para1)

t2 <- bt$all_branch_data[[2]]$t
expr2 <- data[,bt$all_branch_data[[2]]$index]
# plot(t2,expr2['Gfi1',])
filtered_genes2 <- filter_de_genes(t2,as.matrix(expr2),limit_cor = 0.45,limit_bf = 1e17,limit_min_t0 = 0.15,limit_max_t0 = 1,
                                   limit_min_k = 0,limit_min_mu0 = 0,head = 100)

genes2 <- filtered_genes2[['selected_genes']]

para2 <- filtered_genes2$bf
ordered_para2 <- para2[order(para2[,1],decreasing = T),]
head(ordered_para2)

t3 <- bt$all_branch_data[[3]]$t
expr3 <- data[,bt$all_branch_data[[3]]$index]
# plot(t3,expr3['Gfi1',])
filtered_genes3 <- filter_de_genes(t3,as.matrix(expr3),limit_cor = 0.4,limit_bf = 1e14,limit_min_t0 = 0.15,limit_max_t0 = 1,
                                   limit_min_k = 0,limit_min_mu0 = 0)

genes3 <- filtered_genes3[['selected_genes']]

para3 <- filtered_genes3$bf
ordered_para3 <- para3[order(para3[,1],decreasing = T),]
head(ordered_para3)


par(mfrow = c(2,5))
for (gene in rownames(ordered_para1)[1:min(20,nrow(ordered_para1))]) {
  plot(t1, expr1[gene,], main = gene, xlab = "pseudotime", ylab = "expression", pch = 16,col = 'grey20')
  curve(ordered_para1[gene,2] / (1 + exp(-ordered_para1[gene,3] * (x - ordered_para1[gene,4]))), 
        add = TRUE, col = "red", lwd = 2)
}
for (gene in rownames(ordered_para2)[1:min(20,nrow(ordered_para2))]) {
  plot(t2, expr2[gene,], main = gene, xlab = "pseudotime", ylab = "expression", pch = 16,col = 'grey20')
  curve(ordered_para2[gene,2] / (1 + exp(-ordered_para2[gene,3] * (x - ordered_para2[gene,4]))), 
        add = TRUE, col = "red", lwd = 2)
}
for (gene in rownames(ordered_para3)[1:min(20,nrow(ordered_para3))]) {
  plot(t3, expr3[gene,], main = gene, xlab = "pseudotime", ylab = "expression", pch = 16,col = 'grey20')
  curve(ordered_para3[gene,2] / (1 + exp(-ordered_para3[gene,3] * (x - ordered_para3[gene,4]))), 
        add = TRUE, col = "red", lwd = 2)
}


expr1 <- expr1[,order(t1)]
expr2 <- expr2[,order(t2)]
expr3 <- expr3[,order(t3)]
t1 <- t1[order(t1)]
t2 <- t2[order(t2)]
t3 <- t3[order(t3)]

genes1 <- setdiff(genes1, union(genes2,genes3))
genes2 <- setdiff(genes2, union(genes1,genes3))
genes3 <- setdiff(genes3, union(genes1,genes2))

genes_multi1 <- union(rownames(filtered_genes1$cor_matrix[filtered_genes1$cor_matrix[,2] > 0.28,]),
                      intersect(rownames(filtered_genes2$cor_matrix[filtered_genes2$cor_matrix[,2] > 0.2,]),
                                rownames(filtered_genes3$cor_matrix[filtered_genes3$cor_matrix[,2] > 0.25,])))
genes_multi2 <- union(intersect(rownames(filtered_genes1$cor_matrix[filtered_genes1$cor_matrix[,2] > 0.2,]),
                                rownames(filtered_genes2$cor_matrix[filtered_genes2$cor_matrix[,2] > 0.2,])),
                      intersect(rownames(filtered_genes1$cor_matrix[filtered_genes1$cor_matrix[,2] > 0.2,]),
                                rownames(filtered_genes3$cor_matrix[filtered_genes3$cor_matrix[,2] > 0.2,])))

genes1 <- setdiff(genes1,genes_multi2)
genes2 <- setdiff(genes2,genes_multi1)
genes3 <- setdiff(genes3,genes_multi1)

genes <- c(genes1,genes2,genes3)
smooth_expr1 <- process_and_smooth_data(expr1,genes,smooth_k = 8)
smooth_expr2 <- process_and_smooth_data(expr2,genes,smooth_k = 8)
smooth_expr3 <- process_and_smooth_data(expr3,genes,smooth_k = 8)

input_list <- list(smooth_expr1,smooth_expr2,smooth_expr3)
generate_heatmap(input_list)
















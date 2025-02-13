
library(BayesTraj)

setwd('D:/study/2024/cellfate/paper/monocle2/lung/')
data <- read.csv('lung_exprs_data.csv',header = T,row.names = 1)

gene_names <- rownames(data)


data_counts <- t(log1p(data))

AT1_genes <- c('Clic5','Pdpn','Cav1','Ager','Hopx','Vegfa')
AT2_genes <- c('Sftpa1','Lyz2','Sftpb','Sftpc','Slc34a2')

AT1_genes %in% gene_names
AT2_genes %in% gene_names

gene_expression_matrix=data_counts[,c(AT1_genes,AT2_genes)]
gene_expression_matrix <- t(gene_expression_matrix)


cell_labels <- read.csv('cell_labels.csv',header = F)$V1
desired_order <- c('cell_type_1','cell_type_2','cell_type_3')
# 将 cell_labels 设置为有序因子，使得类别按照 desired_order 排序
annotation <- data.frame(Label = factor(cell_labels, levels = desired_order))
rownames(annotation) <- colnames(t(data_counts))

sorted_order <- order(annotation$Label)
gene_expression_matrix <- gene_expression_matrix[, sorted_order]
annotation <- annotation[sorted_order, , drop = FALSE]


unique_labels <- levels(annotation$Label)
# colors <- rainbow(length(unique_labels))  # 生成 9 种不同的颜色
colors <- c('blue','#FFFF00','green')
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
         annotation_col = annotation,       # 列注释
         annotation_colors = ann_colors,    # 标签颜色
         color = color_palette,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F
)


b1_genes <- AT1_genes
b2_genes <- AT2_genes

b1_genes %in% gene_names
b2_genes %in% gene_names

sub_data_counts <- normalized_expression_data(data_counts)
sub_data_counts <- sub_data_counts[,c(b1_genes,b2_genes)]
Y1 <- sub_data_counts[,c(b1_genes),drop = F]
Y2 <- sub_data_counts[,c(b2_genes),drop = F]

AT1_df <- sub_data_counts[cell_labels %in% c('cell_type_2'),]
AT2_df <- sub_data_counts[cell_labels %in% c('cell_type_3'),]
par(mfrow = c(2,1))
boxplot(AT1_df)
boxplot(AT2_df)


quantile_AT1 <- apply(AT1_df,2,function(x){
  quantile(x,probs = c(0.25,0.5,0.75))
})
quantile_AT2 <- apply(AT2_df,2,function(x){
  quantile(x,probs = c(0.25,0.5,0.75))
})

quantile_all <- rbind(quantile_AT1[2,],quantile_AT2[2,])
quantile_all <- quantile_all[,unique(c(AT1_genes,AT2_genes))]

transient_mu <- apply(quantile_all, 2, function(x){
  sort(x, decreasing = TRUE)[2]
})

quantile_all <- rbind(quantile_AT1[3,],quantile_AT2[3,])
quantile_all <- quantile_all[,unique(c(AT1_genes,AT2_genes))]

switch_mu <- apply(quantile_all, 2, function(x){
  sort(x, decreasing = TRUE)[1]
})


b1_mu0_s_mu = switch_mu[AT1_genes]/2
b1_mu0_s_sd = rep(0.05,length(AT1_genes))
b2_mu0_s_mu = switch_mu[c(AT2_genes)]/2
b2_mu0_s_sd = rep(0.05,length(c(AT2_genes)))

num_genes <- length(unique(c(AT1_genes,AT2_genes)))

b1_mu0_t_mu <- transient_mu[!(names(transient_mu) %in% AT1_genes)]/4
b1_mu0_t_sd <- rep(0.02,num_genes-length(AT1_genes))

b2_mu0_t_mu <- transient_mu[!(names(transient_mu) %in% c(AT2_genes))]/4
b2_mu0_t_sd <- rep(0.02,num_genes-length(c(AT2_genes)))

bt <- BayesTraj(iter = 100,n_branch = 2,
                Y1 = Y1,
                Y2 = Y2,
                alpha = 0.1,eplison = 0.04,
                lower_k_s = 10,lower_k_t = 4,lower_t0 = 0.3,upper_b1_t0 = 0.8,lower_b2_t0 = 0.3,
                lower_mu0_s = 0.2,upper_mu0_s = 0.5,lower_mu0_t = 0.05,upper_mu0_t=0.4,
                
                b1_mu0_s_mu = b1_mu0_s_mu,b1_mu0_s_sd = b1_mu0_s_sd,
                b2_mu0_s_mu = b2_mu0_s_mu,b2_mu0_s_sd = b2_mu0_s_sd ,
                b1_mu0_t_mu = b1_mu0_t_mu,b1_mu0_t_sd = b1_mu0_t_sd,
                b2_mu0_t_mu = b2_mu0_t_mu,b2_mu0_t_sd = b2_mu0_t_sd ,
                
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

label <- read.csv('cell_labels.csv',header = F)
rownames(label) <- rownames(data_counts)

all(rownames(Y1) == rownames(label))
label <- label[match(rownames(Y1),rownames(label)),]

plot_branches(bt,branch = 1,ncol = 4,label = label)
plot_branches(bt,branch = 2,ncol = 4,label = label)

t1 <- bt$all_branch_data[[1]]$t
expr1 <- data[,bt$all_branch_data[[1]]$index]
plot(t1,expr1['Clic5',])
filtered_genes1 <- filter_de_genes(t1,as.matrix(expr1),limit_cor = 0.4,limit_bf = 1,limit_min_t0 = 0.1,limit_max_t0 = 1,
                                  limit_min_k = 0,limit_min_mu0 = 0)

genes1 <- filtered_genes1[['selected_genes']]

para1 <- filtered_genes1$bf
ordered_para1 <- para1[order(para1[,1],decreasing = T),]
head(ordered_para1)

t2 <- bt$all_branch_data[[2]]$t
expr2 <- data[,bt$all_branch_data[[2]]$index]
plot(t2,expr2['Clic5',])
filtered_genes2 <- filter_de_genes(t2,as.matrix(expr2),limit_cor = 0.4,limit_bf = 1,limit_min_t0 = 0.1,limit_max_t0 = 1,
                                   limit_min_k = 0,limit_min_mu0 = 0)

genes2 <- filtered_genes2[['selected_genes']]

para2 <- filtered_genes2$bf
ordered_para2 <- para2[order(para2[,1],decreasing = T),]
head(ordered_para2)


par(mfrow = c(2,5))
for (gene in rownames(ordered_para1)[1:20]) {
  plot(t1, expr1[gene,], main = gene, xlab = "pseudotime", ylab = "expression", pch = 16,col = 'grey20')
  curve(ordered_para1[gene,2] / (1 + exp(-ordered_para1[gene,3] * (x - ordered_para1[gene,4]))), 
        add = TRUE, col = "red", lwd = 2)
}
for (gene in rownames(ordered_para2)[1:20]) {
  plot(t2, expr2[gene,], main = gene, xlab = "pseudotime", ylab = "expression", pch = 16,col = 'grey20')
  curve(ordered_para2[gene,2] / (1 + exp(-ordered_para2[gene,3] * (x - ordered_para2[gene,4]))), 
        add = TRUE, col = "red", lwd = 2)
}


expr1 <- expr1[,order(t1)]
expr2 <- expr2[,order(t2)]
t1 <- t1[order(t1)]
t2 <- t2[order(t2)]

genes <- c(genes1,genes2)
smooth_expr1 <- process_and_smooth_data(expr1,genes,smooth_k = 8)
smooth_expr2 <- process_and_smooth_data(expr2,genes,smooth_k = 8)

input_list <- list(smooth_expr1,smooth_expr2)
generate_heatmap(input_list)















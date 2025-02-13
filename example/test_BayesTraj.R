
library(BayesTraj)

setwd('D:/study/2024/cellfate/data/human_1')
data <- read.csv('ExpressionData.csv',header = T,row.names = 1)

gene_names <- rownames(data)

data_counts <- t(data)

AT1_genes <- c('EOMES','CER1', 'GATA4', 'DKK4','MYCT1','PRDM1')
AT2_genes <- c('CDX1','MSX2','T')

AT1_genes %in% gene_names
AT2_genes %in% gene_names

gene_expression_matrix=data_counts[,c(AT1_genes,AT2_genes)]
gene_expression_matrix <- t(gene_expression_matrix)

library(pheatmap)
library(RColorBrewer)
bluecolors <- brewer.pal(9,"Blues")
redcolors <- brewer.pal(9,"Reds")
color_palette <- c(colorRampPalette(colors = c(bluecolors[9], bluecolors[7]))(20),
                   colorRampPalette(colors = c(bluecolors[7], bluecolors[5]))(20),
                   colorRampPalette(colors = c(bluecolors[5], "white"))(10),
                   colorRampPalette(colors = c("white", redcolors[5]))(10),
                   colorRampPalette(colors = c(redcolors[5], redcolors[6]))(20),
                   colorRampPalette(colors = c(redcolors[6], redcolors[7]))(20))
pheatmap(gene_expression_matrix,
         scale = "row",
         color = color_palette,
         cluster_rows = TRUE,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F
)

b1_genes <- AT1_genes
b2_genes <- AT2_genes

b1_genes %in% gene_names
b2_genes %in% gene_names

sub_data_counts <- normalized_expression_data(t(gene_expression_matrix))
Y1 <- sub_data_counts[,c(b1_genes),drop = F]
Y2 <- sub_data_counts[,c(b2_genes),drop = F]

b1_mu0_s_mu = rep(0.5,length(AT1_genes))
b1_mu0_s_sd = rep(0.2,length(AT1_genes))
b2_mu0_s_mu = rep(0.5,length(c(AT2_genes)))
b2_mu0_s_sd = rep(0.2,length(c(AT2_genes)))

num_genes <- length(unique(c(AT1_genes,AT2_genes)))

b1_mu0_t_mu <- rep(0.25,num_genes-length(AT1_genes))
b1_mu0_t_sd <- rep(0.2,num_genes-length(AT1_genes))

b2_mu0_t_mu <- rep(0.25,num_genes-length(c(AT2_genes)))
b2_mu0_t_sd <- rep(0.2,num_genes-length(c(AT2_genes)))

b1_k_s_mu <- rep(0,length(AT1_genes))
b1_k_s_sd <- rep(50,length(AT1_genes))

bt <- BayesTraj(iter = 100,n_branch = 1,
                            Y1 = Y1,
                            Y2 = Y2,
                            alpha = 0.1,eplison = 0.04,
                            lower_k_s = 5,lower_k_t = 1,lower_t0 = 0.1,upper_b1_t0 = 0.8,lower_b2_t0 = 0.1,
                            b1_mu0_s_mu = b1_mu0_s_mu,b1_mu0_s_sd = b1_mu0_s_sd,
                            # b1_k_s_mu = b1_k_s_mu,b1_k_s_sd = b1_k_s_sd,
                            b2_mu0_s_mu = b2_mu0_s_mu,b2_mu0_s_sd = b2_mu0_s_sd ,
                            b1_mu0_t_mu = b1_mu0_t_mu,b1_mu0_t_sd = b1_mu0_t_sd,b2_mu0_t_mu = b2_mu0_t_mu,b2_mu0_t_sd = b2_mu0_t_sd ,
                            lower_mu0_s = 0.2,upper_mu0_s = 0.5,lower_mu0_t = 0.05,upper_mu0_t=0.4,
                            chains = 1)

bt <- infer_pseudotime(bt)
bt <- infer_mu0(bt)
bt <- infer_k(bt)
bt <- infer_t0(bt)
bt <- branch_probability(bt)
bt <- calculate_branches(bt)

label <- rep('esc',nrow(Y1))
color <- c('lightblue')
names(color) <- 'esc'
plot_branches(bt,branch = 1,ncol = 4,label = label,label_color = color)

t <- bt$all_branch_data[[1]]$t
expr <- data[,bt$all_branch_data[[1]]$index]
# debug(filter_de_genes)
filtered_genes <- filter_de_genes(t,as.matrix(expr),limit_cor = 0.4,limit_bf = 1e40,limit_min_t0 = 0.1,limit_max_t0 = 1,
                limit_min_k = 10,limit_min_mu0 = 2.5)

genes <- filtered_genes[['selected_genes']]

para <- filtered_genes$bf
ordered_para <- para[order(para[,1],decreasing = T),]
head(ordered_para)

par(mfrow = c(2,5))
for (gene in rownames(ordered_para)[1:20]) {
  plot(t, expr[gene,], main = gene, xlab = "pseudotime", ylab = "expression", pch = 16,col = 'grey20')
  curve(ordered_para[gene,2] / (1 + exp(-ordered_para[gene,3] * (x - ordered_para[gene,4]))), 
        add = TRUE, col = "red", lwd = 2)
}

# debug(process_and_smooth_data)
expr1 <- process_and_smooth_data(expr,genes,smooth_k = 8)

input_list <- list(expr1)
# debug(generate_heatmap)
generate_heatmap(input_list)















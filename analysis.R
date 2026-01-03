# load libraries
library(survival)
library(survminer)
library(tidyverse)
library(gtsummary)
library(fmsb)
library(caret)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggsankey)
library(tidyestimate)
library(decoderr)
library(survcomp)

# load summarized data
tmp.surv <- read.csv("path_to_clinical_data.csv", na.strings = c("", "NA"))

# binomial test to compare probability of basal-like with exome-capture using whole-transcriptome as reference
binom.test(length(which(tmp.surv$PurIST_tempus=="Basal-like")), 
           nrow(tmp.surv), 
           p = length(which(tmp.surv$PurIST_rnaseq=="Basal-like")) / nrow(tmp.surv),
           alternative = "greater")

## relevel subtypes
tmp.surv$PurIST_rnaseq <- factor(tmp.surv$PurIST_rnaseq, levels = c("Classical", "Basal-like"))
tmp.surv$PurIST_tempus <- factor(tmp.surv$PurIST_tempus, levels = c("Classical", "Basal-like"))

# kappa test (figure 1a)
confusionMatrix(tmp.surv$PurIST_tempus, tmp.surv$PurIST_rnaseq, mode = "prec_recall", positive = "Basal-like")

## to get the confidence interval of kappa
tmp.surv %>%
  {table(RNAseq = .$PurIST_rnaseq, Tempus = .$PurIST_tempus)} %>%
  Kappa.test()

# survival by purist rna-seq (figure 3a)
os.fit.rna <- survfit(Surv(Time.to.os, OS.code) ~ PurIST_rnaseq, conf.type = "plain", data = tmp.surv)
ggsurvplot(os.fit.rna, size = 1.2, pval = T, censor = T, censor.shape = 108, censor.size = 3,
           palette = c("blue", "orange"),
           xlim = c(1.5, 72), break.x.by = 12,
           risk.table = T, risk.table.y.text = F, tables.col = "strata", tables.height = .2,
           tables.theme = list(theme_cleantable(),
                               theme(plot.title=element_text(size=12))),
           legend.labs = c("Classical", "Basal-like"), legend.title = "RNA-seq subtype", legend = c(.9,.9),
           surv.median.line = "hv",
           xlab = "Time (months)")
print(os.fit.rna) # get median with confidence intervals
rna.cox <- coxph(Surv(Time.to.os, OS.code) ~ PurIST_rnaseq, data = tmp.surv) # coxph model (figure s3a)
summary(rna.cox)

# survival by tempus rna-seq (figure 3b)
os.fit.tmp <- survfit(Surv(Time.to.os, OS.code) ~ PurIST_tempus, conf.type = "plain", data = tmp.surv)
ggsurvplot(os.fit.tmp, size = 1.2, pval = T, censor = T, censor.shape = 108, censor.size = 3,
           palette = c("blue", "orange"),
           xlim = c(1.5, 72), break.x.by = 12,
           risk.table = T, risk.table.y.text = F, tables.col = "strata", tables.height = .2,
           tables.theme = list(theme_cleantable(),
                               theme(plot.title=element_text(size=12))),
           legend.labs = c("Classical", "Basal-like"), legend.title = "Tempus subtype", legend = c(.9,.9),
           surv.median.line = "hv",
           xlab = "Time (months)")
print(os.fit.tmp) # get median with confidence intervals
tmp.cox <- coxph(Surv(Time.to.os, OS.code) ~ PurIST_tempus, data = tmp.surv) # coxph model (figure s3a)
summary(tmp.cox)

## compare log-likelihood between methods
## model with higher log-likelihood (less negative value) provides a better fit
rna.loglik <- logLik(rna.cox)
tmp.loglik <- logLik(tmp.cox)
rna.loglik
tmp.loglik
anova(rna.cox, tmp.cox, test = "LRT")

# forest plot of hazard ratios for os (figure s3a)
rna.cox.sum <- data.frame(
  method = "RNA-seq",
  outcome = "OS",
  variable = "Basal-like",
  hazard_ratio = exp(coef(rna.cox)),
  ci_lower = exp(confint(rna.cox))[,1],
  ci_upper = exp(confint(rna.cox))[,2],
  p = summary(rna.cox)$coefficients[,"Pr(>|z|)"],
  n_in_model = nobs(rna.cox)
)

tmp.cox.sum <- data.frame(
  method = "Tempus",
  outcome = "OS",
  variable = "Basal-like",
  hazard_ratio = exp(coef(tmp.cox)),
  ci_lower = exp(confint(tmp.cox))[,1],
  ci_upper = exp(confint(tmp.cox))[,2],
  p = summary(tmp.cox)$coefficients[,"Pr(>|z|)"],
  n_in_model = nobs(tmp.cox)
)

os.cox.sum <- rbind(rna.cox.sum, tmp.cox.sum)
os.cox.sum$method <- factor(os.cox.sum$method, levels = c("Tempus", "RNA-seq"))

print(os.cox.sum[c(1,2,4:8)])

ggplot(os.cox.sum, aes(y = method, x = hazard_ratio, color = method)) +
  geom_point(shape = 18, size = 15) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                 height = .2, linewidth = .5) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .5, color = "red") +
  geom_vline(xintercept = c(.1, 10), linetype = "dotted",linewidth = .5, color = "gray70") +
  scale_x_log10(limits = c(.1, max(os.cox.sum$ci_upper) * 1.2),
                labels = function(x) ifelse(x == .1, "0.1", as.character(x))) +
  labs(x = "HR with 95% CI",
       y = NULL) +
  scale_color_manual(values = c("RNA-seq"="black", "Tempus"="orange")) +
  theme(panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        legend.position = "none",
        legend.key = element_blank(),
        aspect.ratio = .5)

# progression free survival (figures 3c and 3d)
pfs.fit.rna <- survfit(Surv(Time.to.pfs, PFS.code) ~ PurIST_rnaseq, conf.type = "plain", data = tmp.surv)
ggsurvplot(pfs.fit.rna, size = 1.2, pval = T, censor = T, censor.shape = 108, censor.size = 3,
           palette = c("blue", "orange"),
           xlim = c(1.5, 72), break.x.by = 12,
           risk.table = T, risk.table.y.text = F, tables.col = "strata", tables.height = .2,
           tables.theme = list(theme_cleantable(),
                               theme(plot.title=element_text(size=12))),
           legend.labs = c("Classical", "Basal-like"), legend.title = "RNA-seq subtype", legend = c(.9,.9),
           surv.median.line = "hv",
           xlab = "Time (months)")
print(pfs.fit.rna)
rna.cox2 <- coxph(Surv(Time.to.pfs, PFS.code) ~ PurIST_rnaseq, data = tmp.surv)
summary(rna.cox2)

pfs.fit.tmp <- survfit(Surv(Time.to.pfs, PFS.code) ~ PurIST_tempus, conf.type = "plain", data = tmp.surv)
ggsurvplot(pfs.fit.tmp, size = 1.2, pval = T, censor = T, censor.shape = 108, censor.size = 3,
           palette = c("blue", "orange"),
           xlim = c(1.5, 72), break.x.by = 12,
           risk.table = T, risk.table.y.text = F, tables.col = "strata", tables.height = .2,
           tables.theme = list(theme_cleantable(),
                               theme(plot.title=element_text(size=12))),
           legend.labs = c("Classical", "Basal-like"), legend.title = "Tempus subtype", legend = c(.9,.9),
           surv.median.line = "hv",
           xlab = "Time (months)")
print(pfs.fit.tmp)
tmp.cox2 <- coxph(Surv(Time.to.pfs, PFS.code) ~ PurIST_tempus, data = tmp.surv)
summary(tmp.cox2)

# forest plot of hazard ratios for pfs (figure s3b)
rna.cox.sum2 <- data.frame(
  method = "RNA-seq",
  outcome = "PFS",
  variable = "Basal-like",
  hazard_ratio = exp(coef(rna.cox2)),
  ci_lower = exp(confint(rna.cox2))[,1],
  ci_upper = exp(confint(rna.cox2))[,2],
  p = summary(rna.cox2)$coefficients[,"Pr(>|z|)"],
  n_in_model = nobs(rna.cox2)
)

tmp.cox.sum2 <- data.frame(
  method = "Tempus",
  outcome = "PFS",
  variable = "Basal-like",
  hazard_ratio = exp(coef(tmp.cox2)),
  ci_lower = exp(confint(tmp.cox2))[,1],
  ci_upper = exp(confint(tmp.cox2))[,2],
  p = summary(tmp.cox2)$coefficients[,"Pr(>|z|)"],
  n_in_model = nobs(tmp.cox2)
)

pfs.cox.sum <- rbind(rna.cox.sum2, tmp.cox.sum2)
pfs.cox.sum$method <- factor(pfs.cox.sum$method, levels = c("Tempus", "RNA-seq"))

print(pfs.cox.sum[c(1,2,4:8)])

ggplot(pfs.cox.sum, aes(y = method, x = hazard_ratio, color = method)) +
  geom_point(shape = 18, size = 15) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                 height = .2, linewidth = .5) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .5, color = "red") +
  geom_vline(xintercept = c(.1, 10), linetype = "dotted",linewidth = .5, color = "gray70") +
  scale_x_log10(limits = c(.1, max(os.cox.sum$ci_upper) * 1.2),
                labels = function(x) ifelse(x == .1, "0.1", as.character(x))) +
  labs(x = "HR with 95% CI",
       y = NULL) +
  scale_color_manual(values = c("RNA-seq"="black", "Tempus"="orange")) +
  theme(panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        legend.position = "none",
        legend.key = element_blank(),
        aspect.ratio = .5)

# make heatmaps for relative gene expression of purist tsps by method (figure 2b)
rna <- read.csv("path_to_rna_tpm.csv", row.names = 1) # load tpm matrix of rna-seq
tmp <- read.csv("path_to_tempus_tpm.csv", row.names = 1) # load tpm matrix of tempus results
purist_genes <- read.csv("path_to_purist_genes.csv") # load purist genes

## rank normalize genes (entire genes in expression matrix)
rna.ranked <- as.data.frame(apply(rna, 2, percent_rank))
rownames(rna.ranked) <- rownames(rna)

tmp.ranked <- as.data.frame(apply(tmp, 2, percent_rank))
rownames(tmp.ranked) <- rownames(tmp)

## extract purist gene pairs from the rank normalized expression matrix
rna.tsp.ranked <- rna.ranked[rownames(rna.ranked) %in% purist_genes$Gene,]
tmp.tsp.ranked <- tmp.ranked[rownames(tmp.ranked) %in% purist_genes$Gene,]

## run purist on rna-seq and tempus to make annotation tables
## make sure to have purist classifier functions loaded per https://github.com/naimurashid/PurIST
rna.pred <- apply_classifier(rna, classifier)
rna.pred$ptID <- rownames(rna.pred) %>% str_extract("\\d{4}")
rna.pred <- rna.pred %>%
  rownames_to_column("sampID")
tmp.pred <- apply_classifier(tmp, classifier)
tmp.pred$ptID <- rownames(tmp.pred) %>% str_extract("\\d{4}")
tmp.pred <- tmp.pred %>% 
  rownames_to_column("sampID")

## merge rna-seq predictions with tempus predictions
merged.pred <- merge(tmp.pred, rna.pred, by = "ptID") # y will be rna-seq and x will be tempus

## need to remove extraneous characters (periods and parentheses) from names to match
tmp.surv$Sample <- gsub("\\.|\\(|\\)", "", tmp.surv$Sample) 
merged.pred$sampID.y <- gsub("\\.|\\(|\\)", "", merged.pred$sampID.y)

## order predictions by purist probability from rna-seq
rna.pred <- rna.pred[order(rna.pred$PurIST_prob),]
merged.pred <- merged.pred[order(merged.pred$PurIST_prob.y),] # remember y is rna and x is tempus

## fix column names from expression data and sample names from rna prediction to remove special characters
colnames(rna.tsp.ranked) <- gsub("\\(|\\)", "", colnames(rna.tsp.ranked))
rna.pred$sampID <- gsub("\\.", "", rna.pred$sampID)

## order columns of tsp expression matrices by sample ID in prediction
rna.tsp.ranked <- rna.tsp.ranked[,rna.pred$sampID]
tmp.tsp.ranked <- tmp.tsp.ranked[,merged.pred$sampID.x] # lose 1 sample (3029 since it was excluded from freeze)

## data prep for heatmap
## define colors to be used
cont_colors <- viridis::viridis(10)
sample_color = list(PurIST = c("Classical"="blue", "Basal-like"="orange")) # column annotation
sample_color2 = list(Tempus_PurIST = c("Classical"="blue", "Basal-like"="orange"),
                     RNAseq_PurIST = c("Classical"="blue", "Basal-like"="orange")) # for the tempus heatmap
gene_color = structure(c("blue", "orange"), names = c("classical", "basal")) # row heatmap

## combine data frames
rna_purist <- merge(purist_genes, rna.tsp.ranked, by.x = 1, by.y = 0)
rna_purist <- rna_purist[order(rna_purist$TSP_wt, rna_purist$TSP, rna_purist$Subtype),]
rownames(rna_purist) <- rna_purist$Gene

tmp_purist <- merge(purist_genes, tmp.tsp.ranked, by.x = 1, by.y = 0)
tmp_purist <- tmp_purist[order(tmp_purist$TSP_wt, tmp_purist$TSP, tmp_purist$Subtype),]
rownames(tmp_purist) <- tmp_purist$Gene

## prepare heatmaps
## set universal padding between column annotation and heatmap body
ht_opt$COLUMN_ANNO_PADDING = unit(2, "mm")

## rna-seq
rna_column_ha <- HeatmapAnnotation(df = rna.pred["PurIST"],
                                   col = sample_color,
                                   simple_anno_size = unit(4, "mm"),
                                   show_annotation_name = F)
rna_row_hm <- Heatmap(as.matrix(rna_purist$Subtype),
                      row_title = "PurIST gene pairs", row_title_rot = 90,
                      show_row_names = F,
                      cluster_rows = F,
                      row_split = rna_purist$TSP_wt,
                      col = gene_color,
                      show_heatmap_legend = F,
                      use_raster = T,
                      width = unit(4, "mm"),
                      height = unit(7, "cm"))
tsp_row_hm <- Heatmap(as.matrix(rna_purist$TSP_wt), name = "TSP weights",
                      width = unit(4, "mm"), 
                      height = unit(7, "cm"),
                      row_title = "PurIST TSP weights", row_title_rot = 90,
                      show_row_names = F,
                      cluster_rows = F,
                      column_title = "Weights",
                      column_title_side = "bottom",
                      column_title_rot = 90,
                      column_title_gp = gpar(fontsize = 8),
                      col = colorRamp2(c(0,1.5), c("grey90", "black")),
                      use_raster = T,
                      heatmap_legend_param = list(border = "transparent", 
                                                  legend_direction = "horizontal",
                                                  title_position = "topleft",
                                                  legend_width = unit(1.5, "cm"),
                                                  at = c(0,1.5)))
rna_hm <- Heatmap(as.matrix(rna_purist[,7:ncol(rna_purist)]),
                  name = "Percentile rank",
                  cluster_rows = F, cluster_columns = F,
                  show_row_names = T, show_column_names = F, row_names_side = "right",
                  row_names_gp = gpar(fontsize = 9),
                  top_annotation = rna_column_ha,
                  col = cont_colors,
                  use_raster = T,
                  heatmap_legend_param = list(border = "transparent",
                                              legend_direction = "horizontal",
                                              title_position = "topleft",
                                              legend_width = unit(1.5, "cm")))

## tempus
tmp_column_ha <- HeatmapAnnotation(Tempus_PurIST = merged.pred$PurIST.x,
                                   RNAseq_PurIST = merged.pred$PurIST.y,
                                   col = sample_color2,
                                   simple_anno_size = unit(4, "mm"), gap = unit(1, "mm"),
                                   annotation_name_gp = gpar(fontsize = 8),
                                   show_annotation_name = T)
tmp_row_hm <- Heatmap(as.matrix(tmp_purist$Subtype),
                      row_title = "PurIST gene pairs", row_title_rot = 90,
                      show_row_names = F,
                      cluster_rows = F,
                      row_split = tmp_purist$TSP_wt,
                      col = gene_color,
                      show_heatmap_legend = F,
                      use_raster = T,
                      width = unit(4, "mm"),
                      height = unit(7, "cm"))
tmp_hm <- Heatmap(as.matrix(tmp_purist[,7:ncol(tmp_purist)]),
                  name = "Percentile rank",
                  cluster_rows = F, cluster_columns = F,
                  show_row_names = T, show_column_names = F, row_names_side = "right",
                  row_names_gp = gpar(fontsize = 9),
                  top_annotation = tmp_column_ha,
                  col = cont_colors,
                  use_raster = T,
                  heatmap_legend_param = list(border = "transparent",
                                              legend_direction = "horizontal",
                                              title_position = "topleft"))

## draw them side by side
draw(rna_row_hm + tsp_row_hm + rna_hm + tmp_hm,
     adjust_annotation_extension = T,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = T)

# sankey plots (figure 2a)
## extract expression matrix only from ordered data frame
r <- rna_purist[,7:ncol(rna_purist)] 
t <- tmp_purist[,7:ncol(tmp_purist)]

## do for loops to divide basal gene by classical gene for each tsp

## for rna-seq
rr <- nrow(r)
r_ratios <- list()
r_rownames <- list()

for(i in seq(1, rr, by = 2)) {
  res <- r[i,] / r[i+1,]
  r_ratios[[length(r_ratios) + 1]] <- res
  rm(i, res)
}

for(i in seq(1, rr, by = 2)) {
  names <- paste0(rownames(r[i,]), "-", rownames(r[i+1,]))
  r_rownames[[length(r_rownames) + 1]] <- names
  rm(i, names)
}

rna_ratios <- do.call(rbind, r_ratios)
rownames(rna_ratios) <- do.call(rbind, r_rownames)

## for tempus
tt <- nrow(t)
t_ratios <- list()
t_rownames <- list()

for(i in seq(1, tt, by = 2)) {
  res <- t[i,] / t[i+1,]
  t_ratios[[length(t_ratios) + 1]] <- res
  rm(i, res)
}

for(i in seq(1, tt, by = 2)) {
  names <- paste0(rownames(t[i,]), "-", rownames(t[i+1,]))
  t_rownames[[length(t_rownames) + 1]] <- names
  rm(i, names)
}

tmp_ratios <- do.call(rbind, t_ratios)
rownames(tmp_ratios) <- do.call(rbind, t_rownames)

## log2 the ratios so that data is centered around 0
log_rna_ratios <- log2(rna_ratios)
log_tmp_ratios <- log2(tmp_ratios)

## turn log b/c ratios into basal (>0) or classical (<0)
## there is one in rna-seq where the b/c ratio is 1 - will treat this is classical
rna_class <- t(log_rna_ratios)
tmp_class <- t(log_tmp_ratios)

rna_class <- ifelse(rna_class<=0, "Classical", "Basal-like")
tmp_class <- ifelse(tmp_class<0, "Classical", "Basal-like")

rna_class <- data.frame(rna_class)
tmp_class <- data.frame(tmp_class)

colnames(rna_class) <- paste0(colnames(rna_class), "_rna")
colnames(tmp_class) <- paste0(colnames(tmp_class), "_tmp")

merged_tsp <- merge(merged.pred, rna_class, by.x = "sampID.y", by.y = 0)
merged_tsp <- merge(merged_tsp, tmp_class, by.x = "sampID.x", by.y = 0)
rownames(merged_tsp) <- merged_tsp$ptID

## make long tables for each gene-pair
s100a2.long <- merged_tsp %>%
  make_long(S100A2.SLC40A1_rna, S100A2.SLC40A1_tmp) %>%
  mutate(node = factor(node, levels = c("Classical", "Basal-like")),
         next_node = factor(next_node, levels = c("Classical", "Basal-like")))

krt6a.long <- merged_tsp %>%
  make_long(KRT6A.ANXA10_rna, KRT6A.ANXA10_tmp) %>%
  mutate(node = factor(node, levels = c("Classical", "Basal-like")),
         next_node = factor(next_node, levels = c("Classical", "Basal-like")))

gpr87.long <- merged_tsp %>% filter(!is.na(GPR87.REG4_rna)) %>%
  make_long(GPR87.REG4_rna, GPR87.REG4_tmp) %>%
  mutate(node = factor(node, levels = c("Classical", "Basal-like")),
         next_node = factor(next_node, levels = c("Classical", "Basal-like")))

bcar3.long <- merged_tsp %>%
  make_long(BCAR3.GATA6_rna, BCAR3.GATA6_tmp) %>%
  mutate(node = factor(node, levels = c("Classical", "Basal-like")),
         next_node = factor(next_node, levels = c("Classical", "Basal-like")))

krt5.long <- merged_tsp %>%
  make_long(KRT5.CLRN3_rna, KRT5.CLRN3_tmp) %>%
  mutate(node = factor(node, levels = c("Classical", "Basal-like")),
         next_node = factor(next_node, levels = c("Classical", "Basal-like")))

ptges.long <- merged_tsp %>%
  make_long(PTGES.CLDN18_rna, PTGES.CLDN18_tmp) %>%
  mutate(node = factor(node, levels = c("Classical", "Basal-like")),
         next_node = factor(next_node, levels = c("Classical", "Basal-like")))

c16orf.long <- merged_tsp %>% filter(!is.na(C16orf74.DDC_rna)) %>%
  make_long(C16orf74.DDC_rna, C16orf74.DDC_tmp) %>%
  mutate(node = factor(node, levels = c("Classical", "Basal-like")),
         next_node = factor(next_node, levels = c("Classical", "Basal-like")))

itga3.long <- merged_tsp %>%
  make_long(ITGA3.LGALS4_rna, ITGA3.LGALS4_tmp) %>%
  mutate(node = factor(node, levels = c("Classical", "Basal-like")),
         next_node = factor(next_node, levels = c("Classical", "Basal-like")))

## sankey plots
ggplot(s100a2.long, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node))) +
  geom_sankey(node.color=NA, width=.4, smooth=7, space=3, flow.alpha=.5) +
  scale_fill_manual(values=c("Classical"="blue", "Basal-like"="orange")) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(size=10, hjust=.5),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size=12, hjust=.5),
        aspect.ratio=.75) +
  scale_x_discrete(labels = c("RNA-seq", "Tempus")) +
  labs(title = "S100A2 - SLC40A1 (weight: 1.505)")

ggplot(krt6a.long, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node))) +
  geom_sankey(node.color=NA, width=.4, smooth=7, space=3, flow.alpha=.5) +
  scale_fill_manual(values=c("Classical"="blue", "Basal-like"="orange")) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(size=10, hjust=.5),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size=12, hjust=.5),
        aspect.ratio=.75) +
  scale_x_discrete(labels = c("RNA-seq", "Tempus")) +
  labs(title = "KRT6A - ANXA10 (weight: 1.031)")

ggplot(gpr87.long, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node))) +
  geom_sankey(node.color=NA, width=.4, smooth=7, space=3, flow.alpha=.5) +
  scale_fill_manual(values=c("Classical"="blue", "Basal-like"="orange")) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(size=10, hjust=.5),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size=12, hjust=.5),
        aspect.ratio=.75) +
  scale_x_discrete(labels = c("RNA-seq", "Tempus")) +
  labs(title = "GPR87 - REG4 (weight: 0.994)")

ggplot(bcar3.long, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node))) +
  geom_sankey(node.color=NA, width=.4, smooth=7, space=3, flow.alpha=.5) +
  scale_fill_manual(values=c("Classical"="blue", "Basal-like"="orange")) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(size=10, hjust=.5),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size=12, hjust=.5),
        aspect.ratio=.75) +
  scale_x_discrete(labels = c("RNA-seq", "Tempus")) +
  labs(title = "BCAR3 - GATA6 (weight: 0.618)")

ggplot(krt5.long, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node))) +
  geom_sankey(node.color=NA, width=.4, smooth=7, space=3, flow.alpha=.5) +
  scale_fill_manual(values=c("Classical"="blue", "Basal-like"="orange")) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(size=10, hjust=.5),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size=12, hjust=.5),
        aspect.ratio=.75) +
  scale_x_discrete(labels = c("RNA-seq", "Tempus")) +
  labs(title = "KRT5 - CLRN3 (weight: 0.515)")

ggplot(ptges.long, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node))) +
  geom_sankey(node.color=NA, width=.4, smooth=7, space=3, flow.alpha=.5) +
  scale_fill_manual(values=c("Classical"="blue", "Basal-like"="orange")) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(size=10, hjust=.5),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size=12, hjust=.5),
        aspect.ratio=.75) +
  scale_x_discrete(labels = c("RNA-seq", "Tempus")) +
  labs(title = "PTGES - CLDN18 (weight: 0.078)")

ggplot(c16orf.long, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node))) +
  geom_sankey(node.color=NA, width=.4, smooth=7, space=3, flow.alpha=.5) +
  scale_fill_manual(values=c("Classical"="blue", "Basal-like"="orange")) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(size=10, hjust=.5),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size=12, hjust=.5),
        aspect.ratio=.75) +
  scale_x_discrete(labels = c("RNA-seq", "Tempus")) +
  labs(title = "C16orf74 - DDC (weight: 0.071)")

ggplot(itga3.long, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=factor(node))) +
  geom_sankey(node.color=NA, width=.4, smooth=7, space=3, flow.alpha=.5) +
  scale_fill_manual(values=c("Classical"="blue", "Basal-like"="orange")) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(size=10, hjust=.5),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size=12, hjust=.5),
        aspect.ratio=.75) +
  scale_x_discrete(labels = c("RNA-seq", "Tempus")) +
  labs(title = "ITGA3 - LGALS4 (weight = 0.059)")

# draw paired gene expression plots (figure s2)
## create data frame with rna-seq and tempus purist genes side by side
rna.purist.genes <- as.data.frame(t(rna.tsp.ranked))
colnames(rna.purist.genes) <- paste0(colnames(rna.purist.genes), "_rnaseq")
rna.purist.genes <- rna.purist.genes %>% 
  rownames_to_column("sampID")

tmp.purist.genes <- as.data.frame(t(tmp.tsp.ranked))
colnames(tmp.purist.genes) <- paste0(colnames(tmp.purist.genes), "_tempus")
tmp.purist.genes <- tmp.purist.genes %>%
  rownames_to_column("sampID")

tsp.paired <- merged_tsp %>% select(sampID.y, sampID.x, ptID, PurIST.y, PurIST_prob.y, PurIST.x, PurIST_prob.x) %>%
  left_join(rna.purist.genes, by = c("sampID.y" = "sampID")) %>%
  left_join(tmp.purist.genes, by = c("sampID.x" = "sampID"))

## GPR87 (basal)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = GPR87_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = GPR87_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = GPR87_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = GPR87_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = GPR87_rnaseq, yend = GPR87_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = GPR87_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = GPR87_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = GPR87_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = GPR87_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = GPR87_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = GPR87_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "GPR87", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$GPR87_rnaseq, tsp.paired$GPR87_tempus, paired = T)
wilcox.test(tsp.paired$GPR87_rnaseq, tsp.paired$GPR87_tempus, paired = T, alternative = "greater")

## KRT6A (basal)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = KRT6A_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = KRT6A_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = KRT6A_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = KRT6A_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = KRT6A_rnaseq, yend = KRT6A_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = KRT6A_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = KRT6A_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = KRT6A_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = KRT6A_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = KRT6A_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = KRT6A_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "KRT6A", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$KRT6A_rnaseq, tsp.paired$KRT6A_tempus, paired = T)
wilcox.test(tsp.paired$KRT6A_rnaseq, tsp.paired$KRT6A_tempus, paired = T, alternative = "greater")

## BCAR3 (basal)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = BCAR3_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = BCAR3_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = BCAR3_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = BCAR3_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = BCAR3_rnaseq, yend = BCAR3_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = BCAR3_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = BCAR3_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = BCAR3_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = BCAR3_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = BCAR3_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = BCAR3_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "BCAR3", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$BCAR3_rnaseq, tsp.paired$BCAR3_tempus, paired = T)
wilcox.test(tsp.paired$BCAR3_rnaseq, tsp.paired$BCAR3_tempus, paired = T, alternative = "greater")

## PTGES (basal)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = PTGES_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = PTGES_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = PTGES_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = PTGES_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = PTGES_rnaseq, yend = PTGES_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = PTGES_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = PTGES_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = PTGES_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = PTGES_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = PTGES_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = PTGES_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "PTGES", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$PTGES_rnaseq, tsp.paired$PTGES_tempus, paired = T)
wilcox.test(tsp.paired$PTGES_rnaseq, tsp.paired$PTGES_tempus, paired = T, alternative = "greater")

## ITGA3 (basal)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = ITGA3_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = ITGA3_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = ITGA3_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = ITGA3_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = ITGA3_rnaseq, yend = ITGA3_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = ITGA3_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = ITGA3_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = ITGA3_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = ITGA3_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = ITGA3_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = ITGA3_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "ITGA3", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$ITGA3_rnaseq, tsp.paired$ITGA3_tempus, paired = T)
wilcox.test(tsp.paired$ITGA3_rnaseq, tsp.paired$ITGA3_tempus, paired = T, alternative = "greater")

## C16orf74 (basal)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = C16orf74_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = C16orf74_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = C16orf74_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = C16orf74_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = C16orf74_rnaseq, yend = C16orf74_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = C16orf74_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = C16orf74_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = C16orf74_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = C16orf74_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = C16orf74_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = C16orf74_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "C16orf74", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$C16orf74_rnaseq, tsp.paired$C16orf74_tempus, paired = T)
wilcox.test(tsp.paired$C16orf74_rnaseq, tsp.paired$C16orf74_tempus, paired = T, alternative = "greater")

## S100A2 (basal)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = S100A2_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = S100A2_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = S100A2_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = S100A2_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = S100A2_rnaseq, yend = S100A2_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = S100A2_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = S100A2_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = S100A2_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = S100A2_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = S100A2_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = S100A2_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "S100A2", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$S100A2_rnaseq, tsp.paired$S100A2_tempus, paired = T)
wilcox.test(tsp.paired$S100A2_rnaseq, tsp.paired$S100A2_tempus, paired = T, alternative = "greater")

## KRT5 (basal)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = KRT5_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = KRT5_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = KRT5_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = KRT5_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = KRT5_rnaseq, yend = KRT5_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = KRT5_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = KRT5_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = KRT5_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = KRT5_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = KRT5_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = KRT5_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "KRT5", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$KRT5_rnaseq, tsp.paired$KRT5_tempus, paired = T)
wilcox.test(tsp.paired$KRT5_rnaseq, tsp.paired$KRT5_tempus, paired = T, alternative = "greater")

## REG4 (classical)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = REG4_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = REG4_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = REG4_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = REG4_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = REG4_rnaseq, yend = REG4_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = REG4_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = REG4_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = REG4_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = REG4_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = REG4_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = REG4_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "REG4", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$REG4_rnaseq, tsp.paired$REG4_tempus, paired = T)
wilcox.test(tsp.paired$REG4_rnaseq, tsp.paired$REG4_tempus, paired = T, alternative = "greater")

## ANXA10 (classical)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = ANXA10_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = ANXA10_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = ANXA10_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = ANXA10_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = ANXA10_rnaseq, yend = ANXA10_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = ANXA10_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = ANXA10_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = ANXA10_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = ANXA10_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = ANXA10_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = ANXA10_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "ANXA10", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$ANXA10_rnaseq, tsp.paired$ANXA10_tempus, paired = T)
wilcox.test(tsp.paired$ANXA10_rnaseq, tsp.paired$ANXA10_tempus, paired = T, alternative = "greater")

## GATA6 (classical)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = GATA6_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = GATA6_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = GATA6_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = GATA6_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = GATA6_rnaseq, yend = GATA6_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = GATA6_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = GATA6_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = GATA6_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = GATA6_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = GATA6_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = GATA6_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "GATA6", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$GATA6_rnaseq, tsp.paired$GATA6_tempus, paired = T)
wilcox.test(tsp.paired$GATA6_rnaseq, tsp.paired$GATA6_tempus, paired = T, alternative = "greater")

## CLDN18 (classical)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = CLDN18_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = CLDN18_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = CLDN18_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = CLDN18_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = CLDN18_rnaseq, yend = CLDN18_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = CLDN18_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = CLDN18_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = CLDN18_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = CLDN18_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = CLDN18_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = CLDN18_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "CLDN18", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$CLDN18_rnaseq, tsp.paired$CLDN18_tempus, paired = T)
wilcox.test(tsp.paired$CLDN18_rnaseq, tsp.paired$CLDN18_tempus, paired = T, alternative = "greater")

## LGALS4 (classical)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = LGALS4_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = LGALS4_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = LGALS4_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = LGALS4_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = LGALS4_rnaseq, yend = LGALS4_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = LGALS4_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = LGALS4_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = LGALS4_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = LGALS4_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = LGALS4_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = LGALS4_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "LGALS4", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$LGALS4_rnaseq, tsp.paired$LGALS4_tempus, paired = T)
wilcox.test(tsp.paired$LGALS4_rnaseq, tsp.paired$LGALS4_tempus, paired = T, alternative = "greater")

## DDC (classical)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = DDC_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = DDC_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = DDC_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = DDC_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = DDC_rnaseq, yend = DDC_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = DDC_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = DDC_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = DDC_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = DDC_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = DDC_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = DDC_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "DDC", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$DDC_rnaseq, tsp.paired$DDC_tempus, paired = T)
wilcox.test(tsp.paired$DDC_rnaseq, tsp.paired$DDC_tempus, paired = T, alternative = "greater")

## SLC40A1 (classical)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = SLC40A1_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = SLC40A1_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = SLC40A1_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = SLC40A1_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = SLC40A1_rnaseq, yend = SLC40A1_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = SLC40A1_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = SLC40A1_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = SLC40A1_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = SLC40A1_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = SLC40A1_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = SLC40A1_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "SLC40A1", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$SLC40A1_rnaseq, tsp.paired$SLC40A1_tempus, paired = T)
wilcox.test(tsp.paired$SLC40A1_rnaseq, tsp.paired$SLC40A1_tempus, paired = T, alternative = "greater")

## CLRN3 (classical)
ggplot(tsp.paired) +
  stat_boxplot(aes(x = 1, y = CLRN3_rnaseq), geom = "errorbar", color = "black", width = .4) +
  stat_boxplot(aes(x = 2, y = CLRN3_tempus), geom = "errorbar", color = "black", width = .4) +
  geom_boxplot(aes(x = 1, y = CLRN3_rnaseq), color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = CLRN3_tempus), color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = 1, xend = 2, y = CLRN3_rnaseq, yend = CLRN3_tempus), color = "grey50", alpha = .7, linewidth = .75) +
  geom_point(aes(x = 1, y = CLRN3_rnaseq), color = "grey30", shape = 16, size = 3) +
  geom_point(aes(x = 2, y = CLRN3_tempus), color = "grey30", shape = 16, size = 3) +
  stat_summary(aes(x = 1, y = CLRN3_rnaseq), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 2, y = CLRN3_tempus), geom = "crossbar", fun = mean, linewidth = .5, color = "darkred") +
  stat_summary(aes(x = 1, y = CLRN3_rnaseq), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  stat_summary(aes(x = 2, y = CLRN3_tempus), geom = "crossbar", fun = median, linewidth = .5, width = .8, color = "black") +
  scale_x_continuous(breaks = c(1, 2), expand = c(.1, .1), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(title = "CLRN3", x = "Method", y = "Relative expression")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if RNA-seq > Tempus
wilcox.test(tsp.paired$CLRN3_rnaseq, tsp.paired$CLRN3_tempus, paired = T)
wilcox.test(tsp.paired$CLRN3_rnaseq, tsp.paired$CLRN3_tempus, paired = T, alternative = "greater")

# plot purist probabilities (figure 1b)
merged_tsp$prob_change <- ifelse(merged_tsp$PurIST_prob.x > merged_tsp$PurIST_prob.y, "increase", "decrease")
merged_tsp$purist_change <- ifelse(merged_tsp$PurIST_prob.y < .5 & merged_tsp$PurIST_prob.x > .5, "yes", "no")
## add jitter in x dimension (more if purist prob under .2)
merged_tsp <- merged_tsp %>%
  mutate(x1_jitter = ifelse(PurIST_prob.y > .2, jitter(rep(1, n()), amount = .005),
                            jitter(rep(1, n()), amount = .05)),
         x2_jitter = ifelse(PurIST_prob.x > .2, jitter(rep(2, n()), amount = .005),
                            jitter(rep(2, n()), amount = .05)))
## add jitter in y dimension (more if purist prob under .5)
merged_tsp <- merged_tsp %>%
  mutate(y1_jitter = ifelse(PurIST_prob.y > .5, jitter(PurIST_prob.y, amount = .01),
                            jitter(PurIST_prob.y, amount = .02)),
         y2_jitter = ifelse(PurIST_prob.x > .5, jitter(PurIST_prob.x, amount = .01),
                            jitter(PurIST_prob.x, amount = .02)))

ggplot(merged_tsp) +
  stat_boxplot(aes(x = 1, y = PurIST_prob.y), geom = "errorbar", color = "black", width = .2) +
  stat_boxplot(aes(x = 2, y = PurIST_prob.x), geom = "errorbar", color = "black", width = .2) +
  geom_boxplot(aes(x = 1, y = PurIST_prob.y), coef = 0, color = "black", fill = "#7aa6dc", outlier.shape = NA) +
  geom_boxplot(aes(x = 2, y = PurIST_prob.x), coef = 0, color = "black", fill = "#cd534c", outlier.shape = NA) +
  geom_segment(aes(x = x1_jitter, xend = x2_jitter, y = y1_jitter, yend = y2_jitter, color = purist_change), linewidth = .75) +
  geom_point(aes(x = x1_jitter, y = y1_jitter, fill = PurIST.y), color = "black", shape = 21, size = 3, stroke = .5) +
  geom_point(aes(x = x2_jitter, y = y2_jitter, fill = PurIST.x), color = "black", shape = 21, size = 3, stroke = .5) +
  scale_color_manual(values = c("no" = "grey30", "yes" = "orange")) +
  scale_fill_manual(values = c("Classical" = "blue", "Basal-like" = "orange")) +
  scale_x_continuous(breaks = c(1, 2), expand = c(0, .25), labels = c("RNA-seq", "Tempus")) +
  scale_y_continuous(limits = c(-.02,1), expand = c(.02,0)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=10, angle=0, hjust=.5, vjust=.5, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio=1.5) +
  labs(x = "Method", y = "Basal-like probability")

## wilcoxon signed-rank test for paired data: first one is any change, second is specifically if Tempus > RNA-seq
wilcox.test(merged_tsp$PurIST_prob.y, merged_tsp$PurIST_prob.x, paired = T)
wilcox.test(merged_tsp$PurIST_prob.x, merged_tsp$PurIST_prob.y, paired = T, alternative = "greater")

# scatter plot (figure 1c)
ggplot(merged_tsp, aes(x = PurIST_prob.y, y = PurIST_prob.x)) +
  geom_point(aes(color = PurIST.y, fill = PurIST.x), shape = 21, stroke = 1, size = 4) +
  scale_color_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.75)) +
  scale_fill_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.6)) +
  theme(panel.background=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=.75),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=14),
        axis.title=element_text(color="black", size=16),
        legend.position="none",
        aspect.ratio=1) +
  scale_x_continuous(expand = c(.01,.01), limits = c(0,1)) +
  scale_y_continuous(expand = c(.01,.01), limits = c(0,1)) +
  geom_vline(xintercept=.5, linetype="dashed", color="grey30") +
  geom_hline(yintercept=.5, linetype="dashed", color="grey30") +
  labs(x = "RNA-seq basal-like probability", y = "Tempus basal-like probability")

# look at probability difference between rna-seq and tempus (figure 1d)
## difference named residuals here; actual minus predicted
merged_tsp$residuals <- merged_tsp$PurIST_prob.y - merged_tsp$PurIST_prob.x # y is rna-seq and x is tempus

ggplot(merged_tsp, aes(x = PurIST_prob.y, y = residuals)) +
  geom_point(aes(color = PurIST.y, fill = PurIST.x), shape = 21, stroke = 1, size = 4, alpha = .7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_color_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.75)) +
  scale_fill_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.6)) +
  theme(panel.background=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=.75),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=14),
        axis.title=element_text(color="black", size=16),
        legend.position="none",
        aspect.ratio=.6) +
  scale_x_continuous(expand = c(.01, .01), limits = c(0,1)) +
  scale_y_continuous(expand = c(.01, .01), limits = c(-1,.25)) +
  labs(x = "RNA-seq basal-like probability", y = "Probability difference")

# use estimate to infer tumor purity (figure s1a)
## LOWER ESTIMATE SCORE IS HIGHER PURITY
rna.estimate <- rna |> 
  filter_common_genes(id = "hgnc_symbol", tell_missing = F, find_alias = F) |> 
  estimate_score(is_affymetrix = F)
colnames(rna.estimate) <- paste0(colnames(rna.estimate), "_rnaseq")

tmp.estimate <- tmp |>
  filter_common_genes(id = "hgnc_symbol", tell_missing = F, find_alias = F) |>
  estimate_score(is_affymetrix = F)
colnames(tmp.estimate) <- paste0(colnames(tmp.estimate), "_tempus")

merged_tsp <- merged_tsp %>%
  left_join(rna.estimate, by = c("sampID.y" = "sample_rnaseq")) %>%
  left_join(tmp.estimate, by = c("sampID.x" = "sample_tempus"))

## plot estimate score vs probability difference
ggplot(merged_tsp, aes(x = estimate_rnaseq, y = residuals)) +
  geom_smooth(method = "lm", se = T, color = "grey50", fullrange = T) +
  geom_point(aes(color = PurIST.y, fill = PurIST.x), shape = 21, stroke = 1, size = 4, alpha = .7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_color_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.75)) +
  scale_fill_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.6)) +
  theme(panel.background=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=.75),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=14),
        axis.title=element_text(color="black", size=16),
        legend.position="none",
        aspect.ratio=.6) +
  scale_x_reverse(expand = c(0,.01), breaks = seq(5000, -3000, by = -1000), limits = c(5100,-3100)) +
  scale_y_continuous(expand = c(0,.01), limits = c(-1,.25)) +
  labs(x = "ESTIMATE score from RNA-seq", y = "Probability difference")

cat("Pearson's correlation between ESTIMATE score from RNA-seq and probability differences")
cor.test(merged_tsp$estimate_rnaseq, merged_tsp$residuals)

# estimate tumor purity using decoder (figure s1b)
## prep data
log.rna <- log2(rna+1)
refSet <- "TCGA_RNAseq_PAAD"

## run decoder
sampleWeights <- Decon_single_sample(refSet,
                                     log.rna,
                                     "geneSymbol")
sampleWeights <- Norm_PDAC_weights(sampleWeights) # bcSum is tumor purity
sampleWeights <- sampleWeights %>%
  rownames_to_column("sample")

## add decoder purity to data frame
merged_tsp <- merged_tsp %>%
  left_join(sampleWeights %>% dplyr::select(sample, bcSum), by = c("sampID.y" = "sample"))

ggplot(merged_tsp, aes(x = bcSum, y = residuals)) +
  geom_smooth(method = "lm", se = T, color = "grey50", fullrange = T) +
  geom_point(aes(color = PurIST.y, fill = PurIST.x), shape = 21, stroke = 1, size = 4, alpha = .7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_color_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.75)) +
  scale_fill_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.6)) +
  theme(panel.background=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=.75),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=14),
        axis.title=element_text(color="black", size=16),
        legend.position="none",
        aspect.ratio=.6) +
  scale_x_continuous(expand = c(0,0), limits = c(.19,.81)) +
  scale_y_continuous(expand = c(0,.01), limits = c(-1,.25)) +
  labs(x = "DECODER tumor purity", y = "Probability difference")

cat("Pearson's correlation between DECODER score from RNA-seq and probability differences")
cor.test(merged_tsp$bcSum, merged_tsp$residuals)


# incorporate pancreas trial subtype matches
## load clia sheet that has combined nanostring purist probability and rna-seq purist probability for clinical trial samples
clia <- read.csv("path_to_clia_data.csv", na.strings = "", stringsAsFactors = F)

# confusion matrix (figure 4a)
clia <- clia %>%
  filter(!is.na(Profile), !is.na(PurIST)) %>%
  mutate(Profile = factor(Profile, levels = c("Classical", "Basal-like")),
         PurIST = factor(PurIST, levels = c("Classical", "Basal-like")))
confusionMatrix(clia$Profile, clia$PurIST, mode = "prec_recall", positive = "Basal-like")

## to get confidence interval for kappa since caret does not report it
Kappa.test(clia$Profile, clia$PurIST)

# paired basal-like probability plot between rna-seq and nanostring (figure 4b)
clia$purist_change <- ifelse(clia$PurIST==clia$Profile, "no", 
                             ifelse(clia$PurIST=="Classical" & clia$Profile=="Basal-like", "to_basal", "to_classical"))
## add random noise to x dimension (more jitter if basal prob < .05)
clia <- clia %>%
  mutate(x1_jitter = ifelse(PurIST_prob_RNA > .05, jitter(rep(1, n()), amount = .005),
                            jitter(rep(1, n()), amount = .05)),
         x2_jitter = ifelse(PurIST_score_CLIA > .05, jitter(rep(2, n()), amount = .005),
                            jitter(rep(2, n()), amount = .05)))
## add random noise to y dimension (more jitter if basal prob < .05)
clia <- clia %>%
  mutate(y1_jitter = ifelse(PurIST_prob_RNA > .05, jitter(PurIST_prob_RNA, amount = .01), 
                            jitter(PurIST_prob_RNA, amount = 0.025)),
         y2_jitter = ifelse(PurIST_score_CLIA > .05, jitter(PurIST_score_CLIA, amount = .01), 
                            jitter(PurIST_score_CLIA, amount = 0.025)))

ggplot(clia) +
  geom_segment(aes(x = x1_jitter, xend = x2_jitter, y = y1_jitter, yend = y2_jitter, color = purist_change), linewidth = .75) +
  geom_point(aes(x = x1_jitter, y = y1_jitter, fill = PurIST), color = "black", shape = 21, size = 4, stroke = .5) +
  geom_point(aes(x = x2_jitter, y = y2_jitter, fill = Profile), color = "black", shape = 21, size = 4, stroke = .5) +
  stat_summary(aes(x = 1, y = PurIST_prob_RNA), geom = "crossbar", fun = median, width = .5, color = "red2") +
  stat_summary(aes(x = 2, y = PurIST_score_CLIA), geom = "crossbar", fun = median, width = .5, color = "red2") +
  scale_color_manual(values = c("no" = "grey30", "to_basal" = "orange", "to_classical" = "blue")) +
  scale_fill_manual(values = c("Classical" = "blue", "Basal-like" = "orange")) +
  scale_x_continuous(breaks = c(1, 2), expand = c(0,.25), labels = c("RNA-seq", "NanoString")) +
  scale_y_continuous(limits = c(-.03,1.02), expand = c(.01,.01)) +
  theme(legend.position="none",
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(color="black", size=14),
        axis.ticks=element_line(color="black"),
        legend.key=element_blank(),
        aspect.ratio = 1.7) +
  labs(x = "Method", y = "Basal-like probability")

## wilcoxon signed-rank test for paired data
wilcox.test(clia$PurIST_prob_RNA, clia$PurIST_score_CLIA, paired = T)

# probability difference between rna-seq and nanostring vs percent malignancy (figure 4c)
## difference named residuals here; actual minus predicted
clia$residuals <- clia$PurIST_prob_RNA - clia$PurIST_score_CLIA 

ggplot(clia, aes(x = Avg_malig, y = residuals)) +
  geom_smooth(method = "lm", se = T, color = "grey50", fullrange = T) +
  geom_point(aes(color = PurIST, fill = Profile), shape = 21, stroke = 1.2, size = 4.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_color_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.75)) +
  scale_fill_manual(values = alpha(c("Classical"="blue", "Basal-like"="orange"),.5)) +
  theme(panel.background=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=.75),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=14),
        axis.title=element_text(color="black", size=16),
        legend.position="none",
        aspect.ratio=.6) +
  scale_x_continuous(expand = c(0,0), limits = c(-2,102)) +
  scale_y_continuous(expand = c(.01, .01), limits = c(-.75,.75), breaks = seq(-.75,.75,.25)) +
  labs(x = "Percent malignant cells", y = "Probability difference")

## pearson's correlation between percent malignant cells and probability differences
cor.test(clia$Avg_malig, clia$residuals)

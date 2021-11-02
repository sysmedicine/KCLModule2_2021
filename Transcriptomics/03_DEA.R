



# @ Differential Expression Analysis ----
# convert to DESeqDataSet format
dset <- DESeqDataSetFromMatrix(countData = round(mat.count),
                               colData   = mat.design, 
                               design    = ~ Treatment + Time)
# @ DESeq2
ds2 <- DESeq(dset, fitType = "parametric", parallel = T)
# dispersion plot (whether the model fits)
plotDispEsts(ds2)

# @ Results ----
# PCA plot - by use of VST function
dset.transformed <- varianceStabilizingTransformation(ds2, blind = F, fitType = "parametric")
plotPCA(dset.transformed, intgroup = c("Time"), ntop = 1000) + 
  stat_ellipse(linetype=1,level=0.95) + 
  # coord_fixed(xlim=c(-50,50), ylim=c(-40,40)) + 
  # theme(text = element_text(size = 35))
  theme_bw()

# Differentially Expressed Genes (DEGs)
summary(DESeq2::results(ds2, cooksCutoff = F, alpha = 0.01, lfcThreshold = 0, contrast = c("Treatment", "MI", "sham")))
res = DESeq2::results(ds2, cooksCutoff = F, alpha = 0.01, lfcThreshold = 0, contrast = c("Treatment", "MI", "sham"))
res = data.frame(res)

summary(DESeq2::results(ds2, cooksCutoff = F, alpha = 0.01, lfcThreshold = 0, contrast = c("Time", "3dMI", "1dMI")))
res.2 = DESeq2::results(ds2, cooksCutoff = F, alpha = 0.01, lfcThreshold = 0, contrast = c("Time", "3dMI", "1dMI"))
res.2 = data.frame(res.2)

# Volcano plot
res.interest = res.2[res.2$padj < 0.0001,]
res.interest = res.interest[abs(res.interest$log2FoldChange) > 2,]

ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = padj < 0.05)) +
  geom_point() + 
  # geom_text_repel(
  #   data = res.interest,
  #   aes(label = res.interest, 
  #       # fill = ifelse(log2FoldChange < 0, "66FFFF", "FF3300")),
  #       fill = ifelse(padj > 0.05, "CCCCCC",ifelse(log2FoldChange < 0, "66FFFF", "FF3300"))),
  #   # fill = ifelse(padj < 0.05, "CCCCCC",ifelse(log2FoldChange < 0, "66FFFF", "FF3300"))),
  #   fontface = 'bold', color = 'white',
  #   # box.padding = unit(box_padding, "lines"),
  #   # point.padding = unit(point_padding, "lines"),
  #   # nudge_y = nudge_y,
  #   nudge_x = ifelse(res.interest[["log2FoldChange"]] < 0, -1, 1),
  #   segment.color = 'grey50'
  # ) +
  # xlim(c(-2.5,2.5)) +
  # ylim(c(0,15)) +
  guides(fill = F) +
  theme_bw()

# Heatmap 
cor.tpm = cor(mat.tpm)
p.heatmap = pheatmap(cor.tpm, fontsize = 14)
p.heatmap

tiff("pheatmap.tiff", height = 1800, width = 1800, res = 300)
p.heatmap
dev.off()



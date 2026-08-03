###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~ age-related DNA methylation sites (DMS) ~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

library(doParallel)
library(foreach)
library(splines)
library(fastcluster)
library(missMethyl)
library(scales)
library(mixOmics)
library(ggplot2)
library(sva)
library(arm)
library(ggradar)
library(patchwork)
registerDoParallel(cores = 20)


###########################step1: Data Overview
#################PCA
df_pca <- pca(X = t(betaSaliva.m), ncomp = 2, scale = TRUE, center = TRUE)

sampleInfo <- as.data.frame(PhenoTypesSaliva.lv)
sampleInfo$PC1 <- df_pca$x[, "PC1"]
sampleInfo$PC2 <- df_pca$x[, "PC2"]
rownames(sampleInfo) <- colnames(betaSaliva.m)

ggplot(sampleInfo, aes(x = PC1, y = PC2, colour = Population, size = Age)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black", size = 11),
        axis.title = element_text(colour = "black", size = 12),
        axis.ticks.length = unit(0.15, "cm")) +
  labs(x = paste0("PC1: ", percent(df_pca$prop_expl_var$X[1], accuracy = 0.01)),
       y = paste0("PC2: ", percent(df_pca$prop_expl_var$X[2], accuracy = 0.01))) +
  scale_color_manual(values = c("#F7CCC6", "#c2dbcf", "#CACDE8"))


##################PVCA
source("PVCA.R")
PVCA(exprSet = betaSaliva.m, meta = sampleInfo[, 1:8])
dev.off()

############remove YNs_a13(Yi) and YNs_b13 (Yi)
betaSaliva.m.rmYi <- betaSaliva.m[, !colnames(betaSaliva.m) %in% c("YNs_a13", "YNs_b13")]
sampleInfo <- sampleInfo[!rownames(sampleInfo) %in% c("YNs_a13", "YNs_b13"), ]

############batch effect correction related to Sample_Plate
betaSaliva.m.rmYi.logit <- logit(betaSaliva.m.rmYi)
betaSaliva.m.logit.adj <- ComBat(dat = betaSaliva.m.rmYi.logit, batch = sampleInfo$Sample_Plate)
betaSaliva.m.adj <- invlogit(betaSaliva.m.logit.adj)

#########PCA with corrected beta value
df_pca <- pca(X = t(betaSaliva.m.adj), ncomp = 2, scale = TRUE, center = TRUE)

sampleInfo$PC1 <- df_pca$x[, "PC1"]
sampleInfo$PC2 <- df_pca$x[, "PC2"]

ggplot(sampleInfo, aes(x = PC1, y = PC2, fill = Age, size = Age)) +
     geom_point(shape = 21, colour = "#b3c5cf", stroke = 0.6) +
     theme_bw() +
     theme(panel.grid.minor = element_blank(),
           axis.text = element_text(colour = "black", size = 11),
           axis.title = element_text(colour = "black", size = 12),
           axis.ticks.length = unit(0.15, "cm"),
           legend.position = "none") +
     labs(x = paste0("PC1: ", percent(df_pca$prop_expl_var$X[1], accuracy = 0.01)),
          y = paste0("PC2: ", percent(df_pca$prop_expl_var$X[2], accuracy = 0.01))) +
     scale_fill_gradient(low = "#93a2a3", high = "#274546")

#########PVCA with corrected beta value
PVCA(exprSet = betaSaliva.m.adj, meta = sampleInfo[, 1:8])

####ggradar
data <- matrix(c("Saliva", 25.184, 12.31, 2.582, 2.167, 1.717, 1.645, 1.03, 0.116), nrow = 1)
colnames(data) <- c("group", "Age", "Sentrix ID", "Sentrix Position", "Sample Well", "Sex", "Population", "Twins Type", "Sample Plate")
data <- as.data.frame(data)
data[, -1] <- as.numeric(data[, -1])

p1 <- ggradar(data, background.circle.colour = "white",
              values.radar = "",
              group.line.width = 0.7,
              group.point.size = 1.5, 
              group.colours = "#a9adbc",
              fill = TRUE,
              fill.alpha = 0.4,
              axis.label.size = 4,
              axis.label.offset = 1.11,
              axis.line.colour = "grey",
              grid.line.width = 0.2, 
              gridline.min.linetype = "solid",
              gridline.mid.linetype = "solid",
              gridline.max.linetype = "solid",
              gridline.min.colour = "grey",
              gridline.mid.colour = "grey",
              gridline.max.colour = "black") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0,0,0,0))

tmp <- data.frame(x = rep(1, 8), y = rep(1.5, 8), 
                  group = paste0("g", 1:8))

anno <- ggplot()+
  geom_bar(data = tmp, aes(x, y, fill = group), stat = "identity", position = "dodge") +
  scale_fill_manual(values = alpha(c("#503d5c", "#624b6e", "#6f597a", "#8b7991","#af92b5", "#c7b1db", "#cdc0d6", "#fcfbf9"), 0.7))+
  ylim(-6, 1.5) +
  coord_polar(start = -0.43) +
  theme_void() +
  theme(legend.position = "none",
        axis.ticks.length = unit(0, "pt"),
        axis.text = element_blank(),
        plot.margin = margin(0,0,0,0))

anno + inset_element(p1, left = -0.04, bottom = -0.06, right = 1.03, top = 1.05)

sampleInfo <- sampleInfo[, 1:8]

############################step2: Age test
source("function_define.r")

Age <- sampleInfo$Age

age_DMS <- age.test.polynomials.1.f(
  Mat = betaSaliva.m.adj, Ages = Age, verbose = F, pvals = T, PCUTOFF = 0.05)

permut <- foreach(i = 1:1000, .combine = "cbind") %dopar% {
  
  AGE <- Age[sample(1:length(Age), length(Age), replace = FALSE)]
  
  age_polyreg <- age.test.polynomials.1.f(
    Mat = betaSaliva.m.adj, Ages = AGE, verbose = F, pvals = T, PCUTOFF = 0.05)
  
  return(age_polyreg)
}


#### Estimate the false discovery rate
FP <- median(apply(permut[, seq(2, ncol(permut), 2)], 2, function(x){sum(x < 0.01)}))
FDR <- FP/sum(age_DMS[, 2] < 0.01)
 
saliva_age <- betaSaliva.m.adj[rownames(betaSaliva.m.adj) %in% rownames(age_DMS)[age_DMS[, 2] < 0.01], ]
 
m <- matrix(c(1, 1, 2, 3, 4, 5), 2, 3, byrow = T); layout(m);
par(mar = c(3, 3, 3, 1), cex.lab = 1.2, cex.axis = 1.1, mgp = c(1, 0.2, 0), tck = -0.02)
 
#### Hierarchical clustering
dist <- 1 - cor(t(saliva_age), method = "pearson")
clust <- fastcluster::hclust(as.dist(dist), method = "complete")
 
### Figure A: Hierarchical dendrogram
plot(clust, labels = FALSE, hang = -1, lwd = 0.1, main = NA)
 
### Figure B: Number of CpG clusters obtained by cutting the hierarchical clustering
### tree at different heights
height <- seq(1.55, 1.95, 0.1)
 
clusters <- c()
CpG_IN_clusters <- c()
for(i in 1:length(height)){
  cIc = table(cutree(clust, h = height[i]))
  CpG_IN_clusters[[i]] = cIc
  bk = length(cIc)
  clusters = append(clusters, bk)
}
 
barplot(clusters, names.arg = height, col = "#9baba2", border = NA,
         xlab = "Height cutoff", ylab = "Number of clusters", cex.names = 1.1)
 
### Figure C: Number of genes annotated by CpG sites in each of the clusters identified
### using given cutting height cutoff
AnnotatedGene <- list()
for(i in 1:length(height)){
  group = cutree(clust, h = height[i])
  sig.eg = lapply(unique(group), function(grp){
     tryCatch(getMappedEntrezIDs(sig.cpg = names(group[group == grp]),
                                 all.cpg = rownames(saliva_age),
                                 array.type = "EPIC",
                                 genomic.features = c("TSS200", "TSS1500", "1stExon", "5'UTR"))$sig.eg,
              error = function(e) NA)
   })
  AnnotatedGene[[i]] = sig.eg
}
 
df <- lapply(AnnotatedGene, function(AG){unlist(lapply(AG, length))})
boxplot(df, names = height, xlab = "Height cutoff", ylab = "Number of annotated genes",
         col = "#9798ae", border = "#8b7991", lwd = 1)
 
### Figure D: Number of GO terms and KEGG pathways significantly overrepresented
### among genes in each of the CpG clusters
Enrich <- list()
for(i in 1:length(height)){
   group = cutree(clust, h = height[i])
   gst = lapply(unique(group), function(grp){
     gst.fdr = tryCatch(gometh(sig.cpg = names(group[group == grp]),
                               all.cpg = rownames(saliva_age),
                               array.type = "EPIC",
                               genomic.features = c("TSS200", "TSS1500", "1stExon", "5'UTR"),
                               collection = "GO",
                               prior.prob = FALSE),
                        error = function(e) NA)
  })
  Enrich[[i]] = gst
}

df <- lapply(Enrich, function(EN){
   count = unlist(lapply(EN, function(x){sum(x[6] < 0.05)}));
   prop = round(sum(count != 0 & !is.na(count))/length(EN), digits = 3);
   prop
})
 
barplot(unlist(df), names.arg = height, col = "#beb0bf", border = NA, cex.names = 1.1,
        xlab = "Height cutoff", ylab = "Proportion of enriched clusters")
 
term <- c()
for(i in 1:length(Enrich)){
  sig = c()
  for(j in 1:length(Enrich[[i]])){
    df = Enrich[[i]][[j]]
    sig = append(sig, rownames(df)[df[6] < 0.05 & df[1] == "BP"])
  }
  term <- c(term, length(unique(sig)))
}
 
barplot(term, names.arg = height, col = "#d2b6be", border = NA, cex.names = 1.1,
         xlab = "Height cutoff", ylab = "Number of enriched terms")
 
source("A2Rplot.R")
A2Rplot(clust, k = 5, boxes = FALSE, col.up = "grey60", lty.up = 1,
         col.down = c("#547980", "#45ADA8", "#b2ac88", "#FFD3B5", "#E84A5F"),
         show.labels = FALSE, lwd.up = 0.7, lwd.down = 0.7)
 
memb <- cutree(clust, k = 5)
 
ageRanges <- range(Age)
ageGrid <- seq(from = ageRanges[1], to = ageRanges[2], length.out = 30)
 
cl <- 5
 
plot.data <- saliva_age[rownames(saliva_age) %in% names(memb)[memb == cl], ]
plot.data <- scale(t(plot.data), center = TRUE, scale = TRUE)
plot.data <- data.frame(cg = rowMeans(plot.data), age = Age)
 
fit <- lm(cg ~ ns(age, df = 3), plot.data)
pred <- predict(fit, newdata = list(age = ageGrid), interval = "confidence")
 
ggplot() +
  geom_point(data = plot.data, aes(x = age, y = cg), size = 4,
             colour = "#695980", shape = 21, stroke = 0.5) +
  geom_line(aes(x = ageGrid, y = pred[, 1]), linewidth = 1, colour = "#695980") +
  geom_ribbon(aes(ymin = pred[, 2], ymax = pred[, 3], x = ageGrid),
               fill = alpha("#695980", 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 11),
        plot.title = element_text(colour = "black", size = 12,
                                  hjust = 0.5, face = "bold"),
        axis.ticks.length = unit(0.2, "cm")) +
  scale_y_continuous(limits = c(-2.1, 2)) +
  labs(x = "Age (years)", y = "Methylation level",
      title = paste(paste0("C", cl), sum(memb == cl), sep = ":"))

########Aging atlas
ageGene <- read.csv("Aging atlas.csv")
ageGene <- ageGene[!duplicated(ageGene), ]

bg <- getMappedEntrezIDs(sig.cpg = rownames(betaSaliva.m.adj),
                         array.type = "EPIC",
                         genomic.features = c("TSS200", "TSS1500", "1stExon", "5'UTR"))$sig.eg

ageGene <- ageGene[ageGene$Gene_ID %in% bg, ]

annoGene <- getMappedEntrezIDs(sig.cpg = rownames(saliva_age),
                               array.type = "EPIC",
                               genomic.features = c("TSS200", "TSS1500", "1stExon", "5'UTR"))$sig.eg

##########hypergeometric test
N <- length(unique(bg))
M <- length(unique(ageGene$Gene_ID))
n <- length(unique(annoGene))
k <- length(intersect(annoGene, ageGene$Gene_ID))

phyper(k-1, M, N-M, n, lower.tail = FALSE)

enrich <- c()
for(terms in unique(ageGene$Gene_Set)){
  termGene = ageGene$Gene_ID[ageGene$Gene_Set == terms]
  M = length(termGene)
  k = length(intersect(annoGene, termGene))
  enrich = rbind(enrich, c(terms, M, k, round(phyper(k-1, M, N-M, n, lower.tail = FALSE), 4)))
}
colnames(enrich) <- c("Category", "CategoryGene", "hit", "pval")
enrich <- as.data.frame(enrich)
enrich$padj <- p.adjust(enrich$pval, method = "BH")

m <- matrix(1:2, 1, 2, byrow = TRUE);layout(m)
pie(c(nrow(saliva_age), nrow(betaSaliva.m.adj) - nrow(saliva_age)),
    col = c("#a9bcd0", "#d8dbe2"), border = NA, labels = c("Age-related", ""))

bp <- barplot(as.numeric(enrich[, 2]), names.arg = enrich[, 1], border = NA,
              col = "#906a44", las = 2, horiz = TRUE)

for(i in 1:length(bp)){
  rect(ybottom = bp[i] - 0.5, xleft = 0, xright = enrich[i, 3],
       ytop = bp[i] + 0.5, border = NA, col = "#684121")
}

############alternative
ageGene <- read.csv("Aging atlas.csv")
ageGene <- ageGene[!duplicated(ageGene), ]

bg <- getMappedEntrezIDs(sig.cpg = rownames(saliva_age),
                         array.type = "EPIC",
                         genomic.features = c("TSS200", "TSS1500", "1stExon", "5'UTR"))$sig.eg

ageGene <- ageGene[ageGene$Gene_ID %in% bg, ]

annoGene <- AnnotatedGene[[4]]

N <- length(unique(bg))

enrich <- c()
for(i in 1:length(annoGene)){
  for(terms in unique(ageGene$Gene_Set)){
    termGene = ageGene$Gene_ID[ageGene$Gene_Set == terms]
    M = length(termGene)
    k = length(intersect(annoGene[[i]], termGene))
    n = length(unique(annoGene[[i]]))
    enrich = rbind(enrich, c(paste0("C", i), terms, M, k, round(phyper(k-1, M, N-M, n, lower.tail = FALSE), 4)))
  }
}

enrich <- as.data.frame(enrich)
colnames(enrich) <- c("Group", "Term", "Ncount", "Hit", "pval")
enrich$prop <- round(as.numeric(enrich$Hit)/as.numeric(enrich$Ncount), 3)
enrich$sig <- ifelse(enrich$pval < 0.05, "*", "")
enrich$Term <- factor(enrich$Term, levels = unique(ageGene$Gene_Set)[order(unlist(lapply(unique(ageGene$Gene_Set), nchar)), decreasing = FALSE)])
enrich$Group <- factor(enrich$Group, levels = paste0("C", c(1, 4, 3, 2)))

ggplot(enrich, aes(x = Group, y = Term)) +
  geom_tile(aes(fill = prop), colour = NA) +
  scale_fill_gradient2(low = "#BBDEEC",high = "#695980", breaks = c(seq(0, 1, 0.2))) +
  geom_text(aes(label = sig), colour = "black", vjust = 0.5, hjust = 0.5, size = 8) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.ticks = element_blank())


####################################################################
#####~~~~~~~~ Differential methylated position (DMP) ~~~~~~~~~~#####
##############~~~~~~~~~ Monozygotic twins ~~~~~~~~~~~###############
############~~~~~~~~~~~ Dizygotic twins ~~~~~~~~~~~~~###############

library(DSS)
library(bsseq)
library(gridExtra)
library(ggplot2)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(GenomicRanges)
library(karyoploteR)
library(ggbreak)
library(doParallel)
registerDoParallel(cores = 20)


###############Step1:Differentially methylated position (DMP)
###########Shuffling the sample labels within similar age, generating undifferentiated distribution of Δβ 
metadata <- data.frame(SampleID = rownames(sampleInfo),
                       AgeGroup = cut(sampleInfo$Age, breaks = c(0, 20, 40, 60, 80), 
                                      labels = c("0-20", "20-40", "40-60", "60-80")))

permut_deltas <- foreach(i = 1:1000, .combine = "cbind") %dopar% {
  
  permut_beta = matrix(nrow = nrow(betaSaliva.m.adj), ncol = ncol(betaSaliva.m.adj))
  
  for(age_group in unique(metadata$AgeGroup)){
    
    idx = which(metadata$AgeGroup == age_group)
    
    if(length(idx) >= 4){
      
      twin1 = idx[seq(1, length(idx), 2)]
      twin2 = idx[seq(2, length(idx), 2)]
      
      repeat{
        shuffle = sample(twin2, size = length(twin2), replace = FALSE)
        if(all(shuffle != twin2)) break
      }
      
      permut_beta[, twin1] = betaSaliva.m.adj[, twin1]
      permut_beta[, twin2] = betaSaliva.m.adj[, shuffle]
      
    }
    
  }
  
  deltas = abs(permut_beta[, seq(1, ncol(permut_beta), 2)] - permut_beta[, seq(2, ncol(permut_beta), 2)])
  deltas = apply(deltas, 1, quantile, 0.95)
  return(deltas)
  
}

save(permut_deltas, file = "permut_saliva_deltas_within_ageGroup.RData")

cutoff <- apply(permut_deltas, 1, quantile, 0.95)

twin_deltas <- abs(betaSaliva.m.adj[, seq(1, ncol(betaSaliva.m.adj), 2)] - betaSaliva.m.adj[, seq(2, ncol(betaSaliva.m.adj), 2)])

MZ_deltas <- twin_deltas[, colnames(twin_deltas) %in% rownames(sampleInfo)[sampleInfo$Twins_Type == "MZ"]]
DMP_MZ <- apply(MZ_deltas, 2, function(f){f > cutoff})
save(DMP_MZ, file = "Saliva_DMP_MZ.RData")

DZ_deltas <- twin_deltas[, colnames(twin_deltas) %in% rownames(sampleInfo)[sampleInfo$Twins_Type == "DZ"]]
DMP_DZ <- apply(DZ_deltas, 2, function(f){f > cutoff})


###############Step2:Differentially methylated region (DMR)
source("DMRfinder.R")
DMR_MZ <- lapply(1:ncol(DMP_MZ), function(f){
  
  data = data.frame(Probe_ID = rownames(DMP_MZ),
                    Status = as.numeric(DMP_MZ[, f]),
                    stringsAsFactors = FALSE)
  
  dmrs = DMRfinder(data, mismatches = 3, icd = 10000, illumina = TRUE)
  
  dmrs
})

###############Step3:Differentially methylated CpG island (dmCGI)
threshold <- apply(permut_deltas, 1, quantile, 0.95)
twin_CGIs <- abs(CGI_beta[, seq(1, ncol(CGI_beta), 2)] - CGI_beta[, seq(2, ncol(CGI_beta), 2)])

MZ_CGIs <- twin_CGIs[, colnames(twin_CGIs) %in% rownames(sampleInfo)[sampleInfo$Twins_Type == "MZ"]]
dmCGI_MZ <- apply(MZ_CGIs, 2, function(f){f > threshold})

DZ_CGIs <- twin_CGIs[, colnames(twin_CGIs) %in% rownames(sampleInfo)[sampleInfo$Twins_Type == "DZ"]]
dmCGI_DZ <- apply(DZ_CGIs, 2, function(f){f > threshold})

#############CGIs including at least one DMPs
CGI <- probe_islands[probe_islands$Relation_to_Island == "Island" & probe_islands$Probe_ID %in% rownames(betaSaliva.m.adj), ]
dmCGI_in_DMPs <- lapply(1:ncol(dmCGI_MZ), function(f){
  
  dmCGI = rownames(dmCGI_MZ)[dmCGI_MZ[, f] == TRUE]
  dmps = rownames(DMP_MZ)[DMP_MZ[, f] == TRUE]
  
  CID = unique(CGI$Islands_Name[CGI$Islands_Name %in% dmCGI & CGI$Probe_ID %in% dmps])
  
  CID
})


bp <- barplot(colSums(dmCGI_MZ), las = 2, col = alpha("#b2c9c9", 0.3), border = NA,
              cex.axis = 1.1, cex.names = 1.1)
for(i in 1:length(bp)){
  rect(xleft = bp[i] - 0.5, ybottom = 0, xright = bp[i] + 0.5, 
       ytop = unlist(lapply(dmCGI_in_DMPs, length))[i], angle = 45, 
       border = "#a5bcc1", density = 20, col = "#a5bcc1")
}
text(x = bp, y = colSums(dmCGI_MZ), labels = colSums(dmCGI_MZ), adj = 0.5,
     pos = 3, cex  = 1, offset = 0.3)

points(x = bp, y = c(2, 68, 2, 830, 3, 37, 15, 288, 0, 2, 1, 6, 8, 118, 3, 13, 6,
                     97, 109, 78, 0, 4851, 3, 127, 21, 2, 91, 377, 28, 340, 34, 12)
         , type = "b", col = "#b17769", cex = 1.5)


###############Step4:Summarize
###########DMP counts along with Age
withinPair_DMP_counts <- colSums(DMP_MZ)
dat <- data.frame(DMP_count = withinPair_DMP_counts,
                  #DMR_count = unlist(lapply(DMR_MZ, function(f){ifelse(is.null(f), 0, nrow(f))})),
                  Age = sampleInfo$Age[match(names(withinPair_DMP_counts), rownames(sampleInfo))],
                  stringsAsFactors = FALSE)

ggplot(dat, aes(x = Age, y = DMP_count, size = DMP_count, fill = Age)) +
  geom_point(shape = 21, colour = "#cc9b6d", stroke = 0.5) +
  #geom_smooth(method = "gam", colour = "#f1ca89", fill = alpha("#f1ca89", 0.3), 
              #formula = y ~ x + I(x^2), linetype = 2, linewidth = 1) +
  geom_smooth(method = "lm", colour = "#f1ca89", fill = alpha("#f1ca89", 0.3),
              linetype = 2, linewidth = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 11),
        axis.title = element_blank(),
        axis.ticks.length = unit(0.15, "cm"), 
        legend.position = "none") +
  scale_fill_gradient(low = "#f7e9d7", high = "#cc9b6d") +
  labs(y = "Number of DMPs", x = "Age (years)") +
  scale_y_break(breaks = c(40000, 60000), scales = 0.3, space = 0.1)


###########recurrent DMP (rDMP) and sporadic DMP (sDMP)
acrossPair_DMP_counts <- rowSums(DMP_MZ)
acrossPair_DMP_counts <- acrossPair_DMP_counts[acrossPair_DMP_counts > 0]
table(acrossPair_DMP_counts == 1)/length(acrossPair_DMP_counts)*100


dat <- data.frame(group = c("sDMP", "rDMP"),
                  value = c(173308, 44299),
                  percent = c(0.8, 0.2))

ggplot(dat, aes(x = 2, y = percent, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ebd8c3", "#f9f3df"))


CpG_island <- probe_islands[probe_islands$Probe_ID %in% names(acrossPair_DMP_counts)[acrossPair_DMP_counts > 0], ]
CpG_island$Count <- acrossPair_DMP_counts[match(CpG_island$Probe_ID, names(acrossPair_DMP_counts))]
CpG_island$Group <- ifelse(CpG_island$Count == 1, "sporadic DMP", "recurrent DMP")


###########Relation_to_Island
xtab <- table(CpG_island$Relation_to_Island, CpG_island$Group)

bp <- barplot(t(xtab)/rowSums(t(xtab)), beside = TRUE, legend.text = c("rDMP", "sDMP"),
              col = c("#ebd8c3", "#f9f3df"), border = NA, ylab = "Proportion (%)",
              cex.axis = 1.2, cex.names = 1.1, cex.lab = 1.2, xpd = TRUE, las = 2,
              args.legend = list(x = "topleft", border = NA, cex = 1.1, bty = "n"))
text(x = bp[1, ], y = xtab[, 1]/sum(xtab[, 1]), xpd = TRUE,
     labels = percent(xtab[, 1]/sum(xtab[, 1]), accuracy = 0.1),
     cex = 1, adj = 0.5, pos = 3, offset = 0.2)
text(x = bp[2, ], y = xtab[, 2]/sum(xtab[, 2]), xpd = TRUE,
     labels = percent(xtab[, 2]/sum(xtab[, 2]), accuracy = 0.1),
     cex = 1, adj = 0.5, pos = 3, offset = 0.2)


pvalue <- c()
for(i in 1:nrow(xtab)){
  
  col1 = sum(xtab[, 1])
  col2 = sum(xtab[, 2])
  m = matrix(c(xtab[i, 1], xtab[i, 2], col1 - xtab[i, 1], col2 - xtab[i, 2]),
             2, 2, byrow = TRUE)
  
  if(xtab[i, 1]/col1 < xtab[i, 2]/col2){
    pval = fisher.test(m, alternative = "less")$p.value
  } else {
    pval = fisher.test(m, alternative = "greater")$p.value
  }
  
  pvalue = append(pvalue, pval)
}

names(pvalue) <- rownames(xtab)


###########UCSC_RefGene_Group
new.data <- CpG_island[CpG_island$UCSC_RefGene_Name != "", c(1, 7, 8, 10)]
split_data <- lapply(1:nrow(new.data), function(i){
  genes <- unlist(strsplit(new.data[i, "UCSC_RefGene_Name"], ";"))
  region <- unlist(strsplit(new.data[i, "UCSC_RefGene_Group"], ";"))
  data.frame(Probe_ID = rep(new.data[i, "Probe_ID"], length(genes)),
             RefGene_Name = genes, RefGene_Group = region, 
             Group = rep(new.data[i, "Group"], length(genes)))
})

split_data <- do.call(rbind, split_data)
split_data <- split_data[!duplicated(split_data), ]

xtab <- table(split_data$RefGene_Group, split_data$Group)

bp <- barplot(t(xtab)/rowSums(t(xtab)), beside = TRUE, legend.text = c("rDMP", "sDMP"),
              col = c("#ebd8c3", "#f9f3df"), border = NA, ylab = "Proportion (%)",
              cex.axis = 1.2, cex.names = 1.1, cex.lab = 1.2, xpd = TRUE, las = 2,
              args.legend = list(x = "topleft", border = NA, cex = 1.1, bty = "n"))
text(x = bp[1, ], y = xtab[, 1]/sum(xtab[, 1]), xpd = TRUE,
     labels = percent(xtab[, 1]/sum(xtab[, 1]), accuracy = 0.1),
     cex = 1, adj = 0.5, pos = 3, offset = 0.2)
text(x = bp[2, ], y = xtab[, 2]/sum(xtab[, 2]), xpd = TRUE,
     labels = percent(xtab[, 2]/sum(xtab[, 2]), accuracy = 0.1),
     cex = 1, adj = 0.5, pos = 3, offset = 0.2)


pvalue <- c()
for(i in 1:nrow(xtab)){
  
  col1 = sum(xtab[, 1])
  col2 = sum(xtab[, 2])
  m = matrix(c(xtab[i, 1], xtab[i, 2], col1 - xtab[i, 1], col2 - xtab[i, 2]),
             2, 2, byrow = TRUE)
  
  if(xtab[i, 1]/col1 < xtab[i, 2]/col2){
    pval = fisher.test(m, alternative = "less")$p.value
  } else {
    pval = fisher.test(m, alternative = "greater")$p.value
  }
  
  pvalue = append(pvalue, pval)
}

names(pvalue) <- rownames(xtab)



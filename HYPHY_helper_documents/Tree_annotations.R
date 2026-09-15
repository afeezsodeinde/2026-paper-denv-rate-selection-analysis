setwd("~/Desktop/DENV_EVO_2025/DENV4")
library(ape)
tree <- read.tree("capsid_protein_C-nodups.fasta.treefile")
write(tree$tip.label, file = "capsid_tips.csv")
tree <- read.tree("envelope_protein_E-nodups.fasta.treefile")
write(tree$tip.label, file = "envelope_tips.csv")
tree <- read.tree("membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")
write(tree$tip.label, file = "membrane_tips.csv")
tree <- read.tree("nonstructural_protein_NS1-nodups.fasta.treefile")
write(tree$tip.label, file = "NS1.csv")
tree <- read.tree("nonstructural_protein_NS2A-nodups.fasta.treefile")
write(tree$tip.label, file = "NS2A.csv")
tree <- read.tree("nonstructural_protein_NS2B-nodups.fasta.treefile")
write(tree$tip.label, file = "NS2B.csv")
tree <- read.tree("nonstructural_protein_NS3-nodups.fasta.treefile")
write(tree$tip.label, file = "NS3.csv")
tree <- read.tree("nonstructural_protein_NS4A-nodups.fasta.treefile")
write(tree$tip.label, file = "NS4A.csv")
tree <- read.tree("nonstructural_protein_NS4B-nodups.fasta.treefile")
write(tree$tip.label, file = "NS4B.csv")
tree <- read.tree("RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")
write(tree$tip.label, file = "NS5.csv")

setwd("~/Desktop/DENV_EVO_2025/DENV1")
library(ape)
library(ggtree)
library(dplyr)

tree <- read.tree("envelope_protein_E-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/envelope_tips_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_envelope_protein_E-nodups.fasta.treefile")

#Capsid-DENV1
tree <- read.tree("capsid_protein_C-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/capsid_tips_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_capsid_protein_C-nodups.fasta.treefile")

#Membrane-DENV1
tree <- read.tree("membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/membrane_tips_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")

#NS1-DENV1
tree <- read.tree("nonstructural_protein_NS1-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/NS1_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_nonstructural_protein_NS1-nodups.fasta.treefile")

#NS2A-DENV1
tree <- read.tree("nonstructural_protein_NS2A-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/NS2A_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_nonstructural_protein_NS2A-nodups.fasta.treefile")

#NS2B-DENV1
tree <- read.tree("nonstructural_protein_NS2B-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/NS2B_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_nonstructural_protein_NS2B-nodups.fasta.treefile")

#NS3-DENV1
tree <- read.tree("nonstructural_protein_NS3-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/NS3_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_nonstructural_protein_NS3-nodups.fasta.treefile")

#NS4A-DENV1
tree <- read.tree("nonstructural_protein_NS4A-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/NS4A_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_nonstructural_protein_NS4A-nodups.fasta.treefile")

#NS4B-DENV1
tree <- read.tree("nonstructural_protein_NS4B-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/NS4B_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_nonstructural_protein_NS4B-nodups.fasta.treefile")

#NS5-DENV1
tree <- read.tree("RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV1/Tree_Tips_csv/NS5_with_maj_lins.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV1/Annotated_trees/annotated_RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")

###DENV2
#Capsid-DENV2
setwd("~/Desktop/DENV_EVO_2025/DENV2")
tree <- read.tree("capsid_protein_C-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/capsid_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_capsid_protein_C-nodups.fasta.treefile")

#Envelope-DENV2
setwd("~/Desktop/DENV_EVO_2025/DENV2")
tree <- read.tree("envelope_protein_E-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/envelope_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_envelope_protein_E-nodups.fasta.treefile")

#Membrane-DENV2
setwd("~/Desktop/DENV_EVO_2025/DENV2")
tree <- read.tree("membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/membrane_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")

#NS1-DENV2
tree <- read.tree("nonstructural_protein_NS1-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/NS1.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_nonstructural_protein_NS1-nodups.fasta.treefile")

#NS2A-DENV2
tree <- read.tree("nonstructural_protein_NS2A-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/NS2A.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_nonstructural_protein_NS2A-nodups.fasta.treefile")

#NS2B-DENV2
tree <- read.tree("nonstructural_protein_NS2B-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/NS2B.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_nonstructural_protein_NS2B-nodups.fasta.treefile")

#NS4A-DENV2
tree <- read.tree("nonstructural_protein_NS4A-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/NS4A.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_nonstructural_protein_NS4A-nodups.fasta.treefile")

#NS4B-DENV2
tree <- read.tree("nonstructural_protein_NS4B-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/NS4B.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_nonstructural_protein_NS4B-nodups.fasta.treefile")

#NS5-DENV2
tree <- read.tree("RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV2/Tree_Tips_csv/NS5.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV2/Annotated_trees/annotated_RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")


###DENV3
#Capsid-DENV3
setwd("~/Desktop/DENV_EVO_2025/DENV3")
tree <- read.tree("capsid_protein_C-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/capsid_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_capsid_protein_C-nodups.fasta.treefile")

#Envelope-DENV3
tree <- read.tree("envelope_protein_E-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/envelope_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_envelope_protein_E-nodups.fasta.treefile")

#Membrane-DENV3
tree <- read.tree("membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/membrane_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")

#NS1-DENV3
tree <- read.tree("nonstructural_protein_NS1-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/NS1.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_nonstructural_protein_NS1-nodups.fasta.treefile")

#NS2A-DENV3
tree <- read.tree("nonstructural_protein_NS2A-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/NS2A.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_nonstructural_protein_NS2A-nodups.fasta.treefile")

#NS2B-DENV3
tree <- read.tree("nonstructural_protein_NS2B-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/NS2B.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_nonstructural_protein_NS2B-nodups.fasta.treefile")

#NS3-DENV3
tree <- read.tree("nonstructural_protein_NS3-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/NS3.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_nonstructural_protein_NS3-nodups.fasta.treefile")

#NS4A-DENV3
tree <- read.tree("nonstructural_protein_NS4A-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/NS4A.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_nonstructural_protein_NS4A-nodups.fasta.treefile")

#NS4B-DENV3
tree <- read.tree("nonstructural_protein_NS4B-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/NS4B.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_nonstructural_protein_NS4B-nodups.fasta.treefile")

#NS5-DENV3
tree <- read.tree("RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV3/Tree_Tips_csv/NS5.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV3/Annotated_trees/annotated_RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")

###DENV4
#Capsid-DENV4
setwd("~/Desktop/DENV_EVO_2025/DENV4")
tree <- read.tree("capsid_protein_C-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/capsid_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_capsid_protein_C-nodups.fasta.treefile")

#Envelope-DENV4
tree <- read.tree("envelope_protein_E-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/envelope_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_envelope_protein_E-nodups.fasta.treefile")

#Membrane-DENV4
tree <- read.tree("membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/membrane_tips.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_membrane_glycoprotein_precursor_prM-nodups.fasta.treefile")

#NS1-DENV4
tree <- read.tree("nonstructural_protein_NS1-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/NS1.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_nonstructural_protein_NS1-nodups.fasta.treefile")

#NS2A-DENV4
tree <- read.tree("nonstructural_protein_NS2A-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/NS2A.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_nonstructural_protein_NS2A-nodups.fasta.treefile")

#NS2B-DENV4
tree <- read.tree("nonstructural_protein_NS2B-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/NS2B.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_nonstructural_protein_NS2B-nodups.fasta.treefile")

#NS3-DENV4
tree <- read.tree("nonstructural_protein_NS3-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/NS3.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_nonstructural_protein_NS3-nodups.fasta.treefile")

#NS4A-DENV4
tree <- read.tree("nonstructural_protein_NS4A-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/NS4A.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_nonstructural_protein_NS4A-nodups.fasta.treefile")

#NS4B-DENV4
tree <- read.tree("nonstructural_protein_NS4B-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/NS4B.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_nonstructural_protein_NS4B-nodups.fasta.treefile")

#NS5-DENV4
tree <- read.tree("RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")
meta <- read.csv("~/Desktop/DENV_EVO_2025/DENV4/Tree_Tips_csv/NS5.csv", stringsAsFactors = FALSE)
lut  <- setNames(meta$Genotype, meta$Sequence)
lab  <- ifelse(tree$tip.label %in% names(lut),
               paste0(tree$tip.label, " (", lut[tree$tip.label], ")"),
               tree$tip.label)

plot(tree, cex = 0.7)
tiplabels(text = lab, frame = "none", adj = -0.05, cex = 0.6, offset = 0.001)


tree2 <- tree
tree2$tip.label <- ifelse(tree2$tip.label %in% names(lut),
                          paste0(tree2$tip.label, "_", lut[tree2$tip.label]),
                          tree2$tip.label)
write.tree(tree2, "~/Desktop/DENV_EVO_2025/DENV4/Annotated_trees/annotated_RNA-dependent_RNA_polymerase_NS5-nodups.fasta.treefile")



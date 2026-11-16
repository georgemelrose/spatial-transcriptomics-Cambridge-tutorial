### 11.1 Setup# 

# Load libraries
library(Seurat) # single-cell and spatial analysis toolkit
library(sparseMatrixStats) #high-performance summary stats 
#for rows and columns of sparse matrices 
library(paletteer) # colour palettes
library(ggplot2) # plotting
library(dplyr) # data manipulation
library(patchwork) # combining plots

# Load the Seurat object from the previous chapter
visium <- readRDS("precomputed/mouse_brain_visium_deconvolution.rds")

# Load CellChat library
# an R package for inference, analysis, and visualisation of cell-cell
# communication from SC-data
library(CellChat)

### 11.2 Exercise Overview & Purpose ###

#Cells in multi-cellular tissues communicate to coordinate development, immune
#responses and disease progression. A common mechanism is ligands from one cell to
#receptors in another#

#Receptor-mediated pathways can then influence gene expression, cell behaviour, and
#cell fate decisions.#

#Single-cell, spatial transcriptomics let us study these interactions by measuring
#gene expression across many cells at once. Cell-cell communication analysis uses
#curate ligand-receptor pairs to infer which signalling relationships are most likely
#between cell populations. 

# Cells with a ligand = potential senders
# Cells with a receptor = potential receivers

# Therefore, build an interaction network! 

# CCC methods don't discover new signalling pathways - they rely
# on curated databases to ID the interactions most plausible in the data.

# Why use spatial transcriptomics? Because we can then ask how cell
# types communicate with the tissue micro-environment, not just whether 
# they have compatible ligands and receptors. 

# There are lots of different methods for the different datasets - 
# https://link.springer.com/article/10.1186/s13059-026-04063-5





### 11.3 CellChat Analysis ##

#An R package for inferring cell-cell communication networks from cell,
#spatial transcriptomic data. It combines gene expression measurements with a 
#specialist database of ligand-receptor interactions to estimate the probability
#of communication between cell populations.

#CellChat accounts for mulit-subunit receptor complexes, cofactors, pathway level
#organisation - improves performance compared to models only using ligand-receptor pairs

#ANALYSIS STEPS#

#1) Create the CellChat object

#2) Set the ligand-receptor interaction database 

#3) ID variable genes and interactions

#4) Compute communication probabilities and infer signalling pathways 

#5) Visualise and explore the results 

# View scale factors for the image in the Seurat object
# scale factors are numeric ratios that map pixel coordinates from
# a full-resolution histology image to downsample images, for plotting

#e.g. a map legend, zoom ratio that keep data points lined up with 
#your microscope photo.

ScaleFactors(visium[["slice1"]])

# Create spatial factors data frame
spatial_factors <- data.frame(
  ratio = 55 / ScaleFactors(visium[["slice1"]])$spot,
  tol = 55 / 2
)

spatial_factors

# Extract tissue coordinates of the spots
spatial_coords <- GetTissueCoordinates(visium, image = "slice1")[, c("x", "y")]

#AAACAAGTATCTCCCA-1 isa standard barcode/spot identifier 
#found in spatial transcriptomic datasets#
head(spatial_coords)

# compute distances of the spots to each other
spot_dists <- computeCellDistance(
  coordinates = spatial_coords, #each row gives spatial locations of each cell
  ratio = spatial_factors$ratio, #conversion factor when converting spatial coords
                                #to pixels
  tol = spatial_factors$tol #tolerance factor, for increasing the robustness
                            #when comparing center-to-center distance 
)

# minimum distance observed
# remove zero distances (distance of a spot to itself)
min(spot_dists[spot_dists != 0])

# Create CellChat object
cell_chat <- createCellChat(
  object = visium, #normalised data matrix (not count) of genes by cells
  assay = "SCT", #assay to use when input is a Seurat or SingleCellExperiment 
                 #object
  group.by = "first_type", #the object in the metadata that defines cell groups
  datatype = "spatial", #by default its' RNA
  coordinates = spatial_coords, #data matrix in which each row gives spatial locations
  spatial.factors = spatial_factors #
)

# drop levels of cell type labels that are not present in the data
cell_chat@idents <- droplevels(cell_chat@idents)

# Visualise categories of interactions in the mouse database
# shows the different mouse database interactions
showDatabaseCategory(CellChatDB.mouse)

# Subset database to exclude "Non-protein Signalling"
mouse_db <- subsetDB(
  CellChatDB.mouse,
  search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
  non_protein = FALSE
)

# add the database to the CellChat object
cell_chat@DB <- mouse_db

# subset the expression data of signaling genes for saving computation cost
cell_chat <- subsetData(cell_chat)

# Identify overexpressed genes and interactions
# this is essentially a differential expression test between the labels defined above
# similar to cluster marker genes analysis. 

#this is done to ID interactions where both partners are variable 
#the steps focus the analysis on genes, interactions more likely to contribute 
#to communication between cell populations

#these steps focus the analysis on the the genes, interactions more 
#likely to contribute to communication between cell populations
cell_chat <- identifyOverExpressedGenes(cell_chat)
cell_chat <- identifyOverExpressedInteractions(cell_chat)


# Compute communication probabilities
# computes the communication probability between interacting cell groups

cell_chat <- computeCommunProb(
  cell_chat, #object for cell-cell interaction
  contact.range = 100, #interaction range to restrict the 
                       #contact-dependent signaling rate. This matches
                       #the approx.distance between the 2 centres of 2 adjacent 
                       #10X Visium spots. For SC resolution technologies,
                       #it's recommended to use the average cell diameter
  interaction.range = 200 #the max.interaction/diffusion length of ligands
)

# Compute communication probabilities at the signaling pathway level
# this function aggregates probabilities of individual ligand-receptor
# interactions into pathway-level probabilities, using pathway structure
cell_chat <- computeCommunProbPathway(cell_chat)

# Filter interactions to only consider those supported by at least 10 spots
cell_chat <- filterCommunication(cell_chat, min.cells = 10)

# Aggregate the communication network at the cell type level
cell_chat <- aggregateNet(cell_chat)

# calculate number of cells of each type for scaling the node sizes in the plot
group_size <- as.numeric(table(cell_chat@idents))

# Visualise the overall communication network as a circle plot
netVisual_circle(
  cell_chat@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions",
  top = 0.1
)

# Visualise the weights of the interactions
netVisual_circle(
  cell_chat@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength",
  top = 0.1
)

# Visualise the interaction strength between cell types as a heatmap
netVisual_heatmap(cell_chat, measure = "weight")


# Bubble plot of ligand-receptor interactions from Lamp5 to Vip and vice versa
netVisual_bubble(
  cell_chat,
  sources.use = c("Lamp5", "Vip"),
  targets.use = c("Lamp5", "Vip")
)

# barplot of the contribution of each ligand-receptor interaction pair to the selected pathways
p1 <- netAnalysis_contribution(cell_chat, signaling = "COLLAGEN") +
  ggtitle("COLLAGEN pathway")

p2 <- netAnalysis_contribution(cell_chat, signaling = "LAMININ") +
  ggtitle("LAMININ pathway")

p1 +
  p2 +
  plot_annotation(title = "Contribution of ligand-receptor interaction pairs")

# compute centrality scores for the selected pathways
cell_chat <- netAnalysis_computeCentrality(cell_chat)

# visualise the centrality scores for a selected pathway
netAnalysis_signalingRole_network(
  cell_chat,
  signaling = "LAMININ",
  width = 16,
  height = 8
)

# All pathways
netAnalysis_signalingRole_scatter(cell_chat)

# Only the LAMININ pathway
netAnalysis_signalingRole_scatter(cell_chat, signaling = c("LAMININ"))

netVisual_aggregate(
  cell_chat,
  signaling = "LAMININ",
  layout = "spatial",
  edge.width.max = 2,
  vertex.size.max = 1,
  alpha.image = 0.2,
  vertex.label.cex = 3.5
)

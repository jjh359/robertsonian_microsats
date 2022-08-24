#Code to generate population statistics for Robertsonian house mice using 
#microsat data and the adegenet package. 
#Informed by the adegenet manuals by Thibaut Jombart and tutorials by Tom Jenkins
#https://tomjenkins.netlify.app/2020/09/21/r-popgen-getting-started/#5

##############################################################
# NOTE: This "script" will analyze one set of loci at a time. 
# Having it iterate over all sets of loci is a WIP. 
# !!! You will also need to manually adjust how bootstraps
# are displayed - see end of script !!!
# PS Just realised drawSupportOnEdges() may do this automatically
##############################################################

library("ade4") ### Required for adegenet and for PCAs
library("adegenet") ### Used to manipulate genind and genpop files
library("ape") ### Used for drawing neighbour joining trees
library("hierfstat") ### Tests for F-statistics and genetic distances
library("genepop") ### Used for calculating Fst and Rst
library("ggplot2") ### Usd for plotting
library("graph4lg") ### Used to convert genind to genepop format
library("reshape2") ### Used for melt function
library("Rphylip") ### Used for phylogenetics, specifically the hybrid tree
library("stringr") ### Used for manipulating data with str_squish()
library("tidyverse") ### Used for pivot_longer function
library("viridis") ### Used for colour scheme

###############################################################
# Lines you need to change!  ##################################
###############################################################

#Set working directory to wherever you have your input files
setwd("D:/Documents/Research/Robertsonian_microsats/Rb_microsats_final_version/")

#Name CSV file in working directory containing a set of loci to run analyses on
microsat_data <- read.csv("all_loci.csv")

#Choose a suffix that will be added to output file names to indicate which loci are being analyzed
suffix <- "all_loci"

#Choose a phrase to express which loci are represented. Similar to above, but used
#for titles in some plots so you can have spaces. E.g. "loci on chromosome 6"
pca_suffix <- "all loci"

#Path to exe directory in your local phylip installation
path_to_fitch <- "D:/Programs/phylip-3.698/exe"

#Choose number of bootstraps to run during phylogenetic analyses
nboots <- 1000

#Filename for the standardised Dmyu csv file (as edited from RSTCalc output) that matches set of loci being analysed
dmyu_file <- "mouse_standardised_dmyu_all_loci.csv"

#######################################################
# Data preparation  ###################################
#######################################################

#Subset genotype, individual, and site columns in preparation for genind/genpop files
genotypes <- microsat_data[,3:ncol(microsat_data)]
ind <- as.character(microsat_data$ID)
site <- as.character(microsat_data$Site)

#Create genind and genpop variables
mouse_gen <- df2genind(genotypes, ploidy = 2, ind.names = ind, pop = site, sep = "/")
mouse_pop <- genind2genpop(mouse_gen)

#Get lists of loci names and population names
loci_names <- locNames(mouse_gen)
populations <- popNames(mouse_pop)

#Produce table with allelic richness
rich <- allelic.richness(genind2hierfstat(mouse_gen), diploid = TRUE)$Ar
mean_rich <- apply(rich, MARGIN = 2, FUN = mean)
rich <- round(rich, digits = 3)
rich <- rbind(rich, mean_rich)
rownames(rich) <- c(loci_names, "Mean")
write.table(rich, file = paste0("allelic_richness_", suffix, ".txt"), sep = ",", col.names = NA)

#######################################################
# Heterozygosity ######################################
#######################################################

#Get basic stats
basics_mouse <- basic.stats(mouse_gen, diploid = TRUE)

#Observed heterozygosity
Ho_mouse <- apply(basics_mouse$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE)
round(Ho_mouse, digits = 2)
Ho_mouse

#Expected heterozygosity
He_mouse <- apply(basics_mouse$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE)
round(He_mouse, digits = 2)
He_mouse

#Visualize
Het_mouse <- data.frame(Site = names(Ho_mouse), Ho = Ho_mouse, He = He_mouse)
Het_mouse <- melt(Het_mouse, id.vars = "Site")
Het_mouse$Site <- factor(Het_mouse$Site, levels=populations)
ggplot(data = Het_mouse, aes(x = Site, y = value, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge(width = 1), colour = "black")+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.80))+
  scale_fill_viridis(discrete=TRUE)+
  ylab("Heterozygosity")+
  xlab("Chromosomal race")+
  #ggtitle(paste0("Observed (Ho) and expected (He) heterozygosity for ", pca_suffix))+
  theme_bw() +
  theme(legend.title= element_blank(), legend.position = c(0.1,0.9))
ggsave(paste0("heterozygosity_", suffix, ".pdf"), plot = last_plot(), device = "pdf", width = 6, height = 6)

#######################################################
# Fst and Rst #########################################
#######################################################

#Note that this will generate additional files named according to the suffix you selected
#Convert our genind to genepop and use it to calculate Fst and Rst
genind_to_genepop(mouse_gen, output = paste0("mouse_genepop_", suffix, ".txt"))
Fst(paste0("mouse_genepop_", suffix, ".txt"), pairs = TRUE, outputFile = paste0("mouse_genepop_", suffix, "_fst.txt")) #Also produces .MIG
Fst(paste0("mouse_genepop_", suffix, ".txt"), pairs = TRUE, sizes = TRUE, outputFile = paste0("mouse_genepop_", suffix, "_rst.txt")) #Also produces .MIG

#Get a matrix of Fst (lower triangle) against Rst (upper triangle)
#Involves a little bit of finagling to get the data out of the files produced by genepop
#If you change the number of populations you will need to change these parameters accordingly! 
fst_rst_matrix <- matrix(0, 8, 8)
population_names <- c("CHPO", "IBIN", "ICRE", "IGAL", "ILVA", "ST1", "ST2", "STP")
dimnames(fst_rst_matrix) <- list(population_names,population_names)
fst_lines <- readLines(paste0("mouse_genepop_", suffix, "_fst.txt.MIG"))
for (i in 1:7){
  fst_lines[i+3] <- str_squish(fst_lines[i+3])
  fst_value <- as.numeric(strsplit(fst_lines[i+3], split=" ")[[1]])
  fst_rst_matrix[i+1,1:i] <- fst_value
}
rst_lines <- readLines(paste0("mouse_genepop_", suffix, "_rst.txt.MIG"))
for (i in 1:7){
  rst_lines[i+3] <- str_squish(rst_lines[i+3])
  rst_value <- as.numeric(strsplit(rst_lines[i+3], split=" ")[[1]])
  fst_rst_matrix[1:i,i+1] <- rst_value
}

#Convert our now populated matrix into a dataframe for easier visualization
fst_rst_df <- fst_rst_matrix %>% 
  as.data.frame() %>%
  rownames_to_column("site1") %>%
  pivot_longer(-c(site1), names_to = "site2", values_to = "Fst/Rst")
fst_rst_df$`Fst/Rst` <- round(fst_rst_df$`Fst/Rst`, digits = 3)

#Will plot Rst above diagonal and Fst below
fst_label <- expression(italic("R")[ST] / italic("F")[ST])
mid <- max(fst_rst_df$`Fst/Rst`) / 2
ggplot(data = fst_rst_df, aes(x=site1, y=site2, fill=`Fst/Rst`))+
  geom_tile(color = "black")+
  geom_text(aes(label = `Fst/Rst`), color="black", size = 4) +
  scale_fill_viridis(name = fst_label, limits = c(0, max(fst_rst_df$`Fst/Rst`)), breaks = c(0, 0.25, 0.5), direction = 1)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "left")+
  xlab("Chromosomal race 1")+
  ylab("Chromosomal race 2")+
  theme_bw()
ggsave(paste0("Rst_Fst_", suffix, ".pdf"), plot = last_plot(), device = "pdf", width = 8, height = 8)
#######################################################
# PCAs ################################################
#######################################################

#Get allele counts from genind
mouse_pca <- tab(mouse_gen, NA.method = "mean")

# Perform PCA
pca1 <- dudi.pca(mouse_pca, scannf = FALSE, scale = FALSE, nf = 3)

#Analyse what percent of genetic variance is explained by each axis
percent <- pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,14),
        names.arg = round(percent, 1))

#Create data frame of coordinates, and add IDs and sites
ind_coords <- as.data.frame(pca1$li)
colnames(ind_coords) <- c("Axis1","Axis2","Axis3")
ind_coords$Ind <- indNames(mouse_gen)
ind_coords$Site <- mouse_gen$pop

#Calculate centroid position for each population
centroid <- aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

#Add centroid coordinates to ind_coords data frame
ind_coords <- left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

#x and y labels, plus title
xlab <- paste("Principal Component 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab <- paste("Principal Component 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")
pca_title <- paste0("PCA from ", pca_suffix)

#First two Principal Components
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, fontface="bold", color="white", show.legend = FALSE)+
  # coloring
  scale_fill_viridis(discrete = TRUE)+
  scale_colour_viridis(discrete = TRUE)+
  #labels
  labs(x = xlab, y = ylab)+
  #ggtitle(pca_title)+
  theme_bw()
ggsave(paste0("PCA_", suffix, ".pdf"), plot = last_plot(), device = "pdf", width = 8, height = 8)

# #Second and third Principal components
# xlab = paste("Principal Component 3 (", format(round(percent[3], 1), nsmall=1)," %)", sep="")
# ylab = paste("Principal Component 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")
# ggplot(data = ind_coords, aes(x = Axis3, y = Axis2))+
#   geom_hline(yintercept = 0)+
#   geom_vline(xintercept = 0)+
#   # spider segments
#   geom_segment(aes(xend = Axis3.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
#   # points
#   geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
#   # centroids
#   geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, fontface="bold", color="white", show.legend = FALSE)+
#   # coloring
#   scale_fill_viridis(discrete = TRUE)+
#   scale_colour_viridis(discrete = TRUE)+
#   #labels
#   labs(x = xlab, y = ylab)+
#   ggtitle(pca_title)+
#   theme_bw()
# ggsave(paste0("PCA_2", suffix, ".pdf"), plot = last_plot(), device = "pdf", width = 8, height = 8)
# 
# #First and third principal components
# xlab = paste("Principal Component 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
# ylab = paste("Principal Component 3 (", format(round(percent[3], 1), nsmall=1)," %)", sep="")
# ggplot(data = ind_coords, aes(x = Axis1, y = Axis3))+
#   geom_hline(yintercept = 0)+
#   geom_vline(xintercept = 0)+
#   # spider segments
#   geom_segment(aes(xend = Axis1.cen, yend = Axis3.cen, colour = Site), show.legend = FALSE)+
#   # points
#   geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
#   # centroids
#   geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, fontface="bold", color="white", show.legend = FALSE)+
#   # coloring
#   scale_fill_viridis(discrete = TRUE)+
#   scale_colour_viridis(discrete = TRUE)+
#   #labels
#   labs(x = xlab, y = ylab)+
#   ggtitle(pca_title)+
#   theme_bw()
# ggsave(paste0("PCA_2", suffix, ".pdf"), plot = last_plot(), device = "pdf", width = 8, height = 8)

########################################
# Phylogenetics ########################
########################################

#Methodology based on that of White and Searle (2008)
#Combines Cavalli-Sforza and Edwards distance (Dch) with delta mu distance (Dmyu). 
#Former reliably infers topology but not branch lengths, latter reliably infers branch lengths but not topology
#Uses Fitch to create a hybrid tree that enforces the Ddch topology with the Dmyu branch lengths

#Calculate Dch with hierfstat and convert to distance matrix
mouse_distances_Dch <- genet.dist(mouse_gen, diploid = TRUE, method = "Dch")
distance_matrix_Dch <- as.matrix(mouse_distances_Dch)
colnames(distance_matrix_Dch) <- populations
rownames(distance_matrix_Dch) <- populations

#Generate a neighbour-joining tree from the Dch distance matix
nj_tre <- bionj(distance_matrix_Dch)

#Test that NJ tree is most appropriate for the data
x <- as.vector(mouse_distances_Dch)
y <- as.vector(as.dist(cophenetic(nj_tre)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2
#High correlation value suggests good fit. 

#Can also test UPGMA to see if it's more appropriate. Wasn't for any data set I tested. 
tre_upgma <- as.phylo(hclust(mouse_distances_Dch, method = "average"))
y <- as.vector(as.dist(cophenetic(tre_upgma)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree",
     main="Is UPGMA appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2

#Bootstrapping method adapted from function written by Vojtech Zeisek and Emmanuel Paradis
#Original can be found here: https://www.mail-archive.com/r-sig-phylo@r-project.org/msg03123.html
#Take as input a tree, source data in genind format, and number of permutations.
boot.phylo.nj.pop <- function(origtree, genindata, nperm) {
  
  # Custom variant of function sample to correct its behavior when number of individuals within respective population is 1.
  mysample <- function(x){
    if (length(x) == 1) {out <- x}
    else {out <- sample(x, replace=TRUE)}
    return(out)
  }
  
  # Construction of NJ tree from bootstrapped dataset for further comparison with original tree.
  ## Here it is possible to adjust type and/or parameters of distance matrix and/or tree constructing.
  treemaker <- function (indexnum, listgenindobj) {
    distances <- genet.dist(listgenindobj[[indexnum]], diploid = TRUE, method = "Dch")
    distanceMatrix <- as.matrix(distances)
    genpopobj <- genind2genpop(listgenindobj[[indexnum]])
    populations <- popNames(genpopobj)
    colnames(distanceMatrix) <- populations
    rownames(distanceMatrix) <- populations
    bionj(distanceMatrix)
  }
  
  # Extraction of needed information from original data and preparation of temporary variables.
  genindata2 <- genindata
  data <- genindata@tab
  pop <- genindata@pop
  indiv <- 1:nrow(genindata@tab)
  res <- vector(mode="list", length=nperm)
  
  # In every step (from 1 to number of permutations) sample within respective populations of source data and write it into list of genind objects.
  for (index in 1:nperm) {
    struct <- unlist(tapply(indiv, pop, mysample))
    new.tab <- data[struct,]
    rownames(new.tab) <- make.unique(rownames(new.tab))
    genindata2@tab <- new.tab
    # For correct writing of modified data to the list.
    res[[index]] <- genindata2
  }
  # Conversion of all genind objects in the list to genpop objects.
  #population <- lapply(res, genind2genpop)
  #population <- res
  # Construction of trees and saving of them into the list of trees.
  trees <- vector("list", nperm)
  #for (i in 1:nperm) trees[[i]] <- treemaker(i, res)
  trees <- lapply(1:nperm, treemaker, res)
  # Compare every bootstrapped tree with original one and return vector of bootstrap support values.
  for (i in 1:nperm) storage.mode(trees[[i]]$Nnode) <- "integer"
  storage.mode(origtree$Nnode) <- "integer"
  pp <- prop.part(trees)
  pp <- postprocess.prop.part(pp)
  ans <- prop.clades(origtree, part=pp, rooted=FALSE)
  return(ans)
}

#For large values of nboot (e.g. 1000) this is quite slow and inefficient 
#Use this as an opportunity to grab a cup of tea
myboots <- boot.phylo.nj.pop(nj_tre, mouse_gen, nboots)

#Take a look at this tree
plot.phylo(nj_tre, cex=.8)
#drawSupportOnEdges(myboots/nboots*100,adj = c(0,-0.2),frame="none")
nodelabels(myboots/nboots*100,adj = c(1,-0.2),frame="none")
#title(paste0("NJ Tree - ", pca_suffix))
add.scale.bar()

#Write as a nexus file
write.nexus(nj_tre, file = paste0("nj_tree_", suffix, ".nexus"))

#And export as a pdf
pdf(paste0("NJ_tree_(Dch)_", suffix, ".pdf"), 8, 8)
plot.phylo(nj_tre, cex=.8)
nodelabels(myboots/nboots*100,adj = c(1,-0.2),frame="none")
#title(paste0("NJ Tree (Dch) - ", pca_suffix))
add.scale.bar()
dev.off()

#Now need a distance matrix for Dmyu
#Allele frequencies are standardised in RSTCalc, and Dmyu is calculated from them in the same software
#If you want to look at different populations, you'll have to rerun it through RSTCalc
mouse_standardised_dmyu <- read.csv(dmyu_file, row.names = 1)
mouse_standardised_dmyu <- as.matrix(mouse_standardised_dmyu)

#Use Rfitch to produce a hybrid tree. Uses Minimum Evolution model, and assigns SPT as outgroup
#Important to note it prohibits negative branch lengths! Will instead resolve these as a polytomy. 
#Remember this tree is unrooted - can confirm this with is.rooted(hybrid_tree)
hybrid_tree <- Rfitch(mouse_standardised_dmyu, path = path_to_fitch, method = "ME", tree = nj_tre, outgroup = "STP", negative = FALSE)

#For some reason, Rfitch occasionally errors out with no clear reason
#If the tree does not plot, then just rerun the previous command
plot.phylo(hybrid_tree, cex=.8)

###############################################################
# WARNING! Need to remove bootstraps on polytomies! ###########
###############################################################
#This has to be done manually - don't run the R code blindly! 
#First see which nodes are in a polytomy
nodelabels()

#Then to remove labels, place them in the below command separated by commas
#If there are no polytomies, then can make it empty instead
blanks <- c(9,10, 11)

#Can then use this to hide the corresponding bootstraps
blanks <- blanks - length(hybrid_tree$tip.label)
myboots_blanked <- myboots
myboots_blanked[blanks] <- ""
myboots_blanked <- as.numeric(myboots_blanked)

#Now replot to check it worked
plot.phylo(hybrid_tree, cex=.8)
nodelabels(myboots_blanked/nboots*100,adj = c(1,-0.2),frame="none")
#title(paste0("NJ Hybrid Tree - ", pca_suffix))
add.scale.bar()

#Write as a nexus file
write.nexus(hybrid_tree, file = paste0("hybrid_tree_", suffix, ".nexus"))

#And export as a pdf
pdf(paste0("NJ_hybrid_tree_", suffix, ".pdf"), 8, 8)
plot.phylo(hybrid_tree, cex=.8)
nodelabels(myboots_blanked/nboots*100,adj = c(1,-0.2),frame="none")
#title(paste0("NJ Hybrid Tree - ", pca_suffix))
add.scale.bar()
dev.off()


# robertsonian_microsats
code and data for publication

R code, input files, and outputs for analysing microsatellite data

1) SUMMARY
2) DIRECTORY STRUCTURE
3) DESCRIPTION OF R CODE
4) METHODS
5) CITATIONS

1) SUMMARY
The aim of this script is to take as input a set of microsatellite allele data and a distance matrix and output a series of general analyses.
	All input files are provided.
Analyses include estimated and observed heterozygosity; PCA; pairwise Rst and Fst; and two alternative NJ trees.

Primary input of microsatellite data must be in CSV format. Separate CSV files are given for loci on each chromosome, plus one for all loci together. 
	First column (site) contains the ID for sampling site.
	Second column (ID) contains the ID for the sampled individual
	Third columns and beyond (locus name e.g. C1L317) contain pair of alleles at each locus separated by /
Secondary input is a csv file of a distance matrix (Delta mu) derived from the separate program rstcalc
	If, as in my case, your computer is too modern to run rstcalc, it can be simulated with a program like DosBox. While this input wasn't used in the final version of the paper, the code and data has been included here for completeness.

R script requires user to specify a set of variables prior to use (see lines 27-50)
	Path to the working directory containing input files
	Allele data input filename (csv format)
	A suffix to give output files (users choice)
	A suffix for plot titles (users choice)
	Path to a local installation of the Fitch program in Phylip
	A number of bootstraps to perform for NJ tree support (suggest 1000)
	Filename for an rstcalc standardised distance matrix (Delta mu) output converted to CSV format

One set of loci may be analysed at a time! Script can be left mostly unsupervised but for one small plotting adjustment at the end to remove bootstrap labels at polytomies (see lines 374-405)

2) DIRECTORY STRUCTURE
Note that <abbr> indicates the chromosome (1,3,4,6 or all) on which microsat loci are found.

home						Home directory
│   README.md					README - this file!
│   <abbr>_loci.csv				Input csv file containing microsatellite allele data for <abbr>
│   microsat_methods_final.R			R script to perform analyses of microsat loci
│   mouse_standardised_dmyu_<abbr>.csv		Input csv containing distance matrix (Delta mu) manually edited from rstcalc output
│
└───outputs					Dir containing dirs with results
│   │   file011.txt
│   │   file012.txt
│   │
│   └───<abbr>_loci				Dir of results for each subset of loci identified with <abbr>
│       │   allelic_richness_<abbr>.txt	Data frame of allelic richness estimates from Hierfstat	
│       │   heterozygosity_<abbr>.pdf		Plot of observed and expected heterozygosity at each site
│       │   hybrid_tree_<abbr>.nexus		Nexus file containing hybrid phylogeny
│       │   NJ_hybrid_tree_<abbr>.pdf		Plot of hybrid phylogeny using 
│       │   NJ_tree_(Dch)_<abbr>.pdf		Plot of phylogeny Cavalli-Sforza and Edwards distance 
│       │   nj_tree_all_loci.nexus		Nexus file of phylogeny based on Cavalli-Sforza and Edwards distance 
│       │   PCA_<abbr>_.pdf			Plot of PCA (first two components)
│       │   Rst_Fst_<abbr>.pdf			Plot of pairwise Rst and Fst for each site
│       │   
│       └───genepop_files			Dir containing Genepop files for loci <abbr>
│           │	mouse_genepop_<abbr>.txt	Genepop file converted from genind by graph4lg
│           │	mouse_genepop_<abbr>_fst.txt	Genepop file of pairwise Fst at each locus
│           │	mouse_genepop_<abbr>_fst.txt.MIG	Pairwise Fst between sites
│           │	mouse_genepop_<abbr>_rst.txt	Genepop file of pairwise Rst at each locus
│           │	mouse_genepop_<abbr>_rst.txt.MIG	Pairwise Rst between sites
│   
└───rstcalc_intermediate_files			Files as input and output from rstcalc for all subsets of loci
    │   DELTA-MU-<abbr>				Delta-mu distance matrix output from rstcalc
    │   DELTA-MU-<abbr>-edit			Same distance matrix, but manually edited for clarity
    │   MICROSAT-<abbr>				Microsat allele data in rstcalc format
    │   RSTOUT-<abbr>				Standard rstcalc output file
    │   STND<abbr>				Standardised allele size from rstcalc STANDARD function

3) DESCRIPTION OF R CODE

Lines 1-12: Preamble
	Description of code and warnings re bootstrap plotting (see "Lines 248-405: Phylogenetic analyses.")

Lines 13-26: Libraries
	Required libraries for script to function. Please check that all are correctly installed and loaded! 

Lines 27-52: Variable assignment.
	This is where you need to edit your variables in order for the script to work! 

Lines 53-77: Data preparation.
	Reads in your basic input CSV containing microsat allele data (57-59), converts them to genind and genpop data types using functions in adegenet (62-63), and calculates allelic richness using hierfstat (66-75). 

Lines 78-107: Heterozygosity calculations.
	basic.stats() function in hiefstat calculates estimated and observed heterozygosities for each population (81-93), which are then plotted (95-105).

Lines 108-154: Estimation of Fst and Rst. 
	Genind data converted to genepop format using genind_to_genepop() function in graph4lg package (113). Fst() function in genepop then estimates pairwise Fst and Rst (114-115), an produces several output files. Fst (120-128) and Rst (129-134) values are read from these output files and added to a matrix. Matrix is reformatted (137-141) for easy visualization (144-153). 

Lines 155-247: PCA calculations. 
	Allele counts are extracted from genind (159) and used as input for PCA analysis (162) using dudi.pca() function in ade4 package. A plot indicating the percentage of genetic variance explained by each PC axis is output to console for reference but is not saved (165-167). Centroid coordinates for each population in the PCA are calculated (170-179) and the PCA for components 1 and 2 is visualised (182-203). Note that there is commented code (i.e. code that R will not run) to visualise plots for components 2 against 3 (205-224), and 1 against 3 (226-245). These can be uncommented if desired. 

Lines 248-405: Phylogenetic analyses. 
	Cavalli-Sforza and Edwards distance (Dch) is calculated using the genet.dist() function in hierfstat (257-260). This distance matrix is used to create a neighbour joining (NJ) tree (which is therefore unrooted) using the bionj() function in ape (263). Original pairwise distances are plotted against the pairwise distances on the tree to test that NJ is a suitable method (266-270). The same is done using UPGMA (274-279). These two plots are output to console but not saved. For all provided datasets, NJ was more appropriate than UPGMA. A bootstrapping function is then defined (284-335) and used to perform bootstrap analysis for the NJ tree (339). The NJ tree is then visualised (342-345) and saved as both a nexus (348) and pdf (351-356). 
	As Dch does not always accurately infer branch lengths (see METHODS), a second approach is used to generate a hybrid tree such that Dch is used to estimate tree topology, while Goldstein's Delta Mu Squared distance (Dmu) is used to infer branch lengths. Dmu is calculated from standardised alleles using RSTCALC, such that rather than repeat unit number each allele is represented in standard deviations from the mean. RSTCALC must be run beforehand and the appropriate outputs converted into a CSV file - this has been done for all provided datasets. Dmu is read in to R (361-362) and a hybrid tree is produced with Rfitch() in the Rphylip package (367). Note that Rfitch() may randomly throw an error, in which case this line should be rerun. The hybrid tree is then plotted (371). Rfitch() is also set to forbid negative branch lengths, which occur frequently in these datasets. Negative branch lengths will instead be resolved as a polytomy. 
	Before saving the hybrid tree, we need to check for polytomies so that we don't plot multiple bootstraps on a single node. This final part of the script requires user supervision. Node labels are placed on the plotted hybrid tree (378). The user must then edit the variable "blanks" to include any node labels at a polytomy (382). New labels are produced (385-388), the hybrid tree is replotted (391), and then exported as both a nexus (397) and a pdf (400-405). 

4) METHODS

All analyses were performed in R version 4.1.1 (R Core Team, 2021) except where specified, across both all loci and the loci on each individual chromosome. Visualisations make use of the ggplot2 (Wickham, 2016) and viridis (Garnier et al, 2021) packages. 

Alleles for each microsatellite loci were converted to genind and genpop format using the adegenet package (Jombart, 2008; Jombart & Ahmed, 2011). Estimates of allelic richness, as well as the expected and observed heterozygosities of each populaltion, were calculated using the allelic.richness() and basic.stats() functions, respectively, in hierfstat (Goudet & Jombart, 2021). Principle components analysis (PCA) was performed on the allele counts across individual mice using the ade4 (Dray & Dufour, 2007) package's dudi.pca() function. After data was transformed from genind to genepop format using the genind_to_genepop() function in graph4lg (Savary et al., 2020), pairwise Fst (identity based) and Rst (allele-size based) between each population was calculated using the Fst() function in Genepop (Rousset, 2008). 

Phylogenetic relationships between each population were estimated using a neighbour joining (NJ) method based on two sets of distance matrices derived from microsatellite alleles. The first distance method is the chord measure (Dch) of Cavalli-Sforza and Edwards (1967), which was estimated with the genet.dist() function in hierfstat (Goudet & Jombart, 2021). Dch distances were used to estimate an NJ tree using Gascuel's (1997) BIONJ algorithm in the R package ape version 5.6.1 (Paradis and Schliep, 2019). Bootstrapping (1000 replicates) was performed to gauge node stability using a custom function adapted from Zeisek & Paradis (2014).

While Dch is based on the infinite alleles model (IAM) rather than the stepwise mutation model (SMM) of evolution, there is evidence that it outperforms SMM estimates when performing phylogenetic reconstruction between closely related and recently diverged taxa (Goldstein and Pollock, 1997). However, conflicting studies find that Dch will reliably infer tree topology but not branch lengths (Takezaki & Nei, 1996; Angers & Bernatchez, 1998). In contrast, these same analyses find that Goldstein et al's (1995) delta mu squared (Dmu) will accurately infer branch lengths but not topology. 

As such, one may also employ the hybrid tree approach of Angers & Bernatchez (1998) and White & Searle (2008). We include the code and description of such an approach, but ultimately did not include this in the final manuscript as the phylogenies were poorly resolved. Pairwise Dmu was calculated from standardised allele sizes using RSTCALC (Goodman, 1997). A hybrid tree was constructed using the RFitch() function in RPhylip (Felsenstein, 2013; Revell & Chamberlain, 2014) to fix the topology estimated from Dch and estimate branch lengths with Dmu. Negative branch lengths were prohibited, and are instead resolved as polytomies. Both the Dch only trees and the hybrid trees were visualised with functions in ape (Paradis & Schliep, 2019). As they were reconstructed with NJ methods, all phylogenies are unrooted.

5) CITATIONS

Angers B. and Bernatchez L. (1998). Combined use of SMM and non-SMM methods to infer fine structure and evolutionary history of closely related brook charr (Salvelinus fontinalis, Salmonidae) populations from microsatellites. Molecular Biology and Evolution. 15: 143–159.

Cavalli-Sforza L.L., and Edwards A.W.F. (1967). Phylogenetic analysis: models and estimation procedures. Evolution. 32: 550-570. [also published in American Journal of Human Genetics. 19: 233-257 (1967)].

Dray S., Dufour A. (2007). The ade4 Package: Implementing the Duality Diagram for Ecologists. Journal of Statistical
Software. 22(4): 1-20. doi: 10.18637/jss.v022.i04 (URL: https://doi.org/10.18637/jss.v022.i04).

Felsenstein J. (2013). PHYLIP (Phylogeny Inference Package) version 3.695. Distributed by the author. Department of Genome
Sciences, University of Washington, Seattle.

Garnier S., Ross N., Rudis R., Camargo A. P., Sciaini M., and Scherer C. (2021). Rvision - Colorblind-Friendly Color Maps for R. R package version 0.6.2.

Gascuel, O. (1997). BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. Molecular Biology and Evolution. 14: 685–695.

Goldstein D.B., Ruiz Linares A., Cavalli-Sforza L.L., Feldman M.W. (1995). Genetic absolute dating based on microsatellites and the origin of modern humans. Proceedings of the National Academy of Sciences of the USA. 92: 6723–6727

Goldstein D.B. and Pollock D.D. 1997. Launching microsatellites: a review of mutation processes and methods of phylogenetic interference. Journal of Heredity. 88:335-42.

Goodman S.J. (1997). RST CALC: a collection of computer programs for calculating unbiased estimates of genetic differentiation and determining their significance for microsatellite data. Molecular Ecology. 6: 881–885.

Goudet J. and Jombart T. (2021). hierfstat: Estimation and Tests of Hierarchical F-Statistics. R package version 0.5-10. https://CRAN.R-project.org/package=hierfstat

Jombart T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics 24: 1403-1405. doi: 10.1093/bioinformatics/btn129

Jombart T. and Ahmed I. (2011). adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. 27(21): 3070-3071. doi:10.1093/bioinformatics/btr521

Paradis E. and Schliep K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics. 35: 526-528

R Core Team. (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/

Revell, L.J. and Chamberlain, S.A. (2014). Rphylip: An R interface for PHYLIP. Methods in Ecology and Evolution. 5: 976–981. doi: 10.1111/2041-210X.12233

Rousset F. (2008) genepop’007: a complete re-implementation of the genepop software for Windows and Linux. Molecular Ecology Resources. 8: 103-106 http://dx.doi.org/10.1111/j.1471-8286.2007.01931.x

Savary P., Foltête J.C., Moal H., Vuidel G. and Garnier S. (2020). graph4lg: a package for constructing and analysing graphs for landscape genetics in R. Methods in Ecology and Evolution. 00(1): 1-9. URL: https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13530

Takezaki N. and Nei M. (1996). Genetic distances and reconstruction of phylogenetic trees from microsatellite DNA. Genetics. 144: 389–399

White T.A. and Searle J.B. (2008). The colonization of Scottish islands by the common shrew, Sorex araneus (Eulipotyphla: Soricidae). Biological Journal of the Linnean Society. 94: 797–808

Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, New York.

Zeisek V. and Paradis E. (2014). Re: [R-sig-phylo] Problems with bootstraps of NJ tree from SSRs data. The Mail Archive. https://www.mail-archive.com/r-sig-phylo@r-project.org/msg03123.html Last accessed: 7th March 2022.

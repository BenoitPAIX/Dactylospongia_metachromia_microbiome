
####  Project name : Dactylospongia_microbiome ####

####  Description of the project ####
# The sponge Dactylospongia metachromia is known to produce natural products of pharmacological interest, and a ongoing aquaculture program is aiming to grow and cultivate this sponge in French Polynesia
# To study the best conditions of aquaculture of the sponge during the farming, it is essential to consider it as an holobiont system ...
# ... considering that the microbiome could play an important role for the sponge growth, health, resilience to environmental disturbances, and production of natural products of interest
# To this end, the study aims to investigate the spatio-temporal dynamics and stability of the prokaryotic community of the sponge model Dactylospongia metachromia. 

####  Objective of the script ####
# To analyse the 16S rRNA gene metabarcoding reads obtained from the DADA2 pipeline. 
# The main package used for this analysis will be phyloseq (for more details see tutorials such as https://joey711.github.io/phyloseq/index.html)

#### Author of the script : Dr. Benoît PAIX

#### Citation : not yet available (paper submitted)




####  1. Setting the working directory ####

#use setwd, example : setwd("C:/Users/benoit_paix/Google Drive/2021_2023 Postdoc/3. Dactylospongia microbiome project/4. Lab notebook/3. Metabarcoding analyses/3. Phyloseq analysis)

####  2. Installation and loading of the packages ####


# example of few package installations : 
#install.packages("randomcoloR")
#install.packages("ggplot2")
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
#install.packages("stringi")
#install.packages("Rcpp")
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install(version='devel')
#BiocManager::install("microbiome")
#install.packages("hrbrthemes")
#install.packages("ggthemes")
#install.packages("RColorBrewer")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Loading of the packages
library(stringi)
library(vegan)
library(Rcpp)
library(ggplot2)
library(randomcoloR)
library(phyloseq)
library(ape)
library(dplyr)
library(agricolae)
library(RVAideMemoire)
library(microbiome)
library(hrbrthemes)
library(cowplot)
library(RColorBrewer)
library(fossil)
library(reshape2)
library(tidyverse)
library(multcompView)
library(rcompanion)
library(pairwiseAdonis) 
library(ggpubr)
library(scales)
library(microbiomeMarker)


####  3. Reading of the csv files (ASV, Taxonomy and Metadata tables)   ####

ASV_table = read.csv(file = "ASV_table_filtered.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(ASV_table)

TAX_table = read.csv(file = "Taxonomy_table_filtered.csv" , sep = ";" , header = T , row.names = 1)
dim(TAX_table)

META_table = read.csv(file = "Metadata_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(META_table)


# to have a taxonomy table with all taxonomical ranks names concatenated
TAX_table = as.data.frame(TAX_table)
TAX_table

TAX_table.2 <- TAX_table

TAX_table.2$Kingdom <- TAX_table.2$Kingdom
TAX_table.2$Phylum <- paste(TAX_table.2$Kingdom,"|" ,TAX_table$Phylum)
TAX_table.2$Class <- paste(TAX_table.2$Phylum,"|" ,TAX_table$Class)
TAX_table.2$Order <- paste(TAX_table.2$Class ,"|" ,TAX_table$Order)
TAX_table.2$Family <- paste(TAX_table.2$Order,"|" ,TAX_table$Family)
TAX_table.2$Genus <- paste(TAX_table.2$Family,"|" ,TAX_table$Genus)

head(TAX_table.2)
dim(TAX_table.2)


####  4. Creation of the main phyloseq object  ####

ASV_phylo = otu_table(ASV_table, taxa_are_rows = TRUE)
dim(ASV_phylo)

TAX_table.2 = as.matrix(TAX_table.2)
TAX_phylo = tax_table(TAX_table.2)
TAX_phylo
dim(TAX_phylo)

META_phylo = sample_data(META_table)
dim(META_phylo)

physeq = phyloseq(ASV_phylo, TAX_phylo, META_phylo)
physeq


####  5. Verification of the rarefaction curves ####


#to remove negative controle and extraction blank from the analysis
physeq_rarecurve = subset_samples(physeq, Sample_type == "True Sample")
physeq_rarecurve


data_rarecurve = rarecurve(t(otu_table(physeq_rarecurve)), step=50,  tidy = TRUE)
data_rarecurve

plot_rarecurve = ggplot(data_rarecurve, aes(x = Sample   , y = Species, colour = Site))
plot_rarecurve = plot_rarecurve + geom_line(size = 1, alpha = 0.7)
plot_rarecurve = plot_rarecurve + theme_bw(15)
plot_rarecurve = plot_rarecurve + scale_color_manual(values = rep(c("black"),times=52))
plot_rarecurve = plot_rarecurve +  scale_x_continuous(labels = comma_format(big.mark = ".",
                                                                            decimal.mark = ","))
plot_rarecurve = plot_rarecurve + theme(legend.position = "none")
plot_rarecurve = plot_rarecurve + xlab("Sample size (reads number)") + ylab("ASVs number")
plot_rarecurve # Figure saved as Figure S5




####  6.Subsampling for the BioGeographical Study (BGS) ####

physeq_bgs = subset_samples(physeq, Study == "Biogeography")
physeq_bgs


####  7. BGS - Creation of 2 phyloseq objects  ####

# one phyloseq object with the rarefied dataset for alpha div analyses

min(sample_sums(physeq_bgs))

physeq_bgs_rarefied = rarefy_even_depth(physeq_bgs)
physeq_bgs_rarefied

sample_sums(physeq_bgs_rarefied)

# the other one with the compositional dataset (percentages) for compositional analyses

physeq_bgs_compo = transform(physeq_bgs, "compositional")
physeq_bgs_compo


#### 8. BGS - Defining a color palette for the different sites of the BGS ####


colorpal_site = c("01.Red Sea" = "brown1" ,
                  "02.Maldives" = "darkgreen" ,
                  "03.Indonesia" = "darkorange" ,
                  "04.FP_Tetiaroa" = "blueviolet" ,
                  "05.FP_Rangiroa_Avatoru" = "cornflowerblue" ,
                  "06.FP_Rangiroa_Lagon2" = "darkslategray4" ,
                  "07.FP_Rangiroa_Lagon1" = "cyan" ,
                  "08.FP_Rangiroa_Tiputa" = "darkblue" ,
                  "09.FP_Makemo" = "bisque2" ,
                  "10.FP_Raroia" = "lightgreen" ,
                  "11.FP_Tematangi" = "#FFFF33" ,
                  "12.FP_Mangareva" = "#F781BF" )

#### 9. BGS - Analysis of the alpha-diversity : boxplots ####

data_alpha_bgs = estimate_richness(physeq_bgs_rarefied , measures = c("Observed", "Shannon", "Chao1"))
data_alpha_bgs

Pielou = data_alpha_bgs$Shannon / log(data_alpha_bgs$Observed)
Pielou

data_alpha_bgs = cbind(sample_data(physeq_bgs_rarefied), data_alpha_bgs , Pielou)
data_alpha_bgs = data_alpha_bgs[, c("Region_site", "Shannon" ,"Chao1","Pielou")]
data_alpha_bgs

data_alpha_bgs_long = melt(data_alpha_bgs, id.var=c("Region_site"))
data_alpha_bgs_long

plot_alpha_bgs = ggplot(data_alpha_bgs_long, aes(Region_site ,value, fill = Region_site))
plot_alpha_bgs = plot_alpha_bgs + geom_boxplot(alpha = 0.8, size = 1) + facet_grid(  variable ~ ., scales="free") 
plot_alpha_bgs = plot_alpha_bgs + geom_point(size = 2, alpha = 0.8, pch = 21, stroke = 1)
plot_alpha_bgs = plot_alpha_bgs + theme_bw(base_size = 15) 
plot_alpha_bgs = plot_alpha_bgs + theme(legend.position="left")
plot_alpha_bgs = plot_alpha_bgs + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_alpha_bgs = plot_alpha_bgs + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
plot_alpha_bgs = plot_alpha_bgs + scale_fill_manual(values = colorpal_site)
plot_alpha_bgs # Figure saved as Figure 2A


#### 10. BGS - Analysis of the alpha-diversity : statistical tests  ####

shapiro.test(data_alpha_bgs$Shannon)
#From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
#ANOVA test 

shapiro.test(data_alpha_bgs$Chao1)
#From the output, the p-value < 0.05 implying that the distribution of the data are significantly different from normal distribution. In other words, we cannot assume the normality.
#Kruskal-Wallis test 


shapiro.test(data_alpha_bgs$Pielou)
#From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
#ANOVA test 

data_anova_shannon_bgs <- aov(Shannon ~ Region_site, data_alpha_bgs)
data_anova_shannon_bgs
summary(data_anova_shannon_bgs)

data_hsd_shannon_bgs = HSD.test(aov(Shannon ~ Region_site, data_alpha_bgs), "Region_site", group=T)
data_hsd_shannon_bgs


data_kruskal_Chao1_bgs = kruskal.test(Chao1 ~ Region_site, data_alpha_bgs)
data_kruskal_Chao1_bgs


data_wilcox_Chao1_bgs = pairwise.wilcox.test(data_alpha_bgs$Chao1, data_alpha_bgs$Region_site, p.adjust.method = "none")
data_wilcox_Chao1_bgs
data_wilcox_Chao1_bgs$p.value

data_wilcox_Chao1_bgs_full = fullPTable(data_wilcox_Chao1_bgs$p.value)
data_wilcox_Chao1_bgs_full

multcompLetters(data_wilcox_Chao1_bgs_full,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)


data_anova_pielou_bgs <- aov(Pielou ~ Region_site, data_alpha_bgs)
data_anova_pielou_bgs
summary(data_anova_pielou_bgs)

data_hsd_pielou_bgs = HSD.test(aov(Pielou ~ Region_site, data_alpha_bgs), "Region_site", group=T)
data_hsd_pielou_bgs
# Statistical data saved in Table S2 (SI)


#### 11. BGS - BETA DIVERSITY : NMDS analysis ####


nmds_bgs = ordinate(physeq_bgs_compo, "NMDS", "bray")
nmds_bgs$points
nmds_bgs$stress

data_nmds_bgs = plot_ordination(physeq_bgs_compo, nmds_bgs, type="samples", justDF = TRUE)
data_nmds_bgs

plot_nmds_bgs = ggplot(data_nmds_bgs, aes(NMDS1, NMDS2))
plot_nmds_bgs = plot_nmds_bgs + geom_point(aes(fill = Region_site), pch = 21, stroke = 1.5 , size = 7, alpha = 0.9, color = "black" , show.legend = T)
plot_nmds_bgs = plot_nmds_bgs + theme_bw(base_size = 17) 
plot_nmds_bgs = plot_nmds_bgs + theme(legend.position="left")
plot_nmds_bgs = plot_nmds_bgs + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_nmds_bgs = plot_nmds_bgs + scale_fill_manual(values = colorpal_site)
plot_nmds_bgs = plot_nmds_bgs + annotate("text", label="2D stress = 0.117 ", x=1.2, y=0.9, size = 6)
plot_nmds_bgs

legend_bgs = as_ggplot(get_legend(plot_nmds_bgs))
legend_bgs

plot_nmds_bgs = plot_nmds_bgs + theme(legend.position="none")
plot_nmds_bgs  # Figure saved as Figure 3A


#### 12. BGS - BETA DIVERSITY : PERMANOVA ####

metadata_bgs <- as(sample_data(physeq_bgs), "data.frame")
metadata_bgs

data_distbeta_bgs = as.matrix(distance(physeq_bgs_compo, method="bray"))
data_distbeta_bgs

data_permanova_bgs = adonis2(data_distbeta_bgs ~ Region_site, data = metadata_bgs)
data_permanova_bgs # Statistical data saved in Table S3A (SI)

data_pairwiseadonis_bgs = pairwise.adonis(data_distbeta_bgs, metadata_bgs$Region_site)
data_pairwiseadonis_bgs



#### 13. BGS - Distance based RDA analysis using environmental parameters ####

physeq_bgsenv = subset_samples(physeq, Study == "Biogeography" & Env_analysis== "Yes")
physeq_bgsenv

physeq_bgsenv_compo = transform(physeq_bgsenv, "compositional")
physeq_bgsenv_compo

data_distbeta_bgsenv = as.matrix(distance(physeq_bgsenv_compo, method="bray"))
data_distbeta_bgsenv


metadata_bgsenv = as(sample_data(physeq_bgsenv), "data.frame")
metadata_bgsenv

dbRDA_bgsenv <- capscale(data_distbeta_bgsenv ~ Temperature +	Salinity  + PO4 + NOx  +	SiOHx , metadata_bgsenv)
dbRDA_bgsenv

anova.cca(dbRDA_bgsenv, step=999) # Statistical data saved in Table S3B (SI)


data_dbRDA_bgsenv_samples = scores(dbRDA_bgsenv , display = "sites")
data_dbRDA_bgsenv_samples = as.data.frame(data_dbRDA_bgsenv_samples)
data_dbRDA_bgsenv_samples = cbind(data_dbRDA_bgsenv_samples, sample_data(physeq_bgsenv))
head(data_dbRDA_bgsenv_samples)

data_dbRDA_bgsenv_arrow = scores(dbRDA_bgsenv , display = "bp")
data_dbRDA_bgsenv_arrow = as.data.frame(data_dbRDA_bgsenv_arrow)
head(data_dbRDA_bgsenv_arrow)


plot_dbRDA_bgsenv = ggplot(data_dbRDA_bgsenv_samples , aes(CAP1, CAP2))
plot_dbRDA_bgsenv = plot_dbRDA_bgsenv + geom_segment(color = "black", data = data_dbRDA_bgsenv_arrow , mapping=aes(x=0, y=0, xend=CAP1, yend=CAP2), arrow=arrow(length=unit(0.30,"cm"), type = "open" ), size=1)
plot_dbRDA_bgsenv = plot_dbRDA_bgsenv + geom_point( aes(fill = Region_site) , size = 7, stroke = 1.5, alpha = 0.9, pch = 21  , color = "black", show.legend = T)
plot_dbRDA_bgsenv = plot_dbRDA_bgsenv + scale_fill_manual(values = colorpal_site)
plot_dbRDA_bgsenv = plot_dbRDA_bgsenv + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_dbRDA_bgsenv = plot_dbRDA_bgsenv + geom_text(data = data_dbRDA_bgsenv_arrow ,cex = 6, color = "black" ,mapping = aes(CAP1 , CAP2) , label = row.names(data_dbRDA_bgsenv_arrow))
plot_dbRDA_bgsenv = plot_dbRDA_bgsenv + theme_bw(base_size = 17) 
plot_dbRDA_bgsenv = plot_dbRDA_bgsenv + theme(legend.position="none")
plot_dbRDA_bgsenv # Figure saved as Figure 3B



### 14. Variance paritioning analysis using environmental parameters (only) ####

head(data_distbeta_bgsenv)
dim(data_distbeta_bgsenv)

varpart_bgsenv = varpart(data_distbeta_bgsenv , ~ Temperature ,	~Salinity,~ PO4 + NOx  +	SiOHx, data = metadata_bgsenv)
varpart_bgsenv

showvarparts(3, bg=2:4)
plot(varpart_bgsenv, bg=2:4)
# Figure saved as Figure S6 (SI)


#### 15. Mantel test analysis with geographical distances ####

data_gps_bgs = as(sample_data(physeq_bgs), "data.frame")[c("longitude","latitude")]
data_gps_bgs

data_distgps_bgs = earth.dist(data_gps_bgs)
data_distgps_bgs = as.matrix(data_distgps_bgs)
data_distgps_bgs
row.names(data_distgps_bgs) <- row.names(data_gps_bgs)
colnames(data_distgps_bgs) <- row.names(data_gps_bgs) 

head(data_distgps_bgs)
dim(data_distgps_bgs)


data_mantel_gpsbeta = mantel(data_distgps_bgs, data_distbeta_bgs, method = "spearman", permutations = 9999, na.rm = TRUE)
data_mantel_gpsbeta


#### 16. Variance paritioning analysis using environmental parameters and geographical distances ####

data_gps_bgsenv = as(sample_data(physeq_bgsenv), "data.frame")[c("longitude","latitude")]
data_gps_bgsenv

data_distgps_bgsenv = earth.dist(data_gps_bgsenv)
data_distgps_bgsenv = as.matrix(data_distgps_bgsenv)
data_distgps_bgsenv
row.names(data_distgps_bgsenv) <- row.names(data_gps_bgsenv)
colnames(data_distgps_bgsenv) <- row.names(data_gps_bgsenv) 

head(data_distgps_bgsenv)
dim(data_distgps_bgsenv)

data_distgps_bgsenv_pcnm = pcnm(data_distgps_bgsenv)
data_distgps_bgsenv_pcnm$vectors



varpart_bgsenv_gps = varpart(data_distbeta_bgsenv , ~ Temperature +	Salinity +  PO4 + NOx  +	SiOHx,  ~ data_distgps_bgsenv_pcnm$vectors, data = metadata_bgsenv)
varpart_bgsenv_gps

showvarparts(2, bg = c("purple","yellow"))
plot(varpart_bgsenv_gps, bg = c("purple","yellow"))
# Figure saved as Figure S6 (SI)



#### 17. MANTEL test analysis with phylogenetic distances based on 28S marker ####


physeq_bgsphylo = subset_samples(physeq,  Phylogeny_analysis == "Yes")
physeq_bgsphylo


physeq_bgsphylo_compo = transform(physeq_bgsphylo, "compositional")
physeq_bgsphylo_compo

data_distbeta_bgsphylo = distance(physeq_bgsphylo_compo, method="bray") 
data_distbeta_bgsphylo = as.matrix(dist(data_distbeta_bgsphylo))
data_distbeta_bgsphylo

head(data_distbeta_bgsphylo)
dim(data_distbeta_bgsphylo)


data_distphylo = read.csv(file = "Phylogeny_dist_matrix.csv" , header = TRUE , sep = ";" , row.names = 1)
data_distphylo

data_mantel_phylobeta = mantel(data_distbeta_bgsphylo, data_distphylo, method = "spearman", permutations = 9999, na.rm = TRUE)
data_mantel_phylobeta


#### 18 BGS - Barplot analysis  of the totalcommunity at the family level ####


physeq_bgs_compo_family = aggregate_rare(physeq_bgs_compo, level = "Family", detection = 3/100, prevalence = 0/100)
physeq_bgs_compo_family

ntaxa(physeq_bgs_compo_family)

colorpal_family = c("Archaea | Crenarchaeota | Nitrososphaeria | Nitrosopumilales | Nitrosopumilaceae"="darkorchid",
                      "Bacteria | Acidobacteriota | Acidobacteriae | PAUC26f | NA"="darkorange",
                      "Bacteria | Acidobacteriota | Thermoanaerobaculia | Thermoanaerobaculales | Thermoanaerobaculaceae"="darkgoldenrod1",
                      "Bacteria | Acidobacteriota | Vicinamibacteria | Subgroup 9 | NA"="goldenrod3",
                      "Bacteria | Acidobacteriota | Vicinamibacteria | Vicinamibacterales | NA"="gold",
                      "Bacteria | Actinobacteriota | Acidimicrobiia | Microtrichales | Microtrichaceae"="darkslateblue",
                      "Bacteria | Bacteroidota | Rhodothermia | Rhodothermales | Rhodothermaceae"="red",
                      "Bacteria | Chloroflexi | Anaerolineae | Caldilineales | Caldilineaceae"="palegreen4",
                      "Bacteria | Chloroflexi | Anaerolineae | SBR1031 | A4b"="olivedrab1",
                      "Bacteria | Chloroflexi | Dehalococcoidia | SAR202 clade | NA"="seagreen4",
                      "Bacteria | Chloroflexi | TK10 | NA | NA"="green1",
                      "Bacteria | Dadabacteria | Dadabacteriia | Dadabacteriales | NA"="plum1",
                      "Bacteria | Entotheonellaeota | Entotheonellia | Entotheonellales | Entotheonellaceae"="slategray4",
                      "Bacteria | Gemmatimonadota | BD2-11 terrestrial group | NA | NA"="salmon",
                      "Bacteria | Myxococcota | bacteriap25 | NA | NA"="slategray2",
                      "Bacteria | Nitrospirota | Nitrospiria | Nitrospirales | Nitrospiraceae"="aquamarine4",
                      "Bacteria | PAUC34f | NA | NA | NA"="yellow1",
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Defluviicoccales | NA"="mediumpurple1",
                      "Bacteria | Proteobacteria | Alphaproteobacteria | Rhodobacterales | Rhodobacteraceae"="mediumpurple4",
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Nitrosococcales | Nitrosococcaceae"="cyan",
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | Endozoicomonadaceae"="mediumblue",
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | Halieaceae"="dodgerblue",
                      "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | KI89A clade"="dodgerblue4",
                      "Bacteria | Proteobacteria | Gammaproteobacteria | UBA10353 marine group | NA"="cadetblue2",
                      "Bacteria | Spirochaetota | Spirochaetia | Spirochaetales | Spirochaetaceae"="peachpuff",
                      "Bacteria | Verrucomicrobiota | Verrucomicrobiae | Opitutales | Puniceicoccaceae"="sienna",
                      "Other"="gray")



# for a barplot without border the phyloseq function plot_bar need to be changed. 
#Use the plot_bar_2 function as follows

plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


barplot_bgs_family = plot_bar_2(physeq_bgs_compo_family ,"Replicate_name" ,fill = "Family") + theme_bw(base_line_size = 12)
barplot_bgs_family = barplot_bgs_family + geom_bar(stat = "identity"  , position="stack", size = 3)
barplot_bgs_family = barplot_bgs_family + scale_fill_manual(values = colorpal_family)
barplot_bgs_family = barplot_bgs_family + theme(legend.text = element_text(size=9),
                          axis.ticks.y = element_blank(),
                          legend.position="none",
                          axis.text.x = element_blank(),
                          legend.title = element_blank(),
                          axis.ticks.x=element_blank(), 
                          axis.title = element_blank())
barplot_bgs_family = barplot_bgs_family + facet_grid(cols = vars(Region_site), scales = "free", space = "free")
barplot_bgs_family = barplot_bgs_family + guides(fill = guide_legend(ncol = 1)) 
barplot_bgs_family = barplot_bgs_family + scale_y_percent()
# barplot_bgs_family = barplot_bgs_family + coord_flip()
#barplot_bgs_family = barplot_bgs_family + theme_bw()
barplot_bgs_family # Figure saved as Figure 4A


# 19. BGS - LEfSe analysis to identify site-specific biomarker taxa ####


# Sites with n= 1 are removed from the analysis (Makemo and Maldives)

physeq_bgslefse = subset_samples(physeq, Study == "Biogeography" 
                                         & Region_location != "06.FP_Makemo" 
                                         & Region_location != "02.Maldives_Faafu.Atoll")
physeq_bgslefse

physeq_bgslefse_compo = transform(physeq_bgslefse, "compositional")
physeq_bgslefse_compo

lefse_bgs <-run_lefse(physeq_bgslefse ,
                      wilcoxon_cutoff = 0.01,
                      norm = "CPM",
                      group = "Region_location",
                      kw_cutoff = 0.01,
                      multigrp_strat = TRUE,
                      lda_cutoff = 4)
lefse_bgs

data_lefse_bgs = marker_table(lefse_bgs)
data_lefse_bgs

plot_cladogram(lefse_bgs, color = c("brown1","darkorange","blueviolet", "lightgreen","#FFFF33", "#F781BF")) +
  theme(plot.margin = margin(0, 0, 0, 0)) # Figure saved as Figure S7 (SI)


#### 20.Subsampling for the Farming Trial Study (FTS) ####

physeq_fts = subset_samples(physeq, Study == "Aquaculture" | Region_site == "05.FP_Rangiroa_Avatoru")
physeq_fts


#### 21. FTS - Creation of 2 phyloseq objects  ####

# one phyloseq object with the rarefied dataset for alpha div analyses

min(sample_sums(physeq_fts))

physeq_fts_rarefied = rarefy_even_depth(physeq_fts)
physeq_fts_rarefied

sample_sums(physeq_fts_rarefied)

# the other one with the compositional dataset (percentages) for compositional analyses

physeq_fts_compo = transform(physeq_fts, "compositional")
physeq_fts_compo

#### 22. FTS - Analysis of the alpha-diversity : boxplots ####

data_alpha_fts= estimate_richness(physeq_fts_rarefied , measures = c("Observed", "Shannon", "Chao1"))
data_alpha_fts

Pielou = data_alpha_fts$Shannon / log(data_alpha_fts$Observed)
Pielou

data_alpha_fts= cbind(sample_data(physeq_fts_rarefied),data_alpha_fts,Pielou)
data_alpha_fts

data_alpha_fts = data_alpha_fts[, c("Culture_time", "Shannon" ,"Chao1","Pielou")]
data_alpha_fts

data_alpha_fts_long = melt(data_alpha_fts, id.var=c("Culture_time"))
data_alpha_fts_long

plot_alpha_fts = ggplot(data_alpha_fts_long, aes(Culture_time ,value))
plot_alpha_fts = plot_alpha_fts + geom_boxplot(alpha = 0.7, size = 1,  fill = "cornflowerblue") + facet_grid(  variable ~ ., scales="free") 
plot_alpha_fts = plot_alpha_fts + geom_point(aes(pch = Culture_time), fill = "cornflowerblue" , size = 3, alpha = 0.8, stroke = 1)
plot_alpha_fts = plot_alpha_fts + theme_bw(base_size = 15) 
plot_alpha_fts = plot_alpha_fts + theme(legend.position="right")
plot_alpha_fts = plot_alpha_fts + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
# plot_alpha_fts = plot_alpha_fts + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
plot_alpha_fts = plot_alpha_fts + scale_shape_manual(values = c(21, 23, 22,24, 25) )
plot_alpha_fts # Figure saved as Figure 2B


#### 23. FTS - Analysis of the alpha-diversity : statistical tests ####

shapiro.test(data_alpha_fts$Shannon)
#From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
#ANOVA test

shapiro.test(data_alpha_fts$Chao1)
#From the output, the p-value < 0.05 implying that the distribution of the data are significantly different from normal distribution. In other words, we cannot assume the normality.
#Kruskal-Wallis test

shapiro.test(data_alpha_fts$Pielou)
#From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
#ANOVA test



data_anova_shannon_fts <- aov(Shannon ~ Culture_time, data_alpha_fts)
data_anova_shannon_fts
summary(data_anova_shannon_fts)

data_hsd_shannon_fts = HSD.test(aov(Shannon ~ Culture_time, data_alpha_fts), "Culture_time", group=T)
data_hsd_shannon_fts


data_kruskal_Chao1_fts = kruskal.test(Chao1 ~ Culture_time, data_alpha_fts)
data_kruskal_Chao1_fts


data_wilcox_Chao1_fts = pairwise.wilcox.test(data_alpha_fts$Chao1, data_alpha_fts$Culture_time,
                                       p.adjust.method = "none")
data_wilcox_Chao1_fts
data_wilcox_Chao1_fts$p.value

data_wilcox_Chao1_fts_full = fullPTable(data_wilcox_Chao1_fts$p.value)
data_wilcox_Chao1_fts_full

multcompLetters(data_wilcox_Chao1_fts_full,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)



data_anova_pielou_fts <- aov(Pielou ~ Culture_time, data_alpha_fts)
data_anova_pielou_fts
summary(data_anova_pielou_fts)

data_hsd_pielou_fts = HSD.test(aov(Pielou ~ Culture_time, alpha.data.meta), "Culture_time", group=T)
data_hsd_pielou_fts
# Statistical data saved in Table S2 (SI)



#### 24. FTS - BETA DIVERSITY : NMDS analysis ####


nmds_fts = ordinate(physeq_fts_compo, "NMDS", "bray")
nmds_fts$points
nmds_fts$stress

data_nmds_fts = plot_ordination(physeq_fts_compo, nmds_fts, type="samples", justDF = TRUE)
data_nmds_fts

plot_nmds_fts = ggplot(data_nmds_fts, aes(NMDS1, NMDS2))
plot_nmds_fts = plot_nmds_fts + geom_point(aes(pch = Culture_time), fill = "cornflowerblue",  stroke = 2 , size = 7, alpha = 0.9, color = "black" , show.legend = T)
plot_nmds_fts = plot_nmds_fts + theme_bw(base_size = 17) 
plot_nmds_fts = plot_nmds_fts + theme(legend.position="left")
plot_nmds_fts = plot_nmds_fts + guides(fill = guide_legend(override.aes = list(shape = 21)))
plot_nmds_fts = plot_nmds_fts + scale_shape_manual(values = c(21, 23, 22,24, 25))
plot_nmds_fts = plot_nmds_fts + annotate("text", label="2D stress = 0.177", x=0.25, y=0.35, size = 6)
plot_nmds_fts = plot_nmds_fts + theme(legend.position="none")
plot_nmds_fts # Figure saved as figure 3C


#### 25. FTS - BETA DIVERSITY : PERMANOVA ####

metadata_fts <- as(sample_data(physeq_fts_compo), "data.frame")
metadata_fts

data_distbeta_fts = as.matrix(distance(physeq_fts_compo, method="bray"))
data_distbeta_fts


data_permnova_fts = adonis2(data_distbeta_fts ~ Culture_time, data = metadata_fts)
data_permnova_fts # Statistical data saved in Table S3A (SI)

data_pairwiseadonis_fts = pairwise.adonis(data_distbeta_fts, metadata_fts$Culture_time)
data_pairwiseadonis_fts



#### 26 FTS - Distance based RDA analysis using environmental parameters ####

data_distbeta_fts

metadata_fts

dbRDA_fts <- capscale(data_distbeta_fts ~ Temperature +	Salinity  + Chlorophyll_a + PO4 + NOx  +	SiOHx, metadata_fts)
dbRDA_fts

anova.cca(dbRDA_fts, step=999)# Statistical data saved in Table S3B (SI)

data_dbRDA_fts_samples = scores(dbRDA_fts , display = "sites")
data_dbRDA_fts_samples = as.data.frame(data_dbRDA_fts_samples)
data_dbRDA_fts_samples = cbind(data_dbRDA_fts_samples, metadata_fts)
head(data_dbRDA_fts_samples)

data_dbRDA_fts_arrow = scores(dbRDA_fts , display = "bp")
data_dbRDA_fts_arrow = as.data.frame(data_dbRDA_fts_arrow)
data_dbRDA_fts_arrow

plot_dbRDA_fts = ggplot(data_dbRDA_fts_samples , aes(CAP1, CAP2))
plot_dbRDA_fts = plot_dbRDA_fts + geom_segment(color = "black", data = data_dbRDA_fts_arrow , mapping=aes(x=0, y=0, xend=CAP1, yend=CAP2), arrow=arrow(length=unit(0.30,"cm"), type = "open" ), size=1)
plot_dbRDA_fts = plot_dbRDA_fts + geom_point(aes(pch = Culture_time),fill = "cornflowerblue", size =7, stroke = 2, alpha = 0.9,  color = "black", show.legend = T)
plot_dbRDA_fts = plot_dbRDA_fts + scale_shape_manual(values = c(21, 23, 22,24, 25))
plot_dbRDA_fts = plot_dbRDA_fts + geom_text(data = data_dbRDA_fts_arrow ,cex = 6, color = "black" ,mapping = aes(CAP1 , CAP2) , label = row.names(data_dbRDA_fts_arrow))
plot_dbRDA_fts = plot_dbRDA_fts + theme_bw(base_size = 17) 
plot_dbRDA_fts = plot_dbRDA_fts + theme(legend.position="none")
plot_dbRDA_fts # Figure saved as Figure 3D


#### 27. FTS - Variance partioning analysis with env parameters ####

head(data_distbeta_fts)
dim(data_distbeta_fts)

head(metadata_fts)

varpart_fts = varpart(data_distbeta_fts , ~ Temperature ,	~Salinity,~ PO4 + NOx  +	SiOHx , data = metadata_fts)
varpart_fts


showvarparts(3, bg=2:4)
plot(varpart_fts, bg=2:4)
# Figure saved as Figure S6 (SI)


#### 28. FTS - Barplot analysis  of the total community at the family level ####

physeq_fts_compo_family = aggregate_rare(physeq_fts_compo, level = "Family", detection = 3/100, prevalence = 0/100)
physeq_fts_compo_family

ntaxa(physeq_fts_compo_family)

# for a barplot without border the phyloseq function plot_bar need to be changed. 
#Use the plot_bar_2 function as follows

plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

barplot_fts_family = plot_bar_2(physeq_fts_compo_family ,"Replicate_name" ,fill = "Family") + theme_bw(base_line_size = 12)
barplot_fts_family = barplot_fts_family + geom_bar(stat = "identity"  , position="stack", size = 3)
barplot_fts_family = barplot_fts_family + scale_fill_manual(values = colorpal_family)
barplot_fts_family = barplot_fts_family + theme(legend.text = element_text(size=9),
                                                axis.ticks.y = element_blank(),
                                                legend.position="right",
                                                axis.text.x = element_blank(),
                                                legend.title = element_blank(),
                                                axis.ticks.x=element_blank(), 
                                                axis.title = element_blank())
barplot_fts_family = barplot_fts_family + facet_grid(cols = vars(Culture_time), scales = "free", space = "free")
barplot_fts_family = barplot_fts_family + guides(fill = guide_legend(ncol = 1)) 
barplot_fts_family = barplot_fts_family + scale_y_percent()
# barplot_fts_family = barplot_fts_family + coord_flip()
barplot_fts_family = barplot_fts_family + theme(axis.title = element_blank())
barplot_fts_family # Figure saved as Figure 4B




#### 29. Subsampling for the of the study of the CORE community (BGS + FTS) ####


physeq_bgsfts = subset_samples(physeq, Study == "Biogeography" | Study == "Aquaculture")
physeq_bgsfts


#### 30. CORE - Creation of phyloseq objects for the core community  ####


# one phyloseq object with the rarefied dataset for alpha div analyses

physeq_bgsfts_rarefied = rarefy_even_depth(physeq_bgsfts)
physeq_bgsfts_rarefied

physeq_core_rarefied = core(physeq_bgsfts_rarefied, detection = 0.001/100, prevalence = 90/100)
physeq_core_rarefied

# the other one with the compositional dataset (percentages) for compositional analyses

physeq_bgsfts_compo = transform(physeq_bgsfts, "compositional")
physeq_bgsfts_compo

ntaxa(physeq_bgsfts_compo)

physeq_core_compo = core(physeq_bgsfts_compo, detection = 0.001/100, prevalence = 90/100)
physeq_core_compo

data_taxtable_ASVcore  = tax_table(physeq_core_compo)
head(data_taxtable_ASVcore) # Table saved as Table S4

ntaxa(physeq_core_compo)


#### 31. CORE - Barplot analysis  of the CORE community at the GENUS level #########

physeq_core_compo_genus = tax_glom(physeq_core_compo, taxrank = "Genus")
physeq_core_compo_genus

tax_table(physeq_core_compo_genus)

barplot_core_genus = plot_bar_2(physeq_core_compo_genus ,"Replicate_name" ,fill = "Genus") + theme_bw()
barplot_core_genus = barplot_core_genus + geom_bar(stat = "identity" , position="stack")
barplot_core_genus = barplot_core_genus + theme(legend.text = element_text(size=9),
                                                axis.ticks.y = element_blank(),
                                                legend.position="bottom",
                                                axis.text.x = element_blank(),
                                                legend.title = element_blank(),
                                                axis.ticks.x=element_blank(), 
                                                axis.title = element_blank())
barplot_core_genus = barplot_core_genus + facet_grid(cols = vars(Region_site), scales = "free", space = "free")
barplot_core_genus = barplot_core_genus + guides(fill = guide_legend(ncol = 1)) 
barplot_core_genus = barplot_core_genus+  scale_y_percent()
barplot_core_genus = barplot_core_genus + theme(axis.title = element_blank())
barplot_core_genus = barplot_core_genus + scale_fill_manual(values = c("darkorange", 
                                                                        "darkgoldenrod1", 
                                                                        "gold", 
                                                                        "limegreen" ,
                                                                        "seagreen4", 
                                                                        "plum1" , 
                                                                        "salmon" , 
                                                                        "slategray2" , 
                                                                        "aquamarine4" , 
                                                                        "yellow1", 
                                                                        "mediumpurple2", 
                                                                        "orchid", 
                                                                        "cyan"))
barplot_core_genus # Figure saved as Figure 5


#### 32. CORE - Statistical analysis of the percentage of sequences of CORE ASVs #########

ntaxa(physeq_core_compo)

percentages_coresequences = sample_sums(physeq_core_compo)
head(percentages_coresequences)

shapiro.test(percentages_coresequences)

data_coresequences = cbind(percentages_coresequences, sample_data(physeq_core_compo))
head(data_coresequences)


data_anova_sequences_core <- aov(percentages_coresequences ~ Region_site, data_coresequences)
data_anova_sequences_core
summary(data_anova_sequences_core)



#### 33. CORE - Statistical analysis of the percentage of numbers of CORE ASVs #########


ntaxa(physeq_core_rarefied)

ntaxa_core = estimate_richness(physeq_core_rarefied , measures = c("Observed"))
ntaxa_core

ntaxa_all = estimate_richness(physeq_bgsfts_rarefied , measures = c("Observed"))
ntaxa_all


percentage_coreASV = ntaxa_core$Observed / ntaxa_all$Observed
percentage_coreASV

data_coreASV = cbind(sample_data(physeq_core_rarefied), percentage_coreASV)
data_coreASV

min(data_coreASV$percentage_coreASV)
max(data_coreASV$percentage_coreASV)
mean(data_coreASV$percentage_coreASV)
sd(data_coreASV$percentage_coreASV)

data_coreASV_bgs= subset(data_coreASV, Study == "Biogeography")
data_coreASV_bgs

data_coreASV_fts = subset(data_coreASV, Study == "Aquaculture" | Region_site == "05.FP_Rangiroa_Avatoru")
data_coreASV_fts

data_coreASV_bgs = data_coreASV_bgs[, c("Region_site", "percentage_coreASV")]
data_coreASV_bgs

data_coreASV_fts = data_coreASV_fts[, c("Culture_time", "percentage_coreASV")]
data_coreASV_fts

data_coreASV_bgs_long = melt(data_coreASV_bgs, id.var=c("Region_site"))
data_coreASV_bgs_long

data_coreASV_fts_long = melt(data_coreASV_fts, id.var=c("Culture_time"))
data_coreASV_fts_long


plot_coreASV_bgs = ggplot(data_coreASV_bgs_long, aes(Region_site ,value, fill = Region_site), ymin = 0.04, ymax = 0.09)
plot_coreASV_bgs = plot_coreASV_bgs + geom_boxplot(alpha = 0.8, size = 1) #+ facet_grid(  variable ~ ., scales="free") 
plot_coreASV_bgs = plot_coreASV_bgs + geom_point(size = 2, alpha = 0.8, pch = 21, stroke = 1)
plot_coreASV_bgs = plot_coreASV_bgs + theme_bw(base_size = 15) 
plot_coreASV_bgs = plot_coreASV_bgs + theme(legend.position="left")
plot_coreASV_bgs = plot_coreASV_bgs + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_coreASV_bgs = plot_coreASV_bgs + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
plot_coreASV_bgs = plot_coreASV_bgs + scale_fill_manual(values = colorpal_site )
plot_coreASV_bgs #Figure saved as Figure S8A

plot_coreASV_fts = ggplot(data_coreASV_fts_long, aes(Culture_time ,value),ymin = 0.04, ymax = 0.09)
plot_coreASV_fts = plot_coreASV_fts + geom_boxplot(alpha = 0.7, size = 1,  fill = "cornflowerblue") #+ facet_grid(  variable ~ ., scales="free") 
plot_coreASV_fts = plot_coreASV_fts + geom_point(aes(pch = Culture_time), fill = "cornflowerblue" , size = 3, alpha = 0.8, stroke = 1)
plot_coreASV_fts = plot_coreASV_fts + theme_bw(base_size = 15) 
plot_coreASV_fts = plot_coreASV_fts + theme(legend.position="right")
plot_coreASV_fts = plot_coreASV_fts + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_coreASV_fts = plot_coreASV_fts + scale_shape_manual(values = c(21, 23, 22,24, 25) )
plot_coreASV_fts #Figure saved as Figure S8B


shapiro.test(data_coreASV_bgs$percentage_coreASV)

data_anova_coreASV_bgs <- aov(percentage_coreASV ~ Region_site, data_coreASV_bgs)
data_anova_coreASV_bgs
summary(data_anova_coreASV_bgs)

hsd_data_coreASV_bgs = HSD.test(data_anova_coreASV_bgs, "Region_site", group=T)
hsd_data_coreASV_bgs

shapiro.test(hsd_data_coreASV_fts$data_coreASV_fts)


data_anova_coreASV_fts <- aov(percentage_coreASV ~ Culture_time, data_coreASV_fts)
data_anova_coreASV_fts
summary(data_anova_coreASV_fts)

hsd_data_coreASV_fts = HSD.test(data_anova_coreASV_fts, "Culture_time", group=T)
hsd_data_coreASV_fts


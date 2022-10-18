library(qiime2R)
library(phyloseq)
library(vegan)
library(ggplot2)
library(tidyverse)
library(randomcoloR)
library(data.table)
library(reshape2)
library(ggpubr)
library(MiscMetabar)


#functions to extract plots(ggplot from lists)

extract_plots_from_list <- function(x) {
  
  if (!is.list(x)) {
    stop("`x` must be a list.")
  }
  
  if (length(x) < 1) {
    return(x)
  }
  
  more_lists <- sapply(x, function(object_to_check) {
    'list' %in% class(object_to_check) & !(any(c('gg', "gTree", "gtable") %in% class(object_to_check)))
  })
  
  result <- c(x[!more_lists],
              unlist(x[more_lists],
                     recursive = FALSE))
  
  if (sum(more_lists)) { 
    
    Recall(result)
    
  } else {
    
    is_plot_obj <- sapply(result, function(result_object) {
      any(c('gg', "gTree", "gtable") %in% class(result_object))
    })
    
    result <- result[is_plot_obj]
    return(result)
  }
}

#import data
setwd("~/Desktop/lusofona/results")

SVs <- read_qza("~/Desktop/lusofona/results/table.qza")
tree <- read_qza("~/Desktop/lusofona/results/rooted-tree.qza")

#Import metadata
metadata <- read_tsv("~/Desktop/lusofona/meta_tx.txt")
metadata<-metadata[-c(5:8)]
metadata<-as.data.frame(metadata)

#Import taxonomy and manipulate it to be more easy to work with
taxonomy <- read_qza("~/Desktop/lusofona/results/taxonomy.qza")
###dont run this for lefse names extraction
taxonomy <- as.data.frame(taxonomy$data)
taxonomy<-taxonomy[-which(names(taxonomy) %in% c("Confidence"))]
taxonomy <- taxonomy %>% separate(Taxon, sep=";",c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#remove taxa name appendages
taxonomy$Kingdom<-gsub("d__","",taxonomy$Kingdom)
taxonomy$Phylum<-gsub("p__","",taxonomy$Phylum)
taxonomy$Class<-gsub("c__","",taxonomy$Class)
taxonomy$Order<-gsub("o__","",taxonomy$Order)
taxonomy$Family<-gsub("f__","",taxonomy$Family)
taxonomy$Genus<-gsub("g__","",taxonomy$Genus)
taxonomy[taxonomy==" "]<-NA

#crete phyloseq object
phy_otu <- phyloseq(
  otu_table(SVs$data, taxa_are_rows = TRUE),
  phy_tree(tree$data),
  tax_table(taxonomy %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving taxonomy to the way phyloseq wants it.
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample-id")),
  refseqfilename = "rep-seqs.fasta"
)

##filter out mitochondria
phy_otu <- phy_otu %>% subset_taxa( Family!= "Mitochondria" | is.na(Family))

#sequences distribution

##otu histogram
topN <- 100
most_abundant_taxa = sort(taxa_sums(phy_otu), TRUE)[1:topN]

#extract most_abundant_taxa from mouse_data
mouse_100_OTUs <- prune_taxa(names(most_abundant_taxa), phy_otu)

# create a dataframe with the counts per otu
mouse_otu_sums <- data.frame(taxa_sums(mouse_100_OTUs))


pdf("OTUs_distribution.pdf", width = 10, height = 6)
ggplot(mouse_otu_sums,aes(reorder(row.names(mouse_otu_sums),-taxa_sums.mouse_100_OTUs.),taxa_sums.mouse_100_OTUs.)) + 
  geom_histogram(stat="identity",colour="black",fill="darkturquoise")  +
  xlab("ASV Rank") + ylab("Number of Sequences per ASV") +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + theme_classic()+ 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()

#Remove singletons and doubletons, and blanks
ps_phy2<-phy_otu

#singletons
tdt = data.table(tax_table(ps_phy2),
                 TotalCounts = taxa_sums(ps_phy2),
                 OTU = taxa_names(ps_phy2))

#number of singletons
tdt[(TotalCounts <= 0), .N]

#number of doubletons
tdt[(TotalCounts <= 2), .N]

#remove singletons and doubletons
ps_phy2<-prune_taxa(taxa_sums(ps_phy2) >= 2, ps_phy2)

#remove blanks
ps_phy2<-subset_samples(ps_phy2, Treatment !="blank")


#rarefaction
pdf("rarefaction.pdf",height = 7,width = 7)
rarecurve(t(otu_table(ps_phy2)), step = 1000,label= T,ylab="ASVs")
abline(v=(min(rowSums(t(otu_table(ps_phy2))))), col="red", lty=2)
dev.off()

set.seed(99)
ps_rarefied<-rarefy_even_depth(ps_phy2, min(sample_sums(ps_phy2)))

##after rarefaction
rarecurve(t(otu_table(ps_rarefied)), step = 1000, label = T,ylab="ASVs")

#alpha diversity
alpha.diversity <-estimate_richness(ps_rarefied, measures = c("Shannon","InvSimpson","Observed")) 
alpha.diversity


###add pielou eveness
H <- alpha.diversity$Shannon
S1 <- alpha.diversity$Observed
S <- log(S1)
evenness <- H/S
evenness
alpha.diversity$Pielou = evenness
alpha.diversity

data1 <- cbind(sample_data(ps_rarefied), alpha.diversity)
data1$Samples<-rownames(data1)
data1<-as.data.frame(data1)


#subset_data og mouse and pigs
mice_colitis<-data1[data1$Animal=="mouse" & data1$Health_status=="colitis",]
mice_healthy<-data1[data1$Animal == "mouse" & data1$Health_status == "healthy",]
pig<-data1[data1$Animal=="pig",]

status<-list(mice_colitis,mice_healthy,pig)
long_status<-list()

#convert wide format to long format
for (i in 1:length(status)){
  
  long_status[[i]] <- melt(setDT(status[[i]]), 
                           id.vars = c("Samples","Health_status","Treatment","Animal"), 
                           variable.name = "Diversity", 
                           measure.vars= c("Shannon","Pielou","InvSimpson","Observed"))
}


names(long_status) = c("mouse_colitis","mouse_healthy","pig")

ggcenas<-list()

for (i in 1:length(long_status)){
  
  ggcenas[[i]]<-ggplot(data = long_status[[i]], aes(x = Treatment, y = value, fill=Treatment)) +
    geom_point(colour="black",pch=21, size=2)+
    labs(x = "", y = "") +
    facet_wrap(~Diversity, scales = "free", nrow = 2)+
    theme_cowplot(font_size = 10,)+
    theme(strip.background = element_blank())
  
  ggsave(ggcenas[[i]],filename=paste(as.character(names(long_status)[i]),"inv_",".pdf"), width = 8, height = 8)
}

#save table alpha diversity
write.table(data1,"alpha_Diversity_Table.txt")


#beta diversity
#pig data
pigy<-subset_samples(ps_rarefied, Animal=="pig")
m_colitis<-subset_samples(ps_rarefied,Animal=="mouse" & Health_status=="colitis")
m_healthy<-subset_samples(ps_rarefied,Animal=="mouse" & Health_status=="healthy")

mouse<-subset_samples(ps_rarefied,Animal=="mouse")

merged<-list(pigy,m_colitis,m_healthy)
names(merged) = c("pig","mouse_colitis","mouse_healthy")

###use all vegan/phyloseq dissimilarities/distances
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
dist_methods["designdist"]
dist_methods = dist_methods[which(dist_methods=="unifrac" | dist_methods== "wunifrac" | dist_methods=="bray" | dist_methods=="jaccard")]

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods

dist<-list()
final<-list()
disper<-list()
per<-list()
permu<-list()
permutt<-list()
permanova<-list()
f_perma<-list()
NMDS<-list()
p<-list()
plist<-list()
plet<-list()
final_plet<-list()

#loop for distance calculation, betadispersion, permanova, NMDS
set.seed(1)
for( i in 1:length(merged)){
  for (j in dist_methods){
    
    ##distances
    if (j=="jaccard"){
      dist[[j]] <- distance(merged[[i]], method="jaccard", binary = TRUE)
      
    }else{
      dist[[j]]<-distance(merged[[i]], method=j)
    }
    
    final[[i]] <- dist
    names(final)[i]<-paste(names(merged[i]))
    
    #dispersion
    sampledf <- data.frame(sample_data(merged[[i]]))
    disper[[j]]<-betadisper(final[[i]][[j]],sampledf$Treatment)
    per[[i]] <- disper
    
    names(per)[i]<-paste(names(merged[i]))
    
    #permutests
    permu[[j]]<-permutest(per[[i]][[j]])
    permutt[[i]] <- permu
    
    names(permutt)[i]<-paste(names(merged[i]))
    
    #plots
    pdf(paste0("beta_disper_",paste(names(merged[i])),"_",j,".pdf"),width = 7, height =7)
    plot(per[[i]][[j]], main= paste0(j))
    dev.off()
    
    #PERMANOVA
    permanova[[j]]<-adonis(final[[i]][[j]] ~ Treatment, data = sampledf)
    f_perma[[i]]<-permanova
    
    names(f_perma)[i]<-paste(names(merged[i]))
    
    #NMDS
    NMDS[[j]] <-ordinate(merged[[i]], method ="NMDS",distance = j)
    p[[i]]<-NMDS
    
    names(p)[i]<-paste(names(merged[i]))
    
    # create plots
    plet[[j]] <- plot_ordination(merged[[i]], p[[i]][[j]], color="Treatment")
    
    final_plet[[i]]<-plet
    
    names(final_plet)[i]<-paste(names(merged[i]))
    
  }
}

#descrease one level of nested list NMDS plots
cena<-extract_plots_from_list(final_plet)

#stressplots for each NMDS
pdf("stressplot_mice.pdf", width = 7, height =7)
for (i in 1:length(p)){
  for(j in p[[i]]){
    
    stressplot(j)
  }
}
dev.off()

#change and save NMDS plots
pdf("NMDS_all.pdf", width = 7, height =7)
for (i in seq(cena)) {
  gg<-cena[[i]]+
    geom_point(aes(NMDS1,NMDS2, fill=Treatment),size=4, colour="black", shape=21)+
    theme_cowplot(font_size = 13)+
    ggtitle(names(cena[i]))
  print(gg)
}
dev.off()


#print to table all permutation tests together
write.table(capture.output(permutt),"betadisper_tests.txt",quote = F)

#print to table all permanovas togheter
write.table(capture.output(f_perma),"permanova_tests.txt",quote = F)

#print NMDS results
write.table(capture.output(p),"NMDS_results.txt",quote = F)

#empty lists
dist<-list()
final<-list()
disper<-list()
permu<-list()
permanova<-list()
NMDS<-list()
plet<-list()

set.seed(1)

#beta-diversity for mouse data
for (j in dist_methods){
  
  ##distances
  if (j=="jaccard"){
    dist[[j]] <- distance(mouse, method="jaccard", binary = TRUE)
    
  }else{
    dist[[j]]<-distance(mouse, method=j)
  }
  
  final <- dist
  
  
  #dispersion
  sampledf <- data.frame(sample_data(mouse))
  disper[[j]]<-betadisper(final[[j]],sampledf$Health_status)
  
  #permutests
  permu[[j]]<-permutest(disper[[j]])
  
  
  #plots
  pdf(paste0("beta_disper_",j,"_","mouse",".pdf"),width = 7, height =7)
  plot(disper[[j]], main= paste0(j))
  dev.off()
  
  #PERMANOVA
  permanova[[j]]<-adonis(final[[j]] ~ Health_status, data = sampledf)
  
  
  #NMDS
  NMDS[[j]] <-ordinate(mouse, method ="NMDS",distance = j)
  
  # create plots
  plet[[j]] <- plot_ordination(mouse, NMDS[[j]], color="Health_status")
  
}

#descrease one level of nested list NMDS plots
cena1<-extract_plots_from_list(plet)

#stressplots 
pdf("stressplot_micem.pdf", width = 7, height =7)
for (i in 1:length(NMDS)){
  stressplot(NMDS[[i]])
}
dev.off()


#NMDS plots changes and save
pdf("NMDS_all.pdf", width = 8, height =7)
for (i in seq(cena1)) {
  gg<-cena1[[i]]+
    geom_point(aes(NMDS1,NMDS2, fill=Health_status),size=4, colour="black", shape=21)+
    theme_cowplot(font_size = 13)+
    ggtitle(names(cena1[i]))
  print(gg)
}
dev.off()


#print to table all permutation tests together
write.table(capture.output(permu),"betadisper_tests.txt",quote = F)

#print to table all permanovas togheter
write.table(capture.output(permanova),"permanova_tests.txt",quote = F)

#print NMDS results
write.table(capture.output(NMDS),"NMDS_results.txt",quote = F)

#accumulation curves
accu_plot(ps_rarefied, fact = "Treatment", by.fact = T)


#Linear discriminant analysis
set.seed(99)

ps_glom_new1<-tax_glom(mouse,taxrank = "Genus")

#LDA=2
lef1_new1<-run_lefse(ps_glom_new1,
                     wilcoxon_cutoff = 0.05,
                     taxa_rank = "Genus",
                     norm = "CPM",
                     group = "Treatment",
                     kw_cutoff = 0.05,
                     lda_cutoff = 2)

pdf(paste("LDA_2score","_mouse_",".pdf",height = 12,width = 10)
    plot_ef_bar(lef1_new1)+
      geom_bar(stat = "identity", color="black", size=0.5)+
      scale_fill_manual(values = c("grey48","red"))+
      theme_cowplot(font_size = 14)+
      theme(panel.grid.major.x = element_line( size=0.5, color="black", linetype = "dashed" ),
            panel.grid.major.y = element_blank(),
            text = element_text(size=14),
            axis.title=element_text(size=14))+
      guides(fill=guide_legend(title="Health_status"))
    dev.off()
    

#LDA score =4
    
set.seed(99)
lef2_rarefied1<-run_lefse(ps_glom_new1,
                          wilcoxon_cutoff = 0.05,
                          taxa_rank = "Genus",
                          norm = "CPM",
                          group = "Health_status",
                          kw_cutoff = 0.05,
                          lda_cutoff = 4)
    
pdf("LDA_4score_genus.pdf",height = 6,width = 10)
plot_ef_bar(lef2_rarefied1)+
geom_bar(stat = "identity", color="black", size=0.5)+
scale_fill_manual(values = c("grey48","red"))+
theme_cowplot(font_size = 14)+
theme(panel.grid.major.x = element_line( size=0.5, color="black", linetype = "dashed" ),
            panel.grid.major.y = element_blank(),
            text = element_text(size=14),
            axis.title=element_text(size=14))+
guides(fill=guide_legend(title="Health_status"))
dev.off()
    
    
phy_relative_rarefied_new<-ps_glom_new1%>%           
transform_sample_counts(function(x) x/sum(x) )
    
phy_set_new<-psmelt(phy_relative_rarefied_new)
    
phy_set9_new<-subset(phy_set_new, Genus %in% c(" Lachnospiraceae_NK4A136_group",
                                                   " Colidextribacter",
                                                   " Escherichia-Shigella",
                                                   " Prevotellaceae_UCG-001"))
    
#adjust pvalue
stat.test <- phy_set9_new %>%
group_by(Genus) %>%
wilcox_test(Abundance ~ Health_status) %>%
adjust_pvalue(method = "fdr") %>%
add_significance()
    
# Make facet and add p-values
stat.test <- stat.test %>% add_xy_position(x ="Genus")
    
phy_set9_new$Health_status <- as.factor(phy_set9_new$Health_status)

#create graph    
p<-ggboxplot(
  phy_set9_new, x = "Health_status", y = "Abundance", fill ="Health_status", outlier.shape=NA)+  facet_wrap(~Genus, scales = "free")+
  geom_jitter(aes(fill = Health_status), height = 0, width = .2, colour="black",pch=21, size=2)+ 
  labs(x = "", y = "Abundance (%)") +
  scale_fill_manual(values = c("grey48","red"))+
  scale_color_manual(values = c("grey48","red"))+
  theme_cowplot(font_size = 13,)+
  #facet_wrap(~Genus, scales = "free")+
  theme(strip.background = element_blank())

#add significances to graphs and save them  
pdf("lefse_genus_LDA_4_white.pdf", width=8, height=8)  
p + stat_pvalue_manual(stat.test)
dev.off()


#barplots phyla and Genus
#separate pig from mouse
    
pigz<-subset_samples(ps_phy2, Animal == "pig")
    
mouse<-subset_samples(ps_phy2, Animal =="mouse")
    
    
##pigz
#Phyla
pigz_p <- pigz%>%
tax_glom(taxrank = "Phylum") %>%                     
transform_sample_counts(function(x) {x/sum(x)} ) %>% 
psmelt() %>%
arrange(Phylum)
pigz_p$Phylum<- as.character(pigz_p$Phylum)
pigz_p$Phylum[pigz_p$Abundance < 0.01] <- "Others"
    
pigz_p$Sample<-factor(pigz_p$Sample, levels = c("P1F","P2F","P3F","P19F","P20F","P22F"))
    
pigz_p$Treatment<-factor(pigz_p$Treatment, levels = c("control","deflamine"))
    
#Family  
pigz_f <- pigz%>%
  tax_glom(taxrank = "Family") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%
  arrange(Family)
  pigz_f$Family<- as.character(pigz_f$Family)
  pigz_f$Family[pigz_f$Abundance < 0.01] <- "Others"
    
    
pigz_f$Sample<-factor(pigz_f$Sample, levels = c("P1F","P2F","P3F","P19F","P20F","P22F"))
    
pigz_f$Treatment<-factor(pigz_f$Treatment, levels = c("control","deflamine"))
    
#graph phyla
pdf("relative_abundance_Phylum_pigs.pdf",width = 10, height = 7)
    ggplot(pigz_p, aes(x = Sample, y = Abundance, fill = Phylum))+
      scale_y_continuous(labels = c("0","25","50","75","100"))+
      geom_bar(stat="identity", colour="black")+
      scale_fill_manual(values = c("Others" = "black",
                                   " Firmicutes"="tomato1",
                                   " Bacteroidota"="gold2",
                                   " Cyanobacteria"="olivedrab1",
                                   " Proteobacteria"="skyblue"))+
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab("Relative Abundance (%)  \n")+
      theme_minimal(base_size =14)+
      facet_wrap(~Treatment, scales="free_x", ncol = 5, labeller = label_both)+
      theme(axis.text.x = element_text(angle = 0,vjust = 1.2), 
            panel.grid.major = element_line(colour = "white"), 
      )+
      labs(title = NULL)+
      guides(fill=guide_legend(ncol=1))
dev.off()
    
#check families
unique(pigz_f$Family)
    
# no. of colours in the palette
no_of_colors <- 14
    
# sample colors
palette <- distinctColorPalette(no_of_colors)     
    
# hex color codes
palette
    
# colors in the pie chart
pie(rep(1, no_of_colors), col = palette, 
main = "palette using randomcoloR package")
    
#graph family
pdf("relative_abundance_pigs_Family.pdf",width = 10, height = 7)
  ggplot(pigz_f, aes(x = Sample, y = Abundance, fill = Family))+
  scale_y_continuous(labels = c("0","25","50","75","100"))+
  geom_bar(stat="identity", colour="black")+
  scale_fill_manual(values= c("Others"="black",
                                  " Acidaminococcaceae"="#7DDF8C",
                                  " Clostridia_UCG-014"= "#B248E3",
                                  " Gastranaerophilales"="blue",
                                  " Lachnospiraceae"="#D7AAD8",   
                                  " Lactobacillaceae"="#E3A85E",  
                                  " Oscillospiraceae"="brown",
                                  " Prevotellaceae"="skyblue",     
                                  " Rikenellaceae"="#D962B6",      
                                  " Ruminococcaceae"="#D0D799",  
                                  " Selenomonadaceae"="#D96065",   
                                  " Streptococcaceae"="#84EE4E",   
                                  " Succinivibrionaceae"="#D4DE52",
                                  " uncultured"="#CE9F97",         
                                  " Veillonellaceae"="#D5DEDB"))+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (%)  \n")+
  theme_minimal(base_size =14)+
  facet_wrap(~Treatment, scales="free_x", ncol = 5, labeller = label_both)+
  theme(axis.text.x = element_text(angle = 0,vjust = 1.2), 
  panel.grid.major = element_line(colour = "white"), 
      )+
  labs(title = NULL)+
  guides(fill=guide_legend(ncol=2))
dev.off()
    
#graph phyla
pdf("relative_abundance_Phylum_pigs.pdf",width = 10, height = 7)
    ggplot(pigz_p, aes(x = Sample, y = Abundance, fill = Phylum))+
      scale_y_continuous(labels = c("0","25","50","75","100"))+
      geom_bar(stat="identity", colour="black")+
      scale_fill_manual(values = c("Others" = "black",
                                   " Firmicutes"="tomato1",
                                   " Bacteroidota"="gold2",
                                   " Cyanobacteria"="olivedrab1",
                                   " Proteobacteria"="skyblue"))+
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab("Relative Abundance (%)  \n")+
      theme_minimal(base_size =14)+
      facet_wrap(~Treatment, scales="free_x", ncol = 5, labeller = label_both)+
      theme(axis.text.x = element_text(angle = 0,vjust = 1.2), 
            panel.grid.major = element_line(colour = "white"), 
      )+
      labs(title = NULL)+
      guides(fill=guide_legend(ncol=1))
dev.off()
    
#Phyla mouse
mouse_p <- mouse%>%
      tax_glom(taxrank = "Phylum") %>%                     
      transform_sample_counts(function(x) {x/sum(x)} ) %>% 
      psmelt() %>%
      arrange(Phylum)
    mouse_p$Phylum<- as.character(mouse_p$Phylum)
    mouse_p$Phylum[mouse_p$Abundance < 0.01] <- "Others"
    
    
mouse_p$Treatment<-factor(mouse_p$Treatment, levels = c("control","deflamine", "extract","peptide"))
    
mouse_p$Health_status<-factor(mouse_p$Health_status, levels = c("healthy","colitis"))
    
#Family mouse
mouse_f <- mouse%>%
      tax_glom(taxrank = "Family") %>%                     
      transform_sample_counts(function(x) {x/sum(x)} ) %>% 
      psmelt() %>%
      arrange(Family)
    mouse_f$Family<- as.character(mouse_f$Family)
    mouse_f$Family[mouse_f$Abundance < 0.01] <- "Others"
    
mouse_f$Sample<-factor(mouse_f$Sample, levels = c("P1F","P2F","P3F","P19F","P20F","P22F"))
    
mouse_f$Treatment<-factor(mouse_f$Treatment, levels = c("control","deflamine", "extract","peptide"))
    
mouse_f$Health_status<-factor(mouse_f$Health_status, levels = c("healthy","colitis"))
    
unique(mouse_p$Phylum)
    
#graph phyla
pdf("relative_abundance_Phylum_mouse.pdf",width = 17, height = 7)
    ggplot(mouse_p, aes(x = Sample, y = Abundance, fill = Phylum))+
      scale_y_continuous(labels = c("0","25","50","75","100"))+
      geom_bar(stat="identity", colour="black")+
      scale_fill_manual(values = c("Others"="black",             
                                   " Bacteroidota"="moccasin",      
                                   " Cyanobacteria"="olivedrab3",    
                                   " Deferribacterota"="tomato1",  
                                   " Desulfobacterota"="lightseagreen",  
                                   " Firmicutes"="skyblue",       
                                   " Proteobacteria"="red",    
                                   " Verrucomicrobiota"="khaki1"))+
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab("Relative Abundance (%)  \n")+
      theme_minimal(base_size =14)+
      facet_wrap(Treatment~Health_status, scales="free_x", ncol = 8, labeller = label_both)+
      theme(axis.text.x = element_text(angle = 0,vjust = 1.2), 
            panel.grid.major = element_line(colour = "white"), 
      )+
      labs(title = NULL)+
      guides(fill=guide_legend(ncol=1))
dev.off()
    
#check families
unique(mouse_f$Family)
    
# no. of colours in the palette
no_of_colors <- 21
    
# sample colors
palette <- distinctColorPalette(no_of_colors)     
    
# hex color codes
palette
    
# colors in the pie chart
pie(rep(1, no_of_colors), col = palette, 
        main = "palette using randomcoloR package")
    
#graph family
pdf("relative_abundance_mouse_family.pdf",width = 17, height = 7)
   ggplot(mouse_f, aes(x = Sample, y = Abundance, fill = Family))+
   scale_y_continuous(labels = c("0","25","50","75","100"))+
   geom_bar(stat="identity", colour="black")+
   scale_fill_manual(values = c("Others"="black",                      
                                   " Akkermansiaceae"="#7C4A7D",           
                                   " Clostridia_UCG-014"="blue2",    
                                   " Clostridia_vadinBB60_group"="#7E66DE",    
                                   " Desulfovibrionaceae"="#8896DD",        
                                   " Enterobacteriaceae"="#7EE1B3",       
                                   " Gastranaerophilales"="#D6DD98",       
                                   " Lachnospiraceae"="lightseagreen",           
                                   " Lactobacillaceae"="sienna2",          
                                   " Marinifilaceae"="khaki1",           
                                   " Muribaculaceae"="goldenrod",            
                                   " Oscillospiraceae"="#82EF4A",         
                                   " Prevotellaceae"="#D14EE3",            
                                   " RF39"="#E1A853",                      
                                   " Rikenellaceae"="tomato1",
                                   " Deferribacteraceae"="#E06667",    
                                   " Ruminococcaceae"="#70E0E0",          
                                   " Sutterellaceae"="#7DAFD0",
                                   " Tannerellaceae"="#CEE9D8",           
                                   " uncultured"="#E1B096"))+
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance (%)  \n")+
    theme_minimal(base_size =14)+
    facet_wrap(Treatment~Health_status,scales="free_x",ncol = 8, labeller = label_both)+
    theme(axis.text.x = element_text(angle = 0,vjust = 1.2), 
            panel.grid.major = element_line(colour = "white"), 
      )+
    labs(title = NULL)+
    guides(fill=guide_legend(ncol=1))
dev.off()
    
    
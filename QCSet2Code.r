
################################
# Importing required Libraries #
################################

library(pheatmap)
library(RColorBrewer)
library(GUniFrac)
library(vegan)
library(phyloseq)
library(ade4)
library(car)
library(phangorn)
library(multcomp)
library(ggplot2)
library(XLConnect)
library(ggplotify)
library(cowplot)
library(lme4)
library(pbkrtest)

#################################
#Reading in OTU table and  tree #
#converting OTU table           # 
#from json to R dataframe       #
################################



OTU_tab2phylo <- import_biom("jsonOTUtax.biom", parseFunction = parse_taxonomy_greengenes)
OTU_tab2 <- as((OTU_tab2@otu_table),"matrix")
OTU_tab2 <- as.data.frame(t(OTU_tab2))

taxonomyT <- OTU_tab2@tax_table@.Data

##########################
# Upload additional data #
##########################

batdata2 <- readWorksheetFromFile("Batch4to9.xlsx",sheet = 1)

#################
# Cleaning Data #
#################

#Check number of sequences
sum(rowSums(OTU_tab2)) #4830153

#Trim the "020515" from the row.names
rownames(OTU_tab2) <- gsub(pattern="020515",replacement="",rownames(OTU_tab2))

#Load tree
tree2 <- read_tree_greengenes("rep_set.tre")
tree2 <- midpoint(tree2) # midpoint rooted

#Check the concordance between the OTU and tip.label names.

tree2$tip.label
names(OTU_tab)

#Identify OTUs in OTU table, but absent from tree
absent_from_tree <- names(OTU_tab2)[!(names(OTU_tab2) %in% tree2$tip.label)]  #2 OTUs
in_both <- names(OTU_tab2)[names(OTU_tab2) %in% tree2$tip.label]   #5977
#Subset to just those OTUs in the tree, i.e. "in_both"
OTU_trim2 <- subset(OTU_tab2,select=in_both)   #5977



################################
# Convert OTU_tab and OTU_trim #
#    to otu_table objects      #
################################

OTU_tab_object = otu_table(OTU_tab, taxa_are_rows = FALSE)
OTU_trim_object2 = otu_table(OTU_trim2, taxa_are_rows = FALSE)

######################
# Log Transform Data #
#  and make Heat Map #
######################

batOrder <- batdata[order(batdata$Batch),]
sampOrder <- batdata[order(batdata$Cont.ID),]
logOTUtrim2 <- log1p(OTU_trim_object2)


heatdat2 <- batdata2[,c("Batch","Subject")]
heatdat2$Batch <- gsub("ch", "ch ", heatdat2$Batch)
heatdat2$Subject <- gsub("([[:alpha:]])([[:digit:]])", "\\1 \\2", heatdat2$Subject)
heatdat2$Subject <- gsub("([[:digit:]])([[:alpha:]])", "\\1 \\2", heatdat2$Subject)
heatdat2$Batch <- as.factor(heatdat2$Batch)
heatdat2$Subject <- as.factor(heatdat2$Subject)
rownames(heatdat2) <- batdata2$Label
heatcol2 = list(Batch=c("Batch 4" ="deeppink3","Batch 5"="darkgreen",
                        "Batch 6" ="darkorchid","Batch 7"="springgreen3", "Batch 8" = "navy","Batch 9" = "turquoise3"),
                Subject = c("Child 1"= "darkorchid", "Child 2 Extract 1"="navy",
                            "Child 2 Extract 2" ="turquoise3","Child 3" ="deeppink3"))
names(heatcol2$Batch) <- unique(heatdat2$Batch)
names(heatcol2$Subject) <- unique(sort(heatdat2$Subject))

levels(heatdat2$Batch)
levels(heatdat2$Subject)

#plot_heatmap(logOTUtrim, sample.order = batOrder$Label)
#plot_heatmap(logOTUtrim, sample.order = sampOrder$Label)


p1 <- pheatmap(logOTUtrim2, color = colorRampPalette(rev(brewer.pal(n=11, name = "Spectral")))(100), 
                  border_color = NA,
                  cluster_rows = T, cluster_cols = T,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  scale = "none",
                  show_rownames = F, show_colnames = F, 
                  annotation_row = heatdat2, 
                  annotation_colors = heatcol2,
                  legend = FALSE,
                  width = 5,
                  height = 3,
                  fontsize = 10
              
         )


p1


#Create Fig 2 from JCM and QC heatmaps 

heat1 <- as.grob(p1)
heat2 <- as.grob(p2)
Fig2 <- plot_grid(heat2,heat1, ncol=2, labels = c("A","B"),align = "h", rel_widths = c(1,1.175))


tiff("Figure_2.tiff", width = 10, height =4,units = "in", res = 300)
Fig2
dev.off()


#############################
# Calculate alpha diversity #
#############################


#
alpha_div2 <- estimate_richness(OTU_trim_object2, split = TRUE)
rownames(alpha_div2) <- gsub("X","",rownames(alpha_div2))
batdata2 <- batdata2[match(row.names(alpha_div2),batdata2$Label),]
batdata2$Label == rownames(alpha_div2)
color = c("#000093" ,"#0076FF", "#00B8C2", 
          "#49FB25",  "#FEEA02", "#FFC200", "#FF8500",
          "#FF3300")
plot(as.factor(batdata2$Batch), alpha_div2$Chao1,  col = color,  main = "Chao1 index by batch")
plot(as.factor(batdata2$Batch), alpha_div2$Shannon,  col = color,  main = "Shannon index by batch")
plot(as.factor(batdata2$Batch), alpha_div2$Simpson,  col = color,  main = "Simpson index by batch")
summary(o<-lm(alpha_div2$Chao1 ~ as.factor(batdata2$Batch)))
summary(glm(alpha_div2$Shannon ~ as.factor(batdata2$Batch)))
summary(glm(alpha_div2$Simpson ~ as.factor(batdata2$Batch)))

alphadat2 <- cbind(alpha_div2, Batch = batdata2$Batch, Subject = batdata2$Subject)
alphadat2$Subject <- gsub("([[:alpha:]])([[:digit:]])", "\\1 \\2", alphadat2$Subject)
alphadat2$Subject <- gsub("([[:digit:]])([[:alpha:]])", "\\1 \\2", alphadat2$Subject)
alphadat1$Subject <- gsub("([[:alpha:]])([[:digit:]])", "\\1 \\2", alphadat1$Subject)

#increment subject number for combination with set 2
for (i in 5:1){
  alphadat1$Subject <- gsub(i,i+3, alphadat1$Subject )
}


alphadat <- rbind(alphadat1,alphadat2)

summary(o<-glm(alphadat$Simpson ~ as.factor(alphadat$Batch)))
anova(o)

summary(mod<-lmer(alphadat$Simpson ~ alphadat$Batch +(1|factor(alphadat$Subject))))
dfkr<- get_ddf_Lb(mod, fixef(mod))
coefs<- data.frame(coef(summary(mod)))
coefs$pkr <- 2*(1- pt(abs(coefs$t.value), dfkr))
coefs

anova(mod)
set.seed(3103)
rare2 <- rarefy_even_depth(OTU_trim_object2)
alpha2 <- estimate_richness(rare2)

plot(as.factor(batdata$Batch), alpha2$Chao1,  col = color,  main = "Chao1 index by batch post Rarefaction")
plot(as.factor(batdata$Batch), alpha2$Shannon,  col = color,  main = "Shannon index by batch post Rarefaction")
plot(as.factor(batdata$Batch), alpha2$Simpson,  col = color,  main = "Simpson index by batch post Rarefaction")

summary(glm(alpha2$Chao1 ~ as.factor(batdata2$Batch)))
summary(glm(alpha2$Shannon ~ as.factor(batdata2$Batch)))
summary(glm(alpha2$Simpson ~ as.factor(batdata2$Batch)))


eFig2 <- ggplot(data = alphadat, aes(x = Batch, y = Simpson))+
  geom_boxplot(fill = rep(c("#0000FF", "#66CC22", "#FF3399", 
                            "#6699FF", "#990099", "#FF9900", 
                            "darkgreen", "navy", "turquoise3")))+
  theme_bw(base_size = 12)+
  ylab("Simpson Diversity Index")+
  scale_x_discrete(name = "Batch number", labels = c(1:9))

eFig2

tiff("eFigure1.tiff", width = 6, height =4,units = "in", res = 300)
eFig2
dev.off()

############
# GUnifrac #
############
#normalize OTU table 
normOTU2 <- transform_sample_counts(OTU_trim_object2, function(x) x / sum(x))

#calculate GUnifracDistances
unifracs2 <- GUniFrac(normOTU2, tree2, alpha = c(0, 0.5, 1))$unifracs

#extracting distances that used alpha = 0.5
d5_2 <- unifracs2[,,"d_0.5"] 

#Ordering data by subject
SubjOrder <- batdata[order(batdata$Subject),]

#copy for testing
distan2 <- d5_2

#matching order
SubjOrder <- SubjOrder[match(row.names(distan),SubjOrder$Label),]
SubjOrder$Label == rownames(distan)

Subjects <- SubjOrder$Subject


batdata2 <- batdata2[match(row.names(distan2),batdata2$Label),]
batdata2$Label == rownames(distan2)
adonis(as.dist(distan2) ~ batdata2$Batch)


pcoa2<-pcoa(d5_2, correction="cailliez")
PCpercent <- pcoa2$values$Broken_stick[1:2]
pcoa2<-pcoa2$vectors[,1:2]

batdata2$Label == batdata2A$Label
batdata2$Subject <- batdata2A$Subject

batdata2$Subject <- factor(batdata2$Subject, levels = c("Child1" ,"Child3",
                                                        "Child2Extract1", 
                                                        "Child2Extract2"))


B1 <- ggplot(batdata2, aes(pcoa2[,1],pcoa2[,2])) + 
  geom_point(size = 2.5,aes( color = Subject, shape = Subject), alpha=0.75)  + 
  scale_shape_manual(values =c(15,17,18,19), name = "Subject", 
                     labels = c("Child 1" ,"Child 3",
                                "Child 2 Extract 1", 
                                "Child 2 Extract 2")) +
  scale_color_manual(values =c("darkorchid","deeppink3","navy", "turquoise3"), name = "Subject", 
                     labels = c("Child 1" ,"Child 3",
                                "Child 2 Extract 1", 
                                "Child 2 Extract 2"))+
  labs(x = "PC 1 (10%)", y = "PC 2 (8%)")+
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom", legend.title.align = 0.5)+
  guides(color=guide_legend(ncol=2, title.position = "top"))+
  guides(shape=guide_legend(ncol=2, title.position = "top"))
 
B1

B3B <- ggplot(batdata2[batdata2$Subject == "Child2Extract1"|batdata2$Subject == "Child2Extract2",], 
              aes(pcoa2[batdata2$Subject == "Child2Extract1"|batdata2$Subject == "Child2Extract2",1],
                  pcoa2[batdata2$Subject == "Child2Extract1"|batdata2$Subject == "Child2Extract2",2])) + 
  geom_point(size = 5,aes( color = Subject, shape = Subject), alpha = 0.75)  + 
  scale_shape_manual(values =c(17,19), name = "Extraction", 
                     labels = c("Extract 1","Extract 2")) +
  scale_color_manual(values =c("navy", "turquoise3"), name = "Extraction", 
                     labels = c("Extract 1", "Extract 2"))+
  labs(x = "PC 1 (21%)", y = "PC 2 (15%)")+
  theme_bw(base_size = 12) + 
  theme(legend.position = "bottom", legend.title.align = 0.5)+
  guides(color=guide_legend(ncol=2, title.position = "top"))+
  guides(shape=guide_legend(ncol=2, title.position = "top"))

B3B

Fig1 <- plot_grid(B2,B1, ncol=2, labels = c("A","B"),align = "h", rel_widths = c(1,1))

Fig1A <- plot_grid(B2,B1,B3B, ncol=3, labels = c("A","B","C"),align = "h", rel_widths = c(1,1,1))


tiff("Figure_1.tiff", width = 6.25, height =4,units = "in", res = 300)
Fig1
dev.off()


tiff("Figure_1A.tiff", width = 8, height =4,units = "in", res = 300)
Fig1A
dev.off()


tiff("Figure_1B.tiff", width = 8, height =4,units = "in", res = 300)
Fig1B
dev.off()

     
       

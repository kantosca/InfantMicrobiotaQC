
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
library(readxl)


#################################
#Reading in OTU table and  tree #
#converting OTU table           # 
#from json to R dataframe       #
################################


OTU_tab <- import_biom("OTUtable.from_txt_json.biom")
OTU_tab <- as(OTU_tab,"matrix")
OTU_tab <- as.data.frame(t(OTU_tab))

##########################
# Upload additional data #
##########################

batdata <- read_excel("JQC.xlsx")
#correcting colnames since there was an issue with a duplicate column. 
colnames(batdata) <- c("Batch", "Label", "SampleType", "ContID","Labels", "Subject")

#################
# Cleaning Data #
#################

#Check number of sequences
sum(rowSums(OTU_tab)) #6100208

#Trim the "020515" from the row.names
rownames(OTU_tab) <- gsub(pattern="020515",replacement="",rownames(OTU_tab))

#Load tree
tree <- read_tree_greengenes("rep_set.tre")
tree <- midpoint(tree) # midpoint rooted

#Check the concordance between the OTU and tip.label names.

tree$tip.label
names(OTU_tab)

#Identify OTUs in OTU table, but absent from tree
absent_from_tree <- names(OTU_tab)[!(names(OTU_tab) %in% tree$tip.label)]  #2 OTUs
in_both <- names(OTU_tab)[names(OTU_tab) %in% tree$tip.label]   #5977
#Subset to just those OTUs in the tree, i.e. "in_both"
OTU_trim <- subset(OTU_tab,select=in_both)   #5977


  
################################
# Convert OTU_tab and OTU_trim #
#    to otu_table objects      #
################################

OTU_tab_object = otu_table(OTU_tab, taxa_are_rows = FALSE)
OTU_trim_object = otu_table(OTU_trim, taxa_are_rows = FALSE)

######################
# Log Transform Data #
#  and make Heat Map #
######################

batOrder <- batdata[order(batdata$Batch),]
sampOrder <- batdata[order(batdata$Cont.ID),]
logOTUtrim <- log1p(OTU_trim_object)
heatdat <- batdata[,c("Batch","Subject")]
heatdat$Batch <- gsub("ch", "ch ", heatdat$Batch)
heatdat$Subject <- gsub("([[:alpha:]])([[:digit:]])", "\\1 \\2", heatdat$Subject)
heatdat$Subject <- gsub("([[:digit:]])([[:alpha:]])", "\\1 \\2", heatdat$Subject)
heatdat$Batch <- as.factor(heatdat$Batch)
heatdat$Subject <- as.factor(heatdat$Subject)
rownames(heatdat) <- batdata$Label
heatcol = list(Batch=c("Batch 1"="navy", "Batch 2"="turquoise3","Batch 3" ="deeppink3"),
               Subject = c("Child 1" ="deeppink3","Child 2" ="darkgreen",
                           "Child 3" = "darkorchid","Child 4" = "navy","Child 5" = "turquoise3"))
names(heatcol$Batch) <- unique(heatdat$Batch)
names(heatcol$Subject) <- unique(heatdat$Subject)
#plot_heatmap(logOTUtrim, sample.order = batOrder$Label)
#plot_heatmap(logOTUtrim, sample.order = sampOrder$Label)

#generate figure plot

p2 <- pheatmap(logOTUtrim, color = colorRampPalette(rev(brewer.pal(n=11, name = "Spectral")))(100), 
                  border_color = NA,
                  cluster_rows = T, cluster_cols = T,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  scale = "none",
                  show_rownames = F, show_colnames = F, 
                  annotation_row = heatdat, 
                  annotation_colors = heatcol,
                  legend = FALSE,
                  width = 4,
                  height = 3,
                  fontsize = 10
)
 p2


#############################
# Calculate alpha diversity #
#############################

color = c("#000093" ,"#0076FF", "#00B8C2", 
          "#49FB25",  "#FEEA02", "#FFC200", "#FF8500",
          "#FF3300")
alpha_div <- estimate_richness(OTU_trim_object, split = TRUE)
rownames(alpha_div) <- gsub("X","",rownames(alpha_div))
batdata <- batdata[match(row.names(alpha_div),batdata$Label),]
batdata$Label == rownames(alpha_div)

alpha_div$Shannon


plot(as.factor(batdata$Batch), alpha_div$Chao1,  col = color[6:8],  main = "Chao1 index by batch")
plot(as.factor(batdata$Batch), alpha_div$Shannon,  col = color[6:8],  main = "Shannon index by batch")
points(rep(1,5),alpha_div$Shannon[batdata$Batch=="Batch1"], col = c(1:5),pch=19 )
points(rep(2,5),alpha_div$Shannon[batdata$Batch=="Batch2"], col = c(1:5),pch=19 )
points(rep(3,5),alpha_div$Shannon[batdata$Batch=="Batch3"], col = c(1:5),pch=19 )
plot(as.factor(batdata$Batch), alpha_div$Simpson,  col = color[6:8],  main = "Simpson index by batch")
plot(as.factor(batdata$Subject), alpha_div$Shannon,  col = color[6:8],  main = "Shannon index by subject")
summary(glm(alpha_div$Chao1 ~ as.factor(batdata$Batch)))
summary(glm(alpha_div$Shannon ~ as.factor(batdata$Batch)))
summary(glm(alpha_div$Simpson ~ as.factor(batdata$Batch)))

alphadat1 <- cbind(alpha_div,Batch = batdata$Batch, Subject = batdata$Subject)



set.seed(3103)
rare <- rarefy_even_depth(OTU_trim_object)
alpha2 <- estimate_richness(rare)

plot(as.factor(batdata$Batch), alpha2$Chao1,  col = color[6:8],  main = "Chao1 index by batch post Rarefaction")
plot(as.factor(batdata$Batch), alpha2$Shannon,  col = color[6:8],  main = "Shannon index by batch post Rarefaction")
plot(as.factor(batdata$Batch), alpha2$Simpson,  col = color[6:8],  main = "Simpson index by batch post Rarefaction")

summary(glm(alpha2$Chao1 ~ as.factor(batdata$Batch)))
summary(glm(alpha2$Shannon ~ as.factor(batdata$Batch)))
summary(glm(alpha2$Simpson ~ as.factor(batdata$Batch)))


############
# GUnifrac #
############

#normalize OTU table 
normOTU <- transform_sample_counts(OTU_trim_object, function(x) x / sum(x))

#calculate GUnifracDistances
unifracs <- GUniFrac(normOTU, tree, alpha = c(0, 0.5, 1))$unifracs

#extracting distances that used alpha = 0.5
d5 <- unifracs[,,"d_0.5"] 

#Ordering data by subject
SubjOrder <- batdata[order(batdata$Subject),]

#copy for testing
distan <- d5

#matching order
SubjOrder <- SubjOrder[match(row.names(distan),SubjOrder$Label),]
SubjOrder$Label == rownames(distan)

Subjects <- SubjOrder$Subject

batdata<- batdata[match(row.names(d5),batdata$Label),]
batdata$Label == rownames(d5)
adonis(as.dist(d5) ~ batdata$Batch)


#ordination
pcoa1<-pcoa(d5, correction="cailliez")
PCpercent1 <- pcoa1$values$Broken_stick[1:2]
pcoa1<-pcoa1$vectors[,1:2]

#generate figure plot
B2 <- ggplot(batdata, aes(pcoa1[,1],pcoa1[,2])) + 
  geom_point(size = 2.5,aes( color = Subject, shape = Subject), alpha = 0.75)  + 
  scale_shape_manual(values =c(15,19,8,17,18), name = "Subject", 
                     labels = c("Child 1" , "Child 2","Child 3","Child 4", "Child 5")) +
  scale_color_manual(values =c("deeppink3","darkgreen","darkorchid","navy", "turquoise3"), name = "Subject", 
                     labels = c("Child 1" , "Child 2","Child 3","Child 4", "Child 5"))+
  labs(x = "PC 1 (25%)", y = "PC 2 (17%)")+
  theme_bw(base_size = 12)+ 
  theme(legend.position = "bottom", legend.title.align = 0.5)+
  guides(color=guide_legend(ncol=2, title.position = "top"))+
  guides(shape=guide_legend(ncol=2, title.position = "top"))

B2






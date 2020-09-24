##### R codes for data analysis and plotting in "Diversity of gut microbiomes in marine fishes is shaped by host-related factors" ####

##############################################################################################################
####Abundance against the prevalence of bacteria at the genus level in the fish guts among the 115 samples####
##############################################################################################################
design <- read.table("design_nounderline.csv", header=T, row.names= 1, sep=",",stringsAsFactors = F,check.names = FALSE) 
genus_ab<-read.table('Genus_table.txt',row.names=1,header=T, sep="\t",check.names=FALSE,stringsAsFactors = FALSE)
t_genus<-data.frame(t(genus_ab))
genus_ab<-data.frame(cbind(t_genus,'Uni ID'=design[match(rownames(t_genus),rownames(design)),1]),check.names = F) # add the unisample ID into the genus mother table
rownames(genus_ab)<-genus_ab$`Uni ID` #set the unisample name as the rownames
genus_ab<-genus_ab[,-512]
genus_ab<-data.frame(t(genus_ab),check.names = F,stringsAsFactors = F)

genus_ap<-genus_ab
genus_ap[genus_ap>0]<-1

genus_ap$sum<-rowSums(genus_ap) #calculate the sum of each OTU
genus_ap$percent<-genus_ap$sum/115 #calculate the percentage of each OTU

desc_stats = data.frame(Min=min(genus_ap$percent),#minimum
                        Max=max(genus_ap$percent),#maximum
                        Med=median(genus_ap$percent),#median
                        Mean=mean(genus_ap$percent),#mean
                        SD=sd(genus_ap$percent)#Standard deviation
)

desc_stats

nrow(genus_ap[genus_ap$percent>0.5,]) #calculate the number of the dominant genera
a<-genus_ap[genus_ap$percent>0.5,] #show the OTU existed in more than 30% samples
b<-genus_ab[rownames(genus_ab)[genus_ap$percent>0.5],] #show the OTU existed in more than 80% fish species and its relative abundance
c<-cbind(b,genus_ap$percent[match(rownames(b),rownames(genus_ap))]) # combine the relative abundance and a&p percentage together
colnames(c)[116]<-"prevalence"

#relative abundance of genus level
design <- read.table("design_nounderline.csv", header=T, row.names= 1, sep=",",stringsAsFactors = F,check.names = FALSE) 
genus_ab<-read.table('Genus_table.txt',row.names=1,header=T, sep="\t",check.names=FALSE,stringsAsFactors = FALSE)
t_genus<-data.frame(t(genus_ab),check.names = F)
genus_ab<-data.frame(cbind(t_genus,'Uni ID'=design[match(rownames(t_genus),rownames(design)),1]),check.names = F) # add the unisample ID into the genus mother table
rownames(genus_ab)<-genus_ab$`Uni ID` #set the unisample name as the rownames
genus_ab<-genus_ab[,-512]
genus_ab<-data.frame(t(genus_ab),check.names = F,stringsAsFactors = F)

mydescription<-function(x){desc_stats = data.frame(Min=min(x),#minimum
                                                   Max=max(x),#maximum
                                                   Med=median(x),#median
                                                   Mean=mean(x),#mean
                                                   SD=sd(x)); desc_stats}#Standard deviation

mydescription(rowMeans(genus_ab)) #calculate the basic statistical characteristics of average abundance

write.table(mydescription(rowMeans(genus_ab)), file = "Absence_presence.csv", append = TRUE, quote= FALSE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)

nrow(genus_ab[rowMeans(genus_ab)>0.01,])  #calculate the rest number of OTU
a<-genus_ab[rowMeans(genus_ab)>0.01,] 

#create a abundance and prevelance graphic
genus_mean_ab<-data.frame(rowMeans(genus_ab),check.names = F)#create a average abundance martix
ab_ap<-cbind(genus_mean_ab,genus_ap[match(rownames(genus_mean_ab),rownames(genus_ap)),117]) #combine the 
colnames(ab_ap)<-c("Abundance","Prevalence")
ab_ap$genus<-rownames(ab_ap)
nrow(ab_ap[ab_ap$Abundance>0.0125,])


library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggthemes)

windowsFonts(myFont = windowsFont("Times New Roman")) 

p<-ggplot(ab_ap,aes(x=Abundance,y=Prevalence)) + 
  geom_point()+
  geom_point(data=ab_ap[ab_ap$Abundance>0.0125,],aes(fill= genus),size=5,shape=21,colour='black') +
  geom_label_repel(data=ab_ap[ab_ap$Abundance>0.0125,],aes(label=genus),box.padding = 0.2, point.padding =0.3,size=4.5,segment.color='grey50') +
  theme_bw()+theme(legend.position="none")+ 
  theme(axis.title = element_text(size = 15, family="myFont",face = "bold"),axis.text=element_text(size=12,  family="myFont",face="bold", color = "black"))+
  labs(x = "Average abundance of genus in 115 samples", y= "Percentage of samples in which genus is present")
p
ggsave(paste( "genus-abundance-prevalence-2.png", sep=""), p, width = 8.5, height = 5.7, dpi=1200)  

# find the high-abundance genus in each sample
genus_ab[genus_ab>0.1]
cordinate=which(genus_ab>0.1,arr.ind=TRUE) #the cordinate information for the high-abundance genus
fre<-data.frame(table(rownames(cordinate)))
fre$per<-fre$Freq/115
colnames(fre)[1]<-'genus'
rownames(fre)<-fre$genus

genus_mean_ab<-data.frame(rowMeans(genus_ab),check.names = F)#create a average abundance martix
ab_ap_1<-cbind(fre,genus_ap[match(rownames(fre),rownames(genus_ap)),117]) #combine the abundance and prevalence matrix together
colnames(ab_ap_1)<-c('Genus','Frequency',"Abundance","Prevalence")
nrow(ab_ap_1[ab_ap_1$Abundance>0.05,])

library(ggplot2)
library(RColorBrewer)
library(ggrepel)

p<-ggplot(ab_ap_1,aes(x=Abundance,y=Prevalence)) + 
  geom_point()+
  geom_point(data=ab_ap_1[ab_ap_1$Abundance>0.05| ab_ap_1$Prevalence>0.75,],aes(fill= Genus),size=3.5,shape=21,colour='black') +
  geom_label_repel(data=ab_ap_1[ab_ap_1$Abundance>0.05 | ab_ap_1$Prevalence>0.75,],aes(label=Genus),box.padding = 0.2, point.padding =0.3,size=3.2,segment.color='grey50') +
  theme(legend.position="none")+
  labs(x = "Percentage of samples in which genus represents >10% of microbial community", y= "Percentage of samples in which genus is present") +
  xlim(0,0.20)
p
ggsave(paste( "genus-0.1abundance-prevalence.png", sep=""), p, width = 6, height = 4, dpi=1200)  

#add the genus with abundance lower than x%
genus_mean_ab<-data.frame(rowMeans(genus_ab),check.names = F)#create a average abundance martix
ab_ap_2<-cbind(genus_ap,fre[match(rownames(genus_ap),rownames(fre)),]) #combine the abundance and prevalence matrix together
ab_ap_2<-ab_ap_2[,117:120]
colnames(ab_ap_2)<-c("Prevalence",'Genus','Frequency',"Abundance")
ab_ap_2$Abundance[is.na(ab_ap_2$Abundance)] <- 0
nrow(ab_ap_2[ab_ap_2$Abundance>0.05,])

p<-ggplot(ab_ap_2,aes(x=Abundance,y=Prevalence)) + 
  geom_point(size=1)+
  geom_jitter(position=position_jitter(0.001)) +
  geom_point(data=ab_ap_2[ab_ap_2$Abundance>0.05,],aes(fill= Genus),size=3.5,shape=21,colour='black') +
  geom_label_repel(data=ab_ap_2[ab_ap_2$Abundance>0.05,],aes(label=Genus),box.padding = 0.2, point.padding =0.3,size=3.2,segment.color='grey50') +
  theme(legend.position="none")+
  xlim(-0.01,0.22)+
  labs(x = "Percentage of samples in which genus represents >10% of microbial community", y= "Percentage of samples in which genus is present")
p
ggsave(paste( "genus-0.05abundance-prevalence-1.png", sep=""), p, width = 9, height = 7, dpi=1200)  



###########################################################################################################################
####Box plots of the Shannon index of the gut microbiomes in the tested fishes among orders/feeding habit/trophic level####
###########################################################################################################################
library(plotly)
library(RColorBrewer)
library(processx)
design = read.table("design_nounderline.csv", header=T, row.names= 1, sep=",",check.names = F) 
alpha = read.table("alpha.txt", header=T, row.names= 1, sep="\t", check.names = F)

index = cbind(alpha, design[match(rownames(alpha), rownames(design)), ]) 

index$`Feeding habit 1`<-factor(index$`Feeding habit 1`, levels=c("Herbivore/Omnivore", "Zooplanktivore/Zoobenthivore", "Zoobenthivore", "Zoobenthivore/Piscivore", "Piscivore"))

#fish order
p = ggplot(index, aes(x=`Order`, y=shannon, fill=`Order`))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=0.8) +
  labs(x="Fish order", y="Shannon index") + 
  theme_calc()+
  theme(axis.text.y=element_text(hjust=0, size = 18,colour="black"),axis.text.x=element_text(colour="black",size = 11, vjust=1, hjust=1, angle = 30)) + 
  theme(axis.title.y = element_text(vjust=8, size=18), axis.title.x = element_text(size=14),legend.text = element_text(size=18),legend.title = element_text(size=18))+
  guides(fill=guide_legend(title='Order'))+
  scale_fill_manual(values = c( brewer.pal(5, "Set2")))
p
p<-ggplotly(p)
p

#fish feeding habit
index$`Feeding habit 1`<-factor(index$`Feeding habit 1`, levels=c("Herbivore/Omnivore", "Zooplanktivore/Zoobenthivore", "Zoobenthivore", "Zoobenthivore/Piscivore", "Piscivore"))
p = ggplot(index, aes(x=`Feeding habit 1`, y=shannon, fill=`Feeding habit 1`))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=1.2) +
  labs(x="Fish feeding habit", y="Shannon index") + 
  theme_calc()+
  theme(axis.text.y=element_text(hjust=0, size = 18,colour="black"),axis.text.x=element_text(colour="black",size = 11, vjust=1, hjust=1, angle = 30)) + 
  theme(axis.title.y = element_text(vjust=8, size=18), axis.title.x = element_text(size=14),legend.text = element_text(size=18),legend.title = element_text(size=18))+
  guides(fill=guide_legend(title='Feeding habit'))+
  scale_fill_manual(values = c( brewer.pal(5, "Set2")))
p
p<-ggplotly(p)
p


#Trophic level
p = ggplot(index, aes(x=`Feeding habit 3`, y=shannon, fill=`Feeding habit 3`))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=0.8) +
  labs(x="Trophic level group", y="Shannon index") + 
  theme_calc()+
  theme(axis.text.y=element_text(hjust=0, size = 18,colour="black"),axis.text.x=element_text(colour="black",size = 11, vjust=1, hjust=1, angle = 30)) + 
  theme(axis.title.y = element_text(vjust=8, size=18), axis.title.x = element_text(size=14),legend.text = element_text(size=18),legend.title = element_text(size=18))+
  guides(fill=guide_legend(title='Trophic level group'))+
  scale_fill_manual(values = c( brewer.pal(5, "Set2")))
p
p<-ggplotly(p)
p


##################################################################################################################
####Clustering pattern among all fish gut microbiome samples categorized by orders/feeding habit/trophic level####
##################################################################################################################
library("plotly")
library("vegan")

design = read.table("design_nounderline.csv", header=T, row.names= 1, sep=",",stringsAsFactors = F,check.names = F) 

otu_table<-read.table("otu-table.tsv",header=T,sep="\t", row.names= 1, check.names=FALSE,stringsAsFactors = F) # import the OTU data

variability_table = function(cca){
  chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table = cbind(chi, chi/chi[1])
  colnames(variability_table) = c("inertia", "proportion")
  rownames(variability_table) = c("total", "constrained", "unconstrained")
  return(variability_table)
}

# Constrained analysis OTU table by fish orders
capscale.gen = capscale(t(otu_table) ~ `Order` , data=design, add=F, sqrt.dist=T, distance="bray") #Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn

perm_anova.gen = anova.cca(capscale.gen)# ANOVA-like permutation analysis

var_tbl.gen = variability_table(capscale.gen)# generate variability tables and calculate confidence intervals for the variance
eig = capscale.gen$CCA$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]

points = capscale.gen$CCA$wa[, 1:3]# extract the weighted average (sample) scores
points = as.data.frame(points)
colnames(points) = c("x", "y","z")
points = cbind(points, design[match(rownames(points), rownames(design)),])

p <- plot_ly(points, x = ~x, y = ~y, z = ~z, color = ~`Order`)%>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")),
                      yaxis = list(title = paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")),
                      zaxis = list(title = paste("CPCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""))))
p

# Constrained analysis OTU table by feeding habits
capscale.gen = capscale(t(otu_table) ~ `Feeding habit 1` , data=design, add=F, sqrt.dist=T, distance="bray") #Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn

perm_anova.gen = anova.cca(capscale.gen)# ANOVA-like permutation analysis

var_tbl.gen = variability_table(capscale.gen)# generate variability tables and calculate confidence intervals for the variance
eig = capscale.gen$CCA$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]

points = capscale.gen$CCA$wa[, 1:3]# extract the weighted average (sample) scores
points = as.data.frame(points)
colnames(points) = c("x", "y","z")
points = cbind(points, design[match(rownames(points), rownames(design)),])

p <- plot_ly(points, x = ~x, y = ~y, z = ~z, color = ~`Feeding habit 1`)%>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")),
                      yaxis = list(title = paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")),
                      zaxis = list(title = paste("CPCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""))))
p

# Constrained analysis OTU table by trophic level
capscale.gen = capscale(t(otu_table) ~ `Feeding habit 3` , data=design, add=F, sqrt.dist=T, distance="bray") #Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn


perm_anova.gen = anova.cca(capscale.gen)# ANOVA-like permutation analysis


var_tbl.gen = variability_table(capscale.gen)# generate variability tables and calculate confidence intervals for the variance
eig = capscale.gen$CCA$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]


points = capscale.gen$CCA$wa[, 1:3]# extract the weighted average (sample) scores
points = as.data.frame(points)
colnames(points) = c("x", "y","z")
points = cbind(points, design[match(rownames(points), rownames(design)),])

p <- plot_ly(points, x = ~x, y = ~y, z = ~z, color = ~`Feeding habit 3`)%>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep="")),
                      yaxis = list(title = paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")),
                      zaxis = list(title = paste("CPCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""))))
p




##################################################################################################################
####Discriminative gut microbiomes detected among the three distinguishing fish feeding habits/trophic level######
##################################################################################################################
library(labdsv)
library(pheatmap)

#Feeding habit
design <- read.table("design_nounderline.csv", header=T, row.names= 1, sep=",", check.names = F) 
design$`Feeding habit 1` = factor(design$`Feeding habit 1`,levels=c("Herbivore/Omnivore","Zooplanktivore/Zoobenthivore","Zoobenthivore","Zoobenthivore/Piscivore","Piscivore"))

otu_abu<-read.table("otu_abundance.csv", sep=",", header=T, check.names=F,stringsAsFactors = F)#otu abundance table
otu_abu<-data.frame(t(otu_abu), check.names = F)
otu_abu<-as.matrix(otu_abu)
class(otu_abu)<-"numeric"
otu_abu<-otu_abu[rownames(design),]


clust<-design$`Feeding habit 1`
otu_ind<-indval(otu_abu,clust)#calculate the indicator value for feeding habit 1
otu_pval_ind<-cbind(otu_ind$indval,otu_ind$pval)
otu_pval_ind$otuID<-rownames(otu_pval_ind)
colnames(otu_pval_ind)[6]<-"otu_pval"

sel_indotu<-list()
selpval_indotu<-subset(otu_pval_ind,otu_pval<=0.01)

Order<-c("Herbivore/Omnivore","Zooplanktivore/Zoobenthivore","Zoobenthivore","Zoobenthivore/Piscivore","Piscivore")

#select the indicator value for each feeding habit,but "Zoobenthivore" and "Zoobenthivore/Piscivore", 
#didn't get the desired indicator value
selpval_indotu<-selpval_indotu[order(selpval_indotu[,1],decreasing = T),]
sel_indotu_tab<-selpval_indotu[1:15,] #indval>0.28
sel_indotu_tab$ind_cat<-Order[1]
sel_indotu[[Order[1]]]<-sel_indotu_tab

selpval_indotu<-selpval_indotu[order(selpval_indotu[,2],decreasing = T),]
sel_indotu_tab<-selpval_indotu[1:12,] #indval>0.4
sel_indotu_tab$ind_cat<-Order[2]
sel_indotu[[Order[2]]]<-sel_indotu_tab

selpval_indotu<-selpval_indotu[order(selpval_indotu[,5],decreasing = T),]
sel_indotu_tab<-selpval_indotu[1:8,] #indval>0.2
sel_indotu_tab$ind_cat<-Order[5]
sel_indotu[[Order[5]]]<-sel_indotu_tab

sel_indotu_inte <- rbind(sel_indotu[[1]],sel_indotu[[2]],sel_indotu[[3]])

otu_ind_summary<-Reduce(function(...) merge(..., all=T), sel_indotu)
otu_ind_summary$ind_cat<-factor(otu_ind_summary$ind_cat,levels=c("Herbivore/Omnivore","Zooplanktivore/Zoobenthivore","Piscivore"))
otu_ind_summary<-otu_ind_summary[order(otu_ind_summary$ind_cat),]

otu_ind_summary$otuID<-c("e06ff231e225ef1500678035eb3d7655","0a8bbc088b8beec6326c627280b5c0d6","3b1b336cbc81cab353492b4a9955714a",
                         "ba7c1080b4d12da4b8dfe57d480afdf2","04e08696178132addcf37990c127c0d4","16e0520799028e3de103aee82db00d7b",
                         "d110483ba324dac60f260aa8936bfd29","f8a02acc65025ef66f7f1f300b213ede","32d4bbd89c257433688c0d8945fe2366",
                         "54fe3b4d57b7521c32d691759f1b7cf8","69c1dd230e378d44a6833e0ec830527c","c75cf96d0bfc2bfe663b9cedbdda8148",
                         "9f20ec68881e1158cf1eb8c38e6a2926","c6c07cb4cb499e835776969ec52ad82b","d089b08a255484c51d94d5994450e4e3",
                         "eca4da5cf7004dc372f7965a177aa605","885ab07d0f8efa0791b6a6f1eecf3645","490ca131b376b701aef9dd4cdf08faea",
                         "61319d924e9ac98c34564d5a73bd0a19","54cab533f66bd39b30f32ebcda532f08","85380253842e383168f6ebceb4de4f7f",
                         "8f6b02986fd2ce505a3f8714ae7c5f04","1bdaaddfd7436609f8327fe211b52ded","7e6c02c5745a35523da1e73a930b744f",
                         "a1a102064b559c55cf9dc1f50cc36ed2","2df2bb9ba63e1111154657e968941efd","1409fdaeed3befd3064189b22b0b9f97",
                         "54d82fe61718250783d69e6d470674d5","23a464c3e8092ad3acc51493f99ccd68","ac44c690eb5612176b4820aa954e4d11",
                         "158c744a6c24032c5f2a2d6b28f9953f","fb6f934214e1e514652f68f2e6f73d62","9cd9b540c25187cbe9cc44ec80ae0b43",
                         "10fd9c602733f5a2c4cdf1b9e350564b","7f38895c75a4bbdc460259f8c06405b7") #delete the character X in the begining the otu ID
plotdata_ind<-otu_abu[,otu_ind_summary$otuID]
index<-cbind(plotdata_ind,design[match(rownames(plotdata_ind),rownames(design)),])

breaksList = seq(log10(0.000001), log10(0.094923), by = 0.001) #to unify the legend bar range of three heatmap of three orders
png(file=paste("indicator_HVOV.png"), width = 20, height = 8, units="cm",res = 600)
pheatmap(log10(plotdata_HVOV+0.000001), color = colorRampPalette(c("white","#016392"))(50), cluster_rows=FALSE, cluster_cols=FALSE, fontsize=5, fontsize_row=5, cellwidth=10, 
         cellheight=5, legend=TRUE,border_color = NA, show_rownames = F,show_colnames = T)
dev.off()

png(file=paste("indicator_ZPZB.png"), width = 20, height = 8, units="cm",res = 600)
pheatmap(log10(plotdata_ZPZB+0.000001), color = colorRampPalette(c("white","#016392"))(50), cluster_rows=FALSE, cluster_cols=FALSE, fontsize=5, fontsize_row=5, cellwidth=10, 
         cellheight=5, legend=TRUE,border_color = NA, show_rownames = F,show_colnames = T)
dev.off()

png(file=paste("indicator_PV.png"), width = 20, height = 8, units="cm",res = 600)
pheatmap(log10(plotdata_PV+0.000001), color = colorRampPalette(c("white","#016392"))(50), cluster_rows=FALSE, cluster_cols=FALSE, fontsize=5, fontsize_row=5, cellwidth=10, 
         cellheight=5, legend=TRUE,border_color = NA, show_rownames = F,show_colnames = T)
dev.off()



#############################################################################
####The core gut microbiomes in the fish gut of different feeding habit######
#############################################################################
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggthemes)

#Absence and presence of ASV in each feeding habit
core_otu_fh_tax<-read.table("feeding_habit_abu_pre_tax.csv",row.names= 1,  header=T, sep=",",check.names=FALSE)

HVOV_plot <- core_otu_fh_tax[core_otu_fh_tax$`Herbivore/Omnivore`>0.015 | core_otu_fh_tax$HOOV_pre>0.4,] #5 dominant otu
HVOV_plot$label <- HVOV_plot$Family
ZPZB_plot <- core_otu_fh_tax[core_otu_fh_tax$`Zooplanktivore/Zoobenthivore`>0.008 & core_otu_fh_tax$ZPZB_pre>0.6,] #3 dominant otu
ZPZB_plot$label <- c("Ralstonia", "ML635J-21", "Gaiellales")
ZB_plot <- core_otu_fh_tax[core_otu_fh_tax$`Zoobenthivore`>0.007 & core_otu_fh_tax$ZB_pre>0.5,] #4 dominant otu
ZB_plot$label <- c("Ralstonia", "ML635J-21", "Ralstonia", "Ralstonia")
ZBPV_plot <- core_otu_fh_tax[core_otu_fh_tax$`Zoobenthivore/Piscivore`>0.011 | core_otu_fh_tax$ZBPV_pre>0.4,] #9 dominant otu
ZBPV_plot$label <- c("Photobacterium", "ML635J-21", "Photobacterium", "Clostridium perfringens","Photobacterium","Photobacterium","Photobacterium","Photobacterium","Ralstonia")
PV_plot <- core_otu_fh_tax[core_otu_fh_tax$`Piscivore`>0.018 | core_otu_fh_tax$PV_pre>0.4,] #13 dominant otu
PV_plot$label <- c ("Clostridium perfringens", "Cetobacterium", "Photobacterium", "Lactococcus", "Photobacterium", "Acinetobacter  guillouiae",
                    "Acinetobacter  guillouiae","Photobacterium","Clostridium perfringens","Enhydrobacter","Photobacterium","Clostridium perfringens", "Cetobacterium")

p<-ggplot(core_otu_fh_tax,aes(x=`Herbivore/Omnivore`,y=`HOOV_pre`)) + 
  geom_point(alpha=0.5)+ theme_calc()+
  geom_point(data=HVOV_plot,aes(fill= `label`),size=3,shape=21,colour='black') +
  theme(legend.position="none") +
  labs(x = "Average abundance of OTU in Herbivore/Omnivore samples", y= "Percentage of samples in which OTU is present") +
  geom_label_repel(data=HVOV_plot,aes(label=label),box.padding = 0.2, point.padding =0.3,size=3.2,segment.color='grey50')
p


p<-ggplot(core_otu_fh_tax,aes(x=`Zooplanktivore/Zoobenthivore`,y=`ZPZB_pre`)) + 
  geom_point(alpha=0.5)+ theme_calc()+
  geom_point(data=ZPZB_plot,aes(fill=label),size=3,shape=21,colour='black') +
  theme(legend.position="none") +
  labs(x = "Average abundance of OTU in Zooplanktivore/Zoobenthivore samples", y= "Percentage of samples in which OTU is present") +
  geom_label_repel(data=ZPZB_plot,aes(label=label),box.padding = 0.2, point.padding =0.3,size=3.2,segment.color='grey50') 
p


p<-ggplot(core_otu_fh_tax,aes(x=`Zoobenthivore`,y=`ZB_pre`)) + 
  geom_point(alpha=0.5)+ theme_calc()+
  geom_point(data=ZB_plot,aes(fill= `label`),size=3,shape=21,colour='black')+
  theme(legend.position="none") +
  labs(x = "Average abundance of OTU in Zoobenthivore samples", y= "Percentage of samples in which OTU is present") +
  geom_label_repel(data=ZB_plot,aes(label=label),box.padding = 0.2, point.padding =0.3,size=3.2,segment.color='grey50')
p


p<-ggplot(core_otu_fh_tax,aes(x=`Zoobenthivore/Piscivore`,y=`ZBPV_pre`)) + 
  geom_point(alpha=0.5)+ theme_calc()+
  geom_point(data=ZBPV_plot,aes(fill= label),size=3,shape=21,colour='black') +
  theme(legend.position="none") +
  labs(x = "Average abundance of OTU in Zoobenthivore/Piscivore samples", y= "Percentage of samples in which OTU is present") +
  geom_label_repel(data=ZBPV_plot,aes(label=label),box.padding = 0.2, point.padding =0.3,size=3.2,segment.color='grey50')
p


p<-ggplot(core_otu_fh_tax,aes(x=`Piscivore`,y=`PV_pre`)) + 
  geom_point(alpha=0.5)+ theme_calc()+
  geom_point(data=PV_plot,aes(fill= label),size=3,shape=21,colour='black') +
  theme(legend.position="none") +
  labs(x = "Average abundance of OTU in Piscivore samples", y= "Percentage of samples in which OTU is present") +
  geom_label_repel(data=PV_plot,aes(label=label),box.padding = 0.2, point.padding =0.3,size=3.2,segment.color='grey50')
p

#barchat
HVOV_plot$`Feeding habit abbre` <-"HVOV"
ZPZB_plot$`Feeding habit abbre` <-"ZPZB"
ZB_plot$`Feeding habit abbre` <-"ZB"
ZBPV_plot$`Feeding habit abbre` <-"ZBPV"
PV_plot$`Feeding habit abbre` <-"PV"
HVOV_plot$`OTU` <-rownames(HVOV_plot)
ZPZB_plot$`OTU` <-rownames(ZPZB_plot)
ZB_plot$`OTU` <-rownames(ZB_plot)
ZBPV_plot$`OTU` <-rownames(ZBPV_plot)
PV_plot$`OTU` <-rownames(PV_plot)
core_set<-rbind(HVOV_plot,ZPZB_plot,ZB_plot,ZBPV_plot,PV_plot)

library(reshape2)
library(dplyr)
core_set_bar<-melt(core_set[,-6:-10], id.vars = colnames(core_set[,-6:-10])[6:15], variable.name="Feeding habit", value.name = "Relative abundance")
core_set_bar$`Feeding habit`<-as.factor(core_set_bar$`Feeding habit`)
core_set_bar$`Feeding habit abbre`<-factor(core_set_bar$`Feeding habit abbre`,levels=c("HVOV", "ZPZB", "ZB", "ZBPV", "PV"))

OTU_level<-c("d089b08a255484c51d94d5994450e4e3","c75cf96d0bfc2bfe663b9cedbdda8148","2603d8b6951bc0b07fcf06fa7cf22c1b","c6c07cb4cb499e835776969ec52ad82b","9f20ec68881e1158cf1eb8c38e6a2926",
             "e72c95ec6a14134c0df9918d101bd5f0","1409fdaeed3befd3064189b22b0b9f97","d94d61c10c87c8fc5aec8e93fa7bff26",
             "886c51bb637bf84c71a567504746ea51","25644ca7604d1ccbc8f063e9ef2e92fd","1a1bc2aa20394f40d6f578c2b7103cc0","4e8f19d74e2e8ceeefba6e53b4656834",
             "e7fe5bb7a0858494ee3179c6b19b9e49","21cedbbaa6970d8a46382aa6a1e1bec8","6a20acc91daf66a1283fb7fa904ab35e","2cf1bc37656d223762a6a6a7b9eaa482","430fbb3ef20a4ed10cffaf6599ce36a9","066f77b196efbb4237ce44598429fb46",
             "4fc20797dc6c89d073b8a541864a370b","6bdd0a8d021e42f812eed1ee42882254","13791efc306177bf7940f102598a86a0","f5f81981599a6a98c3291954a79a1532","981d7e116af843d490bbe8b53d0f9856","6c94e0d993b850e9cb80e30e238db862","7087e67cb01856a8b5e588f912bc06fc","ace484772ca4436ac0b4450c33c19c1d","c17edb7950b5ccdc32a5303ebf183a1e","fc04efc3c9402be10973a40e4e6a9b9b","aa12cd3df9a3136e19c655b7e6287c6d") 
# 066f77b196efbb4237ce44598429fb46 430fbb3ef20a4ed10cffaf6599ce36a9 (PV;ZB/PV)  1a1bc2aa20394f40d6f578c2b7103cc0  (ZB;ZB/PV)  
# d94d61c10c87c8fc5aec8e93fa7bff26 (ZP/ZB ZB ZB/PV) 


core_set_bar$OTU<-factor(core_set_bar$OTU,levels=OTU_level)
core_set_bar<-na.omit(core_set_bar) #another way: core_set_bar[complete.cases(core_set_bar[,9]),]
core_set_bar<-arrange(core_set_bar,`Feeding habit abbre`, `Relative abundance`)
core_set_bar$`Feeding habit`<-factor(core_set_bar$`Feeding habit`, levels=c("Herbivore/Omnivore", "Zooplanktivore/Zoobenthivore", "Zoobenthivore", "Zoobenthivore/Piscivore","Piscivore"))

p<-ggplot(core_set_bar, aes(fill=`Feeding habit`, y=`Relative abundance`, x=OTU)) + 
  geom_bar(position="dodge", stat="identity", width = 0.8)+ theme_calc() + 
  labs(x = "OTU", y = "Relative abundance", fill = "Feeding habit") + #edit the label
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
  scale_fill_brewer(palette = "Set1") +
  theme(plot.margin = unit(c(0.2,0.2,0.2,2.5), "cm"))



####################################################################################################
####Percentages of unassigned sequences for fish gut microbiomes at different taxonomic levels######
####################################################################################################
library(ggplot2)
library(ggrepel)

index <- read.table("unassigned_list.csv", header=T,  sep=",",stringsAsFactors = FALSE,check.names = F) 

p <- ggplot(index, aes(x=Taxonomy, y=Unassigned_sequence_percentage, group= `Fish Species`)) + 
  geom_line(aes(color=`Fish Species`),size=0.5) + 
  geom_point(size=2.8,shape=21,aes(color=`Fish Species`,fill=`Fish Species`)) +
  guides(col = guide_legend(ncol = 2)) +
  ylab("Percentage of Unassigned Sequences (%)") + 
  xlim("Domain", "Phylum", "Class",  "Order",  "Family" ,"Genus", " ") +
  theme_bw()+
  theme(legend.position=c(0.3,0.7),legend.key.size = unit(0.2, "in"),legend.text=element_text(face="italic"))+
  geom_label_repel(data=subset(index,Taxonomy=="Genus"),aes(label=`Species_abbre2`),fontface="bold.italic", show.legend=F, hjust = 0,nudge_x= 0.5,direction = "y",segment.size = 0.4, point.padding =0.3,size=4.3, color="black",segment.color='grey50')
p


#############################################################################################
####Abundance against the prevalence of bacteria at the ASV level among the 115 samples######
#############################################################################################

#Absence and presence of OTU of all 115 samples
otu<-read.table("otu-table.tsv",row.names= 1,   sep="\t", header= T, check.names=FALSE,stringsAsFactors = FALSE)
otu_ab<-otu
otu_ab[otu_ab>0]<-1

otu_ab$sum<-rowSums(otu_ab) #calculate the sum of each OTU
otu_ab$percent<-otu_ab$sum/115 #calculate the percentage of each OTU

desc_stats = data.frame(Min=min(otu_ab$percent),#minimum
                        Max=max(otu_ab$percent),#maximum
                        Med=median(otu_ab$percent),#median
                        Mean=mean(otu_ab$percent),#mean
                        SD=sd(otu_ab$percent)#Standard deviation
)

desc_stats

otu_RA<-read.table("otu_abundance.csv",row.names= 1, header=T, sep=",",check.names=FALSE,stringsAsFactors = FALSE)

nrow(otu_ab[otu_ab$percent>0.3,]) #calculate the number of the dominant OTU
a<-otu_ab[otu_ab$percent>0.3,] #show the OTU existed in more than 30% samples
b<-otu_RA[rownames(otu_ab)[otu_ab$percent>0.3],] #show the OTU existed in more than 80% fish species and its relative abundance
c<-cbind(b,otu_ab$percent[match(rownames(b),rownames(otu_ab))]) # combine the relative abundance and a&p percentage together
colnames(c)[116]<-"presence"


#create the relative abundance table of all 115 samples
otu_abundance_115_1<-otu_abundance_115
otu_abundance_115_1[otu_abundance_115_1<=0.01]<-0 #set the otu = 0 if the abundance lower than 0.01
otu_abundance_115_2<-otu_abundance_115_1[rowSums(otu_abundance_115_1)!=0,] #create a new mother table with OTU abundance higher than 0.01
nrow(otu_abundance_115_2)  #calculate the rest number of OTU
colSums(otu_abundance_115_2) #calculate the filtered OTU total abundance

which(otu_abundance_115>=0.02,arr.ind=T) # show the row and col of the OTUs with abundance higher than 0.1

#create a abundance and prevelance graphic for otu 

otu_mean_ab<-data.frame(rowMeans(otu_RA),check.names = F)#create a average abundance martix
otu_ab_ap<-cbind(otu_mean_ab,otu_ab[match(rownames(otu_mean_ab),rownames(otu_ab)),117]) #combine the 
colnames(otu_ab_ap)<-c("Abundance","Prevalence")
otu_ab_ap$genus<-rownames(otu_ab_ap)
nrow(otu_ab_ap[otu_ab_ap$Prevalence>0.36,])


library(ggplot2)
library(RColorBrewer)
library(ggrepel)

p<-ggplot(otu_ab_ap,aes(x=Abundance,y=Prevalence)) + 
  geom_point()+
  geom_point(data=otu_ab_ap[otu_ab_ap$Prevalence>0.36,],aes(fill= genus),size=3,shape=21,colour='black') +
  labs(x = "Average abundance of OTU in 115 samples", y= "Percentage of samples in which OTU is present")+
  theme_bw()+
  theme(legend.position="none")
p


##################################################################
####Box plots of the alpha diversity (Shannon index) of ASVs######
##################################################################
library("ggplot2") # load related packages
library("ggthemes") 
library('dplyr')
library("plotly")
library("RColorBrewer")


design = read.table("design_nounderline.csv", header=T, row.names= 1, sep=",", check.names = F, stringsAsFactors = F) 
alpha = read.table("alpha.txt", header=T, row.names= 1, sep="\t", check.names = F, stringsAsFactors = F)
index = cbind(alpha, design[match(rownames(alpha), rownames(design)), ]) 
#shannon index _ ascending order
mean_otu<-tapply(index$shannon, index[,c('Species')], mean)
mean_otu<-data.frame(mean_otu)
mean_otu$Species<-rownames(mean_otu)
mean_otu$Species[1]<-"Acanthopagrus schlegelii schlegelii"
mean_otu$rank<-rank(mean_otu$mean_otu)
mean_otu<-arrange(mean_otu,rank)
index$`Fish Species`<-factor(index$`Fish Species`,levels=mean_otu$`Species`)

index$`Feeding habit 1`<-factor(index$`Feeding habit 1`, levels=c("Herbivore/Omnivore", "Zooplanktivore/Zoobenthivore", "Zoobenthivore", "Zoobenthivore/Piscivore", "Piscivore"))
p = ggplot(index, aes(x=`Fish Species`, y=shannon, fill=`Feeding habit 1`))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.5) +
  labs(x="Fish species", y="Shannon index") + 
  theme_calc()+
  theme(axis.text.y=element_text(hjust=0, size = 18,colour="black"),axis.text.x=element_text(colour="black",size = 11, vjust=1, hjust=1, angle = 30)) + 
  theme(axis.title.y = element_text(vjust=8, size=18), axis.title.x = element_text(size=14),legend.text = element_text(size=18),legend.title = element_text(size=18))+
  guides(fill=guide_legend(title='Feeding habit'))+
  scale_fill_manual(values = c( brewer.pal(5, "Set2")))
p
p<-ggplotly(p)
p


##############################################################################################
####Mean stable isotope C and N  values identified in the dorsal muscle tissue of fishes######
##############################################################################################
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggthemes)

design = read.table("design_nounderline.csv", header=T, row.names= 1, sep=",",stringsAsFactors = FALSE,check.names = F) 
isotope <- aggregate(design[,8:14], by=list('Fish Species'=design$`Fish Species`),mean)
for (i in (1:nrow(isotope))){
  isotope$`Feeding habit 1`[i] = design$`Feeding habit 1`[match(isotope$'Fish Species'[i],design$'Fish Species')]
  isotope$Order[i]=design$Order[match(isotope$'Fish Species'[i],design$'Fish Species')]
  isotope$Suborder[i]=design$Suborder[match(isotope$'Fish Species'[i],design$'Fish Species')]
} #add the Diet information to the isotope

p <- ggplot(isotope,aes(x=d13C,y=d15N)) + 
  geom_label_repel(aes(label=`Fish Species`),box.padding = 0.2, point.padding =0.3,size=5,segment.color='grey50') +
  geom_point(aes(size =`Trophic level`,color=`Trophic level`)) + 
  guides(color=guide_legend(title='Feeding habit'),size=guide_legend(title='Trophic level'))+
  xlim(-21.1,-15.5) +ylim(11,15.2)+theme_bw()+
  theme(axis.title=element_text(size=15),axis.text = element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15))
p


##################################
####Heatmap of KEGG pathways######
##################################
pathway_level<- read.table(file = "pathway_level.txt",header=TRUE, sep="\t", check.names = F,quote="")#merge the pathway level with test_ko_L3
metagenome_L1 <- read.table(file="metagenome_L1.txt",header=TRUE, sep="\t", check.names = F,quote="")
metagenome_L2 <- read.table(file="metagenome_L2.txt",header=TRUE, sep="\t", check.names = F,quote="")
metagenome_L3 <- read.table(file="metagenome_L3.txt",header=TRUE, sep="\t", check.names = F,quote="")

design <- read.table("design_nounderline.csv", header=T, row.names=1, sep=",", stringsAsFactors = FALSE,check.names = F) 

#heatmap for LEVEL 1
L1 <- data.frame(t(metagenome_L1),check.names = F)
L1$sampleID<-rownames(L1)
L1 <- cbind(L1,design[match(rownames(L1),rownames(design)),])
library(pheatmap)

#calculate the mean abundance of each fish species
L1_species = aggregate(L1[, 1:7], list(L1$Species), mean)
rownames(L1_species)<-L1_species$Group.1
L1_species<-L1_species[,-1]
#log transformation
L1_sp <- L1_species[,c('Metabolism','Genetic Information Processing', 'Environmental Information Processing','Cellular Processes','Human Diseases','Organismal Systems')] #reorder the order of different pathways according to its abundances
#draw the heatmap
pheatmap(t(L1_sp), cellwidth=20, cellheight=15, cluster_cols =1, cluster_rows = 1, fontsize_row = 10, fontsize_col = 10, angle_col = 45, filename = "L1_sp_clust.png", width = 10, height = 6)

#heatmap for LEVEL 2
L2 <- data.frame(t(metagenome_L2),check.names = F)
L2$sampleID<-rownames(L2)
L2 <- cbind(L2,design[match(rownames(L2),rownames(design)),])
library(pheatmap)

#calculate the mean abundance of each fish species
L2_species = aggregate(L2[, 1:40], list(L2$Species), mean)
rownames(L2_species)<-L2_species$Group.1
L2_species<-L2_species[,-1]
#log transformation
L2_sp <- L2_species[,-c(6,35)]
#draw the heatmap
pheatmap(t(L2_sp),cellwidth=20, cellheight=15, cluster_cols =0, cluster_rows = 0, fontsize_row = 10, fontsize_col = 10, angle_col = 45, filename = "L2_sp_log.png", width = 12, height = 11)
pheatmap(t(L2_sp),cellwidth=20, cellheight=15, cluster_cols =1, cluster_rows = 1, fontsize_row = 10, fontsize_col = 10, angle_col = 45, filename = "L2_sp_clust.png",width = 12, height = 11)

#heatmap for digestion related gene
design <- read.table("design_nounderline.csv", header=T, row.names=1, sep=",", stringsAsFactors = FALSE,check.names = FALSE) 

digestion_L3 <- read.table("digestion_L3.txt", row.names=1, header=T, sep="\t", quote="", check.names = F, fill=TRUE) 
digestion_L3<-data.frame(t(data.frame(t(digestion_L3),check.names = F)/rowSums(t(digestion_L3))),check.names = F)#relative abundance

digestion_L3 <- data.frame(t(digestion_L3),check.names = FALSE)
digestion_L3$sampleID <- rownames(digestion_L3)
digestion_L3 <- cbind(digestion_L3,design[match(rownames(digestion_L3),rownames(design)),])

digestion_L3_fh = aggregate(digestion_L3[, 1:40], by=list(digestion_L3$`Feeding habit 1`), mean)
rownames(digestion_L3_fh)<-digestion_L3_fh$Group.1
digestion_L3_fh<-digestion_L3_fh[,-1]
digestion_L3_fh <- data.frame(t(digestion_L3_fh),check.names = F)

library(pheatmap)
digestion_L3_fh <- digestion_L3_fh[,c(1,5,3,4,2)]

pheatmap(digestion_L3_fh, cellwidth=28, cellheight=11, cluster_cols =0, cluster_rows = 0, fontsize_row = 10, fontsize_col = 10, angle_col = 45, filename = "di_L3_fh.png", width = 10, height = 8)


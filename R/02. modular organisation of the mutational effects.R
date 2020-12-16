require(ggplot2)
require(dplyr)
require(reshape2)


dir.create("Modular organisation")
path="Modular organisation"

#required data:
load("nscore_df.RData")
load("yeast_controls.RData")


ABseq=c("D","A","E","F","R","H","D","S","G","Y","E","V","H","H","Q","K","L","V","F","F","A","E","D","V","G","S","N","K","G","A","I","I","G","L","M","V","G","G","V","V","I","A")
ABseq_pos=c("D\n1","A\n2","E\n3","F\n4","R\n5","H\n6","D\n7","S\n8","G\n9","Y\n10","E\n11","V\n12","H\n13","H\n14","Q\n15","K\n16","L\n17","V\n18","F\n19","F\n20","A\n21","E\n22","D\n23","V\n24","G\n25","S\n26","N\n27","K\n28","G\n29","A\n30","I\n31","I\n32","G\n33","L\n34","M\n35","V\n36","G\n37","G\n38","V\n39","V\n40","I\n41","A\n42")
vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "*")
disease_mutations<- c("H-6-R", "D-7-N","D-7-H", "E-11-K","K-16-N", "A-21-G", 
                      "E-22-Q", "E-22-K", "E-22-G", "D-23-N", "L-34-V", "A-42-T")
color_axis=c("D1"="#EE0011FF","A2"="darkgrey","E3"="#EE0011FF","F4"="#9A703EFF","R5"="#0C5BB0FF","H6"="#15983DFF",
             "D7"="#EE0011FF","S8"="#15983DFF","G9"="#FA6B09FF","Y10"="#9A703EFF","E11"="#EE0011FF","V12"="darkgrey",
             "H13"="#15983DFF","H14"="#15983DFF","Q15"="#15983DFF","K16"="#0C5BB0FF","L17"="darkgrey","V18"="darkgrey",
             "F19"="#9A703EFF","F20"="#9A703EFF","A21"="darkgrey","E22"="#EE0011FF","D23"="#EE0011FF","V24"="darkgrey",
             "G25"="#FA6B09FF","S26"="#15983DFF","N27"="#15983DFF","K28"="#0C5BB0FF","G29"="#FA6B09FF","A30"="darkgrey",
             "I31"="darkgrey","I32"="darkgrey","G33"="#FA6B09FF","L34"="darkgrey","M35"="darkgrey","V36"="darkgrey",
             "G37"="#FA6B09FF","G38"="#FA6B09FF","V39"="darkgrey","V40"="darkgrey","I41"="darkgrey","A42"="darkgrey")
color_axis_y<-c("black","#FEC10BFF", rep("#15983DFF",6), rep("#EE0011FF", 2), rep("#0C5BB0FF", 2),rep("#9A703EFF", 3),  rep("darkgrey", 5), "#FA6B09FF" )






#####  heatmap (contains stops)

#add syn
positions<-(1:42)
syn.df<-data.frame(
  "WT_AA"= ABseq,
  "Mut"= ABseq,
  "Pos"= positions, 
  "sigma"=0, 
  "nscore_c"=0,
  "ID"="syn"
)

heatmap_df<-rbind(singles_stops[,c("WT_AA", "Mut", "Pos", "sigma", "nscore_c", "ID")], syn.df)


#add info fAD
heatmap_df$box<-"VUS"
heatmap_df[heatmap_df$ID %in% disease_mutations,]$box<-"Dominant"
heatmap_df[heatmap_df$ID=="A-2-V",]$box<-"Recessive"

#add info syn
heatmap_df$label<-""
heatmap_df[heatmap_df$ID=="syn",]$label<-"*"



p_heatmap<-ggplot(heatmap_df)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)),fill=nscore_c), size=0.1)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)),color=box), fill=NA, size=1)+
  theme_minimal()+
  theme()+
  scale_x_continuous(breaks=seq(1:42), labels = ABseq_pos, expand = c(0,0))+
  labs(x="AB(1-42) WT amino acid and position", y="Mutant amino acid", fill="Nucleation score")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x.top = element_line(),       
        plot.title =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=15, face = "bold"),
        legend.text = element_text(size=15), 
        axis.text.x = element_text(color=color_axis, size=12),
        axis.text.y = element_text(color=color_axis_y, size=12),
        axis.title = element_text(size = 20))+
  geom_text(aes(Pos,factor(Mut, levels=rev(vectorAA)),label=label), size=8)+
  
  scale_fill_gradientn(colours=c("darkorange2", "#f5ad66", "lightgrey", "darkblue"), breaks=c(-4, -2, 0, 2),
                       na.value = "white") +
  scale_color_manual(values=c("#EE0011FF","#9A703EFF", NA), 
                     labels=c("Dominant",  "Recessive", ""), name=""  )
p_heatmap

ggsave(p_heatmap,path=path, file="p_heatmap.pdf",width=14, height=8)







####  heatmap FDR

heatmap_fdr<-singles_stops[,c("Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10")]

heatmap_fdr$category<-"WT-like"

heatmap_fdr[(heatmap_fdr$p.adjust<0.25 & heatmap_fdr$nscore_c<0),]$category<- "NS- 25%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.1 & heatmap_fdr$nscore_c<0),]$category<- "NS- 10%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.05 & heatmap_fdr$nscore_c<0),]$category<- "NS- 5%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.01 & heatmap_fdr$nscore_c<0),]$category<- "NS- 1%"

heatmap_fdr[(heatmap_fdr$p.adjust<0.25 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 25%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.1 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 10%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.05 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 5%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.01 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 1%"


#add syn

syn.df<-data.frame(
  "WT_AA"= ABseq,
  "Mut"= ABseq,
  "Pos"= c(1:42), 
  "nscore_c"=0,
  "ID"="syn",
  "p.adjust"=NA,
  "category_10"="WT-like",
  "category"="WT-like"
)


heatmap_fdr<-rbind(heatmap_fdr, syn.df)

#add info syn
heatmap_fdr$label<-""
heatmap_fdr[heatmap_fdr$ID=="syn",]$label<-"*"


levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")
colors<-c("darkorange2", "#f19133", "#f5ad66", "#f8c899","lightgrey", "#9999d1","#6666b9","#3333a2","darkblue")


p_heatmap_fdr<-ggplot(heatmap_fdr)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)),fill=factor(category, levels=levels)), size=0, alpha=0.8, color=NA)+
  theme_minimal()+
  theme()+
  scale_x_continuous(breaks=seq(1:42), labels = ABseq_pos, expand = c(0,0))+
  labs(x="AB(1-42) WT amino acid and position", y="Mutant amino acid", fill="Nucleation score")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x.top = element_line(),       
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=15, face = "bold"),
        legend.text = element_text(size=15), 
        axis.text.x = element_text(color=color_axis, size=12),
        
        axis.text.y = element_text(color=color_axis_y, size=12),
        axis.title = element_text(size = 20)
  )+
  geom_text(aes(Pos,Mut,label=label), size=8)+
  scale_fill_manual("Category (FDR)", values= colors, labels=levels)

p_heatmap_fdr

ggsave(p_heatmap_fdr,path=path, file="p_heatmap_FDR.pdf",width=14, height=8)









#####  stacked barplot nucleation score categories

fdr_categories<-singles[,c("Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10")]

fdr_categories$category<-"WT-like"

fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c<0),]$category<- "NS- 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS- 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c<0),]$category<- "NS- 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c<0),]$category<- "NS- 1%"

fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c>0),]$category<- "NS+ 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+ 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c>0),]$category<- "NS+ 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c>0),]$category<- "NS+ 1%"


categories <- fdr_categories %>% group_by(Pos,category) %>% dplyr::summarise(Freq=n()) 

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")
colors<-c("darkorange2", "#f19133", "#f5ad66", "#f8c899","lightgrey", "#9999d1","#6666b9","#3333a2","darkblue")

p_categories<-ggplot(categories, aes(fill=factor(category, levels=levels), x=factor(Pos), y=Freq)) + 
  theme_bw()+
  theme(legend.title=element_text(size=15, face="bold"),
        
        legend.text = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=25),
        axis.text = element_text(size=20)
  )+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation (%FDR)", values=colors)+
  labs( x= "AB(1-42) amino acid position",y="Counts")+
  scale_x_discrete(breaks=c(10,20,30,40)) +
  scale_y_continuous(expand = c(0,0))
p_categories

ggsave(p_categories, path=path, file="p_fdr_categories.pdf",width=12, height=3)




## defining gatekeepers


gatekeepers=c()

for (i in 1:42){
  
  
  l_wt<-nrow(subset(fdr_categories, Pos==i & category == "WT-like"))
  l_inc_1<-nrow(subset(fdr_categories, Pos==i & category == "NS+ 1%"))
  l_inc_5<-nrow(subset(fdr_categories, Pos==i & category == "NS+ 5%"))
  l_inc_10<-nrow(subset(fdr_categories, Pos==i & category == "NS+ 10%"))
  l_inc_25<-nrow(subset(fdr_categories, Pos==i & category == "NS+ 25%"))
  l_dec_1<-nrow(subset(fdr_categories, Pos==i & category == "NS- 1%"))
  l_dec_5<-nrow(subset(fdr_categories, Pos==i & category == "NS- 5%"))
  l_dec_10<-nrow(subset(fdr_categories, Pos==i & category == "NS- 10%"))
  l_dec_25<-nrow(subset(fdr_categories, Pos==i & category == "NS- 25%"))
  
  
  if((l_inc_1+l_inc_5+l_inc_10)>=(l_dec_1+l_dec_5+l_dec_10)
     &(l_inc_1+l_inc_5+l_inc_10+l_dec_1+l_dec_5+l_dec_10)>= l_wt
     
  ){
    
    gatekeepers<-c(gatekeepers, i)  
    
  }
  
}

gatekeepers








## position vs nucleation boxplots

p_position_NS<-ggplot(singles, aes(x=factor(category_10, levels=c("NS_inc", "WT-like", "NS_dec"), labels=c("+", "WT-like", "-")), y=Pos))+
  geom_violin(size=0.3)+
  theme_bw()+
  geom_jitter(width = 0.2, aes(color=nscore_c), size=2, alpha=0.8, shape=16)+
  geom_boxplot(width=0.15, fill=NA, size=0.4, outlier.shape = NA)+ 
  theme(axis.text.x=element_text(size=15),
        panel.border = element_blank(),
        axis.line.y = element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_gradientn(colours=c("darkorange2", "#f5ad66", "lightgrey", "darkblue"), breaks=c(-4, -2, 0, 2),na.value = "white")+
  labs(x="Nucleation", y="Amino acid position", color="Nucleation score")

p_position_NS
ggsave(p_position_NS, path=path, file="p_position_NS.pdf", width=4, height=2.5)





## NS vs clusters

singles$cluster<-""
singles[singles$Pos %in% c(1,3,7,11,17,22,42),]$cluster<-"gatekeeper"
singles[singles$Pos %in% c(2,4:6,8:10,12:16,18:21,23:26),]$cluster<-"2_26"
singles[singles$Pos %in% c(27:41),]$cluster<-"27_41"


p_clusters<-ggplot(singles, aes(x=factor(cluster, levels=c("gatekeeper", "2_26", "27_41"),
                                    labels=c("Gatekeepers\n(D1, E3, D7, E11,\nL17, E22 & A42)", "Nt (2-26)", "Ct (27-41)" )), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.2)+
  geom_violin(size=0.3)+
    theme_bw()+
  geom_jitter(width = 0.2,  size=2, alpha=0.8, shape=16, color="grey")+
  geom_boxplot(width=0.15, fill=NA, size=0.4, outlier.shape = NA)+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, vjust = 1, hjust=1),
        panel.border = element_blank(),
        axis.line.y = element_line(color="black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(y="Nucleation score")

p_clusters
ggsave(p_clusters,path=path, file="p_NS_clusters.pdf", width=3, height=2.5)




## controls yeast growth AB and SupN

melt_controls<-melt(yeast_controls, id.vars = "variant")

p_controls<-ggplot(melt_controls, aes(x=factor(variant, levels=c("supN", "AB", "supN-AB")), y=value))+

  theme_bw()+
  geom_jitter(width=0.2, size=3, aes(color=variant), show.legend = F)+
  scale_color_manual(values=c("darkblue","grey70","#15983DFF"))+
  theme(axis.title.x=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color="black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(y="Nucleation score")

p_controls
ggsave(p_controls,path=path, file="p_yeast_controls.pdf", width=3, height=2.5)



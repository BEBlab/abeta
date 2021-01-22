require(ggplot2)
require(reshape2)
require(dplyr)

dir.create("Determinants of A-beta nucleation")
path="Determinants of A-beta nucleation"

#required data:
load("nscore_df.RData")
load("truncations")


AA_type<-data.frame("AA"= c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P" ),
                   "name_AA"=c("Glycine", "Alanine","Valine","Leucine","Methionine","Isoleucine","Phenylalanine",
                   "Tyrosine","Tryptophan","Lysine","Arginine","Aspartic acid","Glutamic acid","Serine","Threonine",
                   "Cysteine","Asparagine","Glutamine","Histidine","Proline"),
                    "type"=c("glycine",rep("aliphatic",5),rep("aromatic",3),rep("positive",2),rep("negative",2),rep("polar",6),"proline"))
ABseq=c("D1","A2","E3","F4","R5","H6","D7","S8","G9","Y10","E11","V12","H13","H14","Q15","K16","L17","V18","F19","F20","A21","E22","D23","V24","G25","S26","N27","K28","G29","A30","I31","I32","G33","L34","M35","V36","G37","G38","V39","V40","I41","A42")

labels_ABseq=c("D\n1","A\n2","E\n3","F\n4","R\n5","H\n6","D\n7","S\n8","G\n9","Y\n10",
               "E\n11","V\n12","H\n13","H\n14","Q\n15","K\n16","L\n17","V\n18","F\n19","F\n20",
               "A\n21","E\n22","D\n23","V\n24","G\n25","S\n26","N\n27","K\n28","G\n29","A\n30",
               "I\n31","I\n32","G\n33","L\n34","M\n35","V\n36","G\n37","G\n38","V\n39","V\n40","I\n41","A\n42")
levels=c("aliphatic", "aromatic", "negative", "positive", "polar", "histidine", "glycine", "proline")

disease_mutations<- c("H-6-R", "D-7-N", "D-7-H", "E-11-K","K-16-N", "A-21-G", "E-22-Q", "E-22-K", "E-22-G", "D-23-N", "L-34-V", "A-42-T")

color=c("aliphatic"="darkgrey",
        "aromatic"="#9A703EFF",
        "negative"="#EE0011FF",
        "positive"="#0C5BB0FF", 
        "polar"="#15983DFF", 
        "glycine"="#FA6B09FF", 
        "proline"="#FEC10BFF")

color_axis=c("D1"="#EE0011FF","A2"="darkgrey","E3"="#EE0011FF","F4"="#9A703EFF","R5"="#0C5BB0FF","H6"="#15983DFF",
             "D7"="#EE0011FF","S8"="#15983DFF","G9"="#FA6B09FF","Y10"="#9A703EFF","E11"="#EE0011FF","V12"="darkgrey",
             "H13"="#15983DFF","H14"="#15983DFF","Q15"="#15983DFF","K16"="#0C5BB0FF","L17"="darkgrey","V18"="darkgrey",
             "F19"="#9A703EFF","F20"="#9A703EFF","A21"="darkgrey","E22"="#EE0011FF","D23"="#EE0011FF","V24"="darkgrey",
             "G25"="#FA6B09FF","S26"="#15983DFF","N27"="#15983DFF","K28"="#0C5BB0FF","G29"="#FA6B09FF","A30"="darkgrey",
             "I31"="darkgrey","I32"="darkgrey","G33"="#FA6B09FF","L34"="darkgrey","M35"="darkgrey","V36"="darkgrey",
             "G37"="#FA6B09FF","G38"="#FA6B09FF","V39"="darkgrey","V40"="darkgrey","I41"="darkgrey","A42"="darkgrey")



# add info of AA type for WT and Mut;  AA WT+position; and fAD info
singles$WT_type<-""
singles$Mut_type<-""

for(i in 1:nrow(singles)){
  singles[i,]$WT_type<-as.character(AA_type[AA_type$AA==singles[i,]$WT_AA,]$type)
  singles[i,]$Mut_type<-as.character(AA_type[AA_type$AA==singles[i,]$Mut,]$type)
}


singles$names<-paste0(singles$WT_AA, singles$Pos)


singles$fAD<-"non-fAD"
singles[singles$ID %in% disease_mutations,]$fAD<-"fAD"
singles[singles$ID=="A-2-V",]$fAD<-"recessive"









### master boxplot

p_boxplot<-ggplot(singles, aes(x = factor(names, levels=ABseq, labels=labels_ABseq), y = nscore_c)) +
  geom_hline(yintercept = 0, size=0.5, colour="black", show.legend = F)+
  geom_boxplot(outlier.shape = NA, color="grey", fill="white", size=0.3) +
  
  geom_point(data=singles, aes(y=0, color=factor(WT_type, levels=levels),size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1"))),alpha=0.8, shape=8, size=2.5)+ 
  geom_jitter(data=singles, aes(color=factor(Mut_type, levels=levels),size=factor(sig_10,levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)), width=0.2, alpha=0.8)+
  
  scale_color_manual(values=color)+
  scale_size_manual(values=c(2,4))+
  scale_shape_manual(values=c(17,16,18))+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color="black"),
        axis.text.x = element_text(color=color_axis, size=12),
        axis.title.y = element_text(size=15))+
  guides(fill=F)+
  labs(y="Nucleation score",x="AB(1-42) WT amino acid and position",shape="fAD category", size="FDR", color="AA type")
p_boxplot


ggsave(p_boxplot,path=path, file="master_boxplot.pdf", width = 15, height = 6)





### master boxplot with only Proline and Threonine

singles$outline<-"other"
singles[singles$Mut=="P",]$outline<-"proline"
singles[singles$Mut=="T",]$outline<-"threonine"


p_boxplot<-ggplot(singles, aes(x = factor(names, levels=ABseq, labels=labels_ABseq), y = nscore_c)) +
  geom_hline(yintercept = 0, size=0.5, colour="black", show.legend = F)+
  geom_boxplot( size=0.2, outlier.shape = NA, color="grey90", fill="white") +
  
  geom_jitter(data=singles[singles$outline=="other",],
              aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
              width=0.25, alpha=1, stroke=1, color="grey90", fill="grey90")+
  
  geom_jitter(data=singles[singles$outline %in% "proline",],
              aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
              color="grey30",fill="grey30", width=0.25, alpha=1, stroke=0.1)+
  geom_jitter(data=singles[singles$outline %in% "threonine",],
              aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
              color="grey30", fill="white", width=0.25, alpha=1, stroke=1)+
  

  scale_color_manual(values=c("#A1C720FF", "grey30" ))+
  scale_size_manual(values=c(1.5,3))+
  
  scale_shape_manual(values=c(24,21,23))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color="black"),
        axis.text.x = element_text(color="black", size=12),
        axis.title.y = element_text(size=15))+
  labs(y="Nucleation score",x="AB(1-42) WT amino acid and position",shape="fAD category", size="FDR")
p_boxplot


ggsave(p_boxplot,path=path, file="boxplot_PT.pdf", width = 12, height = 4)







### master boxplot with only Isoleucine and Valine

singles$outline<-"other"
singles[singles$Mut=="V",]$outline<-"valine"
singles[singles$Mut=="I",]$outline<-"isoleucine"


p_boxplot<-ggplot(singles, aes(x = factor(names, levels=ABseq, labels=labels_ABseq), y = nscore_c)) +
  geom_hline(yintercept = 0, size=0.5, colour="black", show.legend = F)+
  geom_boxplot(size=0.2, outlier.shape = NA, color="grey90", fill="white") +
 
  geom_jitter(data=singles[singles$outline=="other",],
              aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
              width=0.25, alpha=1, stroke=1, color="grey90", fill="grey90")+
  
  geom_jitter(data=singles[singles$outline %in% "valine",],
              aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
              color="grey30",fill="grey30", width=0.25, alpha=1, stroke=0.1)+
  geom_jitter(data=singles[singles$outline %in% "isoleucine",],
              aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
              color="grey30", fill="white",  width=0.25, alpha=1, stroke=1)+
  
  geom_point(data=singles[singles$WT_AA=="V",], 
             color="grey30",     aes(y=0),  alpha=0.8, shape=8, size=2.5)+
  geom_point(data=singles[singles$WT_AA=="I",], 
             color="grey30",     aes(y=0),  alpha=0.8, shape=4, size=2.5)+
  
  scale_color_manual(values=c("#EE0011FF", "#0C5BB0FF" ))+
  scale_size_manual(values=c(1.5,3))+
  scale_shape_manual(values=c(24,21,23))+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color="black"),
        axis.text.x = element_text(color="black", size=12),
        axis.title.y = element_text(size=15))+
 
  labs(y="Nucleation score",x="AB(1-42) WT amino acid and position",shape="fAD category", size="FDR")
p_boxplot

ggsave(p_boxplot,path=path, file="boxplot_IV.pdf", width = 12, height = 4)








### master boxplot for each AA


for (AA in unique(singles$Mut)){
  if(length(which(singles$WT_AA==AA))==0){
    

p_AA<-ggplot(singles, aes(x = factor(names, levels=ABseq, labels=labels_ABseq), y = nscore_c)) +
      
      geom_hline(aes(yintercept = 0), size=0.15, colour="black")+
      geom_boxplot( size=0.2, outlier.shape = NA, color="grey90", fill="white") +
      
      geom_jitter(data=singles[singles$Mut!=AA,],
                  aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
                  width=0.2, alpha=1, stroke=1, color="grey90", fill="grey90")+
      
      
      geom_jitter(data=singles[singles$Mut==AA,], 
                  aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
                  color="grey30",fill="grey30", width=0.2, alpha=1, stroke=0.1)+
      
      scale_size_manual(values=c(2.5,4))+
      scale_shape_manual(values=c(17,16, 18))+
      
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.title=element_text(size=25),
            axis.line.y = element_line(color="black", size=0.15),
            axis.text.x = element_text(color="black", size=14),
            axis.text.y = element_text(size=14))+
      labs(title=AA_type[AA_type$AA==AA,]$name_AA, shape="fAD class", size="FDR",y="Nucleation score",x="AB(1-42) WT amino acid and position")
    
    ggsave(p_AA, path=path, file=paste0(AA, "_AA_boxplot.pdf"), width = 12, height = 4)
    
    
  }else{
    
    
    p_AA<-ggplot(singles, aes(x = factor(names, levels=ABseq, labels=labels_ABseq), y = nscore_c)) +
      
      geom_hline(aes(yintercept = 0), size=0.15, colour="black", show.legend = F)+
      geom_boxplot( size=0.2, outlier.shape = NA, color="grey90", fill="white") +
      
      geom_jitter(data=singles[singles$Mut!=AA,],
                  aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),  
                  width=0.2, alpha=1, stroke=1, color="grey90", fill="grey90")+
      
      geom_point(data=singles[singles$WT_AA==AA,],  
                 color="grey30",     aes(y=0),  alpha=1, shape=8, size=3)+
      
      geom_jitter(data=singles[singles$Mut==AA,],    
                  aes(size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")), shape=factor(fAD)),
                  color="grey30",fill="grey30",  width=0.2, alpha=1, stroke=0.1)+
      
      
      scale_size_manual(values=c(2.5,4))+
      scale_shape_manual(values=c(17,16, 18))+
      
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.title=element_text(size=25),
            axis.line.y = element_line(color="black", size=0.15),
            axis.text.x = element_text(color="black", size=14),
            axis.text.y = element_text(size=14)
            )+
      labs(title=AA_type[AA_type$AA==AA,]$name_AA, shape="fAD class", size="FDR", y="Nucleation score",x="AB(1-42) WT amino acid and position")
    
    
    ggsave(p_AA,path=path, file=paste0(AA, "_AA_boxplot.pdf"), width = 12, height = 4)
  }}








### boxplots to each Mut AA

list<-list(
  "Nt (2-26)"=c(2,4:6,8:10,12:16,18:21,24:26),
  "Ct (27-41)"=c(27:41),
  "Negatively charged gatekeepers"=c(1,3,7,11,22),
  "Negatively charged positions"=c(1,3,7,11,22,23)
)



for(region in names(list)){
  
  residues<-list[[region]]
  
  
p<-ggplot(singles[singles$Pos %in% residues,], aes(x= reorder(Mut, -nscore_c), y = nscore_c, fill=Mut_type)) +
    theme_bw()+
    geom_hline(aes(yintercept = 0), size=0.2, colour="black", show.legend = F)+
    geom_boxplot( width=0.5, size=0.2,outlier.shape = NA, alpha=0.8) +
    theme(axis.text.x = element_text( size=15, hjust=0.5), 
          axis.title.y = element_text(size=10),
          title = element_text(size=15),
          axis.line.y = element_line(size=0.2),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(y="Nucleation score", title=paste0(region), fill="Mutant AA type", x="Mutant amino acid")+
    scale_fill_manual(values=color)+
    geom_jitter( position = position_jitter(width = 0.2), size=1, show.legend = T, alpha=0.8, shape=16)
  
p

ggsave(p,path = path, file=paste0(region, "_to_MutAA_boxplot.pdf"), width = 8, height = 3)


  
}





## truncations

melt_trunc<-melt(truncations, id.vars = "variant")

p_trunc<-ggplot(melt_trunc, aes(x=factor(variant, levels=c("supN", "supN-AB","supN-AB(22-42)","supN-AB(24-42)","supN-AB(27-42)" )), y=value))+
  
  theme_bw()+
  geom_jitter(width=0.2, size=3, aes(color=variant), show.legend = F)+
  scale_color_manual(values=c("grey","#15983DFF", rep("darkblue", 3)))+
  theme(axis.title.x=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color="black", size=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=45, vjust = 1, hjust=1))+
  labs(y="% -Ade growth")

p_trunc
ggsave(p_trunc,file="p_yeast_trunc.pdf", width=3, height=3)



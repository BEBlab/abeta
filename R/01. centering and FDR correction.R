require(readxl)
require(ggplot2)
require(ggpubr)
require(hexbin)


dir.create("Centering and FDR correction")
path="Centering and FDR correction"


#required data:
#excel file MS_BL_BB_processed_data.xlsx
#deposited in NCBI's Gene Expression Omnibus (GEO) as record GSE151147 



#import data
silent<-read_excel("MS_BL_BB_processed_data.xlsx", sheet="synonymous")
singles<-read_excel("MS_BL_BB_processed_data.xlsx", sheet="1 aa change")
doubles<-read_excel("MS_BL_BB_processed_data.xlsx", sheet="2 aa changes")


#centering to the weighted mean of synonymous with 1 mut codon
mean_syn_1codon<-weighted.mean(silent[silent$Nmut_codons==1,]$nscore, silent[silent$Nmut_codons==1,]$sigma^-2, na.rm = T)

###silent
silent$nscore_c<-as.numeric(paste(silent$nscore+(-mean_syn_1codon)))
silent$ID<-"silent"
silent$Mut<-"silent"


###singles
singles$nscore_c<-as.numeric(paste(as.numeric(singles$nscore)+(-mean_syn_1codon)))
singles$nscore1_c<-as.numeric(paste(as.numeric(singles$nscore1)+(-mean_syn_1codon)))
singles$nscore2_c<-as.numeric(paste(as.numeric(singles$nscore2)+(-mean_syn_1codon)))
singles$nscore3_c<-as.numeric(paste(as.numeric(singles$nscore3)+(-mean_syn_1codon)))

singles$ID<-paste(singles$WT_AA, singles$Pos, singles$Mut, sep = "-")

# FDR=0.1 correction and assignment into categories

singles$zscore<-singles$nscore_c/singles$sigma
singles$p.adjust<-p.adjust(2*pnorm(-abs(singles$zscore)), method = "BH")

singles$sig_10<-FALSE
singles[singles$p.adjust<0.1,]$sig_10<-TRUE

singles$category_10<-"WT-like"
singles[singles$sig_10==T & singles$nscore_c<0,]$category_10<-"NS_dec"
singles[singles$sig_10==T & singles$nscore_c>0,]$category_10<-"NS_inc"

#with stops
singles_stops<-singles
#remove stops
singles<-singles[singles$Mut!="*",]



### doubles
doubles$nscore_c<-as.numeric(paste(as.numeric(doubles$nscore)+(-mean_syn_1codon)))
doubles$nscore1_c<-as.numeric(paste(as.numeric(doubles$nscore1)+(-mean_syn_1codon)))
doubles$nscore2_c<-as.numeric(paste(as.numeric(doubles$nscore2)+(-mean_syn_1codon)))
doubles$nscore3_c<-as.numeric(paste(as.numeric(doubles$nscore3)+(-mean_syn_1codon)))

doubles$ID_mut1<-paste(doubles$WT_AA1, doubles$Pos1, doubles$Mut1, sep = "-")
doubles$ID_mut2<-paste(doubles$WT_AA2, doubles$Pos2, doubles$Mut2, sep = "-")
doubles$ID<-paste(doubles$ID_mut1, doubles$ID_mut2, sep="_")


# FDR=0.1 correction and assignment into categories

doubles$zscore<-doubles$nscore_c/doubles$sigma
doubles$p.adjust<-p.adjust(2*pnorm(-abs(doubles$zscore)), method = "BH")


doubles$sig_10<-FALSE
doubles[doubles$p.adjust<0.1,]$sig_10<-TRUE

doubles$category_10<-"WT-like"
doubles[doubles$sig_10==T & doubles$nscore_c<0,]$category_10<-"NS_dec"
doubles[doubles$sig_10==T & doubles$nscore_c>0,]$category_10<-"NS_inc"


#remove stops
doubles<-doubles[!doubles$Mut1=="*",]
doubles<-doubles[!doubles$Mut2=="*",]


### singles and doubles
singles_doubles<-singles
colnames(singles_doubles)[c(1,2,3)]<-c("Pos1", "WT_AA1", "Mut1")
singles_doubles[,c("Pos2","WT_AA2", "Mut2", "ID_mut1", "ID_mut2")]<-NA
singles_doubles<-rbind(singles_doubles, doubles)


save(silent, singles, singles_stops, singles_doubles, doubles, file="nscore_df.RData")




#### plots ####

# distribution singles (Missense, Nonsense and Synonymous)

dist<-singles_stops[,c("Mut", "nscore_c")]
dist<-rbind(dist, silent[silent$Nmut_codons==1,c("Mut", "nscore_c")])
dist$type<-"Missense"
dist[dist$Mut=="*",]$type<-"Nonsense"
dist[dist$Mut=="silent",]$type<-"Synonymous"

p_hist<-ggplot(dist, aes(x=nscore_c))+
  geom_histogram(binwidth = 0.2, position="identity", alpha=0.8, color=NA, aes(fill=factor(type, levels=c("Missense", "Nonsense", "Synonymous"))))+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  theme_bw()+
 #scale_fill_manual(values=c("grey70", "grey40", "grey2"))+
  scale_fill_manual(values=c("grey70", "darkblue","#15983DFF"))+
  labs(x="Nucleation score", y="Counts")+
  theme(legend.position = c(0.3,0.8),
        legend.text = element_text(size=18),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 25),
        axis.text = element_text(size=18),
        strip.text = element_text(size=15),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))+
  xlim(-6,3.5)
p_hist


# distribution doubles (don't contain stops)

dist<-doubles
dist$type<-"Missense"

p_hist_doubles<-ggplot(dist, aes(x=nscore_c) )+
  geom_histogram(binwidth = 0.2, position="identity", alpha=0.8, color=NA, aes(fill=factor(type, levels=c("Missense"))))+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  theme_bw()+
  scale_fill_manual(values="grey70")+
  labs(x="Nucleation score", y="Counts")+
  theme(legend.position = c(0.9,0.8),
        legend.text = element_text(size=18),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 25),
        axis.text = element_text(size=18),
        strip.text = element_text(size=15),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))+
  xlim(c(-6,3.5))

p_hist_doubles

#arragne singles and doubles
arr<-ggarrange(p_hist, p_hist_doubles, ncol=1, align = "v")
ggsave(arr, path=path, file="p_hist_singles_doubles.pdf",width=6, height=8, useDingbats=FALSE)



# correlation between replicates

#1 and 2
subset<-singles_doubles[!is.na(singles_doubles$nscore1_c),]
subset<-subset[!is.na(subset$nscore2_c),]
n<-length(subset$ID)

corr<-cor.test(singles_doubles$nscore1_c, singles_doubles$nscore2_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr<-ggplot(singles_doubles, aes(x=nscore1_c, y=nscore2_c) )+
  stat_binhex()+
  theme_bw()+
  labs(x="Replicate 1", y="Replicate 2")+
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'),
          axis.title  = element_text(size = 25),
          axis.text = element_text(size=18))+
  annotate("text", x = -5, y = 3, label = paste0("R=", round(R, 2)), size=8)+
  annotate("text", x = -5, y = 2, label = paste0("p=",format(p, digits = 2, scientific = T)), size=8)+

  scale_fill_gradient(high="grey30", low="grey90")+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")
p_corr


ggsave(p_corr,path=path, file="p_corr reps 1 and 2.pdf",width=5, height=4, useDingbats=FALSE)

#1 and 3
subset<-singles_doubles[!is.na(singles_doubles$nscore1_c),]
subset<-subset[!is.na(subset$nscore3_c),]
n<-length(subset$ID)

corr<-cor.test(singles_doubles$nscore1_c, singles_doubles$nscore3_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr<-ggplot(singles_doubles, aes(x=nscore1_c, y=nscore3_c) )+
  stat_binhex()+
  theme_bw()+
  labs(x="Replicate 1", y="Replicate 3")+
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'),
          axis.title  = element_text(size = 25),
          axis.text = element_text(size=18))+
  annotate("text", x = -5, y = 3, label = paste0("R=", round(R, 2)), size=8)+
  annotate("text", x = -5, y = 2, label = paste0("p=",format(p, digits = 2, scientific = T)), size=8)+
  
  scale_fill_gradient(high="grey30", low="grey90")+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")
p_corr


ggsave(p_corr, path=path, file="p_corr reps 1 and 3.pdf",width=5, height=4, useDingbats=FALSE)


#2 and 3
subset<-singles_doubles[!is.na(singles_doubles$nscore2_c),]
subset<-subset[!is.na(subset$nscore3_c),]
n<-length(subset$ID)

corr<-cor.test(singles_doubles$nscore2_c, singles_doubles$nscore3_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr<-ggplot(singles_doubles, aes(x=nscore2_c, y=nscore3_c) )+
  stat_binhex()+
  theme_bw()+
  labs(x="Replicate 2", y="Replicate 3")+
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'),
          axis.title  = element_text(size = 25),
          axis.text = element_text(size=18))+
  annotate("text", x = -5, y = 3, label = paste0("R=", round(R, 2)), size=8)+
  annotate("text", x = -5, y = 2, label = paste0("p=",format(p, digits = 2, scientific = T)), size=8)+
  
  scale_fill_gradient(high="grey30", low="grey90")+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")
p_corr


ggsave(p_corr,path=path, file="p_corr reps 2 and 3.pdf",width=5, height=4, useDingbats=FALSE)


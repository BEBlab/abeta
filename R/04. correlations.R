require(dplyr)
require(reshape2)
require(weights)
require(ggplot2)
require(hexbin)

dir.create("Correlations")
path="Correlations"

##required data:
load("nscore_df.RData")
load("predictors_scores.RData")
load("loadings_properties.RData")
load("solubility_dataset.RData")
load("kinetics.RData")




### aggregation predictors correlation plots

list<-list(
  "Nt (2_26)"=c(2,4:6,8:10,12:16,18:21,24:26),
  "Ct (27_41)"=c(27:41)
)

predictors<-c("Tango","Waltz","Zyggregator","Camsol" )
singles_doubles_predictors<-left_join(singles_doubles[,c("Pos1", "Pos2", "ID",  "Nmut_codons","nscore_c", "sigma")], predictors_scores[,c("ID", predictors)], by="ID")


for(region in names(list)){
  
  residues<-list[[region]]
  
  subset_region_singles<-singles_doubles_predictors[singles_doubles_predictors$Nmut_codons==1 & singles_doubles_predictors$Pos1 %in% residues,]
  subset_region_doubles<-singles_doubles_predictors[singles_doubles_predictors$Nmut_codons==2 & singles_doubles_predictors$Pos1 %in% residues & singles_doubles_predictors$Pos2 %in% residues ,]
  subset_region<-rbind(subset_region_singles, subset_region_doubles)
  
  melt_subset_region<-melt(subset_region, id.vars = c("Pos1", "Pos2", "ID",  "Nmut_codons","nscore_c", "sigma"))
  
  corr_vector=c()
  pval_vector=c()
  for(i in predictors){
    
    subset<-subset.data.frame(subset_region[,c(i,"nscore_c","sigma")])
    subset<-na.omit(subset)
    
    #pearson
    pearson<-wtd.cor(subset$nscore_c, subset[,i], weight=subset$sigma^-2)
    
    corr<-pearson[1,1]
    corr_vector=c(corr_vector,corr)
    
    p.value<-pearson[1,4]
    pval_vector=c(pval_vector, as.numeric(p.value))
    
    
  }
  
  corr_text <- data.frame(
    label = corr_vector,
    variable=predictors)
  
  pval_text <- data.frame(
    label = as.numeric(pval_vector),
    variable=predictors)
  
p<-ggplot(melt_subset_region,aes(x=value,y=nscore_c))+
    facet_wrap(~variable, nrow = 2, scales = "free_x")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black', size=0.15),
          strip.background = element_rect(fill=NA, color=NA))+
    
    stat_binhex()+
    scale_fill_gradient(high="black", low="grey85")+
    labs(y="Nucleation score",x="Score")+
    ylim(-6.5,3.5)+
    geom_text(data=corr_text, aes(label=paste0("R=",round(label, 2)), x=-Inf, y=Inf),hjust=-0.15,vjust=1.5, size=4, colour="black")+
    geom_text(data=pval_text, aes(label=paste0("p=",format(label, digits = 2, scientific = T)), x=-Inf, y=Inf),hjust=-0.1, vjust=3, size=4, colour="black")+
    geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")

  
ggsave(p, path=path, file=paste0(region, "_predictors_wrap.pdf"), width = 6, height = 5, useDingbats=FALSE)

}






### hydrophobicity (Kyte-Doolittle, 1982) and PC1 (Bolognesi, 2019)

#first calculate the hydrophobicity and PC1 scores for all variants 
#open 'loadings_properties.RData' 

singles_doubles_properties<-singles_doubles[,c("Pos1", "WT_AA1", "Mut1", "Nmut_codons","sigma", "nscore_c", "Pos2", "WT_AA2", "Mut2", "ID")]

singles_doubles_properties$hydrophobicity<-""
singles_doubles_properties$PC1<-""

singles_properties<-singles_doubles_properties[singles_doubles_properties$Nmut_codons==1,]

for(i in 1:nrow(singles_properties)){
  WT_score<-loadings_properties[loadings_properties$amino_acid==singles_properties[i,]$WT_AA1,]$`Hydrophobicity (Kyte-Doolittle)`  
  Mut_score<-loadings_properties[loadings_properties$amino_acid==singles_properties[i,]$Mut1,]$`Hydrophobicity (Kyte-Doolittle)`
  singles_properties[i,]$hydrophobicity<-paste(Mut_score-WT_score)

  WT_score<-loadings_properties[loadings_properties$amino_acid==singles_properties[i,]$WT_AA1,]$PC1  
  Mut_score<-loadings_properties[loadings_properties$amino_acid==singles_properties[i,]$Mut1,]$PC1
  singles_properties[i,]$PC1<-paste(Mut_score-WT_score)
  
}


doubles_properties<-singles_doubles_properties[singles_doubles_properties$Nmut_codons==2,]

for(i in 1:nrow(doubles_properties)){
  WT_score<-paste(as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_properties[i,]$WT_AA1,]$`Hydrophobicity (Kyte-Doolittle)`)+
            as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_properties[i,]$WT_AA2,]$`Hydrophobicity (Kyte-Doolittle)`))
  Mut_score<-paste(as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_properties[i,]$Mut1,]$`Hydrophobicity (Kyte-Doolittle)`)+
                    as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_properties[i,]$Mut2,]$`Hydrophobicity (Kyte-Doolittle)`))
  doubles_properties[i,]$hydrophobicity<-paste(as.numeric(Mut_score)-as.numeric(WT_score))
  
  
  
  WT_score<-paste(as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_properties[i,]$WT_AA1,]$PC1)+
                    as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_properties[i,]$WT_AA2,]$PC1))
  Mut_score<-paste(as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_properties[i,]$Mut1,]$PC1)+
                     as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_properties[i,]$Mut2,]$PC1))
  doubles_properties[i,]$PC1<-paste(as.numeric(Mut_score)-as.numeric(WT_score))
  
}

singles_doubles_properties<-rbind(singles_properties, doubles_properties)




### plot NS vs hydrophobicity (Nt and Ct)

list<-list(
  "Nt (2_26)"=c(2,4:6,8:10,12:16,18:21,24:26),
  "Ct (27_41)"=c(27:41)
)


for(region in names(list)){
  
  residues<-list[[region]]
  
  subset_region_singles<-singles_doubles_properties[singles_doubles_properties$Nmut_codons==1 & singles_doubles_properties$Pos1 %in% residues,]
  subset_region_doubles<-singles_doubles_properties[singles_doubles_properties$Nmut_codons==2 & singles_doubles_properties$Pos1 %in% residues & singles_doubles_properties$Pos2 %in% residues ,]
  subset_region<-rbind(subset_region_singles, subset_region_doubles)
  
  #pearson
  pearson<-wtd.cor(subset_region$nscore_c, as.numeric(subset_region$hydrophobicity), weight=subset_region$sigma^-2)
    
  corr<-pearson[1,1]
  p.value<-pearson[1,4]

    
p_corr<-ggplot(subset_region, aes(x=as.numeric(hydrophobicity), y=nscore_c) )+
  stat_binhex()+
  theme_bw()+
  labs(x="Hydrophobicity (Kyte-Doolittle, 1982)", y="Nucleation score")+
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'))+
  scale_fill_gradient(high="grey20", low="grey90")+
  annotate("text", x = -12, y = 3, label = paste0("R=",format(corr, digits = 2, scientific = T)), size=4)+
  annotate("text", x = -12, y = 2, label =paste0("p=",format(p.value, digits = 2, scientific = T)), size=4)+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")

p_corr

ggsave(p_corr, path=path, file=paste0(region,"_NS_hydrophobicity.pdf"),width=5, height=4, useDingbats=FALSE)

}





### plot NS vs PC1 (Bolognesi 2019)

#PC1 is multiplied by -1 to correlate with hydrophobicity scales- therefore the considered hydrophobic AA have positive scores

#pearson
pearson<-wtd.cor(singles_doubles_properties$nscore_c, -as.numeric(singles_doubles_properties$PC1), weight=singles_doubles_properties$sigma^-2)
corr<-pearson[1,1]
p.value<-pearson[1,4]

p_corr<-ggplot(singles_doubles_properties, aes(x=-as.numeric(PC1), y=nscore_c) )+
  stat_binhex()+
  theme_bw()+
  labs(x="PC1 (Hydrophobicity)", y="Nucleation score")+
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'),
          axis.title  = element_text(size = 25),
          axis.text = element_text(size=18))+
  scale_fill_gradient(high="grey30", low="grey90")+
  annotate("text", x = -50, y = 3, label = paste0("R=",format(corr, digits = 2, scientific = T)), size=4)+
  annotate("text", x = -50, y = 2, label =paste0("p=",format(p.value, digits = 2, scientific = T)), size=4)+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")

p_corr

ggsave(p_corr, path=path, file="p_NS_PC1.pdf",width=5, height=4, useDingbats=FALSE)







### plot Solubility score (Gray 2019) vs PC1

#first calculate the PC1 scores for all single and double AA mutants available in Gray et al, 2019
#open 'loadings_properties.RData' 

singles_doubles_solubility<-solubility_dataset
singles_doubles_solubility$PC1<-""

singles_solubility<-singles_doubles_solubility[singles_doubles_solubility$Nmut_codons==1,]

for(i in 1:nrow(singles_solubility)){
  WT_score<-loadings_properties[loadings_properties$amino_acid==singles_solubility[i,]$WT_AA1,]$PC1  
  Mut_score<-loadings_properties[loadings_properties$amino_acid==singles_solubility[i,]$Mut1,]$PC1
  singles_solubility[i,]$PC1<-paste(Mut_score-WT_score)
  
}


doubles_solubility<-singles_doubles_solubility[singles_doubles_solubility$Nmut_codons==2,]

for(i in 1:nrow(doubles_solubility)){

  WT_score<-paste(as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_solubility[i,]$WT_AA1,]$PC1)+
                    as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_solubility[i,]$WT_AA2,]$PC1))
  Mut_score<-paste(as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_solubility[i,]$Mut1,]$PC1)+
                     as.numeric(loadings_properties[loadings_properties$amino_acid==doubles_solubility[i,]$Mut2,]$PC1))
  doubles_solubility[i,]$PC1<-paste(as.numeric(Mut_score)-as.numeric(WT_score))
  
}

singles_doubles_solubility<-rbind(singles_solubility, doubles_solubility)


### plot

#pearson
pearson<-wtd.cor(as.numeric(singles_doubles_solubility$score), -as.numeric(singles_doubles_solubility$PC1))
corr<-pearson[1,1]
p.value<-pearson[1,4]

p_corr<-ggplot(singles_doubles_solubility, aes(x=-as.numeric(PC1), y=as.numeric(score)))+
  stat_binhex()+
  theme_bw()+
  labs(x="PC1 (Hydrophobicity)", y="Solubility score\n(Gray 2019)")+
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'),
          axis.title  = element_text(size = 25),
          axis.text = element_text(size=18))+
  scale_fill_gradient(high="grey30", low="grey90")+
  annotate("text", x = 20, y = 1.5, label = paste0("R=",round(corr,2)), size=4)+
  annotate("text", x = 20, y = 1, label =paste0("p=",format(p.value, digits = 2, scientific = T)), size=4)+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")

p_corr

ggsave(p_corr,path=path, file="p_solubility_PC1.pdf",width=5.5, height=4, useDingbats=FALSE)







### in vitro measurements of A-beta aggregation

kinetics<-left_join(kinetics, singles[,c("ID", "nscore_c", "sigma")], by="ID")
kinetics<-left_join(kinetics, solubility_dataset[,c("ID", "score")], by="ID")
kinetics$score<-as.numeric(kinetics$score)
kinetics[kinetics$ID=="WT", c("nscore_c", "sigma", "score")]<-2.2e-16

# calculate Kn and K2 parameters
kinetics$kn<-paste(as.numeric(kinetics$primary)/as.numeric(kinetics$elongation))
kinetics$k2<-paste(as.numeric(kinetics$secondary)/as.numeric(kinetics$elongation))



# correlation with nucleation score

corr_vector=c()
pval_vector=c()
measurements<-c("primary","secondary","elongation","nuc_conversion","saturation","kn","k2")

for (i in measurements){
  
  #pearson
  pearson<-wtd.cor(kinetics$nscore_c, log10(as.numeric(kinetics[[i]])), weight = kinetics$sigma^-2)
  
  corr<-pearson[1,1]
  corr_vector=c(corr_vector,corr)
  
  p.value<-pearson[1,4]
  pval_vector=c(pval_vector, as.numeric(p.value))
  
}


corr_text <- data.frame(
  label = corr_vector,
  variable=measurements)

pval_text <- data.frame(
  label = as.numeric(pval_vector),
  variable=measurements)


melt_kinetics<-melt(kinetics, id=c("nscore_c", "sigma", "ID", "score"))

p_kinetics<-ggplot(melt_kinetics, aes(x=as.numeric(value), y=as.numeric(nscore_c)))+
  facet_wrap(vars(factor(variable, levels=c("primary", "secondary", "elongation", "nuc_conversion", "saturation", "kn", "k2"),
                         labels=c("Primary\nnucleation\nk+kn", "Secondary\nnucleation\nk+k2", "Elongation\nk+", "Conversion of\nsecondary nuclei",
                                  "Saturation of\nsecondary nucleation", "primary/elongation\nkn","secondary/elongation\nk2"))),
             scales = "free_x", ncol=5)+
  
  geom_smooth(method='lm',linetype = 2, size=1,  color = "#15983DFF", fill="grey75")+
  scale_x_continuous(trans="log10")+
  geom_point(size=3, color="#15983DFF")+
  labs(y="Nucleation score")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.line = element_line(color="black", size=0.15),
        strip.background=element_rect( color=NA, fill=NA))+
  geom_text(data=corr_text, aes(label=paste0("R=",round(label, 2)), x=0, y=3.5),hjust=-0.15,vjust=1.5, size=4, colour="black")+
  geom_text(data=pval_text, aes(label=paste0("p=",format(label, digits = 2, scientific = T)), x=0, y=3.5),hjust=-0.1, vjust=3, size=4, colour="black")

p_kinetics  


ggsave(p_kinetics,path=path, file="p_kinetics_NS.pdf",width=12, height=6, useDingbats=FALSE)






# correlation with solubility score

corr_vector=c()
pval_vector=c()
measurements<-c("primary","secondary","elongation","nuc_conversion","saturation","kn","k2")

for (i in measurements){
  
  #pearson
  pearson<-wtd.cor(as.numeric(kinetics$score), log10(as.numeric(kinetics[[i]])))
  
  corr<-pearson[1,1]
  corr_vector=c(corr_vector,corr)
  
  p.value<-pearson[1,4]
  pval_vector=c(pval_vector, as.numeric(p.value))
  
}


corr_text <- data.frame(
  label = corr_vector,
  variable=measurements)

pval_text <- data.frame(
  label = as.numeric(pval_vector),
  variable=measurements)


melt_kinetics<-melt(kinetics, id=c("nscore_c", "sigma", "ID", "score"))

p_kinetics_solubility<-ggplot(melt_kinetics, aes(x=as.numeric(value), y=as.numeric(score)))+
  facet_wrap(vars(factor(variable, levels=c("primary", "secondary", "elongation", "nuc_conversion", "saturation", "kn", "k2"),
                         labels=c("Primary\nnucleation\nk+kn", "Secondary\nnucleation\nk+k2", "Elongation\nk+", "Conversion of\nsecondary nuclei",
                                  "Saturation of\nsecondary nucleation", "primary/elongation\nkn","secondary/elongation\nk2"))),
             scales = "free_x", ncol=5)+
  
  geom_smooth(method='lm',linetype = 2, size=1,  color = "black", fill="grey75")+
  scale_x_continuous(trans="log10")+
  geom_point(size=3, color="black")+
  labs(y="Solubility score (Gray 2019)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.line = element_line(color="black", size=0.15),
        strip.background=element_rect( color=NA, fill=NA))+
  geom_text(data=corr_text, aes(label=paste0("R=",round(label, 2)), x=0, y=1),hjust=-0.15,vjust=1.5, size=4, colour="black")+
  geom_text(data=pval_text, aes(label=paste0("p=",format(label, digits = 2, scientific = T)), x=0, y=1),hjust=-0.1, vjust=3, size=4, colour="black")

p_kinetics_solubility  

ggsave(p_kinetics_solubility,path=path, file="p_kinetics_solubility.pdf",width=12, height=6, useDingbats=FALSE)




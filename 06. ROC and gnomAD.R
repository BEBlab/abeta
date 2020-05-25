require(dplyr)
require(pROC)
require(ggplot2)
require(ggrepel)


dir.create("ROC")
path="ROC"

##required data:
load("nscore_df.RData")
load("predictors_scores.RData")
load("solubility_dataset.RData")
load("gnomad.RData")

disease_mutations<- c("H-6-R", "D-7-N","D-7-H","E-11-K", "K-16-N", "A-21-G","E-22-Q", "E-22-K", "E-22-G", "D-23-N", "L-34-V",  "A-42-T")

#test fAD variants for significance
NS_merged<-sum((singles[singles$ID %in% disease_mutations,]$nscore_c)/(singles[singles$ID %in% disease_mutations,]$sigma^2))/sum(1/(singles[singles$ID %in% disease_mutations,]$sigma^2))
sigma_merged<-sqrt(1/sum(1/singles[singles$ID %in% disease_mutations,]$sigma^2))
z_score <- (NS_merged-0)/sigma_merged
p_value <- pnorm(z_score, lower.tail = FALSE)*2



#only single AA mutants

ROC_df<-left_join(singles[,c("ID", "nscore_c", "sigma")], predictors_scores, by="ID")
colnames(ROC_df)[2]<-"Nucleation score"

ROC_df$fAD_12<-0
ROC_df[ROC_df$ID %in% disease_mutations,]$fAD_12<-1



predictors<-c("Nucleation score","Tango","Waltz","Zyggregator","Camsol","Polyphen","CADD")




#build ROC with pROC package

mylist_fAD12 = list()
my_vector_auc = c()

for(i in predictors){
  subset<-ROC_df[,c("fAD_12",i)]
  subset<-na.omit(subset)
  
  glm.fit=glm(subset$fAD_12~subset[[i]], family = binomial)
  roc<-roc(subset$fAD_12, glm.fit$fitted.values)
  mylist_fAD12[[i]] <- roc
  
  auc<-roc$auc
  my_vector_auc=c(my_vector_auc, auc, i)
}


#extract AUC
auc<-as.data.frame(matrix(my_vector_auc, ncol=2, byrow=T))
colnames(auc)<-c("AUC", "Predictor")
auc$AUC<-round(as.numeric(as.character(auc$AUC)), 2)


#prepare labels x plots
labels_auc=c()
for(i in 1:nrow(auc)){
  auc_lab<-paste0(auc[i,]$Predictor," (AUC=", auc[i,]$AUC, ")" )
  pred<-as.character(auc[i,]$Predictor)
  labels_auc=c(labels_auc, auc_lab, pred)
}


label<-as.data.frame(matrix(labels_auc, ncol=2, byrow=T))
colnames(label)<-c("label", "Predictor")




#### plots


##NS, Tango, polyphen
colors=c("#EE0011FF", "#15983DFF", "darkgrey")

labels_plot<-c(as.character(label[label$Predictor=="Nucleation score",]$label),
               as.character(label[label$Predictor=="Tango",]$label),
               as.character(label[label$Predictor=="Polyphen",]$label))



g_fAD12<-ggroc(list(mylist_fAD12$`Nucleation score`, 
                    mylist_fAD12$`Tango`,
                    mylist_fAD12$Polyphen),
               legacy.axes = T, size=1)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color='black'),
    axis.title = element_text(size=20),
    axis.text = element_text(size=15)
  )+
  labs(x="False Positive Rate", y="True Positive Rate")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype="dashed", size=0.15)+
  
  scale_color_manual(labels=labels_plot, values=colors)
g_fAD12


ggsave(g_fAD12, path=path, file="ROC_NS_tango_polyphen.pdf", width = 7, height = 4)





##CADD
colors=c("darkgrey")
labels_plot<-as.character(label[label$Predictor=="CADD",]$label)

g_fAD12<-ggroc(list(mylist_fAD12$CADD),
               legacy.axes = T, size=1)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color='black'),
    axis.title = element_text(size=20),
    axis.text = element_text(size=15)
  )+
  labs(x="False Positive Rate", y="True Positive Rate")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype="dashed")+
  scale_color_manual(labels=labels_plot, values=colors)
g_fAD12

ggsave(g_fAD12, path=path, file="ROC_CADD.pdf", width = 7, height = 4)






##camsol, ZGG, waltz
colors=c("#15983DFF", "#A1C720FF", "#FEC10BFF")
labels_plot<-c(as.character(label[label$Predictor=="Zyggregator",]$label),
               as.character( label[label$Predictor=="Camsol",]$label),
               as.character( label[label$Predictor=="Waltz",]$label))


g_fAD12<-ggroc(list(mylist_fAD12$Camsol, 
                    mylist_fAD12$Zyggregator,
                    mylist_fAD12$Waltz),
               legacy.axes = T, size=1)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color='black'),
    axis.title = element_text(size=20),
    axis.text = element_text(size=15)
  )+
  labs(x="False Positive Rate", y="True Positive Rate")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype="dashed")+
  
  scale_color_manual(labels=labels_plot, values=colors)
g_fAD12

ggsave(g_fAD12, path=path, file="ROC_camsol_ZGG_waltz.pdf", width = 7, height = 4)







### Solubility score

ROC_solubility_df<-solubility_dataset[solubility_dataset$Nmut_codons==1,c("ID", "score")]
colnames(ROC_solubility_df)[2]<-"Solubility score"

ROC_solubility_df$fAD_12<-0
ROC_solubility_df[ROC_solubility_df$ID %in% disease_mutations,]$fAD_12<-1


#build ROC with pROC package

glm.fit=glm(ROC_solubility_df$fAD_12~as.numeric(ROC_solubility_df$`Solubility score`), family = binomial)
roc<-roc(ROC_solubility_df$fAD_12, glm.fit$fitted.values)
mylist_fAD12[["Solubility score"]] <- roc

auc<-roc$auc

auc<-round(auc, 2)

colors=c("black")
labels_plot<-as.character(paste0("Solubility score (AUC=", auc, ")"))


g_fAD12<-ggroc(list(mylist_fAD12$`Solubility score`), legacy.axes = T, size=1)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color='black'),
    axis.title = element_text(size=20),
    axis.text = element_text(size=15)
  )+
  labs(x="False Positive Rate", y="True Positive Rate")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype="dashed")+
  scale_color_manual(labels=labels_plot, values=colors)
g_fAD12

ggsave(g_fAD12, path=path, file="ROC_solubility.pdf", width = 7, height = 4)






####  gnomAD


gnomad<-left_join(gnomad, singles[,c("WT_AA", "Pos", "Mut", "ID", "nscore_c", "sigma", "sig_10")], by="ID")

#create an ID for plot labels
gnomad$label_ID<-paste0(gnomad$WT_AA,gnomad$Pos, gnomad$Mut)



p_gnomad<-ggplot(gnomad)+
              geom_hline(yintercept = 0)+
              
              geom_point(aes(x=freq, y=nscore_c, size=factor(sig_10, levels=c(F,T), labels=c("FDR>0.1", "FDR<0.1")),
                             color=factor(others, levels=c("clinvar- Protective", "clinvar - Likely pathogenic", "fAD", "VUS"), 
                                                  labels=c("Protective (ClinVar)", "Likely Pathogenic (ClinVar)", "Dominant", "VUS"))))+
              
              geom_text_repel(size=6,aes(x=freq, y=nscore_c,label=label_ID, 
                              color=factor(others, levels=c("clinvar- Protective", "clinvar - Likely pathogenic", "fAD", "VUS"), 
                                                  labels=c("Protective (ClinVar)", "Likely Pathogenic (ClinVar)", "Dominant", "VUS"))), show.legend = F)+
              
              labs(x="log10(gnomAD Frequency)", y="Nucleation score")+
              scale_color_manual(values=c("#15983DFF", "#FEC10BFF", "#EE0011FF", "grey60"))+
              theme_bw()+
              theme(axis.title = element_text(size=20),
                    axis.text= element_text(size=16),
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(color='black'),
                    legend.title = element_blank(),
                    legend.text = element_text(size=12),
                    legend.position = c(0.75,0.2))+
              scale_x_continuous(trans = "log10")+
              scale_size_manual(values=c(3,1.5), labels=c("FDR<0.1", "FDR>0.1"))
p_gnomad

ggsave(p_gnomad,path=path, file="p_gnomad.pdf",width=7, height=6, useDingbats=FALSE)


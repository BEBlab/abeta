require(ggplot2)
require(reshape2)
require(dplyr)


dir.create("Charge")
path="Charge"

#required data:
load("nscore_df.RData")


list<-list(
  "Nt (1_26)"=c(2,4:6,8:10,12:16,18:21,24:26),
  "Ct (27_41)"=c(27:41),
  "all (1_42)"=c(1:42),
  "gatekeepers"=c(1,3,7,11,17,22,42)
)


AB.df<-data.frame(seq=c("D","A","E","F","R","H","D","S","G","Y","E","V","H","H","Q","K","L","V","F","F","A",
                        "E","D","V","G","S","N","K","G","A","I","I","G","L","M","V","G","G","V","V","I","A"), 
                  position=c(1:42),
                  charge=0)
AB.df[AB.df$seq %in% c("D", "E"),]$charge<-(-1)
AB.df[AB.df$seq %in% c("R", "K"),]$charge<-1


positive<-c("R", "K")
negative<-c("D", "E")





charge<-singles_doubles

wt_counts=c()

for(region in names(list)){
  
  residues<-list[[region]]
  
  #count the number of positively, negatively charged residues and total net charge of the WT sequence in the specific region
  
  total_pos_wt<-length(which(AB.df[AB.df$position %in% residues,]$charge==1))
  total_neg_wt<-length(which(AB.df[AB.df$position %in% residues,]$charge==(-1)))
  net_charge_wt<-sum(AB.df[AB.df$position %in% residues,]$charge)
  
  wt_counts<-c(wt_counts, region, total_pos_wt, total_neg_wt, net_charge_wt)
  
  tryCatch({
    
    
    #take singles inside the region
    
    singles_charge<-singles_doubles[singles_doubles$Nmut_codons==1,]
    singles_charge<-singles_charge[singles_charge$Pos1 %in% residues,]
    
    
    #count for each mutation the number of positive, negatives residues and net charge
    #this allows to classify each mutation in specific number of positive, negative residues and charge
    all_positive=c()
    all_negative=c()
    all_charge=c()
    
    for( i in 1:nrow(singles_charge)){
      wt_pos<-as.numeric(singles_charge[i,]$WT_AA1 %in% positive)
      mut_pos<-as.numeric(singles_charge[i,]$Mut1 %in% positive)
      total_pos<-as.numeric(total_pos_wt)-wt_pos+mut_pos
      all_positive=c(all_positive, total_pos)
      
      wt_neg<-as.numeric(singles_charge[i,]$WT_AA1 %in% negative)
      mut_neg<-as.numeric(singles_charge[i,]$Mut1 %in% negative)
      total_neg<-as.numeric(total_neg_wt)-wt_neg+mut_neg
      all_negative=c(all_negative, total_neg)
      
      f_charge<-(-abs(total_neg))+total_pos
      all_charge<-c(all_charge, f_charge)
      
    }
    
    singles_charge[[paste0(region, "-positive")]]<-all_positive
    singles_charge[[paste0(region, "-negative")]]<-all_negative
    singles_charge[[paste0(region, "-charge")]]<-all_charge
    
    
    
    #take doubles inside the region and do the same
    
    doubles_charge<-singles_doubles[singles_doubles$Nmut_codons==2,]
    doubles_charge<-doubles_charge[doubles_charge$Pos1 %in% residues &  doubles_charge$Pos2 %in% residues  ,]
    
    
    #count for each mutation the number of positive, negatives residues and net charge
    #this allows to classify each mutation in specific number of positive, negative residues and charge
    
    all_positive=c()
    all_negative=c()
    all_charge=c()
    
    for( i in 1:nrow(doubles_charge)){
      wt_1_pos<-as.numeric(doubles_charge[i,]$WT_AA1 %in% positive)
      wt_2_pos<-as.numeric(doubles_charge[i,]$WT_AA2 %in% positive)
      
      mut_1_pos<-as.numeric(doubles_charge[i,]$Mut1 %in% positive)
      mut_2_pos<-as.numeric(doubles_charge[i,]$Mut2 %in% positive)
      
      total_pos<-as.numeric(total_pos_wt)-wt_1_pos-wt_2_pos+mut_1_pos+mut_2_pos
      all_positive=c(all_positive, total_pos)
      
      wt_1_pos<-as.numeric(doubles_charge[i,]$WT_AA1 %in% negative)
      wt_2_pos<-as.numeric(doubles_charge[i,]$WT_AA2 %in% negative)
      
      mut_1_pos<-as.numeric(doubles_charge[i,]$Mut1 %in% negative)
      mut_2_pos<-as.numeric(doubles_charge[i,]$Mut2 %in% negative)
      
      total_neg<-as.numeric(total_neg_wt)-wt_1_pos-wt_2_pos+mut_1_pos+mut_2_pos
      all_negative=c(all_negative, total_neg)
      
      
      f_charge<-(-abs(total_neg))+total_pos
      all_charge<-c(all_charge, f_charge)
      
      
    }
    
    doubles_charge[[paste0(region, "-positive")]]<-all_positive
    doubles_charge[[paste0(region, "-negative")]]<-all_negative
    doubles_charge[[paste0(region, "-charge")]]<-all_charge
    
    
    subset_region<-rbind(singles_charge, doubles_charge)
    charge<-left_join(charge, subset_region[,c(21,31,32,33)], by="ID")
    
  }, error=function(e){})
}




#melt all added info of charge and split info of region and variants (positive/negative/all)
melt_charge<-melt(charge, id=c(1:30))

region=c()
variants=c()

for (i in 1:nrow(melt_charge)){
  
  region=c(region, unlist(strsplit(as.character(melt_charge[i,]$variable), "-"))[1])
  variants=c(variants, unlist(strsplit(as.character(melt_charge[i,]$variable), "-"))[2])
}

melt_charge$region<-region
melt_charge$variants<-variants

#save(melt_charge, file="melt_charge.RData")





#### arrange table by merging same variants
melt_positive<-melt_charge[melt_charge$variants=="positive",]
colnames(melt_positive)[32]<-"value_positive"
melt_negative<-melt_charge[melt_charge$variants=="negative",]
colnames(melt_negative)[c(32)]<-"value_negative"

melt_pos_neg<-cbind(melt_positive, melt_negative[32])


#remove NA (coming from doubles not assigned)
melt_pos_neg<-melt_pos_neg[!is.na(melt_pos_neg$value_negative),]
melt_pos_neg<-melt_pos_neg[!is.na(melt_pos_neg$value_positive),]



### no aliphatics no aromatics no proline

remove<-c("A", "V", "L", "M", "I", "F", "Y", "W", "P")
subset_melt_pos_neg<-melt_pos_neg[!melt_pos_neg$Mut1%in% remove, ]
subset_melt_pos_neg<-subset_melt_pos_neg[!subset_melt_pos_neg$Mut2%in% remove, ]
subset_melt_pos_neg<-subset_melt_pos_neg[!subset_melt_pos_neg$WT_AA1%in% remove, ]
subset_melt_pos_neg<-subset_melt_pos_neg[!subset_melt_pos_neg$WT_AA2%in% remove, ]





for(region in names(list)){
  
  negative_grids<-as.vector(unique(subset_melt_pos_neg[subset_melt_pos_neg$region==region,]$value_negative))
  negative_grids<-sort(negative_grids, decreasing = F)
  
  positive_grids<-as.vector(unique(subset_melt_pos_neg[subset_melt_pos_neg$region==region,]$value_positive))
  positive_grids<-sort(positive_grids, decreasing = F)
  
  
  my_vector=c()
  for(N in negative_grids[1]:negative_grids[length(negative_grids)]){
    for (P in positive_grids[1]:positive_grids[length(positive_grids)]){
      
      print_charge<-P-N
      my_vector=c(my_vector, print_charge, N, P)
    }}
  
  
  data_text<-as.data.frame(matrix(data=my_vector, ncol=3, byrow = T))
  colnames(data_text)<-c("label","value_negative", "value_positive" )
  data_text$region<-region
  data_text$value_negative<-as.factor(data_text$value_negative)
  data_text$value_positive<-as.factor(data_text$value_positive)
  
  
  subset_melt_pos_neg$color_grid<-paste(subset_melt_pos_neg$value_positive-subset_melt_pos_neg$value_negative)
  subset_melt_pos_neg$total_residues<-paste(subset_melt_pos_neg$value_negative + subset_melt_pos_neg$value_positive)
  
  
  #plot negatively charged residues vs positively charged residues
  
  p1<-ggplot(subset_melt_pos_neg[subset_melt_pos_neg$region==region,], aes(x=factor(region), y=nscore_c))+
    geom_hline(yintercept = 0, color="#15983DFF", linetype="dashed", size=0.3)+
    
    facet_grid(vars(value_negative), vars(value_positive), switch = "both")+
    
    geom_violin(data=subset_melt_pos_neg[subset_melt_pos_neg$region==region,],  
                width=1, size=0.3, aes(fill=as.numeric(color_grid)), show.legend = F)+
    geom_boxplot(width=0.1, outlier.shape = NA, size=0.3)+
    theme_bw()+
    labs(title=paste0(region), y="Negatively charged residues", x="Positively charged residues")+
    theme(axis.text.x = element_blank(),
          panel.grid.major = element_line(size=0.2),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=15),
          legend.title=element_text(),
          panel.border = element_rect(color="grey", size=0.2),
          panel.background = element_rect(fill="white"),
          strip.background = element_rect(fill="white", size=0.2),
          title = element_text(size=15),
          panel.spacing = unit(0, "lines") )+
    scale_y_continuous(position = "right")+
    scale_fill_gradient(high="grey40",low="grey90")+
    scale_color_gradient2(name="total\nnet charge", midpoint=0, low="firebrick", high="darkblue", mid="white")+
    geom_text(data=data_text, aes(label=label, x=-Inf, y=Inf),hjust=-0.15,vjust=1.5, size=4, colour="black")
    
  p1
  
ggsave(p1,path=path, file = paste0(region, "_positive_negative.pdf"), width = 6, height = 5)
  
  
  
  color_grid<-as.vector(unique(subset_melt_pos_neg[subset_melt_pos_neg$region==region,]$color_grid))
  color_grid_order<-as.character(sort(as.numeric(color_grid), decreasing = F))
  
  total_residues<-as.vector(unique(subset_melt_pos_neg[subset_melt_pos_neg$region==region,]$total_residues))
  total_residues_order<-as.character(sort(as.numeric(total_residues), decreasing = F))
  
  text<-as.data.frame(unique(subset_melt_pos_neg[subset_melt_pos_neg$region==region,][,c(32,35,36,37)]))
  text$text<-paste0("+",text$value_positive, "-",text$value_negative)  
  
  
  #plot number of charged resiudes vs total net charge
  
p2<-ggplot(subset_melt_pos_neg[subset_melt_pos_neg$region==region,], aes(x=factor(region), y=nscore_c))+
    geom_hline(yintercept = 0, color="#15983DFF", linetype="dashed", size=0.3)+

    facet_grid(factor(total_residues, levels=total_residues_order)~factor(color_grid, levels=color_grid_order),
               switch = "both", drop=T)+
    
    geom_violin(data=subset_melt_pos_neg[subset_melt_pos_neg$region==region,], width=1,
                aes(fill=as.numeric(total_residues)), size=0.3, show.legend = F)+
    
    
    geom_boxplot(width=0.1, outlier.shape = NA, size=0.3)+
    theme_bw()+
    labs(title=paste0(region), x="Total net charge", y="Number of charged residues")+
    theme(axis.text.x = element_blank(),
          panel.grid.major = element_line(size=0.2),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=15),
          legend.title=element_text(),
          panel.border = element_rect(color="grey", size=0.2),
          panel.background = element_rect(fill="white"),
          strip.background = element_rect(fill="white", size=0.2),
          title = element_text(size=15),
          panel.spacing = unit(0, "lines") )+
    annotate(geom = 'segment', y = Inf, yend = Inf, color = 'grey20', x = -Inf, xend = Inf, size = 0.7)+
    geom_text(data=text, aes(label=text, x=Inf, y=Inf),vjust=1.5, hjust=1,size=3, colour="black")+
    
    scale_fill_gradient(low="grey40",high="grey90")+
    scale_y_continuous(position = "right")
p2
  
  ggsave(p2,path=path, file = paste0(region, "_charged residues vs total charge.pdf"), width = 6, height = 5)

}




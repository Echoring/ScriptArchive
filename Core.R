library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(ComplexHeatmap)


pav <- read.delim("Orthogroups_PAV.tsv", header = T, row.names = 1)
pav <- pav[rowSums(pav)>0,]
samples <- colnames(pav)
indNum <- ncol(pav)
pav[pav > 0] = 1
Ftype <- c("Private" , "Dispensable" ,"Softcore" ,"Core")

Fcut <- c(0, 1 , indNum - 6, indNum - 1,indNum)
#Fcolor <- brewer.pal(length(Ftype), 'Set1')
Fcolor <-c("#BC3C29FF","#0072B5FF","#20854EFF","#E18727FF")

psum <- tibble( ID=rownames(pav), Freq = rowSums(pav)) %>%
mutate( Class  = cut(Freq, breaks = Fcut , labels = Ftype))
pre<- "resultSim0"
write.table(psum,
file=paste(pre, "freq_class.txt", sep = "."),
row.names = F,
sep = "\t",
quote = F)

ggplot(psum,aes( x=Freq) )+geom_bar(aes(fill = Class),width = 0.6) +ylab("Number of Gene/Gene family") +xlab("Frequency") +  scale_fill_manual(values = Fcolor, name = NULL,) +theme_classic()

pfreq <- as.data.frame(table(psum$Class)) %>%
mutate(percent = Freq/sum(Freq)) %>%
dplyr::rename(Class=Var1)
write.table(pfreq,
file=paste(pre, "freq_class.summary.txt", sep = "."),
row.names = F,
sep = "\t",
quote = F)

ggplot(pfreq, aes(x=1 ,y=percent,fill=Class )) +geom_col(color = "white",position = "stack",width = 1) +geom_text( aes(x = 1.3,label = paste0(Class, "\n", round(percent*100,1), "%")) ,position = position_stack(vjust = 0.5)) +scale_fill_manual(values = Fcolor,name = NULL,) +coord_polar(theta = "y",clip = "off")+theme(legend.position = "none") +theme_void()

psum <- psum[order(psum$Freq,decreasing=T),]

names(Fcolor) <- Ftype

pav_sort <- t(pav[psum$ID, ])

Heatmap(pav_sort,cluster_rows = F,cluster_columns = F,col = c("grey","darkorange" ),row_names_side = "left",show_column_names = F ,show_row_names = T,heatmap_legend_param = list(title = " ", at = c(1, 0),labels = c( "Present","Absent"),ncol = 2),top_annotation = HeatmapAnnotation(Type = psum$Class,col = list( Type =  Fcolor)))

pav1 <- read.delim("Orthogroups_PAV.tsv", header = T, row.names = 1)
pav1$type<-psum$Class

typestat<-aggregate(pav1[,1:46],by=list(type=pav1$type),sum)
write.table(typestat,"typestat.txt")
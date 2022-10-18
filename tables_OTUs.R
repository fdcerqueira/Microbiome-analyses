
library(ggplot2)
library(cowplot)
library(randomcoloR)
library(dplyr)
library(scales)

##read tables
dir.create("~/Desktop/tables")
setwd("~/Desktop/tables")

#function for 10^x annotation in scale
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#create function to change rownames
change_rownames <- function(df1) {
  row.names(df1) <- df1[,1]
  df1 <- df1[,-1]
  df1  
}

otu<-readxl::read_xlsx("~/Downloads/Total loads - qPCR and ASVs(1).xlsx", sheet = 2)
dilution<-readxl::read_xlsx("~/Downloads/Total loads - qPCR and ASVs(1).xlsx", sheet = 1)
dilution<-dilution[-c(13:24)]
nn<-dilution[1,]
dilution<-dilution[-1,]
names(dilution)<-nn
dilution<-as.data.frame(dilution)
rownames(dilution)<-dilution[,1]
dilution<-dilution[,-1]

#remove NAs
dilution[is.na(dilution)]<-0

#subset grup SD to get median
dil_SD<-dilution$`Pre-antiB treatment`[grepl("SD*", rownames(dilution))]
dil_SD<-as.numeric(dil_SD)

#use median of SD group to fill the original NA value
dilution$`Pre-antiB treatment`[dilution$`Pre-antiB treatment`==0]<-median(dil_SD)

#separate otu tables
otu<-as.data.frame(otu)

otu_pre<-otu[grepl("Pre_*",otu$index),]
otu_d0<-otu[grepl("D0_*",otu$index),]
otu_d27<-otu[grepl("D27_*",otu$index),]
otu_d34<-otu[grepl("D34_*",otu$index),]
otu_d41<-otu[grepl("D41_*",otu$index),]
otu_d48<-otu[grepl("D48_*",otu$index),]
otu_d55<-otu[grepl("D55_*",otu$index),]
otu_d62<-otu[grepl("D62_*",otu$index),]
otu_d69<-otu[grepl("D69_*",otu$index),]
otu_d76<-otu[grepl("D76_*",otu$index),]
otu_d84<-otu[grepl("D84_*",otu$index),]

#check which sample is missing 
mis<-otu_d62$index
mis<-gsub("D62_","",mis)
mis2<-rownames(dilution)

mis2[!(mis2 %in% mis)]

#index as rowname
otu_data<-list(otu_pre,otu_d0,otu_d27,otu_d34,otu_d41,otu_d48,otu_d55,otu_d62,otu_d69,otu_d76,otu_d84)

otu_names<-c("pre","0","27","34","41","48","55","62","69","76","84")
names(otu_data)<-otu_names

#change colnames dilution
na<-c("pre","0","27","34","41","48","55","62","69","76","84")
names(dilution)<-na

#apply rownames function
otu_data2<-lapply(otu_data, function(x) change_rownames(x))

#apply function transvere
otu_data2<-lapply(otu_data2, function(x) t(x))

#remove extra row
otu_data3<-list()
for (i in 1:length(otu_data2)){
  otu_data3[[i]]<-otu_data2[[i]][-nrow(otu_data2[[i]]),]
}
names(otu_data3)<-names(otu_data2)

#convert to dataframe
otu_data3<-lapply(otu_data3, function(x) as.data.frame(x))

#relative abundances: for loop, could not use lapply
relative<-list()

for (i in names(otu_data3)){                       
  
  relative[[i]]<- lapply(otu_data3[[i]], function(x) x/sum(x, na.rm=TRUE))
}

#convert to dataframe
relative2<-lapply(relative, function(x) as.data.frame(x))

#convert rownames
for (i in 1:length(otu_data3)){
  rownames(relative2[[i]])<-rownames(otu_data3[[i]])
}

#names
names(relative2)<-names(otu_data3)

#remove WD6 sample from day 62 (one missing from otu table day 62)
d61<-dilution$`62`
d61<-as.data.frame(d61)
rownames(d61)<-rownames(dilution)
d61<-d61[rownames(d61) != "WD6",]
dilution<-dilution[-which(names(dilution) %in% "62")]

#extract and remove 62 dataframe from list
d62_new<-relative2[["62"]]
relative3<-relative2[names(relative2) != "62" ]

#empty vector and list
final_list<-list()
finall<-list()
dil_col<-c()

#empty vector
dil_col<-c()
rowss<-c()

#multiply values
for (i in 1:ncol(dilution)){
  for (j in 1:length(relative3)){
    for (k in 1:nrow(relative3[[j]])){
      
      #statement to make sure the correct column of 16s table is going to multiply in with the right 
      if (names(relative3)[j]==names(dilution)[i]){
        
        #save vales of 16s abundance columns and rows from otu data to numeric vectors
        dil_col<-as.numeric(dilution[,i])
        rowss<-as.numeric(relative3[[j]][k,])
        #print(rowss)
        print(paste("multipling 16s day:", names(dilution)[i],",row:",k, ",dataframe:",names(relative3)[j]))
        #multly vector with 16s data (each column- correspond to all samples 1 timepoint ) X relative abundance OTUs (each dataframe-1 time point, per row-all samples) 
        final_list[[k]]<-dil_col*rowss
        #each time point list to be stored in final list
        finall[[j]]<-final_list
        
      }else{break}
    }
  }    
} 
#convert nested.list to dataframe and give names to lists
final_tables <- lapply(finall, function(x) as.data.frame(do.call(rbind, x)))
names(final_tables)<-names(relative3)

#apply to final table rownames and colnames
for (i in 1:length(relative3)){
  names(final_tables[[i]])<-names(relative3[[i]])
  rownames(final_tables[[i]])<-rownames(relative3[[i]])
}

#convert d62 to numeric matrix and 16s values to numeric vector
d62_final<-apply(as.matrix(d62_new), 2, as.numeric)
rownames(d62_final)<-rownames(d62_new)
names(d62_final)<-names(d62_new)
d61<-as.vector(as.numeric(d61))

#multiply 16s values with rows day62
d62_final1<-sweep(d62_final, MARGIN=2, d61, `*`)

#merge d62 to others
final_tables[[length(final_tables)+1]]<-d62_final1
names(final_tables)[11]<-"62"

#write tables
for (i in 1:length(final_tables)){
  write.csv(final_tables[[i]], paste(names(final_tables)[i],".csv"))
}

#create new nested list
table_p<-final_tables[names(final_tables) != "62" ]
d62_g1<-final_tables[names(final_tables) == "62" ]
#remove d62
table_p1<-table_p

#change colnames
cenas1<-names(final_tables[[1]])
cenas1<-gsub("Pre_antiB_treatment_","",cenas1)

#give new colnames each table add otu column, separate otu column
for (i in 1:length(table_p1)) {
  names(table_p1[[i]])<-cenas1
  table_p1[[i]]$otu<-rownames(table_p1[[i]])
  table_p1[[i]] <- table_p1[[i]] %>% separate(otu, sep=";",c("Kingdom","Phylum","Class","Order","Family","Genus"))  
}

#put proper name for genus
for (i in 1:length(table_p1)){
  for ( j in 1:nrow(table_p1[[i]])){
    
    if(str_detect(table_p1[[i]][j,]$Genus,"D_5__unidentified$") ==T)
    {
      table_p1[[i]][j,]$Genus<-table_p1[[i]][j,]$Family
    } 
    else if (str_detect(table_p1[[i]][j,]$Genus,"D_5__uncultured$")==T )
    {
      table_p1[[i]][j,]$Genus<-table_p1[[i]][j,]$Family
    }
    else if (str_detect(table_p1[[i]][j,]$Genus,"D_5__uncultured[:space:]bacterium$")==T )
    {
      table_p1[[i]][j,]$Genus<-table_p1[[i]][j,]$Order 
    }
    else if (str_detect(table_p1[[i]][j,]$Family,".*__$") ==F & 
             str_detect(table_p1[[i]][j,]$Genus,".*__$") ==T)
    { 
      if (str_detect(table_p1[[i]][j,]$Family,"D_4__uncultured[:space:]rumen[:space:]bacterium$") ==T)
      {
        table_p1[[i]][j,]$Genus<-table_p1[[i]][j,]$Order
      }
      else if ((str_detect(table_p1[[i]][j,]$Family,"D_4__uncultured[:space:]bacterium$") ==T))
      {
        table_p1[[i]][j,]$Genus<-table_p1[[i]][j,]$Order
      }
      else  
      {
        table_p1[[i]][j,]$Genus<-table_p1[[i]][j,]$Family
      }
    } 
    else if (str_detect(table_p1[[i]][j,]$Order,".*__$") ==F & 
             str_detect(table_p1[[i]][j,]$Family,".*__$") ==T )
    {
      table_p1[[i]][j,]$Genus<-table_p1[[i]][j,]$Order 
    }
    else if (str_detect(table_p1[[i]][j,]$Class,".*__$") ==F & 
             str_detect(table_p1[[i]][j,]$Order,".*__$") ==T )
    {
      table_p1[[i]][j,]$Genus<-table_p1[[i]][j,]$Class
    }
  }
}

#remove other taxa levels,prefix of Genus, and make each genus unique due to some genus have the same name (taxonomy atribution at the lowest level possible)
for (i in 1:length(table_p1[[i]])){
  table_p1[[i]]<-table_p1[[i]][-which(names(table_p1[[i]]) %in% c("Kingdom","Phylum","Class","Order","Family"))]
  table_p1[[i]]$Genus<-gsub("^D_\\d__", "",table_p1[[i]]$Genus)
  rownames(table_p1[[i]])<-make.unique(as.character(table_p1[[i]]$Genus), sep = "_")
}

#add colour column to each dataframe
# no. of colours in the palette
no_of_colors <- length(table_p1[[1]]$Genus)

#pie(rep(1, no_of_colors), col = palette, main = "palette using randomcoloR package")
# sample colors
palette <- distinctColorPalette(no_of_colors) 

#create new column with colors
for (i in 1:length(table_p1)){
  table_p1[[i]]$colorz<-palette
}

tabl<-list()
#create otu identifier column, just keep name of genus and convert it to long format
for (i in 1:length(table_p1)){
  table_p1[[i]]$Genus<-rownames(table_p1[[i]])
  tabl[[i]]<-table_p1[[i]] %>% gather(sample, abundance, -Genus, -colorz)
}

names(tabl)<-names(table_p1)

tabl1<-list()

#create column for Groups
for ( i in 1:length(tabl)){
  tabl[[i]]$group<-ifelse(grepl("SD*",tabl[[i]]$sample),"SD",
                          ifelse(grepl("WD*",tabl[[i]]$sample),"WD","AD"))
  tabl[[i]]$sample<- factor(tabl[[i]]$sample,levels = 
                              c("SD1","SD2","SD3","SD4","SD5", "SD6","SD7","SD8","SD9","SD10","SD11","SD12","SD13","SD14",
                                "WD1","WD2","WD3","WD4","WD5","WD6","WD7","WD8","WD9","WD10","WD11","WD12","WD13","WD14",
                                "AD1","AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12","AD13","AD14"))
} 

abun_phylum<-list()

#graphs
for (i in 1:length(tabl)){
  
  #get tables and create others group
  abun_phylum[[i]]<- as.data.frame(tabl[[i]])
  names(abun_phylum)[i]<-names(tabl)[i]
  abun_phylum[[i]]$Genus<- as.character(abun_phylum[[i]]$Genus)
  abun_phylum[[i]]$Genus[abun_phylum[[i]]$abundance < 2000000000 ] <- "Others"
  abun_phylum[[i]]$colorz[abun_phylum[[i]]$Genus=="Others" ] <-"black"
  #create a vector with colors, and give each element a name corresponding to genus 
  col <-as.character(abun_phylum[[i]]$colorz)
  names(col) <- as.character(abun_phylum[[i]]$Genus)
  #reorder(sample, -abundance, sum) to reorder by abundance
  #graph
  xx<-ggplot(abun_phylum[[i]], aes(sample, abundance, fill = factor(Genus)))+
    geom_bar(stat="identity", colour="black")+
    scale_fill_manual(values=col)+
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1))+
    ylab("Genera abundance \n")+
    xlab("Samples")+
    facet_wrap(~group, scales = "free")+
    theme_cowplot(font_size = 13)+
    scale_y_continuous(expand = c(0, 0), label=scientific_10)+
    theme(axis.text.x = element_text(angle = 0,vjust = 1.2), 
          panel.grid.major = element_line(colour = "white"),)+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
    theme(plot.margin = margin(20,20,20,20, "points"))+
    theme(strip.background = element_blank())+
    labs(title = NULL)+
    guides(fill=guide_legend("Genera",ncol=1))
  ggsave(paste(names(tabl)[i],".pdf"),height = 9, width = 15)
}

#day 62 graph
d62_g1<-as.data.frame(final_tables[names(final_tables) == "62" ])

#change colnames
cenas2<-names(d62_g1)
cenas2<-gsub("X62.D62_","",cenas2)
names(d62_g1)<-cenas2
d62_g1$otu<-rownames(d62_g1)
d62_g1 <- d62_g1 %>% separate(otu, sep=";",c("Kingdom","Phylum","Class","Order","Family","Genus"))  

#put proper name for genus
for ( j in 1:nrow(d62_g1)){
  
  if(str_detect(d62_g1[j,]$Genus,"D_5__unidentified$") ==T)
  {
    d62_g1[j,]$Genus<-d62_g1[j,]$Family
  } 
  else if (str_detect(d62_g1[j,]$Genus,"D_5__uncultured$")==T )
  {
    d62_g1[j,]$Genus<-d62_g1[j,]$Family
  }
  else if (str_detect(d62_g1[j,]$Genus,"D_5__uncultured[:space:]bacterium$")==T )
  {
    d62_g1[j,]$Genus<-d62_g1[j,]$Order 
  }
  else if (str_detect(d62_g1[j,]$Family,".*__$") ==F & 
           str_detect(d62_g1[j,]$Genus,".*__$") ==T)
  { 
    if (str_detect(d62_g1[j,]$Family,"D_4__uncultured[:space:]rumen[:space:]bacterium$") ==T)
    {
      d62_g1[j,]$Genus<-d62_g1[j,]$Order
    }
    else if ((str_detect(d62_g1[j,]$Family,"D_4__uncultured[:space:]bacterium$") ==T))
    {
      d62_g1[j,]$Genus<-d62_g1[j,]$Order
    }
    else  
    {
      d62_g1[j,]$Genus<-d62_g1[j,]$Family
    }
  } 
  else if (str_detect(d62_g1[j,]$Order,".*__$") ==F & 
           str_detect(d62_g1[j,]$Family,".*__$") ==T )
  {
    d62_g1[j,]$Genus<-d62_g1[j,]$Order 
  }
  else if (str_detect(d62_g1[j,]$Class,".*__$") ==F & 
           str_detect(d62_g1[j,]$Order,".*__$") ==T )
  {
    d62_g1[j,]$Genus<-d62_g1[j,]$Class
  }
}


#remove other taxa levels,prefix of Genus, and make each genus unique due to some genus have the same name (taxonomy atribution at the lowest level possible)
d62_g1<-d62_g1[-which(names(d62_g1) %in% c("Kingdom","Phylum","Class","Order","Family"))]
d62_g1$Genus<-gsub("^D_\\d__", "",d62_g1$Genus)
rownames(d62_g1)<-make.unique(as.character(d62_g1$Genus), sep = "_")
d62_g1$colorz<-palette

#create otu identifier column, just keep name of genus and convert it to long format
d62_g1$Genus<-rownames(d62_g1)
d62_g11<-d62_g1 %>% gather(sample, abundance, -Genus, -colorz)

##create group column
d62_g11$group<-ifelse(grepl("SD*",d62_g11$sample),"SD",
                      ifelse(grepl("WD*",d62_g11$sample),"WD","AD"))

#d62
library(data.table)
long <- melt(setDT(d62_g1), id.vars = c("Genus","colorz"), variable.name = "abundance")

teste1<-long%>%group_by(abundance)%>%
  summarize(Mean = mean(value, na.rm=TRUE))


#get tables and create others group
abun_phylum<- as.data.frame(d62_g11)
abun_phylum$Genus<- as.character(abun_phylum$Genus)
abun_phylum$Genus[abun_phylum$abundance < 2000000000 ] <- "Others"
abun_phylum$colorz[abun_phylum$Genus=="Others" ] <-"black"


abun_phylum$sample<- factor(abun_phylum$sample,levels = 
                              c("SD1","SD2","SD3","SD4","SD5", "SD6","SD7","SD8","SD9","SD10","SD11","SD12","SD13","SD14",
                                "WD1","WD2","WD3","WD4","WD5","WD7","WD8","WD9","WD10","WD11","WD12","WD13","WD14",
                                "AD1","AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12","AD13","AD14"))


#create a vector with colors, and give each element a name corresponding to genus 
col <-as.character(abun_phylum$colorz)
names(col) <- as.character(abun_phylum$Genus)

#reorder(sample, -abundance, sum) ro reorder by abundance
#graph
xx<-ggplot(abun_phylum, aes(sample, abundance, fill = factor(Genus)))+
  geom_bar(stat="identity", colour="black")+
  scale_fill_manual(values=col)+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1))+
  ylab("Genera abundance \n")+
  xlab("Samples")+
  facet_wrap(~group, scales = "free")+
  theme_cowplot(font_size = 13)+
  #scale_y_continuous(expand = c(0, 0), label=scientific_10)+
  theme(axis.text.x = element_text(angle = 0,vjust = 1.2), 
        panel.grid.major = element_line(colour = "white"),)+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
  theme(plot.margin = margin(20,20,20,20, "points"))+
  theme(strip.background = element_blank())+
  labs(title = NULL)+
  guides(fill=guide_legend("Genera",ncol=1))
ggsave("~/Desktop/62.pdf",height = 9, width = 15)

#add d62 to list
tabl2<-tabl
tabl2[["62"]]<-d62_g11
tabl2[["pre"]]<-tabl2[["pre"]][-6]


#create column days for each dataframe 
for (i in 1:length(tabl2)){
  tabl2[[i]]$days<-rep(paste(names(tabl2)[i]),nrow(tabl2[[i]]))
}

time_series<-data.frame()
#merge by row the dataframes
for (i in 1:length(tabl2)){
  time_series<-rbind(time_series,tabl2[[i]])
}

#average by group at each day, and standard deviation and standard error
#df1 <- time_series %>% group_by(group,days,Genus) %>% 
#  summarise(mean = mean(abundance),
#            sd=sd(abundance),
#            n=n(),
#            se= sd(abundance)/sqrt(n()))

df2<-time_series %>% group_by(group,days) %>% 
  summarise(mean = mean(abundance),
            sd=sd(abundance),
            n=n(),
            se= sd(abundance)/sqrt(n()))

#factor order (samples)
df$days<-factor(df$days,levels = c("pre","0","27","34","41","48","55","62","69","76","84"))

#line plot
ggplot(df2, aes(x=days, y=mean, , group=group, fill= group, shape=group))+
  geom_errorbar(aes(ymin = df2$mean-df2$se,
                    ymax = df2$mean+df2$se, width=0.3,color=group))+
  geom_line(aes(color=group))+
  geom_point(size=4)+
  scale_y_continuous(expand = c(0, 0),limits = c(0, NA), label=scientific_10)+
  scale_shape_manual(values=c(22,21,23))+
  xlab("Days")+
  ylab("Mean abundance")+
  facet_wrap(~group, scales = "free")+
  theme_cowplot(font_size = 13)+
  theme(strip.background = element_blank())
ggsave("~/Desktop/lines.pdf",height = 9, width = 15)

#barplot
ggplot(df2, aes(x=days, y=mean,group=group,fill= group))+
  geom_bar(stat="identity", colour="black")+
  geom_errorbar(aes(ymin = df2$mean-df2$se,
                    ymax = df2$mean+df2$se, width=0.3))+
  xlab("Days")+
  ylab("Mean abundance")+
  scale_y_continuous(expand = c(0, 0), label=scientific_10)+
  facet_wrap(~group, scales = "free")+
  theme_cowplot(font_size = 13)+
  theme(strip.background = element_blank())
ggsave("~/Desktop/barplot.pdf",height = 9, width = 15)

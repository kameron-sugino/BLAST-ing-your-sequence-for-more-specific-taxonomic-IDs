#load packages
require(ggplot2)
require(colorspace)
require(vegan)
require(dplyr)

#need to format run 1 so that the runs arent averaged together
#read in data, rename 5A as 3A
tax<-read.table("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/QIIME_gg2_annotations/taxonomy.tsv",sep="\t",header=T)
df.table<-read.table("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/QIIME_gg2_annotations/table.from_biom.txt",header=F)

id<-df.table[1,-1]

test<-merge(tax,df.table,by.x="Feature.ID",by.y="V1")
tax.n<-test$Taxon

otu.id<-data.frame(t(id),t(test[,-c(1:3)]))
colnames(otu.id)<-c("id",tax.n)

#rename 5A to 3A
otu.id[32,1]<-"3A"

#split data for later use
otu<-sapply(otu.id[,-1], as.numeric)
id<-otu.id$id
tax<-colnames(otu)

#sort lists
# otu.ag<-aggregate(x = otu, by = list(gsub("A|B|C|D","",id)), FUN = mean)
ra<-otu/rowSums(otu)*100
ra<-data.frame(id,ra)

#transpose ra to sum columns with the same species designation
ra.t<-t(ra[,-1])
taxa.t<-ifelse(nchar(gsub(".*s__","",rownames(ra.t)))<5,gsub(".*f__","",rownames(ra.t)),gsub(".*s__","",rownames(ra.t)))
taxa.t<-gsub("\\.[0-9]|[0-9]","",taxa.t)

ra.t.ag<-aggregate(ra.t,list(taxa.t),sum)

test<-data.frame(ra[,1],t(ra.t.ag[,-1]))
colnames(test)<-c("id",ra.t.ag[,1])
ra<-test

#reannotating using BLAST results
#need to combine L.helveticus/L.acidophilus (keep name as acidophilus for now, rename later), crispatus also likely acidophilus
ra$Lactobacillus.acidophilus<-ra$Lactobacillus.acidophilus+ra$Lactobacillus.helveticus+ra$Lactobacillus.crispatus
ra<-ra[,-which(colnames(ra)%in%c("Lactobacillus.helveticus","Lactobacillus.crispatus"))]

#these sequences map to reuteri
ra$Limosilactobacillus.reuteri<-ra$Lactobacillaceae..g__Limosilactobacillus..s__+ra$Limosilactobacillus.reuteri+ra$Limosilactobacillus.vaginalis
ra<-ra[,-which(colnames(ra)%in%c("Lactobacillaceae..g__Limosilactobacillus..s__","Limosilactobacillus.vaginalis"))]

#these sequences map to thermophilus
ra$Streptococcus.thermophilus<-ra$Streptococcus.thermophilus+ra$Streptococcaceae..g__Streptococcus..s__
ra<-ra[,-which(colnames(ra)=="Streptococcaceae..g__Streptococcus..s__")]

#refactor bacterial names to first letter of the genus and full species name
taxa.name<-gsub("","",colnames(ra[,-1]))
colnames(ra)<-c("id",paste0(sub("(.)[^ ]*", "\\1", taxa.name),".",gsub(".*\\.","",taxa.name)))

#rename animalis to animalis/lactis, longum/infantis, acidophilus/helveticus
colnames(ra)[which(colnames(ra)%in%c("B.animalis","B.longum","L.acidophilus"))]<-c("B.animalis/B.lactis","B.longum/B.infantis","L.acidophilus/L.helveticus")

df.prob.list<-read.csv("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/221006_Probiotics Master List_table.csv")


ra$unique<-ra$id
ra$id<-gsub("A|B|C|D","",ra$unique)

df.m.list<-merge(df.prob.list,ra,by.x="Number",by.y="id")
#remove non-unique identifier
df.m.list$Number<-df.m.list$unique
df.m.list<-df.m.list[,-which(colnames(df.m.list)=="unique")]

ra.sort<-df.m.list[,-c(1:ncol(df.prob.list))]
list.sort<-df.m.list[,c(1:ncol(df.prob.list))]
#make aggregate df
ra.sort.agg<-aggregate(ra.sort, list(list.sort$Brand),mean)



#########
#write.csv(list.sort,"Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/master_list_sorted.csv",row.names = F)

top.ra<-ra.sort[,apply(ra.sort, 2, function(x) max(x, na.rm = TRUE))>=0.5]
bot.ra<-ra.sort[,apply(ra.sort, 2, function(x) max(x, na.rm = TRUE))<0.5]

ra.new<-data.frame(top.ra,rowSums(bot.ra))
colnames(ra.new)<-c(colnames(top.ra),"Other")
ra.new<-data.frame(list.sort[,c(1:2)],ra.new)

write.csv(ra.new,"Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/20231214_ra_unmerged_gg2_BLAST_names_table.csv",row.names = F)

##########################################################
#since we're here, let's make some plots
l<-list.sort[,-c(1:2)]
l.n<-colnames(l)

onoffsum<-data.frame(paste("No",ra.new$Number),top.ra,rowSums(bot.ra))
colnames(onoffsum)<-c("id",colnames(top.ra),"other")

require(ggplot2)
require(reshape2)
onoffsum.melt<-melt(onoffsum, id.vars = 1)
onoffsum.melt$id<-factor(onoffsum.melt$id,levels=unique(onoffsum.melt$id))
store.melt<-onoffsum.melt
id.list<-unique(onoffsum.melt$id)

#gotta make the column then fill it
onoffsum.melt$label<-NA
for(i in 1:length(id.list)){
  num<-which(onoffsum.melt$id==unique(onoffsum.melt$id)[i])
  onoffsum.melt$label[num]<-ifelse(grepl(gsub("\\.","|",paste(gsub(".*\\.\\.","",l.n[l[i,]==T]),collapse = "|")),onoffsum.melt$variable[num],ignore.case = T)
                                   ,"on label","off label")
}

#create aggregate table base on label
temp<-aggregate(round(value,2) ~ id + label,data=onoffsum.melt,FUN=sum)
test<-recast(temp, id ~ label)

#write.csv(test, "Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/plots/20231129_onoff_label_overview.csv",row.names = F)


tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/plots/20231127_unmerged_gg2_BLAST_onoff_label_overview.tif", res=300, height=7, width=10, units="in", compression="zip")
ggplot(onoffsum.melt, aes(fill=label, y=value, x=id)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x ="", y = "% Abundance",fill="Product \nComposition")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_grid(. ~ factor(gsub("A|B|C|D","",id),levels=c("No 1","No 2","No 3","No 4","No 5","No 6","No 7","No 8","No 9","No 10","No 11","No 12","No 14","No 15","No 16","No 17")), drop=TRUE,scale="free",space="free_x")
dev.off()

#order sample number, variables and create colors
require(pals)
onoffsum.melt$id <- factor(onoffsum.melt$id,levels=unique(as.character(onoffsum.melt$id)))
onoffsum.melt$variable <- factor(onoffsum.melt$variable,levels = unique(as.character(onoffsum.melt$variable)))
onoffsum.melt$col<-onoffsum.melt$variable
on.col<-unname(alphabet())[c(1:length(unique(onoffsum.melt$col)))]

tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/plots/20231214_unmerged_gg2_BLAST_onoff_label_details.tif", res=300, height=7, width=10, units="in", compression="zip")
ggplot(onoffsum.melt, aes(fill=variable, y=value, x=label)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = on.col)+
  facet_wrap(~id)+
  labs(x ="", y = "% Abundance",fill="Product \nComposition")
dev.off()


#####
#heatmap showing present/absent/should be present
pres.table<-ifelse(ra.sort>0.5,T,F)
names.table<-gsub("\\.\\.\\..*","",colnames(pres.table))
names.table<-gsub(".*\\.","",names.table)

onoffsum.melt2<-store.melt

onoffsum.melt2$onlabel<-NA
for(i in 1:length(id.list)){
  num<-which(onoffsum.melt2$id==unique(onoffsum.melt2$id)[i]) #grabs the location of each sample in the melt table
  onoffsum.melt2$onlabel[num]<-grepl(gsub("\\.","|",paste(gsub(".*\\.\\.","",l.n[l[i,]==T]),collapse = "|")),onoffsum.melt2$variable[num],ignore.case = T)
}

onoffsum.melt2$onlabel[num]<-grepl(gsub("\\.","|",paste(gsub(".*\\.\\.","",l.n[l[i,]==T]),collapse = "|")),onoffsum.melt2$variable[num],ignore.case = T)

onoffsum.melt2$insample<-ifelse(onoffsum.melt2$value>0.5,T,F)

onoffsum.melt2$onoff.final<-ifelse(onoffsum.melt2$onlabel==T&onoffsum.melt2$insample==T,"on label in sample",ifelse(onoffsum.melt2$onlabel==F&onoffsum.melt2$insample==T,"not on label in sample",ifelse(onoffsum.melt2$onlabel==T&onoffsum.melt2$insample==F,"on label not in sample","not on label not in sample")))


onoffsum.melt2$id <- factor(onoffsum.melt2$id,levels=unique(as.character(onoffsum.melt2$id)))
onoffsum.melt2$variable <- factor(onoffsum.melt2$variable,levels = unique(as.character(onoffsum.melt2$variable)))

levels(onoffsum.melt2$variable)
tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/plots/20231127_unmerged_gg2_BLAST_heatmap.tif", res=300, height=7, width=22, units="in", compression="zip")
ggplot(onoffsum.melt2, aes(y=forcats::fct_rev(variable), x=id, fill=factor(onoff.final))) + 
  geom_tile()+
  geom_text(aes(label = round(value, 1)))+
  xlab("")+
  ylab("")+
  labs(fill="Label vs Sample")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_grid(. ~ factor(gsub("A|B|C|D","",id),levels=c("No 1","No 2","No 3","No 4","No 5","No 6","No 7","No 8","No 9","No 10","No 11","No 12","No 14","No 15","No 16","No 17")), drop=TRUE,scale="free",space="free_x")
dev.off()




###########################################################
#let's get to merging
run1<-read.csv("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/20231214_ra_unmerged_gg2_BLAST_names_table.csv")

run2<-read.csv("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_2/20231214_ra_gg2_BLAST_names_table.csv")

run1<-data.frame(run1[,c(1:2)],"run1",run1[,-c(1:2)])
colnames(run1)<-c(colnames(run1[,c(1:2)]),"run",colnames(run1[,-c(1:3)]))
run_1_agg<-data.frame(rep("run1",16),aggregate(run1[,-c(1:3)],list(run1$Brand),mean))
colnames(run_1_agg)<-c("run","Brand",colnames(run_1_agg[,-c(1:2)]))

run2<-data.frame(run2[,c(1:2)],"run2",run2[,-c(1:2)])
colnames(run2)<-c(colnames(run2[,c(1:2)]),"run",colnames(run2[,-c(1:3)]))
run2$Number<-as.character(run2$Number)

#merge with dplyr
runs_m<-bind_rows(run1,run2)
runs_m[is.na(runs_m)]<-0

runs_m_agg<-bind_rows(run_1_agg,run2[,-1])
runs_m_agg[is.na(runs_m_agg)]<-0

################################
################################
temp_runs<-runs_m

#temp_runs<-runs_m_agg

#refactor id's to order run 1 then run 2 when run 2 had data
merged_runs<-temp_runs[order(as.numeric(gsub("[A-Z]","",temp_runs$Number))),]
#refactor columns so that other is at the end
temp<-data.frame(merged_runs[,c(1:(which(colnames(merged_runs)=="Other")-1))],merged_runs[,c((which(colnames(merged_runs)=="Other")+1):ncol(merged_runs))],merged_runs[,c(which(colnames(merged_runs)=="Other"))])
names(temp)[ncol(temp)]<-"Other"

merged_runs<-temp

#read in each run's list of on/off label and merge
list.sort.r1<-read.csv("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/master_list_sorted.csv")
list.sort.r2<-read.csv("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_2/master_list_sorted.csv")
list.sort.r2$Number<-as.character(list.sort.r2$Number)

#refactor list to match merged_runs
list.temp<-bind_rows(list.sort.r1,list.sort.r2)
list.sort<-list.temp[order(as.numeric(gsub("[A-Z]","",list.temp$Number))),]
list.sort[is.na(list.sort)]<-FALSE

l<-list.sort[,-c(1:2)]
l.n<-colnames(l)

onoffsum<-data.frame(paste("No",merged_runs$Number),merged_runs[,-c(1:3)])
colnames(onoffsum)<-c("id",colnames(merged_runs[,-c(1:3)]))

#mean of onoffsum
df_agg<-aggregate(onoffsum[,-1],list(list.sort$Brand),mean)

require(ggplot2)
require(reshape2)
onoffsum.melt<-melt(onoffsum, id.vars = 1)
onoffsum.melt$id<-factor(onoffsum.melt$id,levels=unique(onoffsum.melt$id))
store.melt<-onoffsum.melt
id.list<-unique(onoffsum.melt$id)

#gotta make the column then fill it
onoffsum.melt$label<-NA
for(i in 1:length(id.list)){
  num<-which(onoffsum.melt$id==unique(onoffsum.melt$id)[i])
  onoffsum.melt$label[num]<-ifelse(grepl(gsub("\\.","|",paste(gsub(".*\\.\\.","",l.n[l[i,]==T]),collapse = "|")),onoffsum.melt$variable[num],ignore.case = T)
                                   ,"on label","off label")
}

#create aggregate table base on label
temp<-aggregate(round(value,2) ~ id + label,data=onoffsum.melt,FUN=sum)
test<-recast(temp, id ~ label)

#write.csv(test, "Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_2/plots/20231129_onoff_label_overview.csv",row.names = F)


tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20231213_run1_run2_gg2_BLAST_onoff_label_overview.tif", res=300, height=7, width=20, units="in", compression="zip")
ggplot(onoffsum.melt, aes(fill=label, y=value, x=id)) + 
  geom_bar(position="stack",  stat="identity")+
  labs(x ="", y = "% Abundance",fill="Product \nComposition")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  
  facet_grid(. ~ factor(gsub("A|B|C|D","",id),
            levels=c("No 1","No 2","No 3","No 4","No 5",
                     "No 6","No 7","No 8","No 9","No 10",
                     "No 11","No 12", "No 13","No 14",
                     "No 15","No 16","No 17")), 
            drop=TRUE,scale="free",space="free_x")+
  
  stat_summary(aes(x = id, y = value),
               fun.y = sum,
               geom = "col",
               colour = ifelse(grepl("A|B|C|D",unique(onoffsum.melt$id)),0,1),
               fill = ifelse(grepl("A|B|C|D",unique(onoffsum.melt$id)),NA,NA))
  # cant get the boxes to outline or insert a shape parameter because I'm not calling the data in an aes statement correctly
  # geom_segment(aes(xend = id, yend = ifelse(grepl("A|B|C|D",id),0,1),linetype=ifelse(grepl("A|B|C|D",id),"Run 1","Run 2")),linewidth=NA, colour = ifelse(grepl("A|B|C|D",onoffsum.melt$id),0,1))+
  # 
  # scale_linetype_manual(values = c(1,2), name = "Sequencing \nRun")

dev.off()

#order sample number, variables and create colors
require(pals)
onoffsum.melt$id <- factor(onoffsum.melt$id,levels=unique(as.character(onoffsum.melt$id)))
onoffsum.melt$variable <- factor(onoffsum.melt$variable,levels = unique(as.character(onoffsum.melt$variable)))
onoffsum.melt$col<-onoffsum.melt$variable
on.col<-unname(alphabet())[c(1:length(unique(onoffsum.melt$col)))]

tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20231213_run1_run2_gg2_BLAST_onoff_label_details.tif", res=300, height=7, width=15, units="in", compression="zip")
ggplot(onoffsum.melt, aes(fill=variable, y=value, x=label)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = on.col)+
  facet_wrap(~id)+
  labs(x ="", y = "% Abundance",fill="Product \nComposition")
dev.off()


#heatmap showing present/absent/should be present
pres.table<-ifelse(ra.sort>0.5,T,F)
names.table<-gsub("\\.\\.\\..*","",colnames(pres.table))
names.table<-gsub(".*\\.","",names.table)

onoffsum.melt2<-store.melt

onoffsum.melt2$onlabel<-NA
for(i in 1:length(id.list)){
  num<-which(onoffsum.melt2$id==unique(onoffsum.melt2$id)[i]) #grabs the location of each sample in the melt table
  onoffsum.melt2$onlabel[num]<-grepl(gsub("\\.","|",paste(gsub(".*\\.\\.","",l.n[l[i,]==T]),collapse = "|")),onoffsum.melt2$variable[num],ignore.case = T)
}

onoffsum.melt2$onlabel[num]<-grepl(gsub("\\.","|",paste(gsub(".*\\.\\.","",l.n[l[i,]==T]),collapse = "|")),onoffsum.melt2$variable[num],ignore.case = T)

onoffsum.melt2$insample<-ifelse(onoffsum.melt2$value>0.5,T,F)

onoffsum.melt2$onoff.final<-ifelse(onoffsum.melt2$onlabel==T&onoffsum.melt2$insample==T,"on label in sample",ifelse(onoffsum.melt2$onlabel==F&onoffsum.melt2$insample==T,"not on label in sample",ifelse(onoffsum.melt2$onlabel==T&onoffsum.melt2$insample==F,"on label not in sample","not on label not in sample")))


onoffsum.melt2$id <- factor(onoffsum.melt2$id,levels=unique(as.character(onoffsum.melt2$id)))
onoffsum.melt2$variable <- factor(onoffsum.melt2$variable,levels = unique(as.character(onoffsum.melt2$variable)))

levels(onoffsum.melt2$variable)
tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20231213_run1_run2_gg2_BLAST_heatmap.tif", res=300, height=7, width=25, units="in", compression="zip")
ggplot(onoffsum.melt2, aes(y=forcats::fct_rev(variable), x=id, fill=factor(onoff.final))) + 
  geom_tile()+
  geom_text(aes(label = round(value, 1)))+
  xlab("")+
  ylab("")+
  labs(fill="Label vs Sample")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  
  facet_grid(. ~ factor(gsub("A|B|C|D","",id),
                        levels=c("No 1","No 2","No 3","No 4","No 5",
                                 "No 6","No 7","No 8","No 9","No 10",
                                 "No 11","No 12", "No 13","No 14",
                                 "No 15","No 16","No 17")), 
             drop=TRUE,scale="free",space="free_x")+
  
  stat_summary(aes(x = id, y = 1.05),
               fun.y = sum,
               geom = "col",
               colour = ifelse(grepl("A|B|C|D",unique(onoffsum.melt$id)),0,1),
               fill = ifelse(grepl("A|B|C|D",unique(onoffsum.melt$id)),NA,NA))
dev.off()




####################
##heatmap of aggregated data (mean)
#heatmap showing present/absent/should be present
pres.table<-ifelse(df_agg[,-1]>0.5,T,F)
names.table<-gsub("\\.\\.\\..*","",colnames(pres.table))
names.table<-gsub(".*\\.","",names.table)

onoffsum.melt2<-melt(df_agg)
id_list<-df_agg$Group.1
temp.l<-aggregate(list.sort[,-c(1:2)],list(list.sort$Brand),mean)
temp.l<-ifelse(temp.l>0,TRUE,FALSE)
l<-temp.l[,-1]
l.n<-colnames(l[,-1])

temp.l$Group.1==id_list

onoffsum.melt2$onlabel<-NA
for(i in 1:length(id_list)){
  num<-which(onoffsum.melt2$Group.1%in%onoffsum.melt2$Group.1[i]) #grabs the location of each sample in the melt table
  comp<-gsub("\\.","|",paste(gsub(".*\\.\\.","",colnames(l)[l[i,]==T]),collapse = "|"))
  onoffsum.melt2$onlabel[num]<-grepl(comp,onoffsum.melt2$variable[num],ignore.case = T)
}

onoffsum.melt2$insample<-ifelse(onoffsum.melt2$value>0.5,T,F)

onoffsum.melt2$onoff.final<-ifelse(onoffsum.melt2$onlabel==T&onoffsum.melt2$insample==T,"on label in sample",ifelse(onoffsum.melt2$onlabel==F&onoffsum.melt2$insample==T,"not on label in sample",ifelse(onoffsum.melt2$onlabel==T&onoffsum.melt2$insample==F,"on label not in sample","not on label not in sample")))


onoffsum.melt2$Group.1 <- factor(onoffsum.melt2$Group.1,levels=unique(as.character(onoffsum.melt2$Group.1)))
levels(onoffsum.melt2$Group.1)<-c("No.2","No.1","No.3","No.4","No.5","No.7","No.8","No.9","No.6","No.11","No.10","No.12","No.13","No.14","No.15","No.16","No.17")
onoffsum.melt2$Group.1 <- factor(onoffsum.melt2$Group.1, levels=c("No.1","No.2","No.3","No.4","No.5","No.6","No.7","No.8","No.9","No.10","No.11","No.12","No.13","No.14","No.15","No.16","No.17"))


onoffsum.melt2$variable <- factor(onoffsum.melt2$variable,levels = unique(as.character(onoffsum.melt2$variable)))

onoffsum.melt2$cs<-factor(onoffsum.melt2$onoff.final)
levels(onoffsum.melt2$cs)<-c("red","black","lightblue","purple")

levels(onoffsum.melt2$variable)

#hard remove bagri since its not in any samples
onoffsum.melt2<-onoffsum.melt2[onoffsum.melt2$variable!="B.agri",]

#add bacterial category
onoffsum.melt2$category<-onoffsum.melt2$variable
levels(onoffsum.melt2$category)<-c("Bacillus",rep("Bifidobacterium",4),rep("Lactobacillus",8),"Streptococcus","Streptococcus","Weizmannia","Brevibacillus","Lactobacillus","Other")

onoffsum.melt2<-onoffsum.melt2[order(onoffsum.melt2$category),]

require(ggh4x)
levels(onoffsum.melt2$variable)
tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20240129_aggregate_nolabs_gg2_BLAST_heatmap.tif", res=300, height=7, width=12, units="in", compression="zip")
ggplot(onoffsum.melt2, aes(y=forcats::fct_rev(variable), x=Group.1,label=variable)) + 
  geom_tile(aes(fill=onoff.final))+
  geom_text(aes(label = round(value, 1)))+
  scale_fill_manual(values=c("red","black","lightblue","purple"))+
  facet_grid(category~., scales = "free_y") +
  theme(panel.spacing.y=unit(0.2,"lines"))+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+
  xlab("")+
  ylab("")+
  labs(fill="Label vs Sample")+
  force_panelsizes(rows = c(1.7,4,9,3.2,3,1.5),cols = c(1, 1))
#+
#  facet_grid(. ~ factor(gsub("A|B|C|D","",id),levels=c("No 1","No 2","No 3","No 4","No 5","No 6","No 7","No 8","No 9","No 10","No 11","No 12","No 14","No 15","No 16","No 17")), drop=TRUE,scale="free",space="free_x")
dev.off()



tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20240129_aggregate_labs_gg2_BLAST_heatmap.tif", res=300, height=7, width=12, units="in", compression="zip")
ggplot(onoffsum.melt2, aes(y=forcats::fct_rev(variable), x=Group.1,label=variable)) + 
  geom_tile(aes(fill=onoff.final))+
  geom_text(aes(label = round(value, 1)))+
  scale_fill_manual(values=c("red","black","lightblue","purple"))+
  facet_grid(category~., scales = "free_y") +
  theme(panel.spacing.y=unit(0.2,"lines"))+
  #theme(axis.text.y=element_blank(), 
  #      axis.ticks.y=element_blank())+
  xlab("")+
  ylab("")+
  labs(fill="Label vs Sample")+
  force_panelsizes(rows = c(1.7,4,9,3.2,3,1.5),cols = c(1, 1))
#+
#  facet_grid(. ~ factor(gsub("A|B|C|D","",id),levels=c("No 1","No 2","No 3","No 4","No 5","No 6","No 7","No 8","No 9","No 10","No 11","No 12","No 14","No 15","No 16","No 17")), drop=TRUE,scale="free",space="free_x")
dev.off()



#######################################
#######################################
#aggregate run 1 and show next to run 2
temp_runs<-runs_m_agg
 
#refactor id's to order run 1 then run 2 when run 2 had data
merged_runs<-temp_runs[order(temp_runs$Brand,temp_runs$run),]
#refactor columns so that other is at the end
temp<-data.frame(merged_runs[,c(1:(which(colnames(merged_runs)=="Other")-1))],merged_runs[,c((which(colnames(merged_runs)=="Other")+1):ncol(merged_runs))],merged_runs[,c(which(colnames(merged_runs)=="Other"))])
names(temp)[ncol(temp)]<-"Other"

merged_runs<-temp

#read in each run's list of on/off label and merge
list.sort.r1<-read.csv("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_1/master_list_sorted.csv")
list.sort.r2<-read.csv("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_2/master_list_sorted.csv")
list.sort.r2$Number<-as.character(list.sort.r2$Number)

#refactor list to match merged_runs
list.temp<-bind_rows(list.sort.r1,list.sort.r2)
list.sort<-list.temp[order(as.numeric(gsub("[A-Z]","",list.temp$Number))),]
list.sort[is.na(list.sort)]<-FALSE

list.sort<-list.sort[!duplicated(list.sort$Brand),]

list.sort$Number
l<-list.sort[,-c(1:2)]
l.n<-colnames(l)

onoffsum<-data.frame(paste(merged_runs$run,"-",merged_runs$Brand),merged_runs[,-c(1:2)])
colnames(onoffsum)<-c("run - id",colnames(merged_runs[,-c(1:2)]))

require(ggplot2)
require(reshape2)
onoffsum.melt<-melt(onoffsum, id.vars = 1)
onoffsum.melt$id<-factor(onoffsum.melt$`run - id`,levels=unique(onoffsum.melt$`run - id`))
store.melt<-onoffsum.melt
id.list<-unique(onoffsum.melt$id)
# 
# #gotta make the column then fill it
# onoffsum.melt$label<-NA
# for(i in 1:length(id.list)){
#   num<-which(onoffsum.melt$`run - id`==unique(onoffsum.melt$`run - id`)[i])
#   onoffsum.melt$label[num]<-ifelse(grepl(gsub("\\.","|",paste(gsub(".*\\.\\.","",l.n[l[i,]==T]),collapse = "|")),onoffsum.melt$variable[num],ignore.case = T)
#                                    ,"on label","off label")
# }
# 
# #create aggregate table base on label
# temp<-aggregate(round(value,2) ~ id + label,data=onoffsum.melt,FUN=sum)
# test<-recast(temp, id ~ label)

#write.csv(test, "Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/Run_2/plots/20231129_onoff_label_overview.csv",row.names = F)


#tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20231213_run1_run2_gg2_BLAST_onoff_label_overview.tif", res=300, height=7, width=20, units="in", compression="zip")
# ggplot(onoffsum.melt, aes(fill=label, y=value, x=id)) + 
#   geom_bar(position="stack",  stat="identity")+
#   labs(x ="", y = "% Abundance",fill="Product \nComposition")+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   
#   facet_grid(. ~ factor(gsub("A|B|C|D","",id),
#                         levels=c("No 1","No 2","No 3","No 4","No 5",
#                                  "No 6","No 7","No 8","No 9","No 10",
#                                  "No 11","No 12", "No 13","No 14",
#                                  "No 15","No 16","No 17")), 
#              drop=TRUE,scale="free",space="free_x")+
#   
#   stat_summary(aes(x = id, y = value),
#                fun.y = sum,
#                geom = "col",
#                colour = ifelse(grepl("A|B|C|D",unique(onoffsum.melt$id)),0,1),
#                fill = ifelse(grepl("A|B|C|D",unique(onoffsum.melt$id)),NA,NA))
# cant get the boxes to outline or insert a shape parameter because I'm not calling the data in an aes statement correctly
# geom_segment(aes(xend = id, yend = ifelse(grepl("A|B|C|D",id),0,1),linetype=ifelse(grepl("A|B|C|D",id),"Run 1","Run 2")),linewidth=NA, colour = ifelse(grepl("A|B|C|D",onoffsum.melt$id),0,1))+
# 
# scale_linetype_manual(values = c(1,2), name = "Sequencing \nRun")

# dev.off()
# 
# #order sample number, variables and create colors
# require(pals)
# onoffsum.melt$id <- factor(onoffsum.melt$id,levels=unique(as.character(onoffsum.melt$id)))
# onoffsum.melt$variable <- factor(onoffsum.melt$variable,levels = unique(as.character(onoffsum.melt$variable)))
# onoffsum.melt$col<-onoffsum.melt$variable
# on.col<-unname(alphabet())[c(1:length(unique(onoffsum.melt$col)))]
# 
# tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20231213_run1_run2_gg2_BLAST_onoff_label_details.tif", res=300, height=7, width=15, units="in", compression="zip")
# ggplot(onoffsum.melt, aes(fill=variable, y=value, x=label)) + 
#   geom_bar(position="stack", stat="identity")+
#   scale_fill_manual(values = on.col)+
#   facet_wrap(~id)+
#   labs(x ="", y = "% Abundance",fill="Product \nComposition")
# dev.off()


#heatmap showing present/absent/should be present
pres.table<-ifelse(onoffsum[,-1]>0.5,T,F)
names.table<-gsub("\\.\\.\\..*","",colnames(pres.table))
names.table<-gsub(".*\\.","",names.table)

#remove samples with no replicates
onoffsum2<-onoffsum[-c(11,16,17,22,25,26),]
onoffsum.melt<-melt(onoffsum2, id.vars = 1)
onoffsum.melt$id<-factor(onoffsum.melt$`run - id`,levels=unique(onoffsum.melt$`run - id`))
store.melt<-onoffsum.melt

onoffsum.melt2<-store.melt

id_unique<-unique(gsub("run1 - |run2 - ","",onoffsum.melt2$`run - id`))

onoffsum.melt2$onlabel<-NA
for(i in id_unique){
  num<-which(gsub("run1 - |run2 - ","",onoffsum.melt2$id) == i) #grabs the location of each sample in the melt table
  onlabel<-colnames(list.sort)[list.sort[list.sort$Brand==i,]==T]
  temp<-gsub("\\.","|",paste(gsub(".*\\.\\.","",onlabel),collapse = "|"))
  onoffsum.melt2$onlabel[num]<-grepl(temp,onoffsum.melt2$variable[num],ignore.case = T)
}  

onoffsum.melt2$insample<-ifelse(onoffsum.melt2$value>0.5,T,F)

onoffsum.melt2$onoff.final<-ifelse(onoffsum.melt2$onlabel==T&onoffsum.melt2$insample==T,"on label in sample",ifelse(onoffsum.melt2$onlabel==F&onoffsum.melt2$insample==T,"not on label in sample",ifelse(onoffsum.melt2$onlabel==T&onoffsum.melt2$insample==F,"on label not in sample","not on label not in sample")))


onoffsum.melt2$id <- factor(onoffsum.melt2$id,levels=unique(as.character(onoffsum.melt2$id)))
onoffsum.melt2$variable <- factor(onoffsum.melt2$variable,levels = unique(as.character(onoffsum.melt2$variable)))



onoffsum.melt2$`run - id` <- factor(onoffsum.melt2$`run - id`,levels=unique(as.character(onoffsum.melt2$`run - id`)))

levels(onoffsum.melt2$`run - id`)<-c("No.2 - run1",	"No.2 - run2",
                                  "No.1 - run1",	"No.1 - run2",
                                  "No.3 - run1",	"No.3 - run2",
                                  "No.4 - run1",	"No.4 - run2",
                                  "No.5 - run1",	"No.5 - run2",

                                  "No.8 - run1",	"No.8 - run2",
                                  "No.9 - run1",	"No.9 - run2",

                                  
                                  "No.10 - run1",	"No.10 - run2",
                                  "No.12 - run1",	"No.12 - run2",

                                  "No.14 - run1",	"No.14 - run2",

                                  "No.17 - run1",	"No.17 - run2")
onoffsum.melt2$`run - id` <- factor(onoffsum.melt2$`run - id`, 
                                 levels=c("No.1 - run1",	"No.1 - run2",
                                          "No.2 - run1",	"No.2 - run2",
                                          "No.3 - run1",	"No.3 - run2",
                                          "No.4 - run1",	"No.4 - run2",
                                          "No.5 - run1",	"No.5 - run2",

                                          "No.8 - run1",	"No.8 - run2",
                                          "No.9 - run1",	"No.9 - run2",
                                          "No.10 - run1",	"No.10 - run2",

                                          "No.12 - run1",	"No.12 - run2",

                                          "No.14 - run1",	"No.14 - run2",

                                          "No.17 - run1",	"No.17 - run2"))


onoffsum.melt2$variable <- factor(onoffsum.melt2$variable,levels = unique(as.character(onoffsum.melt2$variable)))

onoffsum.melt2$cs<-factor(onoffsum.melt2$onoff.final)
levels(onoffsum.melt2$cs)<-c("red","black","lightblue","purple")

levels(onoffsum.melt2$variable)

#hard remove bagri since its not in any samples
onoffsum.melt2<-onoffsum.melt2[onoffsum.melt2$variable!="B.agri",]

#add bacterial category
onoffsum.melt2$category<-onoffsum.melt2$variable
levels(onoffsum.melt2$category)<-c("Bacillus",rep("Bifidobacterium",4),rep("Lactobacillus",8),"Streptococcus","Streptococcus","Weizmannia","Brevibacillus","Lactobacillus","Other")

onoffsum.melt2<-onoffsum.melt2[order(onoffsum.melt2$category),]

require(ggh4x)
levels(onoffsum.melt2$variable)
tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20240129_aggregate_nolabs_r1_sep_r2_gg2_BLAST_heatmap.tif", res=300, height=7, width=12, units="in", compression="zip")
ggplot(onoffsum.melt2, aes(y=forcats::fct_rev(variable), x=`run - id`,label=variable)) + 
  geom_tile(aes(fill=onoff.final))+
  geom_text(aes(label = round(value, 1)))+
  scale_fill_manual(values=c("red","black","lightblue","purple"))+
  facet_grid(category~., scales = "free_y") +
  theme(panel.spacing.y=unit(0.2,"lines"))+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+
  xlab("")+
  ylab("")+
  labs(fill="Label vs Sample")+
  force_panelsizes(rows = c(1.7,4,9,3.2,3,1.5),cols = c(1, 1))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5), color = "white") 
 #+
#  facet_grid(. ~ factor(gsub("A|B|C|D","",id),levels=c("No 1","No 2","No 3","No 4","No 5","No 6","No 7","No 8","No 9","No 10","No 11","No 12","No 14","No 15","No 16","No 17")), drop=TRUE,scale="free",space="free_x")
dev.off()



tiff("Z:/Research/Friedman-Jonscher Lab/Sugino/Chabaan_Lab_Probiotics_Analysis/20240129_aggregate_labs_r1_sep_r2_gg2_BLAST_heatmap.tif", res=300, height=7, width=12, units="in", compression="zip")
ggplot(onoffsum.melt2, aes(y=forcats::fct_rev(variable), x=`run - id`,label=variable)) + 
  geom_tile(aes(fill=onoff.final))+
  geom_text(aes(label = round(value, 1)))+
  scale_fill_manual(values=c("red","black","lightblue","purple"))+
  facet_grid(category~., scales = "free_y") +
  theme(panel.spacing.y=unit(0.2,"lines"))+
  #theme(axis.text.y=element_blank(), 
  #      axis.ticks.y=element_blank())+
  xlab("")+
  ylab("")+
  labs(fill="Label vs Sample")+
  force_panelsizes(rows = c(1.7,4,9,3.2,3,1.5),cols = c(1, 1))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept = c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5), color = "white") 
#+
#  facet_grid(. ~ factor(gsub("A|B|C|D","",id),levels=c("No 1","No 2","No 3","No 4","No 5","No 6","No 7","No 8","No 9","No 10","No 11","No 12","No 14","No 15","No 16","No 17")), drop=TRUE,scale="free",space="free_x")
dev.off()






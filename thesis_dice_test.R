library(dplyr)
omgetwd()
setwd("C:/Users/rathi/Downloads/DIce Test")


#install.packages('mclustcomp')

BiocManager::install("OmicsMarkeR")



#Reading Cluster X files
file_names_clusterx=list.files('./clutsetrx/')
cluster_values_clusterx=as.numeric(gsub("[^\\d]+", "", file_names_clusterx, perl=TRUE))
df_clusterx=data.frame()
i=1
for(i in 1:length(file_names_clusterx)){
df_clusterx_1=read.csv(paste('./clutsetrx/',file_names_clusterx[i],sep=''))
df_clusterx_1['clusterx_cluster_number']=cluster_values_clusterx[i]
#print(file_names_clusterx[i])
#print(cluster_values_clusterx[i])
#print(i)
#print(ncol(df_clusterx_1))
df_clusterx=dplyr::bind_rows(df_clusterx,df_clusterx_1)
}


#Reading Flowsom X files
file_names_flowsom=list.files('./flowsom/')

cluster_values_flowsom=as.numeric(gsub("[^\\d]+", "", substr(file_names_flowsom,10,100), perl=TRUE))
df_flowsom=data.frame()

for(i in 1:length(file_names_flowsom)){
df_flowsom_1=read.csv(paste('./flowsom/',file_names_flowsom[i],sep=''))
df_flowsom_1['flowsom_cluster_number']=cluster_values_flowsom[i]
#print(i)
#print(file_names_flowsom[i])
#print(cluster_values_flowsom[i])
#print(ncol(df_clusterx_1))
df_flowsom=dplyr::bind_rows(df_flowsom,df_flowsom_1)
}




#Reading Rpheno X files
file_names_Rpheno=list.files('./Rpheno/')
cluster_values_Rpheno=as.numeric(gsub("[^\\d]+", "", file_names_Rpheno, perl=TRUE))
df_Rpheno=data.frame()

for(i in 1:length(file_names_Rpheno)){
df_Rpheno_1=read.csv(paste('./Rpheno/',file_names_Rpheno[i],sep=''))
df_Rpheno_1['rpheno_cluster_number']=cluster_values_Rpheno[i]
#print(i)
#print(file_names_flowsom[i])
#print(cluster_values_flowsom[i])
#print(ncol(df_clusterx_1))
df_Rpheno=dplyr::bind_rows(df_Rpheno,df_Rpheno_1)
}



colnames(df_clusterx)[1]<-'Cell_Sample_Name'
colnames(df_flowsom)[1]<-'Cell_Sample_Name'
colnames(df_Rpheno)[1]<-'Cell_Sample_Name'

####saving the files in rda
save(df_clusterx,df_flowsom,df_Rpheno,file = 'three_cluster_dataset.rda')


#run from below
load('three_cluster_dataset.rda')

df_master=df_clusterx%>%
dplyr::inner_join(df_flowsom%>%
select(Cell_Sample_Name,flowsom_cluster_number),
on='Cell_Sample_Name',
how='inner')%>%
dplyr::inner_join(df_Rpheno%>%
select(Cell_Sample_Name,rpheno_cluster_number),
on='Cell_Sample_Name',
how='inner')



df_master<-df_master%>%mutate(clusterx_signif=ifelse(clusterx_cluster_number==12,1,0),
                   flowsom_signif=ifelse(flowsom_cluster_number==17,1,0),
                   rpheno_signif=ifelse(rpheno_cluster_number==10,1,0))



#RUnning dice test for  
#clusterx vs flowson
OmicsMarkeR::sorensen(as.character(df_master$clusterx_signif), 
                      as.character(df_master$flowsom_signif))

#1

#rpheno vs clusterx
OmicsMarkeR::sorensen(as.character(df_master$clusterx_signif), 
                      as.character(df_master$rpheno_signif))
# 1
#rpheno vs flowson
OmicsMarkeR::sorensen(as.character(df_master$flowsom_signif), 
                      as.character(df_master$rpheno_signif))
# 1

#Grouping marker analysis
#clusterx
markerlist=df_clusterx%>%select(2:50)%>%colnames()
df_clusterx=df_clusterx%>%mutate(clusterx_signif=ifelse(clusterx_cluster_number==12,1,0))
#creating sign and other df
df_clusterx_sign_cluster=df_clusterx%>%filter(clusterx_signif==1)
df_clusterx_other_cluster=df_clusterx%>%filter(clusterx_signif==0)

clusterx_p_value_df=data.frame()
for(i in 1:length(markerlist)){
marker=markerlist[i]
sign_vec=df_clusterx_sign_cluster[,marker]
unsign_vec=df_clusterx_other_cluster[,marker]
t_test<-t.test(sign_vec,unsign_vec)
difference_between_markers=t_test$statistic
p_value_between_markers=t_test$p.value
clusterx_df_temp=data.frame(marker=marker,
                            difference_between_markers=difference_between_markers,
                            p_value_between_markers=p_value_between_markers)
clusterx_p_value_df<-clusterx_p_value_df%>%dplyr::bind_rows(clusterx_df_temp)
}

#flowsom
markerlist=df_flowsom%>%select(2:50)%>%colnames()
df_flowsom=df_flowsom%>%mutate(flowsom_signif=ifelse(flowsom_cluster_number==17,1,0))

#creating sign and other df
df_flowsom_sign_cluster=df_flowsom%>%filter(flowsom_signif==1)
df_flowsom_other_cluster=df_flowsom%>%filter(flowsom_signif==0)

flowsom_p_value_df=data.frame()
for(i in 1:length(markerlist)){
  marker=markerlist[i]
  sign_vec=df_flowsom_sign_cluster[,marker]
  unsign_vec=df_flowsom_other_cluster[,marker]
  t_test<-t.test(sign_vec,unsign_vec)
  difference_between_markers=t_test$statistic
  p_value_between_markers=t_test$p.value
  df_temp=data.frame(marker=marker,
                              difference_between_markers=difference_between_markers,
                              p_value_between_markers=p_value_between_markers)
  flowsom_p_value_df<-flowsom_p_value_df%>%dplyr::bind_rows(df_temp)
}

# rpheno

markerlist=df_Rpheno%>%select(2:50)%>%colnames()
df_Rpheno=df_Rpheno%>%mutate(rpheno_signif=ifelse(rpheno_cluster_number==10,1,0))

#creating sign and other df
df_rpheno_sign_cluster=df_Rpheno%>%filter(rpheno_signif==1)
df_rpheno_other_cluster=df_Rpheno%>%filter(rpheno_signif==0)

Rpheno_p_value_df=data.frame()
for(i in 1:length(markerlist)){
  marker=markerlist[i]
  sign_vec=df_rpheno_sign_cluster[,marker]
  unsign_vec=df_rpheno_other_cluster[,marker]
  t_test<-t.test(sign_vec,unsign_vec)
  difference_between_markers=t_test$statistic
  p_value_between_markers=t_test$p.value
  rpheno_df_temp=data.frame(marker=marker,
                     difference_between_markers=difference_between_markers,
                     p_value_between_markers=p_value_between_markers)
  Rpheno_p_value_df<-flowsom_p_value_df%>%dplyr::bind_rows(rpheno_df_temp)
}






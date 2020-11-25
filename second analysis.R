library(dplyr)
library(stringr)
library(lsr)
setwd("C:/Users/rathi/Downloads/DIce Test")
load('three_cluster_dataset.rda')

df_clusterx_signf<-df_clusterx%>%filter(clusterx_cluster_number==12)
markers<-df_clusterx_signf%>%select(2:50)%>%colnames()

healthy_df_clusterx<-df_clusterx_signf%>%filter(str_detect(Cell_Sample_Name,'Healthy') )
unhealthy_df_clusetrx<-df_clusterx_signf%>%filter(!str_detect(Cell_Sample_Name,'Healthy') )

df_sign_cluster_x_cohen=data.frame()
for(i in 1:length(markers)){
  marker=markers[i]
  v1=healthy_df_clusterx[,marker]
  v2=unhealthy_df_clusetrx[,marker]
  cohen_d=lsr::cohensD(v1,v2)
  p_value=t.test(v1,v2)
  p_value=p_value$p.value
  df_temp=data.frame(marker=marker,cohen_d=cohen_d,p_value=p_value)
  df_sign_cluster_x_cohen<-df_sign_cluster_x_cohen%>%bind_rows(df_temp)
}
df_sign_cluster_x_cohen$sign_p_value=ifelse(df_sign_cluster_x_cohen$p_value<=0.05,1,0)



df_flowsom_sign<-df_flowsom%>%filter(flowsom_cluster_number==17)
markers<-df_flowsom_sign%>%select(2:50)%>%colnames()

healthy_df_flowsom<-df_flowsom_sign%>%filter(str_detect(Cell_Sample_Name,'Healthy') )
unhealthy_df_flowsom<-df_flowsom_sign%>%filter(!str_detect(Cell_Sample_Name,'Healthy') )

df_sign_flowsom_cohen=data.frame()
for(i in 1:length(markers)){
  marker=markers[i]
  v1=healthy_df_flowsom[,marker]
  v2=unhealthy_df_flowsom[,marker]
  cohen_d=lsr::cohensD(v1,v2)
  p_value=t.test(v1,v2)
  p_value=p_value$p.value
  df_temp=data.frame(marker=marker,cohen_d=cohen_d,p_value=p_value)
  df_sign_flowsom_cohen<-df_sign_flowsom_cohen%>%bind_rows(df_temp)
}
df_sign_flowsom_cohen$sign_p_value=ifelse(df_sign_flowsom_cohen$p_value<=0.05,1,0)

#Rpheno
df_Rpheno_sign<-df_Rpheno%>%filter(rpheno_cluster_number==10)
markers<-df_Rpheno_sign%>%select(2:50)%>%colnames()

healthy_df_Rpheno<-df_Rpheno_sign%>%filter(str_detect(Cell_Sample_Name,'Healthy') )
unhealthy_df_Rpheno<-df_Rpheno_sign%>%filter(!str_detect(Cell_Sample_Name,'Healthy') )

df_sign_Rpheno_cohen=data.frame()
for(i in 1:length(markers)){
  marker=markers[i]
  v1=healthy_df_Rpheno[,marker]
  v2=unhealthy_df_Rpheno[,marker]
  cohen_d=lsr::cohensD(v1,v2)
  p_value=t.test(v1,v2)
  p_value=p_value$p.value
  df_temp=data.frame(marker=marker,cohen_d=cohen_d,p_value=p_value)
  df_sign_Rpheno_cohen<-df_sign_Rpheno_cohen%>%bind_rows(df_temp)
}
df_sign_Rpheno_cohen$sign_p_value=ifelse(df_sign_Rpheno_cohen$p_value<=0.05,1,0)

write.csv(df_sign_cluster_x_cohen,"/Users/rathi/Downloads/DIce Test/clusterxcohen.csv", row.names = FALSE)
write.csv(df_sign_flowsom_cohen,"/Users/rathi/Downloads/DIce Test/flowsomcohen.csv", row.names = FALSE)
write.csv(df_sign_Rpheno_cohen,"/Users/rathi/Downloads/DIce Test/rphenoraphcohen.csv", row.names = FALSE)

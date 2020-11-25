#T Test
library(dplyr)
setwd( "Users/rathi/OneDrive/Desktop/flowsom run")
df=read.csv('flowsom#1_FlowSOM_cluster_cell_percentage.csv')
df


#assuming p-value cutoff to be 0.05

df_p_value<-data.frame()
for(i in 1:nrow(df)){
  healthy_values_for_the_cluster<-df[i,]%>%select(contains('Healthy'))
  tb_values_for_the_cluster<-df[i,]%>%select(contains('TB'))
  healthy_values_for_the_cluster<-as.numeric(healthy_values_for_the_cluster)
  tb_values_for_the_cluster<-as.numeric(tb_values_for_the_cluster)
  t_test<-t.test(healthy_values_for_the_cluster,tb_values_for_the_cluster)
  df_p_value_temp<-data.frame(cluster=i,
                              p_value=t_test$p.value)
  df_p_value<-dplyr::bind_rows(df_p_value,df_p_value_temp)
}

p_value_threshold=0.05
df_p_value%>%mutate(p_value_flag=ifelse(p_value<p_value_threshold,'significant difference between healthy and tb',
                                        'not significant'))

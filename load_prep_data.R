load_data = function(){
  ############################################################
  ############################################################
  
  ##prep exp1 data
  data1 = read.csv("dataStudy1.csv")
  data1 = data1[,c(2,8,9,10,11,12,13)]
  data1$version = 1
  ##prep exp2 data
  data2 = read.csv("dataStudy2.csv")
  data2 = data2[,c(2,8,9,10,11,12,13)]
  data2$version = 2
  ##prep exp3 data
  data3 = read.csv("dataStudy3.csv")
  data3 = data3[,c(2,8,9,10,11,12,13)]
  data3$version = 3
  #combine all dataframes
  data = rbind(data1,data2,data3)
  #get rid of bad subjects
  bad_subs = c("3100620","2893759","3658195")
  data = data[!data$Participant.Private.ID%in%bad_subs,]
  
  # make vector for all subj IDs
  S = unique(data$Participant.Private.ID)
  
  versions = c()#version per sub
  effort_selfreport = happy_selfreport = c()
  data$LL = ifelse(data$effortChoice=="high",1,0)#larger-later reward yes/no
  data$T1 = 0 #effort for SS
  data$T2 = data$effortL # effort for LL
  data$X1 = 1 #reward for SS
  data$X2 = data$rewardL #reward for LL
  data$s2_u = 0 # reward variance/uncertainty
  
  data$risk = ifelse(data$gamblingChoice=="risky",1,0)
  data$version[data$version==2]=1
  return(data)
}
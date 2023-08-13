#R codes for data analyses and plotting
#including ARG abundance based sourcetracking, and .
#NOTE: the R script below has been simplified for illurstration

###################
#data preparation #
###################
#data preparation #
###################

source('SourceTracker.r')
> database_str<-read.delim(file="database_str.txt",sep="\t",head=F)
> colnames(database_str)<-c("seqID","subtype","type")
> sample_info<-read.delim(file="sample_info1.0.txt",sep='\t',header=T)
> sample_info$eco_type<-factor(sample_info$eco_type,levels = c("HF","AF","WA","NT","UN"))
> sample_info<-sample_info[order(sample_info$eco_type),]
> sample_info$sampleID<-factor(sample_info$sampleID,levels=sample_info$sampleID)
> ARG_abund<-read.delim(file="ARG_abund1.0.txt",sep='\t',header=T)
> ARG_abundmatrix<-as.matrix(ARG_abund)
> class(ARG_abundmatrix)<-"numeric"
> ARG_abundmatrix<-ARG_abundmatrix[sample_info$sampleID,]
> maxname<-function(data_frame){
  +   envs<-c("HF","AF","WA","NT","UN")
  +   max_name<-envs[which.max(data_frame[2:6])]
  +   return(max_name)}
> maxratio<-function(data_frame){
  +   max_ratio<-data_frame[2:6][which.max(data_frame[2:6])]
  +   return(max_ratio)}

#leave-one-out strategy
> pred_list<-list()
> proptab_list<-list()
> for (i in 1:nrow(ARG_abundmatrix))
  + {
    +   st<-sourcetracker(ARG_abundmatrix[-i,], sample_info$eco_type[-i])
    +   st_predic1<- predict(st,ARG_abundmatrix[i,], alpha1=0.001, alpha2=0.001)
    +   st_predic2 <- as.data.frame(st_predic1[[2]])
    +   st_predic2$sampleID<-sample_info$sampleID[i]
    +   st_predic2<-merge(st_predic2,sample_info,by="sampleID")
    +   pred_list[[i]]<-st_predic1
    +   proptab_list[[i]]<-st_predic2
    + }
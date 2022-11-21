## code to prepare `chr_dat_l` dataset goes here
library(Matrix)
library(data.tree)
library(tidyverse)
library(parallel)
library(caret)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
########################################################################################################
# Data dowloaded from the Juicertool repository
# https://encode-public.s3.amazonaws.com/2019/02/08/fc1d9d5d-8fa0-4e29-9080-3da674d9490d/ENCFF543USQ.hic
# https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic
########################################################################################################
## Util. fn. to load R-objects
get_tbl_in_fn<-function(tmp_file){
  out_tbl<-base::get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

########################################################################################################
dat_folder<-"~/Documents/multires_bhicect/data/HMEC/"
########################################################################################################

chromo<-"chr22"
hub_res_set<-grep("b$",list.files(dat_folder),value=T)
hub_res_set<-names(sort(res_num[hub_res_set],decreasing = T))
#-----------------------
#Extract the interactions of every constitutive child-cluster
chr_dat_l<-lapply(hub_res_set,function(x)read_delim(file = paste0(dat_folder,x,'/',chromo,'.txt'),delim = '\t',col_names = F,col_types = list("i","i","d")))
names(chr_dat_l)<-hub_res_set
chr_dat_l<-lapply(chr_dat_l,function(x){
  x%>%
    filter(!(is.nan(X3)))
})
chr_dat_l<-lapply(chr_dat_l,function(x){
  preprocessParams <- BoxCoxTrans(x$X3,na.rm = T)
  x <- data.frame(x,weight=predict(preprocessParams, x$X3))
  x$weight<-x$weight+(1-min(x$weight,na.rm = T))
  return(x)
})

usethis::use_data(chr_dat_l, overwrite = TRUE,compress = "xz")

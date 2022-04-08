library(dplyr)
library(peaky)
	
        base <-commandArgs(TRUE)[1]
	interactions_file <-commandArgs(TRUE)[2]
	bins_dir <-commandArgs(TRUE)[3]
	fragments_file<-commandArgs(TRUE)[4] 

bin_interactions_fs(interactions_file, fragments_file, output_dir=bins_dir,bins=10)

fits_dir = paste0(base,"/fits")

for(bin_index in 1:10){
  model_bin_fs(bins_dir,bin_index,output_dir=fits_dir,subsample_size=1000)
}

baits_dir = paste0(base,"/baits")
split_baits_fs(bins_dir,residuals_dir = fits_dir, indices=1:10, output_dir = baits_dir)

allpairstat <- list.files(path=baits_dir, pattern=paste0("^[[:print:]]*\\.rds$"), full.names=TRUE)

pky <- do.call(dplyr::bind_rows,lapply(allpairstat,readRDS)) 

write.table(paste(pky$b.chr,pky$b.mid,pky$dist,pky$p.chr,pky$p.mid-as.integer(pky$p.length/2),pky$p.mid+as.integer(pky$p.length/2),-log(pky$fdr.res)),paste0(interactions_file,"_fdr.wig"),quote=FALSE,col.names=FALSE,row.names=FALSE)
rjmcmc_dir = paste0(base,"/rjmcmc")
omega_power = -3.8
baitlist = paste0(baits_dir,"/baitlist.txt")
baitlist2 = paste0(baits_dir,"/baitlist2.txt")
cmd <- paste0(paste0(paste0("cat ", baitlist), "|gawk \'{print \" ls -lt \" $1}\'|sh|gawk \'{if($5>10000){print $NF}}\' > "), baitlist2) 
system(cmd)

list<-read.table(baitlist2)

for(i in 1:dim(list)[1]){
print(i)
peaky_fs(baitlist2,i,output_dir=rjmcmc_dir,omega_power=omega_power)
}

rjmcmc_list(rjmcmc_dir)

rjmcmclist = paste0(rjmcmc_dir,"/rjmcmclist.txt")
baits_rjmcmc_dir = paste0(base,"/baits_rjmcmc")
for(i in 1:dim(list)[1]){
print(i)
interpret_peaky_fs(rjmcmclist,i,baits_dir,baits_rjmcmc_dir,omega_power)
}

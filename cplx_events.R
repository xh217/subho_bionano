# Extract intricate blocks corresponding to each complex event
cplx_events <-function(df,threshold) {
  library(stringr)
  ## initiating the genomic block list
  gblk=df[,1]
  gblk.list=list()
  for(i in 1:nrow(df)) {
    M.grp=strsplit(df[i,2],",")
    gblk.list[df[i,1]]<-M.grp
  }
  
  ## initiating the contig list
  contigs=unique(unlist(gblk.list, use.names=FALSE))
  t=as.matrix(stack(gblk.list))[, c(1,2)]
  ctg.list=split(t[,2],f=t[,1])
  
  ## identify the seeding complex regions with multiple contigs to initiate
  seed.gblk=gblk[which(lengths(gblk.list)>=threshold)]
  
  ## complex sv identification
  complex.sv.t=list()
  cblk=c()
  cctg=c()
  
  for(i in 1:length(seed.gblk)) {
    cplx.blk=seed.gblk[i]
    cplx.ctg=unique(unlist(gblk.list[seed.gblk[i]]))
    cplx.blk.pre=c()
    cplx.ctg.pre=c()
    loop=1
    
    ## iteratively collect all interacting complex blocks and contigs therein
    while(loop) {
      cplx.blk.pre=cplx.blk
      cplx.ctg.pre=cplx.ctg
      
      cplx.blk=unique(c(cplx.blk, unlist(ctg.list[unlist(cplx.ctg)])))
      cplx.ctg=unique(c(cplx.ctg,unlist(gblk.list[unlist(cplx.blk)])))
      
      if(identical(str_sort(cplx.blk.pre,numeric = TRUE),str_sort(cplx.blk,numeric = TRUE)) & identical(str_sort(cplx.ctg.pre,numeric = TRUE),str_sort(cplx.ctg,numeric = TRUE))) { 
      loop=0
      }
    }
    
    ## Redundant complex SV calling
      complex.sv.t=append(complex.sv.t,list(list("blocks"=sort(cplx.blk),"contigs"=sort(cplx.ctg))))
  }
  
  ## Non-redundant complex SV calling
  ctgs.t=unique(sapply(complex.sv.t, function(x)x["contigs"]))
  blks.t=list()
  for(j in 1:length(ctgs.t)) {
    blks.t[[j]]=unique(unlist(ctg.list[ctgs.t[[j]]], use.names=FALSE))
  }
  complex.sv<-list("blocks"=blks.t,"contigs"=ctgs.t)
  return(complex.sv)
}

# 1. Four-group TDI protein comparisons
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))
use('omics')
if(sys.nframe()==0L) {
  file<-file.path(out,'Supp.rds')
  r<-if(file.exists(file))readRDS(file) else list()
  r$omics<-list(Protein=supp_omics())
  r$pattern<-supp_pattern(r)
  a<-r$omics$Protein$contrasts
  features<-unique(a[bonferroni<.05][order(p)][,head(.SD,6),by=comparison]$feature)
  small<-r$pattern$map[,.(n=.N),by=pattern][n<=12,pattern]
  features<-unique(c(features,r$pattern$map[pattern%in%small,feature]))
  mat<-dcast(a[feature%in%features],feature~comparison,value.var='beta')
  M<-as.matrix(mat[,-1]);rownames(M)<-mat$feature
  r$display<-list(proteins=rownames(M)[hclust(dist(M))$order])
  saveRDS(r,file,compress='gzip')
}

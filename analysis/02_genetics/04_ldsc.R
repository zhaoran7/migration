# 1. LDSC preparation and collection
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))


# 2. Prepare signed summary statistics
ldsc_prepare <- function(file,map,output) {
  x<-fread(file);map<-fread(map)
  stopifnot(all(c('CHR','POS','SNP')%in%names(map)),!anyDuplicated(map[,.(CHR,POS)]))
  x[map,on=.(CHR,POS),SNP:=i.SNP]
  x<-x[!is.na(SNP)&nzchar(SNP)&is.finite(BETA)&is.finite(SE)&SE>0&
    is.finite(P)&P>=0&P<=1&is.finite(N)&N>0]
  x[,P:=pmax(P,1e-300)]
  fwrite(x[,.(SNP,EA,NEA,BETA,P,N)],output,sep='\t')
}


# 3. Collect heritability and genetic correlations
ldsc_collect <- function(folder) {
  alias <- function(x) {
    x<-sub('\\.sumstats(\\.gz)?$','',basename(x))
    x[x=='migration_distance_pc40']<-'migration_distance'
    x[x=='delta_tdi_ukb_cond']<-'migration_direction'
    x
  }
  values <- function(lines,pattern) {
    s<-grep(pattern,lines,value=TRUE)
    if(!length(s))return(c(NA_real_,NA_real_))
    s<-sub(pattern,'',tail(s,1))
    v<-regmatches(s,gregexpr('[-+]?[0-9]+(?:\\.[0-9]*)?(?:[eE][-+]?[0-9]+)?',s,perl=TRUE))[[1]]
    c(as.numeric(v),NA_real_)[1:2]
  }
  h2 <- rbindlist(lapply(list.files(file.path(folder,'h2'),'\\.log$',full.names=TRUE),function(file) {
    s<-readLines(file,warn=FALSE)
    h<-values(s,'^Total Observed scale h2:');i<-values(s,'^Intercept:');r<-values(s,'^Ratio:')
    name<-sub('\\.log$','',basename(file))
    data.table(trait=alias(name),h2=h[1],h2_se=h[2],lambda_gc=values(s,'^Lambda GC:')[1],
      mean_chi2=values(s,'^Mean Chi\\^2:')[1],intercept=i[1],intercept_se=i[2],ratio=r[1],ratio_se=r[2],GWAS=name)
  }))
  rg <- rbindlist(lapply(list.files(file.path(folder,'rg'),'\\.log$',full.names=TRUE),function(file) {
    s<-readLines(file,warn=FALSE);i<-grep('^\\s*p1\\s+p2\\s+rg\\s+se\\s+z\\s+p',s)
    if(!length(i))stop('Missing LDSC correlation summary: ',file)
    d<-fread(text=paste(s[i[1]+0:1],collapse='\n'))
    d[,`:=`(GWAS1=sub('\\.sumstats(\\.gz)?$','',basename(p1)),
      GWAS2=sub('\\.sumstats(\\.gz)?$','',basename(p2)),p1=alias(p1),p2=alias(p2))]
    d
  }),fill=TRUE)
  stopifnot(nrow(h2)>0,nrow(rg)>0,all(is.finite(h2$h2)),all(is.finite(rg$rg)))
  fwrite(h2,file.path(folder,'all.h2.res'),sep='\t')
  fwrite(rg,file.path(folder,'all.rg.res'),sep='\t')
  published<-Sys.getenv('PUBLISHED_MD','published_md')
  p<-rg[(GWAS1=='migration_distance_pc40'&GWAS2==published)|
    (GWAS2=='migration_distance_pc40'&GWAS1==published)]
  if(nrow(p))fwrite(p,file.path(folder,'md_published.res'),sep='\t')
}

if(sys.nframe()==0L) {
  a<-commandArgs(TRUE)
  if(a[1]=='prepare')ldsc_prepare(a[2],a[3],a[4]) else ldsc_collect(a[1])
}

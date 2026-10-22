# 1. Merge independent aggregate results
code <- Sys.getenv('MOVE_CODE')
if (!nzchar(code)) code <- dirname(dirname(normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(), value=TRUE)[1]))))
Sys.setenv(MOVE_CODE=code)
source(file.path(code, '00_common/01_base.R'))
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE)
  stopifnot(length(a)>=3L)
  result<-list()
  for(file in a[-1]) {
    x<-readRDS(file);stopifnot(is.list(x),!is.null(names(x)))
    common<-intersect(names(result),names(x))
    if(any(!vapply(common,function(n)isTRUE(all.equal(result[[n]],x[[n]])),TRUE)))
      stop('Conflicting result components: ',file)
    result[names(x)]<-x
  }
  saveRDS(result,a[1],compress='gzip')
}

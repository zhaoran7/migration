# 1. Numerical source-data format
source_write <- function(x, file) {
  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, 'Contents')
  index <- list()
  json <- function(x) as.character(jsonlite::toJSON(x, auto_unbox=TRUE, null='null', digits=17))
  add <- function(x, parent=0L, name='data') {
    id <- length(index)+1L
    item <- list(ID=id, Parent=parent, Name=name, Kind='', Sheet='', Class='', Metadata='')
    index[[id]] <<- item
    if (is.function(x) || is.environment(x) || inherits(x,c('lm','coxph','coxme','merMod','ggplot','grob')))
      stop('Only numerical source data may be exported: ', name)
    item$Class <- paste(setdiff(class(x),'data.table'),collapse='|')
    if (is.null(x)) item$Kind <- 'null'
    else if (inherits(x,'sf')) {
      asset <- paste0(tools::file_path_sans_ext(basename(file)), '_', id, '.geojson')
      sf::st_write(x, file.path(dirname(file),asset), quiet=TRUE, delete_dsn=TRUE)
      item$Kind <- 'sf'; item$Sheet <- asset
    } else if (is.data.frame(x) || is.matrix(x) || is.atomic(x)) {
      kind <- if(is.data.frame(x)) 'table' else if(is.matrix(x)) 'matrix' else 'vector'
      values <- if(kind=='vector') data.frame(Value=x, check.names=FALSE) else as.data.frame(x)
      if (any(vapply(values,is.list,logical(1)))) {
        cols <- names(values)[vapply(values,is.list,logical(1))]
        values[cols] <- lapply(values[cols],function(v)vapply(v,json,''))
      } else cols <- character()
      types <- lapply(as.data.frame(values),function(v) list(type=typeof(v),class=class(v),levels=levels(v)))
      item$Kind <- kind
      item$Sheet <- sprintf('Data%03d',id)
      item$Metadata <- json(list(columns=names(values),types=types,list_columns=cols,
        rownames=if(kind=='matrix')rownames(x) else NULL,names=if(kind=='vector')names(x) else NULL,
        nrow=nrow(values),ncol=ncol(values)))
      openxlsx::addWorksheet(wb,item$Sheet)
      for(n in names(values)) {
        if(is.factor(values[[n]])) values[[n]] <- as.character(values[[n]])
        if(inherits(values[[n]],c('Date','IDate','POSIXct')))values[[n]] <- as.numeric(values[[n]])
        if(is.character(values[[n]]))values[[n]][values[[n]]==''] <- '<EMPTY>'
        if(is.numeric(values[[n]])&&any(is.infinite(values[[n]])|is.nan(values[[n]])))
          values[[n]] <- as.character(values[[n]])
      }
      labels <- c(n='N',event='Events',events='Events',beta='Beta',se='SE',HR='HR',
        lo='95% CI lower',hi='95% CI upper',lower='95% CI lower',upper='95% CI upper',
        p='p',q='q',FDR='q',p_bonf='p Bonferroni',exposure='Exposure',outcome='Outcome')
      old <- names(values);nice <- unname(labels[old]);missing <- is.na(nice)
      nice[missing] <- tools::toTitleCase(gsub('_',' ',old[missing]))
      pcol <- grepl('^[pq]($|_)',old)
      nice[pcol] <- paste0(substr(old[pcol],1,1),substring(nice[pcol],2))
      names(values) <- make.unique(nice,sep=' ')
      openxlsx::writeData(wb,item$Sheet,values,keepNA=TRUE,na.string='<NA>')
      if(any(pcol))openxlsx::addStyle(wb,item$Sheet,
        openxlsx::createStyle(textDecoration=c('bold','italic')),rows=1,cols=which(pcol),gridExpand=TRUE,stack=TRUE)
      openxlsx::freezePane(wb,item$Sheet,firstRow=TRUE)
      openxlsx::setColWidths(wb,item$Sheet,cols=seq_len(ncol(values)),widths=18)
    } else if(is.list(x)) {
      item$Kind <- 'list'
      item$Metadata <- json(list(names=names(x)))
      for(i in seq_along(x))add(x[[i]],id,if(is.null(names(x)))as.character(i) else names(x)[i])
    } else stop('Unsupported source-data type: ',name)
    index[[id]] <<- item
    invisible(id)
  }
  add(x)
  idx <- data.table::rbindlist(index)
  openxlsx::writeData(wb,'Contents',idx)
  openxlsx::setColWidths(wb,'Contents',1:7,c(7,8,28,12,22,20,65))
  header <- openxlsx::createStyle(fontName='Arial',fontSize=10,textDecoration='bold',
    fgFill='#F2F2F2',border='bottom',borderColour='#B9BEC2')
  for(s in names(wb)) {
    openxlsx::showGridLines(wb,s,FALSE)
    openxlsx::addStyle(wb,s,header,rows=1,cols=1:if(s=='Contents')7 else
      jsonlite::fromJSON(idx[Sheet==s,Metadata])$ncol,gridExpand=TRUE,stack=TRUE)
    if(s!='Contents') {
      original <- jsonlite::fromJSON(idx[Sheet==s,Metadata])$columns
      pcol <- grepl('^[pq]($|_)|^FDR$',original)
      if(any(pcol))openxlsx::addStyle(wb,s,
        openxlsx::createStyle(textDecoration=c('bold','italic')),
        rows=1,cols=which(pcol),gridExpand=TRUE,stack=TRUE)
    }
  }
  openxlsx::modifyBaseFont(wb,fontName='Arial',fontSize=10)
  openxlsx::saveWorkbook(wb,file,overwrite=TRUE)
  invisible(file)
}

source_read <- function(file) {
  idx <- data.table::as.data.table(openxlsx::read.xlsx(file,'Contents',check.names=FALSE))
  read <- function(id) {
    a <- idx[ID==id]
    if(a$Kind=='null')return(NULL)
    if(a$Kind=='sf')return(sf::st_read(file.path(dirname(file),a$Sheet),quiet=TRUE))
    meta <- if(is.na(a$Metadata)||!nzchar(a$Metadata))list() else
      jsonlite::fromJSON(a$Metadata,simplifyVector=FALSE)
    if(a$Kind=='list') {
      kids <- idx[Parent==id,ID]
      x <- lapply(kids,read)
      if(!is.null(meta$names))names(x) <- unlist(meta$names)
      if(!is.na(a$Class)&&nzchar(a$Class))class(x) <- strsplit(a$Class,'|',fixed=TRUE)[[1]]
      return(x)
    }
    x <- openxlsx::read.xlsx(file,sheet=a$Sheet,check.names=FALSE,
      skipEmptyRows=FALSE,skipEmptyCols=FALSE,na.strings='<NA>')
    if(nrow(x)<meta$nrow)x <- x[seq_len(meta$nrow),,drop=FALSE]
    names(x) <- unlist(meta$columns)
    for(j in seq_along(x)) {
      info <- meta$types[[j]]
      cl <- unlist(info$class)
      v <- x[[j]]
      if(is.character(v))v[v=='<EMPTY>'] <- ''
      if('factor'%in%cl)v <- factor(v,levels=unlist(info$levels),ordered='ordered'%in%cl)
      else if(any(c('Date','IDate')%in%cl)) {
        v <- as.Date(as.numeric(v),origin='1970-01-01')
        if('IDate'%in%cl)v <- data.table::as.IDate(v)
      } else if(info$type=='integer')v <- as.integer(v)
      else if(info$type=='double')v <- as.numeric(v)
      else if(info$type=='logical')v <- as.logical(v)
      else if(info$type=='character')v <- as.character(v)
      if(names(x)[j]%in%unlist(meta$list_columns))v <- lapply(v,jsonlite::fromJSON)
      x[[j]] <- v
    }
    if(a$Kind=='vector') {
      v <- x[[1]]
      if(!is.null(meta$names))names(v) <- unlist(meta$names)
      return(v)
    }
    if(a$Kind=='matrix') {
      x <- as.matrix(x)
      if(!is.null(meta$rownames))rownames(x) <- unlist(meta$rownames)
      return(x)
    }
    data.table::as.data.table(x)
  }
  read(1L)
}

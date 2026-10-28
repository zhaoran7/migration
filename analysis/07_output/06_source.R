# 1. Panel tables
source_names <- function(d) {
  d <- data.table::copy(data.table::as.data.table(d))
  labels <- c(outcome='Outcome',exposure='Exposure',term='Comparison',group='Group',
    reference='Reference',level='Modifier group',interaction='Interaction',cohort='Cohort',sample='Sample',
    feature='Feature ID',Gene='Gene',comparison='Comparison',profile='Profile',
    prs='PGS ID',prs_label='PGS trait',tdi4='TDI group',
    Observed.Intervals='Observed intervals',Recorded.Fraction='Recorded fraction',Cumulative.Share='Cumulative share',
    precision='Conditional 95% CI width (log HR)',
    origin='Birthplace',destination='Recruitment',variable='Measure',component='Component',
    md='MD',dtdi=paste0(intToUtf8(916),'TDI'),ttwa='Migration cross TTWA',
    camsis='CAMSIS',income='Income',edu.years='Education years',
    TDI='Recruitment TDI','Migration distance'='MD',
    x='X',y='Y',value='Exposure value (SD)',Row='Trait',Column='Trait',MRI_feature='MRI feature',
    md_z='MD',dtdi_z=paste0(intToUtf8(916),'TDI'),density_z='Density change',cross_ttwa='Migration cross TTWA',
    MD_low='MD 25th percentile (km)',MD_high='MD 75th percentile (km)',
    Delta='Delta TDI median',Delta_low='Delta TDI 25th percentile',Delta_high='Delta TDI 75th percentile',
    Density_low='Density change 25th percentile',Density_high='Density change 75th percentile',
    Cross='Cross TTWA',Route='Rural-urban route',SES_n='SES N',
    n='N',N='N',event='Events',events='Events',observations='Observations',
    estimate='Estimate',beta='Beta',se='SE',HR='HR',ratio='RHR',fit='Estimate',
    lo='95% CI lower',hi='95% CI upper',lower='95% CI lower',upper='95% CI upper',
    p='p',q='q (BH)',FDR='q (BH)',p_bonf='p (Bonferroni)',bonferroni='p (Bonferroni)',
    Bonferroni='p (Bonferroni)',p_overall='p overall',p_nonlinear='p non-linear',
    p_interaction='p interaction',interaction_FDR='q interaction',
    context='Location',code='Area code',area='Area',conditional_lower='95% conditional CI lower',
    conditional_upper='95% conditional CI upper',p_approx='p approximate',q_approx='q approximate',
    Trait='Exposure',Cell='Cell type',Cells='Cells',Atlas='Atlas',Tissue='Tissue',
    proportion='Mediation proportion',proportion_lower='Mediation proportion 95% CI lower',
    proportion_upper='Mediation proportion 95% CI upper',indirect='Indirect effect',
    indirect_lower='Indirect effect 95% CI lower',indirect_upper='Indirect effect 95% CI upper',
    a='Path A beta',se_a='Path A SE',b='Path B beta',se_b='Path B SE',
    p_a='p Path A',p_b='p Path B',p_a_bonf='p Path A (Bonferroni)',
    p_b_bonf='p Path B (Bonferroni)',p_indirect='p indirect',
    p_indirect_bonf='p indirect (Bonferroni)',total='Total effect',direct='Direct effect')
  for(v in names(d)) {
    if(is.factor(d[[v]]))data.table::set(d,j=v,value=as.character(d[[v]]))
    if(is.list(d[[v]]))stop('Nested column in submission source data: ',v)
    if(is.numeric(d[[v]]))data.table::set(d,which(!is.finite(d[[v]])),v,NA_real_)
  }
  stopifnot(!any(names(d)%in%c('eid','IID','FID','pidp','HHIDPN')))
  old <- names(d);new <- unname(labels[old]);missing <- is.na(new)
  if('exposure'%in%old)new[old=='x'] <- 'Exposure value'
  new[missing] <- gsub('_',' ',old[missing],fixed=TRUE)
  data.table::setnames(d,make.unique(new,sep=' '))
  readable <- c(md='MD',md_z='MD',dtdi=paste0(intToUtf8(916),'TDI'),dtdi_z=paste0(intToUtf8(916),'TDI'),
    ttwa='Migration cross TTWA',cross_ttwa='Migration cross TTWA',density_z='Density change',
    density_shift='Density change',camsis='CAMSIS',inc='Income',income='Income',
    edu.years='Education years',PC1='PC1',PC2='PC2',PC3='PC3',TDI='Recruitment TDI','Migration distance'='MD')
  for(v in intersect(c('Exposure','Measure','Trait','Component'),names(d)))if(is.character(d[[v]])) {
    hit <- d[[v]]%in%names(readable)
    d[[v]][hit] <- unname(readable[d[[v]][hit]])
  }
  for(v in intersect(c('Outcome','Exposure','Location','Group','Comparison'),names(d)))if(is.character(d[[v]])) {
    d[[v]] <- gsub('CAD/MI','CAD',d[[v]],fixed=TRUE)
    d[[v]] <- gsub('Type 2 diabetes','T2D',d[[v]],fixed=TRUE)
    d[[v]] <- gsub('Heart failure','HF',d[[v]],fixed=TRUE)
    d[[v]] <- gsub('Atrial fibrillation','AF',d[[v]],fixed=TRUE)
    d[[v]] <- gsub('Anxiety disorder','Anxiety',d[[v]],fixed=TRUE)
    if(v=='Location')d[[v]] <- ifelse(d[[v]]=='origin','Birthplace',ifelse(d[[v]]=='destination','Recruitment',d[[v]]))
  }
  d
}

source_block <- function(title,d) list(title=title,data=source_names(d))

source_matrix <- function(title,d,row,column,value) {
  d <- data.table::copy(data.table::as.data.table(d))
  rows <- unique(as.character(d[[row]]));cols <- unique(as.character(d[[column]]))
  stopifnot(!anyDuplicated(d[,c(row,column),with=FALSE]))
  m <- matrix(NA_real_,length(rows),length(cols),dimnames=list(rows,cols))
  m[cbind(match(d[[row]],rows),match(d[[column]],cols))] <- d[[value]]
  table <- data.table::data.table(Label=rows)
  table <- cbind(table,data.table::as.data.table(m))
  data.table::setnames(table,1L,row)
  source_block(title,table)
}

source_curves <- function(r,selected,bounds) {
  data.table::rbindlist(lapply(seq_len(nrow(selected)),function(i) {
    s <- selected[i];d <- data.table::copy(r$Effects[interaction==s$interaction & outcome==s$outcome])
    if(grepl('age$',s$interaction))d <- d[order%in%range(order)]
    b <- bounds[exposure==if(grepl('^MD',s$interaction))'md_z' else 'dtdi_z']
    g <- d[,.(Exposure_SD=seq(b$lower,b$upper,length.out=101)),by=.(level,beta,se)]
    g[,`:=`(estimate=beta*Exposure_SD,lower=beta*Exposure_SD-1.96*se*abs(Exposure_SD),
      upper=beta*Exposure_SD+1.96*se*abs(Exposure_SD))]
    if(s$outcome!='SES')g[,c('estimate','lower','upper'):=lapply(.SD,exp),.SDcols=c('estimate','lower','upper')]
    test <- r$Tests[interaction==s$interaction & outcome==s$outcome][1]
    g[,`:=`(interaction=s$interaction,outcome=s$outcome,p_interaction=test$p,q_interaction=test$FDR)]
    g
  }),fill=TRUE)
}

# 2. Main figures
source_panels <- function(id,x) {
  B <- source_block;M <- source_matrix
  delta <- paste0(intToUtf8(916),'TDI');views <- c('MD',delta)
  take <- function(d,cols)d[,intersect(cols,names(d)),with=FALSE]
  forest <- c('outcome','exposure','term','strat','group','sample','cohort','n','observations','event','events',
    'beta','se','HR','estimate','lo','hi','lower','upper','p','FDR','q')
  association <- function(d)take(d,forest)
  if(id=='Fig1') {
    r <- x$description;p <- x$map$plotting
    d <- r$density[section%in%c('MD RINT','Density change')]
    areas <- unique(p$poly[,.(code,Recruitment_TDI_median=tdi_med)])
    areas <- merge(p$nodes[,.(code,name,in_n,out_n,total)],areas,by='code',all.x=TRUE,sort=FALSE)
    return(list(Fig1b=list(B('Density curves',d[,.(Measure=section,N=n,X=x,Density=density)]),
      B('Group sizes',r$composition[,.(Measure=Variable,Group=Route,N,Proportion=Share)])),
      Fig1c=list(B('TTWA summaries',areas),B('Displayed migration flows',p$edges[,.(Birthplace_TTWA=from,Recruitment_TTWA=to,Birthplace=from_name,Recruitment=to_name,N)]))))
  }
  if(id=='Fig2') {
    r <- x$result;d <- data.table::copy(x$a$d);labs0 <- x$a$labs
    d[,`:=`(Row=labs0[row],Column=labs0[col])]
    symmetric <- function(part,value) {
      z <- d[half==part & row!=col,c('Row','Column',value),with=FALSE]
      back <- data.table::copy(z);data.table::setnames(back,c('Row','Column'),c('Column','Row'))
      z <- unique(data.table::rbindlist(list(z,back),use.names=TRUE))
      z <- merge(data.table::CJ(Row=labs0,Column=labs0),z,by=c('Row','Column'),all.x=TRUE,sort=FALSE)
      z[order(match(Row,labs0),match(Column,labs0))]
    }
    a <- list(M('Phenotypic correlation',symmetric('Phenotypic','val'),'Row','Column','val'),
      M('LDSC genetic correlation',symmetric('LDSC','val'),'Row','Column','val'),
      M('Phenotypic correlation p',symmetric('Phenotypic','p'),'Row','Column','p'),
      M('LDSC genetic correlation p',symmetric('LDSC','p'),'Row','Column','p'),
      B('SNP heritability',x$a$h[,.(Trait=label,Heritability=h2,SE=h2_se,`95% CI lower`=lo,`95% CI upper`=hi)]))
    b <- list(M('Adjusted SES by TDI group',r$ses$tdi,'origin','destination','fit'),
      B('Adjusted SES by TDI group: confidence intervals',take(r$ses$tdi,c('group','origin','destination','fit','lo','hi'))),
      M('Adjusted SES by rural-urban route',r$ses$urban,'origin','destination','fit'),
      B('Adjusted SES by rural-urban route: confidence intervals',take(r$ses$urban,c('group','n','origin','destination','fit','lo','hi'))))
    z <- data.table::rbindlist(list(r$ses$main,r$ses$density),fill=TRUE);z[,FDR:=p.adjust(p,'BH')]
    b <- c(b,list(B('SES associations',association(z)),B('MD-SES curve',r$ses$md[,.(MD_km=x,SES_SD=fit,lower=lo,upper=hi)])))
    c <- list(B('MD and cross-TTWA associations',association(r$forest$main)),
      B('TDI associations and group contrasts',association(r$forest$tdi)),
      B('Density and rural-urban associations',association(r$direction$results)),
      B('MD by occupational physical demand',association(r$strata[strat=='work'])))
    z <- data.table::copy(r$interaction$Terms)
    z[,level:=gsub('md_z|dtdi_z|age_z|sex2|tdi4|work_g|urban_route|:','',term)]
    z[!interaction%in%c('MD x TDI group','MD x physical demand','MD x sex','dTDI x sex','MD x rural-urban route'),level:='']
    z[,key:=paste(interaction,level,sep='|')]
    tests <- r$interaction$Tests[,.(outcome,interaction,omnibus_q=FDR)]
    z <- merge(z,tests,by=c('outcome','interaction'),all.x=TRUE,sort=FALSE)
    z[,eligible:=is.finite(FDR)&FDR<.05&is.finite(omnibus_q)&omnibus_q<.05]
    hkeys <- c('MD x age|','dTDI x age|','MD x sex|Male','dTDI x sex|Male','MD x physical demand|Usually')
    skeys <- c('MD x dTDI|','MD x density change|','MD x rural-urban route|Urban to Rural',
      'MD x rural-urban route|Rural to Urban','MD x rural-urban route|Rural to Rural')
    health <- unique(z[model=='Cox' & eligible==TRUE,outcome])
    z <- z[(model=='Cox' & key%in%hkeys & outcome%in%health)|(outcome=='SES' & key%in%skeys)]
    z[,`:=`(Interaction=key,display=fifelse(eligible,ratio,NA_real_))]
    d <- list(M('Health interaction RHR (blank if not significant)',z[model=='Cox'],'Interaction','outcome','display'),
      M('Health interaction q',z[model=='Cox'],'Interaction','outcome','FDR'),
      M('SES interaction beta (blank if not significant)',z[outcome=='SES'][,display_beta:=fifelse(eligible,beta,NA_real_)],'Interaction','outcome','display_beta'),
      B('SES interaction coefficients',take(z[outcome=='SES'],c('interaction','term','beta','se','lower','upper','p','FDR','omnibus_q'))),
      B('Displayed interaction coefficients',take(z,c('outcome','interaction','term','beta','se','ratio','lower','upper','p','FDR','omnibus_q'))),
      B('Conditional association curves',source_curves(r$interaction,x$selected,x$bounds)))
    return(list(Fig2a=a,Fig2b=b,Fig2c=c,Fig2d=d))
  }
  if(id=='Fig3') {
    r <- x$result;z <- data.table::copy(r$outcomes)
    z[,display:=fifelse(supported,estimate,NA_real_)]
    boxes <- take(r$boxes,c('profile','variable','N','minimum','lower','q1','median','q3','upper','maximum'))
    data.table::setnames(boxes,c('minimum','lower','q1','median','q3','upper','maximum'),
      c('Minimum','Lower whisker','25th percentile','Median','75th percentile','Upper whisker','Maximum'))
    return(list(Fig3a=list(B('Profile definitions and sizes',r$description),
      B('Boxplot summaries',boxes),
      M('Adjusted SES and 10-year risk (%)',z,'outcome','profile','display'),
      B('Adjusted SES and risk: uncertainty',take(z,c('profile','outcome','n','events','supported','estimate','lower','upper','p','q','measure')))),
      Fig3b=list(B('SES and 10-year risk by profile',take(r$points,c('profile','outcome','n','events','ses_n','available','estimate','ses','N',
        'ses_reference','risk_reference','ses_difference','risk_difference','pattern'))))))
  }
  if(id=='Fig4') {
    r <- x$result;t <- r$tests;a <- r$area
    panels <- list(Fig4a=list(B('Regional variance tests',take(t,c('outcome','context','n','events','variance','MHR','LR','df','p','q')))))
    g <- t[q<.05,.(contexts=list(intersect(c('origin','destination'),context)),maps=.N),by=.(id,outcome)]
    g <- data.table::rbindlist(list(g[maps==2][order(id)],g[maps==1][order(-id)]))
    for(i in seq_len(nrow(g)))panels[[paste0('Fig4',letters[i+1])]] <- list(B('Regional HR estimates',take(a[id==g$id[i] & context%in%g$contexts[[i]]],
      c('outcome','context','code','area','n','events','person_years','beta','se','HR','conditional_lower','conditional_upper','p_approx','q_approx'))),
      B('Regional variance test',take(t[id==g$id[i] & context%in%g$contexts[[i]]],c('outcome','context','n','events','variance','MHR','p','q'))))
    return(panels)
  }
  if(id=='Fig5') {
    r <- x$result
    return(list(Fig5a=list(B('Participants by wave',r$Waves[,.(Cohort=cohort,Wave=wave,N=n)]),
      B('Migrant status',r$Status[,.(Cohort=cohort,Group=ifelse(mover==1,'Migrant','Non-migrant'),N=n)]),
      B('Migration times',r$Counts[,.(Cohort=cohort,MT=count,N=n,Proportion=fraction)]),
      B('Cumulative MD density',r$Density[,.(Cohort=cohort,Log_cumulative_MD=x,Density=y)])),
      Fig5b=list(B('SES and mortality associations',association(r$Main))),
      Fig5c=list(B('SES curves',take(r$Curves[outcome=='SES'],c('cohort','exposure','x','reference','estimate','lo','hi','n','observations','p_overall','p_nonlinear')))),
      Fig5d=list(B('Mortality curves',take(r$Curves[outcome=='Mortality'],c('cohort','exposure','x','reference','estimate','lo','hi','n','events','p_overall','p_nonlinear')))),
      Fig5e=list(B('Within-person SES curves',take(r$Within,c('cohort','exposure','x','reference','estimate','lo','hi','n','observations','p_overall','p_nonlinear'))))))
  }
  if(id=='Fig6') {
    r <- x$result;protein <- r$Protein$association[exposure%in%views]
    mri <- r$MRI$association[exposure%in%views];mri <- data.table::copy(mri[feature%in%mri[p_bonf<.05,unique(feature)]])
    mri[,`:=`(Z=beta/se,display=fifelse(p_bonf<.05,beta/se,NA_real_))]
    mri[,MRI_feature:=paste0(`Display name`,' [',feature,']')]
    med <- data.table::copy(r$Protein$mediation[selected==TRUE & exposure%in%views])
    rank <- med[,.(proportion=max(proportion),p=min(p_indirect_bonf)),by=feature][order(-proportion,p,feature)]
    md <- rank[feature%in%med[exposure=='MD',feature],feature]
    ids <- head(c(md,setdiff(rank$feature,md)),40)
    med <- med[feature%in%ids]
    return(list(Fig6a=list(B('Protein associations',take(protein,c('exposure','feature','Display name','Gene','Chromosome','Start','End','n','beta','se','p','p_bonf')))),
      Fig6b=list(M('MRI Z scores (blank if not significant)',mri,'MRI_feature','exposure','display'),
        B('MRI associations',take(mri,c('exposure','feature','Display name','MRI class','Atlas','Measure','Region','Hemisphere','n','beta','se','Z','p','p_bonf')))),
      Fig6c=list(B('Displayed protein mediation paths',take(med,c('exposure','feature','Display name','Gene','outcome','n','event',
        'a','se_a','p_a','p_a_bonf','b','se_b','p_b','p_b_bonf','total','direct','indirect','indirect_lower','indirect_upper','p_indirect','p_indirect_bonf',
        'proportion','proportion_lower','proportion_upper'))))))
  }

# 3. Extended Data
  if(id=='ED1')return(list())
  if(id=='ED2')return(list(ED2c=list(B('Median Recruitment TDI by TTWA',x$result$tdi)),
    ED2d=list(B('Official rural-urban classifications',x$result$urban_source)),
    ED2f=list(B('2011 public population grid',x$result$pop))))
  if(id=='ED3') {
    result <- list(ED3b=list(B('MD Q-Q data',x$qc$qq$md)),ED3c=list(B('Delta TDI Q-Q data',x$qc$qq$dtdi)))
    for(i in 1:2) {
      z <- x$comparison[[c('md','dtdi')[i]]]
      result[[paste0('ED3',letters[3+i])]] <- list(B('Matched SNP Z-score bins',z$bins[,.(Reference_Z=x,Current_Z=y,Variants=variants)]),
        B('Genome-wide significant SNPs',take(z$signals,c('CHR','POS','SNP','EA','NEA','z_current','z_reference','p_current','p_reference'))),
        B('LDSC comparison',take(z$metrics,c('reference','matched','rg','rg_se'))))
    }
    return(result)
  }
  if(id%in%c('ED4','ED5','ED6')) {
    keys <- switch(id,ED4=c('md_q','ttwa2'),ED5=c('dtdi_q','tdi4'),ED6=c('density_shift_q','urban4'))
    result <- list()
    for(i in 1:2) {
      key <- keys[i];d <- data.table::copy(x$prs$bars[grouping==key])
      d[,`:=`(lower=estimate-1.96*se,upper=estimate+1.96*se)]
      result[[paste0(id,letters[2*i-1])]] <- list(B('Mean standardized PGS',take(d,c('prs_label','group','n','estimate','se','lower','upper','p_group','q_group'))))
      a <- data.table::copy(x$prs$associations[grouping==key])
      lookup <- unique(x$prs$bars[,.(prs,prs_label)])
      a[,prs_label:=lookup$prs_label[match(prs,lookup$prs)]]
      result[[paste0(id,letters[2*i])]] <- list(B('PGS associations',take(a,c('prs_label','prs','group','reference','n','beta','se','lower','upper','p','q'))))
    }
    return(result)
  }
  if(id=='ED7') {
    r <- x$result;result <- list()
    for(i in 1:3)result[[paste0(id,letters[i])]] <- list(B('Sensitivity associations',association(r$Effects[exposure==c(views,'Migration cross TTWA')[i]])))
    result$ED7d <- list(B('Cross-TTWA with additional MD adjustment',association(data.table::rbindlist(list(r$Effects[group=='Baseline' & exposure=='Migration cross TTWA'],r$TTWA),fill=TRUE))))
    for(i in 1:2)result[[paste0(id,letters[i+4])]] <- list(B('Spline curve',r$Curves[exposure==views[i]]),
      B('Spline association tests',take(r$Tests[exposure==views[i]],c('outcome','exposure','n','event','p_overall','p_nonlinear','q_overall','q_nonlinear'))))
    return(result)
  }
  if(id=='ED8') {
    r <- x$result;selected <- r$Tests[FDR<.05][!x$selected,on=.(interaction,outcome)]
    return(list(ED8a=list(B('MD by TDI group',association(r$Strata[strat=='tdi']))),
      ED8b=list(B('MD by sex',association(r$Extra[strat=='sex']))),
      ED8c=list(B('MD by age',association(r$Extra[strat=='age']))),
      ED8d=list(B('Interaction tests',take(selected,c('outcome','interaction','n','event','df','p','FDR'))),
        B('Conditional association curves',source_curves(r,selected,x$bounds)))))
  }
  if(id=='ED9') {
    r <- x$result
    return(list(ED9a=list(B('Analysis sample sizes',take(r$Effects[group%in%c('Core','Adjusted')],c('cohort','outcome','exposure','group','n','observations','event')))),
      ED9b=list(B('Adjusted associations',association(r$Effects[group=='Adjusted']))),
      ED9c=list(B('SES curves',take(r$Curves[outcome=='SES'],c('cohort','exposure','x','reference','estimate','lo','hi','n','p_overall','p_nonlinear')))),
      ED9d=list(B('Mortality curves',take(r$Curves[outcome=='Mortality'],c('cohort','exposure','x','reference','estimate','lo','hi','n','events','p_overall','p_nonlinear')))),
      ED9e=list(B('Within-person SES curves',take(r$Within,c('cohort','exposure','x','reference','estimate','lo','hi','n','p_overall','p_nonlinear'))))))
  }
  if(id=='ED10') {
    r <- x$result;d <- r$cells[Display==TRUE];d[,Row:=paste(Atlas,Cell,sep=': ')]
    spatial <- data.table::dcast(r$spatial,spot+x+y+Tissue~Trait,value.var='p')
    return(list(ED10a=list(M('Tissue enrichment: -log10(p)',r$tissue,'Tissue','Trait','score'),
      M('Tissue enrichment: p',r$tissue,'Tissue','Trait','p'),M('Tissue enrichment: Bonferroni p',r$tissue,'Tissue','Trait','Bonferroni'),
      B('Spot-level gsMap associations',spatial)),
      ED10b=list(M('Cell-type scDRS Z',d,'Row','Trait','Z'),M('Cell-type scDRS p',d,'Row','Trait','p'),
        M('Cell-type scDRS Bonferroni p',d,'Row','Trait','Bonferroni'),B('Displayed cell-type counts',r$counts[,.(Atlas,Cell,Cells)]))))
  }

# 4. Supplementary figures
  if(id=='S1') {
    r <- x$result;stages <- unique(r$nodes[,.(axis,stage)])
    links <- merge(r$links[,.(axis,quartile,from,to,N)],stages,by='axis',all.x=TRUE,sort=FALSE)
    links[,To_stage:=stages$stage[match(axis+1,stages$axis)]]
    return(list(S1=list(B('Adjacent-stage flows',links[,.(From_stage=stage,To_stage,MD_quartile=quartile,From_group=from,To_group=to,N)]),
      B('Stage group sizes',r$nodes[,.(Stage=stage,Group=group,N)]),B('Quartile cutoffs',r$cuts))))
  }
  if(id=='S2') {
    r <- x$supp$base$ses;d <- data.table::copy(r$cor);data.table::setnames(d,c('Row','Column','r'))
    return(list(S2a=list(B('SES component availability',r$available)),S2b=list(M('Component Pearson correlation',d,'Row','Column','r')),
      S2c=list(B('PC1 loadings',r$loading)),S2d=list(B('Variance explained',r$variance)),
      S2e=list(M('SES component associations: beta',x$result$ses$comp,'component','exposure','estimate'),
        B('SES component association results',association(x$result$ses$comp)[,Component:=x$result$ses$comp$component]))))
  }
  if(id=='S3') {
    r <- x$supp$base$geo;d <- data.table::copy(r$cor);data.table::setnames(d,c('Row','Column','rho'))
    z <- data.table::copy(r$routes);z[,Proportion:=n/sum(n),by=tdi4]
    return(list(S3a=list(B('TDI matching-distance cumulative distribution',r$nearest[,.(Distance_km=expm1(X),Cumulative_proportion=Y)])),
      S3b=list(B('Joint TDI histogram bins',r$tdi)),S3c=list(M('Spearman correlation',d,'Row','Column','rho')),
      S3d=list(M('Rural-urban proportions within TDI groups',z,'tdi4','urban_route','Proportion')),
      S3e=list(B('TDI matching sensitivity: health',association(x$supp$matching$effects[outcome!='SES' & term==delta]))),
      S3f=list(B('TDI matching sensitivity: SES',take(x$supp$matching$effects[outcome=='SES' & term!='MD'],c('sample','outcome','term','n','beta','se','estimate','lower','upper','p','q'))))))
  }
  if(id=='S4') {
    r <- x$supp$adjustment
    make <- function(z)B('Model comparison',take(z,c('outcome','exposure','adjustment','parameter','n','event','beta','se','estimate','lower','upper','p','q')))
    return(list(S4a=list(make(r[outcome=='SES'])),S4b=list(make(r[outcome!='SES' & exposure=='MD'])),
      S4c=list(make(r[outcome!='SES' & exposure==delta]))))
  }
  if(id=='S5') {
    r <- data.table::copy(x$supp$base$regional$area)
    r[,precision:=2*1.96*se]
    d <- r[,.(Areas=.N,Adequate_areas=sum(n>=100 & events>=10)),by=.(outcome,context)]
    d[,Proportion:=Adequate_areas/Areas]
    return(list(S5a=list(B('Participants per estimated area',take(r,c('outcome','context','code','area','n')))),
      S5b=list(B('Areas meeting support thresholds',d)),
      S5c=list(B('Events and interval widths',take(r,c('outcome','context','code','area','events','precision'))))))
  }
  if(id=='S6') {
    r <- x$result;result <- list()
    for(i in 1:2)result[[paste0('S6',letters[i])]] <- list(B('Regional SES deviations',r$areas[context==c('origin','destination')[i]]),
      B('Regional variance test',r$tests[context==c('origin','destination')[i]]))
    for(i in 1:2)result[[paste0('S6',letters[i+2])]] <- list(B('Regional slope distribution',r$slopes[[c('md_z','dtdi_z')[i]]]$metrics))
    return(result)
  }
  if(id=='S7')return(stats::setNames(lapply(letters[1:4],function(k)list(B('Residential-history coverage',x[[k]]))),paste0('S7',letters[1:4])))
  if(id=='S8') {
    r <- x$supp;d <- data.table::copy(r$omics$Protein$contrasts)
    counts <- d[,.(Negative=sum(bonferroni<.05 & beta<0),Positive=sum(bonferroni<.05 & beta>0)),by=comparison]
    features <- unique(d[bonferroni<.05][order(p)][,head(.SD,6),by=comparison]$feature)
    small <- r$pattern$map[,.(n=.N),by=pattern][n<=12,pattern]
    features <- unique(c(features,r$pattern$map[pattern%in%small,feature]))
    selected <- d[feature%in%features]
    selected[,Protein:=paste0(ifelse(!is.na(Gene)&nzchar(Gene),Gene,`Display name`),' [',feature,']')]
    return(list(S8a=list(B('Significant protein counts',counts)),
      S8b=list(B('Four-group protein contrasts',take(d,c('comparison','feature','Display name','Gene','Chromosome','Start','End','n','beta','se','lower','upper','p','bonferroni')))),
      S8c=list(M('Scaled group means',r$pattern$profiles,'feature','group','shape'),B('Protein pattern membership',r$pattern$map)),
      S8d=list(M('Displayed protein contrast beta',selected,'Protein','comparison','beta'),
        M('Displayed protein contrast Bonferroni p',selected,'Protein','comparison','bonferroni'),
        B('Displayed contrast results',take(selected,c('feature','Gene','comparison','n','beta','se','lower','upper','p','bonferroni'))))))
  }
  if(id=='S9') {
    result <- list()
    for(i in 1:2) {
      atlas <- c('Heart','Liver')[i];d <- x$cells[Atlas==atlas]
      result[[paste0('S9',letters[2*i-1])]] <- list(B('Reference cell-type counts',unique(d[,.(Atlas,Cell,Cells)])))
      result[[paste0('S9',letters[2*i])]] <- list(M('scDRS Z',d,'Cell','Trait','Z'),M('scDRS p',d,'Cell','Trait','p'),
        M('scDRS Bonferroni p',d,'Cell','Trait','Bonferroni'))
    }
    return(result)
  }
  if(id=='S10') {
    result <- list()
    for(i in 1:2) {
      atlas <- c('cs8','dlpfc')[i]
      d <- data.table::rbindlist(lapply(c('MD','DS'),function(e) {
        z <- data.table::copy(x$tissues[[atlas]][[e]]);z[,Trait:=if(e=='MD')'MD' else delta];z
      }))
      d[,`:=`(Bonferroni=p.adjust(p_cauchy,'bonferroni'),Score=pmin(15,-log10(p_cauchy))),by=Trait]
      spots <- data.table::rbindlist(lapply(c('MD','DS'),function(e) {
        z <- merge(x$coordinates[[atlas]],x$spatial[[atlas]][[e]],by='spot',sort=FALSE)
        z[,Trait:=if(e=='MD')'MD' else delta];z
      }))
      result[[paste0('S10',letters[i])]] <- list(M('Annotation enrichment: -log10(p)',d,'annotation','Trait','Score'),
        M('Annotation enrichment: p',d,'annotation','Trait','p_cauchy'),
        M('Annotation enrichment: Bonferroni p',d,'annotation','Trait','Bonferroni'),
        B('Spot-level associations',data.table::dcast(spots,spot+x+y+annotation~Trait,value.var=c('beta','se','z','p'))))
    }
    return(result)
  }
  stop('Unknown source-data figure: ',id)
}

# 5. Submission workbooks
source_workbook <- function(panels,file) {
  if(!length(panels)){if(file.exists(file))unlink(file);return(invisible(NULL))}
  dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
  wb <- openxlsx::createWorkbook(creator='MOVE')
  openxlsx::modifyBaseFont(wb,fontName='Arial',fontSize=10)
  header <- openxlsx::createStyle(textDecoration='bold',fgFill='#F2F2F2',border='bottom',borderColour='#A9AFB4',wrapText=TRUE)
  title <- openxlsx::createStyle(textDecoration='bold')
  for(panel in names(panels)) {
    openxlsx::addWorksheet(wb,panel,gridLines=FALSE)
    row <- 1L
    for(block in panels[[panel]]) {
      d <- as.data.frame(block$data)
      openxlsx::writeData(wb,panel,block$title,startRow=row,colNames=FALSE)
      openxlsx::addStyle(wb,panel,title,rows=row,cols=1)
      openxlsx::writeData(wb,panel,d,startRow=row+1,headerStyle=header)
      openxlsx::setRowHeights(wb,panel,row+1,30)
      pcols <- which(grepl('^p($| )|^q($| )',names(d)))
      probability <- pcols
      if(grepl('correlation p$|scDRS p$|Bonferroni p$|enrichment: p$|interaction q$',block$title))probability <- which(vapply(d,is.numeric,logical(1)))
      if(grepl('Spot-level gsMap',block$title))probability <- which(vapply(d,is.numeric,logical(1)) & seq_along(d)>4)
      if(length(pcols))openxlsx::addStyle(wb,panel,openxlsx::createStyle(textDecoration=c('bold','italic')),rows=row+1,cols=pcols,gridExpand=TRUE,stack=TRUE)
      for(j in which(vapply(d,is.numeric,logical(1))))if(nrow(d))openxlsx::addStyle(wb,panel,
        openxlsx::createStyle(numFmt=if(j%in%probability)'0.00E+00' else if(all(is.na(d[[j]])|d[[j]]==trunc(d[[j]])))'0' else '0.0000'),
        rows=row+2+seq_len(nrow(d))-1,cols=j,stack=TRUE)
      row <- row+nrow(d)+5L
    }
    openxlsx::freezePane(wb,panel,firstActiveRow=3)
    openxlsx::setColWidths(wb,panel,1:40,18)
    openxlsx::setColWidths(wb,panel,1,30)
  }
  openxlsx::saveWorkbook(wb,file,overwrite=TRUE)
  invisible(file)
}

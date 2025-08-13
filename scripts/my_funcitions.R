theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Times"),
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 14),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          # plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 14, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

mg_limma_DEG_use=function(exp,group,ulab,dlab=NULL,paired=FALSE){
  library(limma)
  ind1=which(group==ulab)
  if(is.null(dlab)){
    ind2=which(group!=ulab)
  }else{
    ind2=which(group==dlab)
  }
  sml <- c(rep('G1',length(ind1)),rep('G0',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  if(paired){
    pairinfo=factor(rep(1:length(ind1),each=2))
    design <- model.matrix(~fl+0+pairinfo)
    colnames(design)[1:length(levels(fl))] <- levels(fl)
  }else{
    design <- model.matrix(~fl+0)
    colnames(design)[1:length(levels(fl))] <- levels(fl)
  }
  cont.matrix<-limma::makeContrasts(contrasts='G1-G0',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", coef=1,sort.by="B", number=nrow(eset))
  regulated=ifelse(tT$logFC>0,'Up','Down')
  lfcs=c(log2(1.2),log2(1.3),log2(1.5),1)
  all.deg.cnt=cbind()
  for(lfc in lfcs){
    deg1=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.05)]
    deg2=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.01)]
    deg3=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.05)]
    deg4=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.01)]
    all.deg.cnt=cbind(all.deg.cnt,c(paste0(sum(deg1=='Up'),'|',sum(deg1=='Down'))
                                    ,paste0(sum(deg2=='Up'),'|',sum(deg2=='Down'))
                                    ,paste0(sum(deg3=='Up'),'|',sum(deg3=='Down'))
                                    ,paste0(sum(deg4=='Up'),'|',sum(deg4=='Down'))))
  }
  row.names(all.deg.cnt)=c('p<0.05','p<0.01','FDR<0.05','FDR<0.01')
  colnames(all.deg.cnt)=paste0(c('1.2','1.3','1.5','2'),'-fold')
  return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=all.deg.cnt))
}

get_DEG=function(df_deg,p.cutoff=0.05,logfc.cutoff=1){
  df.deg.res=df_deg$DEG
  df.deg.sig=df.deg.res[which(df.deg.res$adj.P.Val<p.cutoff & abs(df.deg.res$logFC)>logfc.cutoff),]
}

ggplotKMCox=function(dat,title='Groups',labs=NULL,add_text=NULL,pal=NULL){
  library(ggplot2)
  library(ggsci)
  colnames(dat)=c('time','status','groups')
  #sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  #print((sdf))
  #summary(sdf)
  #p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  if(is.null(pal)){
    surv=survminer::ggsurvplot(sf, data = dat
                               , palette =pal_lancet()(9)[c(2,4,3,1,5:6)]
                               ,pval = TRUE, surv.median.line='hv'
                               # ,conf.int = T
                               , xlab="Time(years)"
                               ,linetype = "strata"
                               ,conf.int.style ='step'
                               , pval.coord=c(0, 0.2) #Add p-value 
                               , risk.table.title=""
                               , risk.table = TRUE
                               , legend.title = title
                               , legend.labs = labs
    )
    
  }else{
    if(pal=='npg'){
      surv=survminer::ggsurvplot(sf, data = dat
                                 , palette = pal_npg()(9)[c(1,2,3,4:9)] #jco palette 
                                 ,pval = TRUE, surv.median.line='hv'
                                 # ,conf.int = T
                                 , xlab="Time(years)"
                                 ,linetype = "strata"
                                 ,conf.int.style ='step'
                                 , pval.coord=c(0, 0.2), #Add p-value 
                                 risk.table = TRUE, 
                                 legend.title = title
                                 ,legend.labs = labs
      )
    }else{
      surv=survminer::ggsurvplot(sf, data = dat
                                 , palette = pal #jco palette 
                                 ,pval = TRUE, surv.median.line='hv'
                                 # ,conf.int = T
                                 ,linetype = "strata"
                                 , xlab="Time(years)"
                                 ,conf.int.style ='step'
                                 , pval.coord=c(0, 0.2), #Add p-value 
                                 risk.table = TRUE, 
                                 legend.title = title
                                 ,legend.labs = labs
      )
    }
    
  }
  
  p1=surv$plot+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                ,axis.text.x=element_blank()
                                ,axis.title.x=element_blank()
                                ,plot.margin=unit(c(0, 0, 0, 0), "inches")
                                #,axis.title.y=element_blank()
                                ,legend.position=c(1,1), legend.justification=c(1,1)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,legend.title = element_text(family="Times",face="plain",size = 10)
                                ,legend.text = element_text(family="Times",face="plain",size = 10))
  #p1=p1+text()
  #tms=data.frame(Group=tms.gp,value=tms.tps,Attribute=rep(data_m[1,1],length(tms.gp))
  #               ,ymax=rep(max(ylim),length(tms.gp)))
  #p4=p4+geom_text(data=tms,aes(x=Group, y=ymax, label=value),color="yellow")
  if(!is.null(add_text)){
    text.tb=surv$data.survplot[1,]
    text.tb[1,1]=0
    text.tb[1,5]=0
    text.tb$Text=add_text
    p1=p1+geom_text(data=text.tb,aes(x=time, y=surv, label=Text),color="yellow",hjust =0)
  }
  
  p2=surv$table+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                 #,axis.text.x=element_blank()
                                 #,axis.title.x=element_blank()
                                 #,axis.title.y=element_blank()
                                 ,plot.margin=unit(c(0, 0, 0, 0), "inches")
                                 ,plot.title=element_blank()
                                 ,legend.position=c(1,1), legend.justification=c(1,1)
                                 #,legend.background = element_rect(fill = NA, colour = NA)
                                 ,legend.title = element_text(family="Times",face="plain",size = 10)
                                 ,legend.text = element_text(family="Times",face="plain",size = 10))
  p1=surv$plot
  p2=surv$table
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(0.9,0.3),align = "v")
  return(g2)
}

plotMutiBar_tmp=function(dat,ist=F,margin=T,xlb='',ylb='',fill.color=mycolors1,lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T,legend.nrow=0){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol,position=position_fill()) 
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_d3()+theme(legend.position = "bottom",legend.key.size = unit(0.2, 'cm'))
  pg=pg+ggsci::scale_fill_d3()+scale_fill_discrete(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  pg=pg+scale_fill_manual(values = fill.color,breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  pg=pg+scale_y_continuous(expand=c(0,0)) ## 新增
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  if(legend.nrow>0){
    print('YES')
    print(legend.nrow)
    pg=pg+guides(fill = guide_legend(nrow = legend.nrow, byrow = TRUE))
  }
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        # g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
        g.tb[i,j]=as.numeric(format(chisq.test(bk_dat[,c(i,j)])$p.value,digits=3))
        
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  # g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]<0.05,'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") +xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(anova p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

useMyCol <- function(platte = NULL,
                     n = NULL,
                     showAll = FALSE){
  #---------------------------------------------------------------
  # Primarily Discrete Palettes
  #---------------------------------------------------------------
  # ============================================================================
  #20-colors
  stallion = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
               "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
               "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
               "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D")
  
  stallion2 = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                "#6E4B9E","#0C727C", "#7E1416","#D8A767")
  
  calm = c("#7DD06F", "#844081","#688EC1", "#C17E73", "#484125",
           "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
           "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D",
           "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736")
  
  kelly = c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020",
            "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
            "#FF7A5C", "#53377A", "#FF8E00", "#B32851", "#F4C800",
            "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16")
  
  #16-colors
  bear = c("#faa818", "#41a30d","#fbdf72", "#367d7d", "#d33502",
           "#6ebcbc", "#37526d","#916848", "#f5b390", "#342739",
           "#bed678","#a6d9ee", "#0d74b6",
           "#60824f","#725ca5", "#e0598b")
  
  #15-colors
  ironMan = c('#371377','#7700FF','#9E0142','#FF0080', '#DC494C',
              "#F88D51","#FAD510","#FFFF5F",'#88CFA4','#238B45',
              "#02401B", "#0AD7D3","#046C9A", "#A2A475", 'grey35')
  
  circus = c("#D52126","#88CCEE", "#FEE52C", "#117733", "#CC61B0",
             "#99C945", "#2F8AC4", "#332288","#E68316", "#661101",
             "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74")
  
  #12-colors
  paired = c("#A6CDE2","#1E78B4","#74C476","#34A047","#F59899","#E11E26",
             "#FCBF6E","#F47E1F","#CAB2D6","#6A3E98","#FAF39B","#B15928")
  
  #11-colors
  grove = c("#1a1334","#01545a","#017351","#03c383","#aad962",
            "#fbbf45","#ef6a32","#ed0345","#a12a5e","#710162","#3B9AB2")
  
  #7-colors
  summerNight = c("#2a7185","#a64027","#fbdf72","#60824f","#9cdff0",
                  "#022336","#725ca5")
  
  #5-colors
  zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
  darjeeling = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
  rushmore = c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" , "#F2300F")
  captain = c("grey","#A1CDE1","#12477C","#EC9274","#67001E")
  
  #---------------------------------------------------------------
  # Primarily Continuous Palettes
  #---------------------------------------------------------------
  #10-colors
  horizon = c('#000075','#2E00FF', '#9408F7', '#C729D6', '#FA4AB5',
              '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60')
  
  #9-colors
  horizonExtra =c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5",
                  "#FD619D","#FF9965","#FFD32B","#FFFC5A")
  
  blueYellow = c("#352A86","#343DAE","#0262E0","#1389D2","#2DB7A3",
                 "#A5BE6A","#F8BA43","#F6DA23","#F8FA0D")
  
  sambaNight = c('#1873CC','#1798E5','#00BFFF','#4AC596','#00CC00',
                 '#A2E700','#FFFF00','#FFD200','#FFA500')
  
  solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC',
                 '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D')
  
  whitePurple = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                  '#8c6bb1','#88419d','#810f7c','#4d004b')
  
  whiteBlue = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                '#3690c0','#0570b0','#045a8d','#023858')
  
  
  comet = c("#E6E7E8","#3A97FF","#8816A7","black")
  
  #7-colors
  greenBlue = c('#e0f3db','#ccebc5','#a8ddb5','#4eb3d3','#2b8cbe',
                '#0868ac','#084081')
  
  #6-colors
  beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30",
            "#F7962E","#FCEE2B")
  
  #5-colors
  coolwarm = c("#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29")
  fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A")
  greyMagma = c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF")
  fireworks2 = c("black", "#2488F0","#7F3F98","#E22929","#FCB31A")
  purpleOrange = c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F")
  
  # ============================================================================
  if(showAll == FALSE){
    if(platte == "stallion"){
      col <- stallion[1:n]
      print(paste0('This palatte have ',length(stallion),' colors!'))
    }else if(platte == "stallion2"){
      col <- stallion2[1:n]
      print(paste0('This palatte have ',length(stallion2),' colors!'))
    }else if(platte == "calm"){
      col <- calm[1:n]
      print(paste0('This palatte have ',length(calm),' colors!'))
    }else if(platte == "kelly"){
      col <- kelly[1:n]
      print(paste0('This palatte have ',length(kelly),' colors!'))
    }else if(platte == "bear"){
      col <- bear[1:n]
      print(paste0('This palatte have ',length(bear),' colors!'))
    }else if(platte == "ironMan"){
      col <- ironMan[1:n]
      print(paste0('This palatte have ',length(ironMan),' colors!'))
    }else if(platte == "circus"){
      col <- circus[1:n]
      print(paste0('This palatte have ',length(circus),' colors!'))
    }else if(platte == "paired"){
      col <- paired[1:n]
      print(paste0('This palatte have ',length(paired),' colors!'))
    }else if(platte == "grove"){
      col <- grove[1:n]
      print(paste0('This palatte have ',length(grove),' colors!'))
    }else if(platte == "summerNight"){
      col <- summerNight[1:n]
      print(paste0('This palatte have ',length(summerNight),' colors!'))
    }else if(platte == "zissou"){
      col <- zissou[1:n]
      print(paste0('This palatte have ',length(zissou),' colors!'))
    }else if(platte == "darjeeling"){
      col <- darjeeling[1:n]
      print(paste0('This palatte have ',length(darjeeling),' colors!'))
    }else if(platte == "rushmore"){
      col <- rushmore[1:n]
      print(paste0('This palatte have ',length(rushmore),' colors!'))
    }else if(platte == "captain"){
      col <- captain[1:n]
      print(paste0('This palatte have ',length(captain),' colors!'))
    }else if(platte == "horizon"){
      col <- horizon[1:n] # continues colors
      print(paste0('This palatte have ',length(horizon),' colors!'))
    }else if(platte == "horizonExtra"){
      col <- horizonExtra[1:n]
      print(paste0('This palatte have ',length(horizonExtra),' colors!'))
    }else if(platte == "blueYellow"){
      col <- blueYellow[1:n]
      print(paste0('This palatte have ',length(blueYellow),' colors!'))
    }else if(platte == "sambaNight"){
      col <- sambaNight[1:n]
      print(paste0('This palatte have ',length(sambaNight),' colors!'))
    }else if(platte == "solarExtra"){
      col <- solarExtra[1:n]
      print(paste0('This palatte have ',length(solarExtra),' colors!'))
    }else if(platte == "whitePurple"){
      col <- whitePurple[1:n]
      print(paste0('This palatte have ',length(whitePurple),' colors!'))
    }else if(platte == "whiteBlue"){
      col <- whiteBlue[1:n]
      print(paste0('This palatte have ',length(whiteBlue),' colors!'))
    }else if(platte == "comet"){
      col <- comet[1:n]
      print(paste0('This palatte have ',length(comet),' colors!'))
    }else if(platte == "greenBlue"){
      col <- greenBlue[1:n]
      print(paste0('This palatte have ',length(greenBlue),' colors!'))
    }else if(platte == "beach"){
      col <- beach[1:n]
      print(paste0('This palatte have ',length(beach),' colors!'))
    }else if(platte == "coolwarm"){
      col <- coolwarm[1:n]
      print(paste0('This palatte have ',length(coolwarm),' colors!'))
    }else if(platte == "fireworks"){
      col <- fireworks[1:n]
      print(paste0('This palatte have ',length(fireworks),' colors!'))
    }else if(platte == "greyMagma"){
      col <- greyMagma[1:n]
      print(paste0('This palatte have ',length(greyMagma),' colors!'))
    }else if(platte == "fireworks2"){
      col <- fireworks2[1:n]
      print(paste0('This palatte have ',length(fireworks2),' colors!'))
    }else if(platte == "purpleOrange"){
      col <- purpleOrange[1:n]
      print(paste0('This palatte have ',length(purpleOrange),' colors!'))
    }else{
      print('Please give the correct name!')
    }
    return(col)
  }else{
    dicrete_col <- c('stallion', 'stallion2', 'calm', 'kelly', 'bear' ,
                     'ironMan', 'circus', 'paired', 'grove', 'summerNight' ,
                     'zissou', 'darjeeling', 'rushmore', 'captain')
    
    continues_col <- c('horizon', 'horizonExtra', 'blueYellow', 'sambaNight', 'solarExtra',
                       'whitePurple', 'whiteBlue', 'comet', 'greenBlue', 'beach' ,
                       'coolwarm', 'fireworks', 'greyMagma', 'fireworks2', 'purpleOrange')
    return(c(dicrete_col,continues_col))
  }
}

mg_volcano_custom<- function(diffData = NULL,
                             myMarkers = NULL,
                             order.by = c("logFC"), # c("logFC","P.Value")
                             log2FC.cutoff = 1,
                             pvalue.cutoff = 0.05,
                             adjustP.cutoff = 0.05,
                             topGeneN = 5,
                             col.type = "updown",
                             back.col = 'grey93',
                             pSize = 0.75,
                             aesCol = c('#0099CC','#CC3333'),
                             legend.position = c(0.7,0.9),
                             base_size = 14,
                             tile.col = useMyCol("paired",n = 9),
                             Group.order = NULL,
                             polar = FALSE,
                             expand = c(-1,1),
                             flip = FALSE,
                             ...){
  
  # filter data
  library(dplyr)
  diff.marker <- diffData %>%
    dplyr::filter(abs(logFC) >= log2FC.cutoff & P.Value < pvalue.cutoff)
  
  # assign type,
  diff.marker <- diff.marker %>%
    dplyr::mutate(type = ifelse(logFC >= log2FC.cutoff,"Up","Down")) %>%
    dplyr::mutate(type2 = ifelse(adj.P.Val < adjustP.cutoff,
                                 paste("adjust Pvalue < ",adjustP.cutoff,sep = ''),
                                 paste("adjust Pvalue >= ",adjustP.cutoff,sep = '')))
  head(diff.marker)
  
  # Group orders
  if(!is.null(Group.order)){
    diff.marker$Group <- factor(diff.marker$Group,
                                levels = Group.order)
  }
  
  diff.marker %>% group_by(Group) %>% summarise_at(vars(logFC),list(min,max))
  # get background cols
  back.data<-purrr::map_df(unique(diff.marker$Group),function(x){
    tmp <- diff.marker %>%
      dplyr::filter(Group == x)
    
    new.tmp <- data.frame(Group = x,
                          min = min(tmp$logFC) - 1,
                          max = max(tmp$logFC) + 1)
    return(new.tmp)
  })
  
  # get top gene
  top.marker.tmp <- diff.marker %>%
    dplyr::group_by(Group)
  
  # order
  # if(length(order.by) == 1){
  #   top.marker.max <- top.marker.tmp %>%
  #     dplyr::slice_max(n = topGeneN,order_by = get(order.by))
  #
  #   top.marker.min <- top.marker.tmp %>%
  #     dplyr::group_by(Group) %>%
  #     dplyr::slice_min(n = topGeneN,order_by = get(order.by))
  #
  # }else{
  #   top.marker.max <- top.marker.tmp %>%
  #     dplyr::arrange(dplyr::desc(get(order.by[1])),get(order.by[2])) %>%
  #     dplyr::slice_head(n = topGeneN)
  #
  #   top.marker.min <- top.marker.tmp %>%
  #     dplyr::arrange(dplyr::desc(get(order.by[1])),get(order.by[2])) %>%
  #     dplyr::slice_tail(n = topGeneN)
  # }
  
  top.marker.max <- top.marker.tmp %>%
    dplyr::slice_max(n = topGeneN,order_by = get(order.by))
  
  top.marker.min <- top.marker.tmp %>%
    dplyr::slice_min(n = topGeneN,order_by = get(order.by))
  
  # combine
  top.marker <- rbind(top.marker.max,top.marker.min)
  
  # whether supply own genes
  if(!is.null(myMarkers)){
    top.marker <- diff.marker %>%
      dplyr::filter(gene %in% myMarkers)
  }else{
    top.marker <- top.marker
  }
  
  # ====================================================================
  # plot
  p1 <- ggplot2::ggplot(diff.marker,
                        ggplot2::aes(x = Group,y = logFC)) +
    # add back cols
    ggplot2::geom_col(data = back.data,
                      ggplot2::aes(x = Group,y = min),fill = back.col) +
    ggplot2::geom_col(data = back.data,
                      ggplot2::aes(x = Group,y = max),fill = back.col)
  
  # ap1 <- paste("adjust Pvalue >= ",adjustP.cutoff,sep = '')
  # ap2 <- paste("adjust Pvalue < ",adjustP.cutoff,sep = '')
  
  # color type
  if(col.type == "updown"){
    p2 <- p1 +
      # add point
      ggplot2::geom_jitter(ggplot2::aes(color = type),size = pSize) +
      ggplot2::scale_color_manual(values = c("Down" = aesCol[1],"Up" = aesCol[2]))
  }else if(col.type == "adjustP"){
    p2 <- p1 +
      # add point
      ggplot2::geom_jitter(ggplot2::aes(color = type2),size = pSize) +
      ggplot2::scale_color_manual(values = c(aesCol[2],aesCol[1]))
  }
  
  # theme details
  p3 <- p2 +
    ggplot2::scale_y_continuous(n.breaks = 6) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = legend.position,
                   legend.title = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank()) +
    ggplot2::xlab('Clusters') + ggplot2::ylab('log2FoldChange') +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
  
  # add tile
  p4 <- p3 +
    ggplot2::geom_tile(ggplot2::aes(x = Group,y = 0,fill = Group),
                       color = 'black',
                       height = log2FC.cutoff*2,
                       alpha = 0.3,
                       show.legend = F) +
    ggplot2::scale_fill_manual(values = tile.col) +
    # add gene label
    ggrepel::geom_text_repel(data = top.marker,
                             ggplot2::aes(x = Group,y = logFC,label = gene),
                             max.overlaps = 50,
                             ...)
  
  # whether coord_plolar
  if(polar == TRUE){
    p5 <- p4 +
      geom_textpath(ggplot2::aes(x = Group,y = 0,label = Group)) +
      ggplot2::scale_y_continuous(n.breaks = 6,
                                  expand = ggplot2::expansion(mult = expand)) +
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(legend.position = legend.position,
                     legend.title = ggplot2::element_blank()) +
      ggplot2::coord_polar(clip = 'off',theta = 'x')
  }else{
    # whether flip plot
    if(flip == TRUE){
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_label(ggplot2::aes(x = Group,y = 0,label = Group)) +
        ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank()) +
        ggplot2::coord_flip()
    }else{
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_text(ggplot2::aes(x = Group,y = 0,label = Group)) +
        ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }
  }
  return(p5)
}

get.IOBR.immu.format=function(tcga.t.exp.cibersort){
  tcga.t.exp.cibersort = data.frame(tcga.t.exp.cibersort)
  rownames(tcga.t.exp.cibersort) = tcga.t.exp.cibersort$ID
  tcga.t.exp.cibersort = tcga.t.exp.cibersort[, -1]
  colnames(tcga.t.exp.cibersort) = gsub('(.*)_.*', "\\1", colnames(tcga.t.exp.cibersort))
  return(tcga.t.exp.cibersort)
}

mg_lasso_cox_use=function(dat,time,event,nfolds=3,lambda.min=T,show_text=T,figLabels=c('A','B')){
  library("glmnet") 
  library('survival')
  t.inds=which(!is.na(time)&!is.na(event)&time>0)
  dat=dat[t.inds,]
  time=as.numeric(time[t.inds])
  event=as.numeric(event[t.inds])
  y=Surv(time,event)
  set.seed(123455)
  # set.seed(123456)
  fit1_cv = cv.glmnet(as.matrix(dat), y, family = "cox", nfolds=nfolds
                      ,nlambda=100, alpha=1
  )
  fit<-glmnet(dat, y, family = "cox")
  if(lambda.min){
    lambda=fit1_cv$lambda.min
  }else{
    lambda=fit1_cv$lambda.1se
  }
  coefficients<-coef(fit,s=lambda)
  Active.Index<-which(coefficients[,1]!=0)
  genes=row.names(coefficients)[Active.Index]
  Active.coefficients<-coefficients[Active.Index]  
  g=mg_plot_lasso(fit,fit1_cv,lambda = lambda,show_text=show_text,figLabels=figLabels)
  return(list(Mode1=fit,Model2=fit1_cv,Genes=genes,Coef=Active.coefficients,lambda=lambda,plot=g))
}

createCoxModel_use=function(dat,time,event,isStep=F,direction=c("both", "backward", "forward")[1],check=T){
  cls=colnames(dat)
  dat1=cbind(dat,time,event)
  colnames(dat1)=c(paste0('g',1:ncol(dat)),'time','status')
  dat1=as.data.frame(dat1)
  if(ncol(dat)>nrow(dat)&check){
    print('gene count > sample count')
    return(NULL)
  }
  #nas=apply(dat1, 1, function(x){
  #  return(sum(is.na(x)))
  #})
  #dat1=dat1[which(nas==0),]
  
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(dat1)[1:ncol(dat)],collapse = '+')))
  library(survival)
  cox <- coxph(fmla, data = dat1)
  
  if(isStep){
    tryCatch({
      cox=step(cox,direction =direction,steps = 10000)
    },error = function(e) {
      print(conditionMessage(e))
      return(NULL)
    })
  }
  #score=predict(cox,data=dat1)
  sig.genes=cls[as.numeric(gsub('g','',names(cox$coefficients)))]
  fls=c('RiskScore=')
  for(i in 1:length(sig.genes)){
    if(cox$coefficients[i]>0){
      fls=c(fls,'+',round(cox$coefficients[i],3),'*',sig.genes[i])
    }else{
      fls=c(fls,round(cox$coefficients[i],3),'*',sig.genes[i])
    }
  }
  score=predictRiskScore(cox$coefficients,sig.genes,dat)
  return(list(Cox=cox,Score=score,Genes=sig.genes,Coef=cox$coefficients,fmla=paste0(fls,collapse = '')))
}


ggplotTimeROC_use=function(time,status,score,mks=c(1,3,5)){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    #roc.tm=ROC.DSST$times[which(ROC.DSST$times>0)]
    
    #p.dat=rbind()
    #for(i in which(ROC.DSST$times>0)){
    #los=lowess(ROC.DSST$FP[,i], y=ROC.DSST$TP[,i], f = 1/3, iter = 100)
    #los$x=c(0,los$x,1)
    #los$y=c(0,los$y,1)
    # p.dat=rbind(p.dat,data.frame(los$x, y=los$y,rep(lbs[i],length(los$y)),stringsAsFactors = F))
    #}
    
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction')
    # p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction')
    
    p1=p1+scale_colour_manual(values = c(pal_npg(alpha =0.8)(9)[c(8,4,3,5)]))
    
    p1=p1+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_text(family="Times",face="plain")
                ,axis.title.x=element_text(family="Times",face="plain"),axis.title.y=element_text(family="Times",face="plain")
                ,plot.title=element_blank()
                ,plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches")
                ,legend.position=c(1,0)
                ,legend.justification=c(1,0)
                ,legend.background = element_rect(fill = NA, colour = NA)
                ,legend.title = element_text(family="Times",face="plain")
                ,legend.text = element_text(family="Times",face="plain"))
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}

################### 
mg_Forestplot=function(df_m,outFile,width=6,height=3){
  colnames(df_m)=c('HR','HR.95L','HR.95H','pvalue')
  gene=rownames(df_m)
  # hr=sprintf("%.3f",df_m$"HR")
  # hrLow=sprintf("%.3f",df_m$"HR.95L")
  # hrHigh=sprintf("%.3f",df_m$"HR.95H")
  hr=format(df_m$"HR", digits = 3)
  hrLow=format(df_m$"HR.95L", digits = 3)
  hrHigh=format(df_m$"HR.95H",digits =3)
  Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
  # pVal=ifelse(df_m$pvalue<1e-8, "<1e-8", sprintf("%.3f", df_m$pvalue))
  pVal=format(df_m$pvalue, digits = 3)
  #########
  pdf(file=outFile, width = width, height =height,onefile = FALSE)
  n=nrow(df_m)
  nRow=n+1
  ylim=c(1,nRow)
  layout.show(layout(matrix(c(1,2),nc=2),width=c(2,1.2)))
  xlim = c(0,3)
  par(mar=c(4,1,2,1),mpg=c(2,0.5,0))
  # par(mar=c(3,2,1.5,1.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  

  par(mar=c(4,1,2,1),mpg=c(2,0.5,0))
  # par(mar=c(3,1,1.5,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="black",lwd=2.5,lty=5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
  points(as.numeric(hr), n:1, pch = 20, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
unicox<-function(vars=c("T","N"),time=NULL,event=NULL,data=LUAD_clinical){
  library(survival)
  require(survminer)
  y<-Surv(as.numeric(time),
          as.numeric(event))
  
  ## for
  pvalue<-c()
  HR<-c()
  lower<-c()
  upper<-c()
  varname<-c()
  
  for(i in vars ){
    cox.fit_uni<-coxph(y~data[[i]])
    uni_res <-summary(cox.fit_uni)
    pvalue[[i]]<-uni_res$waldtest[[3]]#pvalue
    HR[[i]]<-uni_res$conf.int[[1]]# HR提取
    lower[[i]] <-uni_res$conf.int[[3]]#lower
    upper[[i]] <-uni_res$conf.int[[4]]# upper
  }

  univar_res<-data.frame(
    # varname<-vars,
    HR<-as.numeric(HR),
    # CI=paste(as.numeric(lower),"-",as.numeric(upper),sep = ""),
    HR.95L<-as.numeric(lower),
    HR.95H<-as.numeric(upper),
    pvalue<-as.numeric(pvalue)
  )
  colnames(univar_res)<-c("HR","HR.95L","HR.95H","pvalue")
  rownames(univar_res)=vars
  univar_res## 
}
multicox<-function(vars=c("T","N","M","age"),time=NULL,event=NULL,data=LUAD_clinical,forest=T){
  library(survival)
  require(survminer)
  y<-Surv(as.numeric(time),
          as.numeric(event))
 
  FM<-as.formula(paste0("y~",paste(vars,collapse = "+")))
  cox.fit_multi <- coxph(FM,data = data)
  munivar_res<-summary(cox.fit_multi)#cox
  pvalue<-munivar_res$coefficients[,"Pr(>|z|)"]#pvalue
  HR<-munivar_res$coefficients[,"exp(coef)"]# HR
  lower<-munivar_res$conf.int[,3]
  upper<-munivar_res$conf.int[,4]
  
  munivar_res<-data.frame(
    # varname<-vars,
    HR<-as.numeric(HR),
    # CI=paste(as.numeric(lower),"-",as.numeric(upper),sep = ""),
    HR.95L<-as.numeric(lower),
    HR.95H<-as.numeric(upper),
    pvalue<-as.numeric(pvalue)
  )
  colnames(munivar_res)<-c("HR","HR.95L","HR.95H","pvalue")
  rownames(munivar_res)=vars
  

  if (forest==T){
    ggforest(cox.fit_multi,data = data)
    ggsave(filename = "mult_cox.pdf")}
  munivar_res## 
  
}
mg_plotDCA=function(status,fmlas,modelNames,data,col=NULL){
  set.seed(123)
  all.mod=list()
  for(i in 1:length(fmlas)){
    fmla <- as.formula(paste0("status~",fmlas[i]))
    model<-rmda::decision_curve(fmla,
                                data=data,
                                bootstraps=500)
    all.mod=c(all.mod,list(model))
  }
  rmda::plot_decision_curve(all.mod,
                            curve.names=modelNames,
                            # col=mg_colors[c(1,10:12,4,5,7:8)],
                            col=col,
                            xlim=c(0,1),legend.position="topright",
                            lwd=1,
                            confidence.intervals=FALSE)
}

plot_cor_point=function(x,y,method=c('pearson','spearman')[1],top_col='#D55E00',right_col='#009E73'
                        ,ylab='y expression',xlab='x expression',title=NULL
                        ,marginal.type=c("histogram", "boxplot", "density", "violin", "densigram")[1]){
  library(ggplot2)
  library(ggpubr)
  library(ggExtra)
  corT=cor.test(x,y,method=method)
  cor=corT$estimate
  pValue=corT$p.value
  df=data.frame(x,y)
  p=ggplot(df, aes(x, y)) + 
    xlab(xlab)+ylab(ylab)+
    geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = method, aes(x =x, y =y))
  p=ggMarginal(p, type = marginal.type, xparams = list(fill = top_col),yparams = list(fill = right_col))
  return(p)
}


# RSF --------------------------------------------------------------
RSF.fs <- function(train) {
  library(randomForestSRC)
  mod.RSF <- rfsrc(Surv(time, status) ~ .,
                   data = train,
                   ntree = 1000,
                   splitrule = "logrank",
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = 123456
  )
  vi <- data.frame(imp=vimp.rfsrc(mod.RSF)$importance)
  vi$imp <- (vi$imp-min(vi$imp))/(max(vi$imp)-min(vi$imp))
  vi$ID <- rownames(vi)
  rid <- rownames(vi)[vi$imp>0.2] ## 0.4 can be adjust
  return(rid)
}

# CoxBoost --------------------------------------------------------------
CoxBoost.fs <- function(train) {
  library(CoxBoost)
  time <- train$time
  status <- train$status
  x <- train %>%
    dplyr::select(-c(1:2)) %>%
    as.matrix()
  # determine penalty parameter
  set.seed(123456)
  optim.CoxBoost <- optimCoxBoostPenalty(
    time = time, status = status, x = x,
    trace = T, start.penalty = 500
  )
  # Fit with obtained penalty parameter and optimal number of boosting
  # steps obtained by cross-validation
  mod.CoxBoost <- CoxBoost(
    time = time, status = status, x = x,
    stepno = optim.CoxBoost$cv.res$optimal.step,
    penalty = optim.CoxBoost$penalty
  )
  rid <- names(coef(mod.CoxBoost)[which(coef(mod.CoxBoost)!=0)])
  return(rid)
}

# stepwiseCox.both --------------------------------------------------------------
stepwiseCox.both.fs <- function(train) {
  fit0 <- coxph(Surv(time, status) ~ ., data = train)
  fit_both <- MASS::stepAIC(fit0, direction = "both")
  rid <- names(coef(fit_both))
  return(rid)
}

# Lasso --------------------------------------------------------------
Lasso.fs <- function(train) {
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      nfolds = 5,
                      family = "cox", # cox
                      grouped = FALSE, 
                      alpha = 1, 
                      # type.measure = "mse"
  )
  coef.min = coef(cv.fit, s = "lambda.min") 
  active.min = which(as.numeric(coef.min)!=0)
  rid <- colnames(x2)[active.min]
  return(rid)
}


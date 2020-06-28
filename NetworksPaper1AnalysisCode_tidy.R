##Load packages
require(lme4)
require(viridis)
require(car)
require(boot)
require(igraph)
require(scatterplot3d)
require(sjPlot)

##Define functions used in generation of dataframe

min_finder<-function(object){
  return(apply(object,1,min))
}

whichmin_finder<-function(object){
  return(apply(object,1,which.min))
}

start_finder<-function(object){
  return(as.numeric(object[,1]>0))
}

max_finder<-function(object){
  return(apply(object,1,max))
}

whichmax_finder<-function(object){
  return(apply(object,1,which.max))
}

ep_start<-function(input,thresh=5){
  return(min(which(input>(thresh-1))))
}

es_finder<-function(object){
  return(apply(object,1,ep_start))
}

##################################################################################

##Set path to parameter files
sourcefolder <- "."
paramsfolder <- "params"

##read in parameter files
params1<-read.csv(file.path(paramsfolder,"model_params.csv"))
params2<-read.csv(file.path(paramsfolder,"model_params2.csv"))

##set path to results files
outputfolder<-file.path("output","results4/")

##################################################################################

##read in results
conc_res<-list()
exp_res<-list()
inf_res<-list()
hosp_res<-list()
he<-list()

c<-1
for(nt in c(1,2,3,4,5,6,7,8,9)){
  for(md in seq(3,77,1)){
    for(r in 1:5){
      tmp_in<-readRDS(paste0(outputfolder,"nets",params1$NetSelect[nt],"mods",params2$ModSelect[md],"rep",r,".RDS"))
      conc_res[[c]]<-tmp_in[[1]]
      exp_res[[c]]<-tmp_in[[2]]
      inf_res[[c]]<-tmp_in[[3]]
      hosp_res[[c]]<-tmp_in[[4]]
      he[[c]]<-tmp_in[[5]]
      c<-c+1
    }
  }
}

##################################################################################

#set up dataframe of summary results
net<-seq(1,9,1)
params<-seq(3,77,1)
reps<-seq(1,5,1)
mods<-seq(1,10,1)

net2<-rep(net,each=length(params)*length(reps)*length(mods))
params2<-rep(rep(params,each=length(reps)*length(mods)),length(net))
reps2<-rep(rep(reps,each=length(mods)),length(params)*length(net))
mods2<-rep(mods,length(params)*length(reps)*length(net))

dat<-data.frame(net2,params2,reps2,mods2)

dat$MinConc<-unlist(lapply(conc_res,min_finder))
dat$WhichMinConc<-unlist(lapply(conc_res,whichmin_finder))
dat$Seeded<-unlist(lapply(exp_res,start_finder))
dat$EpStartExp<-unlist(lapply(exp_res,es_finder))
dat$MaxExp<-unlist(lapply(exp_res,max_finder))
dat$WhichMaxExp<-unlist(lapply(exp_res,whichmax_finder))
dat$EpStartInf<-unlist(lapply(inf_res,es_finder))
dat$MaxInf<-unlist(lapply(inf_res,max_finder))
dat$WhichMaxInf<-unlist(lapply(inf_res,whichmax_finder))
dat$MaxHosp<-unlist(lapply(hosp_res,max_finder))
dat$WhichMaxHosp<-unlist(lapply(hosp_res,whichmax_finder))

##names of the first 4 dataframe columns
names(dat)[1:4]<-c("NetworkID","ParameterID","Rep","Community")

##Add reassurance effect to the dataframe
dat$he<-rep(unlist(he),each=10)

##Set situations where epidemics failed to start to NA
dat$EpStartExp[dat$EpStartExp==Inf]<-NA
dat$EpStartInf[dat$EpStartInf==Inf]<-NA

##Add column to use an identifier of each simulation run
dat$r_eff<-paste0("N",dat$NetworkID,"P",dat$ParameterID,"R",dat$Rep)

##Add column and related information to have colours associated with the reassurance effect 
dat$he_col<-rep(NA,nrow(dat))
he_pal<-plasma(20)
he_qs<-quantile(dat$he,probs=seq(0,1,length.out=21))
for(i in 1:nrow(dat)){
  dat$he_col[i]<-sum(dat$he[i]<=he_qs)
}

#Add columns for mean-centred reassurance effect and epidemic start day
dat$mche<-(dat$he-mean(dat$he))/sd(dat$he)
dat$mcEpStartInf<-dat$EpStartInf-mean(dat$EpStartInf,na.rm=T)

##################################################################################

##need to read in one of the parameter files again due to naming issues
##I need both because of how I coded stuff below
p2<-read.csv(paste0(paramsfolder,"model_params2.csv"))
params2<-read.csv(paste0(paramsfolder,"model_params2.csv"))

##################################################################################

##Create separate dataframes for high, intermediate and low starting levels of concern/adherence
dat_high<-dat[p2[dat$ParameterID,2]==0.5,]
dat_mid<-dat[p2[dat$ParameterID,2]==0.2,]
dat_low<-dat[p2[dat$ParameterID,2]==0.05,]

##################################################################################

##Add columns for severity and logit-transformed severity to the high dataset
dat_high$sev<-(dat_high$MaxInf/max(dat_high$MaxInf))*(1-(dat_high$WhichMaxInf-dat_high$EpStartInf)/max(dat_high$WhichMaxInf-dat_high$EpStartInf,na.rm=TRUE))
dat_high$lt_sev<-car::logit(dat_high$sev)

##################################################################################

#Fit statistical model to logit severity as used in the main text
mod_nohe<-lmer(lt_sev~as.factor(NetworkID)*EpStartInf+as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE,])

##################################################################################

###Code to generate Figure 3

tdat<-dat_high[is.na(dat_high$sev)==FALSE,]
the<-aggregate(tdat$he,by=list(tdat$r_eff),mean)

########################

par(mfrow=c(1,2),mar=c(5,5,2,2))
plot(NULL,xlim=c(0,300),ylim=c(0,1),xlab="Time that epidemic starts",ylab="Predicted severity",las=1,cex.lab=1.8,cex.axis=1.3)

randoms<-ranef(mod_nohe)$r_eff

randomsA<-randoms[rownames(randoms)%in%dat_high$r_eff[p2[dat_high$ParameterID,8]<0.5],]
randomsB<-randoms[rownames(randoms)%in%dat_high$r_eff[p2[dat_high$ParameterID,8]>0.49],]

theA<-the[rownames(randoms)%in%dat_high$r_eff[p2[dat_high$ParameterID,8]<0.5],]
theB<-the[rownames(randoms)%in%dat_high$r_eff[p2[dat_high$ParameterID,8]>0.49],]

x<-seq(1,300,0.01)

for(i in 1:nrow(randomsA)){
  x<-seq(min(dat_high$EpStartInf[dat_high$r_eff==rownames(randomsA)[i]],na.rm=T),max(dat_high$EpStartInf[dat_high$r_eff==rownames(randomsA)[i]],na.rm=T),length.out=200)
  y_n<-(fixef(mod_nohe)[1]+fixef(mod_nohe)[11]+randomsA[i,1]+(fixef(mod_nohe)[10]+fixef(mod_nohe)[27]+randomsA[i,2])*x)
  lines(x,boot::inv.logit(y_n),col=he_pal[dat$he_col[dat$he==theA[i,2]]])
}


plot(NULL,xlim=c(0,300),ylim=c(0,1),xlab="Time that epidemic starts",ylab="Predicted severity",las=1,cex.lab=1.8,cex.axis=1.3)

for(i in 1:nrow(randomsB)){
  x<-seq(min(dat_high$EpStartInf[dat_high$r_eff==rownames(randomsB)[i]],na.rm=T),max(dat_high$EpStartInf[dat_high$r_eff==rownames(randomsB)[i]],na.rm=T),length.out=200)
  y_n<-(fixef(mod_nohe)[1]+fixef(mod_nohe)[14]+randomsB[i,1]+(fixef(mod_nohe)[10]+fixef(mod_nohe)[30]+randomsB[i,2])*x)
  lines(x,boot::inv.logit(y_n),col=he_pal[dat$he_col[dat$he==theB[i,2]]])
}

##################################################################################


#5x5 plots for supplement

socs<-c(0,0.1,0.2,0.5,1)
infs<-c(0,0.2,0.4,0.8,1.6)

#Code for Figure S1
par(mfrow=c(5,5),mar=c(2,2,2,2))

for(i in socs){
  for(j in infs){
    plot(dat_high$MaxInf[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]~dat_high$WhichMaxInf[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j],col=adjustcolor(he_pal[dat_high$he_col[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]],0.4),pch=16,xlim=c(0,300),ylim=c(0,100),las=1)
  }
}

#Code for Figure S2
par(mfrow=c(5,5),mar=c(0,0,0,0))

for(i in socs){
  for(j in infs){
    plot(dat_high$sev[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]~dat_high$EpidemicStartTime1[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j],col=adjustcolor(he_pal[dat_high$he_col[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]],0.4),pch=16,xlim=c(0,300),ylim=c(0,1),yaxt="n",xaxt="n")
  }
}

#Code for Figure S5
par(mfrow=c(5,5),mar=c(0,0,0,0))
for(i in socs){
  for(j in infs){
    plot(dat_high$sev[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]~dat_high$MinConc[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j],col=adjustcolor(he_pal[dat_high$he_col[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]],0.4),pch=16,xlim=c(0,1),ylim=c(0,1),yaxt="n",xaxt="n")
  }
}

##################################################################################

##Code for Figure 2

socs2<-c(0,0.2,1)
infs2<-c(0,0.4,1.6)

par(mfrow=c(3,3),mar=c(5,5,2,2))
for(i in socs2){
  for(j in infs2){
    ms<-c(5,5,2,2)
    par(mar=ms)
    if(j==0){
      ms[c(2,4)]<-c(5,0)
      par(mar=ms)
    }
    if(j==0.4){
      ms[c(2,4)]<-c(4,1)
      par(mar=ms)
    }
    if(j==1.6){
      ms[c(2,4)]<-c(3,2)
      par(mar=ms)
    }
    if(i==0){
      ms[c(1,3)]<-c(3,2)
      par(mar=ms)
    }
    if(i==0.2){
      ms[c(1,3)]<-c(4,1)
      par(mar=ms)
    }
    if(i==1){
      ms[c(1,3)]<-c(5,0)
      par(mar=ms)
    }
    plot(dat_high$sev[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]~dat_high$EpStartInf[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j],col=adjustcolor(he_pal[dat_high$he_col[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]],0.4),pch=16,xlim=c(0,300),ylim=c(0,1),ylab=ifelse(j==0,"Epidemic Severity",""),xlab=ifelse(i==1,"Time that epidemic starts",""),cex.lab=1.8,cex.axis=1.1,las=1)
  }
}

##################################################################################

##Code for Figure 4

socs2<-c(0,0.2,1)
infs2<-c(0,0.4,1.6)

par(mfrow=c(3,3),mar=c(5,5,2,2))
for(i in socs2){
  for(j in infs2){
    ms<-c(5,5,2,2)
    par(mar=ms)
    if(j==0){
      ms[c(2,4)]<-c(5,0)
      par(mar=ms)
    }
    if(j==0.4){
      ms[c(2,4)]<-c(4,1)
      par(mar=ms)
    }
    if(j==1.6){
      ms[c(2,4)]<-c(3,2)
      par(mar=ms)
    }
    if(i==0){
      ms[c(1,3)]<-c(3,2)
      par(mar=ms)
    }
    if(i==0.2){
      ms[c(1,3)]<-c(4,1)
      par(mar=ms)
    }
    if(i==1){
      ms[c(1,3)]<-c(5,0)
      par(mar=ms)
    }
    plot(dat_high$sev[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]~dat_high$MinConc[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j],col=adjustcolor(he_pal[dat_high$he_col[p2[dat_high$ParameterID,8]==i&p2[dat_high$ParameterID,10]==j]],0.4),pch=16,xlim=c(0,1),ylim=c(0,1),ylab=ifelse(j==0,"Epidemic Severity",""),xlab=ifelse(i==1,"Minimum concern",""),cex.lab=1.8,cex.axis=1.1,las=1)
  }
}

##################################################################################

##Code for Figure 1c and d

par(mfrow=c(2,1),mar=c(4.5,5,1,2))
col10<-c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")

plot(NULL,xlim=c(0,150),ylim=c(0,60),xlab="Days",ylab="Number infected",cex.lab=1.6,cex.axis=1.3,las=1)
for(i in 1:nrow(inf_res[[800]])){
  lines(inf_res[[3305]][i,],col=col10[i],lwd=3)
}

plot(NULL,xlim=c(0,150),ylim=c(0,1),xlab="Days",ylab="Proportion concerned",cex.lab=1.6,cex.axis=1.3,las=1)
for(i in 1:nrow(inf_res[[800]])){
  lines(conc_res[[3305]][i,],col=col10[i],lwd=3)
}

mtext(paste0("Reassurance effect is ",round(he[[3305]],3)),side=1,adj=0.95,line=-1.5,cex=1.2)
mtext(paste0("Weak Social Construction"),side=1,adj=0.95,line=-3,cex=1)
mtext(paste0("Weak Awareness"),side=1,adj=0.95,line=-4.5,cex=1)

##################################################################################

##Code for Figure 1a and b

#read in network and filter so that it is only young adults of one predisposition
n_plot<-readRDS(paste0(networksfolder,"88net_and_parents.RDS"))
n_plotB<-n_plot[[1]][c(1:240,481:1110,1741:1870),c(1:240,481:1110,1741:1870)]
n_plotBy<-n_plot[[1]][481:1110,481:1110]

n_plot2<-graph.adjacency(n_plotBy,mode="undirected")

#Plot the network
par(mfrow=c(1,1),mar=c(0,0,0,0))                
plot(n_plot2,vertex.label=NA,vertex.size=6,vertex.color=rep(col10,100),edge.color="grey")


###Plot multiplex structure for dark blue community

#Communication layer
n_plotBy1<-n_plot[[1]][seq(481,1101,10),seq(481,1101,10)]
n_plot3<-graph.adjacency(n_plotBy1,mode="undirected")
lo<-layout.fruchterman.reingold(n_plot3)
plot(n_plot3,layout=lo,vertex.label=NA,vertex.size=6,vertex.color=col10[1],edge.color="grey")

#Infection layer
m_plot<-readRDS(paste0(networksfolder,"8net_and_parents.RDS"))
m_plotBy1<-m_plot[[1]][seq(481,1101,10),seq(481,1101,10)]
m_plot3<-graph.adjacency(m_plotBy1,mode="undirected")
plot(m_plot3,layout=lo,vertex.label=NA,vertex.size=6,vertex.color=col10[1],edge.color="grey")

#3D scatterplot
Dheight<-1
Iheight<-9

cD<-cbind(lo,rep(Dheight,63))
cI<-cbind(lo,rep(Iheight,63))

C<-scatterplot3d(cD,xlim=c(-6,6),ylim=c(-6,6),zlim=c(0,12),color=col10[1],pch=16,box=F,grid=F,cex.symbols=1.5,angle=70,axis=F,scale.y=0.5)
C$points3d(cI[,1],cI[,2],cI[,3],col=col10[1],pch=18,cex=1.5)

theta<-seq(0,2*pi,length=1000)
alpha<-pi/10
ell.top.x<- -0.25+7*cos(theta)*cos(alpha)-7*sin(theta)*sin(alpha)
ell.top.y<- +0.25+7*cos(theta)*sin(alpha)+7*sin(theta)*cos(alpha)
C$points3d(ell.top.x,ell.top.y,rep(Dheight,length(ell.top.x)),type="l",col=adjustcolor(col10[1],alpha=1))
C$points3d(ell.top.x,ell.top.y,rep(Iheight,length(ell.top.x)),type="l",col=adjustcolor(col10[1],alpha=1))
polygon(C$xyz.convert(ell.top.x,ell.top.y,rep(Dheight,length(ell.top.x)))$x,C$xyz.convert(ell.top.x,ell.top.y,rep(Dheight,length(ell.top.x)))$y,col=adjustcolor(col10[1],0.08),border=NA)
polygon(C$xyz.convert(ell.top.x,ell.top.y,rep(Iheight,length(ell.top.x)))$x,C$xyz.convert(ell.top.x,ell.top.y,rep(Iheight,length(ell.top.x)))$y,col=adjustcolor(col10[1],0.08),border=NA)

for(i in 1:63){
  C$points3d(c(cI[i,1],cD[i,1]),c(cI[i,2],cD[i,2]),c(cI[i,3],cD[i,3]),type="l",col="grey",lty=2)
}

for(i in 1:(ncol(n_plotBy1)-1)){
  for(j in (i+1):ncol(n_plotBy1)){
    if(n_plotBy1[i,j]>0){
      C$points3d(c(cI[i,1],cI[j,1]),c(cI[i,2],cI[j,2]),c(cI[i,3],cI[j,3]),type="l",col="dark grey",lwd=1,lty=1)
    }
  }
}

for(i in 1:(ncol(m_plotBy1)-1)){
  for(j in (i+1):ncol(m_plotBy1)){
    if(m_plotBy1[i,j]>0){
      C$points3d(c(cD[i,1],cD[j,1]),c(cD[i,2],cD[j,2]),c(cD[i,3],cD[j,3]),type="l",col="dark grey",lwd=1,lty=1)
    }
  }
}

C$points3d(cI[,1],cI[,2],cI[,3],col=col10[1],pch=18,cex=1.5)
C$points3d(cD[,1],cD[,2],cD[,3],col=col10[1],pch=16,cex=1.5)

##################################################################################

##Separate models for each multiplex network

mod_nohe1<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==1,])

mod_nohe2<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==2,])

mod_nohe3<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==3,])

mod_nohe4<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==4,])

mod_nohe5<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==5,])

mod_nohe6<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==6,])

mod_nohe7<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==7,])

mod_nohe8<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==8,])

mod_nohe9<-lmer(lt_sev~as.factor(p2[ParameterID,8])*EpStartInf*as.factor(p2[ParameterID,10])+(EpStartInf|r_eff),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==9,])

##################################################################################

##Code for Figures S6 and S7

randoms1<-ranef(mod_nohe1)$r_eff
the1<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==1],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==1]),mean)

randoms2<-ranef(mod_nohe2)$r_eff
the2<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==2],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==2]),mean)

randoms3<-ranef(mod_nohe3)$r_eff
the3<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==3],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==3]),mean)

randoms4<-ranef(mod_nohe4)$r_eff
the4<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==4],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==4]),mean)

randoms5<-ranef(mod_nohe5)$r_eff
the5<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==5],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==5]),mean)

randoms6<-ranef(mod_nohe6)$r_eff
the6<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==6],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==6]),mean)

randoms7<-ranef(mod_nohe7)$r_eff
the7<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==7],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==7]),mean)

randoms8<-ranef(mod_nohe8)$r_eff
the8<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==8],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==8]),mean)

randoms9<-ranef(mod_nohe9)$r_eff
the9<-aggregate(dat_high$he[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==9],by=list(dat_high$r_eff[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==9]),mean)

par(mar=c(5,5,2,2))
par(mfrow=c(3,3))

for(i in 1:9){
  xt<-get(paste0("the",i))
  yt<-get(paste0("randoms",i))
  plot(yt[,2]~abs(xt$x),pch=16,ylab="Random Slope Estimates",xlab="Strength of Reassurance Effect",las=1,cex.axis=0.8,cex.lab=1.2)
  mtext(text=paste("NetworkID",i),side=3,adj=0.5,line=0.5)
  mtext(text=round(cor.test(abs(xt$x),yt[,2])$estimate,3),side=1,adj=1,line=-1)
}

for(i in 1:9){
  xt<-get(paste0("the",i))
  yt<-get(paste0("randoms",i))
  plot(yt[,1]~abs(xt$x),pch=16,ylab="Random Intercept Estimates",xlab="Strength of Reassurance Effect",las=1,cex.axis=0.8,cex.lab=1.2)
  mtext(text=paste("NetworkID",i),side=3,adj=0.5,line=0.5)
  mtext(text=round(cor.test(abs(xt$x),yt[,1])$estimate,3),side=1,adj=1,line=-1)
}

##################################################################################

##Code for Figures S3 and S4

peak_vars<-aggregate(dat$WhichMaxInf,by=list(dat$r_eff),var,na.rm=T)

ses<-aggregate(dat$ParameterID,by=list(dat$r_eff),unique)[,2]
soc_effs<-rep(NA,length(ses))
dis_effs<-rep(NA,length(ses))

for(i in 1:length(ses)){
  soc_effs[i]<-params2[ses[i],8]
  dis_effs[i]<-params2[ses[i],10]
}

hea_effs<-aggregate(dat$he,by=list(dat$r_eff),unique)[,2]


par(mar=c(5,5,2,2))
par(mfrow=c(1,1))

#Figure S3
boxplot(log(peak_vars[,2])~soc_effs,range=0,lty=1,xlab="Social Construction Effect Size",ylab="Log Variation in in Outbreak Peak Timings",cex.lab=1.5,cex.axis=1.5,las=1)

#Figure S4
plot(log(peak_vars[,2])~abs(hea_effs),las=1,cex.axis=1.5,cex.lab=1.5,ylab="Log Variation in Outbreak Peak Timings",xlab="Strength of Reasurrance Effect",pch=16,col=adjustcolor("black",0.5))

##################################################################################
##################################################################################
##################################################################################

##In this section of code we produce the Supplementary Tables
#We first make versions of the models with more informative variable names, refit them and produce table using sjPlot

##First for dat_high

names(dat_high)[22]<-"LogitSeverity"
names(dat_high)[11]<-"EpidemicStartTime1"
names(dat_high)[17]<-"SimulationNumber"

dat_high$NetworkID<-as.factor(dat_high$NetworkID)
dat_high$EpidemicStartTime<-dat_high$EpidemicStartTime1-60
dat_high$SocialConstruction<-rep(NA,nrow(dat_high))
dat_high$Awareness<-rep(NA,nrow(dat_high))

for(i in 1:nrow(dat_high)){
  dat_high$SocialConstruction[i]<-p2[dat_high$ParameterID[i],8]
  dat_high$Awareness[i]<-p2[dat_high$ParameterID[i],10]
}

dat_high$SocialConstruction<-as.factor(dat_high$SocialConstruction)
dat_high$Awareness<-as.factor(dat_high$Awareness)

Mod_nohe<-lmer(LogitSeverity~NetworkID*EpidemicStartTime+SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE,])

tab_model(Mod_nohe)

##################################################################################

##Now for dat_mid

dat_mid$sev<-(dat_mid$MaxInf/max(dat_mid$MaxInf))*(1-(dat_mid$WhichMaxInf-dat_mid$EpStartInf)/max(dat_mid$WhichMaxInf-dat_mid$EpStartInf,na.rm=TRUE))
dat_mid$lt_sev<-car::logit(dat_mid$sev)

names(dat_mid)[22]<-"LogitSeverity"
names(dat_mid)[11]<-"EpidemicStartTime1"
names(dat_mid)[17]<-"SimulationNumber"

dat_mid$NetworkID<-as.factor(dat_mid$NetworkID)
dat_mid$EpidemicStartTime<-dat_mid$EpidemicStartTime1-60
dat_mid$SocialConstruction<-rep(NA,nrow(dat_mid))
dat_mid$Awareness<-rep(NA,nrow(dat_mid))

for(i in 1:nrow(dat_mid)){
  dat_mid$SocialConstruction[i]<-p2[dat_mid$ParameterID[i],8]
  dat_mid$Awareness[i]<-p2[dat_mid$ParameterID[i],10]
}

dat_mid$SocialConstruction<-as.factor(dat_mid$SocialConstruction)
dat_mid$Awareness<-as.factor(dat_mid$Awareness)

Mod_noheM<-lmer(LogitSeverity~NetworkID*EpidemicStartTime+SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_mid[is.na(dat_mid$sev)==FALSE,])

tab_model(Mod_noheM)

##################################################################################

##Now for separate networks and dat_high

Mod_nohe1<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==1,])

Mod_nohe2<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==2,])

Mod_nohe3<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==3,])

Mod_nohe4<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==4,])

Mod_nohe5<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==5,])

Mod_nohe6<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==6,])

Mod_nohe7<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==7,])

Mod_nohe8<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==8,])

Mod_nohe9<-lmer(LogitSeverity~SocialConstruction*EpidemicStartTime*Awareness+(EpidemicStartTime|SimulationNumber),data=dat_high[is.na(dat_high$sev)==FALSE&dat_high$NetworkID==9,])

##For these tables we aggregate by homophily
tab_model(Mod_nohe1,Mod_nohe2,Mod_nohe3,collapse.ci=TRUE)
tab_model(Mod_nohe4,Mod_nohe5,Mod_nohe6,collapse.ci=TRUE)
tab_model(Mod_nohe7,Mod_nohe8,Mod_nohe9,collapse.ci=TRUE)


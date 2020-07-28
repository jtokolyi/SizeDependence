library(ggplot2); library(lme4); library(nlme); library(EnvStats)

x <- read.csv("~/hidra/2019/SizeDependence/M28seasonal_full_dataset2.csv",sep="\t")
x$Time <- factor(as.character(x$Time), levels=c("2018spring","2018autumn","2019spring","2019autumn"))
x$Season <- factor(as.character(x$Season), levels=c("Spring","Autumn"))
x$Year <- factor(as.character(x$Year), levels=c("2018","2019"))

x$ReprMode <- ifelse(x$Sex=="ASEX","asex","sex")
x$ReprModeBinomial <- ifelse(x$Sex=="ASEX",0,1)
x$Strain <- as.character(x$Strain)
x$Survival <- ifelse(x$Survival==1, "Yes", "No"); x$Survival <- as.factor(x$Survival)
levels(x$Sex) <- c("Asexual","Female","Male")

x$BodySizeScaled <- scale(x$Body_size)
x$PolypAgeScaled <- scale(x$PolypAge)

### Body size vs. reproductive mode
m1 <- glmer(ReprModeBinomial~scale(Body_size)+scale(PolypAge)+Time+(1|Strain),family="binomial",data=x)
summary(m1)

### Body size vs. fecundity
m <- x[x$Sex=="Male",]
f <- x[x$Sex=="Female",]
m2 <- glmer(MaxTestis~scale(Body_size)+scale(PolypAge)+Time+(1|Strain),family="poisson",data=m)
summary(m2)
f$obs <- 1:nrow(f)
m3 <- glmer(SumEggs~scale(Body_size)+scale(PolypAge)+Time+(1|Strain)+(1|obs),family="poisson",data=f)
summary(m3)

### Body size vs sex start
m4 <- lme(log(SexStart)~scale(Body_size)+scale(PolypAge)+Time, random=~1|Strain,data=m[which(!is.na(m$Body_size)),])
summary(m4)
m5 <- lme(log(SexStart)~scale(Body_size)+scale(PolypAge)+Time, random=~1|Strain,data=f[which(!is.na(f$Body_size)),])
summary(m5)

library(GLMMadaptive)
### Body size vs. survival
m6 <- glmer(Survival~scale(Body_size)+scale(PolypAge)+Time+(1|Strain),family="binomial",data=x[which(!is.na(x$Survival)),])
summary(m6)
m7 <- glmer(Survival~scale(Body_size)+scale(PolypAge)+Time+(1|Strain),family="binomial",data=m[which(!is.na(m$Survival)),])
summary(m7)
m8 <- glmer(Survival~scale(Body_size)+scale(PolypAge)+Time+(1|Strain),family="binomial",data=f[which(!is.na(f$Survival)),])
summary(m8)


tiff("population2.tif",compression="lzw", height=8, width=7, res=300, units="in")
layout(matrix(1:6,ncol=2,byrow=T))
par(mar=c(2.5,4,0,1),oma=c(3,0,2,0),xpd=NA,cex.lab=1.2)

plot(x=m$Body_size, y=log(m$SexStart),pch=19, col=rgb(35/255,139/255,34/255,0.5),cex=1.5,xlab="",ylab="Time to gonadogenensis",ylim=c(2.2,5),xlim=c(0,4.5),axes=F)
axis(1)
axis(2,at=seq(2.5,5,0.5),labels=round(exp(seq(2.5,5,0.5))))
box()
mtext(side=3,text="Males",font=2)
text("a",y=5,x=-0.8,xpd=NA,cex=2,font=2)
m.mod <- lme(log(SexStart)~Body_size, random=~1|Strain,data=m[which(!is.na(m$Body_size)),])
lines(y=predict(m.mod, newdata=data.frame(Body_size=seq(0,4.4,0.1)), level=0, type="response"), x=seq(0,4.4,0.1),lwd=3,col=grey(0.3))
plot(x=f$Body_size, y=log(f$SexStart),pch=19, col=rgb(1,0,0,0.5),cex=1.5,xlab="",ylab="Time to gonadogenesis",ylim=c(2.2,5),xlim=c(0,4.5),axes=F)
mtext(side=3,text="Females",font=2)
axis(1)
axis(2,at=seq(2.5,5,0.5),labels=round(exp(seq(2.5,5,0.5))))
box()
text("b",y=5,x=-0.8,xpd=NA,cex=2,font=2)
f.mod <- lme(log(SexStart)~Body_size, random=~1|Strain,data=f[which(!is.na(f$Body_size)),])
lines(y=predict(f.mod, newdata=data.frame(Body_size=seq(0.1,3.8,0.1)), level=0, type="response"), x=seq(0.1,3.8,0.1),lwd=3,col=grey(0.3))

plot(x=m$Body_size, y=(m$MaxTestis),pch=19, col=rgb(34/255,139/255,34/255,0.5),cex=1.5,xlab="",ylab="No. testes",xlim=c(0,4.5))
text("c",y=25,x=-0.8,xpd=NA,cex=2,font=2)
m.mod <- glmer(MaxTestis~Body_size+(1|Strain),family="poisson",data=m)
lines(y=predict(m.mod, newdata=data.frame(Body_size=seq(0,4.4,0.1)), re.form=NA, type="response"), x=seq(0,4.4,0.1),lwd=3,col=grey(0.3))
plot(x=f$Body_size, y=(f$SumEggs),pch=19, col=rgb(1,0,0,0.5),cex=1.5,xlab="",ylab="No. eggs",xlim=c(0,4.5))
text("d",y=67.5,x=-0.8,xpd=NA,cex=2,font=2)
f.mod <- glmer(SumEggs~Body_size+(1|Strain)+(1|obs),family="poisson",data=f)
lines(y=predict(f.mod, newdata=data.frame(Body_size=seq(0.1,3.8,0.1)), re.form=NA, type="response"), x=seq(0.1,3.8,0.1),lwd=3,col=grey(0.3))

plot(x=m$Body_size, y=as.numeric(m$Survival)-1,col=rgb(34/255,139/255,34/255,0.5),cex=1.5,pch=19,
     xlab=expression(paste("Body size (mm"^2,")",sep="")),xlim=c(0,4.5),
     ylab="Survived",axes=F)
axis(1)
axis(2, at=c(0,1), labels=c("No","Yes"))
box()
text("e",y=1,x=-0.8,xpd=NA,cex=2,font=2)
m.mod <- glmer(Survival~Body_size+(1|Strain),family="binomial",data=m[which(!is.na(m$Body_size) & !is.na(m$Survival)),])
lines(y=predict(m.mod, newdata=data.frame(Body_size=seq(0,4.4,0.1)), re.form=NA, type="response"), x=seq(0,4.4,0.1), lwd=3, col=grey(0.3))

plot(x=f$Body_size, y=as.numeric(f$Survival)-1,col=rgb(1,0,0,0.5),cex=1.5,xlab=expression(paste("Body size (mm"^2,")",sep="")),
                                                                                          xlim=c(0,4.5),pch=19,
     ylab="Survived",axes=F)
axis(1)
axis(2, at=c(0,1), labels=c("No","Yes"))
box()
text("f",y=1,x=-0.8,xpd=NA,cex=2,font=2)
f.mod <- glmer(Survival~Body_size+(1|Strain),family="binomial",data=f[which(!is.na(f$Body_size) & !is.na(f$Survival)),])
lines(y=predict(f.mod, newdata=data.frame(Body_size=seq(0.1,3.8,0.1)), re.form=NA, type="response"), x=seq(0.1,3.8,0.1), lwd=3, col=grey(0.3))

dev.off()

x$RelGonads <- NA
x$RelGonads[which(x$Sex=="Male")] <- log(x$MaxTestis[which(x$Sex=="Male")]/x$Body_size[which(x$Sex=="Male")]+1)
x$RelGonads[which(x$Sex=="Female")] <- log(x$SumEggs[which(x$Sex=="Female")]/x$Body_size[which(x$Sex=="Female")]+1)
m <- x[which(x$Sex=="Male"),]
f <- x[which(x$Sex=="Female"),]

m.model <- lme(RelGonads~scale(Body_size)+scale(PolypAge)+Time, random=~1|Strain,data=m[which(!is.na(m$Body_size)),])
summary(m.model)
f.model <- lme(RelGonads~scale(Body_size)+scale(PolypAge)+Time, random=~1|Strain,data=f[which(!is.na(f$Body_size)),])
summary(f.model)

### Sex interaction
x$Gonads <- NA
x$Gonads[which(x$Sex=="Male")] <- scale(x$MaxTestis[which(x$Sex=="Male")])
x$Gonads[which(x$Sex=="Female")] <- scale(x$SumEggs[which(x$Sex=="Female")])
x$RelGonads[which(x$Sex=="Male")] <- scale(x$RelGonads[which(x$Sex=="Male")])
x$RelGonads[which(x$Sex=="Female")] <- scale(x$RelGonads[which(x$Sex=="Female")])

m1 <- lme(log(SexStart)~Sex*scale(Body_size)+scale(PolypAge)+Time, random=~1|Strain,data=x[-which(is.na(x$Sex)|is.na(x$Body_size)|is.na(x$SexStart)),])
summary(m1)
m2 <- lme(Gonads~Sex*scale(Body_size)+scale(PolypAge)+Time, random=~1|Strain,data=x[-which(is.na(x$Sex)|is.na(x$Body_size)|is.na(x$Gonads)),])
summary(m2)
m3 <- lme(RelGonads~Sex*scale(Body_size)+scale(PolypAge)+Time, random=~1|Strain,data=x[-which(is.na(x$Sex)|is.na(x$Body_size)|is.na(x$RelGonads)),])
summary(m3)
m4 <- glmer(Survival~as.character(Sex)*scale(Body_size)+scale(PolypAge)+Time+(1|Strain),family="binomial",data=x[-which((is.na(x$Sex)|is.na(x$Body_size)|is.na(x$Survival)|x$Sex=="Asexual")),])
summary(m4)
## convergence error. Refitting with GLMMadaptive
m4 <- mixed_model(Survival~as.character(Sex)*scale(Body_size)+scale(PolypAge)+Time, random=~1|Strain,family="binomial",data=x[-which((is.na(x$Sex)|is.na(x$Body_size)|is.na(x$Survival)|x$Sex=="Asexual")),])
summary(m4) ## same result!

setwd("~/hidra/2019/SizeDependence/experiment/")

library(readxl); library(RColorBrewer); library(lme4); library(nlme); library(multcomp)

sm_data <- as.data.frame(read_excel("sm_data.xlsx"))

sm_data$Start_date<-as.Date(paste("2019_",sm_data$Start_date,collapse=NULL),format="%Y_%m_%d")
sm_data$Cooling_date<-as.Date(paste("2019_",sm_data$Cooling_date,collapse=NULL),format="%Y_%m_%d")
sm_data$Start_time<-as.numeric(sm_data$Start_date - sm_data$Cooling_date)
sm_data$PairID <- paste(sm_data$Plate_number,
                        ifelse(sm_data$Exp_group%in%c("reduced","enlarged"), "RE", "CC"), sep="_")

sm_data <- sm_data[-which(sm_data$PairID%in%unique(sm_data$PairID[which(is.na(sm_data$post_area))])),]

## Exp group modified: combining the two control groups
sm_data$Exp_group2<-sm_data$Exp_group
sm_data$Exp_group2[sm_data$Exp_group2=="ctrl1"]<-"ctrl"
sm_data$Exp_group2[sm_data$Exp_group2=="ctrl2"]<-"ctrl"
sm_data$Exp_group2 <- factor(sm_data$Exp_group2,levels=c("reduced","ctrl","enlarged"))

## Defining data to be read, excluding N/A
sm_data_good<-sm_data[-which(sm_data$Fate=="Fail"),]
sm_data_good$Fate[-which(sm_data_good$Fate=="Regenerated")]<-"Not regenerated"

males<-sm_data_good[sm_data_good$Sex=="C2/7",]
females<-sm_data_good[sm_data_good$Sex=="X11/14",]

males <- males[which(!is.na(males$pre_area)),]

f.cols<-brewer.pal(3,"Reds")
m.cols <- brewer.pal(3,"Greens")
names <- c("Reduced", "Control", "Enlarged")

mnull <- lme(pre_area~1, random=~1|PairID, data=males, method="ML")
fnull <- lme(pre_area~1, random=~1|PairID, data=females, method="ML")

m0<-lme(pre_area~Exp_group2, random=~1|PairID, data=males,method="ML")
anova(m0, mnull)
ps <- summary(glht(m0,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
f0<-lme(pre_area~Exp_group2, random=~1|PairID, data=females,method="ML")
anova(f0, fnull)
ps <- summary(glht(f0,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))

tiff("size_manipulation_results.tif", height=8, width=7, units="in",res=300, compression="lzw")
layout(matrix(1:8, ncol=2, byrow=T))
par(mar=c(1,4,0,1),oma=c(2,0,2,0))

mnull<-lme(post_area~1, random=~1|PairID, data=males, method="ML")
fnull<-lme(post_area~1, random=~1|PairID, data=females, method="ML")

m1<-lme(post_area~Exp_group2, random=~1|PairID, data=males, method="ML")
anova(m1, mnull)
ps <- summary(glht(m1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
range.dep <- range(males$post_area, na.rm=T)
xx <- boxplot(males$post_area~males$Exp_group2, col=m.cols, names=NA, outline=F,
              ylim=c(range.dep[1] - 0.1*(diff(range.dep)), 1.2 * range.dep[2]),
              xlab="",ylab="Post-treatment size")
points(y=males$post_area, x=jitter(as.numeric(as.factor(males$Exp_group2))),pch=19, col=rgb(0,0,0,0.5))
if(ps[[10]]$pvalues[1] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=1.05, x1=1.95,length=0)
    text(ifelse(ps[[10]]$pvalues[1]>0.01, "*", ifelse(ps[[10]]$pvalues[1]>0.01, "**", "***")),x=1.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[2] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=2.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[2]>0.01, "*", ifelse(ps[[10]]$pvalues[2]>0.01, "**", "***")),x=2.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[3] < 0.05){
    arrows(y0=1.15*range.dep[2], y1=1.15*range.dep[2], x0=1.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[3]>0.01, "*", ifelse(ps[[10]]$pvalues[3]>0.01, "**", "***")),x=2, y=1.175*range.dep[2])
}
mtext(side=3,text="Male strain",font=2)
text("a", x=-0.1, y=4.1,xpd=NA,cex=2, font=2)

f1<-lme(post_area~Exp_group2, random=~1|PairID, data=females, method="ML")
anova(f1, fnull)
ps <- summary(glht(f1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
range.dep <- range(females$post_area, na.rm=T)
xx <- boxplot(females$post_area~females$Exp_group2, col=f.cols, names=NA, outline=F,
              ylim=c(range.dep[1] - 0.1*(diff(range.dep)), 1.2 * range.dep[2]),
              xlab="",ylab="Post-treatment size")
points(y=females$post_area, x=jitter(as.numeric(as.factor(females$Exp_group2))),pch=19, col=rgb(0,0,0,0.5))
if(ps[[10]]$pvalues[1] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=1.05, x1=1.95,length=0)
    text(ifelse(ps[[10]]$pvalues[1]>0.01, "*", ifelse(ps[[10]]$pvalues[1]>0.01, "**", "***")),x=1.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[2] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=2.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[2]>0.01, "*", ifelse(ps[[10]]$pvalues[2]>0.01, "**", "***")),x=2.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[3] < 0.05){
    arrows(y0=1.15*range.dep[2], y1=1.15*range.dep[2], x0=1.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[3]>0.01, "*", ifelse(ps[[10]]$pvalues[3]>0.01, "**", "***")),x=2, y=1.175*range.dep[2])
}
mtext(side=3,text="Female strain",font=2)
text("b", x=-0.1, y=4.1,xpd=NA,cex=2, font=2)

mnull<-lme(Start_time~1, random=~1|PairID, data=males, method="ML")
fnull<-lme(Start_time~1, random=~1|PairID, data=females, method="ML")

m1<-lme(Start_time~Exp_group2, random=~1|PairID, data=males, method="ML")
anova(m1, mnull)
ps <- summary(glht(m1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
range.dep <- c(4,36)#range(males$Start_time, na.rm=T)
xx <- boxplot(males$Start_time~males$Exp_group2, col=m.cols, names=NA, outline=F,
              ylim=c(range.dep[1] - 0.1*(diff(range.dep)), 1.2 * range.dep[2]),
              xlab="",ylab="Time to gonadogenesis (days)")
points(y=males$Start_time, x=jitter(as.numeric(as.factor(males$Exp_group2))),pch=19, col=rgb(0,0,0,0.5))
if(ps[[10]]$pvalues[1] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=1.05, x1=1.95,length=0)
    text(ifelse(ps[[10]]$pvalues[1]>0.01, "*", ifelse(ps[[10]]$pvalues[1]>0.01, "**", "***")),x=1.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[2] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=2.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[2]>0.01, "*", ifelse(ps[[10]]$pvalues[2]>0.01, "**", "***")),x=2.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[3] < 0.05){
    arrows(y0=1.15*range.dep[2], y1=1.15*range.dep[2], x0=1.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[3]>0.01, "*", ifelse(ps[[10]]$pvalues[3]>0.01, "**", "***")),x=2, y=1.175*range.dep[2])
}
text("c", x=-0.1, y=45,xpd=NA,cex=2, font=2)

f1<-lme(Start_time~Exp_group2, random=~1|PairID, data=females, method="ML")
anova(f1, fnull)
ps <- summary(glht(f1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
range.dep <- c(4,36)#range(females$Start_time, na.rm=T)
xx <- boxplot(females$Start_time~females$Exp_group2, col=f.cols, names=NA, outline=F,
              ylim=c(range.dep[1] - 0.1*(diff(range.dep)), 1.2 * range.dep[2]),
              xlab="",ylab="Time to gonadogenesis (days)")
points(y=females$Start_time, x=jitter(as.numeric(as.factor(females$Exp_group2))),pch=19, col=rgb(0,0,0,0.5))
if(ps[[10]]$pvalues[1] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=1.05, x1=1.95,length=0)
    text(ifelse(ps[[10]]$pvalues[1]>0.01, "*", ifelse(ps[[10]]$pvalues[1]>0.01, "**", "***")),x=1.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[2] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=2.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[2]>0.01, "*", ifelse(ps[[10]]$pvalues[2]>0.01, "**", "***")),x=2.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[3] < 0.05){
    arrows(y0=1.15*range.dep[2], y1=1.15*range.dep[2], x0=1.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[3]>0.01, "*", ifelse(ps[[10]]$pvalues[3]>0.01, "**", "***")),x=2, y=1.175*range.dep[2])
}
text("d", x=-0.1, y=45,xpd=NA,cex=2, font=2)

mnull<-glmer(Gonads~1+(1|PairID), data=males, family="poisson")
fnull<-glmer(Gonads~1+(1|PairID), data=females, family="poisson")

m1<-glmer(Gonads~Exp_group2+(1|PairID), data=males, family="poisson")
anova(m1, mnull)

ps <- summary(glht(m1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
range.dep <- range(males$Gonads, na.rm=T)
xx <- boxplot(males$Gonads~males$Exp_group2, col=m.cols, names=NA, outline=F,
              ylim=c(range.dep[1] - 0.1*(diff(range.dep)), 1.2 * range.dep[2]),
              xlab="",ylab="No. testes")
points(y=males$Gonads, x=jitter(as.numeric(as.factor(males$Exp_group2))),pch=19, col=rgb(0,0,0,0.5))
if(ps[[10]]$pvalues[1] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=1.05, x1=1.95,length=0)
    text(ifelse(ps[[10]]$pvalues[1]>0.01, "*", ifelse(ps[[10]]$pvalues[1]>0.01, "**", "***")),x=1.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[2] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=2.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[2]>0.01, "*", ifelse(ps[[10]]$pvalues[2]>0.01, "**", "***")),x=2.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[3] < 0.05){
    arrows(y0=1.15*range.dep[2], y1=1.15*range.dep[2], x0=1.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[3]>0.01, "*", ifelse(ps[[10]]$pvalues[3]>0.01, "**", "***")),x=2, y=1.175*range.dep[2])
}
text("e", x=-0.1, y=30,xpd=NA,cex=2, font=2)

f1<-glmer(Gonads~Exp_group2+(1|PairID), data=females, family="poisson")
anova(f1,fnull)
ps <- summary(glht(f1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
range.dep <- range(females$Gonads, na.rm=T)
xx <- boxplot(females$Gonads~females$Exp_group2, col=f.cols, names=NA, outline=F,
              ylim=c(range.dep[1] - 0.1*(diff(range.dep)), 1.2 * range.dep[2]),
              xlab="",ylab="No. eggs")
points(y=females$Gonads, x=jitter(as.numeric(as.factor(females$Exp_group2))),pch=19, col=rgb(0,0,0,0.5))
if(ps[[10]]$pvalues[1] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=1.05, x1=1.95,length=0)
    text(ifelse(ps[[10]]$pvalues[1]>0.01, "*", ifelse(ps[[10]]$pvalues[1]>0.01, "**", "***")),x=1.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[2] < 0.05){
    arrows(y0=1.05*range.dep[2], y1=1.05*range.dep[2], x0=2.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[2]>0.01, "*", ifelse(ps[[10]]$pvalues[2]>0.01, "**", "***")),x=2.5, y=1.075*range.dep[2])
}
if(ps[[10]]$pvalues[3] < 0.05){
    arrows(y0=1.15*range.dep[2], y1=1.15*range.dep[2], x0=1.05, x1=2.95,length=0)
    text(ifelse(ps[[10]]$pvalues[3]>0.01, "*", ifelse(ps[[10]]$pvalues[3]>0.01, "**", "***")),x=2, y=1.175*range.dep[2])
}
text("f", x=-0.1, y=20,xpd=NA,cex=2, font=2)

mnull<-glmer(as.factor(Fate)~1+(1|PairID), data=males, family="binomial")
fnull<-glmer(as.factor(Fate)~1+(1|PairID), data=females, family="binomial")

m1<-glmer(as.factor(Fate)~Exp_group2+(1|PairID), data=males, family="binomial")

summary(m1) ## boundary singular message because the random effect is estimated to be 0
### run the same model without a random effect to check for consistency
m1.mod <- glm(as.factor(Fate)~Exp_group2, data=males, family="binomial")
summary(m1.mod)

ps <- summary(glht(m1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
m.tab <- table(males$Fate, males$Exp_group2)
x.pos <- barplot(prop.table(m.tab,2)[2,], col=m.cols,names.arg=names,ylab="Prop. surviving",ylim=c(-0.1,1.1))
axis(1,labels=FALSE,at=x.pos)
box()

max.height <- max(prop.table(m.tab,2)[2,])
if(ps[[10]]$pvalues[1] < 0.05){
    arrows(y0=1.05*max.height, y1=1.05*max.height, x0=0.75, x1=1.85,length=0)
    text(ifelse(ps[[10]]$pvalues[1]>0.01, "*", ifelse(ps[[10]]$pvalues[1]>0.01, "**", "***")),x=1.3, y=1.075*max.height)
}
if(ps[[10]]$pvalues[2] < 0.05){
    arrows(y0=1.05*max.height, y1=1.05*max.height, x0=1.95, x1=3.05,length=0)
    text(ifelse(ps[[10]]$pvalues[2]>0.01, "*", ifelse(ps[[10]]$pvalues[2]>0.01, "**", "***")),x=2.5, y=1.075*max.height)
}
if(ps[[10]]$pvalues[3] < 0.05){
    arrows(y0=1.15*max.height, y1=1.15*max.height, x0=0.75, x1=3.05,length=0)
    text(ifelse(ps[[10]]$pvalues[3]>0.01, "*", ifelse(ps[[10]]$pvalues[3]>0.01, "**", "***")),x=2, y=1.175*max.height)
}
text("g", x=-0.45, y=1.05,xpd=NA,cex=2, font=2)

f1<-glmer(as.factor(Fate)~Exp_group2+(1|PairID), data=females, family="binomial")
anova(f1, fnull)
summary(f1) ## boundary singular message because the random effect is estimated to be 0
### run the same model without a random effect to check for consistency
f1.mod <- glm(as.factor(Fate)~Exp_group2, data=females, family="binomial")
summary(f1.mod)

ps <- summary(glht(f1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))
f.tab <- table(females$Fate, females$Exp_group2)
x.pos <- barplot(prop.table(f.tab,2)[2,], col=f.cols,names.arg=names,ylab="Prop. surviving",ylim=c(-0.1,1.1))
axis(1,labels=FALSE,at=x.pos)
box()

if(ps[[10]]$pvalues[1] < 0.05){
    arrows(y0=1.05*max.height, y1=1.05*max.height, x0=0.75, x1=1.85,length=0)
    text(ifelse(ps[[10]]$pvalues[1]>0.01, "*", ifelse(ps[[10]]$pvalues[1]>0.01, "**", "***")),x=1.3, y=1.075*max.height)
}
if(ps[[10]]$pvalues[2] < 0.05){
    arrows(y0=1.05*max.height, y1=1.05*max.height, x0=1.95, x1=3.05,length=0)
    text(ifelse(ps[[10]]$pvalues[2]>0.01, "*", ifelse(ps[[10]]$pvalues[2]>0.01, "**", "***")),x=2.5, y=1.075*max.height)
}
if(ps[[10]]$pvalues[3] < 0.05){
    arrows(y0=1.15*max.height, y1=1.15*max.height, x0=0.75, x1=3.05,length=0)
    text(ifelse(ps[[10]]$pvalues[3]>0.01, "*", ifelse(ps[[10]]$pvalues[3]>0.01, "**", "***")),x=2, y=1.175*max.height)
}
text("h", x=-0.45, y=1.05,xpd=NA,cex=2, font=2)

dev.off()

#system("eog size_manipulation_results.tif")

males$RelGonads <- log(males$Gonads/males$post_area)
females$RelGonads <- log(females$Gonads/females$post_area)

mnull <- lme(RelGonads~1, random=~1|PairID, data=males,method="ML")
fnull <- lme(RelGonads~1, random=~1|PairID, data=females,method="ML")

m1<-lme(RelGonads~Exp_group2, random=~1|PairID, data=males,method="ML")
anova(m1, mnull)
ps <- summary(glht(m1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))

f1<-lme(RelGonads~Exp_group2, random=~1|PairID, data=females,method="ML")
anova(f1, fnull)
summary(f1)
ps <- summary(glht(f1,linfct=mcp(Exp_group2 = c("ctrl - reduced = 0", "enlarged - ctrl = 0", "enlarged - reduced = 0"))),
              test=adjusted("hochberg"))

mf <- rbind(females, males)
mf$Sex <- as.factor(mf$Sex)
levels(mf$Sex) <- c("MaleStrain","FemaleStrain")
mf$Sex <- relevel(mf$Sex, ref="FemaleStrain")
mf$Gonads_scaled <- NA
mf$Gonads_scaled[mf$Sex=="MaleStrain"] <- scale(log(mf$Gonads[mf$Sex=="MaleStrain"]))
mf$Gonads_scaled[mf$Sex=="FemaleStrain"] <- scale(log(mf$Gonads[mf$Sex=="FemaleStrain"]))

m1 <- lme(post_area~Sex*Exp_group2, random=~1|PairID, data=mf, method="ML")
m2 <- lme(Start_time~Sex*Exp_group2, random=~1|PairID, data=mf, method="ML")
m3 <- lme(Gonads_scaled~Sex*Exp_group2, random=~1|PairID, data=mf, method="ML")
m4 <- lme(RelGonads~Sex*Exp_group2, random=~1|PairID, data=mf, method="ML")
m5 <- glmer(as.factor(Fate)~Sex*Exp_group2+(1|PairID), data=mf, family="binomial")
## boundary singular warning message because random effect is estimated to be close to zero
## refitting without random effect
m5 <- glm(as.factor(Fate)~Sex*Exp_group2, data=mf, family="binomial")

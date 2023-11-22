#This code reproduces results in the paper titled:
#"How do trade-offs emerge? A test of the antagonistic pleiotropy hypothesis using a replicated artificial selection experiment"
#by Wilson, Wolak, Hiramatsu, Garland, and Careau

#For questions, please contact the corresponding author of the article:
#Vincent Careau, University of Ottawa, vcareau@uottawa.ca

#to run this code, you need the following files:
# "QG-AP-test_data.txt"
# "QG-AP-test_pedigree.txt""

#overview of the main sections of the code:
#section #1 - make Figure 1 -> requires data as .txt file
#section #2 - make Figure 2 -> requires data as .txt file
#section #3A- make the A-inverse to run animal models -> requires pedigree and data (.txt files)
#                                                     -> saves data + A-inverses as RData objects
#section #3B- run MCMCglmm models -> load .RData produced in section #3A
#                                 -> check model output as shown in Table S1
#                                 -> saves MCMCglmm models as RData
#section #3C- extract genetic correlations -> load MCMCglmm models produced in section #3B
#                                          -> saves posteriors as RData
#section #3D- make Figure 3 -> requires posteriors estimates for breeding values and genetic correlations
#section #4 - make Figure 4
#section #5 - make Figure S1
#section #6 - make Table S1
#section #7 - make Table S2



#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
#.################ section #1 - make figure 1 adaptive landscapes ###############
dat 	   <- read.table("QG-AP-test_data.txt", header = TRUE)

#show missing data for generations:
  plot(RPM56l~GEN,subset(dat,linetype==1), col=rgb(1,0,0,0.05),pch=16)
points(RPM56l~GEN,subset(dat,linetype==0), col=rgb(0,0,1,0.05),pch=16)
abline(v=32);abline(v=33);abline(v=34);abline(v=35);abline(v=52);abline(v=63);abline(v=64);abline(v=67)
GEN.list<-c(0:31,37:51,53:62,65,66,68:77)#list of GENs without running gdata

#calculate SELECTION GRADIENTS ###########################
DATA.C.WIS<-subset(dat,linetype=="0" & UNI=="WIS")
DATA.S.WIS<-subset(dat,linetype=="1" & UNI=="WIS")
#calculate relative fitness using the number of pups produced
DATA.C.WIS$relFIT<-DATA.C.WIS$pups/mean(DATA.C.WIS$pups)
DATA.S.WIS$relFIT<-DATA.S.WIS$pups/mean(DATA.S.WIS$pups)
#fit linear model
SEL.C.WIS<-lm(relFIT~RPM56l+INT56l+GEN+WHLSTAGE+sex+Fcoeff+line,DATA.C.WIS)
SEL.S.WIS<-lm(relFIT~RPM56l+INT56l+GEN+WHLSTAGE+sex+Fcoeff+line,DATA.S.WIS)
#extract betas OVERALL
BETA.RPM.C.WIS<-summary(SEL.C.WIS)$coefficients[2,1:2]
BETA.INT.C.WIS<-summary(SEL.C.WIS)$coefficients[3,1:2]
BETA.RPM.S.WIS<-summary(SEL.S.WIS)$coefficients[2,1:2]
BETA.INT.S.WIS<-summary(SEL.S.WIS)$coefficients[3,1:2]
############## estimate BETAS in the UCR generations #########################
DATA.C.UCR<-subset(dat,linetype=="0" & UNI=="UCR")
DATA.S.UCR<-subset(dat,linetype=="1" & UNI=="UCR")
#calculate relative fitness using the number of pups produced
DATA.C.UCR$relFIT<-DATA.C.UCR$pups/mean(DATA.C.UCR$pups)
DATA.S.UCR$relFIT<-DATA.S.UCR$pups/mean(DATA.S.UCR$pups)
#fit linear model
SEL.C.UCR<-lm(relFIT~RPM56l+INT56l+GEN+WHLSTAGE+sex+Fcoeff+line,DATA.C.UCR)
SEL.S.UCR<-lm(relFIT~RPM56l+INT56l+GEN+WHLSTAGE+sex+Fcoeff+line,DATA.S.UCR)
#extract betas OVERALL
BETA.RPM.C.UCR<-summary(SEL.C.UCR)$coefficients[2,1:2]
BETA.INT.C.UCR<-summary(SEL.C.UCR)$coefficients[3,1:2]
BETA.RPM.S.UCR<-summary(SEL.S.UCR)$coefficients[2,1:2]
BETA.INT.S.UCR<-summary(SEL.S.UCR)$coefficients[3,1:2]
################### estimate BETAS in each generation ##########################
DATA.C<-subset(dat,linetype=="0")
DATA.S<-subset(dat,linetype=="1")
summary(DATA.C)
BETAS<-data.frame(BLOCK=NA,GEN=0:77,
B.RPM.C=0,se.RPM.C=0,B.INT.C=0,se.INT.C=0,
B.RPM.S=0,se.RPM.S=0,B.INT.S=0,se.INT.S=0,
B.RUN.C=0,se.RUN.C=0,B.RUN.S=0,se.RUN.S=0)
#
for(i in GEN.list) {
#estimate gradients for the first block
BLOCK.C<-subset(dat,GEN==i & linetype=="0")
BLOCK.S<-subset(dat,GEN==i & linetype=="1")
#calculate relative fitness using the number of pups produced
BLOCK.C$relFIT<-BLOCK.C$pups/mean(BLOCK.C$pups)
BLOCK.S$relFIT<-BLOCK.S$pups/mean(BLOCK.S$pups)
#fit linear model
SEL.C<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+line,BLOCK.C)
SEL.S<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+line,BLOCK.S)
RUN.C<-lm(relFIT~RUN56l+WHLSTAGE+sex+line,BLOCK.C)
RUN.S<-lm(relFIT~RUN56l+WHLSTAGE+sex+line,BLOCK.S)
#extract betas
BETAS[which(BETAS$GEN==i),3:4] <-summary(SEL.C)$coefficients[2,1:2]
BETAS[which(BETAS$GEN==i),5:6] <-summary(SEL.C)$coefficients[3,1:2]
BETAS[which(BETAS$GEN==i),7:8] <-summary(SEL.S)$coefficients[2,1:2]
BETAS[which(BETAS$GEN==i),9:10]<-summary(SEL.S)$coefficients[3,1:2]
BETAS[which(BETAS$GEN==i),11:12]<-summary(RUN.C)$coefficients[2,1:2]
BETAS[which(BETAS$GEN==i),13:14]<-summary(RUN.S)$coefficients[2,1:2]
}
#
#selection gradients throughout the experiment
BETAS$n<-1
BETAS[!BETAS$GEN%in%GEN.list,3:11]<-NA
BETAS$GEN.C<-BETAS$GEN-0.05
BETAS$GEN.S<-BETAS$GEN+0.05
#
#calculate MEAN + SE BETAS
BETAS$UNI<-"UCR"
BETAS$UNI[which(BETAS$GEN<32)]<-"WIS"
BETAS[!BETAS$GEN%in%GEN.list,]<-NA
BETAS.WIS<-subset(BETAS,UNI=="WIS")
BETAS.UCR<-subset(BETAS,UNI=="UCR")
#
BETA.RUN.C.WIS<-cbind(mean(BETAS.WIS$B.RUN.C,na.rm=T),sd(BETAS.WIS$B.RUN.C,na.rm=T)/sqrt(sum(BETAS.WIS$n,na.rm=T)))
BETA.RUN.S.WIS<-cbind(mean(BETAS.WIS$B.RUN.S,na.rm=T),sd(BETAS.WIS$B.RUN.S,na.rm=T)/sqrt(sum(BETAS.WIS$n,na.rm=T)))
BETA.RPM.C.WIS<-cbind(mean(BETAS.WIS$B.RPM.C,na.rm=T),sd(BETAS.WIS$B.RPM.C,na.rm=T)/sqrt(sum(BETAS.WIS$n,na.rm=T)))
BETA.INT.C.WIS<-cbind(mean(BETAS.WIS$B.INT.C,na.rm=T),sd(BETAS.WIS$B.INT.C,na.rm=T)/sqrt(sum(BETAS.WIS$n,na.rm=T)))
BETA.RPM.S.WIS<-cbind(mean(BETAS.WIS$B.RPM.S,na.rm=T),sd(BETAS.WIS$B.RPM.S,na.rm=T)/sqrt(sum(BETAS.WIS$n,na.rm=T)))
BETA.INT.S.WIS<-cbind(mean(BETAS.WIS$B.INT.S,na.rm=T),sd(BETAS.WIS$B.INT.S,na.rm=T)/sqrt(sum(BETAS.WIS$n,na.rm=T)))
#
BETA.RUN.C.UCR<-cbind(mean(BETAS.UCR$B.RUN.C,na.rm=T),sd(BETAS.UCR$B.RUN.C,na.rm=T)/sqrt(sum(BETAS.UCR$n,na.rm=T)))
BETA.RUN.S.UCR<-cbind(mean(BETAS.UCR$B.RUN.S,na.rm=T),sd(BETAS.UCR$B.RUN.S,na.rm=T)/sqrt(sum(BETAS.UCR$n,na.rm=T)))
BETA.RPM.C.UCR<-cbind(mean(BETAS.UCR$B.RPM.C,na.rm=T),sd(BETAS.UCR$B.RPM.C,na.rm=T)/sqrt(sum(BETAS.UCR$n,na.rm=T)))
BETA.INT.C.UCR<-cbind(mean(BETAS.UCR$B.INT.C,na.rm=T),sd(BETAS.UCR$B.INT.C,na.rm=T)/sqrt(sum(BETAS.UCR$n,na.rm=T)))
BETA.RPM.S.UCR<-cbind(mean(BETAS.UCR$B.RPM.S,na.rm=T),sd(BETAS.UCR$B.RPM.S,na.rm=T)/sqrt(sum(BETAS.UCR$n,na.rm=T)))
BETA.INT.S.UCR<-cbind(mean(BETAS.UCR$B.INT.S,na.rm=T),sd(BETAS.UCR$B.INT.S,na.rm=T)/sqrt(sum(BETAS.UCR$n,na.rm=T)))
#
############ estimate BETAS in each generation for each LINE ###################
BETAS.LINE<-data.frame(BLOCK=NA,
                  GEN=0:77,
                  B.RPM.1=0,se.RPM.1=0,B.INT.1=0,se.INT.1=0,B.RUN.1=0,se.RUN.1=0,
                  B.RPM.2=0,se.RPM.2=0,B.INT.2=0,se.INT.2=0,B.RUN.2=0,se.RUN.2=0,
                  B.RPM.3=0,se.RPM.3=0,B.INT.3=0,se.INT.3=0,B.RUN.3=0,se.RUN.3=0,
                  B.RPM.4=0,se.RPM.4=0,B.INT.4=0,se.INT.4=0,B.RUN.4=0,se.RUN.4=0,
                  B.RPM.5=0,se.RPM.5=0,B.INT.5=0,se.INT.5=0,B.RUN.5=0,se.RUN.5=0,
                  B.RPM.6=0,se.RPM.6=0,B.INT.6=0,se.INT.6=0,B.RUN.6=0,se.RUN.6=0,
                  B.RPM.7=0,se.RPM.7=0,B.INT.7=0,se.INT.7=0,B.RUN.7=0,se.RUN.7=0,
                  B.RPM.8=0,se.RPM.8=0,B.INT.8=0,se.INT.8=0,B.RUN.8=0,se.RUN.8=0)
for(i in GEN.list) {
#estimate gradients for the first block
BLOCK.1<-subset(DATA.C,GEN==i & line==1)
BLOCK.2<-subset(DATA.C,GEN==i & line==2)
BLOCK.3<-subset(DATA.S,GEN==i & line==3)
BLOCK.4<-subset(DATA.C,GEN==i & line==4)
BLOCK.5<-subset(DATA.C,GEN==i & line==5)
BLOCK.6<-subset(DATA.S,GEN==i & line==6)
BLOCK.7<-subset(DATA.S,GEN==i & line==7)
BLOCK.8<-subset(DATA.S,GEN==i & line==8)
#calculate relative fitness using the number of pups produced
BLOCK.1$relFIT<-BLOCK.1$pups/mean(BLOCK.1$pups)
BLOCK.2$relFIT<-BLOCK.2$pups/mean(BLOCK.2$pups)
BLOCK.3$relFIT<-BLOCK.3$pups/mean(BLOCK.3$pups)
BLOCK.4$relFIT<-BLOCK.4$pups/mean(BLOCK.4$pups)
BLOCK.5$relFIT<-BLOCK.5$pups/mean(BLOCK.5$pups)
BLOCK.6$relFIT<-BLOCK.6$pups/mean(BLOCK.6$pups)
BLOCK.7$relFIT<-BLOCK.7$pups/mean(BLOCK.7$pups)
BLOCK.8$relFIT<-BLOCK.8$pups/mean(BLOCK.8$pups)
#fit linear model
SEL.1<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+Fcoeff,BLOCK.1)
SEL.2<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+Fcoeff,BLOCK.2)
SEL.3<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+Fcoeff,BLOCK.3)
SEL.4<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+Fcoeff,BLOCK.4)
SEL.5<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+Fcoeff,BLOCK.5)
SEL.6<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+Fcoeff,BLOCK.6)
SEL.7<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+Fcoeff,BLOCK.7)
SEL.8<-lm(relFIT~RPM56l+INT56l+WHLSTAGE+sex+Fcoeff,BLOCK.8)
SEL.1.RUN<-lm(relFIT~RUN56l+WHLSTAGE+sex+Fcoeff,BLOCK.1)
SEL.2.RUN<-lm(relFIT~RUN56l+WHLSTAGE+sex+Fcoeff,BLOCK.2)
SEL.3.RUN<-lm(relFIT~RUN56l+WHLSTAGE+sex+Fcoeff,BLOCK.3)
SEL.4.RUN<-lm(relFIT~RUN56l+WHLSTAGE+sex+Fcoeff,BLOCK.4)
SEL.5.RUN<-lm(relFIT~RUN56l+WHLSTAGE+sex+Fcoeff,BLOCK.5)
SEL.6.RUN<-lm(relFIT~RUN56l+WHLSTAGE+sex+Fcoeff,BLOCK.6)
SEL.7.RUN<-lm(relFIT~RUN56l+WHLSTAGE+sex+Fcoeff,BLOCK.7)
SEL.8.RUN<-lm(relFIT~RUN56l+WHLSTAGE+sex+Fcoeff,BLOCK.8)
#extract BETAS.line
BETAS.LINE[which(BETAS.LINE$GEN==i), 3:4 ] <-summary(SEL.1)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i), 5:6 ] <-summary(SEL.1)$coefficients[3,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i), 9:10] <-summary(SEL.2)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),11:12] <-summary(SEL.2)$coefficients[3,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),15:16] <-summary(SEL.3)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),17:18] <-summary(SEL.3)$coefficients[3,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),21:22] <-summary(SEL.4)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),23:24] <-summary(SEL.4)$coefficients[3,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),27:28] <-summary(SEL.5)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),29:30] <-summary(SEL.5)$coefficients[3,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),33:34] <-summary(SEL.6)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),35:36] <-summary(SEL.6)$coefficients[3,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),39:40] <-summary(SEL.7)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),41:42] <-summary(SEL.7)$coefficients[3,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),45:46] <-summary(SEL.8)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),47:48] <-summary(SEL.8)$coefficients[3,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),7:8] <-summary(SEL.1.RUN)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),13:14] <-summary(SEL.2.RUN)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),19:20] <-summary(SEL.3.RUN)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),25:26] <-summary(SEL.4.RUN)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),31:32] <-summary(SEL.5.RUN)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),37:38] <-summary(SEL.6.RUN)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),43:44] <-summary(SEL.7.RUN)$coefficients[2,1:2]
BETAS.LINE[which(BETAS.LINE$GEN==i),49:50] <-summary(SEL.8.RUN)$coefficients[2,1:2]
                  }
#
#selection gradients throughout the experiment
BETAS.LINE$n<-1
BETAS.LINE[!BETAS.LINE$GEN%in%GEN.list,]<-NA
BETAS.LINE$GEN.1<-BETAS.LINE$GEN-0.1
BETAS.LINE$GEN.2<-BETAS.LINE$GEN-0.075
BETAS.LINE$GEN.4<-BETAS.LINE$GEN-0.05
BETAS.LINE$GEN.5<-BETAS.LINE$GEN-0.025
BETAS.LINE$GEN.3<-BETAS.LINE$GEN+0.1
BETAS.LINE$GEN.6<-BETAS.LINE$GEN+0.075
BETAS.LINE$GEN.7<-BETAS.LINE$GEN+0.05
BETAS.LINE$GEN.8<-BETAS.LINE$GEN+0.025
#
#calculate MEAN + SE BETAS.LINE
BETA.RPM.1<-cbind(mean(BETAS.LINE$B.RPM.1,na.rm=T),sd(BETAS.LINE$B.RPM.1,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.INT.1<-cbind(mean(BETAS.LINE$B.INT.1,na.rm=T),sd(BETAS.LINE$B.INT.1,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.RPM.2<-cbind(mean(BETAS.LINE$B.RPM.2,na.rm=T),sd(BETAS.LINE$B.RPM.2,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.INT.2<-cbind(mean(BETAS.LINE$B.INT.2,na.rm=T),sd(BETAS.LINE$B.INT.2,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.RPM.3<-cbind(mean(BETAS.LINE$B.RPM.3,na.rm=T),sd(BETAS.LINE$B.RPM.3,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.INT.3<-cbind(mean(BETAS.LINE$B.INT.3,na.rm=T),sd(BETAS.LINE$B.INT.3,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.RPM.4<-cbind(mean(BETAS.LINE$B.RPM.4,na.rm=T),sd(BETAS.LINE$B.RPM.4,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.INT.4<-cbind(mean(BETAS.LINE$B.INT.4,na.rm=T),sd(BETAS.LINE$B.INT.4,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.RPM.5<-cbind(mean(BETAS.LINE$B.RPM.5,na.rm=T),sd(BETAS.LINE$B.RPM.5,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.INT.5<-cbind(mean(BETAS.LINE$B.INT.5,na.rm=T),sd(BETAS.LINE$B.INT.5,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.RPM.6<-cbind(mean(BETAS.LINE$B.RPM.6,na.rm=T),sd(BETAS.LINE$B.RPM.6,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.INT.6<-cbind(mean(BETAS.LINE$B.INT.6,na.rm=T),sd(BETAS.LINE$B.INT.6,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.RPM.7<-cbind(mean(BETAS.LINE$B.RPM.7,na.rm=T),sd(BETAS.LINE$B.RPM.7,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.INT.7<-cbind(mean(BETAS.LINE$B.INT.7,na.rm=T),sd(BETAS.LINE$B.INT.7,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.RPM.8<-cbind(mean(BETAS.LINE$B.RPM.8,na.rm=T),sd(BETAS.LINE$B.RPM.8,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))
BETA.INT.8<-cbind(mean(BETAS.LINE$B.INT.8,na.rm=T),sd(BETAS.LINE$B.INT.8,na.rm=T)/sqrt(sum(BETAS.LINE$n,na.rm=T)))

# use thin-plate splines to display the adaptive landscape
library(fields)
LINE.3<-subset(dat,line==3)
LINE.6<-subset(dat,line==6)
LINE.7<-subset(dat,line==7)
LINE.8<-subset(dat,line==8)
LINE.3$relFIT<-LINE.3$pups/mean(LINE.3$pups)
LINE.6$relFIT<-LINE.6$pups/mean(LINE.6$pups)
LINE.7$relFIT<-LINE.7$pups/mean(LINE.7$pups)
LINE.8$relFIT<-LINE.8$pups/mean(LINE.8$pups)

     TPS.3<-Tps(data.matrix(data.frame(LINE.3$INT56,LINE.3$RPM56)),LINE.3$relFIT)
save(TPS.3,file="TPS.3.raw.RData")
     TPS.6<-Tps(data.matrix(data.frame(LINE.6$INT56,LINE.6$RPM56)),LINE.6$relFIT)
save(TPS.6,file="TPS.6.raw.RData")
     TPS.7<-Tps(data.matrix(data.frame(LINE.7$INT56,LINE.7$RPM56)),LINE.7$relFIT)
save(TPS.7,file="TPS.7.raw.RData")
     TPS.8<-Tps(data.matrix(data.frame(LINE.8$INT56,LINE.8$RPM56)),LINE.8$relFIT)
save(TPS.8,file="TPS.8.raw.RData")

load(file="TPS.3.raw.RData")
load(file="TPS.6.raw.RData")
load(file="TPS.7.raw.RData")
load(file="TPS.8.raw.RData")
sp.3<-predictSurface(TPS.3)
sp.6<-predictSurface(TPS.6)
sp.7<-predictSurface(TPS.7)
sp.8<-predictSurface(TPS.8)


#make Figure 1
dev.new(width=9,height=4, units = "cm")
par(mfrow=c(6,5),las=1, oma=c(4,2,0.5,0), mar=c(0,2,2,2))
layout(matrix(c(1,1,1,2,3,
                1,1,1,2,3,
                4,4,4,2,3,
                4,4,4,5,6,
                7,7,7,5,6,
                7,7,7,5,6), 6, 5, byrow = TRUE))
#layout.show(7)
  plot(B.RUN.S~GEN, BETAS, pch=16, col="red", cex=1,type="o", ylim=c(-6,14),xlim=c(1,77),xlab="",xaxt="n",ylab="")
points(B.RUN.C~GEN.C, BETAS, pch=16, col="blue", cex=1,type="o")
rect(0, BETA.RUN.S.WIS[1]-1.96*BETA.RUN.S.WIS[2],31,BETA.RUN.S.WIS[1]+1.96*BETA.RUN.S.WIS[2], col=rgb(1,0,0,0.25),border=F)
rect(35,BETA.RUN.S.UCR[1]-1.96*BETA.RUN.S.UCR[2],78,BETA.RUN.S.UCR[1]+1.96*BETA.RUN.S.UCR[2], col=rgb(1,0,0,0.25),border=F)
rect(0, BETA.RUN.C.WIS[1]-1.96*BETA.RUN.C.WIS[2],31,BETA.RUN.C.WIS[1]+1.96*BETA.RUN.C.WIS[2], col=rgb(0,0,1,0.25),border=F)
rect(35,BETA.RUN.C.UCR[1]-1.96*BETA.RUN.C.UCR[2],78,BETA.RUN.C.UCR[1]+1.96*BETA.RUN.C.UCR[2], col=rgb(0,0,1,0.25),border=F)
arrows(x0=BETAS$GEN.C,y0=BETAS$B.RUN.C-BETAS$se.RUN.C,x1=BETAS$GEN.C,y1=BETAS$B.RUN.C+BETAS$se.RUN.C,col="blue",code=3,angle=90,length=0)
arrows(x0=BETAS$GEN.S,y0=BETAS$B.RUN.S-BETAS$se.RUN.S,x1=BETAS$GEN.S,y1=BETAS$B.RUN.S+BETAS$se.RUN.S,col="red", code=3,angle=90,length=0)
points(B.RUN.1~GEN, BETAS.LINE, pch="1", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.RUN.2~GEN, BETAS.LINE, pch="2", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.RUN.4~GEN, BETAS.LINE, pch="4", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.RUN.5~GEN, BETAS.LINE, pch="5", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.RUN.3~GEN, BETAS.LINE, pch="3", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.RUN.6~GEN, BETAS.LINE, pch="6", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.RUN.7~GEN, BETAS.LINE, pch="7", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.RUN.8~GEN, BETAS.LINE, pch="8", col=rgb(1,0,0,0.25), cex=1,type="l")
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
abline(h=0, lty=2)
mtext("A", side=3, adj=0.025, line=-1.5)
legend(47,14, c("control","selected"),pch=c(16,17), col=c("blue","red"),bty="n")
mtext("Wisconsin generations", line=0.4, adj=0.10, side=3, cex=1)
mtext("Riverside generations", line=0.4, adj=0.80, side=3, cex=1)
#
image(sp.3,col=terrain.colors(1000),zlim=c(-0,3))
mtext("D", side=3, adj=0.025, line=-1.5)
mtext("line S3", side=1,cex=0.75,line=-1)
image(sp.6,col=terrain.colors(1000),zlim=c(-0,3))
mtext("E", side=3, adj=0.025, line=-1.5)
mtext("line S6", side=1,cex=0.75,line=-1)
#
  plot(B.RPM.S~GEN.C, BETAS, pch=16, col="red", cex=1,type="o", ylim=c(-6,14),xlim=c(1,77),xlab="",xaxt="n",ylab="")
points(B.RPM.C~GEN.C, BETAS, pch=16, col="blue", cex=1,type="o")
rect(0, BETA.RPM.S.WIS[1]-1.96*BETA.RPM.S.WIS[2],31,BETA.RPM.S.WIS[1]+1.96*BETA.RPM.S.WIS[2], col=rgb(1,0,0,0.25),border=F)
rect(35,BETA.RPM.S.UCR[1]-1.96*BETA.RPM.S.UCR[2],78,BETA.RPM.S.UCR[1]+1.96*BETA.RPM.S.UCR[2], col=rgb(1,0,0,0.25),border=F)
rect(0, BETA.RPM.C.WIS[1]-1.96*BETA.RPM.C.WIS[2],31,BETA.RPM.C.WIS[1]+1.96*BETA.RPM.C.WIS[2], col=rgb(0,0,1,0.25),border=F)
rect(35,BETA.RPM.C.UCR[1]-1.96*BETA.RPM.C.UCR[2],78,BETA.RPM.C.UCR[1]+1.96*BETA.RPM.C.UCR[2], col=rgb(0,0,1,0.25),border=F)
arrows(x0=BETAS$GEN.C,y0=BETAS$B.RPM.C-BETAS$se.RPM.C,x1=BETAS$GEN.C,y1=BETAS$B.RPM.C+BETAS$se.RPM.C,col="blue",code=3,angle=90,length=0)
arrows(x0=BETAS$GEN.S,y0=BETAS$B.RPM.S-BETAS$se.RPM.S,x1=BETAS$GEN.S,y1=BETAS$B.RPM.S+BETAS$se.RPM.S,col="red", code=3,angle=90,length=0)
points(B.RPM.1~GEN, BETAS.LINE, pch="1", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.RPM.2~GEN, BETAS.LINE, pch="2", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.RPM.4~GEN, BETAS.LINE, pch="4", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.RPM.5~GEN, BETAS.LINE, pch="5", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.RPM.3~GEN, BETAS.LINE, pch="3", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.RPM.6~GEN, BETAS.LINE, pch="6", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.RPM.7~GEN, BETAS.LINE, pch="7", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.RPM.8~GEN, BETAS.LINE, pch="8", col=rgb(1,0,0,0.25), cex=1,type="l")
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
abline(h=0, lty=2)
mtext("B", side=3, adj=0.025, line=-1.5)
#
image(sp.7,col=terrain.colors(1000),zlim=c(-0,3))
mtext("F", side=3, adj=0.025, line=-1.5)
mtext("line S7", side=1,cex=0.75,line=-1)
mtext("Running speed (revs/min)",  side=2, line=2.5,cex=1.25,adj=-0.5,las=3)
image(sp.8,col=terrain.colors(1000),zlim=c(-0,3))
mtext("G", side=3, adj=0.025, line=-1.5)
mtext("line S8", side=1,cex=0.75,line=-1)
#
  plot(B.INT.C~GEN.C, BETAS, pch=16, col="white", cex=1,type="o", ylim=c(-6,14),xlim=c(1,77),xlab="",xaxt="n",ylab="")
rect(0, BETA.INT.S.WIS[1]-1.96*BETA.INT.S.WIS[2],31,BETA.INT.S.WIS[1]+1.96*BETA.INT.S.WIS[2], col=rgb(1,0,0,0.25),border=F)
rect(35,BETA.INT.S.UCR[1]-1.96*BETA.INT.S.UCR[2],78,BETA.INT.S.UCR[1]+1.96*BETA.INT.S.UCR[2], col=rgb(1,0,0,0.25),border=F)
rect(0, BETA.INT.C.WIS[1]-1.96*BETA.INT.C.WIS[2],31,BETA.INT.C.WIS[1]+1.96*BETA.INT.C.WIS[2], col=rgb(0,0,1,0.25),border=F)
rect(35,BETA.INT.C.UCR[1]-1.96*BETA.INT.C.UCR[2],78,BETA.INT.C.UCR[1]+1.96*BETA.INT.C.UCR[2], col=rgb(0,0,1,0.25),border=F)
points(B.INT.C~GEN.C, BETAS, pch=16, col="blue", cex=1,type="o")
points(B.INT.S~GEN.S, BETAS, pch=17, col="red", cex=1,type="o")
arrows(x0=BETAS$GEN.C,y0=BETAS$B.INT.C-BETAS$se.INT.C,x1=BETAS$GEN.C,y1=BETAS$B.INT.C+BETAS$se.INT.C,col="blue",code=3,angle=90,length=0)
arrows(x0=BETAS$GEN.S,y0=BETAS$B.INT.S-BETAS$se.INT.S,x1=BETAS$GEN.S,y1=BETAS$B.INT.S+BETAS$se.INT.S,col="red", code=3,angle=90,length=0)
points(B.INT.1~GEN, BETAS.LINE, pch="1", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.INT.2~GEN, BETAS.LINE, pch="2", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.INT.4~GEN, BETAS.LINE, pch="4", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.INT.5~GEN, BETAS.LINE, pch="5", col=rgb(0,0,1,0.25), cex=1,type="l")
points(B.INT.3~GEN, BETAS.LINE, pch="3", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.INT.6~GEN, BETAS.LINE, pch="6", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.INT.7~GEN, BETAS.LINE, pch="7", col=rgb(1,0,0,0.25), cex=1,type="l")
points(B.INT.8~GEN, BETAS.LINE, pch="8", col=rgb(1,0,0,0.25), cex=1,type="l")
abline(h=0, lty=2)
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
mtext("C", side=3, adj=0.025, line=-1.5)
#
mtext("Generation",side=1, line=2.5,cex=1.25)
mtext("Selection gradients (±se)",side=2, line=0.5,cex=1.25,outer=T,las=3)
mtext("Running duration (min/day)",side=1, line=2.5,cex=1.25,adj=2.5)

rm(list=ls())  # delete all objects, to start fresh in new section below
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
#.################ section #2 - make figure 2 responses to selection ###############
dat 	   <- read.table("QG-AP-test_data.txt", header = TRUE)
dat$n<-1
MEANS       <-aggregate(RPM56~GEN+line,dat,FUN=mean)
MEANS$SD.RPM<-aggregate(RPM56~GEN+line,dat,FUN=sd)[3]
MEANS$n     <-aggregate(n~GEN+line,dat,FUN=sum)$n
MEANS$SE.RPM<-MEANS$SD.RPM/sqrt(MEANS$n)
MEANS$INT   <-aggregate(INT56~GEN+line,dat,FUN=mean)$INT
MEANS$SD.INT<-aggregate(INT56~GEN+line,dat,FUN=sd)[3]
MEANS$n     <-aggregate(n~GEN+line,dat,FUN=sum)$n
MEANS$SE.INT<-MEANS$SD.INT/sqrt(MEANS$n)
colnames(MEANS)<-c("GEN","LINE","RPM","SD.RPM","N","SE.RPM","INT","SD.INT","SE.INT")
head(MEANS)
#MEANS<-MEANS[-which(MEANS$GEN=="-1"),]
MEANS$linetype<-"S"
MEANS$linetype[which(MEANS$LINE=="1")]<-"C"
MEANS$linetype[which(MEANS$LINE=="2")]<-"C"
MEANS$linetype[which(MEANS$LINE=="4")]<-"C"
MEANS$linetype[which(MEANS$LINE=="5")]<-"C"
#
MEANS       <-aggregate(RPM56~GEN+line,dat,FUN=mean)
MEANS$SD.RPM<-aggregate(RPM56~GEN+line,dat,FUN=sd)[3]
MEANS$n     <-aggregate(n~GEN+line,dat,FUN=sum)$n
MEANS$SE.RPM<-MEANS$SD.RPM/sqrt(MEANS$n)
MEANS$INT   <-aggregate(INT56~GEN+line,dat,FUN=mean)$INT56
MEANS$SD.INT <-aggregate(INT56~GEN+line,dat,FUN=sd)[3]
MEANS$n      <-aggregate(n~GEN+line,dat,FUN=sum)$n
MEANS$SE.INT<-MEANS$SD.INT/sqrt(MEANS$n)
#
MEANS$RUN     <-aggregate(RUN56~GEN+line,dat,FUN=mean)$RUN56
MEANS$SD.RUN56<-aggregate(RUN56~GEN+line,dat,FUN=sd)[3]
MEANS$n       <-aggregate(n~GEN+line,dat,FUN=sum)$n
MEANS$SE.RUN56<-MEANS$SD.RUN56/sqrt(MEANS$n)
colnames(MEANS)<-c("GEN","LINE","RPM","SD.RPM","N","SE.RPM","INT","SD.INT","SE.INT","RUN","SD.RUN","SE.RUN")
head(MEANS)
#MEANS<-MEANS[-which(MEANS$GEN=="-1"),]
MEANS$linetype<-"S"
MEANS$linetype[which(MEANS$LINE=="1")]<-"C"
MEANS$linetype[which(MEANS$LINE=="2")]<-"C"
MEANS$linetype[which(MEANS$LINE=="4")]<-"C"
MEANS$linetype[which(MEANS$LINE=="5")]<-"C"
#
#add RATIO
MEANS.C       <-aggregate(RPM56~GEN,subset(dat,linetype==0),FUN=mean)
MEANS.C$INT56 <-aggregate(INT56~GEN,subset(dat,linetype==0),FUN=mean)$INT56
MEANS.C$RUN56 <-aggregate(RUN56~GEN,subset(dat,linetype==0),FUN=mean)$RUN56
MEANS.S       <-aggregate(RPM56~GEN,subset(dat,linetype==1),FUN=mean)
MEANS.S$INT56 <-aggregate(INT56~GEN,subset(dat,linetype==1),FUN=mean)$INT56
MEANS.S$RUN56 <-aggregate(RUN56~GEN,subset(dat,linetype==1),FUN=mean)$RUN56
RATIOS<-data.frame(GEN=MEANS$GEN[1:71], RPM=MEANS.S$RPM56/MEANS.C$RPM56,INT=MEANS.S$INT56/MEANS.C$INT56,RUN=MEANS.S$RUN56/MEANS.C$RUN56)

dev.new(width=9,height=6, units = "cm")
par(mfrow=c(3,1),las=1, mar=c(2,3,0,0), oma=c(2,4,2.5,5))
  plot(RUN~GEN, subset(MEANS, LINE=="1"), pch=16,lty=1,col="blue",type="l",ylim=c(6,15000),xlim=c(0,78), xaxt="n")
points(RUN~GEN, subset(MEANS, LINE=="2"), pch=16,lty=1,col="blue",type="l")
points(RUN~GEN, subset(MEANS, LINE=="4"), pch=16,lty=1,col="blue",type="l")
points(RUN~GEN, subset(MEANS, LINE=="5"), pch=16,lty=1,col="blue",type="l")
points(RUN~GEN, subset(MEANS, LINE=="3"), pch=17,lty=1,col="red",type="l")
points(RUN~GEN, subset(MEANS, LINE=="6"), pch=17,lty=1,col="red",type="l")
points(RUN~GEN, subset(MEANS, LINE=="7"), pch=17,lty=1,col="red",type="l")
points(RUN~GEN, subset(MEANS, LINE=="8"), pch=17,lty=1,col="red",type="l")
  mtext("Wisconsin generations", line = 0.75, adj = 0.12, side = 3, cex = 1)
  mtext("Riverside generations", line =  0.75, adj = 0.80, side = 3, cex = 1)
legend(-4,14000,c("Control lines","Selected lines","ratio"), cex=1,lty=c(1,1,1),lwd=c(1,1,2), col=c("blue","red","black"),bty="n")
#
par(new=T)
  plot(RUN~GEN, RATIOS,axes=F,type="l",lwd=2); axis(side = 4); mtext("Ratio", side = 4, line = 3,las=3)
#
rect(xleft=31.1,ybottom=min(RATIOS$RUN),xright=35.9,ytop=max(RATIOS$RUN), col = "white", border = "white")
rect(xleft=51.1,ybottom=min(RATIOS$RUN),xright=52.9,ytop=max(RATIOS$RUN), col = "white", border = "white")
rect(xleft=62.1,ybottom=min(RATIOS$RUN),xright=64.9,ytop=max(RATIOS$RUN), col = "white", border = "white")
rect(xleft=66.1,ybottom=min(RATIOS$RUN),xright=67.9,ytop=max(RATIOS$RUN), col = "white", border = "white")
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
mtext("Distance run",side=2,las=3,line=5)
mtext("(revs/day)",side=2,las=3,line=3.5)
mtext("A",side=3,adj=0.015,line=-1.5)
#
  plot(RPM~GEN, subset(MEANS, LINE=="1"), pch=16,lty=1,col="blue",type="l",ylim=c(0,30),xlim=c(0,78), xaxt="n")
points(RPM~GEN, subset(MEANS, LINE=="2"), pch=16,lty=1,col="blue",type="l")
points(RPM~GEN, subset(MEANS, LINE=="4"), pch=16,lty=1,col="blue",type="l")
points(RPM~GEN, subset(MEANS, LINE=="5"), pch=16,lty=1,col="blue",type="l")
points(RPM~GEN, subset(MEANS, LINE=="3"), pch=17,lty=1,col="red",type="l")
points(RPM~GEN, subset(MEANS, LINE=="6"), pch=17,lty=1,col="red",type="l")
points(RPM~GEN, subset(MEANS, LINE=="7"), pch=17,lty=1,col="red",type="l")
points(RPM~GEN, subset(MEANS, LINE=="8"), pch=17,lty=1,col="red",type="l")
#
par(new=T)
plot(RPM~GEN, RATIOS,axes=F,type="l",lwd=2); axis(side = 4); mtext("Ratio", side = 4, line = 3,las=3)
#
rect(xleft=31.1,ybottom=min(RATIOS$RPM),xright=35.9,ytop=max(RATIOS$RPM), col = "white", border = "white")
rect(xleft=51.1,ybottom=min(RATIOS$RPM),xright=52.9,ytop=max(RATIOS$RPM), col = "white", border = "white")
rect(xleft=62.1,ybottom=min(RATIOS$RPM),xright=64.9,ytop=max(RATIOS$RPM), col = "white", border = "white")
rect(xleft=66.1,ybottom=min(RATIOS$RPM),xright=67.9,ytop=max(RATIOS$RPM), col = "white", border = "white")
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
mtext("Running speed",side=2,las=3,line=5)
mtext("(revs/min)",side=2,las=3,line=3.5)
mtext("B",side=3,adj=0.015,line=-1.5)
#
  plot(INT~GEN, subset(MEANS, LINE=="1"), pch=16,lty=1,col="blue",type="l", ylim=c(0,650),xlim=c(0,78),xaxt="n")
points(INT~GEN, subset(MEANS, LINE=="2"), pch=16,lty=1,col="blue",type="l")
points(INT~GEN, subset(MEANS, LINE=="4"), pch=16,lty=1,col="blue",type="l")
points(INT~GEN, subset(MEANS, LINE=="5"), pch=16,lty=1,col="blue",type="l")
points(INT~GEN, subset(MEANS, LINE=="3"), pch=17,lty=1,col="red",type="l")
points(INT~GEN, subset(MEANS, LINE=="6"), pch=17,lty=1,col="red",type="l")
points(INT~GEN, subset(MEANS, LINE=="7"), pch=17,lty=1,col="red",type="l")
points(INT~GEN, subset(MEANS, LINE=="8"), pch=17,lty=1,col="red",type="l")
#
par(new=T)
plot(INT~GEN, RATIOS,axes=F,type="l",lwd=2) ; axis(side = 4); mtext("Ratio", side = 4, line = 3,las=3)
#
rect(xleft=31.1,ybottom=min(RATIOS$INT),xright=35.9,ytop=max(RATIOS$INT), col = "white", border = "white")
rect(xleft=51.1,ybottom=min(RATIOS$INT),xright=52.9,ytop=max(RATIOS$INT), col = "white", border = "white")
rect(xleft=62.1,ybottom=min(RATIOS$INT),xright=64.9,ytop=max(RATIOS$INT), col = "white", border = "white")
rect(xleft=66.1,ybottom=min(RATIOS$INT),xright=67.9,ytop=max(RATIOS$INT), col = "white", border = "white")
mtext("C",side=3,adj=0.015,line=-1.5)
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
mtext("Running duration",side=2,las=3,line=5)
mtext("(min/day)",side=2,las=3,line=3.5)
#
#legend("bottomleft",c("Control lines","Selected lines"), cex=1,lty=c(1,1), col=c("blue","red"),bty="n")
mtext("Generation", side=1,  line=2.5)
#

#calculate average S/C ratios after generation 25
mean(RATIOS$RUN[which(RATIOS$GEN>24)])
mean(RATIOS$RPM[which(RATIOS$GEN>24)])
mean(RATIOS$INT[which(RATIOS$GEN>24)])

#calculate average S/C ratios for UCR gens
mean(RATIOS$RUN[which(RATIOS$GEN>30)])
mean(RATIOS$RPM[which(RATIOS$GEN>30)])
mean(RATIOS$INT[which(RATIOS$GEN>30)])

rm(list=ls())  # delete all objects, to start fresh in new section below
#.############# section #3A - make A-inverse to run animal models#################
#.############# section #3A - make A-inverse to run animal models#################
#.############# section #3A - make A-inverse to run animal models#################
#.############# section #3A - make A-inverse to run animal models#################
#.############# section #3A - make A-inverse to run animal models#################
#.############# section #3A - make A-inverse to run animal models#################
#.############# section #3A - make A-inverse to run animal models#################
#.############# section #3A - make A-inverse to run animal models#################
#.############# section #3A - make A-inverse to run animal models#################
library(MCMCglmm)
library(nadiv)
PED_pruned <-read.table("QG-AP-test_pedigree.txt", header = TRUE)
dat 	   <- read.table("QG-AP-test_data.txt", header = TRUE)

#CONTROL MICE
BLOCK<-subset(dat,GEN >= -1 & linetype=="0")
BLOCK$animal  <-factor(BLOCK$animal)
BLOCK$line    <-factor(BLOCK$line)
BLOCK$damid   <-factor(BLOCK$damid)
BLOCK$GENfac  <-factor(BLOCK$GEN)
PED.tmp.2  <-subset(PED_pruned, GEN > -2)
PED.tmp.1 <-prepPed(PED.tmp.2)
PED      <-prunePed(PED.tmp.1, unique(BLOCK$animal))
AINVout <- inverseA(PED[, c("animal", "dam", "sire")])
fout <- AINVout$inbreeding
AINV <- AINVout$Ainv
stopifnot(all(abs(fout[match(BLOCK$animal, PED$animal)] - BLOCK$Fcoeff) < 1e-6))
C0Ainv <- AINV
C0datf <- BLOCK
C0ped <- PED
rm(list=c("AINV","AINVout","PED","PED.tmp.1","PED.tmp.2","BLOCK","fout"))

#SELECTED MICE
BLOCK<-subset(dat,GEN >= -1 & linetype=="1")
BLOCK$animal <-factor(BLOCK$animal)
BLOCK$line   <-factor(BLOCK$line)
BLOCK$damid  <-factor(BLOCK$damid)
BLOCK$GENfac <-factor(BLOCK$GEN)
PED.tmp.2  <-subset(PED_pruned, GEN > -2)
PED.tmp.1 <-prepPed(PED.tmp.2)
PED      <-prunePed(PED.tmp.1, unique(BLOCK$animal))
AINVout <- inverseA(PED[, c("animal", "dam", "sire")])
fout <- AINVout$inbreeding
AINV <- AINVout$Ainv
stopifnot(all(abs(fout[match(BLOCK$animal, PED$animal)] - BLOCK$Fcoeff) < 1e-6))
HR1Ainv <- AINV
HR1datf <- BLOCK
HR1ped <- PED
rm(list=c("AINV","AINVout","PED","PED.tmp.1","PED.tmp.2","BLOCK","fout","PED_pruned"))

nrow(HR1datf)
nrow(C0datf)
nrow(HR1datf)+nrow(C0datf)


#save as RData for section below
save(HR1datf, C0datf,HR1Ainv, C0Ainv, file="QG-AP-test_Rdata.RData")
rm(list=ls())  # delete all objects, to start fresh in section below

#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#load data and pedigree objects prepared in section #1
load(file="QG-AP-test_Rdata.RData")

#get data ready for running the models
C0datf$WSTRTymd <-as.Date(C0datf$WSTRTymd, "%Y-%m-%d")
HR1datf$WSTRTymd<-as.Date(HR1datf$WSTRTymd, "%Y-%m-%d")
C0datf$Mnth     <- factor(format(C0datf$WSTRTymd, "%m")) #make a month variable as factor
HR1datf$Mnth    <- factor(format(HR1datf$WSTRTymd, "%m")) #make a month variable as factor


library(MCMCglmm)
#Priors
#multivariate extension of parameter expanded/central-scaled F - gives FLAT prior on the correlations
k <- 2
kPEflatR <- list(R = list(V = diag(k), nu = 0),
  G = list(G1 = list(V = diag(k)*0.02, nu = k+1, alpha.mu = rep(0, k), alpha.V = diag(k)*1000),
    G2 = list(V = diag(k)*0.02, nu = k+1, alpha.mu = rep(0, k), alpha.V = diag(k)*1000)))

#set number of iterations
nsamp <- 3000
THIN <- 250
BURN <- 10000
(NITT <- BURN + nsamp*THIN)

#CONTROL MICE - run bivariate model
C_RPM_INT_mnth <- MCMCglmm(cbind(RPM56l, INT56l) ~ trait +
                                                   trait:sex +
                                                   trait:line +
                                                   trait:GENfac +
                                                   trait:WHLSTAGE +
                                                   trait:Fcoeff +
                                                   trait:Mnth,
                 random = ~ us(trait):animal + us(trait):damid,
                 rcov = ~ us(trait):units,
                 data = C0datf,
                 ginverse = list(animal = C0Ainv),
                 prior = kPEflatR,
                 family = rep("gaussian", 2),
                 nitt = NITT, thin = THIN, burnin = BURN,
                 pr = TRUE, saveX = TRUE, saveZ = TRUE)
summary(C_RPM_INT_mnth)    #these numbers are in Table S1

#SELECTED MICE  - run bivariate model
#set number of iterations
nsamp <- 3000
THIN <- 500
BURN <- 10000
(NITT <- BURN + nsamp*THIN)

HR_RPM_INT_mnth<- MCMCglmm(cbind(RPM56l, INT56l) ~ trait +
                                                   trait:sex +
                                                   trait:line +
                                                   trait:GENfac +
                                                   trait:WHLSTAGE +
                                                   trait:Fcoeff +
                                                   trait:Mnth,
                 random = ~ us(trait):animal + us(trait):damid,
                 rcov = ~ us(trait):units,
                 data = HR1datf,
                 ginverse = list(animal = HR1Ainv),
                 prior = kPEflatR,
                 family = rep("gaussian", 2),
                 nitt = NITT, thin = THIN, burnin = BURN,
                 pr = TRUE, saveX = TRUE, saveZ = TRUE)
summary(HR_RPM_INT_mnth)  #these numbers are in Table S1



#also need to run univariate models with RUN56l as response variable
# Priors
## parameter expanded/central-scaled F
univPEflatR <- list(R = list(V = diag(1), nu = 0),
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
    G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))


nsamp <- 3000
THIN <- 250
BURN <- 10000
(NITT <- BURN + nsamp*THIN)
C_RUN_mnth <- MCMCglmm(RUN56l ~ 1 + sex + line + GENfac + WHLSTAGE + Fcoeff + Mnth,
  random = ~ animal + damid,
  data = C0datf,
  ginverse = list(animal = C0Ainv),
  prior = univPEflatR,
  family = rep("gaussian", 1),
  nitt = NITT, thin = THIN, burnin = BURN,
  pr = TRUE, saveX = TRUE, saveZ = TRUE)


nsamp <- 3000
THIN <- 400
BURN <- 10000
(NITT <- BURN + nsamp*THIN)
HR_RUN_mnth <- MCMCglmm(RUN56l ~ 1 + sex + line + GENfac + WHLSTAGE + Fcoeff + Mnth,
  random = ~ animal + damid,
  data = HR1datf,
  ginverse = list(animal = HR1Ainv),
  prior = univPEflatR,
  family = rep("gaussian", 1),
  nitt = NITT, thin = THIN, burnin = BURN,
  pr = TRUE, saveX = TRUE, saveZ = TRUE)


#save models as RData
save(C_RPM_INT_mnth, HR_RPM_INT_mnth,
  file="QG-AP-test_bivariate_models.RData")
save(C_RUN_mnth, HR_RUN_mnth,      
  file="QG-AP-test_univariate_models.RData")
rm(list=ls())  # delete all objects, to start fresh section below

#.######### section #3C - breeding values and genetic correlations ################
#.######### section #3C - breeding values and genetic correlations ################
#.######### section #3C - breeding values and genetic correlations ################
#.######### section #3C - breeding values and genetic correlations ################
#.######### section #3C - breeding values and genetic correlations ################
#.######### section #3C - breeding values and genetic correlations ################
#.######### section #3C - breeding values and genetic correlations ################
#.######### section #3C - breeding values and genetic correlations ################
#.######### section #3C - breeding values and genetic correlations ################

library(MCMCglmm)
#load models fitted in section #2
load(file="QG-AP-test_bivariate_models.RData")
#load data (needed to assign generations and lines and breeding values)
load(file="QG-AP-test_Rdata.RData")


traits <- c("RPM56l", "INT56l")

# Create a function to calculate variance in breeding values == VA
## `var()` uses "n-1" in the denominator. Sorensen et al. 2001, eqn. 8 uses "n"
## `genVar` for GENERATION-by-line Variance, NOT specifically genetic variance
genVarLine <- function(gen, indices, datfxSub, iline, lindices){
  apply(datfxSub[which(indices == gen & lindices == iline), ], MARGIN = 2,
    FUN = function(x) { (1/length(x)) * sum((x^2) - mean(x)^2)} )
}
# Covariance among breeding values for 3-D ARRAY, for a single line
genCovLine <- function(gen, indices, datArr, iline, lindices){
  apply(datArr[, which(indices == gen & lindices == iline), ], MARGIN = 3,
    FUN = function(x) { (1/ncol(x)) * sum( (x[1, ] - mean(x[1, ])) *
                                           (x[2, ] - mean(x[2, ])) )} )
}

#extract genetic correlations for each line in each generation
#CONTROL MICE
#CONTROL MICE
#CONTROL MICE
#CONTROL MICE
#CONTROL MICE
#CONTROL MICE
#CONTROL MICE
LINES<-c(1,2,4,5)
  # Setup random effects posteriors into array (datfx)
  ## Get names from Z/random effects design matrix so can split up BLUPs
  Zcnms <- C_RPM_INT_mnth$Z@Dimnames[[2L]]
    ZcnmsLst <- strsplit(Zcnms, split = ".", fixed = TRUE)
    ZcnmTrt <- sapply(ZcnmsLst, FUN = "[[", 1)
    ZcnmTrm <- sapply(ZcnmsLst, FUN = "[[", 2)

  ## get names of random effect terms from formula
  rfmla <- C_RPM_INT_mnth$Random$formula
    if(length(rfmla) > 1) rfmla <- rfmla[[2]]
  rtrms <- strsplit(deparse(rfmla), split = " \\+ ")[[1]]
  rfacNms <- sapply(strsplit(rtrms, "\\:"), FUN = tail, 1)  #<-- XXX ASSUME only ever 1 ":" used


  nffx <- C_RPM_INT_mnth$Fixed$nfl  #<-- number of fixed effects (first n solutions in `Sol`)

  ## Z is C0datf diagonal for traits - so has 0s for all rows involved with other trait
  ### Create an index to avoid Z row zeroes.
  traitZrowInd <- rbind(seq(from = 1, by = C_RPM_INT_mnth$Z@Dim[1L] / length(traits),
        length.out = length(traits)),
    seq(from = C_RPM_INT_mnth$Z@Dim[1L] / length(traits), by = C_RPM_INT_mnth$Z@Dim[1L] / length(traits),
        length.out = length(traits)))

  ncomp <- 1  #keep track of how many components have been worked through
  subZlst <- vector("list", sum(C_RPM_INT_mnth$Random$nfl))
  datfx <- array(NA, dim = c(sum(C_RPM_INT_mnth$Random$nfl),
		nrow(C_RPM_INT_mnth$X) / length(traits),
		nrow(C_RPM_INT_mnth$Sol)))
    dimnames(datfx) <- list(rep(NA, dim(datfx)[[1L]]), NULL, NULL)

  for(g in 1:length(rfacNms)){
   cat("\n g=", g, "of", length(rfacNms))
    gfacNm <- rfacNms[g]
    for(t in 1:length(traits)){  # loop through TRAITs of gth factor/term in model
     cat("\n\t t=", t, "of", length(traits))
      tname <- traits[t]
      gtZind <- which(ZcnmTrt == paste0("trait", tname) & ZcnmTrm == gfacNm)
      tmpZsub <- C_RPM_INT_mnth$Z[, gtZind]
        subZlst[[ncomp]] <- tmpZsub
      tmpDatfx <- sapply(seq(nrow(C_RPM_INT_mnth$Sol)),
          FUN = function(i) (tmpZsub %*% C_RPM_INT_mnth$Sol[i, nffx + gtZind])@x)
      # Z is filled with zeroes for rows of other trait(s) so above gives BLUPs=0
      ## First check, then drop 0.0s
      nott <- seq(length(traits))[-t] #<--
      stopifnot(all(tmpDatfx[c(sapply(nott, FUN = function(x){
            seq(traitZrowInd[1,x], traitZrowInd[2,x], by=1)})), ] == 0.0))
      datfx[ncomp, , ] <- tmpDatfx[seq(traitZrowInd[1,t], traitZrowInd[2,t], by=1), ]
        dimnames(datfx)[[1L]][ncomp] <- paste(tname, gfacNm, sep = ".")
      ncomp <- ncomp + 1
    }  # end for t
  }  # end for g
  cat("\n")

########## make check
estPts <- sort(unique(C0datf$GEN)) #seq(0, max(C0DATF$GEN), 1)
# (check to make sure order in MCMCglmm model is same as in C0DATF)
## pull out GEN for first trait from X/fixed effects design matrix and subtract
###  from C0DATF$GEN <--> expect no difference/0s (What done in code below)

# BEGIN check
  # Make index of X col. names for a single trait that correspond to generation in data
  ## (exclude first generation, as this lumped in with intercept)
  tmpXcnmsTrtGenFac <- match(paste0("trait", traits[1], ":GENfac", estPts[-1]),
      C_RPM_INT_mnth$X@Dimnames[[2]])
  # Grab first trait rows and columns corresponding to generation factor
  t1GenFacX <- C_RPM_INT_mnth$X[seq(1, C_RPM_INT_mnth$Z@Dim[1] / length(traits), 1), tmpXcnmsTrtGenFac]
  # Column with non-zero indicates generation, so if multiply this sub-X matrix
  ## by vector of GENfac as integer will get generation for every individual from X
  ## Gives 0 to GEN0 because no non-zero entry in this sub-X for such rows
  ### (intercept taken out)
  gensFromX <- t1GenFacX %*% matrix(estPts[-1], ncol = 1) # XXX excludes GEN=0/intercept
stopifnot((gensFromX@x - C0datf$GEN) == 0)
# END check

##XXX Only do for additive genetic (co)variances and NOT dam (co)variances
postVarLines <- lapply(1:2, FUN = function(t){
      tmpArr <- array(NA, dim = c(dim(datfx)[3], length(estPts), length(LINES)),
        dimnames = list(NULL, as.character(estPts), as.character(LINES)))
      for(l in 1:length(LINES)){
        tmpArr[, , l] <- structure(sapply(estPts,
          FUN = genVarLine, indices = C0datf$GEN, datfxSub = datfx[t, , ],
                iline = LINES[l], lindices = C0datf$line))
      }
      return(structure(tmpArr, mcpar = attr(C_RPM_INT_mnth$Sol, "mcpar")))})
  names(postVarLines) <- dimnames(datfx)[[1L]][1:2]

postCovLines <- array(NA, dim = c(dim(datfx)[3], length(estPts), length(LINES)),
  dimnames = list(NULL, as.character(estPts), as.character(LINES)))
for(l in 1:length(LINES)){
  postCovLines[, , l] <- structure(sapply(estPts,
                  FUN = genCovLine, indices = C0datf$GEN, datArr = datfx[1:2, , ],
                        iline = LINES[l], lindices = C0datf$line))
}  #<-- end for l
  attr(postCovLines, "mcpar") <- attr(C_RPM_INT_mnth$Sol, "mcpar")
  postCovLines <- list(postCovLines)
    names(postCovLines) <- paste(paste0(traits, collapse = ":"), "animal",
      sep = ".")

## calculate correlation over posterior distribution
postRa <- list(rA = postCovLines[[1L]] /sqrt(postVarLines[["RPM56l.animal"]] * postVarLines[["INT56l.animal"]]))
POST.C.LINES <- c(postVarLines, postCovLines, postRa)

#SELECTED MICE
#SELECTED MICE
#SELECTED MICE
#SELECTED MICE
#SELECTED MICE
#SELECTED MICE
#SELECTED MICE
#SELECTED MICE
LINES <- c(3,6,7,8)
      # Setup random effects posteriors into array (datfx)
  ## Get names from Z/random effects design matrix so can split up BLUPs
  Zcnms <- HR_RPM_INT_mnth$Z@Dimnames[[2L]]
    ZcnmsLst <- strsplit(Zcnms, split = ".", fixed = TRUE)
    ZcnmTrt <- sapply(ZcnmsLst, FUN = "[[", 1)
    ZcnmTrm <- sapply(ZcnmsLst, FUN = "[[", 2)

  ## get names of random effect terms from formula
  rfmla <- HR_RPM_INT_mnth$Random$formula
    if(length(rfmla) > 1) rfmla <- rfmla[[2]]
  rtrms <- strsplit(deparse(rfmla), split = " \\+ ")[[1]]
  rfacNms <- sapply(strsplit(rtrms, "\\:"), FUN = tail, 1)  #<-- XXX ASSUME only ever 1 ":" used


  nffx <- HR_RPM_INT_mnth$Fixed$nfl  #<-- number of fixed effects (first n solutions in `Sol`)

  ## Z is C0datf diagonal for traits - so has 0s for all rows involved with other trait
  ### Create an index to avoid Z row zeroes.
  traitZrowInd <- rbind(seq(from = 1, by = HR_RPM_INT_mnth$Z@Dim[1L] / length(traits),
        length.out = length(traits)),
    seq(from = HR_RPM_INT_mnth$Z@Dim[1L] / length(traits), by = HR_RPM_INT_mnth$Z@Dim[1L] / length(traits),
        length.out = length(traits)))

  ncomp <- 1  #keep track of how many components have been worked through
  subZlst <- vector("list", sum(HR_RPM_INT_mnth$Random$nfl))
  datfx <- array(NA, dim = c(sum(HR_RPM_INT_mnth$Random$nfl),
		nrow(HR_RPM_INT_mnth$X) / length(traits),
		nrow(HR_RPM_INT_mnth$Sol)))
    dimnames(datfx) <- list(rep(NA, dim(datfx)[[1L]]), NULL, NULL)

  for(g in 1:length(rfacNms)){
   cat("\n g=", g, "of", length(rfacNms))
    gfacNm <- rfacNms[g]
    for(t in 1:length(traits)){  # loop through TRAITs of gth factor/term in model
     cat("\n\t t=", t, "of", length(traits))
      tname <- traits[t]
      gtZind <- which(ZcnmTrt == paste0("trait", tname) & ZcnmTrm == gfacNm)
      tmpZsub <- HR_RPM_INT_mnth$Z[, gtZind]
        subZlst[[ncomp]] <- tmpZsub
      tmpDatfx <- sapply(seq(nrow(HR_RPM_INT_mnth$Sol)),
          FUN = function(i) (tmpZsub %*% HR_RPM_INT_mnth$Sol[i, nffx + gtZind])@x)
      # Z is filled with zeroes for rows of other trait(s) so above gives BLUPs=0
      ## First check, then drop 0.0s
      nott <- seq(length(traits))[-t] #<--
      stopifnot(all(tmpDatfx[c(sapply(nott, FUN = function(x){
            seq(traitZrowInd[1,x], traitZrowInd[2,x], by=1)})), ] == 0.0))
      datfx[ncomp, , ] <- tmpDatfx[seq(traitZrowInd[1,t], traitZrowInd[2,t], by=1), ]
        dimnames(datfx)[[1L]][ncomp] <- paste(tname, gfacNm, sep = ".")
      ncomp <- ncomp + 1
    }  # end for t
  }  # end for g
  cat("\n")

########## make check
estPts <- sort(unique(HR1datf$GEN)) #seq(0, max(C0DATF$GEN), 1)
# (check to make sure order in MCMCglmm model is same as in C0DATF)
## pull out GEN for first trait from X/fixed effects design matrix and subtract
###  from C0DATF$GEN <--> expect no difference/0s (What done in code below)

# BEGIN check
  # Make index of X col. names for a single trait that correspond to generation in data
  ## (exclude first generation, as this lumped in with intercept)
  tmpXcnmsTrtGenFac <- match(paste0("trait", traits[1], ":GENfac", estPts[-1]),
      HR_RPM_INT_mnth$X@Dimnames[[2]])
  # Grab first trait rows and columns corresponding to generation factor
  t1GenFacX <- HR_RPM_INT_mnth$X[seq(1, HR_RPM_INT_mnth$Z@Dim[1] / length(traits), 1), tmpXcnmsTrtGenFac]
  # Column with non-zero indicates generation, so if multiply this sub-X matrix
  ## by vector of GENfac as integer will get generation for every individual from X
  ## Gives 0 to GEN0 because no non-zero entry in this sub-X for such rows
  ### (intercept taken out)
  gensFromX <- t1GenFacX %*% matrix(estPts[-1], ncol = 1) # XXX excludes GEN=0/intercept
stopifnot((gensFromX@x - HR1datf$GEN) == 0)
# END check

# LINE-specific (Co)Variances
##XXX Only do for additive genetic (co)variances and NOT dam (co)variances
postVarLines <- lapply(1:2, FUN = function(t){
      tmpArr <- array(NA, dim = c(dim(datfx)[3], length(estPts), length(LINES)),
        dimnames = list(NULL, as.character(estPts), as.character(LINES)))
      for(l in 1:length(LINES)){
        tmpArr[, , l] <- structure(sapply(estPts,
          FUN = genVarLine, indices = HR1datf$GEN, datfxSub = datfx[t, , ],
                iline = LINES[l], lindices = HR1datf$line))
      }
      return(structure(tmpArr, mcpar = attr(HR_RPM_INT_mnth$Sol, "mcpar")))})
  names(postVarLines) <- dimnames(datfx)[[1L]][1:2]

postCovLines <- array(NA, dim = c(dim(datfx)[3], length(estPts), length(LINES)),
  dimnames = list(NULL, as.character(estPts), as.character(LINES)))
for(l in 1:length(LINES)){
  postCovLines[, , l] <- structure(sapply(estPts,
                  FUN = genCovLine, indices = HR1datf$GEN, datArr = datfx[1:2, , ],
                        iline = LINES[l], lindices = HR1datf$line))
}  #<-- end for l
  attr(postCovLines, "mcpar") <- attr(HR_RPM_INT_mnth$Sol, "mcpar")
  postCovLines <- list(postCovLines)
    names(postCovLines) <- paste(paste0(traits, collapse = ":"), "animal",
      sep = ".")

## calculate correlation over posterior distribution
postRa <- list(rA = postCovLines[[1L]] /sqrt(postVarLines[["RPM56l.animal"]] * postVarLines[["INT56l.animal"]]))
POST.S.LINES <- c(postVarLines, postCovLines, postRa)

save("POST.C.LINES", "POST.S.LINES", file = "QG-AP-test_covar_POST.LINES.RData")
rm(list=ls())  # delete all objects, to start fresh ion section below



#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################

library(MCMCglmm)

GEN.list<-c(0:31,36:51,53:62,65,66,68:78)#list of GENs without data
#load(file = "covar_POST.LINES.RData")
load(file = "QG-AP-test_covar_POST.LINES.RData")

rA.C.LINES<-POST.C.LINES$rA
rA.S.LINES<-POST.S.LINES$rA
rA.LINE.1<-rA.C.LINES[,,1]
rA.LINE.2<-rA.C.LINES[,,2]
rA.LINE.4<-rA.C.LINES[,,3]
rA.LINE.5<-rA.C.LINES[,,4]
rA.LINE.3<-rA.S.LINES[,,1]
rA.LINE.6<-rA.S.LINES[,,2]
rA.LINE.7<-rA.S.LINES[,,3]
rA.LINE.8<-rA.S.LINES[,,4]
#
Ra.1<-data.frame (GEN=GEN.list,r=NA)
Ra.2<-data.frame (GEN=GEN.list,r=NA)
Ra.3<-data.frame (GEN=GEN.list,r=NA)
Ra.4<-data.frame (GEN=GEN.list,r=NA)
Ra.5<-data.frame (GEN=GEN.list,r=NA)
Ra.6<-data.frame (GEN=GEN.list,r=NA)
Ra.7<-data.frame (GEN=GEN.list,r=NA)
Ra.8<-data.frame (GEN=GEN.list,r=NA)
#
Ra.1$r<-posterior.mode(as.mcmc(rA.LINE.1))
Ra.2$r<-posterior.mode(as.mcmc(rA.LINE.2))
Ra.3$r<-posterior.mode(as.mcmc(rA.LINE.3))
Ra.4$r<-posterior.mode(as.mcmc(rA.LINE.4))
Ra.5$r<-posterior.mode(as.mcmc(rA.LINE.5))
Ra.6$r<-posterior.mode(as.mcmc(rA.LINE.6))
Ra.7$r<-posterior.mode(as.mcmc(rA.LINE.7))
Ra.8$r<-posterior.mode(as.mcmc(rA.LINE.8))

Ra.C<-data.frame (GEN=GEN.list,r=NA, lower=NA,upper=NA)
Ra.S<-data.frame (GEN=GEN.list,r=NA, lower=NA,upper=NA)
acrossLinePost.C <- as.mcmc(apply(POST.C.LINES$rA, MARGIN = 2,FUN = function(v){ c(v)}))
acrossLinePost.S <- as.mcmc(apply(POST.S.LINES$rA, MARGIN = 2,FUN = function(v){ c(v)}))
Ra.C$r<-posterior.mode(acrossLinePost.C)
Ra.S$r<-posterior.mode(acrossLinePost.S)
Ra.C[,3:4]<-HPDinterval(acrossLinePost.C)
Ra.S[,3:4]<-HPDinterval(acrossLinePost.S)
new_rows<- data.frame(GEN=c(32,33,34,35,52,63,64,67), r=NA,lower=NA,upper=NA)
Ra.C<-rbind(Ra.C,new_rows)
Ra.S<-rbind(Ra.S,new_rows)
Ra.C$GEN<-Ra.C$GEN-0.1
Ra.S$GEN<-Ra.S$GEN+0.1
#
Ra.C<-Ra.C[order(Ra.C$GEN),]
Ra.S<-Ra.S[order(Ra.S$GEN),]
#
new_rows<- data.frame(GEN=c(32,33,34,35,52,63,64,67),r=NA)
Ra.1<-rbind(Ra.1,new_rows)
Ra.2<-rbind(Ra.2,new_rows)
Ra.3<-rbind(Ra.3,new_rows)
Ra.4<-rbind(Ra.4,new_rows)
Ra.5<-rbind(Ra.5,new_rows)
Ra.6<-rbind(Ra.6,new_rows)
Ra.7<-rbind(Ra.7,new_rows)
Ra.8<-rbind(Ra.8,new_rows)
#
Ra.1<-Ra.1[order(Ra.1$GEN),]
Ra.2<-Ra.2[order(Ra.2$GEN),]
Ra.3<-Ra.3[order(Ra.3$GEN),]
Ra.4<-Ra.4[order(Ra.4$GEN),]
Ra.5<-Ra.5[order(Ra.5$GEN),]
Ra.6<-Ra.6[order(Ra.6$GEN),]
Ra.7<-Ra.7[order(Ra.7$GEN),]
Ra.8<-Ra.8[order(Ra.8$GEN),]

# get breeding values for running speed and distance
load(file="QG-AP-test_Rdata.RData")
nrow(C0datf)     #N = 10,878
nrow(HR1datf)    #N = 25,124

load(file="QG-AP-test_bivariate_models.RData")
# Speed and Duration models
#extract breeding values
SOL.C<-C_RPM_INT_mnth$Sol
SOL.S<-HR_RPM_INT_mnth$Sol
colnames(SOL.S)
#
MODE.C<-posterior.mode(SOL.C)
CI.C<-data.frame(HPDinterval(SOL.C))
TEMP.C<-data.frame(names(MODE.C),MODE.C,CI.C)
BLUPS.C <- TEMP.C[grep("animal", TEMP.C$names.MODE.C.), ]
nrow(BLUPS.C)
colnames(BLUPS.C)<-c("name","BLUP","lower","upper")
head(BLUPS.C,25)
BLUPS.C$TRAIT<-substr(BLUPS.C$name,6,8)
x<-strsplit(as.character(BLUPS.C$name),".animal.")
x2<-data.frame(Reduce(rbind, x))
BLUPS.C$animal<-x2$X2
head(BLUPS.C)
rownames(BLUPS.C)<-NULL
BLUPS.C$name<-NULL
BLUPS.C.wide<-data.frame(reshape(BLUPS.C, v.names = c("BLUP","lower","upper"), idvar = "animal",timevar = "TRAIT", direction = "wide"))
nrow(BLUPS.C)
nrow(BLUPS.C.wide)
DATA.C<-merge(C0datf,BLUPS.C.wide,by="animal",all.x = T)
#
MODE.S<-posterior.mode(SOL.S)
CI.S<-data.frame(HPDinterval(SOL.S))
TEMP.S<-data.frame(names(MODE.S),MODE.S,CI.S)
BLUPS.S <- TEMP.S[grep("animal", TEMP.S$names.MODE.S.), ]
nrow(BLUPS.S)
colnames(BLUPS.S)<-c("name","BLUP","lower","upper")
head(BLUPS.S,25)
BLUPS.S$TRAIT<-substr(BLUPS.S$name,6,8)
x<-strsplit(as.character(BLUPS.S$name),".animal.")
x2<-data.frame(Reduce(rbind, x))
BLUPS.S$animal<-x2$X2
head(BLUPS.S)
rownames(BLUPS.S)<-NULL
BLUPS.S$name<-NULL
BLUPS.S.wide<-data.frame(reshape(BLUPS.S, v.names = c("BLUP","lower","upper"), idvar = "animal",timevar = "TRAIT", direction = "wide"))
nrow(BLUPS.S)
nrow(BLUPS.S.wide)
DATA.S<-merge(HR1datf,BLUPS.S.wide,by="animal",all.x = T)
head(DATA.S)
DATA<-rbind(DATA.C,DATA.S)
DATA<-droplevels(DATA)
summary(DATA)
#BLUPs RPM trajectories over generations
LINE.C<-aggregate(cbind(BLUP.RPM,BLUP.INT)~GEN+line,data=DATA.C,FUN=mean)
LINE.S<-aggregate(cbind(BLUP.RPM,BLUP.INT)~GEN+line,data=DATA.S,FUN=mean)


#make figure 3 BREEDING VALUES and rA
dev.new(width=9,height=6, units = "cm")
par(mfrow=c(2,2),las=1, oma=c(4,4,1,1),mar=c(1,1,1,1))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
layout.show(3)
#
  plot(BLUP.RPM~GEN,subset(LINE.C,line==1),type="l", col="blue",ylim=c(-0.1,0.6),ylab="",xlab="n",xaxt="n")
points(BLUP.RPM~GEN,subset(LINE.C,line==2),type="l", col="blue")
points(BLUP.RPM~GEN,subset(LINE.C,line==4),type="l", col="blue")
points(BLUP.RPM~GEN,subset(LINE.C,line==5),type="l", col="blue")
points(BLUP.RPM~GEN,subset(LINE.S,line==3),type="l", col="red")
points(BLUP.RPM~GEN,subset(LINE.S,line==6),type="l", col="red")
points(BLUP.RPM~GEN,subset(LINE.S,line==7),type="l", col="red")
points(BLUP.RPM~GEN,subset(LINE.S,line==8),type="l", col="red")
abline(h=0,lty=3)
mtext("Genetic merit (breeding values)",side=2,las=3, line=3)
mtext("A running speed",side=3,adj=0.015,line=-1.5)
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
#BLUPs INT trajectories over generations
  plot(BLUP.INT~GEN,subset(LINE.C,line==1),type="l", col="blue",ylim=c(-0.1,0.6),ylab="",xlab="",yaxt="n",xaxt="n")
points(BLUP.INT~GEN,subset(LINE.C,line==2),type="l", col="blue")
points(BLUP.INT~GEN,subset(LINE.C,line==4),type="l", col="blue")
points(BLUP.INT~GEN,subset(LINE.C,line==5),type="l", col="blue")
points(BLUP.INT~GEN,subset(LINE.S,line==3),type="l", col="red")
points(BLUP.INT~GEN,subset(LINE.S,line==6),type="l", col="red")
points(BLUP.INT~GEN,subset(LINE.S,line==7),type="l", col="red")
points(BLUP.INT~GEN,subset(LINE.S,line==8),type="l", col="red")
abline(h=0,lty=3)
axis(2,labels=F)
mtext("B running duration",side=3,adj=0.015,line=-1.5)
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
#
  plot(r~GEN, Ra.C, col="blue",pch=16, cex=1,type="o", ylim=c(-0.35,1),xlim=c(1,77),xlab="",xaxt="n",ylab="")
points(r~GEN, Ra.S, col="red" ,pch=17, cex=1,type="o")
arrows(x0=Ra.C$GEN,y0=Ra.C$lower,x1=Ra.C$GEN,y1=Ra.C$upper,col="black",code=3,angle=90,length=0)
arrows(x0=Ra.S$GEN,y0=Ra.S$lower,x1=Ra.S$GEN,y1=Ra.S$upper,col="black" ,code=3,angle=90,length=0)
points(r~GEN, Ra.1, col=rgb(0,0,1,0.75), cex=1,type="l")
points(r~GEN, Ra.2, col=rgb(0,0,1,0.75), cex=1,type="l")
points(r~GEN, Ra.3, col=rgb(1,0,0,0.75), cex=1,type="l")
points(r~GEN, Ra.4, col=rgb(0,0,1,0.75), cex=1,type="l")
points(r~GEN, Ra.5, col=rgb(0,0,1,0.75), cex=1,type="l")
points(r~GEN, Ra.6, col=rgb(1,0,0,0.75), cex=1,type="l")
points(r~GEN, Ra.7, col=rgb(1,0,0,0.75), cex=1,type="l")
points(r~GEN, Ra.8, col=rgb(1,0,0,0.75), cex=1,type="l")
points(r~GEN, Ra.C, col="black" ,pch=16, cex=1,type="o")
points(r~GEN, Ra.S, col="black" ,pch=17, cex=1,type="o")
axis(1,at=0:78, labels=F)
axis(1, at = seq(0, 78, 6),padj=-0.8)
abline(h=0, lty=2)
mtext(expression(paste(italic(r)[A]," (±se)")), side=2, line=2.5, las=3,cex=2)
mtext("Generation",side=1, line=2.5,cex=1.5)
#
#mtext("Wisconsin generations", line=0.5, adj=0.13, side=3, cex=1)
#mtext("Riverside generations", line=0.5, adj=0.75, side=3, cex=1)
legend(0,1.1, c("control (average)","selected (average)"),          pch=c(16,17), col=c("black","black"),bty="n")
legend(60,1.1, c("replicate control lines","replicate selected lines"),lty=1, col=c(rgb(0,0,1,0.75),rgb(1,0,0,0.75)),bty="n")
arrows(26,1,26,0.9,length=0.1,lwd=3)
arrows(60,1,60,0.9,length=0.1,lwd=3)
mtext("C",side=3,adj=0.015,line=-1.5)

#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#to make figure 4, you need objects DATA.C and DATA.S produced in section 3D above
dev.new(width=9,height=5, units = "cm")
par(mfrow=c(2,4),las=1, oma=c(5,5,0,0),mar=c(1,1,2,2))
layout(matrix(c(1,1,2,3,1,1,4,5), 2, 4, byrow = TRUE))
layout.show(5)
#
  plot(BLUP.RPM~BLUP.INT,DATA,col="white",ylab="",xlab="",cex.lab=2)
points(BLUP.RPM~BLUP.INT,DATA.C,col=rgb(0,0,1,0.1),pch=16,cex=1.5)
points(BLUP.RPM~BLUP.INT,DATA.S,col=rgb(1,0,0,0.1),pch=17,cex=1.5)
#add trajectories
points(BLUP.RPM~BLUP.INT,subset(LINE.C,line==1),type="l",lwd=2, col="grey")
points(BLUP.RPM~BLUP.INT,subset(LINE.C,line==2),type="l",lwd=2, col="grey")
points(BLUP.RPM~BLUP.INT,subset(LINE.C,line==4),type="l",lwd=2, col="grey")
points(BLUP.RPM~BLUP.INT,subset(LINE.C,line==5),type="l",lwd=2, col="grey")
points(BLUP.RPM~BLUP.INT,subset(LINE.S,line==3),type="l",lwd=2)
points(BLUP.RPM~BLUP.INT,subset(LINE.S,line==6),type="l",lwd=2)
points(BLUP.RPM~BLUP.INT,subset(LINE.S,line==7),type="l",lwd=2)
points(BLUP.RPM~BLUP.INT,subset(LINE.S,line==8),type="l",lwd=2)
#
#add ellipses as specific points in time
library(heplots)
covEllipses(DATA[which(DATA$GEN==15&DATA$line==3),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S3",cex=1.2)
covEllipses(DATA[which(DATA$GEN==15&DATA$line==6),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S6",cex=1.2)
covEllipses(DATA[which(DATA$GEN==15&DATA$line==7),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S7",cex=1.2)
covEllipses(DATA[which(DATA$GEN==15&DATA$line==8),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S8",cex=1.2)
#
covEllipses(DATA[which(DATA$GEN==40&DATA$line==3),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S3",cex=1.2)
covEllipses(DATA[which(DATA$GEN==40&DATA$line==6),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S6",cex=1.2)
covEllipses(DATA[which(DATA$GEN==40&DATA$line==7),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S7",cex=1.2)
covEllipses(DATA[which(DATA$GEN==40&DATA$line==8),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S8",cex=1.2)
#
covEllipses(DATA[which(DATA$GEN==78&DATA$line==1),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="C1",cex=1.2)
covEllipses(DATA[which(DATA$GEN==78&DATA$line==2),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="C2",cex=1.2)
covEllipses(DATA[which(DATA$GEN==78&DATA$line==4),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="C4",cex=1.2)
covEllipses(DATA[which(DATA$GEN==78&DATA$line==5),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="C5",cex=1.2)
#
covEllipses(DATA[which(DATA$GEN==78&DATA$line==3),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S3",cex=1.2)
covEllipses(DATA[which(DATA$GEN==78&DATA$line==6),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S6",cex=1.2)
covEllipses(DATA[which(DATA$GEN==78&DATA$line==7),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S7",cex=1.2)
covEllipses(DATA[which(DATA$GEN==78&DATA$line==8),c("BLUP.INT","BLUP.RPM")],add=T,col="white",center=F,center.pch="",labels="S8",cex=1.2)
#
legend( x="bottomright",
        legend=c("Control mice","Control G trajectories","Selected mice","Selected G trajectories"),
        col=c("blue","grey","red","black"), lwd=2, lty=c(NA,1,NA,1),
        pch=c(16,NA,17,NA), merge=T, bty="n")
text(-0.05,0.35,"generation 15")
arrows(-0.05,0.33,0.025,0.31,length=0.1)
text(0.24,0.25,"generation 40")
arrows(0.24,0.26,0.215,0.32,length=0.1)
text(0.085,0.625,"generation 78")
arrows(0.1,0.61,0.13,0.59,length=0.1)
mtext("Genetic merit for running speed",   side=2,line=2,outer=T,las=3,cex=1.5)
mtext("Genetic merit for running duration",side=1,line=3,outer=T,las=1,cex=1.5)
mtext("A) 78 generations",side=3,adj=0.015,line=-1.5)

#show breeding values for 4 generation
LAST.1<-subset(DATA.S,GEN==1)
LAST.2<-subset(DATA.S,GEN==15)
LAST.3<-subset(DATA.S,GEN==40)
LAST.4<-subset(DATA.S,GEN==78)
ylim.1<-c(min(LAST.1$lower.RPM),max(LAST.1$upper.RPM))
ylim.2<-c(min(LAST.2$lower.RPM),max(LAST.2$upper.RPM))
ylim.3<-c(min(LAST.3$lower.RPM),max(LAST.3$upper.RPM))
ylim.4<-c(min(LAST.4$lower.RPM),max(LAST.4$upper.RPM))
xlim.1<-c(min(LAST.1$lower.INT),max(LAST.1$upper.INT))
xlim.2<-c(min(LAST.2$lower.INT),max(LAST.2$upper.INT))
xlim.3<-c(min(LAST.3$lower.INT),max(LAST.3$upper.INT))
xlim.4<-c(min(LAST.4$lower.INT),max(LAST.4$upper.INT))
#
  plot(BLUP.RPM~BLUP.INT,LAST.1,ylim=ylim.1,xlim=xlim.1,ylab="",xlab="")
 arrows(LAST.1$BLUP.INT, LAST.1$lower.RPM,LAST.1$BLUP.INT, LAST.1$upper.RPM,col=rgb(0,0,0,0.05),code=3,angle=90,length=0)
 arrows(LAST.1$lower.INT,LAST.1$BLUP.RPM, LAST.1$upper.INT,LAST.1$BLUP.RPM, col=rgb(0,0,0,0.05),code=3,angle=90,length=0)
points(BLUP.RPM~BLUP.INT,LAST.1,col=line,pch=16,cex=1.5)
mtext("B) generation 1",side=3,adj=0.015,line=-1.5, cex=0.8)
#
  plot(BLUP.RPM~BLUP.INT,LAST.2,ylim=ylim.2,xlim=xlim.2,ylab="",xlab="")
 arrows(LAST.2$BLUP.INT, LAST.2$lower.RPM,LAST.2$BLUP.INT, LAST.2$upper.RPM,col=rgb(0,0,0,0.05),code=3,angle=90,length=0)
 arrows(LAST.2$lower.INT,LAST.2$BLUP.RPM, LAST.2$upper.INT,LAST.2$BLUP.RPM, col=rgb(0,0,0,0.05),code=3,angle=90,length=0)
points(BLUP.RPM~BLUP.INT,LAST.2,col=line,pch=16,cex=1.5)
mtext("C) generation 15",side=3,adj=0.015,line=-1.5, cex=0.8)
#
  plot(BLUP.RPM~BLUP.INT,LAST.3,ylim=ylim.3,xlim=xlim.3,ylab="",xlab="")
 arrows(LAST.3$BLUP.INT, LAST.3$lower.RPM,LAST.3$BLUP.INT, LAST.3$upper.RPM,col=rgb(0,0,0,0.05),code=3,angle=90,length=0)
 arrows(LAST.3$lower.INT,LAST.3$BLUP.RPM, LAST.3$upper.INT,LAST.3$BLUP.RPM, col=rgb(0,0,0,0.05),code=3,angle=90,length=0)
points(BLUP.RPM~BLUP.INT,LAST.3,col=line,pch=16,cex=1.5)
mtext("D) generation 40",side=3,adj=0.015,line=-1.5, cex=0.8)
#
  plot(BLUP.RPM~BLUP.INT,LAST.4,ylim=ylim.4,xlim=xlim.4,ylab="",xlab="")
 arrows(LAST.4$BLUP.INT, LAST.4$lower.RPM,LAST.4$BLUP.INT, LAST.4$upper.RPM,col=rgb(0,0,0,0.05),code=3,angle=90,length=0)
 arrows(LAST.4$lower.INT,LAST.4$BLUP.RPM, LAST.4$upper.INT,LAST.4$BLUP.RPM, col=rgb(0,0,0,0.05),code=3,angle=90,length=0)
points(BLUP.RPM~BLUP.INT,LAST.4,col=line,pch=16,cex=1.5)
mtext("E) generation 78",side=3,adj=0.015,line=-1.5, cex=0.8)
#
legend( x="topright",legend=c("line S3","line S6","line S7","line S8"),col=c("black","red","blue","green"),pch=16,bty="n")


cor.test(DATA$BLUP.INT[which(DATA$GEN==1&DATA$linetype==1)] ,DATA$BLUP.RPM[which(DATA$GEN==1&DATA$linetype==1)])
cor.test(DATA$BLUP.INT[which(DATA$GEN==15&DATA$linetype==1)],DATA$BLUP.RPM[which(DATA$GEN==15&DATA$linetype==1)])
cor.test(DATA$BLUP.INT[which(DATA$GEN==40&DATA$linetype==1)],DATA$BLUP.RPM[which(DATA$GEN==40&DATA$linetype==1)])
cor.test(DATA$BLUP.INT[which(DATA$GEN==78&DATA$linetype==1)],DATA$BLUP.RPM[which(DATA$GEN==78&DATA$linetype==1)])










#.######################## section #5 ANIMATED GIF VERSION FIGURE S1 ###################
#.######################## section #5 ANIMATED GIF VERSION FIGURE S1 ###################
#.######################## section #5 ANIMATED GIF VERSION FIGURE S1 ###################
#.######################## section #5 ANIMATED GIF VERSION FIGURE S1 ###################
#.######################## section #5 ANIMATED GIF VERSION FIGURE S1 ###################
#.######################## section #5 ANIMATED GIF VERSION FIGURE S1 ###################
#.######################## section #5 ANIMATED GIF VERSION FIGURE S1 ###################
#.######################## section #5 ANIMATED GIF VERSION FIGURE S1 ###################
library(transformr)
library(gganimate)
library(ggplot2)

theme_set(theme_bw())
   theme_update(text = element_text(size=22),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_blank())

PLOT<-ggplot(DATA, aes(BLUP.INT, BLUP.RPM, colour = line)) +
   geom_point(alpha = 0.5, show.legend = F) +
   scale_colour_manual(values=c("#000099","#000099","#000099","#000099", "#990000", "#990000", "#990000", "#990000")) +
   stat_ellipse()+
   labs(title = 'Generation: {frame_time}', x = 'Breeding value for running duration', y = 'Breeding value for running speed')+
   #theme_bw(base_size = 22) +
   transition_time(GEN) +
   ease_aes('linear')
PLOT

animate(PLOT)
anim_save('Figure_S1.gif')

rm(list=ls())
#.######################## section 6 MAKE Table S1 ###############################
#.######################## section 6 MAKE Table S1 ###############################
#.######################## section 6 MAKE Table S1 ###############################
#.######################## section 6 MAKE Table S1 ###############################
#.######################## section 6 MAKE Table S1 ###############################
#.######################## section 6 MAKE Table S1 ###############################
load(file="QG-AP-test_bivariate_models.RData")
#make table S1 showing the estimates for the MCMCglmm models
Va.C<-rbind(c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.animal"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.animal"])),
            c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.animal"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.animal"])),
            c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.animal"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.animal"])))
Va.S<-rbind(c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.animal"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.animal"])),
            c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.animal"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.animal"])),
            c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.animal"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.animal"])))
Vc.C<-rbind(c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.damid"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.damid"])),
            c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.damid"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.damid"])),
            c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.damid"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.damid"])))
Vc.S<-rbind(c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.damid"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.damid"])),
            c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.damid"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.damid"])),
            c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.damid"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.damid"])))
Ve.C<-rbind(c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.units"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.units"])),
            c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.units"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.units"])),
            c(posterior.mode(C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.units"]),HPDinterval(C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.units"])))
Ve.S<-rbind(c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.units"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.units"])),
            c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.units"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.units"])),
            c(posterior.mode(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.units"]),HPDinterval(HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.units"])))
#Make Table 1
rbind(cbind(Va.C,NA,Va.S),NA,
      cbind(Vc.C,NA,Vc.S),NA,
      cbind(Ve.C,NA,Ve.S),NA)


#.######################## section 7 MAKE Table S2 ###############################
#.######################## section 7 MAKE Table S2 ###############################
#.######################## section 7 MAKE Table S2 ###############################
#.######################## section 7 MAKE Table S2 ###############################
#.######################## section 7 MAKE Table S2 ###############################
#.######################## section 7 MAKE Table S2 ###############################
#.######################## section 7 MAKE Table S2 ###############################
#.######################## section 7 MAKE Table S2 ###############################
library(MCMCglmm)
load("QG-AP-test_bivariate_models.RData")
load("QG-AP-test_univariate_models.RData")
summary(C_RUN_mnth)
summary(HR_RUN_mnth)

#make table Ss showing the estimates for RUN56 aversus calculated from the bivariate models
Va.C.RUN.est<-c(posterior.mode(C_RUN_mnth$VCV[,"animal"]), HPDinterval(C_RUN_mnth$VCV[,"animal"]))
Va.S.RUN.est<-c(posterior.mode(HR_RUN_mnth$VCV[,"animal"]),HPDinterval(HR_RUN_mnth$VCV[,"animal"]))
post.Va.C.RUN<-C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.animal"] +C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.animal"] +2*C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.animal"]
post.Va.S.RUN<-HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.animal"]+HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.animal"]+2*HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.animal"]
Va.C.RUN.calc<-c(posterior.mode(post.Va.C.RUN),HPDinterval(post.Va.C.RUN))
Va.S.RUN.calc<-c(posterior.mode(post.Va.S.RUN),HPDinterval(post.Va.S.RUN))
#
Vc.C.RUN.est<-c(posterior.mode(C_RUN_mnth$VCV[,"damid"]),HPDinterval(C_RUN_mnth$VCV[,"damid"]))
Vc.S.RUN.est<-c(posterior.mode(HR_RUN_mnth$VCV[,"damid"]),HPDinterval(HR_RUN_mnth$VCV[,"damid"]))
post.Vc.C.RUN<-C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.damid"] +C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.damid"] +2*C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.damid"]
post.Vc.S.RUN<-HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.damid"]+HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.damid"]+2*HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.damid"]
Vc.C.RUN.calc<-c(posterior.mode(post.Vc.C.RUN),HPDinterval(post.Vc.C.RUN))
Vc.S.RUN.calc<-c(posterior.mode(post.Vc.S.RUN),HPDinterval(post.Vc.S.RUN))
#
Ve.C.RUN.est<-c(posterior.mode(C_RUN_mnth$VCV[,"units"]),HPDinterval(C_RUN_mnth$VCV[,"units"]))
Ve.S.RUN.est<-c(posterior.mode(HR_RUN_mnth$VCV[,"units"]),HPDinterval(HR_RUN_mnth$VCV[,"units"]))
post.Ve.C.RUN<-C_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.units"]+C_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.units"]+2*C_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.units"]
post.Ve.S.RUN<-HR_RPM_INT_mnth$VCV[,"traitINT56l:traitINT56l.units"] +HR_RPM_INT_mnth$VCV[,"traitRPM56l:traitRPM56l.units"] +2*HR_RPM_INT_mnth$VCV[,"traitINT56l:traitRPM56l.units"]
Ve.C.RUN.calc<-c(posterior.mode(post.Ve.C.RUN),HPDinterval(post.Ve.C.RUN))
Ve.S.RUN.calc<-c(posterior.mode(post.Ve.S.RUN),HPDinterval(post.Ve.S.RUN))

#Make Table S2
C.RUN<-rbind(Va.C.RUN.est,Va.C.RUN.calc,NA,
             Vc.C.RUN.est,Vc.C.RUN.calc,NA,
             Ve.C.RUN.est,Ve.C.RUN.calc)
S.RUN<-rbind(Va.S.RUN.est,Va.S.RUN.calc,NA,
             Vc.S.RUN.est,Vc.S.RUN.calc,NA,
             Ve.S.RUN.est,Ve.S.RUN.calc)
cbind(C.RUN,NA,S.RUN)




################################################################################
################################################################################
################################################################################
################################################################################

# Collect R and package version information
R.version.string
packageVersion("fields")
packageVersion("nadiv")
packageVersion("MCMCglmm")
packageVersion("heplots")
packageVersion("ggplot2")
packageVersion("gganimate")
packageVersion("transformr")



#This code reproduces results in the paper titled:
#"How do trade-offs emerge across populations while being absent within? A replicated selection experiment suggests all four evolutionary mechanisms are at play"
#by Wolak, Wilson, Hiramatsu, Garland, and Careau

#For questions, please contact the corresponding author of the article:
#Vincent Careau, University of Ottawa, vcareau@uottawa.ca

#to run this code, you need the following files:
# "QG-mouse-trade-off_data.txt"
# "QG-mouse-trade-off_pedigree.txt""

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
library(fields)  #<-- needed for thin plate splines `Tps()`
dat       <- read.table("QG-mouse-trade-off_data.txt", header = TRUE)
GEN.lst <- c(0:31, 37:51, 53:62, 65, 66, 68:77) #<-- GENs with running data

#calculate SELECTION GRADIENTS ###########################
DATA.C.WIS <- subset(dat, linetype == "0" & UNI == "WIS")
DATA.S.WIS <- subset(dat, linetype == "1" & UNI == "WIS")
#calculate relative fitness using the number of pups produced
DATA.C.WIS$relFIT <- DATA.C.WIS$pups / mean(DATA.C.WIS$pups)
DATA.S.WIS$relFIT <- DATA.S.WIS$pups / mean(DATA.S.WIS$pups)
#fit linear model
SEL.C.WIS <- lm(relFIT ~ scale(RPM56) + scale(INT56) + GEN + 
    WHLSTAGE + sex + Fcoeff + line,
  data = DATA.C.WIS)
SEL.S.WIS <- lm(relFIT ~ scale(RPM56) + scale(INT56) + GEN +
    WHLSTAGE + sex + Fcoeff + line,
  data = DATA.S.WIS)
RUN.S <- lm(relFIT ~ scale(RUN56)+ GEN +
    WHLSTAGE + sex + Fcoeff + line,
  data = DATA.S.WIS)

#extract betas OVERALL
BETA.RPM.C.WIS <- summary(SEL.C.WIS)$coefficients[2,1:2]
BETA.INT.C.WIS <- summary(SEL.C.WIS)$coefficients[3,1:2]
BETA.RPM.S.WIS <- summary(SEL.S.WIS)$coefficients[2,1:2]
BETA.INT.S.WIS <- summary(SEL.S.WIS)$coefficients[3,1:2]
############## estimate BETAS in the UCR generations #########################
DATA.C.UCR <- subset(dat, linetype == "0" & UNI == "UCR")
DATA.S.UCR <- subset(dat, linetype == "1" & UNI == "UCR")
#calculate relative fitness using the number of pups produced
DATA.C.UCR$relFIT <-DATA.C.UCR$pups / mean(DATA.C.UCR$pups)
DATA.S.UCR$relFIT <-DATA.S.UCR$pups / mean(DATA.S.UCR$pups)
#fit linear model
SEL.C.UCR<-lm(relFIT ~ scale(RPM56) + scale(INT56) + GEN +
    WHLSTAGE + sex + Fcoeff + line,
  data = DATA.C.UCR)
SEL.S.UCR<-lm(relFIT ~ scale(RPM56) + scale(INT56) + GEN +
    WHLSTAGE + sex + Fcoeff + line,
  data = DATA.S.UCR)
#extract betas OVERALL
BETA.RPM.C.UCR <- summary(SEL.C.UCR)$coefficients[2,1:2]
BETA.INT.C.UCR <- summary(SEL.C.UCR)$coefficients[3,1:2]
BETA.RPM.S.UCR <- summary(SEL.S.UCR)$coefficients[2,1:2]
BETA.INT.S.UCR <- summary(SEL.S.UCR)$coefficients[3,1:2]
################### estimate BETAS in each generation ##########################
DATA.C<-subset(dat, linetype == "0")
DATA.S<-subset(dat, linetype == "1")
#summary(DATA.C)
BETAS <- data.frame(BLOCK = NA, GEN = 0:77,
  B.RPM.C = 0, se.RPM.C = 0, B.INT.C = 0, se.INT.C = 0,
  B.RPM.S = 0, se.RPM.S = 0, B.INT.S = 0, se.INT.S = 0,
  B.RUN.C = 0, se.RUN.C = 0, B.RUN.S = 0, se.RUN.S = 0,
  X.RUN.S = 0, X.RPM.S = 0, X.INT.S = 0,
  sd.RUN.S = 0,sd.RPM.S = 0, sd.INT.S = 0)
for(g in GEN.lst) {
  #estimate gradients for the first block
  BLOCK.C <- subset(dat, GEN == g & linetype == "0")
  BLOCK.S <- subset(dat, GEN == g & linetype == "1")
  #calculate relative fitness using the number of pups produced
  BLOCK.C$relFIT <- BLOCK.C$pups / mean(BLOCK.C$pups)
  BLOCK.S$relFIT <- BLOCK.S$pups / mean(BLOCK.S$pups)
  #fit linear model
  SEL.C <- lm(relFIT ~ scale(RPM56) + scale(INT56) + WHLSTAGE + sex + line,
    data = BLOCK.C)
  SEL.S <- lm(relFIT ~ scale(RPM56) + scale(INT56) + WHLSTAGE + sex + line,
    data = BLOCK.S)
  RUN.C <- lm(relFIT ~ scale(RUN56) + WHLSTAGE + sex + line, data = BLOCK.C)
  RUN.S <- lm(relFIT ~ scale(RUN56) + WHLSTAGE + sex + line, data = BLOCK.S)
  #extract betas
  betasGENind <- which(BETAS$GEN == g) 
  BETAS[betasGENind, 3:4] <- summary(SEL.C)$coefficients[2, 1:2]
  BETAS[betasGENind, 5:6] <- summary(SEL.C)$coefficients[3, 1:2]
  BETAS[betasGENind, 7:8] <- summary(SEL.S)$coefficients[2, 1:2]
  BETAS[betasGENind, 9:10] <- summary(SEL.S)$coefficients[3, 1:2]
  BETAS[betasGENind, 11:12] <- summary(RUN.C)$coefficients[2, 1:2]
  BETAS[betasGENind, 13:14] <- summary(RUN.S)$coefficients[2, 1:2]
  #
  BETAS$X.RUN.S[betasGENind] <- mean(BLOCK.S$RUN56)
  BETAS$X.RPM.S[betasGENind] <- mean(BLOCK.S$RPM56)
  BETAS$X.INT.S[betasGENind] <- mean(BLOCK.S$INT56)
  BETAS$sd.RUN.S[betasGENind] <- sd(BLOCK.S$RUN56)
  BETAS$sd.RPM.S[betasGENind] <- sd(BLOCK.S$RPM56)
  BETAS$sd.INT.S[betasGENind] <- sd(BLOCK.S$INT56)
}

#selection gradients throughout the experiment
BETAS$n <- 1
BETAS[!BETAS$GEN %in% GEN.lst, 3:11] <- NA
BETAS$GEN.C <- BETAS$GEN - 0.05
BETAS$GEN.S <- BETAS$GEN + 0.05
#
#calculate MEAN + SE BETAS
BETAS$UNI <- "UCR"
BETAS$UNI[which(BETAS$GEN < 32)] <- "WIS"
BETAS[!BETAS$GEN %in% GEN.lst, ] <- NA
BETAS.WIS <- subset(BETAS, UNI == "WIS")
BETAS.UCR <- subset(BETAS, UNI == "UCR")
#
BETA.RUN.C.WIS <- cbind(mean(BETAS.WIS$B.RUN.C, na.rm = TRUE),
  sd(BETAS.WIS$B.RUN.C, na.rm = TRUE) / sqrt(sum(BETAS.WIS$n, na.rm = TRUE)))
BETA.RUN.S.WIS <- cbind(mean(BETAS.WIS$B.RUN.S, na.rm = TRUE),
  sd(BETAS.WIS$B.RUN.S, na.rm = TRUE) / sqrt(sum(BETAS.WIS$n, na.rm = TRUE)))
BETA.RPM.C.WIS<-cbind(mean(BETAS.WIS$B.RPM.C,na.rm = TRUE),
  sd(BETAS.WIS$B.RPM.C, na.rm = TRUE) / sqrt(sum(BETAS.WIS$n, na.rm = TRUE)))
BETA.INT.C.WIS<-cbind(mean(BETAS.WIS$B.INT.C, na.rm = TRUE),
  sd(BETAS.WIS$B.INT.C, na.rm = TRUE) / sqrt(sum(BETAS.WIS$n, na.rm = TRUE)))
BETA.RPM.S.WIS<-cbind(mean(BETAS.WIS$B.RPM.S, na.rm = TRUE),
  sd(BETAS.WIS$B.RPM.S, na.rm = TRUE) / sqrt(sum(BETAS.WIS$n, na.rm = TRUE)))
BETA.INT.S.WIS<-cbind(mean(BETAS.WIS$B.INT.S, na.rm = TRUE),
  sd(BETAS.WIS$B.INT.S, na.rm = TRUE) / sqrt(sum(BETAS.WIS$n, na.rm = TRUE)))
#
BETA.RUN.C.UCR <- cbind(mean(BETAS.UCR$B.RUN.C, na.rm = TRUE),
  sd(BETAS.UCR$B.RUN.C, na.rm = TRUE) / sqrt(sum(BETAS.UCR$n, na.rm = TRUE)))
BETA.RUN.S.UCR <- cbind(mean(BETAS.UCR$B.RUN.S, na.rm = TRUE),
  sd(BETAS.UCR$B.RUN.S, na.rm = TRUE) / sqrt(sum(BETAS.UCR$n, na.rm = TRUE)))
BETA.RPM.C.UCR <- cbind(mean(BETAS.UCR$B.RPM.C, na.rm = TRUE),
  sd(BETAS.UCR$B.RPM.C, na.rm = TRUE) / sqrt(sum(BETAS.UCR$n, na.rm = TRUE)))
BETA.INT.C.UCR <- cbind(mean(BETAS.UCR$B.INT.C, na.rm = TRUE),
  sd(BETAS.UCR$B.INT.C, na.rm = TRUE) / sqrt(sum(BETAS.UCR$n, na.rm = TRUE)))
BETA.RPM.S.UCR <- cbind(mean(BETAS.UCR$B.RPM.S, na.rm = TRUE),
  sd(BETAS.UCR$B.RPM.S, na.rm = TRUE) / sqrt(sum(BETAS.UCR$n, na.rm = TRUE)))
BETA.INT.S.UCR <- cbind(mean(BETAS.UCR$B.INT.S, na.rm = TRUE),
  sd(BETAS.UCR$B.INT.S, na.rm = TRUE) / sqrt(sum(BETAS.UCR$n, na.rm = TRUE)))

#numbers provided in the main text
#"with a noticeable reduction in the average standardized selection gradients (se) in the Riverside generations
BETA.RUN.S.UCR
BETA.RPM.S.UCR
BETA.INT.S.UCR
#compared to Wisconsin generations"
BETA.RUN.S.WIS
BETA.RPM.S.WIS
BETA.INT.S.WIS

############ estimate BETAS in each generation for each LINE ###################
BETAS.LINE <- data.frame(BLOCK = NA, GEN = 0:77)
nms <- paste(rep(paste(rep(c("B", "se"), 3), rep(c("RPM", "INT", "RUN"), each = 2),
                sep = "."), 8),
      rep(seq.int(8), each = 6), sep = ".")
BETAS.LINE[, nms] <- 0
                 
for(g in GEN.lst) {
  for(l in 1:8){
    #estimate gradients for the first block
    if(l %in% c(1, 2, 4, 5)){
      BLOCK <- subset(DATA.C, subset = GEN == g & line == l)
    } else{
        BLOCK <- subset(DATA.S, subset = GEN == g & line == l)
      }
    #calculate relative fitness using the number of pups produced
    BLOCK$relFIT <- BLOCK$pups / mean(BLOCK$pups)
    #fit linear model
    SEL <- lm(relFIT ~ scale(RPM56) + scale(INT56) + WHLSTAGE + sex + Fcoeff,
      data = BLOCK)
    SEL.RUN <-lm(relFIT ~ scale(RUN56l) + WHLSTAGE + sex + Fcoeff,
      data = BLOCK)
    #extract BETAS.line
    BETAS.LINE[which(BETAS.LINE$GEN == g), 
               paste(c("B", "se"), "RPM", l, sep = ".")] <-
                                               summary(SEL)$coefficients[2, 1:2]
    BETAS.LINE[which(BETAS.LINE$GEN == g),
               paste(c("B", "se"), "INT", l, sep = ".")] <-
                                               summary(SEL)$coefficients[3, 1:2]
    BETAS.LINE[which(BETAS.LINE$GEN == g),
               paste(c("B", "se"), "RUN", l, sep = ".")]   <-
                                           summary(SEL.RUN)$coefficients[2, 1:2]
  }  #<-- end for l
}  #<-- end for g



# selection gradients throughout the experiment
BETAS.LINE$n <- 1
BETAS.LINE[!BETAS.LINE$GEN %in% GEN.lst,] <- NA
BETAS.LINE$GEN.1 <- BETAS.LINE$GEN - 0.1
BETAS.LINE$GEN.2 <- BETAS.LINE$GEN - 0.075
BETAS.LINE$GEN.4 <- BETAS.LINE$GEN - 0.05
BETAS.LINE$GEN.5 <- BETAS.LINE$GEN - 0.025
BETAS.LINE$GEN.3 <- BETAS.LINE$GEN + 0.1
BETAS.LINE$GEN.6 <- BETAS.LINE$GEN + 0.075
BETAS.LINE$GEN.7 <- BETAS.LINE$GEN + 0.05
BETAS.LINE$GEN.8 <- BETAS.LINE$GEN + 0.025

# calculate MEAN + SE BETAS.LINE
for(l in 1:8){
  assign(paste0("tmpBETA.RPM.", l),
    cbind(mean(BETAS.LINE[, paste0("B.RPM.", l)], na.rm = TRUE),
          sd(BETAS.LINE[, paste0("B.RPM.", l)], na.rm = TRUE) /
                                         sqrt(sum(BETAS.LINE$n, na.rm = TRUE))))
  assign(paste0("tmpBETA.INT.", l),
    cbind(mean(BETAS.LINE[, paste0("B.INT.", l)], na.rm = TRUE),
          sd(BETAS.LINE[, paste0("B.INT.", l)], na.rm = TRUE) /
                                         sqrt(sum(BETAS.LINE$n, na.rm = TRUE))))
}







# use thin-plate splines to display the fitness surfaces
LINE.3 <- subset(dat, line == 3 & GEN < 78)
LINE.6 <- subset(dat, line == 6 & GEN < 78)
LINE.7 <- subset(dat, line == 7 & GEN < 78)
LINE.8 <- subset(dat, line == 8 & GEN < 78)
LINE.3[, c("RPM56z", "INT56z", "relFIT")] <- NA
LINE.6[, c("RPM56z", "INT56z", "relFIT")] <- NA
LINE.7[, c("RPM56z", "INT56z", "relFIT")] <- NA
LINE.8[, c("RPM56z", "INT56z", "relFIT")] <- NA

for(l in c(3, 6, 7, 8)){
  tmpL <- get(paste0("LINE.", l))
  for(g in GEN.lst) {
    genINDl <- which(tmpL$GEN == g)
    tmpL$RPM56z[genINDl] <- scale(tmpL$RPM56[genINDl])
    tmpL$INT56z[genINDl] <- scale(tmpL$INT56[genINDl])
    tmpL$relFIT[genINDl] <- tmpL$pups[genINDl] / mean(tmpL$pups[genINDl],
                                                                   na.rm = TRUE)
  }  #<-- end for g
  assign(paste0("LINE.", l), tmpL)
}  #<-- end for l




# fields::Tps()
#XXX WARNING: since each line has a lot of data (>6,000 rows) the call to `Tps()`
## takes a while. These calls are commented out and objects created in a previous
### session are loaded from the hard drive
#TPS.3 <- Tps(data.matrix(data.frame(LINE.3$INT56z,LINE.3$RPM56z)), LINE.3$relFIT)
#  save(TPS.3,file = "TPS.3.scale.RData")
#TPS.6 <- Tps(data.matrix(data.frame(LINE.6$INT56z,LINE.6$RPM56z)), LINE.6$relFIT)
#  save(TPS.6,file = "TPS.6.scale.RData")
#TPS.7 <- Tps(data.matrix(data.frame(LINE.7$INT56z,LINE.7$RPM56z)), LINE.7$relFIT)
#  save(TPS.7,file = "TPS.7.scale.RData")
#TPS.8 <- Tps(data.matrix(data.frame(LINE.8$INT56z,LINE.8$RPM56z)), LINE.8$relFIT)
#  save(TPS.8,file = "TPS.8.scale.RData")
load(file = "TPS.3.scale.RData")
load(file = "TPS.6.scale.RData")
load(file = "TPS.7.scale.RData")
load(file = "TPS.8.scale.RData")

sp.3 <- predictSurface(TPS.3)
sp.6 <- predictSurface(TPS.6)
sp.7 <- predictSurface(TPS.7)
sp.8 <- predictSurface(TPS.8)
ZLIM <- range(c(sp.3$z, sp.6$z, sp.7$z, sp.8$z), na.rm = TRUE)
LIM <- range(c(scale(LINE.3[, c("INT56", "RPM56")]),
    scale(LINE.6[, c("INT56", "RPM56")]),
    scale(LINE.7[, c("INT56", "RPM56")]),
    scale(LINE.8[, c("INT56", "RPM56")])),
  na.rm = TRUE)*1.45


# Figure 1
dev.new(width=10,height=4, units = "cm")
par(mfrow=c(6,5),las = 1, oma = c(4,3,0.5,0), mar = c(0,2,2,2))
layout(matrix(c(1,1,1,2,3,
                1,1,1,2,3,
                4,4,4,2,3,
                4,4,4,5,6,
                7,7,7,5,6,
                7,7,7,5,6), 6, 5, byrow = TRUE))
#layout.show(7)
plot(B.RUN.S ~ GEN, data = BETAS, pch = 16, col = "red", cex = 1, type="o",
  ylim = c(-0.675, 1.5), xlim = c(1, 77), xlab = "", xaxt = "n", ylab = "")
  points(B.RUN.C ~ GEN.C, data = BETAS,
    pch = 16, col = "blue", cex = 1, type = "o")
  rect(xleft = 0, ybottom = BETA.RUN.S.WIS[1]-1.96*BETA.RUN.S.WIS[2],
    xright = 31, ytop = BETA.RUN.S.WIS[1]+1.96*BETA.RUN.S.WIS[2],
    col = rgb(1,0,0,0.25), border = FALSE)
  rect(xleft = 35, ybottom = BETA.RUN.S.UCR[1]-1.96*BETA.RUN.S.UCR[2],
    xright = 78, ytop = BETA.RUN.S.UCR[1]+1.96*BETA.RUN.S.UCR[2],
    col = rgb(1,0,0,0.25), border = FALSE)
  rect(xleft = 0, ybottom = BETA.RUN.C.WIS[1]-1.96*BETA.RUN.C.WIS[2],
    xright = 31, ytop = BETA.RUN.C.WIS[1]+1.96*BETA.RUN.C.WIS[2],
    col = rgb(0,0,1,0.25), border = FALSE)
  rect(xleft = 35, ybottom = BETA.RUN.C.UCR[1]-1.96*BETA.RUN.C.UCR[2],
    xright = 78, ytop = BETA.RUN.C.UCR[1]+1.96*BETA.RUN.C.UCR[2],
    col = rgb(0,0,1,0.25), border = FALSE)
  arrows(x0 = BETAS$GEN.C, y0 = BETAS$B.RUN.C-BETAS$se.RUN.C,
    x1 = BETAS$GEN.C, y1 = BETAS$B.RUN.C+BETAS$se.RUN.C,
    col = "blue", code = 3, angle = 90, length = 0)
  arrows(x0 = BETAS$GEN.S, y0 = BETAS$B.RUN.S-BETAS$se.RUN.S,
    x1 = BETAS$GEN.S, y1 = BETAS$B.RUN.S+BETAS$se.RUN.S,
    col = "red", code = 3, angle = 90, length = 0)
  points(B.RUN.1~GEN, BETAS.LINE, pch = "1", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.RUN.2~GEN, BETAS.LINE, pch = "2", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.RUN.4~GEN, BETAS.LINE, pch = "4", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.RUN.5~GEN, BETAS.LINE, pch = "5", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.RUN.3~GEN, BETAS.LINE, pch = "3", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.RUN.6~GEN, BETAS.LINE, pch = "6", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.RUN.7~GEN, BETAS.LINE, pch = "7", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.RUN.8~GEN, BETAS.LINE, pch = "8", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
  abline(h = 0, lty = 2)
mtext("A", side = 3, adj = 0.025, line = -1.5)
legend(x = 47, y = 14, legend = c("control", "selected"),
  pch = c(16, 17), col = c("blue", "red"), bty = "n")
mtext("Wisconsin generations", line = 0.4, adj = 0.05, side = 3, cex = 1)
mtext("Riverside generations", line = 0.4, adj = 0.80, side = 3, cex = 1)
mtext(side = 2, expression(italic(beta)[distance]), las = 3, line = 3)
  
#
image(sp.3, col = heat.colors(n = 30, rev = TRUE),
  zlim = ZLIM, ylim = LIM, xlim = LIM)
  abline(v = 0, lty = 3)
  abline(h = 0, lty = 3)
contour(sp.3, add = TRUE)
mtext("D", side = 3, adj = 0.025, line = -1.5)
mtext("line S3", side = 1, cex = 0.75, line = -1, adj = 0.95)
image(sp.6, col = heat.colors(n = 30, rev = TRUE),
  zlim = ZLIM, ylim = LIM, xlim = LIM)
  abline(v = 0, lty = 3)
  abline(h = 0, lty = 3)
contour(sp.6, add = TRUE)
mtext("E", side = 3, adj = 0.025, line = -1.5)
mtext("line S6", side = 1, cex = 0.75, line = -1, adj = 0.95)

#
plot(B.RPM.S ~ GEN.C, BETAS, pch = 16, col = "red", cex = 1, type = "o",
  ylim = c(-0.675,1.5), xlim = c(1,77), xlab = "", xaxt = "n", ylab = "")
  points(B.RPM.C ~ GEN.C, BETAS, pch = 16, col = "blue", cex = 1, type = "o")
  rect(xleft = 0, ybottom = BETA.RPM.S.WIS[1]-1.96*BETA.RPM.S.WIS[2],
    xright = 31, ytop = BETA.RPM.S.WIS[1]+1.96*BETA.RPM.S.WIS[2],
    col = rgb(1,0,0,0.25), border = FALSE)
  rect(xleft = 35, ybottom = BETA.RPM.S.UCR[1]-1.96*BETA.RPM.S.UCR[2],
    xright = 78, ytop = BETA.RPM.S.UCR[1]+1.96*BETA.RPM.S.UCR[2],
    col = rgb(1,0,0,0.25), border = FALSE)
  rect(xleft = 0, ybottom = BETA.RPM.C.WIS[1]-1.96*BETA.RPM.C.WIS[2],
    xright = 31, ytop = BETA.RPM.C.WIS[1]+1.96*BETA.RPM.C.WIS[2],
    col = rgb(0,0,1,0.25), border = FALSE)
  rect(xleft = 35, ybottom = BETA.RPM.C.UCR[1]-1.96*BETA.RPM.C.UCR[2],
    xright = 78, ytop = BETA.RPM.C.UCR[1]+1.96*BETA.RPM.C.UCR[2],
    col = rgb(0,0,1,0.25), border = FALSE)
  arrows(x0 = BETAS$GEN.C, y0 = BETAS$B.RPM.C-BETAS$se.RPM.C,
    x1 = BETAS$GEN.C, y1 = BETAS$B.RPM.C+BETAS$se.RPM.C,
    col ="blue", code = 3, angle = 90, length = 0)
  arrows(x0 = BETAS$GEN.S, y0 = BETAS$B.RPM.S-BETAS$se.RPM.S,
    x1 = BETAS$GEN.S, y1 = BETAS$B.RPM.S+BETAS$se.RPM.S,
    col = "red", code = 3, angle = 90, length = 0)
  points(B.RPM.1~GEN, BETAS.LINE, pch = "1", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.RPM.2~GEN, BETAS.LINE, pch = "2", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.RPM.4~GEN, BETAS.LINE, pch = "4", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.RPM.5~GEN, BETAS.LINE, pch = "5", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.RPM.3~GEN, BETAS.LINE, pch = "3", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.RPM.6~GEN, BETAS.LINE, pch = "6", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.RPM.7~GEN, BETAS.LINE, pch = "7", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.RPM.8~GEN, BETAS.LINE, pch = "8", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
  abline(h = 0, lty = 2)
mtext("B", side = 3, adj = 0.025, line = -1.5)
mtext(side = 2, expression(italic(beta)[speed]), las = 3, line = 3)
  
#
image(sp.7, col = heat.colors(n=30,rev=T),
  zlim = ZLIM, ylim = LIM, xlim = LIM)
  abline(v = 0, lty = 3)
  abline(h = 0, lty = 3)
contour(sp.7, add = TRUE)
mtext("F", side = 3, adj = 0.025, line = -1.5)
mtext("line S7", side = 1, cex = 0.75, line = -1, adj = 0.95)
mtext("Running speed (sd units)", side = 2, line = 2, cex = 1, adj = -2, las = 3)
image(sp.8, col = heat.colors(n = 30, rev = TRUE),
  zlim = ZLIM, ylim = LIM, xlim = LIM)
  abline(v = 0, lty = 3)
  abline(h = 0, lty = 3)
contour(sp.8, add = TRUE)
mtext("G", side = 3, adj = 0.025, line = -1.5)
mtext("line S8", side = 1, cex = 0.75, line = -1, adj = 0.95)

#
plot(B.INT.C~GEN.C, BETAS, pch = 16, col = "white", cex = 1, type = "o",
  ylim = c(-1, 1.5), xlim = c(1, 77), xlab = "", xaxt = "n", ylab = "")
  rect(xleft = 0, ybottom = BETA.INT.S.WIS[1]-1.96*BETA.INT.S.WIS[2],
    xright = 31, ytop = BETA.INT.S.WIS[1]+1.96*BETA.INT.S.WIS[2], 
    col = rgb(1,0,0,0.25), border = FALSE)
  rect(xleft = 35, ybottom = BETA.INT.S.UCR[1]-1.96*BETA.INT.S.UCR[2],
    xright = 78, ytop = BETA.INT.S.UCR[1]+1.96*BETA.INT.S.UCR[2], 
    col = rgb(1,0,0,0.25), border = FALSE)
  rect(xleft = 0, ybottom = BETA.INT.C.WIS[1]-1.96*BETA.INT.C.WIS[2],
    xright = 31, ytop = BETA.INT.C.WIS[1]+1.96*BETA.INT.C.WIS[2], 
    col = rgb(0,0,1,0.25), border = FALSE)
  rect(xleft = 35, ybottom = BETA.INT.C.UCR[1]-1.96*BETA.INT.C.UCR[2],
    xright = 78, ytop = BETA.INT.C.UCR[1]+1.96*BETA.INT.C.UCR[2], 
    col = rgb(0,0,1,0.25), border = FALSE)
  points(B.INT.C~GEN.C, BETAS, pch = 16, col = "blue", cex = 1, type = "o")
  points(B.INT.S~GEN.S, BETAS, pch = 17, col = "red", cex = 1, type = "o")
  arrows(x0 = BETAS$GEN.C, y0 = BETAS$B.INT.C-BETAS$se.INT.C,
    x1=BETAS$GEN.C, y1 = BETAS$B.INT.C+BETAS$se.INT.C,
    col = "blue", code = 3, angle = 90, length = 0)
  arrows(x0 = BETAS$GEN.S, y0 = BETAS$B.INT.S-BETAS$se.INT.S,
    x1 = BETAS$GEN.S, y1 = BETAS$B.INT.S+BETAS$se.INT.S,
    col = "red", code = 3, angle = 90, length = 0)
  points(B.INT.1~GEN, BETAS.LINE, pch = "1", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.INT.2~GEN, BETAS.LINE, pch = "2", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.INT.4~GEN, BETAS.LINE, pch = "4", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.INT.5~GEN, BETAS.LINE, pch = "5", col = rgb(0,0,1,0.25), cex = 1, type = "l")
  points(B.INT.3~GEN, BETAS.LINE, pch = "3", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.INT.6~GEN, BETAS.LINE, pch = "6", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.INT.7~GEN, BETAS.LINE, pch = "7", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  points(B.INT.8~GEN, BETAS.LINE, pch = "8", col = rgb(1,0,0,0.25), cex = 1, type = "l")
  abline(h = 0, lty = 2)
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
mtext("C", side = 3, adj = 0.025, line = -1.5)
mtext(side = 2, expression(italic(beta)[duration]), las = 3, line = 3)

#
mtext("Generation", side = 1, line = 2.5, cex = 1)
mtext("Running duration (sd units)", side = 1, line = 2.5, cex = 1, adj = 2)















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
dat 	   <- read.table("QG-mouse-trade-off_data.txt", header = TRUE)
library(visreg)
library(splines)
##################################
dat$n <- 1
MEANS        <- aggregate(RPM56 ~ GEN + line, data = dat, FUN = mean)
MEANS$SD.RPM <- aggregate(RPM56 ~ GEN + line, data = dat, FUN = sd)[3]
MEANS$N      <- aggregate(n ~ GEN + line, data = dat, FUN = sum)$n
MEANS$SE.RPM <- MEANS$SD.RPM / sqrt(MEANS$N)
MEANS$INT    <- aggregate(INT56 ~ GEN + line, data = dat, FUN = mean)$INT
MEANS$SD.INT <- aggregate(INT56 ~ GEN + line, data = dat, FUN = sd)[3]
MEANS$SE.INT <- MEANS$SD.INT / sqrt(MEANS$N)
MEANS$RUN    <- aggregate(RUN56 ~ GEN + line, data = dat, FUN = mean)$RUN56
MEANS$SD.RUN <- aggregate(RUN56 ~ GEN + line, data = dat, FUN = sd)[3]
MEANS$SE.RUN <- MEANS$SD.RUN / sqrt(MEANS$N)
  names(MEANS)[c(2:3)] <- c("LINE", "RPM")

MEANS$linetype <- "S"  #<-- Selected Lines
MEANS$linetype[which(MEANS$LINE == "1")] <- "C"  #<-- Control line
MEANS$linetype[which(MEANS$LINE == "2")] <- "C"  #<-- Control line
MEANS$linetype[which(MEANS$LINE == "4")] <- "C"  #<-- Control line
MEANS$linetype[which(MEANS$LINE == "5")] <- "C"  #<-- Control line


# add RATIO
MEANS.C       <- aggregate(RPM56 ~ GEN,
			data = subset(dat, linetype == 0), FUN = mean)
MEANS.C$INT56 <- aggregate(INT56 ~ GEN,
			data = subset(dat, linetype == 0), FUN = mean)$INT56
MEANS.C$RUN56 <- aggregate(RUN56 ~ GEN,
			data = subset(dat, linetype == 0), FUN = mean)$RUN56
MEANS.S       <- aggregate(RPM56 ~ GEN,
			data = subset(dat, linetype == 1), FUN = mean)
MEANS.S$INT56 <- aggregate(INT56 ~ GEN,
			data = subset(dat, linetype == 1), FUN = mean)$INT56
MEANS.S$RUN56 <- aggregate(RUN56 ~ GEN,
			data = subset(dat, linetype == 1), FUN = mean)$RUN56
RATIOS <- data.frame(GEN = MEANS$GEN[1:71],
  RPM = MEANS.S$RPM56 /MEANS.C$RPM56,
  INT = MEANS.S$INT56 / MEANS.C$INT56,
  RUN = MEANS.S$RUN56 / MEANS.C$RUN56)




#.################ section #2D - selection GAIN vs differential ###############
#.################ section #2D - selection GAIN vs differential ###############
#.################ section #2D - selection GAIN vs differential ###############
#.################ section #2D - selection GAIN vs differential ###############
#.################ section #2D - selection GAIN vs differential ###############
#.################ section #2D - selection GAIN vs differential ###############
#.################ section #2D - selection GAIN vs differential ###############
dat <- dat[order(dat$GEN), ]
DAT <- dat[which(dat$GEN < 78), ]
DAT <- DAT[which(dat$GEN > -1), ]
DAT$BRED <- 0  #<-- column to indicate individuals chosen as breeders
DAT$BRED[which(DAT$pups > 0)] <- 1  #<-- assign 1 for breeders
DATb <- DAT[which(DAT$BRED > 0), ]  #<-- object which contains only BREEDERS
SEL   <- aggregate(RUN56 ~ GEN + line, data = DAT, FUN = mean) #RUN56 average per generation
colnames(SEL) <- c("GEN", "LINE", "Xp")  #(Xp) = mean population
SEL$SD <- aggregate(RUN56 ~ GEN + line, data = DAT, FUN = sd)$RUN56
SEL$Xb <- aggregate(RUN56 ~ GEN + line, data = DATb, FUN = mean)$RUN56 #(Xb) = mean breeders
SEL$S.abs <- SEL$Xb - SEL$Xp # RUN56 absolute selection differential

CONTROL <- aggregate(RUN56 ~ GEN + linetype, data = DAT, FUN = mean)
CONTROL <- subset(CONTROL, linetype == 0)

SEL <- merge(SEL, CONTROL, by = "GEN")
SEL$linetype <- NULL
  names(SEL)[match("RUN56", names(SEL))] <- "Xp.C"
SEL$S.std <- SEL$S.abs / SEL$SD
SEL.3 <- subset(SEL, LINE == 3)
SEL.6 <- subset(SEL, LINE == 6)
SEL.7 <- subset(SEL, LINE == 7)
SEL.8 <- subset(SEL, LINE == 8)

# split by CONTROL vs SELECTED lines
## Selective gain in the S line
SEL.3$GAIN <- with(SEL.3, Xp - Xp.C)
SEL.6$GAIN <- with(SEL.6, Xp - Xp.C)
SEL.7$GAIN <- with(SEL.7, Xp - Xp.C)
SEL.8$GAIN <- with(SEL.8, Xp - Xp.C)
## cummulative selection differential
SEL.3$S.cum <- with(SEL.3, cumsum(S.std) - S.std)
SEL.6$S.cum <- with(SEL.6, cumsum(S.std) - S.std)
SEL.7$S.cum <- with(SEL.7, cumsum(S.std) - S.std)
SEL.8$S.cum <- with(SEL.8, cumsum(S.std) - S.std)

SEL.3$GAIN.std <- SEL.3$GAIN / 2050
SEL.6$GAIN.std <- SEL.6$GAIN / 2050
SEL.7$GAIN.std <- SEL.7$GAIN / 2050
SEL.8$GAIN.std <- SEL.8$GAIN / 2050


DF <- 6
spline_model.3 <- lm(GAIN.std ~ ns(S.cum, df = DF), data = SEL.3)
spline_model.6 <- lm(GAIN.std ~ ns(S.cum, df = DF), data = SEL.6)
spline_model.7 <- lm(GAIN.std ~ ns(S.cum, df = DF), data = SEL.7)
spline_model.8 <- lm(GAIN.std ~ ns(S.cum, df = DF), data = SEL.8)

# estimate realized heritability in UCR era
SEL.3$ERA <- SEL.6$ERA <- SEL.7$ERA <- SEL.8$ERA <- 16
SEL.3$ERA[which(SEL.3$GEN > 31)] <- 17
SEL.6$ERA[which(SEL.6$GEN > 31)] <- 17
SEL.7$ERA[which(SEL.7$GEN > 31)] <- 17
SEL.8$ERA[which(SEL.8$GEN > 31)] <- 17

(lm.3 <- summary(lm(GAIN.std ~ S.cum, data = subset(SEL.3, ERA == 17))))
(lm.6 <- summary(lm(GAIN.std ~ S.cum, data = subset(SEL.6, ERA == 17))))
(lm.7 <- summary(lm(GAIN.std ~ S.cum, data = subset(SEL.7, ERA == 17))))
(lm.8 <- summary(lm(GAIN.std ~ S.cum, data = subset(SEL.8, ERA == 17))))

# adjusting for the within-family selection index, using the following formula:
# h2 = h2R(1-t)/1-r)
# where h2 is the overall heritability expected from mass selection
# or a randomly mating population, h2, is the realized heritability obtained for this particular
# index, r is the coefficient of relationship of full sibs (0.5)
# and t is the intraclass correlation of full sibs for nesting scores
IC <- 0.246 #as in Careau et al. (2013; Evolution)
(h2w.3 <- lm.3$coefficients[2, 1] * (1-0.5) / (1-IC))
(h2w.6 <- lm.6$coefficients[2, 1] * (1-0.5) / (1-IC))
(h2w.7 <- lm.7$coefficients[2, 1] * (1-0.5) / (1-IC))
(h2w.8 <- lm.8$coefficients[2, 1] * (1-0.5) / (1-IC))
mean(h2w.3, h2w.6, h2w.7, h2w.8)                 # = average h2w
sd(c(h2w.3, h2w.6, h2w.7, h2w.8)) / sqrt(4)        # = se for h2w with N = 4

LIMITS <- data.frame(line = c(3, 6, 7, 8),
  limit = c(17.78, 16.39, 19.60, 25.25),
  CI = c(3.55, 1.92, 2.54, 3.57))
LIMITS$lcl <- LIMITS$limit - LIMITS$CI
LIMITS$ucl <- LIMITS$limit + LIMITS$CI




YLIM <- c(-0.7, 5.15)
# make Figure 2
dev.new(width = 8, height = 4, units = "cm")
par(mfrow = c(6, 3), las = 1, oma = c(3, 5, 1, 0), mar = c(0, 0, 1, 1))
layout(matrix(c(1,2,3,
                1,2,3,
                4,2,3,
                4,5,6,
                7,5,6,
                7,5,6), 6, 3, byrow = TRUE))
# Distance
plot(RUN ~ GEN, data = subset(MEANS, LINE == "1"),
  pch = 16, lty = 1, col = "blue", type = "l",
  xlim = c(0, 78), xaxt = "n", ylim = c(6, 15000))
  # draw lines for control lines
  for(l in c("2", "4", "5")){
    lines(RUN ~ GEN, data = subset(MEANS, LINE == l), lty = 1, col = "blue")
  }
  # draw lines for selected lines  
  for(l in c("3", "6", "7", "8")){
    lines(RUN ~ GEN, data = subset(MEANS, LINE == l), lty = 1, col = "red")
  }
legend(-4, 3000, "Control lines",
  cex = 1, lty = 1, lwd = 1, col = "blue", bty = "n")
legend(35, 3000, "Selected lines",
  cex = 1, lty = 1, lwd = 1, col = "red", bty = "n")

  rect(xleft = 31.1, ybottom = min(MEANS$RUN),
    xright = 35.9, ytop = max(MEANS$RUN), col = "white", border = "white")
  rect(xleft = 51.1, ybottom = min(MEANS$RUN),
    xright = 52.9, ytop = max(MEANS$RUN), col = "white", border = "white")
  rect(xleft = 62.1, ybottom = min(MEANS$RUN),
    xright = 64.9, ytop = max(MEANS$RUN), col = "white", border = "white")
  rect(xleft = 66.1, ybottom = min(MEANS$RUN),
    xright = 67.9, ytop = max(MEANS$RUN), col = "white", border = "white")
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
  mtext("revs/day", side = 2, las = 3, line = 3.5)
  mtext("A", side = 3, adj = 0.015, line = -1.5)
  mtext("Wisconsin", line = 0.4, adj = 0.1, side = 3, cex = 1)
  mtext("Riverside", line = 0.4, adj = 0.85, side = 3, cex = 1)

# Selection gain (2 lines, see below for other 2)
par(mar = c(1, 3, 1, 1))
visreg(spline_model.3, ylab = "", xlab = "",
  points = list(cex = 1, pch = 16, col = rgb(0,0,0,0)),
  line = list(col = rgb(1,0,0,0.5)),
  xlim = c(0, 50), ylim = YLIM)
  arrows(x0 = LIMITS$ucl[1], y0 = -0.5, x1 = LIMITS$lcl[1], y1 = -0.5,
    code = 3, length = 0, col = 16)
  points(LIMITS$limit[1], y = -0.5, pch = 15, col = 16)
  points(GAIN.std ~ S.cum, data = SEL.3, pch = SEL.3$ERA, col = SEL.3$ERA)
  text(35, 1, "slope = 0.072\u00B10.027**")
  mtext("line S3", side = 1, cex = 0.75, line = -1, adj = 0.95)
  mtext("D", side = 3, adj = 0.025, line = -1.5)
  legend(-1,4.5, c("UCR","WIS"), col = c(17, 16), pch = c(17, 16), bty = "n")
     clip(x1 = SEL.3$S.cum[which(SEL.3$GEN == 36)],
       x2 = SEL.3$S.cum[which(SEL.3$GEN == 77)], y1 = 0, y2 = 10)
  abline(lm(GAIN.std ~ S.cum, data = subset(SEL.3, ERA == 17)), lwd = 2)
visreg(spline_model.6, ylab = "", xlab = "",
  points = list(cex = 1, pch = 16, col = rgb(0,0,0,0)),
  line = list(col = rgb(1,0,0,0.5)),
  xlim = c(0, 50), ylim = YLIM)
  arrows(x0 = LIMITS$ucl[2], y0 = -0.5, x1 = LIMITS$lcl[2], y1 = -0.5,
    code = 3, length = 0, col = 16)
  points(LIMITS$limit[2], y = -0.5, pch = 15,col = 16)
  points(GAIN.std ~ S.cum, data = SEL.6, pch = SEL.6$ERA, col = SEL.6$ERA)
  text(35, 1, "slope = 0.037\u00B10.017*")
  mtext("line S6", side = 1, cex = 0.75, line = -1, adj = 0.95)
  mtext("E", side = 3, adj = 0.025, line = -1.5)
     clip(x1 = SEL.6$S.cum[which(SEL.6$GEN == 36)],
       x2 = SEL.6$S.cum[which(SEL.6$GEN == 77)], y1 = 0, y2 = 10)
  abline(lm(GAIN.std ~ S.cum, data = subset(SEL.6, ERA == 17)), lwd = 2)
#


# Speed
par(mar = c(0,0,2,1))
plot(RPM ~ GEN, data = subset(MEANS, LINE == "1"),
  pch = 16, lty = 1, col = "blue", type = "l",
  xlim = c(0, 78), xaxt = "n", ylim = c(0, 30))
  # draw lines for control lines  
  for(l in c("2", "4", "5")){
    lines(RPM ~ GEN, data = subset(MEANS, LINE == l), lty = 1, col = "blue")
  }
  # draw lines for selected lines  
  for(l in c("3", "6", "7", "8")){
    lines(RPM ~ GEN, data = subset(MEANS, LINE == l), lty = 1, col = "red")
  }

  rect(xleft = 31.1, ybottom = min(MEANS$RPM),
    xright = 35.9, ytop = max(MEANS$RPM), col = "white", border = "white")
  rect(xleft = 51.1, ybottom = min(MEANS$RPM),
    xright = 52.9, ytop = max(MEANS$RPM), col = "white", border = "white")
  rect(xleft = 62.1, ybottom = min(MEANS$RPM),
    xright = 64.9, ytop = max(MEANS$RPM), col = "white", border = "white")
  rect(xleft = 66.1, ybottom = min(MEANS$RPM),
    xright = 67.9, ytop = max(MEANS$RPM), col = "white", border = "white")
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
  mtext("revs/min", side = 2, las = 3, line = 3.5)
  mtext("B", side = 3, adj = 0.015, line = -1.5)

# Selection gain (2 lines, see above for other 2)
par(mar = c(1,3,1,1))
visreg(spline_model.7, ylab = "", xlab = "",
  points = list(cex = 1, pch = 16, col = rgb(0,0,0,0)),
  line = list(col = rgb(1,0,0,0.5)),
  xlim = c(0,50), ylim = YLIM)
  arrows(x0 = LIMITS$ucl[3], y0 = -0.5, x1 = LIMITS$lcl[3], y1 = -0.5,
    code = 3, length = 0, col = 16)
  points(LIMITS$limit[3], y = -0.5, pch = 15, col = 16)
  points(GAIN.std ~ S.cum, data = SEL.7, pch = SEL.7$ERA, col = SEL.7$ERA)
  text(36,1,"slope = 0.041\u00B10.019*")
  mtext("line S7", side = 1, cex = 0.75, line = -1, adj = 0.95)
  mtext("F", side = 3, adj = 0.025, line = -1.5)
     clip(x1 = SEL.7$S.cum[which(SEL.7$GEN==36)],
       x2 = SEL.7$S.cum[which(SEL.7$GEN==77)], y1 = 0, y2 = 10)
  abline(lm(GAIN.std ~ S.cum, data = subset(SEL.7, ERA == 17)), lwd = 2)

mtext("Selective gain (revs/day; sd units)",
    side = 2, line = 2, adj = -0.5, las = 3)

visreg(spline_model.8, ylab = "", xlab = "",
  points = list(cex = 1, pch = 16, col = rgb(0,0,0,0)),
  line = list(col = rgb(1,0,0,0.5)),
  xlim = c(0, 50), ylim = YLIM)
  arrows(x0 = LIMITS$ucl[4], y0 = -0.5, x1 = LIMITS$lcl[4], y1 = -0.5,
    code = 3,length = 0, col = 16)
  points(LIMITS$limit[4], y = -0.5, pch = 15, col = 16)
  points(GAIN.std ~ S.cum, data = SEL.8, pch = SEL.8$ERA, col = SEL.8$ERA)
  text(36, 1, "slope = 0.066\u00B10.018**")
  mtext("line S8", side = 1, cex = 0.75, line = -1, adj = 0.95)
  mtext("G", side = 3, adj = 0.025, line = -1.5)
     clip(x1 = SEL.8$S.cum[which(SEL.8$GEN == 36)],
       x2 = SEL.8$S.cum[which(SEL.8$GEN == 77)], y1 = 0, y2 = 10)
  abline(lm(GAIN.std ~ S.cum, data = subset(SEL.8, ERA == 17)), lwd = 2)

mtext("Cummulative differential (sd units)", side = 1, line = 2.5, adj = 2.5)
#
# Duration
par(mar = c(0,0,2,1))
plot(INT ~ GEN, data = subset(MEANS, LINE == "1"),
  pch = 16, lty = 1, col = "blue", type = "l",
  xlim = c(0, 78), xaxt = "n", ylim = c(0, 650))
  # draw lines for control lines  
  for(l in c("2", "4", "5")){
    lines(INT ~ GEN, data = subset(MEANS, LINE == l), lty = 1, col = "blue")
  }
  # draw lines for selected lines  
  for(l in c("3", "6", "7", "8")){
    lines(INT ~ GEN, data = subset(MEANS, LINE == l), lty = 1, col = "red")
  }

  rect(xleft = 31.1, ybottom = min(MEANS$INT),
    xright = 35.9, ytop = max(MEANS$INT), col = "white", border = "white")
  rect(xleft = 51.1, ybottom = min(MEANS$INT),
    xright = 52.9, ytop = max(MEANS$INT), col = "white", border = "white")
  rect(xleft = 62.1, ybottom = min(MEANS$INT),
    xright = 64.9, ytop = max(MEANS$INT), col = "white", border = "white")
  rect(xleft = 66.1, ybottom = min(MEANS$INT),
    xright = 67.9, ytop = max(MEANS$INT), col = "white", border = "white")
  mtext("C", side = 3, adj = 0.015, line = -1.5)
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
  mtext("min/day", side = 2, las = 3, line = 3.5)

mtext("Generation", side = 1,  line = 1.5)
#

# END of figure






#calculate average S/C ratios after generation 25
mean(RATIOS$RUN[which(RATIOS$GEN > 24)])
mean(RATIOS$RPM[which(RATIOS$GEN > 24)])
mean(RATIOS$INT[which(RATIOS$GEN > 24)])

#calculate average S/C ratios for UCR gens
mean(RATIOS$RUN[which(RATIOS$GEN > 30)])
mean(RATIOS$RPM[which(RATIOS$GEN > 30)])
mean(RATIOS$INT[which(RATIOS$GEN > 30)])


















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
PED_pruned <- read.table("QG-mouse-trade-off_pedigree.txt", header = TRUE)
dat 	   <- read.table("QG-mouse-trade-off_data.txt", header = TRUE)

# prepare objects separately for Control and Selected mice 
## (data are going into separate models)
for(ltype in c("0", "1")){
  BLOCK <- subset(dat, GEN >= -1 & linetype == ltype)
  BLOCK$animal  <- factor(BLOCK$animal)
  BLOCK$line    <- factor(BLOCK$line)
  BLOCK$damid   <- factor(BLOCK$damid)
  BLOCK$GENfac  <- factor(BLOCK$GEN)
  PED.tmp.2  <- subset(PED_pruned, GEN > -2)
  PED.tmp.1 <- prepPed(PED.tmp.2)
  PED      <- prunePed(PED.tmp.1, unique(BLOCK$animal))
  AINVout <- inverseA(PED[, c("animal", "dam", "sire")])
  fout <- AINVout$inbreeding
  AINV <- AINVout$Ainv
  stopifnot(all(abs(fout[match(BLOCK$animal, PED$animal)] - BLOCK$Fcoeff) < 1e-6))
  # CONTROL MICE
  if(ltype == "0"){
    C0Ainv <- AINV
    C0datf <- BLOCK
    C0ped <- PED
  }
  # SELECTED MICE
  if(ltype == "1"){
    HR1Ainv <- AINV
    HR1datf <- BLOCK
    HR1ped <- PED
  }
}  #<-- end for ltype  

nrow(HR1datf)
nrow(C0datf)
nrow(HR1datf) + nrow(C0datf)

# save as RData for section below
save(HR1datf, C0datf,HR1Ainv, C0Ainv, file = "QG-mouse-trade-off.RData")



rm(list=ls())  # delete all objects, to start fresh in section below
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
#.##################### section #3B - run MCMCglmm models ############################
library(MCMCglmm)

#load data and pedigree objects prepared in section 3A
load(file="QG-mouse-trade-off.RData")

#get data ready for running the models
C0datf$WSTRTymd <- as.Date(C0datf$WSTRTymd, "%Y-%m-%d")
HR1datf$WSTRTymd<- as.Date(HR1datf$WSTRTymd, "%Y-%m-%d")
C0datf$Mnth     <- factor(format(C0datf$WSTRTymd, "%m")) #month variable as factor
HR1datf$Mnth    <- factor(format(HR1datf$WSTRTymd, "%m")) #month variable as factor


# Priors
# multivariate extension of parameter expanded/central-scaled F
## gives FLAT prior on the correlations
k <- 2
kPEflatR <- list(R = list(V = diag(k), nu = 0),
  G = list(G1 = list(V = diag(k)*0.02, nu = k+1,
                     alpha.mu = rep(0, k), alpha.V = diag(k)*1000),
           G2 = list(V = diag(k)*0.02, nu = k+1,
                     alpha.mu = rep(0, k), alpha.V = diag(k)*1000)))

# set number of iterations
nsamp <- 3000
THIN <- 250
BURN <- 10000
(NITT <- BURN + nsamp*THIN)

# CONTROL MICE - run bivariate model
C_RPM_INT_mnth <- MCMCglmm(cbind(RPM56l, INT56l) ~ trait +
    trait:sex + trait:line + trait:GENfac + trait:WHLSTAGE + trait:Fcoeff +
    trait:Mnth,
  random = ~ us(trait):animal + us(trait):damid,
  rcov = ~ us(trait):units,
  data = C0datf,
  ginverse = list(animal = C0Ainv),
  prior = kPEflatR,
  family = rep("gaussian", 2),
  nitt = NITT, thin = THIN, burnin = BURN,
  pr = TRUE, saveX = TRUE, saveZ = TRUE)
summary(C_RPM_INT_mnth)    # see Table S2

# SELECTED MICE  - run bivariate model
# set number of iterations
nsamp <- 3000
THIN <- 500
BURN <- 10000
(NITT <- BURN + nsamp*THIN)

HR_RPM_INT_mnth<- MCMCglmm(cbind(RPM56l, INT56l) ~ trait +
    trait:sex + trait:line + trait:GENfac + trait:WHLSTAGE + trait:Fcoeff +
    trait:Mnth,
  random = ~ us(trait):animal + us(trait):damid,
  rcov = ~ us(trait):units,
  data = HR1datf,
  ginverse = list(animal = HR1Ainv),
  prior = kPEflatR,
  family = rep("gaussian", 2),
  nitt = NITT, thin = THIN, burnin = BURN,
  pr = TRUE, saveX = TRUE, saveZ = TRUE)
summary(HR_RPM_INT_mnth)  # see Table S2



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


# save models as RData
save(C_RPM_INT_mnth, HR_RPM_INT_mnth,
  file = "QG-mouse-trade-off_bivariate_models.RData")
save(C_RUN_mnth, HR_RUN_mnth,      
  file = "QG-mouse-trade-off_univariate_models.RData")







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
load(file="QG-mouse-trade-off_bivariate_models.RData")
#load data (needed to assign generations and lines and breeding values)
load(file="QG-mouse-trade-off.RData")


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
LINES <- c(1, 2, 4, 5)
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
  rfacNms <- sapply(strsplit(rtrms, "\\:"), FUN = tail, 1) #<-- XXX ASSUME only ever 1 ":" used


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
            seq(traitZrowInd[1,x], traitZrowInd[2,x], by = 1)})), ] == 0.0))
      datfx[ncomp, , ] <- tmpDatfx[seq(traitZrowInd[1,t], traitZrowInd[2,t], by = 1), ]
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
            seq(traitZrowInd[1,x], traitZrowInd[2,x], by = 1)})), ] == 0.0))
      datfx[ncomp, , ] <- tmpDatfx[seq(traitZrowInd[1,t], traitZrowInd[2,t], by = 1), ]
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

save("POST.C.LINES", "POST.S.LINES",
  file = "QG-mouse-trade-off_covar_POST.LINES.RData")









rm(list=ls())  # delete all objects, to start fresh section below
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################
#.######### section #3D - make Figure 3 breeding values and genetic correlations ################

library(MCMCglmm)

load(file = "QG-mouse-trade-off.RData")
load(file = "QG-mouse-trade-off_bivariate_models.RData")
load(file = "QG-mouse-trade-off_covar_POST.LINES.RData")



GEN.lst <- c(0:31, 36:51, 53:62, 65, 66, 68:78)  #<-- GENs with data
 # (GEN.lst should match column names of arrays in POST.C.LINES or POST.S.LINES)

#
Ra <- vector("list", length = 8)
for(l in 1:8){
  # Control lines
  if(l < 5){
    Ra[[l]] <- data.frame(GEN = GEN.lst,
      r = posterior.mode(mcmc(POST.C.LINES$rA[, , l])))
  }
  if(l >= 5){
    Ra[[l]] <- data.frame(GEN = GEN.lst,
      r = posterior.mode(mcmc(POST.S.LINES$rA[, , l - 4])))
  }
}  #<-- end for l   


acrossLinePost.C <- as.mcmc(apply(POST.C.LINES$rA, MARGIN = 2,
                                  FUN = function(v){ c(v)}))
acrossLinePost.S <- as.mcmc(apply(POST.S.LINES$rA, MARGIN = 2,
                                  FUN = function(v){ c(v)}))
Ra.C <- data.frame(GEN=GEN.lst, r = posterior.mode(acrossLinePost.C))
  Ra.C[, c("lower", "upper")] <- HPDinterval(acrossLinePost.C)
Ra.S <- data.frame(GEN=GEN.lst, r = posterior.mode(acrossLinePost.S))
  Ra.S[, c("lower", "upper")] <- HPDinterval(acrossLinePost.S)
# Fill in missing generations with blank info  
new_rows <- data.frame(GEN = c(32,33,34,35,52,63,64,67),
  r = NA, lower = NA, upper = NA)
Ra.C <- rbind(Ra.C, new_rows)
Ra.S <- rbind(Ra.S, new_rows)
Ra.C$GEN <- Ra.C$GEN-0.1
Ra.S$GEN <- Ra.S$GEN+0.1
#
Ra.C <- Ra.C[order(Ra.C$GEN), ]
Ra.S <- Ra.S[order(Ra.S$GEN), ]
#
for(l in 1:8){
  Ra[[l]] <- rbind(Ra[[l]], new_rows[, c("GEN", "r")])
  Ra[[l]] <- Ra[[l]][order(Ra[[l]][, "GEN"]),]
}







# get breeding values for running speed and distance
nrow(C0datf)     #N = 10,878
nrow(HR1datf)    #N = 25,124

# Speed and Duration models: extract breeding values

# Control breeding value summary and processing
breedvalCnms <- with(C_RPM_INT_mnth, grep("animal", colnames(Sol)))
#XXX NOTE: posterior.mode and HPDinterval functions can take several seconds
## AND require a decent amount of RAM
MODE.C <- with(C_RPM_INT_mnth, posterior.mode(Sol[, breedvalCnms]))
BLUPS.C <- data.frame(name = names(MODE.C), BLUP = MODE.C)
  BLUPS.C[, c("lower", "upper")] <- with(C_RPM_INT_mnth,
                                         HPDinterval(Sol[, breedvalCnms]))
  nrow(BLUPS.C)
  head(BLUPS.C)
BLUPS.C$TRAIT <- substr(BLUPS.C$name, 6, 8)
BLUPS.C$animal <- sapply(strsplit(as.character(BLUPS.C$name), ".animal."),
  FUN = "[", i = 2)
rownames(BLUPS.C) <- NULL
BLUPS.C$name <- NULL

BLUPS.C.wide <- data.frame(reshape(BLUPS.C, v.names = c("BLUP","lower","upper"),
  idvar = "animal",  timevar = "TRAIT", direction = "wide"))
nrow(BLUPS.C)
nrow(BLUPS.C.wide); nrow(BLUPS.C.wide)*2
DATA.C <- merge(C0datf, BLUPS.C.wide, by = "animal", all.x = TRUE)

# free up some memory
rm(list = c("breedvalCnms", "MODE.C"))


# Selected breeding value summary and processing
breedvalSnms <- with(HR_RPM_INT_mnth, grep("animal", colnames(Sol)))
#XXX NOTE: posterior.mode and HPDinterval functions can take several seconds
## AND require a decent amount of RAM (more so for Selected)
MODE.S <- with(HR_RPM_INT_mnth, posterior.mode(Sol[, breedvalSnms]))
BLUPS.S <- data.frame(name = names(MODE.S), BLUP = MODE.S)
  BLUPS.S[, c("lower", "upper")] <- with(HR_RPM_INT_mnth,
                                          HPDinterval(Sol[, breedvalSnms]))
  nrow(BLUPS.S)
  head(BLUPS.S)
BLUPS.S$TRAIT <- substr(BLUPS.S$name, 6, 8)
BLUPS.S$animal <- sapply(strsplit(as.character(BLUPS.S$name), ".animal."),
  FUN = "[", i = 2)
rownames(BLUPS.S) <- NULL
BLUPS.S$name <- NULL

BLUPS.S.wide <- data.frame(reshape(BLUPS.S, v.names = c("BLUP","lower","upper"),
  idvar = "animal", timevar = "TRAIT", direction = "wide"))
nrow(BLUPS.S)
nrow(BLUPS.S.wide); nrow(BLUPS.S.wide)*2
DATA.S <- merge(HR1datf, BLUPS.S.wide, by = "animal", all.x = TRUE)

# free up some memory
rm(list = c("breedvalSnms", "MODE.S"))


DATA <- rbind(DATA.C,DATA.S)
DATA <- droplevels(DATA)
summary(DATA)


#BLUPs RPM trajectories over generations
LINE.C <- aggregate(cbind(BLUP.RPM, BLUP.INT) ~ GEN + line,
  data = DATA.C, FUN = mean)
LINE.S <- aggregate(cbind(BLUP.RPM, BLUP.INT) ~ GEN + line,
  data = DATA.S, FUN = mean)


# save data objects needed for plots since above can be computationally intensive
save("acrossLinePost.C", "acrossLinePost.S",
    "DATA", "DATA.C", "DATA.S", "LINE.C", "LINE.S", "Ra", "Ra.C", "Ra.S", 
  file = "QG-mouse-trade-off_BVs.RData")




################################################################
# if starting new session or importing objects from just above
load(file = "QG-mouse-trade-off_BVs.RData")
################################################################





# make figure 3 BREEDING VALUES and rA
dev.new(width = 9, height = 6, units = "cm")
par(mfrow = c(2, 2), las = 1, oma = c(4, 4, 1, 1), mar = c(1, 1, 1, 1))
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
layout.show(3)
# BLUPs RPM trajectories over generations
plot(BLUP.RPM ~ GEN, data = subset(LINE.C, line == 1),
  type = "l", col = "blue",
  xlab = "n", xaxt = "n", ylab = "", ylim = c(-0.1, 0.6))
  for(l in c(2, 4, 5)){
    lines(BLUP.RPM ~ GEN, data = subset(LINE.C, line == l), col = "blue")
  }
  for(l in c(3, 6, 7, 8)){
    lines(BLUP.RPM ~ GEN, data = subset(LINE.S, line == l), col = "red")
  }
  abline(h = 0, lty = 3)
  arrows(x0 = 26, y0 = 0.2, x1 = 26, y1 = 0.3, length = 0.1, lwd = 2)
  arrows(x0 = 60, y0 = 0.35, x1 = 60, y1 = 0.45, length = 0.1, lwd = 2)
  mtext("Genetic merit (breeding values)", side = 2, las = 3, line = 3)
  mtext("A running speed", side = 3, adj = 0.015, line = -1.5)
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
  
# BLUPs INT trajectories over generations
  plot(BLUP.INT ~ GEN, data = subset(LINE.C, line == 1),
    type = "l", col = "blue",
    xlab = "" , xaxt = "n", ylab = "", yaxt = "n", ylim = c(-0.1, 0.6))
  for(l in c(2, 4, 5)){
    lines(BLUP.INT ~ GEN, data = subset(LINE.C, line == l), col = "blue")
  }
  for(l in c(3, 6, 7, 8)){
    lines(BLUP.INT ~ GEN, data = subset(LINE.S, line == l), col = "red")
  }
  abline(h = 0, lty = 3)
  axis(2, labels = FALSE)
  mtext("B running duration", side = 3, adj = 0.015, line = -1.5)
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
  
# cross-trait correlation between RPM and INT over generations
plot(r ~ GEN, data = Ra.C,
  col = "blue", pch = 16, cex = 1, type = "o",
  xlab = "", xaxt = "n", xlim = c(1, 77), ylab = "", ylim = c(-0.35, 1))
  points(r ~ GEN, data = Ra.S, col = "red" , pch = 17, cex = 1, type = "o")
  arrows(x0 = Ra.C$GEN, y0 = Ra.C$lower, x1 = Ra.C$GEN, y1 = Ra.C$upper,
    col = "black", code = 3, angle = 90, length = 0)
  arrows(x0 = Ra.S$GEN, y0 = Ra.S$lower, x1 = Ra.S$GEN, y1 = Ra.S$upper,
    col = "black", code = 3, angle = 90, length = 0)
  for(l in c(1, 2, 4, 5)){
    lines(r~GEN, data = Ra[[l]], col = rgb(0, 0, 1, 0.9))
  }
  for(l in c(3, 6, 7, 8)){
    lines(r ~ GEN, data = Ra[[l]], col = rgb(1, 0, 0, 0.9))
  }
  points(r ~ GEN, data = Ra.C, col = "black", pch = 16, cex = 1, type = "o")
  points(r ~ GEN, data = Ra.S, col = "black", pch = 17, cex = 1, type = "o")
  axis(1, at = 0:78, labels = FALSE)
  axis(1, at = seq(0, 78, 6), padj = -0.8)
  abline(h = 0, lty = 2)
  mtext(expression(paste(italic(r)[A]," (±se)")), side = 2, line = 2.5, las = 3,cex = 2)
  mtext("Generation", side = 1, line = 2.5, cex = 1.5)

legend(0, 1.1, legend = c("control (average)", "selected (average)"),
  pch = c(16, 17), col = c("black", "black"), bty = "n")
legend(60,1.1, legend = c("replicate control lines", "replicate selected lines"),
  lty = 1, col = c(rgb(0,0,1,0.9),rgb(1,0,0,0.9)), bty = "n")
arrows(x0 = 26, y0 = 1, x1 = 26, y1 = 0.9, length = 0.1, lwd = 3)
arrows(x0 = 60, y0 = 1, x1 = 60, y1 = 0.9, length = 0.1, lwd = 3)
mtext("C", side = 3, adj = 0.015, line = -1.5)




rm(list = ls())
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
#.########################### section #4 emergence of a trade-off ##################################
library(MCMCglmm)
load(file = "QG-mouse-trade-off_BVs.RData")


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

legend( x="bottomright",
        legend=c("Control mice","Control G trajectories","Selected mice","Selected G trajectories"),
        col = c("blue","grey","red","black"), lwd=2, lty = c(NA,1,NA,1),
        pch = c(16,NA,17,NA), merge=T, bty = "n")
text(-0.05,0.35,"generation 15")
arrows(-0.05,0.33,0.025,0.31,length = 0.1)
text(0.24,0.25,"generation 40")
arrows(0.24,0.26,0.215,0.32,length = 0.1)
text(0.085,0.625,"generation 78")
arrows(0.1,0.61,0.13,0.59,length = 0.1)
mtext("Genetic merit for running speed",   side = 2,line = 2,outer=T,las = 3,cex = 1.5)
mtext("Genetic merit for running duration",side = 1,line = 3,outer=T,las = 1,cex = 1.5)
mtext("A) 78 generations",side = 3,adj = 0.015,line = -1.5)

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
  plot(BLUP.RPM~BLUP.INT,LAST.1,ylim = ylim.1,xlim = xlim.1,ylab = "",xlab = "")
 arrows(LAST.1$BLUP.INT, LAST.1$lower.RPM,LAST.1$BLUP.INT, LAST.1$upper.RPM, col = rgb(0,0,0,0.05),code = 3,angle=90,length = 0)
 arrows(LAST.1$lower.INT,LAST.1$BLUP.RPM, LAST.1$upper.INT,LAST.1$BLUP.RPM, col = rgb(0,0,0,0.05),code = 3,angle=90,length = 0)
points(BLUP.RPM~BLUP.INT,LAST.1, col = line,pch = 16,cex = 1.5)
mtext("B) generation 1",side = 3,adj = 0.015,line = -1.5, cex = 0.8)
#
  plot(BLUP.RPM~BLUP.INT,LAST.2,ylim = ylim.2,xlim = xlim.2,ylab = "",xlab = "")
 arrows(LAST.2$BLUP.INT, LAST.2$lower.RPM,LAST.2$BLUP.INT, LAST.2$upper.RPM, col = rgb(0,0,0,0.05),code = 3,angle=90,length = 0)
 arrows(LAST.2$lower.INT,LAST.2$BLUP.RPM, LAST.2$upper.INT,LAST.2$BLUP.RPM, col = rgb(0,0,0,0.05),code = 3,angle=90,length = 0)
points(BLUP.RPM~BLUP.INT,LAST.2, col = line,pch = 16,cex = 1.5)
mtext("C) generation 15",side = 3,adj = 0.015,line = -1.5, cex = 0.8)
#
  plot(BLUP.RPM~BLUP.INT,LAST.3,ylim = ylim.3,xlim = xlim.3,ylab = "",xlab = "")
 arrows(LAST.3$BLUP.INT, LAST.3$lower.RPM,LAST.3$BLUP.INT, LAST.3$upper.RPM, col = rgb(0,0,0,0.05),code = 3,angle=90,length = 0)
 arrows(LAST.3$lower.INT,LAST.3$BLUP.RPM, LAST.3$upper.INT,LAST.3$BLUP.RPM, col = rgb(0,0,0,0.05),code = 3,angle=90,length = 0)
points(BLUP.RPM~BLUP.INT,LAST.3, col = line,pch = 16,cex = 1.5)
mtext("D) generation 40",side = 3,adj = 0.015,line = -1.5, cex = 0.8)
#
  plot(BLUP.RPM~BLUP.INT,LAST.4,ylim = ylim.4,xlim = xlim.4,ylab = "",xlab = "")
 arrows(LAST.4$BLUP.INT, LAST.4$lower.RPM,LAST.4$BLUP.INT, LAST.4$upper.RPM, col = rgb(0,0,0,0.05),code = 3,angle=90,length = 0)
 arrows(LAST.4$lower.INT,LAST.4$BLUP.RPM, LAST.4$upper.INT,LAST.4$BLUP.RPM, col = rgb(0,0,0,0.05),code = 3,angle=90,length = 0)
points(BLUP.RPM~BLUP.INT,LAST.4, col = line,pch = 16,cex = 1.5)
mtext("E) generation 78",side = 3,adj = 0.015,line = -1.5, cex = 0.8)
#
legend( x="topright",legend=c("line S3","line S6","line S7","line S8"), col = c("black","red","blue","green"),pch = 16,bty = "n")


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
load(file="QG-mouse-trade-off_bivariate_models.RData")
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
load("QG-mouse-trade-off_bivariate_models.RData")
load("QG-mouse-trade-off_univariate_models.RData")
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



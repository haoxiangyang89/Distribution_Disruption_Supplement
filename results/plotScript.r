LBList <- read.csv("./Desktop/Git/disruptionN-1/test/LBList.csv")
png(file = "./Desktop/Git/disruptionN-1/test/LBList.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(LBList$LB96,type = "l",ylim=range(c(0,40000)),xlab = "Iterations",ylab = "Obj Value", 
     main = "Lower Bound", col = "#377EB8", lwd = 3, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(LBList$LB72, col = "#E41A1C", lwd = 3)
lines(LBList$LB48, col = "#4DAF4A", lwd = 3)
lines(LBList$LB24, col = "#984EA3", lwd = 3)
lines(LBList$LB12, col = "#FF7F00", lwd = 3)
legend("topleft",c("T = 96","T = 72","T = 48","T = 24","T = 12"), col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5)
dev.off();

tauList <- read.csv("./Desktop/Git/disruptionN-1/test/tauList.csv")
png(file = "./Desktop/Git/disruptionN-1/test/tauList.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(tauList$tau,tauList$LB,type = "l",ylim=range(c(1580,2800)),xlab = "Recovery Length",ylab = "Obj Value", 
     main = "Cost vs. Recovery Length", col = "#377EB8", lwd = 5, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(tauList$tau, tauList$X1000_UB, col = "#E41A1C", lwd = 5)
legend("topleft",c("Lower bound","Simulated Estimation"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5)

dev.off();

#=============================================================================================
# plot 1: GenAll vs. DOnly
dOnlylb1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_1.csv", header = FALSE)
dOnlylb2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_2.csv", header = FALSE)
dOnlylb3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_3.csv", header = FALSE)

dOnlyt1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_1.csv", header = FALSE)
dOnlyt2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_2.csv", header = FALSE)
dOnlyt3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_3.csv", header = FALSE)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/GenAll.png"
png(file = outString, width= 13, height = 8, units = 'in',res = 300);
par(mfrow=c(2,3));
par(mar = c(5,5,2.5,2.5));
plot(1:20, dOnlylb1$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,30000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
      col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, dOnlylb1$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:20, dOnlylb2$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,38000)), xlab = "Iteration", ylab = "LB", main = "Case 33",
      col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, dOnlylb2$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:20, dOnlylb3$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,50000)), xlab = "Iteration", ylab = "LB", main = "Case 123",
      col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, dOnlylb3$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);

par(mar = c(5,5,1.5,2.5));
plot(1:19, dOnlyt1$V1[1:19], type = "l", xlim = range(c(1,100)), ylim=range(c(0,10000)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, dOnlyt1$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
legend("topleft",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:20, dOnlyt2$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,28000)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, dOnlyt2$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
legend("topleft",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:19, dOnlyt3$V1[1:19], type = "l", xlim = range(c(1,100)), ylim=range(c(0,63000)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, dOnlyt3$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
legend("topleft",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
dev.off();

# plot 2: preGen
outString = "./Desktop/Git/disruptionN-1/test/csvOut/pgFig.png"
png(file = outString, width= 13, height = 8, units = 'in',res = 300);
par(mfrow=c(2,3));
par(mar = c(5,5,2.5,2.5));
plot(1:100, dOnlylb1$V3[1:100], type = "l", ylim=range(c(0,30000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(1:100, dOnlylb1$V2[1:100], col = "#E41A1C", lwd = 2);
legend("bottomright",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:100, dOnlylb2$V3[1:100], type = "l", ylim=range(c(0,38000)), xlab = "Iteration", ylab = "LB", main = "Case 33",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(1:100, dOnlylb2$V2[1:100], col = "#E41A1C", lwd = 2);
legend("bottomright",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:100, dOnlylb3$V3[1:100], type = "l", ylim=range(c(0,50000)), xlab = "Iteration", ylab = "LB", main = "Case 123",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(1:100, dOnlylb3$V2[1:100], col = "#E41A1C", lwd = 2);
legend("bottomright",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);

par(mar = c(5,5,1.5,2.5));
plot(dOnlyt1$V3[1:100], type = "l", ylim=range(c(0,2500)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(dOnlyt1$V2[1:100], col = "#E41A1C", lwd = 2);
legend("topleft",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(dOnlyt2$V3[1:100], type = "l", ylim=range(c(0,5500)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(dOnlyt2$V2[1:100], col = "#E41A1C", lwd = 2);
legend("topleft",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(dOnlyt3$V3[1:100], type = "l", ylim=range(c(0,20000)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(dOnlyt3$V2[1:100], col = "#E41A1C", lwd = 2);
legend("topleft",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
dev.off();

# plot 3: NTest
Nlb1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_1.csv", header = FALSE)
Nlb2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_2.csv", header = FALSE)
Nlb3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_3.csv", header = FALSE)
Nlb4 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_4.csv", header = FALSE)
Nlb5 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_5.csv", header = FALSE)
#Nlb6 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_6.csv", header = FALSE)
#Nlb7 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_7.csv", header = FALSE)

Nt1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_1.csv", header = FALSE)
Nt2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_2.csv", header = FALSE)
Nt3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_3.csv", header = FALSE)
Nt4 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_4.csv", header = FALSE)
Nt5 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_5.csv", header = FALSE)
#Nt6 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_6.csv", header = FALSE)
#Nt7 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_7.csv", header = FALSE)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/NFig.png"
png(file = outString, width= 13, height = 8, units = 'in',res = 300);
par(mfrow=c(2,3));
par(mar = c(5,5,2.5,2.5));
plot(1:500, Nlb1$V1[1:500], type = "l", ylim=range(c(0,30000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, Nlb2$V1[1:100], col = "#E41A1C", lwd = 2, log="x");
lines(1:50, Nlb3$V1[1:50], col = "#4DAF4A", lwd = 2, log="x");
lines(1:33, Nlb4$V1[1:33], col = "#984EA3", lwd = 2, log="x");
lines(1:25, Nlb5$V1[1:25], col = "#3CAEA3", lwd = 2, log="x");
legend("bottomright",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(1:500, Nlb1$V2[1:500], type = "l", ylim=range(c(0,38000)), xlab = "Iteration", ylab = "LB", main = "Case 33",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, Nlb2$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
lines(1:50, Nlb3$V2[1:50], col = "#4DAF4A", lwd = 2, log="x");
lines(1:33, Nlb4$V2[1:33], col = "#984EA3", lwd = 2, log="x");
lines(1:25, Nlb5$V2[1:25], col = "#3CAEA3", lwd = 2, log="x");
legend("bottomright",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(1:500, Nlb1$V3[1:500], type = "l", ylim=range(c(0,50000)), xlab = "Iteration", ylab = "LB", main = "Case 123",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(1:100, Nlb2$V3[1:100], col = "#E41A1C", lwd = 2, log="x");
lines(1:50, Nlb3$V3[1:50], col = "#4DAF4A", lwd = 2, log="x");
lines(1:33, Nlb4$V3[1:33], col = "#984EA3", lwd = 2, log="x");
lines(1:25, Nlb5$V3[1:25], col = "#3CAEA3", lwd = 2, log="x");
legend("bottomright",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(0:500, Nt1$V1[1:501], type = "l", ylim=range(c(0,6000)), xlab = "Iteration", ylab = "Time (sec.)",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(0:100, Nt2$V1[1:101], col = "#E41A1C", lwd = 2, log="x");
lines(0:50, Nt3$V1[1:51], col = "#4DAF4A", lwd = 2, log="x");
lines(0:33, Nt4$V1[1:34], col = "#984EA3", lwd = 2, log="x");
lines(0:25, Nt5$V1[1:26], col = "#3CAEA3", lwd = 2, log="x");
legend("topleft",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(0:500, Nt1$V2[1:501], type = "l", ylim=range(c(0,12000)), xlab = "Iteration", ylab = "Time (sec.)",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(0:100, Nt2$V2[1:101], col = "#E41A1C", lwd = 2, log="x");
lines(0:50, Nt3$V2[1:51], col = "#4DAF4A", lwd = 2, log="x");
lines(0:33, Nt4$V2[1:34], col = "#984EA3", lwd = 2, log="x");
lines(0:25, Nt5$V2[1:26], col = "#3CAEA3", lwd = 2, log="x");
legend("topleft",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(0:500, Nt1$V3[1:501], type = "l", ylim=range(c(0,28000)), xlab = "Iteration", ylab = "Time (sec.)",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(0:100, Nt2$V3[1:101], col = "#E41A1C", lwd = 2, log="x");
lines(0:50, Nt3$V3[1:51], col = "#4DAF4A", lwd = 2, log="x");
lines(0:33, Nt4$V3[1:34], col = "#984EA3", lwd = 2, log="x");
lines(0:25, Nt5$V3[1:26], col = "#3CAEA3", lwd = 2, log="x");
legend("topleft",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);
dev.off()

#=============================================================================================
# plot: GenAll vs. DOnly, LB against time
dOnlylb1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_1.csv", header = FALSE)
dOnlylb2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_2.csv", header = FALSE)
dOnlylb3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_3.csv", header = FALSE)

dOnlyt1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_1.csv", header = FALSE)
dOnlyt2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_2.csv", header = FALSE)
dOnlyt3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_3.csv", header = FALSE)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/GenAll_LB.png"

png(file = outString, width= 13, height = 4.5, units = 'in',res = 300);
par(mfrow=c(1,3));
par(mar = c(5,5,2,2.5));
plot(dOnlyt1$V1[1:20],dOnlylb1$V1[1:20], type = "l", xlim = range(c(80,10000)), ylim=range(c(0,30000)), xlab = "Time (sec.)", ylab = "LB", main = "Case 13", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log = "x");
lines(dOnlyt1$V2[1:100],dOnlylb1$V2[1:100], col = "#E41A1C", lwd = 2, log = "x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(dOnlyt2$V1[1:20],dOnlylb2$V1[1:20], type = "l", xlim = range(c(30,28000)), ylim=range(c(0,38000)), xlab = "Time (sec.)", ylab = "LB", main = "Case 33", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(dOnlyt2$V2[1:100],dOnlylb2$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(dOnlyt3$V1[1:20],dOnlylb3$V1[1:20], type = "l", xlim = range(c(50,63000)), ylim=range(c(0,50000)), xlab = "Time (sec.)", ylab = "LB", main = "Case 123", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(dOnlyt3$V2[1:99],dOnlylb3$V2[1:99], col = "#E41A1C", lwd = 2, log="x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
dev.off();

# plot 2: preGen
outString = "./Desktop/Git/disruptionN-1/test/csvOut/pgFig_LB.png"
png(file = outString, width= 13, height = 4, units = 'in',res = 300);
par(mfrow=c(1,3));

par(mar = c(5,5,2,2.5));
plot(dOnlyt1$V3[1:100], dOnlylb1$V3[1:100], type = "l", xlim=range(c(0,2500)), xlab = "Time (sec.)", ylab = "LB", main = "Case 13", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(dOnlyt1$V2[1:100], dOnlylb1$V2[1:100], col = "#E41A1C", lwd = 2);
legend("bottomright",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(dOnlyt2$V3[1:100], dOnlylb2$V3[1:100], type = "l", xlim=range(c(0,5500)), xlab = "Time (sec.)", ylab = "LB", main = "Case 33", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(dOnlyt2$V2[1:100], dOnlylb2$V2[1:100], col = "#E41A1C", lwd = 2);
legend("bottomright",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(dOnlyt3$V3[1:100], dOnlylb3$V3[1:100], type = "l", xlim=range(c(0,20000)), xlab = "Time (sec.)", ylab = "LB", main = "Case 123", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
lines(dOnlyt3$V2[1:100], dOnlylb3$V2[1:100], col = "#E41A1C", lwd = 2);
legend("bottomright",c("Pre-generated cuts","No pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
dev.off();

# plot 3: NTest
Nlb1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_1.csv", header = FALSE)
Nlb2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_2.csv", header = FALSE)
Nlb3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_3.csv", header = FALSE)
Nlb4 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_4.csv", header = FALSE)
Nlb5 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_5.csv", header = FALSE)
#Nlb6 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_6.csv", header = FALSE)
#Nlb7 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_7.csv", header = FALSE)

Nt1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_1.csv", header = FALSE)
Nt2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_2.csv", header = FALSE)
Nt3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_3.csv", header = FALSE)
Nt4 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_4.csv", header = FALSE)
Nt5 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_5.csv", header = FALSE)
#Nt6 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_6.csv", header = FALSE)
#Nt7 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_7.csv", header = FALSE)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/NFig_LB.png"
png(file = outString, width= 13, height = 8, units = 'in',res = 300);
par(mfrow=c(2,3));
par(mar = c(5,5,2.5,2.5));
plot(0:500, Nt1$V1[1:501], type = "l", ylim=range(c(0,6000)), xlab = "Iteration", ylab = "Time (sec.)", main = "Case 13", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(0:100, Nt2$V1[1:101], col = "#E41A1C", lwd = 2, log="x");
lines(0:50, Nt3$V1[1:51], col = "#4DAF4A", lwd = 2, log="x");
lines(0:33, Nt4$V1[1:34], col = "#984EA3", lwd = 2, log="x");
lines(0:25, Nt5$V1[1:26], col = "#3CAEA3", lwd = 2, log="x");
legend("topleft",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(0:500, Nt1$V2[1:501], type = "l", ylim=range(c(0,12000)), xlab = "Iteration", ylab = "Time (sec.)", main = "Case 33", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(0:100, Nt2$V2[1:101], col = "#E41A1C", lwd = 2, log="x");
lines(0:50, Nt3$V2[1:51], col = "#4DAF4A", lwd = 2, log="x");
lines(0:33, Nt4$V2[1:34], col = "#984EA3", lwd = 2, log="x");
lines(0:25, Nt5$V2[1:26], col = "#3CAEA3", lwd = 2, log="x");
legend("topleft",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(0:500, Nt1$V3[1:501], type = "l", ylim=range(c(0,28000)), xlab = "Iteration", ylab = "Time (sec.)", main = "Case 123", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(0:100, Nt2$V3[1:101], col = "#E41A1C", lwd = 2, log="x");
lines(0:50, Nt3$V3[1:51], col = "#4DAF4A", lwd = 2, log="x");
lines(0:33, Nt4$V3[1:34], col = "#984EA3", lwd = 2, log="x");
lines(0:25, Nt5$V3[1:26], col = "#3CAEA3", lwd = 2, log="x");
legend("topleft",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(Nt1$V1[1:500], Nlb1$V1[1:500], type = "l", ylim=range(c(0,30000)), xlab = "Time (sec.)", ylab = "LB",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(Nt2$V1[1:100], Nlb2$V1[1:100], col = "#E41A1C", lwd = 2, log="x");
lines(Nt3$V1[1:50], Nlb3$V1[1:50], col = "#4DAF4A", lwd = 2, log="x");
lines(Nt4$V1[1:33], Nlb4$V1[1:33], col = "#984EA3", lwd = 2, log="x");
lines(Nt5$V1[1:25], Nlb5$V1[1:25], col = "#3CAEA3", lwd = 2, log="x");
legend("bottomright",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(Nt1$V2[1:500], Nlb1$V2[1:500], type = "l", ylim=range(c(0,38000)), xlab = "Time (sec.)", ylab = "LB",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(Nt2$V2[1:100], Nlb2$V2[1:100], col = "#E41A1C", lwd = 2, log="x");
lines(Nt3$V2[1:50], Nlb3$V2[1:50], col = "#4DAF4A", lwd = 2, log="x");
lines(Nt4$V2[1:33], Nlb4$V2[1:33], col = "#984EA3", lwd = 2, log="x");
lines(Nt5$V2[1:25], Nlb5$V2[1:25], col = "#3CAEA3", lwd = 2, log="x");
legend("bottomright",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(Nt1$V3[1:500], Nlb1$V3[1:500], type = "l", ylim=range(c(0,50000)), xlab = "Time (sec.)", ylab = "LB",
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2, log="x");
lines(Nt2$V3[1:100], Nlb2$V3[1:100], col = "#E41A1C", lwd = 2, log="x");
lines(Nt3$V3[1:50], Nlb3$V3[1:50], col = "#4DAF4A", lwd = 2, log="x");
lines(Nt4$V3[1:33], Nlb4$V3[1:33], col = "#984EA3", lwd = 2, log="x");
lines(Nt5$V3[1:25], Nlb5$V3[1:25], col = "#3CAEA3", lwd = 2, log="x");
legend("bottomright",c(expression(~N^p~" = 1",~N^p~" = 5",~N^p~" = 10",~N^p~" = 15",~N^p~" = 20")), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);
dev.off()

# plot 4: UB plot
dOnlyub1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/ub_1.csv", header = FALSE)
dOnlyub2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/ub_2.csv", header = FALSE)
dOnlyub3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/ub_3.csv", header = FALSE)
redRGB <- col2rgb("#E41A1C");
red_1 <- rgb(redRGB[1]/255, redRGB[2]/255, redRGB[3]/255, alpha=0.7)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/LBUB.png"
png(file = outString, width= 13, height = 4.5, units = 'in',res = 300);
par(mfrow=c(1,3));
par(mar = c(5,5,2,2.5));

plot(1:100,dOnlyub1$V1[1:100], type = "l", xlim = range(c(0,100)), ylim=range(c(0,50000)), xlab = "Iterations", ylab = "Bounds", main = "Case 13", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
points(1:100,dOnlyub1$V2[1:100], col = "#E41A1C", lwd = 2, pch = 16);
legend("bottomright",c("Lower Bound","StatisticalUpper Bound"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);

plot(1:100,dOnlyub2$V1[1:100], type = "l", xlim = range(c(0,100)), ylim=range(c(0,60000)), xlab = "Iterations", ylab = "Bounds", main = "Case 33", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
points(1:100,dOnlyub2$V2[1:100], col = "#E41A1C", lwd = 2, pch = 16);
legend("bottomright",c("Lower Bound","StatisticalUpper Bound"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);

plot(1:100,dOnlyub3$V1[1:100], type = "l", xlim = range(c(0,100)), ylim=range(c(0,60000)), xlab = "Iterations", ylab = "Bounds", main = "Case 123", 
     col = "#377EB8", lwd = 2, cex.main = 2, cex.lab = 2, cex.axis = 2);
points(1:100,dOnlyub3$V2[1:100], col = "#E41A1C", lwd = 2, pch = 16);
legend("bottomright",c("Lower Bound","Statistical Upper Bound"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
dev.off()

#polygon(c(1:100,rev(1:100)), c(dOnlyub1$V3[1:100], rev(dOnlyub1$V4[1:100])), col = red_1, border = NA)

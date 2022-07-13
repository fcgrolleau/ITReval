##### Make Figure 2 from original paper

mu1 <- 50
sig1 <- 15

mu2 <- 20
sig2 <- 5

mu3 <- 40
sig3 <- 5

mu4 <- 50
sig4 <- 20

mu5 <- 10
sig5 <- 30

f1 <- function(x) {exp(-(x-mu1)^2/(2*sig1^2)) /(sig1*sqrt(2*pi)) }
f2 <- function(x) {exp(-(x-mu2)^2/(2*sig2^2)) /(sig2*sqrt(2*pi)) }
f3 <- function(x) {exp(-(x-mu3)^2/(2*sig3^2)) /(sig3*sqrt(2*pi)) }
f4 <- function(x) {exp(-(x-mu4)^2/(2*sig4^2)) /(sig4*sqrt(2*pi)) }
f5 <- function(x) {exp(-(x-mu5)^2/(2*sig5^2)) /(sig5*sqrt(2*pi)) }

f <- function(x) f1(x)
e <- function(x) 70*f5(x)

lwst <- .039
expan <- .00165
expit2 <- function(x) abs(1/(1+exp(-x))-.5+lwst)
se <- function(x) expit2(expan*(50-x)^2)
tau <- function(x) { 10*(f3(x-30)-f2(x-30)) }

r_e <- function(x) { 10*(f2(x)-f4(x)) }
r <- function(x) { ifelse(tau(x)<0,1,0) }
i <- function(x) {.5*log((x+1)/(1-x)) }

p <- function(x) {gamma*(1-abs(r_e(x)))^i(alpha) }

p1 <- function(x) {gam(x)*abs(tau(x))^i(alpha)*(1-se(x))^i(beta) }


asre_f <- function(x) { p(x)*f(x)*tau(x)*r_e(x) }
asre_f0 <- function(x) { p0(x)*f(x)*tau(x)*(ifelse(tau(x)<0,1,0)-e(x)) }
asre_f1 <- function(x) { p1(x)*f(x)*tau(x)*(ifelse(tau(x)<0,1,0)-e(x)) }

are_f1 <- function(x) { f(x)*tau(x)*(ifelse(tau(x)<0,1,0)-e(x)) }

p2 <- function(x) {gamma*abs(tau(x))^i(alpha) }
asre_f2 <- function(x) { p2(x)*f(x)*tau(x)*r_e(x) }

sim<- seq(10,100,by=.1)

wl <- 7
dev.new(width=wl, height=wl*1.1, pointsize=7, noRStudioGD = TRUE)
par(mar=c(4,6,1,1))
par(mfcol=c(4,2))
par(xpd=FALSE)

### plot 1

# A
plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(-1,1), xlab="", ylab="", las=1, cex.axis=1.75)
axis(side=1, at = seq(10,100, by=15), cex.axis=1.75)
abline(h=0)
points(sim, 35*sapply(sim, f), type="l", lwd=2, col="#80796BFF")

points(sim, sapply(sim, e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, tau), type="l", lwd=2, col="#DF8F44FF")

points(c(10,60,60,100), c(1,1,0,0), type="l", lwd=2,, col="#6A6599FF")
points(sim, sapply(sim, function(x) tau(x)+se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)-se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")

par(xpd=TRUE)
are_val <- integrate(are_f1, lower=10, upper=100)$value
legend(85,1, box.lty = 0, bg="white", title=paste0("ARE = ", format(round(are_val, 2), nsmall = 2)),
       legend=c("r(x)", expression(pi(x)), expression(tau(x)), expression(p[X](x)~35)),
       lwd=2, col=c("#6A6599FF", "#00A1D5FF","#DF8F44FF", "#80796BFF"), cex=1.2)

mtext('A', side=3, line=-1.1, at=-5, cex=2)

# B
plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab=expression(rho[paste("rd,", alpha)]^"*"*(x)), las=1, cex.lab=2, cex.axis=1.75)
axis(side=1, at = seq(10,100, by=15), cex.axis=1.75)

jamacol1 <- c("#79AF97FF", "#79C897")[1]
jamacol2 <- c("#B24745FF", "#C84745")[1]
jamacol <- c(jamacol2, jamacol1)

n <- 4
it <- 0
p0 <- function(x) {alpha }

for (k in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- k
        it <- it+1
        lines(sim, sapply(sim, p0), lwd=2, col=jamacol[it])
}

n <- 4
it <- 0
asre_val <- list()
for (k in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- k
        it <- it+1
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
}
asre_val <- sapply(asre_val,c)

legend(10,1.00, box.lty = 0,
       legend=rev(c(expression(paste(alpha, "=1/3 ")),
                    expression(paste(alpha, "=2/3 ")))),
       lwd=2, col=rev(jamacol), cex=1.2)

text(29, 1.03, "AIE", pos=4, cex=1.2)
legend(25,1.001, box.lty = 0, 
       legend=rev(format(round(asre_val,2), nsmall = 2)), xjust = 0, cex=1.2)

mtext('C', side=3, line=-1.1, at=-5, cex=2)

# C

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab=expression(rho[paste("cb,", alpha)]^"*"*(x)), las=1, cex.lab=2, cex.axis=1.75)
axis(side=1, at = seq(10,100, by=15), cex.axis=1.75)

n <- 4
it <- 0
p0 <- function(x) (1-abs(ifelse(tau(x)<0,1,0)-e(x)))^i(alpha)

for (k in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- k
        it <- it+1
        lines(sim, sapply(sim, p0), lwd=2, col=rev(jamacol)[it])
}

n <- 4
it <- 0
asre_val <- list()
for (k in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- k
        it <- it+1
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
}
asre_val <- sapply(asre_val,c)

legend(10,0.30, box.lty = 0,
       legend=c(expression(paste(alpha, "=1/3 ")),
                expression(paste(alpha, "=2/3 "))),
       lwd=2, col=rev(jamacol), cex=1.2)

text(29, 0.33, "AIE", pos=4, cex=1.2)
legend(25,0.305, box.lty = 0,
       legend=format(round(asre_val,2), nsmall = 2), xjust=0, cex=1.2)

mtext('E', side=3, line=-1.1, at=-5, cex=2)

# D

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="x", ylab=expression(rho[paste("cl,", alpha)]^"*"*(x)), las=1, cex.lab=2, cex.axis=1.75)
axis(side=1, at = seq(10,100, by=15), cex.axis=1.75)

n <- 4
it <- 0
p0 <- function(x) as.numeric(sign((tau(x)+qt*se(x))*(tau(x)-qt*se(x)))==1)

it <- 0
for (j in c(0.45, .05) ) {
        qt <- qnorm(1-j/2)
        it <- it+1
        lines(sim, sapply(sim, p0)-(it-1)*.01, lwd=2, col=rev(jamacol)[it])
}

it <- 0
asre_val <- list()
for (k in c(0.45, .05) ) {
        alpha <- k
        it <- it+1
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
}
asre_val <- sapply(asre_val,c)

legend(10,1.00, box.lty = 0,
       legend=c(expression(paste(alpha, "=0.45 ")),
                expression(paste(alpha, "=0.05 "))),
       lwd=2, col=rev(jamacol), cex=1.2)

text(29, 1.03, "AIE", pos=4, cex=1.2)
legend(25,1.001, box.lty = 0,
       legend=format(round(asre_val,2), nsmall = 2), xjust=0, cex=1.2)

mtext('G', side=3, line=-1.1, at=-5, cex=2)


### RIGHT PANELS


mu1 <- 50
sig1 <- 15

mu2 <- 20
sig2 <- 5

mu3 <- 40
sig3 <- 5

mu4 <- 50
sig4 <- 20

mu5 <- 100
sig5 <- 30

f1 <- function(x) {exp(-(x-mu1)^2/(2*sig1^2)) /(sig1*sqrt(2*pi)) }
f2 <- function(x) {exp(-(x-mu2)^2/(2*sig2^2)) /(sig2*sqrt(2*pi)) }
f3 <- function(x) {exp(-(x-mu3)^2/(2*sig3^2)) /(sig3*sqrt(2*pi)) }
f4 <- function(x) {exp(-(x-mu4)^2/(2*sig4^2)) /(sig4*sqrt(2*pi)) }
f5 <- function(x) {exp(-(x-mu5)^2/(2*sig5^2)) /(sig5*sqrt(2*pi)) }

f <- function(x) f1(x)
e <- function(x) 70*f5(x)

lwst <- .039
expan <- .00165
expit2 <- function(x) abs(1/(1+exp(-x))-.5+lwst)
se <- function(x) expit2(expan*(50-x)^2)
tau <- function(x) { 10*(f3(x-30)-f2(x-30)) }

r_e <- function(x) { 10*(f2(x)-f4(x)) }
r <- function(x) { ifelse(tau(x)<0,1,0) }
i <- function(x) {.5*log((x+1)/(1-x)) }

p <- function(x) {gamma*(1-abs(r_e(x)))^i(alpha) }

p1 <- function(x) {gam(x)*abs(tau(x))^i(alpha)*(1-se(x))^i(beta) }


asre_f <- function(x) { p(x)*f(x)*tau(x)*r_e(x) }
asre_f0 <- function(x) { p0(x)*f(x)*tau(x)*(ifelse(tau(x)<0,1,0)-e(x)) }
asre_f1 <- function(x) { p1(x)*f(x)*tau(x)*(ifelse(tau(x)<0,1,0)-e(x)) }

are_f1 <- function(x) { f(x)*tau(x)*(ifelse(tau(x)<0,1,0)-e(x)) }

p2 <- function(x) {gamma*abs(tau(x))^i(alpha) }
asre_f2 <- function(x) { p2(x)*f(x)*tau(x)*r_e(x) }

sim<- seq(10,100,by=.1)

wl <- 7


size = 1

# E
plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(-1,1), xlab="", ylab="", las=1, cex.axis=1.75)
axis(side=1, at = seq(10,100, by=15), cex.axis=1.75)
lines(c(5,100), c(0,0))
points(sim, 35*sapply(sim, f), type="l", lwd=2, col="#80796BFF")
points(sim, sapply(sim, e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, tau), type="l", lwd=2, col="#DF8F44FF")
points(c(10,60,60,100), c(1,1,0,0), type="l", lwd=2,, col="#6A6599FF")
points(sim, sapply(sim, function(x) tau(x)+se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)-se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")

par(xpd=TRUE)
are_val <- integrate(are_f1, lower=10, upper=100)$value
legend(85,1, box.lty = 0, bg="white", title=paste0("ARE = ", format(round(are_val, 2), nsmall = 2)),
       legend=c("r(x)", expression(pi(x)), expression(tau(x)), expression(p[X](x)~35)),
       lwd=2, col=c("#6A6599FF", "#00A1D5FF","#DF8F44FF", "#80796BFF"), cex= 1.2)

mtext('B', side=3, line=-1.1, at=-5, cex=2)

# F
plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab="", las=1, cex.axis=1.75)
axis(side=1, at = seq(10,100, by=15), cex.axis=1.75)

jamacol1 <- c("#79AF97FF", "#79C897")[1]
jamacol2 <- c("#B24745FF", "#C84745")[1]
jamacol <- c(jamacol2, jamacol1)

n <- 4
it <- 0
p0 <- function(x) {alpha }

for (k in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- k
        it <- it+1
        lines(sim, sapply(sim, p0), lwd=2, col=jamacol[it])
}

n <- 4
it <- 0
asre_val <- list()
for (k in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- k
        it <- it+1
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
}
asre_val <- sapply(asre_val,c)

legend(10,1.00, box.lty = 0,
       legend=rev(c(expression(paste(alpha, "=1/3 ")),
                    expression(paste(alpha, "=2/3 ")))),
       lwd=2, col=rev(jamacol), cex=1.2)

text(29, 1.03, "AIE", pos=4, cex=1.2)
legend(25,1.001, box.lty = 0, 
       legend=rev(format(round(asre_val,2), nsmall = 2)), xjust = 0,cex=1.2)

mtext('D', side=3, line=-1.1, at=-5, cex=2)

# G

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab="", las=1, cex.axis=1.75)
axis(side=1, at = seq(10,100, by=15), cex.axis=1.75)

n <- 4
it <- 0
p0 <- function(x) (1-abs(ifelse(tau(x)<0,1,0)-e(x)))^i(alpha)

for (k in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- k
        it <- it+1
        lines(sim, sapply(sim, p0), lwd=2, col=rev(jamacol)[it])
}

n <- 4
it <- 0
asre_val <- list()
for (k in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- k
        it <- it+1
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
}
asre_val <- sapply(asre_val,c)

legend(10,1.00, box.lty = 0,
       legend=c(expression(paste(alpha, "=1/3 ")),
                expression(paste(alpha, "=2/3 "))),
       lwd=2, col=rev(jamacol), cex=1.2)

text(29, 1.03, "AIE", pos=4, cex=1.2)
legend(25,1.001, box.lty = 0,
       legend=format(round(asre_val,2), nsmall = 2), xjust=0, cex=1.2)

mtext('F', side=3, line=-1.1, at=-5, cex=2)

# H

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="x", ylab="", las=1, cex.lab=2, cex.axis=1.75)
axis(side=1, at = seq(10,100, by=15), cex.axis=1.75)

n <- 4
it <- 0
p0 <- function(x) as.numeric(sign((tau(x)+qt*se(x))*(tau(x)-qt*se(x)))==1)

it <- 0
for (j in c(0.45, .05) ) {
        qt <- qnorm(1-j/2)
        it <- it+1
        lines(sim, sapply(sim, p0)-(it-1)*.01, lwd=2, col=rev(jamacol)[it])
}

it <- 0
asre_val <- list()
for (k in c(0.45, .05) ) {
        alpha <- k
        it <- it+1
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
}
asre_val <- sapply(asre_val,c)

legend(10,1.00, box.lty = 0,
       legend=c(expression(paste(alpha, "=0.45 ")),
                expression(paste(alpha, "=0.05 "))),
       lwd=2, col=rev(jamacol), cex=1.2)

text(29, 1.03, "AIE", pos=4, cex=1.2)
legend(25,1.001, box.lty = 0,
       legend=format(round(asre_val,2), nsmall = 2), xjust=0, cex=1.2)

mtext('H', side=3, line=-1.1, at=-5, cex=2)


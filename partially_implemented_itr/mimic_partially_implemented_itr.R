source("algo1.R")

set.seed(4321)

# Load data
mimic_si <- read.csv("~/Desktop/github repos/Emulated-ITR/Python/data/mimic_si_preds.csv")

# Create r(X)
mimic_si$r_old <- as.numeric(mimic_si$SOFA_24hours >11 )

# Add intercept
mimic_si$icpt <- 1

# Specify the variables used in the different models
a_var <- "a1"
r_var <- "r_old"
y_var <- "d60d"
mu_vars <- c("admission_age", "SOFA_24hours")
psi_vars <- c("icpt", "admission_age", "SOFA_24hours", "weight", "bun_k1", "ph_k1", "pot_k1")
rho_vars <- c("icpt", "admission_age", "bun_k1", "ph_k1", "pot_k1")

# Estimate MIG, ARE and AIE along their bootstrap standard errors
res <- boot(mimic_si, estimators_boot,
            a_var=a_var, r_var=r_var, y_var=y_var,
            mu_vars=mu_vars,
            psi_vars=psi_vars,
            rho_vars=rho_vars,
            maxit=15, epsi = 10^-2, verbose=TRUE, 
            R=2, parallel="multicore", ncpus=parallel::detectCores(logical = FALSE))

# Store estimates
estimates <- res$t0

# Store confidence intervals
estimates_ci <- emp_boot_ci(res$t0, res$t, alpha = .05)

##### Make Figure 6 from original paper
library(latex2exp)
dev.new(width=6,height=3,pointsize=9, noRStudioGD = TRUE)
par(mar=c(3,5,2,2)+.1, mgp=c(1.75,0.5,0), tcl=-.4)
xl <- c(-.04, .08)
yl <- c(0.5,3.5)
#
plot(xl,yl,axes=F,type="n",xlab="",ylab="", xaxs="i")
segments(0,par('usr')[3],0,3.5, lty="dotted")
segments(estimates_ci[1,],length(estimates):1,estimates_ci[2,],length(estimates):1, lwd=1.5)
points(estimates,length(estimates):1, pch=15, col="#00A1D5FF", cex=3)
axis(1, at=seq(from=xl[1], to=xl[2], by=.02), labels=seq(from=xl[1], to=xl[2], by=.02), cex.axis=.8, lwd=0.75, lwd.tick=0.75)
mtext("Absolute Risk Difference", side=1, line=1.75, cex=1.2)
#
arrows(c(-.005/2,.005/2), 3.5, 15.5*c(-.005/2,.005/2), 3.5, length = 0.07, xpd=TRUE)
mtext("Favors ITR implementation",side=3, line=-.2, at=-.005/2, font=1, las=1, adj=1, cex=.9)
mtext("Favors no ITR implementation",side=3, line=-.2, at=.005/2, font=1, las=1, adj=0, cex=.9)

mtext(c("MIG",
        "ARE",
        "AIE"), 
      side=2, line=3, at=length(estimates):1, las=1, adj=0, cex=1.2)



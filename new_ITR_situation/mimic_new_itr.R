library(boot)
library(gplots)
library(latex2exp)

source("boot_func_new_itr.R")

set.seed(324)

# Load data
mimic_si <-read.csv("/Users/francois/Desktop/github repos/ITR-evaluation/ITR-evaluation/mimic_si_preds.csv")

# Fit propensity score model
ps_mod <- glm(a1 ~ admission_age + weight + bun_k1 + ph_k1 + pot_k1 + SOFA_24hours + immunosuppressant, 
              data=mimic_si, family = "binomial")

# Get propensity score predictions
mimic_si$ps_hat <-  ps_mod$fitted.values

resamples <- 10 # Bootstrap resamples

# cognitive bias scenario
Lambda_cb <- c()
Alpha_seq_cb <- seq(0, 1-10^-100, length=10)
mean_imp_cb <- c()
boot_ci_cb <- list()
it <- 1
it_boot <- 0

for (i in Alpha_seq_cb ) {
  alpha_cb <- i
  res <- boot(mimic_si, cb_new, R=resamples )
  Lambda_cb <- c(Lambda_cb, res$t0 )
  boot_ci_cb[[it]] <- quantile(res$t, probs = c(.025, .975))
  mimic_si$rho <- with(mimic_si, (1- abs(r - ps_hat) )^(.5*log((alpha_cb+1)/(1-alpha_cb))) )
  mean_imp_cb <- c(mean_imp_cb, mean(mimic_si$rho))
  it_boot <- 0
  it <- it+1
}
boot_ci_cb <- sapply(boot_ci_cb, c)

# confidence level scenario
Lambda_cl <- c()
Alpha_seq_cl <- seq(0, 1, length=10)
mean_imp_cl <- c()
boot_ci_cl <- list()
it <- 1
it_boot <- 0

for (i in Alpha_seq_cl ) {
  alpha_cl <- i
  res <- boot(mimic_si, cl_new, R=resamples )
  Lambda_cl <- c(Lambda_cl, res$t0 )
  boot_ci_cl[[it]] <- quantile(res$t, probs = c(.025, .975))
  mimic_si$rho <- with(mimic_si, (iard - qnorm(1 - alpha_cl/2)*ise)*(iard + qnorm(1 - alpha_cl/2))*ise > 0 )
  mean_imp_cl <- c(mean_imp_cl, mean(mimic_si$rho))
  it_boot <- 0
  it <- it+1
}     
boot_ci_cl <- sapply(boot_ci_cl, c)

ci_width <- 0.004
diamond_size <- 2
lab_size <- 2
xlab_pos <- 3.7
ylab_pos <- 3.4
wl <- 8

##### Make Figure 5 from original paper
dev.new(width=wl, height=wl/2, pointsize=7, noRStudioGD = TRUE)
par(mfcol=c(1,2), mar = c(5.5, 6.2, 2, 2))
plot(NULL, xlim=c(0,1), ylim=round(c(min(c(boot_ci_cl, boot_ci_cb)), max(c(boot_ci_cl, boot_ci_cb))), 3), bty="n", las=1,
     xlab="", ylab="")
mtext("A",side=3, line=0, at=-.15, font=1, las=1, adj=1, cex=2)
title(xlab="", line=xlab_pos, cex.lab=lab_size)
abline(h = 0, lty=1, lwd = 2)

title(xlab=TeX('$\\alpha$'), line=xlab_pos, cex.lab=lab_size)
title(ylab=TeX('$\\widehat{\\Lambda}_{ITE}(r,\\rho^{*}_{"•",\\alpha})$'), line=ylab_pos, cex.lab=lab_size)

plotCI(Alpha_seq_cl+.01, Lambda_cl, ui=boot_ci_cl[2,], li=boot_ci_cl[1,], pch=15, gap=0, cex=diamond_size, sfrac=ci_width, col="#df8f44ff", barcol="black", add=TRUE)
plotCI(Alpha_seq_cb-.02, Lambda_cb, ui=boot_ci_cb[2,], li=boot_ci_cb[1,], pch=16, gap=0, cex=diamond_size, sfrac=ci_width, col="#00a1d5ff", barcol="black", add=TRUE)

legend(.027, -0.024, legend=c(TeX('$\\rho^{*}_{cb,\\alpha}$'), TeX('$\\rho^{*}_{cl,\\alpha}$')),
       col=c("#00a1d5ff", "#df8f44ff"), pch=c(16,15), cex=1.5, pt.cex = 1.5, horiz = TRUE)

plot(NULL, xlim=c(0,1), ylim=round(c(min(c(boot_ci_cl, boot_ci_cb)), max(c(boot_ci_cl, boot_ci_cb))), 3), bty="n", las=1,
     xlab="", ylab="")
mtext("B",side=3, line=0, at=-.15, font=1, las=1, adj=1, cex=2)
title(xlab="", line=xlab_pos, cex.lab=lab_size*.9)
abline(h = 0, lty=1, lwd = 2)

title(xlab="Proportion of patients implementing the new ITR", line=xlab_pos, cex.lab=lab_size*.75)
plotCI(mean_imp_cl-.01, Lambda_cl, ui=boot_ci_cl[2,], li=boot_ci_cl[1,], pch=15, gap=0, cex=diamond_size, sfrac=ci_width, col="#df8f44ff", barcol="black", add=TRUE)
plotCI(mean_imp_cb+.02, Lambda_cb, ui=boot_ci_cb[2,], li=boot_ci_cb[1,], pch=16, gap=0, cex=diamond_size, sfrac=ci_width, col="#00a1d5ff", barcol="black", add=TRUE)


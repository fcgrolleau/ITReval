library(reshape2)
library(ggplot2)
library(hrbrthemes)

setwd("/Users/francois/Desktop/github repos/ITR-evaluation/ITR-evaluation/modfied hme/perhaps_maybe/simuls/results_04_22_FP")
load("end_third_scenario_04_25.RData")

simul_to_results <- function(
true=TrueEstimands_psi_different,
e=estimates_psi_different,
lbs=lb_ci_psi_different,
ubs=ub_ci_psi_different) {
  
mse_ss_row <- matrix(NA, nrow=3, ncol=3)
variances_ss_row <- matrix(NA, nrow=3, ncol=3)
bias_ss_row <- matrix(NA, nrow=3, ncol=3)
relative_bias_ss_row <- matrix(NA, nrow=3, ncol=3)
bias_sq_ss_row <- matrix(NA, nrow=3, ncol=3)
coverage_ss_row <- matrix(NA, nrow=3, ncol=3)
ci_width_ss_row <- matrix(NA, nrow=3, ncol=3)

colnames(mse_ss_row) <- c("MIG", "ARE", "AIE")
colnames(variances_ss_row) <- c("MIG", "ARE", "AIE")
colnames(bias_ss_row) <- c("MIG", "ARE", "AIE")
colnames(relative_bias_ss_row) <- c("MIG", "ARE", "AIE")
colnames(bias_sq_ss_row) <- c("MIG", "ARE", "AIE")
colnames(coverage_ss_row) <- c("MIG", "ARE", "AIE")
colnames(ci_width_ss_row) <- c("MIG", "ARE", "AIE")

for(i in 1:length(e))
{
  est <- e[[i]]
  lb <- lbs[[i]]
  ub <- ubs[[i]]
    
mse <- list(MIG=NA, ARE=NA, AIE=NA)
variances <- list(MIG=NA, ARE=NA, AIE=NA)
bias <- list(MIG=NA, ARE=NA, AIE=NA)
relative_bias <- list(MIG=NA, ARE=NA, AIE=NA)
bias_sq <- list(MIG=NA, ARE=NA, AIE=NA)
coverage <- list(MIG=NA, ARE=NA, AIE=NA)
ci_width <- list(MIG=NA, ARE=NA, AIE=NA)

mse$MIG <- mean((est[,"MIG"] - true[,"True_MIG"] )^2)
variances$MIG <- mean((est[,"MIG"] - mean(est[,"MIG"]))^2)
bias$MIG <- mean(est[,"MIG"]) - true[,"True_MIG"]
relative_bias$MIG <- bias$MIG/true[,"True_MIG"]
bias_sq$MIG <- bias$MIG^2

mse$ARE <- mean((est[,"ARE"] - true[,"True_ARE"] )^2)
variances$ARE <- mean((est[,"ARE"] - mean(est[,"ARE"]))^2)
bias$ARE <- mean(est[,"ARE"]) - true[,"True_ARE"]
relative_bias$ARE <- bias$ARE/true[,"True_ARE"]
bias_sq$ARE <- bias$ARE^2

mse$AIE <- mean((est[,"AIE"] - true[,"True_AIE"] )^2)
variances$AIE <- mean((est[,"AIE"] - mean(est[,"AIE"]))^2)
relative_bias$AIE <- bias$ARE/true[,"True_AIE"]
bias$AIE <- mean(est[,"AIE"]) - true[,"True_AIE"]
bias_sq$AIE <- bias$AIE^2


coverage$MIG <- mean(lb[,"MIG"] <=  true[,"True_MIG"] & true[,"True_MIG"] <= ub[,"MIG"])
coverage$ARE <- mean(lb[,"ARE"] <=  true[,"True_ARE"] & true[,"True_ARE"] <= ub[,"ARE"])
coverage$AIE <- mean(lb[,"AIE"] <=  true[,"True_AIE"] & true[,"True_AIE"] <= ub[,"AIE"])

ci_width$MIG <- mean(ub[,"MIG"] - lb[,"MIG"])
ci_width$ARE <- mean(ub[,"ARE"] - lb[,"ARE"])
ci_width$AIE <- mean(ub[,"AIE"] - lb[,"AIE"])

mse_ss_row[i,] <- sapply(mse, c)
variances_ss_row[i,] <- sapply(variances, c)
bias_ss_row[i,] <- sapply(bias, c)
relative_bias_ss_row[i,] <- sapply(relative_bias, c)
bias_sq_ss_row[i,] <- sapply(bias_sq, c)
coverage_ss_row[i,] <- sapply(coverage, c)
ci_width_ss_row[i,] <- sapply(ci_width, c)

}

return(list(orig_mse=mse, mse=mse_ss_row, variances=variances_ss_row, bias=bias_ss_row, relative_bias=relative_bias_ss_row, bias_sq=bias_sq_ss_row, coverage=coverage_ss_row, ci_width=ci_width_ss_row))

} 

### test table building 

simul_to_results(
  true=TrueEstimands_psi_random,
  e=estimates_psi_random,
  lbs=lb_ci_psi_random,
  ubs=ub_ci_psi_random
)

write.csv(format(round(simul_to_results(
  true=TrueEstimands_psi_similar,
  e=estimates_psi_similar,
  lbs=lb_ci_psi_similar,
  ubs=ub_ci_psi_similar
)$ci_width, 3), nsmall=3), file = "temp.csv")

###

long_formating <- function(TrueEstimands=TrueEstimands_psi_different,
         e=estimates_psi_different,
         lbs=lb_ci_psi_different,
         ubs=ub_ci_psi_different,
         ss=c(200,800,2000),
         scenario="test scenario")
{
temp <- simul_to_results(
  true=TrueEstimands,
  e=e,
  lbs=lbs,
  ubs=ubs) 


plot_dat <- list()

for(i in 1:3)
{  
res <- rbind(
  temp$mse[i,],
  temp$variances[i,],
  temp$bias_sq[i,],
  temp$coverage[i,],
  temp$ci_width[i,])

rownames(res) <- c("mse", "variance", "bias_sq", "coverage", "ci_width")

res <- data.frame(res)
res$metric <- c("mse", "variance", "bias_sq", "coverage", "ci_width")

if(i==1){
res$n <- paste0("n=", ss[1]) }

if(i==2){
  res$n <- paste0("n=", ss[2]) }

if(i==3){
  res$n <- paste0("n=", ss[3]) }

plot_dat[[i]] <- res 
}

plot_dat_long <- data.frame()

for(i in 1:3)
{  
temp <- melt(plot_dat[[i]][c("mse", "bias_sq"), c("MIG", "ARE", "AIE", "metric")],
                id.vars       = c("metric"),
                variable.name = "variable",
                value.name    = "value")

temp$sqrt_value <- sqrt(temp$value)
temp$n <- plot_dat[[i]]$n[1]

plot_dat_long <- rbind(plot_dat_long, temp)
}

plot_dat_long$n = factor(plot_dat_long$n, levels = c(paste0("n=", ss[1]), paste0("n=", ss[2]), paste0("n=", ss[3])) )
plot_dat_long$scenario <- scenario

return(plot_dat_long)
}


final_long <- function(TrueEstimands=list(TrueEstimands_psi_different, TrueEstimands_psi_random, TrueEstimands_psi_similar),
         e=list(estimates_psi_different_narm, estimates_psi_random_narm, estimates_psi_similar_narm),
         lbs=list(lb_ci_psi_different_narm, lb_ci_psi_random_narm, lb_ci_psi_similar_narm),
         ubs=list(ub_ci_psi_different_narm, ub_ci_psi_random_narm, ub_ci_psi_similar_narm),
         ss=c(200,800,2000),
         scenario=list("psi/r different", "psi random", "psi/r similar"))
{
  lon_df <- data.frame()
  for(i in 1:length(e)){
  temp <- long_formating(TrueEstimands[[i]], e[[i]], lbs[[i]], ubs[[i]], ss, scenario[[i]])
  lon_df <- rbind(lon_df, temp)
  }
  return(lon_df)
}

plot_dat_long <- final_long(TrueEstimands=list(TrueEstimands_psi_different, TrueEstimands_psi_random, TrueEstimands_psi_similar),
                            e=list(estimates_psi_different, estimates_psi_random, estimates_psi_similar),
                            lbs=list(lb_ci_psi_different, lb_ci_psi_random, lb_ci_psi_similar),
                            ubs=list(ub_ci_psi_different, ub_ci_psi_random, ub_ci_psi_similar))


plot.simul <- function(plot_dat_long)
{  
ggplot(plot_dat_long, aes(fill=metric, y=sqrt_value, x=variable)) +
  facet_grid(scenario ~ n) + 
  scale_fill_manual("legend", values = c("mse" = "#00A1D5FF", "bias_sq" = "#00468BFF")) + 
  geom_bar(position="identity", stat="identity") +
  theme_light(base_size = 22) +
  theme(legend.position="none", strip.text = element_text(colour = 'black'), strip.background =element_rect(fill="gray90"),
        axis.text.x=element_text(angle=45,hjust=1),
        axis.title.y = element_text(face="bold")) +
  xlab("") +
  ylab("Bias/RMSE")
}
                            
plot_dat_long$scenario <- factor(plot_dat_long$scenario, levels = c("psi random", "psi/r similar", "psi/r different") )
plot_dat_long$scenario <- factor(plot_dat_long$scenario, levels = c("psi/r different", "psi random", "psi/r similar") )
plot_dat_long$variable <- factor(gsub("_", "-", plot_dat_long$variable), level=c("MIG", "ARE", "AIE"))

# Relabel scenarios  
levels(plot_dat_long$scenario) <- c("Scenario A", "Scenario B", "Scenario C")

##### Make Figure 4 from original paper
fig <- plot.simul(plot_dat_long)
fig
ggsave("plot_simul.pdf",  width = 30, height = 30, units = "cm")

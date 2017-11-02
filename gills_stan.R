setwd("C:/Users/tyatabe/OneDrive/Docs/Projects/Jamie's")

# Laptop
setwd('/home/tyatabe/Onedrive/Docs/Projects/Jamies')

##############################################################################################
################################## Descriptive stats #########################################
##############################################################################################
# Read in data
ct <- read.csv("ctvalues.csv")

d <- ct[,c("Sample.no.","Week", "Temp", "time_w", "Cage.ID", "Fish_ID", 
           "mort", "at_risk", "Histopathology.score", "AGD", "PRV", 
           "cyst", "pox", "Ten", "Des")]

summary(d)

d.complete <- d[complete.cases(d$Histopathology.score),]

# Getting the standardized (by two sd) temp, Ct (times -1) and mort values
# NEED TO LOG-TRANSFORM MORT BEFORE SCALING
d.cent <- scale(d.complete[,c(3:4, 10:15)], scale=apply(d.complete[,c(3:4, 10:15)], 2, sd)*2)
d.cent <- scale(d.cent, center=F, scale=c(1, 1, rep(-1, 6)))
d.cent <- as.data.frame(d.cent)

d.cent <- cbind(d.complete[,c(1:2, 5:9)], d.cent)
colnames(d.cent)[7] <- "histo"
d.cent$cage <- ifelse(d.cent$Cage.ID==1, 1,2)
d.cent$sample.no <- as.integer(as.factor(d.cent$Week))

# Creatig a pairwise difference in time matrix (time distance b/w obs)
dist.u <- dist(unique(d.cent$Week), method = "maximum", diag = T, upper = T, p = 2)
dist.m.u <- as.matrix(dist.u)

##############################################################################################
#################### Models for fish gill scores ########################
#######################################################################
# Setting up data
data=list(histo=d.cent$histo+1, AGD=d.cent$AGD, Des=d.cent$Des,
          cyst=d.cent$cyst, PRV=d.cent$PRV, pox=d.cent$pox, Ten=d.cent$Ten, cage=d.cent$cage, 
          sample=d.cent$sample.no, Dmat=dist.m.u, Temp = d.cent$Temp, 
          time_w=d.cent$time_w, N=nrow(d.cent), N_cage=length(unique(d.cent$Cage.ID)), 
          N_sample=length(unique(d.cent$Sample.no.)))

# Stan code is stored in another rcode named gills stancode

# Full model: 45 parameters
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0, bagd=0, bdes=0, bcyst=0, 
              bpox=0, bten=0, bagd_des=0, bagd_cyst=0, bagd_pox=0, bagd_ten=0, 
              a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start





m1.1 <- stan(model_code = stancode1, iter=20000, chains=1, cores=1, init=init,
            warmup=4000, control = list(adapt_delta = 0.95),
            data=data)

print(m1.1, probs=c(0.025, 0.975), pars=c("cutpoints", "bt", "bwater_t", "bagd", "bdes", "bcyst", "bpox", "bten", 
                                          "bagd_des", "bagd_cyst", "bagd_pox", "bagd_ten"))
print(m1.1, probs=c(0.025, 0.975), pars=c("a_sample", "sigma", "etasq", "rhosq"))

plot(m1.1, pars=c("bt", "bwater_t", "bagd", "bdes", "bcyst", "bpox", "bten", 
                  "bagd_des", "bagd_cyst", "bagd_pox", "bagd_ten"))

traceplot(m1.1, pars=c("cutpoints", "bt", "bwater_t", "bagd", "bdes", "bcyst", "bpox", "bten", 
                       "bagd_des", "bagd_cyst", "bagd_pox", "bagd_ten"))
traceplot(m1.1, pars=c("a_sample"))
traceplot(m1.1, pars=c("etasq", "rhosq", "sigma"))

log_lik_m1.1 <- extract_log_lik(m1.1)
(m1.1_loo <- loo(log_lik_m1.1))
(m1.1_waic <- waic(log_lik_m1.1))


#Null model
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.2 <- stan(model_code = stancode2, iter=20000, chains=1, cores=1, init=init,
             warmup=4000, control = list(adapt_delta = 0.95),
             data=data)

log_lik_m1.2 <- extract_log_lik(m1.2)
(m1.2_loo <- loo(log_lik_m1.2))
(m1.2_waic <- waic(log_lik_m1.2))


#AGD, temp, time since FW
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0, bagd=0, 
              a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.3 <- stan(model_code = stancode3, iter=20000, chains=1, cores=1, init=init,
             warmup=4000, control = list(adapt_delta = 0.95),
             data=data)

print(m1.3, probs=c(0.025, 0.975), 
      pars=c("cutpoints", "bt", "bwater_t", "bagd"))

plot(m1.3, pars=c("cutpoints", "bt", "bwater_t", "bagd"))


log_lik_m1.3 <- extract_log_lik(m1.3)
(m1.3_loo <- loo(log_lik_m1.3))
(m1.3_waic <- waic(log_lik_m1.3))

# Desmozoan
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0, bagd=0, bdes=0, 
              a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.4 <- stan(model_code = stancode4, iter=20000, chains=1, cores=1, init=init,
             warmup=4000, control = list(adapt_delta = 0.95),
             data=data)

print(m1.4, probs=c(0.025, 0.975), 
      pars=c("cutpoints", "bt", "bwater_t", "bagd", "bdes"))

plot(m1.4, pars=c("bt", "bwater_t", "bagd", "bdes"))


log_lik_m1.4 <- extract_log_lik(m1.4)
(m1.4_loo <- loo(log_lik_m1.4))
(m1.4_waic <- waic(log_lik_m1.4))


# Model with temp and FW treatment
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0,
              a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.5 <- stan(model_code = stancode5, iter=20000, chains=1, cores=1, init=init,
             warmup=4000, control = list(adapt_delta = 0.95),
             data=data)

print(m1.5, probs=c(0.025, 0.975), 
      pars=c("cutpoints", "bt", "bwater_t"))

plot(m1.5, pars=c("bt", "bwater_t"))


log_lik_m1.5 <- extract_log_lik(m1.5)
(m1.5_loo <- loo(log_lik_m1.5))
(m1.5_waic <- waic(log_lik_m1.5))


# Model with AGD b. c. cysticola, temp, and time since fw
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0, bagd=0, bcyst=0,
              a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.6 <- stan(model_code = stancode6, iter=20000, chains=1, cores=1, init=init,
             warmup=4000, control = list(adapt_delta = 0.95),
             data=data)

print(m1.6, probs=c(0.025, 0.975), 
      pars=c("cutpoints", "bt", "bwater_t", "bagd", "bcyst"))

plot(m1.6, pars=c("bt", "bwater_t", "bagd", "bcyst"))


log_lik_m1.6 <- extract_log_lik(m1.6)
(m1.6_loo <- loo(log_lik_m1.6))
(m1.6_waic <- waic(log_lik_m1.6))

# Main effects model (no interactions between pathogens)
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0, bagd=0, bdes=0, bcyst=0, 
              bpox=0, bten=0, a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.7 <- stan(model_code = stancode7, iter=20000, chains=1, cores=1, init=init,
             warmup=4000, control = list(adapt_delta = 0.95),
             data=data)

print(m1.7, probs=c(0.025, 0.975), 
      pars=c("cutpoints", "bt", "bwater_t", "bagd", "bdes", "bcyst", "bpox", "bten"))

plot(m1.7, pars=c("cutpoints", "bt", "bwater_t", "bagd", "bdes", "bcyst", "bpox", "bten"))


log_lik_m1.7 <- extract_log_lik(m1.7)
(m1.7_loo <- loo(log_lik_m1.7))
(m1.7_waic <- waic(log_lik_m1.7))


# Model with AGD, temp, and time since fw. Just one long chain for inference
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0, bagd=0, bpox=0, bten=0, 
              a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.8 <- stan(model_code = stancode3, iter=16000, chains=1, cores=1, init=init,
             warmup=8000, control = list(adapt_delta = 0.95),
             data=data)

print(m1.8, probs=c(0.025, 0.975), 
      pars=c("cutpoints", "bt", "bwater_t", "bagd"))

plot(m1.8, pars=c("bt", "bwater_t", "bagd"))

log_lik_m1.8 <- extract_log_lik(m1.8)
(m1.8_loo <- loo(log_lik_m1.8))
(m1.8_waic <- waic(log_lik_m1.8))



# Model with AGD, pox, temp, and time since fw
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0, bagd=0, bpox=0, 
              a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.9 <- stan(model_code = stancode9, iter=20000, chains=1, cores=1, init=init,
             warmup=4000, control = list(adapt_delta = 0.95),
             data=data)

print(m1.9, probs=c(0.025, 0.975), 
      pars=c("cutpoints", "bagd", "bpox"))

plot(m1.9, pars=c("bagd", "bpox"))


log_lik_m1.9 <- extract_log_lik(m1.9)
(m1.9_loo <- loo(log_lik_m1.9))
(m1.9_waic <- waic(log_lik_m1.9))


# Model with AGD, ten, temp, and time since fw
n_chains <- 1
start <- list(cutpoints=c(-1, 1, 2), bt=0, bwater_t=0, bagd=0, bten=0, 
              a_sample= rep(0, length(unique(d.cent$Sample.no.))))

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m1.10 <- stan(model_code = stancode10, iter=20000, chains=1, cores=1, init=init,
             warmup=4000, control = list(adapt_delta = 0.95),
             data=data)

print(m1.10, probs=c(0.025, 0.975), 
      pars=c("cutpoints", "bagd", "bten"))

plot(m1.10, pars=c("bagd", "bten"))


log_lik_m1.10 <- extract_log_lik(m1.10)
(m1.10_loo <- loo(log_lik_m1.10))
(m1.10_waic <- waic(log_lik_m1.10))


# Renaming for table (and sanity) purposes
full_model <- m1.1_loo
null_model <- m1.2_loo
agd_model <- m1.3_loo
des_model <- m1.4_loo
temp_fw_model <- m1.5_loo
bcyst_model <- m1.6_loo
main_effects_model <- m1.7_loo
prv_model <- m1.8_loo
pox_model <- m1.9_loo
ten_model <- m1.10_loo

print(compare(ten_model, pox_model, main_effects_model, bcyst_model, temp_fw_model, des_model, agd_model, null_model, 
        full_model), digits=1)

## Posterior predictive checking from top ranking model
# Extract posterior samples
postm1.6 <- extract(m1.6)
phi <- postm1.6$phi
a <- postm1.6$cutpoints
cyst <- postm1.6$bcyst
agd <- postm1.6$bagd
mean(cyst>0)
mean(cyst<0)
mean(agd<0)
mean(agd>0)
# Simulate
library(rethinking)
simul <- matrix(data= rep(NA, 16000*300), ncol=300, nrow=16000)

for(i in 1:nrow(phi)){
  for(k in 1:ncol(phi)){
simul[i,k] <- rordlogit(1, phi[i,k], a[i,])
  }
}
# Subtract 1, as it was added before for modeling
simul <- simul-1
histo <- d.cent$histo

# Histograms
simplehist(simul)
simplehist(histo)
# Descriptive stats
mean(as.vector(simul))
mean(histo)
var(as.vector(simul))
var(histo)
# Proportion of each score
mean(histo==0)
mean(as.vector(simul)==0)
mean(histo==1)
mean(as.vector(simul)==1)
mean(histo==2)
mean(as.vector(simul)==2)
mean(histo==3)
mean(as.vector(simul)==3)

# Using Bayesplot
library(bayesplot)
# Density and histogram
ppc_dens_overlay(histo, simul[1:50, ])
ppc_hist(histo, simul[1:11, ])
# Distribution of test statistics
# Prop 1, 2, 3, and 4
prop_0 <- function(x) mean(x == 0)
prop_1 <- function(x) mean(x == 1)
prop_2 <- function(x) mean(x == 2)
prop_3 <- function(x) mean(x == 3)

# Getting summary statistics for simulated data
prop0 <- apply(simul, 1, prop_0)
prop1 <- apply(simul, 1, prop_1)
prop2 <- apply(simul, 1, prop_2)
prop3 <- apply(simul, 1, prop_3)
mean_score <- apply(simul, 1, mean)
var_score <- apply(simul, 1, var)

# Mean and 95% PI
mean(histo==0)
mean(prop0)
PI(prop0, prob=.95)

mean(histo==1)
mean(prop1)
PI(prop1, prob=.95)

mean(histo==2)
mean(prop2)
PI(prop2, prob=.95)

mean(histo==3)
mean(prop3)
PI(prop3, prob=.95)

# Tweaking with ppc_stat function to add y axis text and ticks

p1 <- ppc_stat(histo, simul, stat = "prop_0", binwidth = 0.001) + yaxis_text() + labs(x = "Probability", y="No of simulations")
p2 <- ppc_stat(histo, simul, stat = "prop_1", binwidth = 0.001) + yaxis_text() + labs(x = "Probability", y="No of simulations")
p3 <- ppc_stat(histo, simul, stat = "prop_2", binwidth = 0.001) + yaxis_text() + labs(x = "Probability", y="No of simulations")
p4 <- ppc_stat(histo, simul, stat = "prop_3", binwidth = 0.001) + yaxis_text() + labs(x = "Probability", y="No of simulations")

# Mean and variance
p5 <- ppc_stat(histo, simul, stat = "mean", binwidth = 0.002) + yaxis_text() + labs(x = "Mean", y="No of simulations")
p6 <- ppc_stat(histo, simul, stat = "var", binwidth = 0.004) + yaxis_text() + labs(x = "Variance", y="No of simulations")


# Using bayesplot_grid
bayesplot_grid(p1, p2, p3, p4, p5, p6, titles=c("Probability of score = 0", "Probability of score = 1", 
                                                "Probability of score = 2", "Probability of score = 3",
                                                "Score mean", "Score variance"), legends=F)
# Storing the plot as a tiff file
tiff(file = "fig4.tiff", width=1500, height=1500, res=300)  # create tiff device
# do the plot
bayesplot_grid(p1, p2, p3, p4, p5, p6, titles=c("Probability of score = 0", "Probability of score = 1", 
                                                "Probability of score = 2", "Probability of score = 3",
                                                "Score mean", "Score variance"), legends=F)

dev.off()                 # return to default device (X11)




## Parameter estimates (dot plot)
postm1.1 <- extract(m1.1)
param <- as.data.frame(postm1.1)
param <- param[,c("bt", "bwater_t", "bagd", "bdes", "bcyst",
                                                    "bpox", "bten", "bagd_des", "bagd_cyst", 
                                                    "bagd_pox", "bagd_ten")]
colnames(param) <- c("Temperature", "Weeks since last FW Tx", "N. perurans", 
                     "D. lepeophtherii", "C. B. cysticola", "Pox virus", "T. maritimum", 
                     "N. perurans X D. lepeophtherii", "N. perurans X C. B. cysticola",
                     "N. perurans X Pox virus", "N. perurans X T. maritimum")
var.labels <- names(param)
med <- apply(param, 2, median)
pi.95 <- apply(param, 2, PI, prob=0.95)
lb <- pi.95[1,]; ub <- pi.95[2,]
pi.80 <- apply(param, 2, PI, prob=0.80)
lb.80 <- pi.80[1,]; ub.80 <- pi.80[2,]

sum.param <- data.frame(factor(var.labels, levels=names(param)),med, lb, ub, lb.80, ub.80)
colnames(sum.param) <- c("param", "Estimate", "lb", "ub", "lb.80", "ub.80")
row.names(sum.param) <- NULL

# Invert order of parameters...IT DOES NOT WORK!!!
sum.param$order <- seq(1:nrow(sum.param))

# Plot
library(lattice)
plot1 <- dotplot(reorder(param, -order)~Estimate, data=sum.param, xlim=c(-max(sum.param$ub)-0.02, max(sum.param$ub)+0.02), 
                 scales=list(y=list(cex=8/12), x=list(cex=8/12)), xlab="", 
                 panel=function(x,y){
                   panel.xyplot(x, y, pch=16, cex=1, col="black", grid=T)
                   panel.segments(sum.param$lb, as.numeric(y), sum.param$ub, as.numeric(y), lty=1, 
                                  col="black", lwd=1)
                   panel.segments(sum.param$lb.80, as.numeric(y), sum.param$ub.80, as.numeric(y), lty=1, 
                                  lwd=2, col="red3")
                   panel.xyplot(x, y, pch=16, cex=1, col="black")
                   panel.abline(v=0, col=1, lty=2)
                 })

# Making fig3
tiff(file = "fig3.tiff", width=1500, height=1500, res=300)  # create tiff device
# do the plot
print(plot1)
dev.off()                 # return to default device (X11)


## Counterfactual predictions
# Posterior of parameter estimates
a <- postm1.1$cutpoints
bt <- postm1.1$bt
bwater_t <- postm1.1$bwater_t
btemp_timew <- postm1.1$btemp_timew
bagd <- postm1.1$bagd 
bdes <- postm1.1$bdes 
bcyst <- postm1.1$bcyst
bprv <- postm1.1$bprv
bpox <- postm1.1$bpox
bten <- postm1.1$bten
bagd_des <- postm1.1$bagd_des
bagd_cyst <- postm1.1$bagd_cyst
bagd_prv <- postm1.1$bagd_prv
bagd_pox <- postm1.1$bagd_pox
bagd_ten <- postm1.1$bagd_ten

# Create sequences of values for gill pathogens' concentrations
ct_agd <- seq(from = min(data$AGD), to=max(data$AGD), length.out = 80)
ct_des.l <- rep(min(data$Des), 80)
ct_des.h <- rep(max(data$Des), 80)

ct_cyst <- seq(from = min(data$cyst), to=max(data$cyst), length.out = 80)
ct_prv <- seq(from = min(data$PRV), to=max(data$PRV), length.out = 80)
ct_pox <- seq(from = min(data$pox), to=max(data$pox), length.out = 80)
ct_ten <- seq(from = min(data$Ten), to=max(data$Ten), length.out = 80)
# Since the other covariates are centered at 0, I just won't include them (equivalent of 
# working with their mean value)
# Sequence of mean week effect (I'll use the mean effect)
week <- rep(mean(as.vector(postm1.1$a_sample)), 80)


# Estimate cumulative log-odds of histo score for AGD at varying levels of other pathogens
# AGD X Des interaction

# Low (min) levels of Des
phi.des.low <- matrix(rep(NA, 80*8000), ncol=80)

for ( i in 1:length(ct_agd) ) {
  phi.des.low[,i] = bagd * ct_agd[i] + bdes * ct_des.l[i] + bagd_des * ct_agd[i] * ct_des.l[i]
  + week[i]
}

# Cumulative probability calculation

prob.des.low <- array(rep(NA, 3*8000*80), dim=c(8000,3,80))


for (i in 1:length(ct_agd)){
  
prob.des.low[,,i] <- pordlogit( 1:3 , phi.des.low[,i] , a , log=FALSE )

}

# Getting the median, 66, and 95% PI

prob.des.low.med <- apply(prob.des.low,c(2,3), median)
prob.des.low.66pi <- apply(prob.des.low,c(2,3), PI, prob=0.66)
prob.des.low.95pi <- apply(prob.des.low,c(2,3), PI, prob=0.95)


# Mean (zero) levels of Des
phi.des.med <- matrix(rep(NA, 80*8000), ncol=80)

for ( i in 1:length(ct_agd) ) {
  phi.des.med[,i] = bagd * ct_agd[i] + week[i]
}

# Cumulative probability calculation

prob.des.med <- array(rep(NA, 3*8000*80), dim=c(8000,3,80))


for (i in 1:length(ct_agd)){
  
  prob.des.med[,,i] <- pordlogit( 1:3 , phi.des.med[,i] , a , log=FALSE )
  
}

# Getting the median, 66, and 95% PI

prob.des.med.med <- apply(prob.des.med,c(2,3), median)
prob.des.med.66pi <- apply(prob.des.med,c(2,3), PI, prob=0.66)
prob.des.med.95pi <- apply(prob.des.med,c(2,3), PI, prob=0.95)


# High (max) levels of Des
phi.des.high <- matrix(rep(NA, 80*8000), ncol=80)

for ( i in 1:length(ct_agd) ) {
  phi.des.high[,i] = bagd * ct_agd[i] + bdes * ct_des.h[i] + bagd_des * ct_agd[i] * ct_des.h[i]
  + week[i]
}

# Cumulative probability calculation

prob.des.high <- array(rep(NA, 3*8000*80), dim=c(8000,3,80))


for (i in 1:length(ct_agd)){
  
  prob.des.high[,,i] <- pordlogit( 1:3 , phi.des.high[,i] , a , log=FALSE )
  
}

# Getting the median, 66, and 95% PI

prob.des.high.med <- apply(prob.des.high,c(2,3), median)
prob.des.high.66pi <- apply(prob.des.high,c(2,3), PI, prob=0.66)
prob.des.high.95pi <- apply(prob.des.high,c(2,3), PI, prob=0.95)


# Plotting

yprob <- rep(1, length(ct_agd))

par(mfrow=c(1,3)) 
plot(ct_agd , yprob , type="n" , xlab="Ct N Perurans" , ylab="probability",
      ylim=c(0,0.5)) 
lines(ct_agd, prob.des.low.med[1,], col=col.alpha(rangi2,1))
lines(ct_agd, prob.des.low.med[2,], col=col.alpha(rangi2,1))
lines(ct_agd, prob.des.low.med[3,], col=col.alpha(rangi2,1))

shade(prob.des.low.95pi[,1,], ct_agd, col=col.alpha("pink",0.2))
shade(prob.des.low.95pi[,2,], ct_agd, col=col.alpha(rangi2,0.2))
shade(prob.des.low.95pi[,3,], ct_agd, col=col.alpha("green",0.2))

shade(prob.des.low.66pi[,1,], ct_agd, col=col.alpha("pink",0.2))
shade(prob.des.low.66pi[,2,], ct_agd, col=col.alpha(rangi2,0.2))
shade(prob.des.low.66pi[,3,], ct_agd, col=col.alpha("green",0.2))


plot(ct_agd , yprob , type="n" , xlab="Ct N Perurans" , ylab="probability",
     ylim=c(0,0.5)) 
lines(ct_agd, prob.des.med.med[1,], col=col.alpha(rangi2,1))
lines(ct_agd, prob.des.med.med[2,], col=col.alpha(rangi2,1))
lines(ct_agd, prob.des.med.med[3,], col=col.alpha(rangi2,1))

shade(prob.des.med.95pi[,1,], ct_agd, col=col.alpha("pink",0.2))
shade(prob.des.med.95pi[,2,], ct_agd, col=col.alpha(rangi2,0.2))
shade(prob.des.med.95pi[,3,], ct_agd, col=col.alpha("green",0.2))

plot(ct_agd , yprob , type="n" , xlab="Ct N Perurans" , ylab="probability",
     ylim=c(0,0.5)) 
lines(ct_agd, prob.des.high.med[1,], col=col.alpha(rangi2,1))
lines(ct_agd, prob.des.high.med[2,], col=col.alpha(rangi2,1))
lines(ct_agd, prob.des.high.med[3,], col=col.alpha(rangi2,1))

shade(prob.des.high.95pi[,1,], ct_agd, col=col.alpha("pink",0.2))
shade(prob.des.high.95pi[,2,], ct_agd, col=col.alpha(rangi2,0.2))
shade(prob.des.high.95pi[,3,], ct_agd, col=col.alpha("green",0.2))


##############################################################################################
#################### Models for cage mortality rate ########################
#######################################################################
# Setting up data for stan model of mortality as a function of Ct values
# Read in data
ct <- read.csv("ctvalues.csv")

d <- ct[,c("Sample.no.","Week", "Temp", "time_w", "Cage.ID", "Fish_ID", 
           "mort", "at_risk", "Histopathology.score", "AGD", "PRV", 
           "cyst", "pox", "Ten", "Des")]

summary(d)

d.complete <- d[complete.cases(d$mort),]

# Getting the standardized (by two sd) temp, Ct (times -1) and mort values
d.cent <- scale(d.complete[,c(3:4, 10:15)], scale=apply(d.complete[,c(3:4, 10:15)], 2, sd)*2)
d.cent <- scale(d.cent, center=F, scale=c(1, 1, rep(-1, 6)))
d.cent <- as.data.frame(d.cent)

d.cent <- cbind(d.complete[,c(1:2, 5:9)], d.cent)
colnames(d.cent)[7] <- "histo"
d.cent$cage <- ifelse(d.cent$Cage.ID==1, 1,2)
d.cent$sample.no <- as.integer(as.factor(d.cent$Week))

# Creatig a pairwise difference in time matrix (time distance b/w obs)
dist.u <- dist(unique(d.cent$Week), method = "maximum", diag = T, upper = T, p = 2)
dist.m.u <- as.matrix(dist.u)


#Creating offset
d.cent$logpop <- log(d.cent$at_risk)
# Creating data set with mean Ct and SE
library(doBy)
serror <- function(x){ sd(x)/sqrt(length(x))}
# Mean and SE of repeated samples
d.mean <- summaryBy(AGD + PRV + cyst + pox + Ten + Des + histo ~ cage + sample.no, 
                    FUN=c(mean, serror), data=d.cent)
colnames(d.mean) <- c("cage", "sample", "AGD", "PRV", "cyst", "pox", "Ten", "Des", "histo", "AGDse", "PRVse", "cystse", "poxse", "Tense", "Desse", "histose")

# Value of non-repeated variables
d.one <- summaryBy(Temp + time_w + mort + at_risk + logpop ~ cage + sample.no, 
                   FUN=min, data=d.cent)
colnames(d.one) <- c("cage", "sample", "Temp", "time_w", "mort", "at_risk", "logpop")

d.merge <- cbind(d.mean, d.one[,-1:-2])


# Creating list for stan

data=list(logpop = d.merge$logpop, mort=d.merge$mort, AGD=d.merge$AGD, Des=d.merge$Des,
          cyst=d.merge$cyst, PRV=d.merge$PRV, pox=d.merge$pox, Ten=d.merge$Ten, cage=d.merge$cage, sample=d.merge$sample,
          Dmat=dist.m.u, Temp = d.merge$Temp, time_w=d.merge$time_w, N = nrow(d.merge), N_sample = length(unique(d.merge$sample)),
          N_cage = length(unique(d.merge$cage)), N_AGD=nrow(d.merge), N_Des=nrow(d.merge),N_cyst=nrow(d.merge),
          N_PRV= nrow(d.merge), N_pox=nrow(d.merge), N_Ten =nrow(d.merge), 
          AGDse = d.merge$AGDse+0.01, PRVse = d.merge$PRVse+0.01, 
          cystse = d.merge$cystse+0.01, poxse = d.merge$poxse+0.01, 
          Tense = d.merge$Tense+0.01, Desse = d.merge$Desse+0.01)
  
# Negative binomial: Full model
n_chains <- 4

start <- list( btemp=0, btime_w=0, bagd=0, bdes=0, bcyst=0, a=0, theta=0.1, etasq=0.1, rhosq=0.1, 
               sigma=0.1, bprv=0, bpox=0, bten=0, bagd_des=0, bagd_cyst=0, bagd_prv=0, bagd_pox=0,
               bagd_ten=0, a_sample= rep(0, length(unique(d.merge$sample))), AGD_est = d.merge$AGD, 
              Des_est = d.merge$Des, cyst_est = d.merge$cyst, PRV_est=d.merge$PRV, 
              pox_est=d.merge$pox, Ten_est=d.merge$Ten)

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start



m.11 <- stan(model_code = stancode1, iter=4000, chains=4, cores=4, init=init,
             warmup=2000, control = list(adapt_delta = 0.95),
            data=data)



traceplot(m.11, pars=c("a", "bagd", "bdes", "bcyst", "bprv", "bpox", 
                       "bten", "bagd_des", "bagd_cyst", "bagd_prv", "bagd_pox", 
                      "bagd_ten","btemp", "btime_w"))
traceplot(m.11, pars=c("a_sample", "etasq", "rhosq", "sigma", "theta"))
traceplot(m.11, pars=c("AGD_est"))
traceplot(m.11, pars=c("Des_est"))
traceplot(m.11, pars=c("cyst_est"))
traceplot(m.11, pars=c("PRV_est"))
traceplot(m.11, pars=c("pox_est"))
traceplot(m.11, pars=c("Ten_est"))


print(m.11, pars=c("a", "bagd", "bdes", "bcyst", "bprv", "bpox", "bten", "bagd_des", 
                   "bagd_cyst", "bagd_prv", "bagd_pox", "bagd_ten","btemp", "btime_w"), 
      probs=c(0.025, 0.975))
print(m.11, pars=c("a_sample", "etasq", "rhosq", "sigma", "theta"), 
      probs=c(0.025, 0.975))

stan_plot(m.11, pars=c("bagd", "bdes", "bcyst", "bprv", "bpox", "bten", "bagd_des", "bagd_cyst", 
                       "bagd_prv", "bagd_pox", "bagd_ten","btemp", "btime_w"))


# Posterior predictive checks
post <- extract(m.11)
simul <- post$y_pred

# mean and variance
mean_mort <- apply(simul, 1, mean)
mean(d.merge$mort)
var_mort <- apply(simul, 1, var)
var(d.merge$mort)
# Prop zeroes
prop_0 <- function(x) mean(x == 0)
prop_0(d.merge$mort)

# Plots
library(bayesplot)
# Descriptive stats
pmean <- ppc_stat(d.merge$mort, simul, stat = "mean", binwidth = 0.002) + yaxis_text() + labs(x = "Mean", y="No of simulations")
pvar <- ppc_stat(d.merge$mort, simul[1:1000,], stat = "var") + yaxis_text() + labs(x = "Variance", y="No of simulations")
prop0 <- ppc_stat(d.merge$mort, simul, stat = "prop_0") + yaxis_text() + labs(x = "Prop 0", y="No of simulations")

# Histograms and density
ppc_hist(d.merge$mort, simul[1:5,])
ppc_dens_overlay(d.merge$mort, simul[1:50, ])


# Zero-inflated NB model...converges like shit. Improve model (rewrite in stan)
n_chains <- 4

start <- list( a=0, a_sample= rep(0, length(unique(d.merge$sample))), b=0, 
               b_sample= rep(0, length(unique(d.merge$sample))), theta=0.1, etasq_a=0.1, 
               etasq_b=0.1, rhosq_a=0.1, rhosq_b=0.1, sigma_a=0.1, sigma_b=0.1)

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m.12 <- stan(model_code = stancode2, iter=4000, chains=4, cores=4, init=init,
             warmup=2000, control = list(adapt_delta = 0.95),
             data=data)

print(m.12,  pars=c("a", "a_sample", "b", "b_sample", "theta", "etasq_a", "etasq_b", "rhosq_a", "rhosq_b",
                    "sigma_a", "sigma_b"), probs=c(0.025, 0.975))

traceplot(m.12, pars=c("a", "a_sample", "sigma_a", "rhosq_a", "etasq_a"))
traceplot(m.12, pars=c("b", "b_sample", "sigma_b", "rhosq_b", "etasq_b"))
pairs(m.112, pars=c("a", "a_cage", "sigma_cage"))



# extract samples from the posterior

postm.12 <- extract(m.12)
# Mental note: this pred is just the rate!!! (expected value), not a random sample from
# ZINB
simul <- postm.12$pred
lambda <- postm.12$lambda
theta <- postm.12$theta
p <- postm.12$p

# Simulating observations from a ZINB
library(ZIM)
kk <- matrix(rep(NA,8000*53), ncol=53)

for (i in 1:53){   
  kk[,i] <- rzinb(1, theta, lambda[,i], omega=p[,i])
}


# mean and variance
mean_mort <- apply(simul, 1, mean)
mean(d.merge$mort)
var_mort <- apply(simul, 1, var)
var(d.merge$mort)
# Prop zeroes
prop_0 <- function(x) mean(x == 0)
prop_0(d.merge$mort)
prop_0(as.vector(simul))
pmean <- ppc_stat(d.merge$mort, simul, stat = "mean") + yaxis_text() + labs(x = "Mean", y="No of simulations")
pvar <- ppc_stat(d.merge$mort, simul[-which.maxn(var_mort, length(var_mort[var_mort>2e4])),], stat = "var") + yaxis_text() + labs(x = "Variance", y="No of simulations")
prop0 <- ppc_stat(d.merge$mort, simul, stat = "prop_0") + yaxis_text() + labs(x = "Prop 0", y="No of simulations")


# Zero-inflated NB model, no Gaussian process
n_chains <- 4

start <- list( a=0, a_sample_raw= rep(0, length(unique(d.merge$sample))), b=0, 
               b_sample_raw= rep(0, length(unique(d.merge$sample))), theta=0.1,
               sigma_a=0.1, sigma_b=0.1)

init <- list()
for ( i in 1:n_chains ) init[[i]] <- start


m.13 <- stan(model_code = stancode3, iter=4000, chains=4, cores=4, init=init,
             warmup=2000, control = list(adapt_delta = 0.95),
             data=data)

print(m.13,  pars=c("a", "a_sample", "b", "b_sample", "theta",
                    "sigma_a", "sigma_b"), probs=c(0.025, 0.975))

traceplot(m.13, pars=c("a", "a_sample", "sigma_a"))
traceplot(m.13, pars=c("b", "b_sample", "sigma_b"))





# Comparing models
log_lik1 <- extract_log_lik(m.11)
(waic1 <- waic(log_lik1))
(loo1 <- loo(log_lik1))

log_lik2 <- extract_log_lik(m.12)
(waic2 <- waic(log_lik2))
(loo2 <- loo(log_lik2))

log_lik3 <- extract_log_lik(m.13)
(waic3 <- waic(log_lik3))
(loo3 <- loo(log_lik3))

print(compare(waic1, waic2, waic3), digits=3)
print(compare(loo1, loo2, loo3), digits=3)

plot(loo1)
plot(loo2)


# extract samples from the posterior

postm.13 <- extract(m.13)
simul <- postm.13$pred


# mean and variance
mean_mort <- apply(simul, 1, mean)
mean(d.merge$mort)
var_mort <- apply(simul, 1, var)
var(d.merge$mort)
# Prop zeroes
prop_0 <- function(x) mean(x == 0)
prop_0(d.merge$mort)
prop_0(as.vector(simul))
pmean <- ppc_stat(d.merge$mort, simul, stat = "mean") + yaxis_text() + labs(x = "Mean", y="No of simulations")
pvar <- ppc_stat(d.merge$mort, simul[-which.maxn(var_mort, length(var_mort[var_mort>2e4])),], stat = "var") + yaxis_text() + labs(x = "Variance", y="No of simulations")
prop0 <- ppc_stat(d.merge$mort, simul, stat = "prop_0") + yaxis_text() + labs(x = "Prop 0", y="No of simulations")



# extracting samples from posterior: weekly rate of mortality
post <- extract(m.112)
rate <- post$lambda

# Plotting
rate.med <- apply(rate, 2, FUN=median)
rate.90 <- apply(rate, 2,FUN=PI ,prob=0.9)

plot(rate.med ~ seq(1:nrow(d.merge)), xlab="Sampling week", ylab="Cage mortality",
     ylim=c(0, max(rate.90)), xaxt = "n")
axis(side = 1, at = seq(1:nrow(d.merge)), lwd=0.5,
     labels = c(unique(d.cent$Week),unique(d.cent$Week)[-26]))
abline(v=26.5)

for ( i in 1:nrow(d.merge) ) {
  ci <- rate.90[,i]
  x <- seq(1:nrow(d.merge))[i]
  lines( c(x,x) , ci)
  points(c(x,x) , ci, cex=0.7, pch=3)
}

points(seq(1:nrow(d.merge)), d.merge$mort, cex=1,pch=20, col=rangi2)


# extracting samples from posterior: AGD Ct
Ten <- post$Ten_est
Ten.med <- apply(Ten, 2, FUN=median)
Ten.90 <- apply(Ten, 2, FUN=PI, prob=0.90)

plot(Ten.med ~ seq(1:nrow(d.merge)), xlab="Sampling week", ylab="Standardized Ct",
     ylim=c(min(Ten.90), max(Ten.90)), xaxt = "n")
axis(side = 1, at = seq(1:nrow(d.merge)), lwd=0.5,
     labels = c(unique(d.cent$Week),unique(d.cent$Week)[-26]))
abline(v=26.5)

for ( i in 1:nrow(d.merge) ) {
  ci <- Ten.90[,i]
  x <- seq(1:nrow(d.merge))[i]
  lines( c(x,x) , ci)
  points(c(x,x) , ci, cex=0.7, pch=3)
}
points(seq(1:nrow(d.merge)), d.merge$Ten, cex=1,pch=20, col=rangi2)

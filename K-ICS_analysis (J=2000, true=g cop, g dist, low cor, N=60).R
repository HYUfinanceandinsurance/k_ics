#install.packages("copula")
library(copula)
#install.packages("agop")
library(agop)
#install.packages("VineCopula")
library(VineCopula)

#### Generation of hypothetical population data for dependent risks   ####

set.seed(10000)
obj.cop <- normalCopula(param = c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10),
                        dim=4, dispstr = "un")
samples.cop <- rCopula(120000, obj.cop)
corr_std    <-  getSigma(obj.cop)

mktlsd <- sqrt(2*(log(500)-6))

pop.risk_life <- qnorm(    samples.cop[,1], mean=500, sd=1000) - 500
pop.risk_pnc  <- qnorm(    samples.cop[,2], mean=500, sd=1000) - 500
pop.risk_cred <- qnorm(    samples.cop[,3], mean=500, sd=1000) - 500
pop.risk_mkt  <- qnorm(    samples.cop[,4], mean=500, sd=1000) - 500
pop.risk_total <- pop.risk_life + pop.risk_pnc + pop.risk_cred + pop.risk_mkt

pop.var_true  <- quantile(pop.risk_total, probs=0.995)
hist(pop.risk_total, breaks = 100, freq=FALSE)
lines(density(pop.risk_total),col="red")
abline(v=pop.var_true , col="blue")

J=120000 / 60 # To reduce randomness in calculation of VaR due to random seeds, samples were extracted to compute VaR 

#### VaR estimation with sample size of 60 (one can replace 60 with any number, say, 60, 240, or 480)####

prd60.vstd <- rep(NA, J)
prd60.vemp <- rep(NA, J)
prd60.varc <- rep(NA, J)
prd60.velp <- rep(NA, J)
prd60.cvine <- rep(NA, J)

prd60.var_std <- 0
prd60.var_emp <- 0
prd60.var_arc <- 0
prd60.var_elp <- 0
prd60.var_CVine <- 0

for (i in 1:J) {
  set.seed(i)
  sam60.cop <- rCopula(60, obj.cop)
  
  sam60.risk_life <- qnorm(    sam60.cop[,1], mean=500, sd=1000) - 500
  sam60.risk_pnc  <- qnorm(    sam60.cop[,2], mean=500, sd=1000) - 500
  sam60.risk_cred <- qnorm(    sam60.cop[,3], mean=500, sd=1000) - 500
  sam60.risk_mkt  <- qnorm(    sam60.cop[,4], mean=500, sd=1000) - 500
  sam60.risk_total <- sam60.risk_life + sam60.risk_pnc + sam60.risk_cred + sam60.risk_mkt
  
  est60.parm_life  <- c(0, sd(sam60.risk_life))
  life_trans       <- -mean(sam60.risk_life)
  trs60.risk_life  <- sam60.risk_life+life_trans 
  
  est60.parm_pnc  <- c(0, sd(sam60.risk_pnc))
  pnc_trans       <- -mean(sam60.risk_pnc)
  trs60.risk_pnc  <- sam60.risk_pnc+pnc_trans 
  
  est60.parm_cred  <- c(0, sd(sam60.risk_cred))
  cred_trans       <- -mean(sam60.risk_cred)
  trs60.risk_cred  <- sam60.risk_cred+cred_trans 
  
  est60.parm_mkt  <- c(0, sd(sam60.risk_mkt))
  mkt_trans       <- -mean(sam60.risk_mkt)
  trs60.risk_mkt  <- sam60.risk_mkt+mkt_trans 
  
  
  prd60.var_marginal <- c(
    qnorm(    0.995, sd   =est60.parm_life[2], mean=est60.parm_life[1]),
    qnorm(    0.995, sd   =est60.parm_pnc[ 2], mean=est60.parm_pnc[ 1]),
    qnorm(    0.995, sd   =est60.parm_cred[2], mean=est60.parm_cred[1]),
    qnorm(    0.995, sd   =est60.parm_mkt[ 2], mean=est60.parm_mkt[ 1]))
  
  psd60.cop <- cbind(
    pnorm(      trs60.risk_life, mean=est60.parm_life[1], sd =est60.parm_life[2]),
    pnorm(      trs60.risk_pnc , mean=est60.parm_pnc[ 1], sd =est60.parm_pnc[ 2]),
    pnorm(      trs60.risk_cred, mean=est60.parm_cred[1], sd =est60.parm_cred[2]),
    pnorm(      trs60.risk_mkt , mean=est60.parm_mkt[ 1], sd =est60.parm_mkt[ 2]))
  
  init.arcparm <- 1/(1-mean(c(corKendall(psd60.cop)[1,2:4], corKendall(psd60.cop)[2,3:4], corKendall(psd60.cop)[3,4])))
  
  prd60.arccopfit <- fitCopula(gumbelCopula(dim=4), psd60.cop, method="mpl", start=init.arcparm)
  
  init.elpparm <- c(cor(psd60.cop)[1,2:4], cor(psd60.cop)[2,3:4], cor(psd60.cop)[3,4])
  
  cvm <- RVineStructureSelect(psd60.cop,
                              c(1,2,3,4,5,6),
                              rotations = F, 
                              selectioncrit = "AIC",
                              indeptest = TRUE, level = 0.05, 
                              type = "RVine")
  
  prd60.elpcopfit <- fitCopula(tCopula(dim=4, dispstr = "un", df=4, df.fixed = TRUE), psd60.cop, method="mpl",
                                start=init.elpparm)
  pld60.cop <- rCopula(20000, prd60.elpcopfit@copula)
  
  pld60.risk_life <- qnorm(    pld60.cop[,1], sd   =est60.parm_life[2], mean=est60.parm_life[1])
  pld60.risk_life <- pld60.risk_life - mean(pld60.risk_life)
  pld60.risk_pnc  <- qnorm(    pld60.cop[,2], sd   =est60.parm_pnc[2], mean=est60.parm_pnc[1])
  pld60.risk_pnc  <- pld60.risk_pnc  - mean(pld60.risk_pnc)
  pld60.risk_cred <- qnorm(    pld60.cop[,3], sd   =est60.parm_cred[2], mean=est60.parm_cred[1])
  pld60.risk_cred <- pld60.risk_cred  - mean(pld60.risk_cred)
  pld60.risk_mkt  <- qnorm(    pld60.cop[,4], sd   =est60.parm_mkt[2], mean=est60.parm_mkt[1])
  pld60.risk_mkt  <- pld60.risk_mkt  - mean(pld60.risk_mkt)
  pld60.risk_total <- pld60.risk_life + pld60.risk_pnc + pld60.risk_cred + pld60.risk_mkt
  
  pad60.cop <- rCopula(20000, prd60.arccopfit@copula)
  
  pad60.risk_life <- qnorm(    pad60.cop[,1], sd   =est60.parm_life[2], mean=est60.parm_life[1])
  pad60.risk_life <- pad60.risk_life - mean(pad60.risk_life)
  pad60.risk_pnc  <- qnorm(    pad60.cop[,2], sd   =est60.parm_pnc[2], mean=est60.parm_pnc[1])
  pad60.risk_pnc  <- pad60.risk_pnc  - mean(pad60.risk_pnc)
  pad60.risk_cred <- qnorm(    pad60.cop[,3], sd   =est60.parm_cred[2], mean=est60.parm_cred[1])
  pad60.risk_cred <- pad60.risk_cred  - mean(pad60.risk_cred)
  pad60.risk_mkt  <- qnorm(    pad60.cop[,4], sd   =est60.parm_mkt[2], mean=est60.parm_mkt[1])
  pad60.risk_mkt  <- pad60.risk_mkt  - mean(pad60.risk_mkt)
  pad60.risk_total <- pad60.risk_life + pad60.risk_pnc + pad60.risk_cred + pad60.risk_mkt
  
  pcv60.cop <- RVineSim(20000, cvm)
  
  pcv60.risk_life <- qnorm(    pcv60.cop[,1], sd   =est60.parm_life[2], mean=est60.parm_life[1])
  pcv60.risk_life <- pcv60.risk_life - mean(pad60.risk_life)
  pcv60.risk_pnc  <- qnorm(    pcv60.cop[,2], sd   =est60.parm_pnc[2], mean=est60.parm_pnc[1])
  pcv60.risk_pnc  <- pcv60.risk_pnc  - mean(pad60.risk_pnc)
  pcv60.risk_cred <- qnorm(    pcv60.cop[,3], sd   =est60.parm_cred[2], mean=est60.parm_cred[1])
  pcv60.risk_cred <- pcv60.risk_cred  - mean(pad60.risk_cred)
  pcv60.risk_mkt  <- qnorm(    pcv60.cop[,4], sd   =est60.parm_mkt[2], mean=est60.parm_mkt[1])
  pcv60.risk_mkt  <- pcv60.risk_mkt  - mean(pcv60.risk_mkt)
  pcv60.risk_total <- pcv60.risk_life + pcv60.risk_pnc + pcv60.risk_cred + pcv60.risk_mkt
  
  
  prd60.vstd[i]  <- sqrt(prd60.var_marginal %*% (diag(4)*0.75+0.25) %*% prd60.var_marginal)
  prd60.vemp[i]  <- quantile(sam60.risk_total, probs=0.995)
  prd60.varc[i]  <- quantile(pad60.risk_total, probs=0.995)
  prd60.velp[i]  <- quantile(pld60.risk_total, probs=0.995)
  prd60.cvine[i] <- quantile(pcv60.risk_total, probs=0.995)
  
  prd60.var_std  <- prd60.var_std + sqrt(prd60.var_marginal %*% (diag(4)*0.75+0.25) %*% prd60.var_marginal)/J
  prd60.var_emp  <- prd60.var_emp + quantile(sam60.risk_total, probs=0.995)/J
  prd60.var_arc  <- prd60.var_arc + quantile(pad60.risk_total, probs=0.995)/J
  prd60.var_elp  <- prd60.var_elp + quantile(pld60.risk_total, probs=0.995)/J
  prd60.var_CVine <- prd60.var_CVine + quantile(pcv60.risk_total, probs=0.995)/J
}

c(prd60.var_std, prd60.var_emp, prd60.var_arc, prd60.var_elp, prd60.var_CVine)
c(prd60.var_std, prd60.var_emp, prd60.var_arc, prd60.var_elp, prd60.var_CVine) - pop.var_true

summary(prd60.vstd - pop.var_true)
summary(prd60.vemp - pop.var_true)
summary(prd60.varc - pop.var_true)
summary(prd60.velp - pop.var_true)
summary(prd60.cvine - pop.var_true)

save.image("K-ICS_analysis (J=2000, true=g cop, g dist, low cor, N=60).RData")

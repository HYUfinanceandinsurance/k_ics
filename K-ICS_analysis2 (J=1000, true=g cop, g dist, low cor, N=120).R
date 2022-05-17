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

J=120000 / 120 # To reduce randomness in calculation of VaR due to random seeds, samples were extracted to compute VaR 

#### VaR estimation with sample size of 120 (one can replace 120 with any number, say, 60, 240, or 480)####

prd120.vstd <- rep(NA, J)
prd120.vemp <- rep(NA, J)
prd120.varc <- rep(NA, J)
prd120.velp <- rep(NA, J)
prd120.cvine <- rep(NA, J)
prd120.c2vine <- rep(NA, J)

prd120.var_std <- 0
prd120.var_emp <- 0
prd120.var_arc <- 0
prd120.var_elp <- 0
prd120.var_CVine <- 0
prd120.var_C2Vine <- 0

for (i in 1:J) {
  set.seed(i)
  sam120.cop <- rCopula(120, obj.cop)
  
  sam120.risk_life <- qnorm(    sam120.cop[,1], mean=500, sd=1000) - 500
  sam120.risk_pnc  <- qnorm(    sam120.cop[,2], mean=500, sd=1000) - 500
  sam120.risk_cred <- qnorm(    sam120.cop[,3], mean=500, sd=1000) - 500
  sam120.risk_mkt  <- qnorm(    sam120.cop[,4], mean=500, sd=1000) - 500
  sam120.risk_total <- sam120.risk_life + sam120.risk_pnc + sam120.risk_cred + sam120.risk_mkt
  
  est120.parm_life  <- c(0, sd(sam120.risk_life))
  life_trans       <- -mean(sam120.risk_life)
  trs120.risk_life  <- sam120.risk_life+life_trans 
  
  est120.parm_pnc  <- c(0, sd(sam120.risk_pnc))
  pnc_trans       <- -mean(sam120.risk_pnc)
  trs120.risk_pnc  <- sam120.risk_pnc+pnc_trans 
  
  est120.parm_cred  <- c(0, sd(sam120.risk_cred))
  cred_trans       <- -mean(sam120.risk_cred)
  trs120.risk_cred  <- sam120.risk_cred+cred_trans 
  
  est120.parm_mkt  <- c(0, sd(sam120.risk_mkt))
  mkt_trans       <- -mean(sam120.risk_mkt)
  trs120.risk_mkt  <- sam120.risk_mkt+mkt_trans 
  
  
  prd120.var_marginal <- c(
    qnorm(    0.995, sd   =est120.parm_life[2], mean=est120.parm_life[1]),
    qnorm(    0.995, sd   =est120.parm_pnc[ 2], mean=est120.parm_pnc[ 1]),
    qnorm(    0.995, sd   =est120.parm_cred[2], mean=est120.parm_cred[1]),
    qnorm(    0.995, sd   =est120.parm_mkt[ 2], mean=est120.parm_mkt[ 1]))
  
  psd120.cop <- cbind(
    pnorm(      trs120.risk_life, mean=est120.parm_life[1], sd =est120.parm_life[2]),
    pnorm(      trs120.risk_pnc , mean=est120.parm_pnc[ 1], sd =est120.parm_pnc[ 2]),
    pnorm(      trs120.risk_cred, mean=est120.parm_cred[1], sd =est120.parm_cred[2]),
    pnorm(      trs120.risk_mkt , mean=est120.parm_mkt[ 1], sd =est120.parm_mkt[ 2]))
  
  init.arcparm <- 1/(1-mean(c(corKendall(psd120.cop)[1,2:4], corKendall(psd120.cop)[2,3:4], corKendall(psd120.cop)[3,4])))
  
  prd120.arccopfit <- fitCopula(gumbelCopula(dim=4), psd120.cop, method="mpl", start=init.arcparm)
  
  init.elpparm <- c(cor(psd120.cop)[1,2:4], cor(psd120.cop)[2,3:4], cor(psd120.cop)[3,4])
  
  cvm <- RVineStructureSelect(psd120.cop,
                              c(1,2,3,4,5,6),
                              rotations = F, 
                              selectioncrit = "AIC",
                              indeptest = TRUE, level = 0.05, 
                              type = "RVine")
  
  c2vm <- RVineStructureSelect(psd120.cop,
                               rotations = T, 
                               selectioncrit = "AIC",
                               indeptest = TRUE, level = 0.05, 
                               type = "RVine")
  
  
  prd120.elpcopfit <- fitCopula(tCopula(dim=4, dispstr = "un", df=4, df.fixed = TRUE), psd120.cop, method="mpl",
                                start=init.elpparm)
  pld120.cop <- rCopula(20000, prd120.elpcopfit@copula)
  
  pld120.risk_life <- qnorm(    pld120.cop[,1], sd   =est120.parm_life[2], mean=est120.parm_life[1])
  pld120.risk_life <- pld120.risk_life - mean(pld120.risk_life)
  pld120.risk_pnc  <- qnorm(    pld120.cop[,2], sd   =est120.parm_pnc[2], mean=est120.parm_pnc[1])
  pld120.risk_pnc  <- pld120.risk_pnc  - mean(pld120.risk_pnc)
  pld120.risk_cred <- qnorm(    pld120.cop[,3], sd   =est120.parm_cred[2], mean=est120.parm_cred[1])
  pld120.risk_cred <- pld120.risk_cred  - mean(pld120.risk_cred)
  pld120.risk_mkt  <- qnorm(    pld120.cop[,4], sd   =est120.parm_mkt[2], mean=est120.parm_mkt[1])
  pld120.risk_mkt  <- pld120.risk_mkt  - mean(pld120.risk_mkt)
  pld120.risk_total <- pld120.risk_life + pld120.risk_pnc + pld120.risk_cred + pld120.risk_mkt
  
  pad120.cop <- rCopula(20000, prd120.arccopfit@copula)
  
  pad120.risk_life <- qnorm(    pad120.cop[,1], sd   =est120.parm_life[2], mean=est120.parm_life[1])
  pad120.risk_life <- pad120.risk_life - mean(pad120.risk_life)
  pad120.risk_pnc  <- qnorm(    pad120.cop[,2], sd   =est120.parm_pnc[2], mean=est120.parm_pnc[1])
  pad120.risk_pnc  <- pad120.risk_pnc  - mean(pad120.risk_pnc)
  pad120.risk_cred <- qnorm(    pad120.cop[,3], sd   =est120.parm_cred[2], mean=est120.parm_cred[1])
  pad120.risk_cred <- pad120.risk_cred  - mean(pad120.risk_cred)
  pad120.risk_mkt  <- qnorm(    pad120.cop[,4], sd   =est120.parm_mkt[2], mean=est120.parm_mkt[1])
  pad120.risk_mkt  <- pad120.risk_mkt  - mean(pad120.risk_mkt)
  pad120.risk_total <- pad120.risk_life + pad120.risk_pnc + pad120.risk_cred + pad120.risk_mkt
  
  # Internal 3 (CVine)
  
  pcv120.cop <- RVineSim(20000, cvm)
  
  pcv120.risk_life <- qnorm(    pcv120.cop[,1], sd   =est120.parm_life[2], mean=est120.parm_life[1])
  pcv120.risk_life <- pcv120.risk_life - mean(pad120.risk_life)
  pcv120.risk_pnc  <- qnorm(    pcv120.cop[,2], sd   =est120.parm_pnc[2], mean=est120.parm_pnc[1])
  pcv120.risk_pnc  <- pcv120.risk_pnc  - mean(pad120.risk_pnc)
  pcv120.risk_cred <- qnorm(    pcv120.cop[,3], sd   =est120.parm_cred[2], mean=est120.parm_cred[1])
  pcv120.risk_cred <- pcv120.risk_cred  - mean(pad120.risk_cred)
  pcv120.risk_mkt  <- qnorm(    pcv120.cop[,4], sd   =est120.parm_mkt[2], mean=est120.parm_mkt[1])
  pcv120.risk_mkt  <- pcv120.risk_mkt  - mean(pcv120.risk_mkt)
  pcv120.risk_total <- pcv120.risk_life + pcv120.risk_pnc + pcv120.risk_cred + pcv120.risk_mkt
  
  # Internal 4 (C2Vine)
  
  pc2v120.cop <- RVineSim(20000, c2vm)
  
  pc2v120.risk_life <- qnorm(    pc2v120.cop[,1], sd   =est120.parm_life[2], mean=est120.parm_life[1])
  pc2v120.risk_life <- pc2v120.risk_life - mean(pad120.risk_life)
  pc2v120.risk_pnc  <- qnorm(    pc2v120.cop[,2], sd   =est120.parm_pnc[2], mean=est120.parm_pnc[1])
  pc2v120.risk_pnc  <- pc2v120.risk_pnc  - mean(pad120.risk_pnc)
  pc2v120.risk_cred <- qnorm(    pc2v120.cop[,3], sd   =est120.parm_cred[2], mean=est120.parm_cred[1])
  pc2v120.risk_cred <- pc2v120.risk_cred  - mean(pad120.risk_cred)
  pc2v120.risk_mkt  <- qnorm(    pc2v120.cop[,4], sd   =est120.parm_mkt[2], mean=est120.parm_mkt[1])
  pc2v120.risk_mkt  <- pc2v120.risk_mkt  - mean(pc2v120.risk_mkt)
  pc2v120.risk_total <- pc2v120.risk_life + pc2v120.risk_pnc + pc2v120.risk_cred + pc2v120.risk_mkt
  
  
  prd120.vstd[i]  <- sqrt(prd120.var_marginal %*% (diag(4)*0.75+0.25) %*% prd120.var_marginal)
  prd120.vemp[i]  <- quantile(sam120.risk_total, probs=0.995)
  prd120.varc[i]  <- quantile(pad120.risk_total, probs=0.995)
  prd120.velp[i]  <- quantile(pld120.risk_total, probs=0.995)
  prd120.cvine[i] <- quantile(pcv120.risk_total, probs=0.995)
  prd120.c2vine[i] <- quantile(pcv120.risk_total, probs=0.995)
  
  prd120.var_std  <- prd120.var_std + sqrt(prd120.var_marginal %*% (diag(4)*0.75+0.25) %*% prd120.var_marginal)/J
  prd120.var_emp  <- prd120.var_emp + quantile(sam120.risk_total, probs=0.995)/J
  prd120.var_arc  <- prd120.var_arc + quantile(pad120.risk_total, probs=0.995)/J
  prd120.var_elp  <- prd120.var_elp + quantile(pld120.risk_total, probs=0.995)/J
  prd120.var_CVine <- prd120.var_CVine + quantile(pcv120.risk_total, probs=0.995)/J
  prd120.var_CVine <- prd120.var_C2Vine + quantile(pcv120.risk_total, probs=0.995)/J
}

c(prd120.var_std, prd120.var_emp, prd120.var_arc, prd120.var_elp, prd120.var_CVine, prd120.var_C2Vine)
c(prd120.var_std, prd120.var_emp, prd120.var_arc, prd120.var_elp, prd120.var_CVine, prd120.var_C2Vine) - pop.var_true

summary(prd120.vstd - pop.var_true)
summary(prd120.vemp - pop.var_true)
summary(prd120.varc - pop.var_true)
summary(prd120.velp - pop.var_true)
summary(prd120.cvine - pop.var_true)
summary(prd120.c2vine - pop.var_true)

save.image("K-ICS_analysis2 (J=1000, true=g cop, g dist, low cor, N=120).RData")

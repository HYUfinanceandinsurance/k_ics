#install.packages("copula")
library(copula)
#install.packages("agop")
library(agop)
#install.packages("VineCopula")
library(VineCopula)

#### Generation of hypothetical population data for dependent risks   ####

set.seed(10000)
obj.cop <- gumbelCopula(param = 1.38, dim=4)
samples.cop <- rCopula(120000, obj.cop)
cor(samples.cop)
mktlsd <- sqrt(2*(log(500)-6))

pop.risk_life <- qnorm(    samples.cop[,1], mean=500, sd=1000) - 500
pop.risk_pnc  <- qpareto2( samples.cop[,2], k=4     , s=1500)  - 500
pop.risk_cred <- qt(       samples.cop[,3], df=4)*1000+500     - 500
pop.risk_mkt  <- exp(qnorm(samples.cop[,4], mean=6, sd=mktlsd))- 500
pop.risk_total <- pop.risk_life + pop.risk_pnc + pop.risk_cred + pop.risk_mkt

pop.var_true  <- quantile(pop.risk_total, probs=0.995)
hist(pop.risk_total, breaks = 100, freq=FALSE)
lines(density(pop.risk_total),col="red")
abline(v=pop.var_true , col="blue")

J=120000 / 240 # To reduce randomness in calculation of VaR due to random seeds, samples were extracted to compute VaR 

#### VaR estimation with sample size of 240 (one can replace 240 with any number, say, 60, 240, or 480)####

prd240.vstd <- rep(NA, J)
prd240.vemp <- rep(NA, J)
prd240.varc <- rep(NA, J)
prd240.velp <- rep(NA, J)
prd240.cvine <- rep(NA, J)

prd240.var_std <- 0
prd240.var_emp <- 0
prd240.var_arc <- 0
prd240.var_elp <- 0
prd240.var_CVine <- 0

for (i in 1:J) {
  set.seed(i)
  sam240.cop <- rCopula(240, obj.cop)
  
  sam240.risk_life <- qnorm(    sam240.cop[,1], mean=500, sd=1000) - 500
  sam240.risk_pnc  <- qpareto2( sam240.cop[,2], k=4     , s=1500)  - 500
  sam240.risk_cred <- qt(       sam240.cop[,3], df=4)*1000+500     - 500
  sam240.risk_mkt  <- exp(qnorm(sam240.cop[,4], mean=6, sd=mktlsd))- 500
  sam240.risk_total <- sam240.risk_life + sam240.risk_pnc + sam240.risk_cred + sam240.risk_mkt
  
  est240.parm_life  <- c(0, sd(sam240.risk_life))
  life_trans       <- -mean(sam240.risk_life)
  trs240.risk_life  <- sam240.risk_life+life_trans 
  
  pnc_trans       <- 1e-16 - min(sam240.risk_pnc)
  trs240.risk_pnc  <- sam240.risk_pnc+pnc_trans 
  pnc_alpha       <- max(3, 2/(1-mean(trs240.risk_pnc)^2/
                                 var(trs240.risk_pnc)))
  est240.parm_pnc  <- c(pnc_alpha, mean(trs240.risk_pnc)*(pnc_alpha-1))
  
  est240.parm_cred  <- c(0, sd(sam240.risk_cred)/sqrt(2))
  cred_trans       <- -mean(sam240.risk_cred)
  trs240.risk_cred  <- sam240.risk_cred+cred_trans 
  
  mkt_trans       <- 1e-16 - min(sam240.risk_mkt)
  trs240.risk_mkt  <- sam240.risk_mkt+mkt_trans 
  est240.parm_mkt  <- c(log(mean(trs240.risk_mkt))-0.5*log(var(trs240.risk_mkt)/mean(trs240.risk_mkt)^2+1),
                        sqrt(log(var( trs240.risk_mkt)/        mean(trs240.risk_mkt)^2+1)))
  
  prd240.var_marginal <- c(
    qnorm(    0.995, sd   =est240.parm_life[2], mean=est240.parm_life[1]),
    qpareto2( 0.995, s    =est240.parm_pnc[ 2], k   =est240.parm_pnc[ 1]) - pnc_trans ,
    qt(       0.995, df=4)*est240.parm_cred[2]      +est240.parm_cred[1] ,
    exp(qnorm(0.995, sd   =est240.parm_mkt[ 2], mean=est240.parm_mkt[ 1]))- mkt_trans )
  
  psd240.cop <- cbind(
    pnorm(      trs240.risk_life, mean=est240.parm_life[1], sd =est240.parm_life[2]),
    ppareto2(   trs240.risk_pnc , k   =est240.parm_pnc[ 1], s  =est240.parm_pnc[ 2]),
    pt(df=4, q=(trs240.risk_cred     - est240.parm_cred[1])  /  est240.parm_cred[2]),
    pnorm( (log(trs240.risk_mkt)     - est240.parm_mkt[ 1])  /  est240.parm_mkt[ 2]))
  
  init.arcparm <- 1/(1-mean(c(corKendall(psd240.cop)[1,2:4], corKendall(psd240.cop)[2,3:4], corKendall(psd240.cop)[3,4])))
  
  prd240.arccopfit <- fitCopula(gumbelCopula(dim=4), psd240.cop, method="mpl", start=init.arcparm)
  
  init.elpparm <- c(cor(psd240.cop)[1,2:4], cor(psd240.cop)[2,3:4], cor(psd240.cop)[3,4],4)
  
  cvm <- RVineStructureSelect(psd240.cop,
                              c(1,2,3,4,5,6),
                              rotations = F, 
                              selectioncrit = "AIC",
                              indeptest = TRUE, level = 0.05, 
                              type = "RVine")
  
  prd240.elpcopfit <- fitCopula(tCopula(dim=4, dispstr = "un"), psd240.cop, method="mpl",
                                start=init.elpparm)
  pld240.cop <- rCopula(20000, prd240.elpcopfit@copula)
  
  pld240.risk_life <- qnorm(    pld240.cop[,1], sd   =est240.parm_life[2], mean=est240.parm_life[1])
  pld240.risk_life <- pld240.risk_life - mean(pld240.risk_life)
  pld240.risk_pnc  <- qpareto2( pld240.cop[,2], s    =est240.parm_pnc[ 2], k   =est240.parm_pnc[ 1]) 
  pld240.risk_pnc  <- pld240.risk_pnc  - mean(pld240.risk_pnc)
  pld240.risk_cred <- qt(       pld240.cop[,3], df=4)*est240.parm_cred[2]      +est240.parm_cred[1] 
  pld240.risk_cred <- pld240.risk_cred  - mean(pld240.risk_cred)
  pld240.risk_mkt  <- exp(qnorm(pld240.cop[,4], sd   =est240.parm_mkt[ 2], mean=est240.parm_mkt[ 1]))
  pld240.risk_mkt  <- pld240.risk_mkt  - mean(pld240.risk_mkt)
  pld240.risk_total <- pld240.risk_life + pld240.risk_pnc + pld240.risk_cred + pld240.risk_mkt
  
  pad240.cop <- rCopula(20000, prd240.arccopfit@copula)
  
  pad240.risk_life <- qnorm(    pad240.cop[,1], sd   =est240.parm_life[2], mean=est240.parm_life[1])
  pad240.risk_life <- pad240.risk_life - mean(pad240.risk_life)
  pad240.risk_pnc  <- qpareto2( pad240.cop[,2], s    =est240.parm_pnc[ 2], k   =est240.parm_pnc[ 1]) 
  pad240.risk_pnc  <- pad240.risk_pnc  - mean(pad240.risk_pnc)
  pad240.risk_cred <- qt(       pad240.cop[,3], df=4)*est240.parm_cred[2]      +est240.parm_cred[1] 
  pad240.risk_cred <- pad240.risk_cred  - mean(pad240.risk_cred)
  pad240.risk_mkt  <- exp(qnorm(pad240.cop[,4], sd   =est240.parm_mkt[ 2], mean=est240.parm_mkt[ 1]))
  pad240.risk_mkt  <- pad240.risk_mkt  - mean(pad240.risk_mkt)
  pad240.risk_total <- pad240.risk_life + pad240.risk_pnc + pad240.risk_cred + pad240.risk_mkt
  
  
  pcv240.cop <- RVineSim(20000, cvm)
  
  pcv240.risk_life <- qnorm(    pcv240.cop[,1], sd   =est240.parm_life[2], mean=est240.parm_life[1])
  pcv240.risk_life <- pcv240.risk_life - mean(pcv240.risk_life)
  pcv240.risk_pnc  <- qpareto2( pcv240.cop[,2], s    =est240.parm_pnc[ 2], k   =est240.parm_pnc[ 1]) 
  pcv240.risk_pnc  <- pcv240.risk_pnc  - mean(pcv240.risk_pnc)
  pcv240.risk_cred <- qt(       pcv240.cop[,3], df=4)*est240.parm_cred[2]      +est240.parm_cred[1] 
  pcv240.risk_cred <- pcv240.risk_cred  - mean(pcv240.risk_cred)
  pcv240.risk_mkt  <- exp(qnorm(pcv240.cop[,4], sd   =est240.parm_mkt[ 2], mean=est240.parm_mkt[ 1]))
  pcv240.risk_mkt <- pcv240.risk_mkt - mean(pcv240.risk_mkt)
  pcv240.risk_total <- pcv240.risk_life + pcv240.risk_pnc + pcv240.risk_cred + pcv240.risk_mkt
  
  prd240.vstd[i]  <- sqrt(prd240.var_marginal %*% (diag(4)*0.75+0.25) %*% prd240.var_marginal)
  prd240.vemp[i]  <- quantile(sam240.risk_total, probs=0.995)
  prd240.varc[i]  <- quantile(pad240.risk_total, probs=0.995)
  prd240.velp[i]  <- quantile(pld240.risk_total, probs=0.995)
  prd240.cvine[i] <- quantile(pcv240.risk_total, probs=0.995)
  
  prd240.var_std  <- prd240.var_std + sqrt(prd240.var_marginal %*% (diag(4)*0.75+0.25) %*% prd240.var_marginal)/J
  prd240.var_emp  <- prd240.var_emp + quantile(sam240.risk_total, probs=0.995)/J
  prd240.var_arc  <- prd240.var_arc + quantile(pad240.risk_total, probs=0.995)/J
  prd240.var_elp  <- prd240.var_elp + quantile(pld240.risk_total, probs=0.995)/J
  prd240.var_CVine <- prd240.var_CVine + quantile(pcv240.risk_total, probs=0.995)/J
}

c(prd240.var_std, prd240.var_emp, prd240.var_arc, prd240.var_elp, prd240.var_CVine)
c(prd240.var_std, prd240.var_emp, prd240.var_arc, prd240.var_elp, prd240.var_CVine) - pop.var_true


summary(prd240.vstd - pop.var_true)
summary(prd240.vemp - pop.var_true)
summary(prd240.varc - pop.var_true)
summary(prd240.velp - pop.var_true)
summary(prd240.cvine - pop.var_true)

save.image("K-ICS_analysis (J=500, true=Gumbel cop, 4 dist, high cor, N=240).RData")

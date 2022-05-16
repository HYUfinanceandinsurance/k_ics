load("K-ICS_analysis (J=1000, true=g cop, g dist, high cor, N=60).RData")
prd60_ggh  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

load("K-ICS_analysis (J=1000, true=g cop, g dist, mid cor, N=60).RData")
prd60_ggm  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

load("K-ICS_analysis (J=1000, true=g cop, g dist, low cor, N=60).RData")
prd60_ggl  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

load("K-ICS_analysis (J=1000, true=Gumbel cop, 4 dist, high cor, N=60).RData")
prd60_a4h  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

load("K-ICS_analysis (J=1000, true=Gumbel cop, 4 dist, mid cor, N=60).RData")
prd60_a4m  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

load("K-ICS_analysis (J=1000, true=Gumbel cop, 4 dist, low cor, N=60).RData")
prd60_a4l  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

load("K-ICS_analysis (J=1000, true=t cop, 4 dist, high cor, N=60).RData")
prd60_t4h  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

load("K-ICS_analysis (J=1000, true=t cop, 4 dist, mid cor, N=60).RData")
prd60_t4m  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

load("K-ICS_analysis (J=1000, true=t cop, 4 dist, low cor, N=60).RData")
prd60_t4l  <- cbind(prd60.vstd, prd60.vemp, prd60.varc, prd60.velp, prd60.cvine, pop.var_true) 

rm(list = ls()[!ls() %in% c("prd60_t4h", "prd60_t4m", "prd60_t4l",
                            "prd60_a4h", "prd60_a4m", "prd60_a4l",
                            "prd60_ggh", "prd60_ggm", "prd60_ggl")])

stable <-
  rbind(c(colMeans(prd60_ggl), colMeans(prd60_ggl)[1:5]-colMeans(prd60_ggl)[6]),
        c(colMeans(prd60_ggm), colMeans(prd60_ggm)[1:5]-colMeans(prd60_ggm)[6]),
        c(colMeans(prd60_ggh), colMeans(prd60_ggh)[1:5]-colMeans(prd60_ggh)[6]),
        c(colMeans(prd60_a4l), colMeans(prd60_a4l)[1:5]-colMeans(prd60_a4l)[6]),
        c(colMeans(prd60_a4m), colMeans(prd60_a4m)[1:5]-colMeans(prd60_a4m)[6]),
        c(colMeans(prd60_a4h), colMeans(prd60_a4h)[1:5]-colMeans(prd60_a4h)[6]),
        c(colMeans(prd60_t4l), colMeans(prd60_t4l)[1:5]-colMeans(prd60_t4l)[6]),
        c(colMeans(prd60_t4m), colMeans(prd60_t4m)[1:5]-colMeans(prd60_t4m)[6]),
        c(colMeans(prd60_t4h), colMeans(prd60_t4h)[1:5]-colMeans(prd60_t4h)[6]))

rownames(stable) <- c("Low", "Medium", "High",
                      "Low", "Medium", "High",
                      "Low", "Medium", "High")

colnames(stable) <- c("Standard", "Empirical", "Archimedian", "Elliptical", "CVine", "True",
                      "Standard", "Empirical", "Archimedian", "Elliptical" "CVine")

round(stable)

labels <- c(rep("Standard"  ,1000), rep("Empirical" ,1000),
            rep("Internal 1",1000), rep("Internal 2",1000),
            rep("Internal 3",1000))
labels <- factor(labels, levels = c("Standard", "Empirical", "Internal 1", "Internal 2", "Internal 3"))

library(ggplot2)

prt60_ggh     <- as.data.frame(cbind(as.vector(prd60_ggh[,1:5]-prd60_ggh[,6]),
                                      rep(1, 1000)))
prt60_ggh[,2] <- labels
colnames(prt60_ggh) <- c("Diff", "Label")

plt60_ggh <- ggplot(prt60_ggh, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 3")

plt60_ggh


prt60_ggm     <- as.data.frame(cbind(as.vector(prd60_ggm[,1:5]-prd60_ggm[,6]),
                                      rep(1, 1000)))
prt60_ggm[,2] <- labels
colnames(prt60_ggm) <- c("Diff", "Label")

plt60_ggm <- ggplot(prt60_ggm, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 2")

plt60_ggm


prt60_ggl     <- as.data.frame(cbind(as.vector(prd60_ggl[,1:5]-prd60_ggl[,6]),
                                      rep(1, 1000)))
prt60_ggl[,2] <- labels
colnames(prt60_ggl) <- c("Diff", "Label")

plt60_ggl <- ggplot(prt60_ggl, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 1")

plt60_ggl


prt60_a4h     <- as.data.frame(cbind(as.vector(prd60_a4h[,1:5]-prd60_a4h[,6]),
                                      rep(1, 1000)))
prt60_a4h[,2] <- labels
colnames(prt60_a4h) <- c("Diff", "Label")

plt60_a4h <- ggplot(prt60_a4h, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 6")

plt60_a4h

prt60_a4m     <- as.data.frame(cbind(as.vector(prd60_a4m[,1:5]-prd60_a4m[,6]),
                                      rep(1, 1000)))
prt60_a4m[,2] <- labels
colnames(prt60_a4m) <- c("Diff", "Label")

plt60_a4m <- ggplot(prt60_a4m, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 5")

plt60_a4m


prt60_a4l     <- as.data.frame(cbind(as.vector(prd60_a4l[,1:5]-prd60_a4l[,6]),
                                      rep(1, 1000)))
prt60_a4l[,2] <- labels
colnames(prt60_a4l) <- c("Diff", "Label")

plt60_a4l <- ggplot(prt60_a4l, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 4")

plt60_a4l

prt60_t4h     <- as.data.frame(cbind(as.vector(prd60_t4h[,1:5]-prd60_t4h[,6]),
                                      rep(1, 1000)))
prt60_t4h[,2] <- labels
colnames(prt60_t4h) <- c("Diff", "Label")

plt60_t4h <- ggplot(prt60_t4h, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 9")

plt60_t4h

prt60_t4m     <- as.data.frame(cbind(as.vector(prd60_t4m[,1:5]-prd60_t4m[,6]),
                                      rep(1, 1000)))
prt60_t4m[,2] <- labels
colnames(prt60_t4m) <- c("Diff", "Label")

plt60_t4m <- ggplot(prt60_t4m, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 8")

plt60_t4m


prt60_t4l     <- as.data.frame(cbind(as.vector(prd60_t4l[,1:5]-prd60_t4l[,6]),
                                      rep(1, 1000)))
prt60_t4l[,2] <- labels
colnames(prt60_t4l) <- c("Diff", "Label")

plt60_t4l <- ggplot(prt60_t4l, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 7")

plt60_t4l

library(gridExtra)

grid.arrange(plt60_ggl, plt60_ggm, plt60_ggh,
             plt60_a4l, plt60_a4m, plt60_a4h,
             plt60_t4l, plt60_t4m, plt60_t4h, ncol=3)

save.image("K-ICS_summary (N=60).RData")

load("K-ICS_analysis2 (J=1000, true=g cop, g dist, high cor, N=120).RData")
prd120_ggh  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=1000, true=g cop, g dist, mid cor, N=120).RData")
prd120_ggm  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=1000, true=g cop, g dist, low cor, N=120).RData")
prd120_ggl  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=1000, true=Gumbel cop, 4 dist, high cor, N=120).RData")
prd120_a4h  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=1000, true=Gumbel cop, 4 dist, mid cor, N=120).RData")
prd120_a4m  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=1000, true=Gumbel cop, 4 dist, low cor, N=120).RData")
prd120_a4l  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=1000, true=t cop, 4 dist, high cor, N=120).RData")
prd120_t4h  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=1000, true=t cop, 4 dist, mid cor, N=120).RData")
prd120_t4m  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=1000, true=t cop, 4 dist, low cor, N=120).RData")
prd120_t4l  <- cbind(prd120.vstd, prd120.vemp, prd120.varc, prd120.velp, prd120.cvine, prd120.c2vine, pop.var_true) 


rm(list = ls()[!ls() %in% c("prd120_t4h", "prd120_t4m", "prd120_t4l",
                            "prd120_a4h", "prd120_a4m", "prd120_a4l",
                            "prd120_ggh", "prd120_ggm", "prd120_ggl")])
stable <-
  rbind(c(colMeans(prd120_ggl), colMeans(prd120_ggl)[1:6]-colMeans(prd120_ggl)[7]),
        c(colMeans(prd120_ggm), colMeans(prd120_ggm)[1:6]-colMeans(prd120_ggm)[7]),
        c(colMeans(prd120_ggh), colMeans(prd120_ggh)[1:6]-colMeans(prd120_ggh)[7]),
        c(colMeans(prd120_a4l), colMeans(prd120_a4l)[1:6]-colMeans(prd120_a4l)[7]),
        c(colMeans(prd120_a4m), colMeans(prd120_a4m)[1:6]-colMeans(prd120_a4m)[7]),
        c(colMeans(prd120_a4h), colMeans(prd120_a4h)[1:6]-colMeans(prd120_a4h)[7]),
        c(colMeans(prd120_t4l), colMeans(prd120_t4l)[1:6]-colMeans(prd120_t4l)[7]),
        c(colMeans(prd120_t4m), colMeans(prd120_t4m)[1:6]-colMeans(prd120_t4m)[7]),
        c(colMeans(prd120_t4h), colMeans(prd120_t4h)[1:6]-colMeans(prd120_t4h)[7]))

rownames(stable) <- c("1", "2", "3",
                      "4", "5", "6",
                      "7", "8", "9")
stable
colnames(stable) <- c("Standard", "Empirical", "Archimedian", "Elliptical", "RVine",  "R2Vine","True",
                      "Standard", "Empirical", "Archimedian", "Elliptical", "RVine","R2Vine")

tableA <- stable[,1:7]
tableB <- stable[,8:13]

round(stable)

labels <- c(rep("Standard"  ,1000), rep("Empirical" ,1000),
            rep("Internal 1",1000), rep("Internal 2",1000),
            rep("Internal 3",1000), rep("Internal 4",1000))
labels <- factor(labels, levels = c("Standard", "Empirical", "Internal 1", "Internal 2", "Internal 3", "Internal 4"))

library(ggplot2)

prt120_ggh     <- as.data.frame(cbind(as.vector(prd120_ggh[,1:6]-prd120_ggh[,7]),
                                      rep(1, 1000)))
prt120_ggh[,2] <- labels
colnames(prt120_ggh) <- c("Diff", "Label")

plt120_ggh <- ggplot(prt120_ggh, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 3")

plt120_ggh


prt120_ggm     <- as.data.frame(cbind(as.vector(prd120_ggm[,1:6]-prd120_ggm[,7]),
                                      rep(1, 1000)))
prt120_ggm[,2] <- labels
colnames(prt120_ggm) <- c("Diff", "Label")

plt120_ggm <- ggplot(prt120_ggm, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 2")

plt120_ggm


prt120_ggl     <- as.data.frame(cbind(as.vector(prd120_ggl[,1:6]-prd120_ggl[,7]),
                                      rep(1, 1000)))
prt120_ggl[,2] <- labels
colnames(prt120_ggl) <- c("Diff", "Label")

plt120_ggl <- ggplot(prt120_ggl, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 1")

plt120_ggl


prt120_a4h     <- as.data.frame(cbind(as.vector(prd120_a4h[,1:6]-prd120_a4h[,7]),
                                      rep(1, 1000)))
prt120_a4h[,2] <- labels
colnames(prt120_a4h) <- c("Diff", "Label")

plt120_a4h <- ggplot(prt120_a4h, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 6")

plt120_a4h

prt120_a4m     <- as.data.frame(cbind(as.vector(prd120_a4m[,1:6]-prd120_a4m[,7]),
                                      rep(1, 1000)))
prt120_a4m[,2] <- labels
colnames(prt120_a4m) <- c("Diff", "Label")

plt120_a4m <- ggplot(prt120_a4m, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 5")

plt120_a4m


prt120_a4l     <- as.data.frame(cbind(as.vector(prd120_a4l[,1:6]-prd120_a4l[,7]),
                                      rep(1, 1000)))
prt120_a4l[,2] <- labels
colnames(prt120_a4l) <- c("Diff", "Label")

plt120_a4l <- ggplot(prt120_a4l, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 4")

plt120_a4l

prt120_t4h     <- as.data.frame(cbind(as.vector(prd120_t4h[,1:6]-prd120_t4h[,7]),
                                      rep(1, 1000)))
prt120_t4h[,2] <- labels
colnames(prt120_t4h) <- c("Diff", "Label")

plt120_t4h <- ggplot(prt120_t4h, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 9")

plt120_t4h

prt120_t4m     <- as.data.frame(cbind(as.vector(prd120_t4m[,1:6]-prd120_t4m[,7]),
                                      rep(1, 1000)))
prt120_t4m[,2] <- labels
colnames(prt120_t4m) <- c("Diff", "Label")

plt120_t4m <- ggplot(prt120_t4m, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 8")

plt120_t4m


prt120_t4l     <- as.data.frame(cbind(as.vector(prd120_t4l[,1:6]-prd120_t4l[,7]),
                                      rep(1, 1000)))
prt120_t4l[,2] <- labels
colnames(prt120_t4l) <- c("Diff", "Label")

plt120_t4l <- ggplot(prt120_t4l, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 7")

plt120_t4l

library(gridExtra)

grid.arrange(plt120_ggl, plt120_ggm, plt120_ggh,
             plt120_a4l, plt120_a4m, plt120_a4h,
             plt120_t4l, plt120_t4m, plt120_t4h, ncol=3)

save.image("K-ICS_summary2 (N=120).RData")

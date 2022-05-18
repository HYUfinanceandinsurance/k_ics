load("K-ICS_analysis2 (J=500, true=g cop, g dist, high cor, N=240).RData")
prd240_ggh  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=500, true=g cop, g dist, mid cor, N=240).RData")
prd240_ggm  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=500, true=g cop, g dist, low cor, N=240).RData")
prd240_ggl  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=500, true=Gumbel cop, 4 dist, high cor, N=240).RData")
prd240_a4h  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=500, true=Gumbel cop, 4 dist, mid cor, N=240).RData")
prd240_a4m  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=500, true=Gumbel cop, 4 dist, low cor, N=240).RData")
prd240_a4l  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=500, true=t cop, 4 dist, high cor, N=240).RData")
prd240_t4h  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=500, true=t cop, 4 dist, mid cor, N=240).RData")
prd240_t4m  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 

load("K-ICS_analysis2 (J=500, true=t cop, 4 dist, low cor, N=240).RData")
prd240_t4l  <- cbind(prd240.vstd, prd240.vemp, prd240.varc, prd240.velp, prd240.cvine, prd240.c2vine, pop.var_true) 


rm(list = ls()[!ls() %in% c("prd240_t4h", "prd240_t4m", "prd240_t4l",
                            "prd240_a4h", "prd240_a4m", "prd240_a4l",
                            "prd240_ggh", "prd240_ggm", "prd240_ggl")])
stable <-
  rbind(c(colMeans(prd240_ggl), colMeans(prd240_ggl)[1:6]-colMeans(prd240_ggl)[7]),
        c(colMeans(prd240_ggm), colMeans(prd240_ggm)[1:6]-colMeans(prd240_ggm)[7]),
        c(colMeans(prd240_ggh), colMeans(prd240_ggh)[1:6]-colMeans(prd240_ggh)[7]),
        c(colMeans(prd240_a4l), colMeans(prd240_a4l)[1:6]-colMeans(prd240_a4l)[7]),
        c(colMeans(prd240_a4m), colMeans(prd240_a4m)[1:6]-colMeans(prd240_a4m)[7]),
        c(colMeans(prd240_a4h), colMeans(prd240_a4h)[1:6]-colMeans(prd240_a4h)[7]),
        c(colMeans(prd240_t4l), colMeans(prd240_t4l)[1:6]-colMeans(prd240_t4l)[7]),
        c(colMeans(prd240_t4m), colMeans(prd240_t4m)[1:6]-colMeans(prd240_t4m)[7]),
        c(colMeans(prd240_t4h), colMeans(prd240_t4h)[1:6]-colMeans(prd240_t4h)[7]))

rownames(stable) <- c("1", "2", "3",
                      "4", "5", "6",
                      "7", "8", "9")
stable
colnames(stable) <- c("Standard", "Empirical", "Archimedian", "Elliptical", "RVine",  "R2Vine","True",
                      "Standard", "Empirical", "Archimedian", "Elliptical", "RVine","R2Vine")

tableA <- stable[,1:7]
tableB <- stable[,8:13]

round(stable)

labels <- c(rep("Standard"  ,500), rep("Empirical" ,500),
            rep("Internal 1",500), rep("Internal 2",500),
            rep("Internal 3",500), rep("Internal 4",500))
labels <- factor(labels, levels = c("Standard", "Empirical", "Internal 1", "Internal 2", "Internal 3", "Internal 4"))

library(ggplot2)

prt240_ggh     <- as.data.frame(cbind(as.vector(prd240_ggh[,1:6]-prd240_ggh[,7]),
                                      rep(1, 500)))
prt240_ggh[,2] <- labels
colnames(prt240_ggh) <- c("Diff", "Label")

plt240_ggh <- ggplot(prt240_ggh, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 3")

plt240_ggh


prt240_ggm     <- as.data.frame(cbind(as.vector(prd240_ggm[,1:6]-prd240_ggm[,7]),
                                      rep(1, 500)))
prt240_ggm[,2] <- labels
colnames(prt240_ggm) <- c("Diff", "Label")

plt240_ggm <- ggplot(prt240_ggm, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 2")

plt240_ggm


prt240_ggl     <- as.data.frame(cbind(as.vector(prd240_ggl[,1:6]-prd240_ggl[,7]),
                                      rep(1, 500)))
prt240_ggl[,2] <- labels
colnames(prt240_ggl) <- c("Diff", "Label")

plt240_ggl <- ggplot(prt240_ggl, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 1")

plt240_ggl


prt240_a4h     <- as.data.frame(cbind(as.vector(prd240_a4h[,1:6]-prd240_a4h[,7]),
                                      rep(1, 500)))
prt240_a4h[,2] <- labels
colnames(prt240_a4h) <- c("Diff", "Label")

plt240_a4h <- ggplot(prt240_a4h, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 6")

plt240_a4h

prt240_a4m     <- as.data.frame(cbind(as.vector(prd240_a4m[,1:6]-prd240_a4m[,7]),
                                      rep(1, 500)))
prt240_a4m[,2] <- labels
colnames(prt240_a4m) <- c("Diff", "Label")

plt240_a4m <- ggplot(prt240_a4m, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 5")

plt240_a4m


prt240_a4l     <- as.data.frame(cbind(as.vector(prd240_a4l[,1:6]-prd240_a4l[,7]),
                                      rep(1, 500)))
prt240_a4l[,2] <- labels
colnames(prt240_a4l) <- c("Diff", "Label")

plt240_a4l <- ggplot(prt240_a4l, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 4")

plt240_a4l

prt240_t4h     <- as.data.frame(cbind(as.vector(prd240_t4h[,1:6]-prd240_t4h[,7]),
                                      rep(1, 500)))
prt240_t4h[,2] <- labels
colnames(prt240_t4h) <- c("Diff", "Label")

plt240_t4h <- ggplot(prt240_t4h, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 9")

plt240_t4h

prt240_t4m     <- as.data.frame(cbind(as.vector(prd240_t4m[,1:6]-prd240_t4m[,7]),
                                      rep(1, 500)))
prt240_t4m[,2] <- labels
colnames(prt240_t4m) <- c("Diff", "Label")

plt240_t4m <- ggplot(prt240_t4m, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 8")

plt240_t4m


prt240_t4l     <- as.data.frame(cbind(as.vector(prd240_t4l[,1:6]-prd240_t4l[,7]),
                                      rep(1, 500)))
prt240_t4l[,2] <- labels
colnames(prt240_t4l) <- c("Diff", "Label")

plt240_t4l <- ggplot(prt240_t4l, aes(x=Label, y=Diff, fill=Label)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Scenario 7")

plt240_t4l

library(gridExtra)

grid.arrange(plt240_ggl, plt240_ggm, plt240_ggh,
             plt240_a4l, plt240_a4m, plt240_a4h,
             plt240_t4l, plt240_t4m, plt240_t4h, ncol=3)

save.image("K-ICS_summary2 (N=240).RData")

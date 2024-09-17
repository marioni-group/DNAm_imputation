library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggpubr)


###########################
## Figure 1
###########################

d = read.csv("elovl2_cg16867657_Hannum_measured_and_predicted.csv")

d2 <- read.csv("elovl2_cg16867657_GSE246337_measured_and_predicted.csv")

d2 <- d2[!is.na(d2$elovl2_episcore),]
d2 <- d2[!is.na(d2$age),]

smk = read.csv("ahrr_cg05575921_LBC_measured_and_predicted.csv")


par(mfrow=c(3,3))
plot(scale(smk$ahrr_DNAm), scale(smk$ahrr_episcore))
boxplot(scale(smk$ahrr_DNAm) ~ smk$smoke)
boxplot(scale(smk$ahrr_episcore) ~ smk$smoke)

plot(scale(d2$elovl2_DNAm), scale(d2$elovl2_episcore))
plot(scale(d2$elovl2_DNAm) ~ d2$age)
plot(scale(d2$elovl2_episcore) ~ d2$age)

plot(scale(d3$elovl2_DNAm), scale(d3$elovl2_episcore))
plot(scale(d3$elovl2_DNAm) ~ d3$age)
plot(scale(d3$elovl2_episcore) ~ d3$age)



smk1 <- ggplot(smk, aes(x = scale(ahrr_DNAm), y = scale(ahrr_episcore))) +
  geom_point(aes(color=factor(smoke))) +
  scale_color_discrete(labels=c("never", "former", "current"), name="smoking status") + 
  stat_smooth(method="loess") +
  xlab("CpG") +
  ylab("EpiScore") +
annotate("text", x = max(scale(smk$ahrr_DNAm)*0.95), y = min(scale(smk$ahrr_episcore)*1.05), size = 6, label = paste("r = 0.87, P = 4.8x10-275"), vjust = 0, hjust = 1) +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14), legend.position="bottom")

smk2 <- ggplot(smk, aes(y = scale(ahrr_DNAm), x = factor(smoke))) +
  geom_violin() +
  geom_boxplot(width=0.2) +
  scale_x_discrete(name="Self-reported smoking", 
                   labels=c("0" = "Never (n=410)", '1'= 'Former (n=376)', '2'= 'Current (n=109)')) +
  ylab("CpG") +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))

smk3 <- ggplot(smk, aes(y = scale(ahrr_episcore), x = factor(smoke))) +
  geom_violin() +
  geom_boxplot(width=0.2) +
  scale_x_discrete(name="Self-reported smoking", 
                   labels=c("0" = "Never (n=410)", '1'= 'Former (n=376)', '2'= 'Current (n=109)')) +
                   ylab("EpiScore") +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))

p1 <- ggarrange(smk1, smk2, smk3, ncol=3) 

p1 <- annotate_figure(p1, top = text_grob("LBC1936, cg05575921 and smoking", face = "bold", size = 18))


elo1 <- ggplot(d, aes(x = scale(elovl2_DNAm), 
                       y = scale(elovl2_episcore))) + 
  geom_point() +
  stat_smooth(method="loess") +
  xlab("CpG") +
  ylab("EpiScore") +
annotate("text", x = max(scale(d$elovl2_DNAm)*0.95), y = min(scale(d$elovl2_episcore)*1.05), size = 6, label = paste("r = 0.90, P = 8.1x10-238"), vjust = 0, hjust = 1) +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))


elo2 <- ggplot(d, aes(x = age, 
                       y = scale(elovl2_DNAm))) + 
  geom_point() +
  stat_smooth(method="loess") +
  xlab("Age (years)") +
  ylab("CpG") + 
annotate("text", x = max(d$age*0.95), y = min(scale(d$elovl2_DNAm)*1.05), size = 6, label = paste("r = 0.83, P = 7.4x10-167"), vjust = 0, hjust = 1) +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))


elo3 <- ggplot(d, aes(x = age, 
                       y = scale(elovl2_episcore))) + 
  geom_point() +
  stat_smooth(method="loess") +
  xlab("Age (years)") +
  ylab("EpiScore") +
annotate("text", x = max(d$age*0.95), y = min(scale(d$elovl2_episcore)*1.05), size = 6, label = paste("r = 0.88, P = 8.4x10-214"), vjust = 0, hjust = 1) +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))

p2 <- ggarrange(elo1, elo2, elo3, ncol=3) 

p2 <- annotate_figure(p2, top = text_grob("Hannum/GSE40279, cg16867657 and Age", face = "bold", size = 18))


elo4 <- ggplot(d2, aes(x = scale(elovl2_DNAm), 
                       y = scale(elovl2_episcore))) + 
  geom_point() +
  stat_smooth(method="loess") +
  xlab("CpG") +
  ylab("EpiScore") +
annotate("text", x = max(scale(d2$elovl2_DNAm)*0.95), y = min(scale(d2$elovl2_episcore)*1.05), size = 6, label = paste("r = 0.96, P = 2.0x10-239"), vjust = 0, hjust = 1) +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))


elo5 <- ggplot(d2, aes(x = age, 
                       y = scale(elovl2_DNAm))) + 
  geom_point() +
  stat_smooth(method="loess") +
  xlab("Age (years)") +
  ylab("CpG") +
annotate("text", x = max(d2$age*0.95), y = min(scale(d2$elovl2_DNAm)*1.05), size = 6, label = paste("r = 0.91, P = 1.5x10-185"), vjust = 0, hjust = 1) +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))


elo6 <- ggplot(d2, aes(x = age, 
                       y = scale(elovl2_episcore))) + 
  geom_point() +
  stat_smooth(method="loess") +
  xlab("Age (years)") +
  ylab("EpiScore") +
annotate("text", x = max(d2$age*0.95), y = min(scale(d2$elovl2_episcore)*1.05), size = 6, label = paste("r = 0.94, P = 3.7x10-204"), vjust = 0, hjust = 1) +
theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))


p3 <- ggarrange(elo4, elo5, elo6, ncol=3) 

p3 <- annotate_figure(p3, top = text_grob("GSE246337, cg16867657 and Age", face = "bold", size = 18))


p <- ggarrange(p1, p2, p3, nrow=3) 

jpeg("DNAm_imputation_figure1.jpeg", width=1000, height=900)
p
dev.off()




###########################
## Supplementary Figure 1
###########################

library(RMariaDB)

smoking <- dbReadTable(gwasp, "smoking")
agesex <- dbReadTable(gwasp, "agesex")

meth_ids <- read.table("methylation_ids.txt", header=TRUE)

smoking <- smoking[which(smoking$id %in% meth_ids$x),]

agesex <- agesex[,1:3]

smoking <- merge(smoking, agesex, by="id", all.x=TRUE)
smoking$smoke[smoking$status=="smoker"] <- "current"
smoking$smoke[smoking$status=="ex-smoker"] <- "former"
smoking$smoke[smoking$status=="non-smoker"] <- "never"

s1 <- ggplot(smoking, aes(x=age.y, color=sex, fill=sex)) +
geom_density(alpha=0.4) +
xlab("Age (years)") +
theme(axis.text=element_text(size=12))

smk_status <- data.frame(smoke=c("never", "former", "current"), value=c(9632, 5400, 3228))
smk_status$percent <- (smk_status$value / 18260) * 100
smk_status$percent <- round(smk_status$percent, 1)
smk_status$label <- paste0(smk_status$value, " (", smk_status$percent, "%)")

s2 <- ggplot(smk_status, aes(x="", y=value, fill=factor(smoke))) + 
geom_bar(width=1, stat="identity") + 
coord_polar("y", start=0) + 
guides(fill=guide_legend(title="smoking status")) + 
theme_minimal() + 
theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
geom_text(aes(label=label), position=position_stack(vjust = 0.5))

smk_packy <- smoking[which(smoking$smoke=="former" | smoking$smoke=="current"),]

sd(smk_packy$pack_years, na.rm=TRUE) * 6

smk_packy$pack_years[smk_packy$pack_years>139.2] <- NA

s3 <- ggplot(smk_packy, aes(x=pack_years, y=age.y, color=factor(smoke))) + 
geom_point() + 
geom_density_2d(color='black') + 
scale_color_discrete(labels=c("current", "former"), name="smoking status") + 
xlab("Pack Years") + 
ylab("Age (years)") +
theme(axis.text=element_text(size=12)) +
facet_wrap(~ smoke)

s <- ggarrange(s1, s2, s3, nrow=3)

jpeg("DNAm_imputation_suppfigure1.jpeg", width=500, height=1000)
s
dev.off()


rm(list=ls())

getwd() # LINE 250
# setwd("/Users/valentinamelica/R/Grey-Whale")

suppressMessages(library(cowplot))
suppressMessages(library(plotrix))
suppressMessages(library(patchwork))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(lattice))
suppressMessages(library(ggmap))
suppressMessages(library(ggbeeswarm))

# data manipulation
suppressMessages(library(openxlsx))
suppressMessages(library(plyr))
suppressMessages(library(here))
suppressMessages(library(dplyr)) # don't need dplyr if you have tidyverse 
suppressMessages(library(tidyverse))

library(MASS)
library(mixtools)
library(MuMIn)
library(gapminder) # package with a gapminder dataset
library(emmeans)
library(gplots) # to plot means
library(rstatix)
library(stringr)
library(patchwork)
library(Hmisc)
library(ggforce)
library(lme4)
library(here)
library(ggsignif) 

install.packages("userfriendlyscience")
citation("stats")
options(na.action = "na.omit")
?qqnorm
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

ErB<- read.csv(here("./data/GW_Biol.csv"))
str(ErB)
colnames(ErB)
colnames(ErB) <- c("vial", "MMSNW" , "ID","BlubberDEPTH","Age.class","UME","Yr","Mo","Da","Month",
                  "Latitude.N","Longitude.W","State", "CarCode" ,  "BodyCond","Adjusted_BC",  "Sex", "COD","Lipid",      
                  "Cts", "Cts CV", "Ctn" ,  "Ctn CV",  "Aldo", "Aldo CV")
head(ErB)
# view(ErB)
ErB$MMSNW = as.factor(ErB$MMSNW)
ErB$ID = as.factor(ErB$ID)
ErB$BlubberDEPTH = as.factor(ErB$BlubberDEPTH)
levels(ErB$BlubberDEPTH)
summary(ErB$BlubberDEPTH)
#35 samples all from surface
ErB$Age.class= as.factor(ErB$Age.class)
summary(ErB$Age.class)
#19 adults and 14 immature. 2 unknown
ErB$Sex = as.factor(ErB$Sex)
summary(ErB$Sex)
#22 female and 13 male
ErB$State = as.factor(ErB$State)
summary(ErB$State)
ErB$Adjusted_BC =as.factor(ErB$Adjusted_BC)
levels(ErB$Adjusted_BC)
summary(ErB$Adjusted_BC)
ErB$Adjusted_BC<- ordered(ErB$Adjusted_BC, levels= c("Good","Fair", "Poor"))
ErB$CarCode = as.factor(ErB$CarCode)
summary(ErB$CarCode)
ErB$COD= as.factor(ErB$COD)
levels(ErB$COD)
#ErB$COD<- ordered(ErB$COD, levels= c("Trauma", "Nutritional stress", "Entanglement", "Predation", "Undertermined"))
summary(ErB$COD)

ErB$Yr = as.factor(ErB$Yr)
summary(ErB$Yr)

##### correlations
cor.test(ErB$Cts, ErB$Ctn)
# t = 1.0388, df = 31, p-value = 0.3069
# r= 0.18
cor.test(ErB$Aldo, ErB$Ctn)
# t = 1.8926, df = 31, p-value = 0.06778
# r= 0.32
cor.test(ErB$Aldo, ErB$Cts)
# t = -0.23633, df = 30, p-value = 0.8148
# r= -0.04
cor.test(ErB$Aldo, ErB$Lipid)
# t = 0.25365, df = 23, p-value = 0.802
# r= 0.05
cor.test(ErB$Cts, ErB$Lipid)
# t = 0.84068, df = 24, p-value = 0.4088
# r= 0.17
cor.test(ErB$Ctn, ErB$Lipid)
# t = 1.9794, df = 23, p-value = 0.05986
# r= 0.38

#create a table of ourliers
ErB_outliers <-ErB %>%
  filter(Cts > 30 | Ctn >20 | Aldo >10) %>%
  unite(date, Yr, Mo, Da, sep= "-")%>%
  select(-c(vial, `Cts CV`, `Ctn CV`, `Aldo CV`))

write.csv(ErB_outliers, here("./output/Table_outliers.csv"), row.names = FALSE)

#create a table of ourliers
# ErB_table <-ErB %>%
#   unite(date, Yr, Mo, Da, sep= "-")%>%
#   select(-c(vial,BlubberDEPTH, `Cts CV`, `Ctn CV`, `Aldo CV`, BodyCond, Month)) 
#  
# 
# write.csv(ErB_table, here("./output/Table_1.csv"), row.names = FALSE)


ErB %>%
  group_by(COD) %>%  
  dplyr::summarize(n = n())
ErB %>%
  group_by(COD) %>%  
  dplyr::summarize(n = n(),
                   meanCTS = mean(Cts, na.rm=T),
                   sdCTS = sd(Cts, na.rm=T),
                   rangeCTS = max(Cts, na.rm=T) - min(Cts,na.rm=T),
                   meanCTN = mean(Ctn, na.rm=T),
                   sdCTN = sd(Ctn, na.rm=T),
                   rangeCTN = max(Ctn, na.rm=T) - min(Ctn, na.rm=T),
                   meanALDO = mean(Aldo, na.rm=T),
                   sdALDO = sd(Aldo, na.rm=T),
                   rangeALDO = max(Aldo, na.rm=T) - min(Aldo, na.rm=T))




#### Section 1. COD
##### compare pre ume and ume samples 
# test if UME is a factor
ErBCTS <- ErB[!is.na(ErB$Cts),]
tapply(ErBCTS$Cts, FUN = mean , INDEX = ErBCTS$COD)
tapply(ErBCTS$Cts, FUN = sd , INDEX = ErBCTS$COD)
summary(ErBCTS)

performance::check_outliers(ErBCTS$Cts)

# test if UME is a factor
hist(ErBCTS$Cts)
shapiro.test(ErBCTS$Cts) # p < 0.001
fit <- lm(Cts ~ UME, data= ErBCTS)
# plot(fit)
# shapiro.test(resid(fit)) # p < 0.001
# #apply log tranformation
# shapiro.test(log(ErBCTS$Cts)) # p = 0.6
# fit <- lm(log(Cts) ~ UME, data= ErBCTS)
# #plot(fit) # no heteroschedasticity
# shapiro.test(resid(fit)) # p= 0.7
# t.test(log(Cts) ~ UME, data= ErBCTS)
#t = 1.6646, df = 31.991, p-value = 0.1058
colnames(ErBCTS)
?wilcox_test
wilcox_test(Cts ~ UME,  p.adjust.method = "bonferroni", data= ErBCTS)
#same results with this p= 0.12
Cts_UME <- ggplot(data = ErBCTS, aes(x = UME, y=Cts))+
  geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))

#conclusion: for cortisol UME samples can be divided in the remaining COD groups

ErBCTN <- ErB[!is.na(ErB$Ctn),]

performance::check_outliers(ErBCTN$Ctn)

# test if UME is a factor
hist(ErBCTN$Ctn)
shapiro.test(ErBCTN$Ctn) # p < 0.001
# fit <- lm(Ctn ~ UME, data= ErBCTN)
# plot(fit)
# shapiro.test(resid(fit)) # p < 0.001
# #apply log tranformation
# shapiro.test(log(ErBCTN$Ctn)) # p = 0.2
# fit <- lm(log(Ctn) ~ UME, data= ErBCTN)
# plot(fit) # variance is good
# shapiro.test(resid(fit)) # p= 0.4
# t.test(log(Ctn) ~ UME, data= ErBCTN)
#t = 0.97649, df = 28.862, p-value = 0.3369
wilcox_test(Ctn ~ UME,  p.adjust.method = "bonferroni", data= ErBCTN)
#same results with this p= 0.34
Ctn_UME <-ggplot(data = ErBCTN, aes(x = UME, y=Ctn))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))
# conclusion : no significant difference in corticosterone concentrations between samples 
# collected before or during the UME


ErBALDO <- ErB[!is.na(ErB$Aldo),]
hist(ErBALDO$Aldo)
performance::check_outliers(ErBALDO$Aldo)
ErBALDO[c(24, 32),] # there is one big outlier , a female with concentration of aldosterone around 15 ng.g
shapiro.test(ErBALDO$Aldo) # p < 0.001
fit <- lm(Aldo ~ UME, data= ErBALDO)
plot(fit)
shapiro.test(resid(fit)) # p < 0.001
#apply log tranformation
shapiro.test(log(ErBALDO$Aldo)) # p = 0.02
fit <- lm(log(Aldo) ~ UME, data= ErBALDO)
plot(fit) # variance is good
shapiro.test(resid(fit)) # p= 0.005
qqnorm(resid(fit)) # its honestly not too bad

t.test(log(Aldo) ~ UME, data= ErBALDO)
# t = 2.2232, df = 27.509, p-value = 0.0346

wilcox_test(Aldo ~ UME,  p.adjust.method = "bonferroni", data= ErBALDO)
#same results with this p= 0.01
# there is a significant difference however it looks like it is driven by that outlier
aldo_f <-subset(ErBALDO, ErBALDO$Aldo < 10)
aldo_f <-subset(ErBALDO, ErBALDO$Aldo < 1)
shapiro.test(aldo_f$Aldo) #p= 0.00001
wilcox_test(Aldo ~ UME,  p.adjust.method = "bonferroni", data= aldo_f)
#still significant 
aldo_f[c(26, 33),] 

# fit <- lm(Aldo ~ UME, data= aldo_f)
# 
# plot(fit)
# shapiro.test(resid(fit)) # p < 0.001
#apply log tranformation
# shapiro.test(log(aldo_f$Aldo)) # p = 0.355
# fit <- lm(log(Aldo) ~ UME, data= aldo_f)
# plot(fit) # variance is good
# shapiro.test(resid(fit)) # p= 0.12
# qqnorm(resid(fit)) # its honestly not too bad
# t.test(log(Aldo) ~ UME, data= aldo_f)
#t = 1.9634, df = 29.915, p-value = 0.05896

Aldo_UME <-ggplot(data = aldo_f, aes(x = UME, y=Aldo))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))
# conclusion : if we remove the outlier there is no significant difference in aldosterone concentrations between samples 
# collected before or during the UME

####aldosterone validation 
ErBALDO %>%
  group_by(Sex, Age.class) %>%  
  dplyr::summarize(n = n(),
                   meanALDO = mean(Aldo, na.rm=T),
                   sdALDO = sd(Aldo, na.rm=T),
                   rangeALDO = max(Aldo, na.rm=T) - min(Aldo, na.rm=T))

ggplot(data = subset(ErBALDO, ErBALDO$Age.class!= "Unknown") , aes(x = Sex, y=Aldo, fill= Age.class))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=8, face= "bold"))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))

UME<- ggarrange(Cts_UME,Ctn_UME,Aldo_UME, nrow = 1,ncol = 3)
UME
# save figure 
ggsave(here("./output/UME.png"), UME, width = 10, height = 6, dpi = 300)


#####END HERE
######AGE CLASS and SEX
ErBAC <-ErB %>%
  filter(Age.class!= "Unknown")
str(ErBAC)
ErBAC %>%
  #filter(Cts < 30)  %>%
  group_by(Age.class, Sex) %>%
  get_summary_stats(Cts, type = "mean_sd")

ErBAC %>%
  # filter(Aldo < 10)  %>%
  group_by(Age.class, Sex) %>%
  get_summary_stats(Aldo, type = "mean_sd")

ErBAC %>%
  # filter(Ctn < 20)  %>%
  group_by(Age.class, Sex) %>%
  get_summary_stats(Ctn, type = "mean_sd")
ErBAC %>%
  filter(Ctn < 20)  %>%
  group_by(Age.class) %>%
  get_summary_stats(Ctn, type = "mean_sd")
# higher in immature
ErBAC %>%
  filter(Ctn < 20)  %>%
  group_by(Sex) %>%
  get_summary_stats(Ctn, type = "mean_sd")
# higher in females with and w/o the outlier

ErBAC %>%
  group_by(Age.class, Sex) %>%
  get_summary_stats(Lipid, type = "mean_sd")

# Cortisol
shapiro.test(resid(fit1)) # p=0.1847
qqnorm(resid(fit1))
# plot(fit1,1)
bartlett.test(log(ErBAC$Cts), ErBAC$Age.class)
# ok
# kruskal.test(Cts ~ Age.class, data= ErBAC)
#Kruskal-Wallis chi-squared = 3.1236, df = 1, p-value = 0.07717
performance::check_outliers(ErBAC$Cts)
ErBAC[c(8,18),] #samples with cortisol > 30 ng.g
model_GC_null <- lm(log(Cts) ~ -1, data= ErBAC )
fit1= lm(log(Cts) ~ Age.class -1, data= ErBAC)
fit2 = lm(log(Cts) ~ Age.class*Sex -1, data= ErBAC)
fit3 = lm(log(Cts) ~Sex -1, data= ErBAC)
fit4 = lm(log(Cts) ~Age.class+Sex -1, data= ErBAC)
# fit5 = lm(log(Cts) ~Adjusted_BC -1, data= ErBAC)

MuMIn::AICc(model_GC_null,fit1, fit2,fit3, fit4)
# fit 1 the model with Age class 
anova(model_GC_null,fit1, fit2,fit3, fit4)
anova(model_GC_null,fit1, fit4)
#the fit1 model is significantly different from the null model

# the posthoc test is not really needed becayse we only have two levels in both age classes and sex
# posthoc_CTS <- emmeans(fit1, ~ Age.class)
# summary(posthoc_CTS)
# posthoc_CTS %>% 
#   pairs() 
# # sifgnificant difference 

#Immature females have significantly different CTS from Adult female and Adult males
# 
# fit1= glm(Cts ~ Age.class -1,data= subset(ErBAC,ErBAC$Cts < 30))
# fit2 = glm(Cts ~ Age.class*Sex -1, data= subset(ErBAC,ErBAC$Cts < 30))
# fit3 = glm(Cts ~Sex -1, data=subset(ErBAC,ErBAC$Cts < 30))
# 
# model_GC_null <- glm(Cts ~ -1, data= subset(ErBAC,ErBAC$Cts < 30))
# MuMIn::AICc(model_GC_null,fit1, fit2,fit3) # fit1
# 
# wilcox.test(Cts ~ Age.class, data=subset(ErBAC,ErBAC$Cts < 30))
# posthoc_CTS <- emmeans(fit2, ~ Age.class*Sex)
# summary(posthoc_CTS)
# posthoc_CTS %>% 
#   pairs() 

#corticosterone
shapiro.test(log(ErBAC$Ctn)) # p = 0.17
fit1 = lm(log(Ctn) ~ Age.class, data= ErBAC)
summary(fit1)
shapiro.test(resid(fit1)) # p = 0.35
qqnorm(resid(fit1))
bartlett.test(log(ErBAC$Ctn), ErBAC$Age.class) # p = 0.93 
bartlett.test(log(ErBAC$Ctn), ErBAC$Sex) # p= 0.54
# meet the requirements for parametric test

performance::check_outliers(ErBAC$Ctn)

ErBAC[c(3, 9,25,33),]
fit1= lm(log(Ctn) ~ Age.class -1, data= ErBAC)
fit2 = lm(log(Ctn) ~ Age.class*Sex -1, data= ErBAC)
fit3 = lm(log(Ctn) ~Sex -1, data= ErBAC)
fit4 = lm(log(Ctn) ~ Age.class+Sex -1, data= ErBAC)
model_GC_null <- lm(log(Ctn) ~ -1, data= ErBAC)
MuMIn::AICc(model_GC_null,fit1, fit2,fit3, fit4)
#               df     AICc
# model_GC_null  1 110.6472
# fit1           3 109.0164 
# fit2           5 109.4163
# fit3           3 111.2480
# fit4           4 106.7711 <- this is the lowest
anova(model_GC_null,fit1,fit4)
#the addition of sex 
# t.test(log(Ctn) ~ Age.class, data= ErBAC)
# 
# posthoc_CTN <- emmeans(fit1, ~ Age.class)
# summary(posthoc_CTN)
# posthoc_CTN %>% 
#   pairs() 
# posthoc_CTN <- emmeans(fit4, ~ Sex)
# summary(posthoc_CTN)
# posthoc_CTN %>% 
#   pairs() 

library(car)
fit <- lm(log(Ctn) ~ Age.class + Sex, data = ErBAC)
Anova(fit, type = 2)

# 
# fit5 = lm(log(Ctn) ~ Age.class:Sex -1, data= ErBAC)
# summary(fit5)
# posthoc_CTN <- emmeans(fit2, ~ Age.class*Sex)
# summary(posthoc_CTN)
# posthoc_CTN %>% 
#   pairs() 

# fit1= glm(Ctn ~ Age.class - 1, data= subset(ErBAC,ErBAC$Ctn < 20))
# fit2 = glm(Ctn ~ Age.class*Sex - 1, data= subset(ErBAC,ErBAC$Ctn < 20))
# fit3 = glm(Ctn ~Sex - 1, data=subset(ErBAC,ErBAC$Ctn < 20))
# model_GC_null <- glm(Ctn ~ - 1, data= subset(ErBAC,ErBAC$Ctn < 20))
# MuMIn::AICc(model_GC_null,fit1, fit2,fit3)
# ErBAC %>%
#   filter(Ctn < 20)  %>%
#   group_by( Sex) %>%
#   get_summary_stats(Ctn, type = "mean_sd")
# 
# wilcox.test(Ctn ~ Sex, data= subset(ErBAC,ErBAC$Ctn < 20))
# posthoc_CTN <- emmeans(fit3, ~ Sex)
# summary(posthoc_CTN)
# posthoc_CTN %>% 
#   pairs() 

#Aldosterone with outlier
performance::check_outliers(ErBAC$Aldo)
ErBAC[c(13,22,11),]
shapiro.test(log(ErBAC$Aldo))
qqnorm(log(ErBAC$Aldo))

fit1 = lm(log(Aldo) ~ Sex, data= ErBAC)
summary(fit1)
shapiro.test(resid(fit1))
qqnorm(resid(fit1))
bartlett.test(log(ErBAC$Aldo), ErBAC$Sex)

#ok
fit1= lm(log(Aldo) ~ Age.class-1, data= ErBAC)
fit2 = lm(log(Aldo) ~ Age.class*Sex-1, data= ErBAC)
fit3 = lm(log(Aldo) ~Sex -1, data= ErBAC)
fit4 = lm(log(Aldo) ~ Age.class+Sex-1, data= ErBAC)
model_GC_null <- lm(log(Aldo) ~ -1, data= ErBAC)
MuMIn::AICc(model_GC_null,fit1, fit2,fit3, fit4) # fit 3
anova(model_GC_null,fit3, fit4)
# 
# # wilcox.test(Aldo ~ Age.class, data= ErBAC)
# 
# posthoc_Aldo <- emmeans(fit1, ~ Age.class)
# summary(posthoc_Aldo)
# posthoc_Aldo %>% 
#   pairs() 
# posthoc_Aldo <- emmeans(fit3, ~ Sex)
# summary(posthoc_Aldo)
# posthoc_Aldo %>% 
#   pairs() 

# #without outlier
# fit1= glm(Aldo ~ Age.class -1, data= subset(ErBAC,ErBAC$Aldo <10))
# fit2 = glm(Aldo ~ Age.class*Sex -1, data= subset(ErBAC,ErBAC$Aldo <10))
# fit3 = glm(Aldo ~Sex -1, data= subset(ErBAC,ErBAC$Aldo <10))
# 
# model_GC_null <- glm(Aldo ~ -1, data= subset(ErBAC,ErBAC$Aldo <10))
# MuMIn::AICc(model_GC_null,fit1, fit2,fit3)
# posthoc_Aldo <- emmeans(fit2, ~ Age.class*Sex)
# summary(posthoc_Aldo)
# posthoc_Aldo %>% 
#   pairs() 
# # significan for females

shapiro.test(log(ErBAC$Lipid))
#log transformed 

fit1 = lm(log(Lipid) ~ Age.class, data= ErBAC)
summary(fit1)
shapiro.test(resid(fit1))
qqnorm(resid(fit1))
bartlett.test(log(ErBAC$Lipid), ErBAC$Age.class)
#ok
fit1= lm(log(Lipid) ~ Age.class, data= ErBAC)
fit2 = lm(log(Lipid) ~ Age.class*Sex, data= ErBAC)
fit3 = lm(log(Lipid) ~Sex, data= ErBAC)
fit4 = lm(log(Lipid) ~ Age.class+Sex, data= ErBAC)
model_GC_null <- lm(log(Lipid) ~ -1, data= ErBAC )
MuMIn::AICc(model_GC_null,fit1, fit2,fit3, fit4)
anova(model_GC_null,fit1)
# 
# posthoc_Lipid <- emmeans(fit1, ~ Age.class)
# summary(posthoc_Lipid)
# posthoc_Lipid %>% 
#   pairs() 
# 

# showSignificance( c(1,2), 450, -0.05, "ns") + 
#   showSignificance( c(1,3), 500, -0.05, "ns") +
#   showSignificance( c(2,3), 425, -0.05, "ns") +

#####blar plot with sd
# figure
df1= ErBAC %>%
  #filter(Cts < 30)  %>%
  group_by(Age.class, Sex) %>%
  get_summary_stats(Cts, type = "mean_sd")

Fig2a <- ggplot(df1, aes(x=Age.class, y=mean, fill=Sex)) + 
  geom_bar(stat="identity",color= "black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ylab("Cortisol ng/g")+ xlab("")+
  theme_minimal()
  # geom_signif(comparisons = list(c("Adult", "Immature")),
  #             map_signif_level = TRUE,
  #             annotations = c("*"))
Fig2a

df2= ErBAC %>%
  #filter(Cts < 30)  %>%
  group_by(Age.class, Sex) %>%
  get_summary_stats(Ctn, type = "mean_sd")
Fig2b <- ggplot(df2, aes(x=Age.class, y=mean, fill=Sex)) + 
  geom_bar(stat="identity",color= "black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ylab("Corticosterone ng/g")+ xlab("")+
  theme_minimal()
  # geom_signif(comparisons = list(c("Adult", "Immature")),
  #             map_signif_level = TRUE,
  #             annotations = c("*"))
Fig2b 
df3= ErBAC %>%
  #filter(Cts < 30)  %>%
  group_by(Age.class, Sex) %>%
  get_summary_stats(Aldo, type = "mean_sd")
Fig2c <- ggplot(df3, aes(x=Age.class, y=mean, fill=Sex)) + 
  geom_bar(stat="identity",color= "black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ylab("Aldosterone ng/g")+ xlab("")+
  theme_minimal()
Fig2c
df4= ErBAC %>%
  #filter(Cts < 30)  %>%
  group_by(Age.class, Sex) %>%
  get_summary_stats(Lipid, type = "mean_sd")
Fig2d <- ggplot(df4, aes(x=Age.class, y=mean, fill=Sex)) + 
  geom_bar(stat="identity",color= "black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ylab("Lipid content %")+ xlab("")+
  theme_minimal()
Fig2d

Ageclass <- ggarrange(Fig2a, Fig2b, Fig2c,Fig2d, nrow=2, ncol = 2)

here()
ggsave(here("./output/Ageclass.png"),Ageclass, width=8, height=8, dpi=600)
###########
# Looking here at a different approach since we have everything log transformed
fig2a <- ggplot(data = ErBAC, aes(x = Age.class, y=Cts, fill= Sex))+
  geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g, log scale")+
  # geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
  #              width = .75, linetype = "dashed")+
  scale_y_continuous(
    trans = "log",
    breaks = scales::breaks_log(n = 4),
    labels = scales::label_number(accuracy = 0.01)
  ) + annotation_logticks(sides = "l")+
  theme_minimal()+ #scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text= element_text(size=10, face= "bold"))+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))+
  labs(tag = "A")


fig2b <- ggplot(data = ErBAC, aes(x = Age.class, y=Ctn, fill= Sex))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g, log scale")+
  # geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
  #              width = .75, linetype = "dashed")+
  scale_y_continuous(
    trans = "log",
    breaks = scales::breaks_log(n = 4),
    labels = scales::label_number(accuracy = 0.01)
  ) + annotation_logticks(sides = "l")+
  theme_minimal()+ #scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text= element_text(size=10, face= "bold"))+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))+
  labs(tag = "B")

fig2c <- ggplot(data = ErBAC, aes(x = Age.class, y=Aldo, fill= Sex))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g, log scale")+
  # geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
  #              width = .75, linetype = "dashed")+
  scale_y_continuous(
    trans = "log",
    breaks = scales::breaks_log(n = 4),
    labels = scales::label_number(accuracy = 0.01)
  ) + annotation_logticks(sides = "l")+
  theme_minimal()+ #scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text= element_text(size=10, face= "bold"))+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))+
  labs(tag = "C")

fig2d<- ggplot(data = ErBAC, aes(x = Age.class, y=Lipid, fill= Sex))+
  geom_boxplot()+ xlab("")+ ylab("Lipid content, log scale")+
  # geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
  #              width = .75, linetype = "dashed")+
  scale_y_continuous(
    trans = "log",
    breaks = scales::breaks_log(n = 4),
    labels = scales::label_number(accuracy = 0.01)
  ) + annotation_logticks(sides = "l")+
  theme_minimal()+ #scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text= element_text(size=10, face= "bold"))+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))+
  labs(tag = "D")

Ageclass_log <- ggarrange(fig2a, fig2b, fig2c,fig2d, nrow=2, ncol = 2)

here()
ggsave(here("./output/Fig2_Log_Age&Sex.png"),Ageclass_log, width=8, height=8, dpi=600)





Fig5b <- ggplot(data = df3, aes(x = Age.class, y=Cts, fill=Sex))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  scale_fill_manual(values=c("grey43", "ghostwhite"))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))
Fig5b

Fig5c <- ggplot(data = ErBAC, aes(x = Age.class, y=Aldo, fill=Sex))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  scale_fill_manual(values=c("grey43", "ghostwhite"))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))
Fig5c 
Fig5d <- ggplot(data = ErBAC, aes(x = Age.class, y=Lipid, fill=Sex))+
  geom_boxplot()+ xlab("")+ ylab("Lipid content %")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  scale_fill_manual(values=c("grey43", "ghostwhite"))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))
Fig5_with <- ggarrange(Fig5a, Fig5b, Fig5c, Fig5d, nrow=2, ncol = 2)

ggsave(here("./output/Fig_with_outlier/Age_class.png"), Fig5_with, width=6, height=6, dpi=300 )

# # WITHOUT
# 
# Fig5a <- ggplot(data = subset(ErBAC,ErBAC$Cts < 30), aes(x = Age.class, y=Cts, fill=Sex))+
#   geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g")+
#   theme_bw()+scale_shape_manual(values=c(8, 1))+
#   scale_fill_manual(values=c("grey43", "ghostwhite"))+
#   theme(panel.border= element_blank(), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top",
#         axis.line = element_line(linewidth = 0.7, linetype = "solid",
#                                  colour = "black"))+
#   theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
#   
#   theme(axis.title.y=element_text(size=12, face= "bold")) +
#   theme(axis.title.x=element_text(size=10))
# 
# Fig5a
# Fig5b <- ggplot(data = subset(ErBAC, ErBAC$Ctn < 20), aes(x = Age.class, y=Ctn, fill=Sex))+
#   geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
#   theme_bw()+scale_shape_manual(values=c(8, 1))+
#   scale_fill_manual(values=c("grey43", "ghostwhite"))+
#   theme(panel.border= element_blank(), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top",
#         axis.line = element_line(linewidth = 0.7, linetype = "solid",
#                                  colour = "black"))+
#   theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
#   
#   theme(axis.title.y=element_text(size=12, face= "bold")) +
#   theme(axis.title.x=element_text(size=10))
# Fig5b
# Fig5c <- ggplot(data = subset(ErBAC, ErBAC$Aldo < 10), aes(x = Age.class, y=Aldo, fill=Sex))+
#   geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
#   theme_bw()+scale_shape_manual(values=c(8, 1))+
#   scale_fill_manual(values=c("grey43", "ghostwhite"))+
#   theme(panel.border= element_blank(), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top",
#         axis.line = element_line(linewidth = 0.7, linetype = "solid",
#                                  colour = "black"))+
#   theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
#   theme(axis.title.y=element_text(size=12, face= "bold")) +
#   theme(axis.title.x=element_text(size=10))
# Fig5c 
# Fig5d <- ggplot(data = ErBAC, aes(x = Age.class, y=Lipid, fill=Sex))+
#   geom_boxplot()+ xlab("")+ ylab("Lipid content %")+
#   theme_bw()+scale_shape_manual(values=c(8, 1))+
#   scale_fill_manual(values=c("grey43", "ghostwhite"))+
#   theme(panel.border= element_blank(), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top",
#         axis.line = element_line(linewidth = 0.7, linetype = "solid",
#                                  colour = "black"))+
#   theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
#   theme(axis.title.y=element_text(size=12, face= "bold")) +
#   theme(axis.title.x=element_text(size=10))
# Fig5_wo<- ggarrange(Fig5a, Fig5b, Fig5c, Fig5d, nrow=2, ncol = 2)
# ggsave(here("./output/Fig_wo_outlier/Age_class.png"),Fig5_wo, width=6, height=6, dpi=300)
# #####Lipid content stats


######## body condition # Adjusted 
ErB$Adjusted_BC
ErB<-ErB %>% drop_na(Adjusted_BC)
plot(ErB$Adjusted_BC, ErB$Cts)
plot(ErB$Adjusted_BC, ErB$Ctn)
plot(ErB$Adjusted_BC, ErB$Aldo)
plot(ErB$Adjusted_BC, ErB$Lipid)

ErB %>%
  # filter(Cts < 30)%>% 
  group_by(Adjusted_BC) %>%
  get_summary_stats(Cts, type = "mean_sd")

ErB %>%
  # filter(Ctn < 20)%>% 
  group_by(Adjusted_BC) %>%
  get_summary_stats(Ctn, type = "mean_sd")

ErB %>%
  # filter(Aldo < 10)%>% 
  group_by(Adjusted_BC) %>%
  get_summary_stats(Aldo, type = "mean_sd")
shapiro.test(log(ErB$Cts)) # 
fit1 <-lm(log(Cts) ~ Adjusted_BC, data= ErB)
summary(fit1)

kruskal.test(Cts ~ Adjusted_BC, data= ErB)
#Kruskal-Wallis chi-squared = 1.3817, df = 2, p-value = 0.5011
kruskal.test(Ctn ~ Adjusted_BC, data= ErB)
#Kruskal-Wallis chi-squared = 1.2974, df = 2, p-value = 0.5227
kruskal.test(Aldo ~ Adjusted_BC, data= ErB)
#Kruskal-Wallis chi-squared = 0.64407, df = 2, p-value = 0.7247
kruskal.test(Lipid ~ Adjusted_BC, data= ErB)
#Kruskal-Wallis chi-squared = 8.1923, df = 2, p-value = 0.01664
?dunn_test
ErB%>% dunn_test(Lipid ~ Adjusted_BC)


#####blar plot with sd
# figure
df1= ErB %>%
  # filter(Cts < 30)%>% 
  mutate(Adjusted_BC = factor(Adjusted_BC,
      levels = c("Poor", "Fair", "Good")))%>%
  group_by(Adjusted_BC) %>%
  get_summary_stats(Cts, type = "mean_sd")

Fig3a <- ggplot(df1, aes(x=Adjusted_BC, y=mean)) + 
  geom_bar(stat="identity",color= "black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ylab("Cortisol ng/g")+ xlab("")+
  theme_minimal()+
  
Fig3a

df2= ErB %>%
  # filter(Cts < 30)%>% 
  mutate(Adjusted_BC = factor(Adjusted_BC,
  levels = c("Poor", "Fair", "Good")))%>%
  group_by(Adjusted_BC) %>%
  get_summary_stats(Ctn, type = "mean_sd")
Fig3b <- ggplot(df2, aes(x=Adjusted_BC, y=mean)) + 
  geom_bar(stat="identity",color= "black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ylab("Corticosterone ng/g")+ xlab("")+
  theme_minimal()
# geom_signif(comparisons = list(c("Adult", "Immature")),
#             map_signif_level = TRUE,
#             annotations = c("*"))
Fig3b 
df3= ErB %>%
  # filter(Cts < 30)%>% 
  mutate(Adjusted_BC = factor(Adjusted_BC,
  levels = c("Poor", "Fair", "Good")))%>%
  group_by(Adjusted_BC) %>%
  get_summary_stats(Aldo, type = "mean_sd")
Fig3c <- ggplot(df3, aes(x=Adjusted_BC, y=mean)) + 
  geom_bar(stat="identity",color= "black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ylab("Aldosterone ng/g")+ xlab("")+
  theme_minimal()
Fig3c
df4= ErB %>%
  # filter(Cts < 30)%>% 
  mutate(Adjusted_BC = factor(Adjusted_BC,
                              levels = c("Poor", "Fair", "Good")))%>%
  group_by(Adjusted_BC) %>%
  get_summary_stats(Lipid, type = "mean_sd")
Fig4d <- ggplot(df4, aes(x=Adjusted_BC, y=mean)) + 
  geom_bar(stat="identity",color= "black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  ylab("Lipid content %")+ xlab("")+
  theme_minimal()
Fig4d

Ageclass <- ggarrange(Fig2a, Fig2b, Fig2c,Fig2d, nrow=2, ncol = 2)

here()
ggsave(here("./output/Ageclass.png"),Ageclass, width=8, height=8, dpi=600)

Fig4a <- ggplot(data = ErB, aes(x = Adjusted_BC, y=Cts, fill= Sex))+
  geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g, log scale")+
  # geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
  #              width = .75, linetype = "dashed")+
  scale_y_continuous(trans='log10')+ annotation_logticks(sides= "l")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))

Fig4a
Fig4b <- ggplot(data = ErB1, aes(x = Adjusted_BC, y=Ctn, fill= Sex))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g, log scale")+
  # geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
  #              width = .75, linetype = "dashed")+
  scale_y_continuous(trans='log10')+ annotation_logticks(sides= "l")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))
Fig4b

Fig4c <- ggplot(data = ErB1, aes(x = Adjusted_BC, y=Aldo, fill= Sex))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g, log scale")+
  # geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
  #              width = .75, linetype = "dashed")+
  scale_y_continuous(trans='log10')+ annotation_logticks(sides= "l")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))
Fig4c 
Fig4d <- ggplot(data = ErB1, aes(x = Adjusted_BC, y=Lipid, fill= Sex))+
  geom_boxplot()+ xlab("")+ ylab("Lipid content %")+
  # geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
  #              width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))
Fig4d
Body_Condition_with <- ggarrange(Fig4a, Fig4b, Fig4c, Fig4d, nrow=2, ncol = 2)

ggsave(here("./output/Fig_with_outlier/BodyCondition.png"),Body_Condition_with, width=6, height=6, dpi=300)


# #### without outlier
# Fig4a <- ggplot(data = subset(ErB, ErB$Cts < 30), aes(x = Adjusted_BC, y=Cts))+
# geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g")+
#   geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
#   stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
#                width = .75, linetype = "dashed")+
#   theme_bw()+scale_shape_manual(values=c(8, 1))+
#   theme(panel.border= element_blank(), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top",
#         axis.line = element_line(linewidth = 0.7, linetype = "solid",
#                                  colour = "black"))+
#   theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
#   
#   theme(axis.title.y=element_text(size=12, face= "bold")) +
#   theme(axis.title.x=element_text(size=10))
# 
# Fig4a
# Fig4b <- ggplot(data = subset(ErB, ErB$Ctn <20), aes(x = Adjusted_BC, y=Ctn))+
#   geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
#   geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
#   stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
#                width = .75, linetype = "dashed")+
#   theme_bw()+scale_shape_manual(values=c(8, 1))+
#   theme(panel.border= element_blank(), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top",
#         axis.line = element_line(linewidth = 0.7, linetype = "solid",
#                                  colour = "black"))+
#   theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
#   
#   theme(axis.title.y=element_text(size=12, face= "bold")) +
#   theme(axis.title.x=element_text(size=10))
# Fig4b
# Fig4c <- ggplot(data = subset(ErB, ErB$Aldo <10), aes(x = Adjusted_BC, y=Aldo))+
#   geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
#   geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
#   stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
#                width = .75, linetype = "dashed")+
#   theme_bw()+scale_shape_manual(values=c(8, 1))+
#   theme(panel.border= element_blank(), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top",
#         axis.line = element_line(linewidth = 0.7, linetype = "solid",
#                                  colour = "black"))+
#   theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
#   
#   theme(axis.title.y=element_text(size=12, face= "bold")) +
#   theme(axis.title.x=element_text(size=10))
# Fig4c 
# Fig4d <- ggplot(data = ErB, aes(x = Adjusted_BC, y=Lipid))+
#   geom_boxplot()+ xlab("")+ ylab("Lipid content %")+
#   geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
#   stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
#                width = .75, linetype = "dashed")+
#   theme_bw()+scale_shape_manual(values=c(8, 1))+
#   theme(panel.border= element_blank(), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "top",
#         axis.line = element_line(linewidth = 0.7, linetype = "solid",
#                                  colour = "black"))+
#   theme(axis.text.x= element_text(size=10, face= "bold", angle = 45, hjust = 1))+ 
#   
#   theme(axis.title.y=element_text(size=12, face= "bold")) +
#   theme(axis.title.x=element_text(size=10))
# Fig4d
# Body_Condition_wo <- ggarrange(Fig4a, Fig4b, Fig4c, Fig4d, nrow=2, ncol = 2)
# 
# ggsave(here("./output/Fig_wo_outlier/BodyCondition.png"),Body_Condition_wo, width=6, height=6, dpi=300)
# 
# kruskal.test(Cts ~ Adjusted_BC, data= subset(ErB, ErB$Cts <30))
# # Kruskal-Wallis chi-squared = 0.79957, df = 2, p-value = 0.6705
# kruskal.test(Ctn ~ Adjusted_BC, data= subset(ErB, ErB$Ctn <20))
# # Kruskal-Wallis chi-squared = 0.81357, df = 2, p-value = 0.6658
# kruskal.test(Aldo ~ Adjusted_BC, data= subset(ErB, ErB$Aldo <10))
# # Kruskal-Wallis chi-squared = 1.2571, df = 2, p-value = 0.5334
# kruskal.test(Lipid ~ Adjusted_BC, data= ErB)
# #Kruskal-Wallis chi-squared = 8.1923, df = 2, p-value = 0.01664
# ?pairwise_wilcox.test
# pairwise.wilcox.test(ErB$Lipid, ErB$Adjusted_BC, p.adjust.method= "bonf")
# 
# 
# 

#COD
#cortisol

ErBCTS <- ErB[!is.na(ErB$Cts),]
tapply(ErBCTS$Cts, FUN = mean , INDEX = ErBCTS$COD)
tapply(ErBCTS$Cts, FUN = sd , INDEX = ErBCTS$COD)
summary(ErBCTS)
performance::check_outliers(ErBCTS$Cts)
ErBCTS[c(4,5,20,21),]

hist(ErBCTS$Cts)
hist(log(ErBCTS$Cts))
performance::check_outliers(log(ErBCTS$Cts))
shapiro.test(log(ErBCTS$Cts)) # p=0.62
test<- aov(log(Cts) ~ COD, data= ErBCTS)
plot(test)
summary(test)
shapiro.test(resid(test)) # 0.48
summary(test) # no significant difference detected
kruskal.test(Cts ~ COD, data= ErBCTS)
###Kruskal-Wallis chi-squared = 5.1399, df = 4, p-value = 0.2732
kruskal.test(Cts ~ COD, data= subset(ErBCTS, ErBCTS$Cts <30))
#Kruskal-Wallis chi-squared = 5.1809, df = 4, p-value = 0.2692
#corticosterone

ErBCTN <- ErB[!is.na(ErB$Ctn),]
tapply(ErBCTN$Ctn, FUN = mean , INDEX = ErBCTN$COD)
tapply(ErBCTN$Ctn, FUN = sd , INDEX = ErBCTN$COD)

hist(ErBCTN$Ctn)
hist(log(ErBCTN$Ctn))
shapiro.test(log(ErBCTN$Ctn))
#log transformed give you a normal distribution
test1<- lm(log(Ctn) ~ COD, data= ErBCTN)
summary(test1) #no significant difference
shapiro.test(resid(test1)) # p= 0.9

 # no significant difference detected
performance::check_outliers(ErBCTN$Ctn)
ErBCTN[c(3,9,25,33),]


hist(ErBCTN$Ctn)
hist(log(ErBCTN$Ctn))
shapiro.test(log(ErBCTN$Ctn))
#log transformed give you a normal distribution
levene_test(ErBCTN$COD, ErBCTN$Ctn)
test1<- lm(log(Ctn) ~ COD, data= ErBCTN)
summary(test1)
shapiro.test(resid(test1)) #yes
fit1<- aov(log(Ctn) ~ COD, data= ErBCTN)
plot(fit1)
summary(fit1)
kruskal.test(Ctn ~ COD, data= ErBCTN)
#Kruskal-Wallis chi-squared = 5.498, df = 4, p-value = 0.2399
no0 <-subset(ErBCTN, ErBCTN$Ctn < 20)
tapply(no0$Ctn, FUN = mean , INDEX = no0$COD)
tapply(no0$Ctn, FUN = sd , INDEX = no0$COD)
levene_test(no0$COD,no0$Ctn) 
test1<- lm(log(Ctn) ~ COD, data= no0)
summary(test1)
shapiro.test(resid(test1)) #yes
kruskal.test(Ctn ~ COD, data= no0)
#Kruskal-Wallis chi-squared = 6.0413, df = 4, p-value = 0.1961


ErBALDO <- ErB[!is.na(ErB$Aldo),]
performance::check_outliers(ErBALDO$Aldo)
ErBALDO[c(3,18,24,31,34),]


hist(ErBALDO$Aldo)
hist(log(ErBALDO$Aldo))
shapiro.test(log(ErBALDO$Aldo))
qqnorm(log(ErBALDO$Aldo))
#log transformed give you a normal distribution
levene_test(ErBALDO$COD, log(ErBALDO$Aldo))
kruskal.test(Aldo ~ COD, data= ErBALDO)
#Kruskal-Wallis chi-squared = 4.4635, df = 4, p-value = 0.3469
test2<- lm(log(Aldo) ~ COD, data= subset(ErBALDO, ErBALDO$Aldo <10))
summary(test2)
shapiro.test(resid(test2))

kruskal.test(Aldo ~ COD, data= subset(ErBALDO, ErBALDO$Aldo <10))

# ErBALDO <-subset(ErBALDO, ErBALDO$Aldo <10)
tapply(ErBALDO$Aldo, FUN = mean , INDEX = ErBALDO$COD)
tapply(ErBALDO$Aldo, FUN = sd , INDEX = ErBALDO$COD)




#no difference

Fig3a <-ggplot(data = ErBCTS, aes(x = COD, y=Cts))+
  geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_blank())+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))


Fig3b <-ggplot(data = ErBCTN, aes(x = COD, y=Ctn))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_blank())+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))
Fig3c <-ggplot(data = ErBALDO, aes(x = COD, y=Aldo))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=8, face= "bold"))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))


COD_with <- ggarrange(Fig3a, Fig3b, Fig3c, nrow=3, ncol = 1)

ggsave(here("./output/Fig_with_outlier/COD.png"),COD_with, width=9, height=6, dpi=300)

##########wihtout outlier

Fig3a <-ggplot(data = subset(ErBCTS, ErBCTS$Cts <30), aes(x = COD, y=Cts))+
  geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_blank())+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))


Fig3b <-ggplot(data = subset(ErBCTN, ErBCTN$Ctn < 20), aes(x = COD, y=Ctn))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_blank())+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))

Fig3c <-ggplot(data = subset(ErBALDO, ErBALDO$Aldo <10), aes(x = COD, y=Aldo))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=8, face= "bold"))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))


COD_wo <- ggarrange(Fig3a, Fig3b, Fig3c, nrow=3, ncol = 1)
COD_wo 
ggsave(here("./output/Fig_wo_outlier/COD.png"),COD_wo, width=9, height=6, dpi=300)


######
########## PRESENTATION MADEIRA 12/04/2026
# redoing this figure to a barplot
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


df_CTS_no <- ErBCTS %>%
  # filter(ErBCTS$Cts <30 ) %>%
  group_by(COD) %>%
  summarise(meanCTS =mean(Cts), 
            sdCTS = sd(Cts), 
            n= n())


p <- ggplot(df_CTS_no , aes(x=COD, y=meanCTS, fill= COD) )+ 
  geom_bar(position=position_dodge(.9), stat = "identity",
           color="black", show_guide=FALSE)+
  geom_errorbar(aes(ymin=meanCTS, ymax=meanCTS+sdCTS), position=position_dodge(.9), width = 0.25)
p

final <- p+ theme_bw()+
  scale_x_discrete(labels=c("Entanglement (n=2)", "Nutr stress (n=12)", "Predation (n=3)", 
                            "Trauma (n=10)", "Undetermined (n=7)")) + ylab("Glucocorticoids (ng/g)")+
  theme(axis.text.x= element_text(size=12))+
  theme(axis.text.y= element_text(size=12))+
  theme(axis.title.y=element_text(size=14, face= "bold")) +
  theme(axis.title.x=element_blank())+
  scale_fill_brewer(palette="Blues") 

ggsave(here("./output/Fig_with_outlier/COD_madeira.png"),final, width=9, height=6, dpi=600)



Fig3a <-ggplot(data = subset(ErBCTS, ErBCTS$Cts <30), aes(x = COD, y=Cts))+
  geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "top",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_blank())+ 
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))


Fig3b <-ggplot(data = subset(ErBCTN, ErBCTN$Ctn < 20), aes(x = COD, y=Ctn))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_blank())+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))

Fig3c <-ggplot(data = subset(ErBALDO, ErBALDO$Aldo <10), aes(x = COD, y=Aldo))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=8, face= "bold"))+ 
  
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))


COD_wo <- ggarrange(Fig3a, Fig3b, Fig3c, nrow=3, ncol = 1)
COD_wo 
ggsave(here("./output/Fig_wo_outlier/COD.png"),COD_wo, width=9, height=6, dpi=300)




###### carcass condition


plot(ErB$CarCode, ErB$Cts)
plot(ErB$CarCode, ErB$Ctn)
plot(ErB$CarCode, ErB$Aldo)

ErB %>%
  group_by(CarCode) %>%  
  dplyr::summarize(n = n(),
                   meanCTS = mean(Cts, na.rm=T),
                   sdCTS = sd(Cts, na.rm=T),
                   rangeCTS = max(Cts, na.rm=T) - min(Cts,na.rm=T),
                   meanCTN = mean(Ctn, na.rm=T),
                   sdCTN = sd(Ctn, na.rm=T),
                   rangeCTN = max(Ctn, na.rm=T) - min(Ctn, na.rm=T),
                   meanALDO = mean(Aldo, na.rm=T),
                   sdALDO = sd(Aldo, na.rm=T),
                   rangeALDO = max(Aldo, na.rm=T) - min(Aldo, na.rm=T))

hist(ErB$Cts)
hist(log(ErB$Cts))
performance::check_outliers(log(ErB$Cts))
shapiro.test(ErB$Cts) # signficant
bartlett.test(ErB$CarCode,ErB$Cts)
kruskal.test(Cts ~ CarCode, data= ErB)
wilcox_test(Cts ~ CarCode, data= ErB)
# no significant 
kruskal.test(Ctn ~ CarCode, data= ErB)
# no significant 
kruskal.test(Aldo ~ CarCode, data= ErB)
# no significant 
kruskal.test(Lipid ~ CarCode, data= ErB)
# no significant 












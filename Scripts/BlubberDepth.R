getwd() # LINE 250
#this file is just for the analysis of hormone and lipid concentrations in response to blubber depth. 
#I have created a dataset that have samples that have data for all blubber depth (surface, medium and bottom) for either hormones. 19 individuals total
rm(list=ls())
  

# validation after gray whale analysis
library(plotrix)
library(MASS)
library(plyr)
library(lattice)
library(mixtools)
library(MuMIn)
library(dplyr) # data manipulation package
library(ggplot2) # data visualization package
library(gapminder) # package with a gapminder dataset
library(forcats) # package for working with categorical data
library(plotly) # package for interactive data visualizations
library(tidyr)
library(emmeans)
library(ggpubr)
library(gplots) # to plot means
library(tidyverse)
library(ggpubr)
library(rstatix)
library(stringr)
library(patchwork)
library(Hmisc)
library(ggforce)
library(lme4)
library(here)
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

#Gray whale analysis in relation to blubber depth 
Er<- read.csv(here ("./data/GW_BLB.csv"), header=TRUE)
str(Er)
colnames(Er) <- c("vial", "MMSNW" , "ID","BlubberDEPTH","Sex", "Age.class", "BoCo" ,"CarCode"   , "BlubberWeight" ,      
                     "Cts", "Cts CV", "Ctn" ,  "Ctn CV",  "Aldo", "Aldo CV", "LIPID")

Er$MMSNW = as.factor(Er$MMSNW)
Er$ID = as.factor(Er$ID)
Er$BlubberDEPTH = as.factor(Er$BlubberDEPTH)
summary(Er$BlubberDEPTH)
Er$Age.class= as.factor(Er$Age.class)
summary(Er$Age.class)
# Er$BlubberDEPTH <- ordered(Er$BlubberDEPTH, levels=c("Surface", "Medium", "Bottom"))
Er$BoCo =as.factor(Er$BoCo)
levels(Er$BoCo)

Er$CarCode = as.factor(Er$CarCode)
Er_table <- Er %>%
  select(-c(vial, `Cts CV`, `Ctn CV`, `Aldo CV`))

length(unique(Er$ID))

Sum_ER <-Er %>%
  group_by(Sex, BlubberDEPTH) %>%  
  dplyr::summarize(n = n(),
                   meanCTS = mean(Cts, na.rm=T),
                   sdCTS = sd(Cts, na.rm=T),
                   maxCTS = max(Cts, na.rm=T), 
                   minCTS = min(Cts,na.rm=T),
                   meanCTN = mean(Ctn, na.rm=T),
                   sdCTN = sd(Ctn, na.rm=T),
                   maxCTN = max(Ctn, na.rm=T), 
                   minCTN = min(Ctn,na.rm=T),
                   meanALDO = mean(Aldo, na.rm=T),
                   sdALDO = sd(Aldo, na.rm=T),
                   maxALDO = max(Aldo, na.rm=T),
                   minALDO = min(Aldo,na.rm=T),
                   meanLip = mean(LIPID, na.rm=T),
                   sdLip = sd(LIPID, na.rm=T),
                   maxLip= max(LIPID, na.rm=T),
                   minLip =min(LIPID, na.rm=T))

view(Sum_ER)
write.csv(Sum_ER, here("./output/Sum_ER.csv"), row.names = FALSE)
range(Er$BlubberWeight) #0.096-0.57
mean(Er$BlubberWeight) #0.22
sd(Er$BlubberWeight) #0.10
is.na(Er$Cts)
is.na(Er$Ctn)
#EDA
plot.new()
par(mfrow=c(1,1))
plot(Er$BlubberDEPTH, Er$Cts, xlab="Blubber depth", ylab="Cortisol ng/g")
plot(Er$BlubberDEPTH, Er$Ctn, xlab= "Blubber depth", ylab="Corticosterone ng/g")
plot(Er$BlubberDEPTH, Er$LIPID, xlab="Blubber depth", ylab="Lipid %")

##### correlations
#lipid
cor.test(Er$LIPID,Er$Cts) #0.12
cor.test(Er$LIPID,Er$Ctn)#0.15
cor.test(Er$LIPID,Er$Aldo) #0.28
#amongst hormones
cor.test(Er$Ctn,Er$Cts) # 0.49
cor.test(Er$Aldo,Er$Cts) # 0.05
cor.test(Er$Ctn,Er$Aldo) # 0.65
Er1 <- subset(Er, Er$Ctn <3) #0.55
cor.test(Er1$Ctn,Er1$Aldo)
plot(Er$Ctn, Er$Cts, xlab="Corticosterone", ylab="Cortisol")
plot(Er1$Ctn, Er1$Aldo, xlab="Corticosterone", ylab="Aldosterone")
plot(Er$Cts, Er$Aldo, xlab="Cortisol", ylab="Aldosterone")

#EDA hormone concentrations across blubber depth

#Cortisol and figure 2a

ErCTS <- Er[!is.na(Er$Cts),]
summary(ErCTS$BlubberDEPTH)
ggplot(ErCTS) + geom_histogram(aes(log(Cts)))
bartlett.test(log(ErCTS$Cts), ErCTS$BlubberDEPTH)

fit2 <- lm(log(Cts) ~ BlubberDEPTH, data = ErCTS)
plot(fit2)
# 
# performance::check_outliers(ErCTS$Cts)
# ErCTS[c(12:14,20,42),]
# summary(ErCTS$Sex)
model_Cts1 <- lmer(log(Cts) ~ BlubberDEPTH + (1|ID) - 1, data = ErCTS)
model_Cts2 <- lmer(log(Cts) ~ BlubberDEPTH:Sex + (1|ID) - 1, data = ErCTS)
model_GC_null <- lmer(log(Cts)~ (1|ID) - 1, data = ErCTS)
MuMIn::AICc(model_GC_null,model_Cts1,model_Cts2)
anova(model_GC_null,model_Cts1,model_Cts2)
posthoc_Cts1 <- emmeans(model_Cts1,  pairwise~ BlubberDEPTH)
summary(posthoc_Cts1) # significantly different between surface and medium layer
# surface and medium sign diff
ErCTS$BlubberDEPTH<- ordered(ErCTS$BlubberDEPTH, levels= c("Surface", "Medium", "Bottom"))
shapes = c("F" = 24, "M" = 21)
cols = c("F" = '#999999', "M" = '#E69F00')

library(scales)
Fig2a <- ggplot(data=ErCTS, aes(y=Cts, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g")+
  geom_jitter(aes(shape=Sex, fill=Sex),position=position_jitter(0.2), size=2)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  scale_shape_manual(values = shapes, name= "Sex")+
  scale_fill_manual(values = cols, name= "Sex")+
  theme_bw()+
  scale_y_continuous(trans='log10')+ annotation_logticks(sides= "l")+
  annotate("text",
           x = -Inf,
           y = Inf,
           label = "A",
           hjust = 2.5,
           vjust = 1,
           size = 4,
           fontface = "bold")+
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "top",
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 10, 10, 20)  # extra space on left
  )

Fig2a

###########CONTINUE HERE


# no significant difference; however the concentrations in blubber cortisol between the surface 
# and medium layer is barely non significant
#filtered for outliers samples with Cts concentrations lowet than 20 ng/g


# filtered <- ErCTS |> filter(Cts<25)
# #filtered <- ErCTS |> filter(ID !="1524")
# #filtered <- ErCTS |> filter(ID !="SKMMR2019-5-Er")
# model_Cts_filtered <- lmer(Cts ~ BlubberDEPTH + (1|ID) - 1, data = filtered)
# summary(model_Cts_filtered)
# model_GC_null <- lmer(Cts~ (1|ID) - 1, data = filtered)
# MuMIn::AICc(model_GC_null,model_Cts_filtered)
# posthoc_Cts_filtered <- emmeans(model_Cts_filtered, pairwise~BlubberDEPTH )
# summary(posthoc_Cts_filtered)

#no significant difference



#For cortisol:
#17 individuals with complete profile (i.e. hormone concentrations for all three layers) + 2 individuals with concentrations from 2 out of three layers. 
#To test the effects of blubber depth on cortisol concentrations, I fitted a mixed effects model where individual are counted as random effects ( repeated measures ) and blubber depth (and sex) where tested as a fixed effects. 
#The output of the models indicated that the model including the interaction of blubber depth and sex to be the best model (lowest AICc). 
#However, the post hoc test Emmeans did not indicate any significant difference in any combination of blubber depth and sex. Note: what there are two major outliers in the surface layers (both females).
CTS_sum <- ErCTS %>%
  group_by(BlubberDEPTH, Sex) %>%  
  dplyr::summarize(n = n(),
                   mean = mean(Cts),
                   sd = sd(Cts))

# view(CTS_sum)
# write.csv(CTS_sum, "/Users/valentinamelica/R/GW/data", row.names = FALSE)

#Corticosterone and fig 2b
ErCTN <- Er[!is.na(Er$Ctn),]
ErCTN <- ErCTN |> filter(ID !="180413") #remove ID that has only surface
#ErCTN <- subset(ErCTN, ErCTN$CTN<20.0) #remove outlayer
summary(ErCTN$BlubberDEPTH)
ggplot(ErCTN) + geom_histogram(aes(log(Ctn))) 
bartlett.test(log(ErCTN$Ctn), ErCTN$BlubberDEPTH)
fit2 <- lm(log(Ctn) ~ BlubberDEPTH, data = ErCTN)
plot(fit2)

# performance::check_outliers(ErCTN$Ctn)
# ErCTN[c(13,19,34,46),]
tapply(ErCTN$Ctn, FUN=mean, INDEX= ErCTN$BlubberDEPTH)
#Surface   Medium  Bottom 
#1.08      0.78    0.66 
tapply(ErCTN$Ctn, FUN=range, INDEX= ErCTN$BlubberDEPTH)
#Bottom# #0.07 1.7 #Medium #0.15 1.87 #Surface #0.16 3.17
model_Ctn1<- lmer(log(Ctn) ~ BlubberDEPTH + (1|ID) - 1, data = ErCTN) 
model_Ctn2<- lmer(log(Ctn) ~ BlubberDEPTH:Sex+ (1|ID) - 1, data = ErCTN)
#model_Ctn2<- lmer(Ctn ~ BlubberDEPTH*Sex+ (1|ID) - 1, data = ErCTN)
model_CTN_null <- lmer(log(Ctn)~ (1|ID) - 1, data = ErCTN)
MuMIn::AICc(model_CTN_null,model_Ctn1, model_Ctn2)
anova(model_CTN_null,model_Ctn1, model_Ctn2)
#the null model win
summary(model_Ctn2)
posthoc_Ctn2 <- emmeans(model_Ctn2, pairwise~BlubberDEPTH:Sex )
summary(posthoc_Ctn2) # no statistical difference


# model_Ctn<- lmer(Ctn ~ BlubberDEPTH + (1|ID) - 1, data = subset(ErCTN, ErCTN$Ctn < 2.2)) 
# #model_Ctn1<- lmer(Ctn ~ BlubberDEPTH + Sex+ (1|ID) - 1, data = ErCTN)
# #model_Ctn2<- lmer(Ctn ~ BlubberDEPTH*Sex+ (1|ID) - 1, data = ErCTN)
# model_CTN_null <- lmer(Ctn~ (1|ID) - 1, data = subset(ErCTN, ErCTN$Ctn < 2.2))
# MuMIn::AICc(model_CTN_null,model_Ctn)
# 
# #the null model win
# summary(model_Ctn)
# posthoc_Ctn <- emmeans(model_Ctn, pairwise~BlubberDEPTH, adjust= "bonferroni" )
# summary(posthoc_Ctn) # no statistical difference

levels(ErCTS$BlubberDEPTH) 

ErCTN$BlubberDEPTH<- ordered(ErCTN$BlubberDEPTH, levels= c("Surface", "Medium", "Bottom"))

Fig2b <- ggplot(data=ErCTN, aes(y=Ctn, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  geom_jitter(aes(shape=Sex, fill=Sex),position=position_jitter(0.2), size=2)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  scale_shape_manual(values = shapes, name= "Sex")+
  scale_fill_manual(values = cols, name= "Sex")+
  theme_bw()+
  scale_y_continuous(trans='log10')+ annotation_logticks(sides= "l")+
  annotate("text",
           x = -Inf,
           y = Inf,
           label = "B",
           hjust = 2.5,
           vjust = 1,
           size = 4,
           fontface = "bold")+
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "top",
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 10, 10, 20)  # extra space on left
  )

Fig2b

#Aldosterone and Fig 2c
ErALDO <- Er[!is.na(Er$Aldo),]
ErALDO<- ErALDO |> filter(ID !="1712") #remove ID that has only surface
#ErCTN <- subset(ErCTN, ErCTN$CTN<20.0) #remove outlayer
summary(ErALDO$BlubberDEPTH)
ggplot(ErALDO) + geom_histogram(aes(log(Aldo))) 
qqnorm(log(ErALDO$Aldo))
bartlett.test(log(ErALDO$Aldo), ErALDO$BlubberDEPTH)
fit2 <- lm(log(Aldo) ~ BlubberDEPTH, data = ErALDO)
plot(fit2)
# performance::check_outliers(ErALDO$Aldo)
# ErALDO[c(1,8,9,31),]

tapply(ErALDO$Aldo, FUN=mean, INDEX= ErALDO$BlubberDEPTH)
#Surface   Medium  Bottom
#0.23      0.06    0.08
tapply(ErALDO$Aldo, FUN=range, INDEX= ErALDO$BlubberDEPTH)
#Bottom#
#0.003 0.43
#Medium
#0.01 0.26
#Surface
#0.04 1.52
model_Aldo1<- lmer(log(Aldo) ~ BlubberDEPTH + (1|ID) - 1, data = ErALDO)
model_Aldo2<- lmer(Aldo ~ BlubberDEPTH:Sex+ (1|ID) - 1, data = ErALDO)
model_ALDO_null <- lmer(log(Aldo)~ (1|ID) - 1, data = ErALDO)
MuMIn::AICc(model_ALDO_null,model_Aldo1, model_Aldo2)
anova(model_ALDO_null, model_Aldo1, model_Aldo2)
summary(model_Aldo1)
posthoc_Aldo1 <- emmeans(model_Aldo1, pairwise~BlubberDEPTH )
summary(posthoc_Aldo1) # no statistical difference

ErALDO$BlubberDEPTH<- ordered(ErALDO$BlubberDEPTH, levels= c("Surface", "Medium", "Bottom"))

shapes = c("F" = 24, "M" = 21)
cols = c("F" = '#999999', "M" = '#E69F00')

Fig2c <- ggplot(data=ErALDO, aes(y=Aldo, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  geom_jitter(aes(shape=Sex, fill=Sex),position=position_jitter(0.2), size=2)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  scale_shape_manual(values = shapes, name= "Sex")+
  scale_fill_manual(values = cols, name= "Sex")+
  theme_bw()+
  scale_y_continuous(trans='log10')+ annotation_logticks(sides= "l")+
  annotate("text",
           x = -Inf,
           y = Inf,
           label = "C",
           hjust = 2.5,
           vjust = 1,
           size = 4,
           fontface = "bold")+
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "top",
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 10, 10, 20)  # extra space on left
  )
Fig2c

# # remove of outlier
# filteredALDO <- ErALDO |> filter(Aldo<1)
# model_Aldofiltered<- lmer(Aldo ~ BlubberDEPTH + (1|ID) - 1, data = filteredALDO)
# model_ALDO_filtered_null <- lmer(Aldo~ (1|ID) - 1, data = filteredALDO)
# MuMIn::AICc(model_ALDO_filtered_null,model_Aldofiltered)
# summary(model_Aldo)
# posthoc_Aldo <- emmeans(model_Aldo, pairwise~BlubberDEPTH, adjust= "bonferroni" )
# summary(posthoc_Aldo)


###########lipid content

ErLipid <- Er[!is.na(Er$LIPID),]
summary(ErLipid$BlubberDEPTH)
ggplot(ErLipid) + geom_histogram(aes(log(LIPID))) 
fit2 <- lm(log(LIPID)~ BlubberDEPTH, data= ErLipid)
plot(fit2)

performance::check_outliers(ErLipid$LIPID)
ErLipid[c(17:19, 28, 29, 31),] #outliers above 50%
bartlett.test(log(ErLipid$LIPID), ErLipid$BlubberDEPTH)

model_Lipid<- lmer(log(LIPID) ~ BlubberDEPTH + (1|ID) - 1, data = ErLipid)
# model_Lipid2<- lmer(log(LIPID) ~ BlubberDEPTH:Sex+ (1|ID) - 1, data = ErLipid)
model_Lipid_null <- lmer(log(LIPID)~ (1|ID) - 1, data = ErLipid)
MuMIn::AICc(model_Lipid_null,model_Lipid)
anova(model_Lipid_null,model_Lipid)

summary(model_Lipid)
posthoc_Lipid <- emmeans(model_Lipid, pairwise~BlubberDEPTH )
summary(posthoc_Lipid) # surface layer is significantly different from the other 2

# model_Lipid<- lmer(LIPID ~ BlubberDEPTH + (1|ID) - 1, data = subset(ErLipid, ErLipid$LIPID < 50))
# #model_Lipid2<- lmer(LIPID ~ BlubberDEPTH:Sex+ (1|ID) - 1, data = ErLipid)
# model_Lipid_null <- lmer(LIPID~ (1|ID) - 1, data = subset(ErLipid, ErLipid$LIPID < 50))
# MuMIn::AICc(model_Lipid_null,model_Lipid)
# summary(model_Lipid)
# posthoc_Lipid <- emmeans(model_Lipid, pairwise~BlubberDEPTH, adjust= "bonferroni" )
# summary(posthoc_Lipid)
ErLipid$BlubberDEPTH<- ordered(ErLipid$BlubberDEPTH, levels= c("Surface", "Medium", "Bottom"))
shapes = c("F" = 24, "M" = 21)
cols = c("F" = '#999999', "M" = '#E69F00')

Fig2d <- ggplot(data=ErLipid, aes(y=LIPID, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Lipid content %")+
  geom_jitter(aes(shape=Sex, fill=Sex),position=position_jitter(0.2), size=2)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  scale_shape_manual(values = shapes, name= "Sex")+
  scale_fill_manual(values = cols, name= "Sex")+
  theme_bw()+
  scale_y_continuous(trans='log10')+ annotation_logticks(sides= "l")+
  annotate("text",
           x = -Inf,
           y = Inf,
           label = "D",
           hjust = 2.5,
           vjust = 1,
           size = 4,
           fontface = "bold")+
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "top",
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 10, 10, 20)  # extra space on left
  )

Fig2d

BlbDepth_with <- ggarrange(Fig2a, Fig2b, Fig2c,Fig2d, nrow=2, ncol = 2)
BlbDepth_with
ggsave(here("./output/Blubberdepth.png"),BlbDepth_with, width=10, height=10, dpi=600)

####without outlier
fig2a <- ggplot(data=subset(ErCTS, ErCTS$Cts <25), aes(y=Cts, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Cortisol ng/g")+
  geom_jitter(aes(shape=Sex),position=position_jitter(0.2), size=2)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+
  scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth= 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=12, face= "bold"))+
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))
Fig2a<-fig2a +coord_flip()
Fig2a

ErCTN$BlubberDEPTH<- ordered(ErCTN$BlubberDEPTH, levels= c("Bottom", "Medium", "Surface"))


Fig2b <- ggplot(data=ErCTN, aes(y=Ctn, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  geom_jitter(aes(shape=Sex),position=position_jitter(0.2), size=2)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+
  scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth= 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=12, face= "bold"))+
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))

Fig2b <-Fig2b +coord_flip()

Fig2b

# filteredALDO$BlubberDEPTH<- ordered(filteredALDO$BlubberDEPTH, levels= c("Bottom", "Medium", "Surface"))

ErALDO$BlubberDEPTH<- ordered(ErALDO$BlubberDEPTH, levels= c("Bottom", "Medium", "Surface"))

Fig2c<- ggplot(data=subset(ErALDO,ErALDO$Aldo < 1), aes(y=Aldo, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  geom_jitter(aes(shape=Sex),position=position_jitter(0.2), size=2)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+
  scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "noe",
        axis.line = element_line(linewidth= 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=12, face= "bold"))+
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))

Fig2c <-Fig2c +coord_flip()
Fig2c




ErLipid$BlubberDEPTH<- ordered(ErLipid$BlubberDEPTH, levels= c("Bottom", "Medium", "Surface"))
fig2d<- ggplot(data=ErLipid, aes(y=LIPID, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Lipid %")+
  geom_jitter(aes(shape=Sex),position=position_jitter(0.2), size=2)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed")+
  theme_bw()+
  scale_shape_manual(values=c(8, 1))+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth= 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=12, face= "bold"))+
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))

Fig2d <-fig2d +coord_flip()
Fig2d

BlbDepth_wo <- ggarrange(Fig2a, Fig2b, Fig2c,Fig2d, nrow=2, ncol = 2)
BlbDepth_wo
ggsave(here("./output/Fig_wo_outlier/Blubberdepth.png"),BlbDepth_wo, width=6, height=6, dpi=300)



#############old code

hist(ErCTS$CTN) #left skewed
hist(log(ErCTN$CTN))
range(ErCTN$CTN) #0.07-21.66
shapiro.test(ErCTN$CTN) #p<0.001
shapiro.test(log(ErCTN$CTN)) #not significant
qqnorm(log(ErCTN$CTN))

bartlett.test(ErCTN$CTN, ErCTN$Sex) #not homogeneus
summary(ErCTN$Age.class)
bartlett.test(ErCTN$CTN, ErCTN$Age.class) #not homogeneus variances
summary(ErCTN$BlubberDEPTH)
bartlett.test(ErCTN$CTN, ErCTN$BlubberDEPTH) #non equal variances

kruskal.test(Ctn ~ BlubberDEPTH, data=ErCTN)

#Kruskal-Wallis chi-squared = 5.855, df = 2, p-value = 0.05353 considering all data

kruskal.test(CTN ~ BlubberDEPTH, data=subset(ErCTN, ErCTN$CTN<20.0))
#Kruskal-Wallis chi-squared = 5.2589, df = 2,p-value = 0.07212

(fitCTN <-aov(log(CTN) ~ BlubberDEPTH, data=subset(ErCTN, ErCTN$CTN<20.0)))
summary(fitCTN) #p =0.041
TukeyHSD(fitCTN)
#surface and bottom are significantly different p=0.03
tapply(ErCTN$CTN, FUN= mean, INDEX= ErCTN$BlubberDEPTH)
tapply(ErCTN$CTN, FUN= range, INDEX= ErCTN$BlubberDEPTH)

#Bottom    Medium   Surface 
#0.564     0.781    1.56
# 0.07-2.10 #0.15-1.87 #0.11-7.85
plot(ErCTN$BlubberDEPTH, log(ErCTN$CTN))

fig2b <- ggplot(data=ErCTN, aes(y=CTN, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  theme_bw()+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=12, face= "bold"))+ 
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))++ facet_wrap(~Sex)
fig2b<- fig2b+coord_flip()
fig2b


levels(ErCTN$Keep)

prova2 <-subset(ErCTN, ErCTN$Keep != "N") 
kruskal.test(CTN ~ BlubberDEPTH, data=prova2)

(prova2fit <-aov(log(CTN) ~ BlubberDEPTH, data=ErCTN))
summary(prova2fit)

ggplot(data=prova2, aes(y=CTN, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Corticosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  theme_bw()+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=12, face= "bold"))+ 
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))+coord_flip()

#for corticosterone if we remove samples that only have a surface representative then we lose the significance but not the trend

#Aldosterone
ErALDO <- ErTOT[!is.na(ErTOT$ALDO),]
#remove outlayer 
ErALDO <- subset(ErALDO, ErALDO$ALDO < 14.0)
hist(ErALDO$ALDO) #left skewed
hist(log(ErALDO$ALDO))
range(ErALDO$ALDO) #00.01 -2.35
shapiro.test(ErALDO$ALDO) #p<0.001
shapiro.test(log(ErALDO$ALDO)) # not significant when the outlayer is eliminated
qqnorm(log(ErALDO$ALDO)) #good

bartlett.test(ErALDO$ALDO, ErALDO$Sex) #not homogeneus
summary(ErALDO$Age.class)
bartlett.test(ErALDO$ALDO, ErALDO$Age.class) #not homogeneus variances
summary(ErALDO$BlubberDEPTH)
bartlett.test(ErALDO$ALDO, ErALDO$BlubberDEPTH) #non equal variances
kruskal.test(ALDO ~ BlubberDEPTH, data=ErALDO)
#Kruskal-Wallis chi-squared = 17.147, df = 2, p-value = 0.000189

(fitALDO <-aov(log(ALDO) ~ BlubberDEPTH, data=ErALDO))
summary(fitALDO) #p <0.001
TukeyHSD(fitALDO) #surface and bottom are significantly different 
#surface and bottom are significantly different p=0.02
tapply(ErALDO$ALDO, FUN= mean, INDEX= ErALDO$BlubberDEPTH)
tapply(ErALDO$ALDO, FUN= range, INDEX= ErALDO$BlubberDEPTH)
#Bottom    Medium   Surface 
#0.076    0.076     0.75
# 0.01-0.43 #0.01-0.26 #0.03-2.35
plot(ErALDO$BlubberDEPTH, ErALDO$ALDO)

fig2c <- ggplot(data=ErALDO, aes(y=ALDO, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  theme_bw()+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=12, face= "bold")) + 
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))+ facet_wrap(~Sex)

fig2c= fig2c + coord_flip()
fig2c


prova3 <-subset(ErALDO, ErALDO$Keep != "N") 
qqnorm(log(prova3$ALDO))
(fitALDO_prova <-aov(log(ALDO) ~ BlubberDEPTH, data=prova3))
summary(fitALDO_prova) #p <0.001
TukeyHSD(fitALDO_prova) 


ggplot(data=prova3, aes(y=ALDO, x=BlubberDEPTH))+
  geom_boxplot()+ xlab("")+ ylab("Aldosterone ng/g")+
  geom_jitter(aes(shape=Sex), position=position_jitter(0.2), size=2.5)+
  theme_bw()+
  theme(panel.border= element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black"))+
  theme(axis.text.x= element_text(size=12, face= "bold")) + 
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12, face= "bold")) +
  theme(axis.title.x=element_text(size=10))+ coord_flip()

multiplot(fig2a, fig2b, fig2c, cols=3)

#aldosterone validation GW

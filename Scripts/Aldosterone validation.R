library(ggplot2)
library(plotrix)
library(MASS)
library(plyr)
library(lattice)
library(mixtools)
library(MuMIn)
options(na.action = "na.omit")
library(emmeans)
library(drc)
library(here)
library(car)
library(ggplot2)
library(ggpmisc)
suppressMessages(library(ggbeeswarm))
suppressMessages(library(ggpubr))

library(tidyr)
library(dplyr)

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
#VALIDATION
# 1. ALDOSTERONE
# 1.1 Parallelism

ErobALDO <-read.csv(here("./data/ErobALDO.csv"), header=TRUE)
# GWF_ALDO<- drm(RESPONSE~DOSE, CURVE,
                fct=LL.4 (names=c("Slope", "Lower", "Upper", "ED50")),data=ErobALDO)
# plot(GWF_ALDO)


# DOSE : standard concentration
# RESPONSE : OD, % binding
# CURVE : factor to identify the curves you are comparing

ErobALDO$CURVE= as.factor(ErobALDO$CURVE)

# ErobALDO_clean <- subset(
#   ErobALDO,
#   DOSE > 0 & is.finite(log10(DOSE)) &
#     !is.na(CURVE) & !is.na(RESPONSE)
# )

m2 <- lm(RESPONSE ~ CURVE/log10(DOSE), data = ErobALDO)

summary(m2)  # bottom three rows are the slopes


?linearHypothesis
# STD vs FEMALE
linearHypothesis(m2, "CURVESTD:log10(DOSE) - CURVEFEMALE:log10(DOSE) = 0")
linearHypothesis(m2, "CURVESTD:log10(DOSE) - CURVEMALE:log10(DOSE) = 0")

# Save coefficients so you can extract needed quantities:
cf <- summary(m2)$coef   
cf
str(cf)
# comparing the STD slope to Female  

(d.slope <- cf[6,1] - cf[4,1])    # Difference in slopes
(se.d <- sqrt(cf[6,2]^2 + cf[4,2]^2))    # Standard error of difference

# Note that these are identical to the differences in the model 'm.interactions'
# and in the output from 'lstrends'

# Null hypothesis: d.slope = 0
# Compute t-statistic:
(t.stat <- abs(d.slope / se.d)) 


# To formally evaluate significance we compare 't.stat' to a t-distribution 
# with n-p d.f., where n = 10 observations and p = 4 parameters
# Probability of getting a t-value as large as the observed one:
2*(1 - pt(t.stat, df=6))
# Same p-value as above.
# p value is 0.09


# comparing the STD slope to Male 


(d.slope <- cf[6,1] - cf[5,1])    # Difference in slopes
(se.d <- sqrt(cf[6,2]^2 + cf[5,2]^2))    # Standard error of difference
(t.stat<- abs(d.slope / se.d))
# df= n of observation total (STD + DILUTION) - p (parameters , for each logistic curve so 4 *2)
# df= 5
2*(1 - pt(t.stat, df=5 )) # 0.6


FigS1A <- ggplot(ErobALDO, aes(x = DOSE, y = RESPONSE, color = CURVE)) +
  geom_point(aes(shape=CURVE), size = 2) +
  # geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)+
  # geom_smooth(method = "lm", se = FALSE)+
  geom_line(aes(group = CURVE)) +  # connects observed values
  scale_color_manual(values=c("#999999", "#E69F00", "red"))+
  scale_x_log10(expand = expansion(mult = c(0.1, 0.05))) +
  labs(x = "Standard concentrations, pg/ml (log scale)", y = "Optical density") +
  theme_minimal()+
  annotate("text",
           x = min(ErobALDO$DOSE)/2,
           y = 1,
           label = "A",
           hjust = 1.2,
           size= 6,
           fontface= "bold")+
  coord_cartesian(clip = "off")


Erob_long <- ErobALDO %>%
  pivot_longer(
    cols = c(OBSERVED.F, OBSERVED.M),
    names_to = "GROUP",
    values_to = "OBSERVED",
    values_drop_na = TRUE
  )


FigS1B <- ggplot(Erob_long, aes(x = EXPECTED, y = OBSERVED, color = GROUP)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  
  scale_color_manual(values = c("OBSERVED.F" = "#999999",
                                "OBSERVED.M" = "#E69F00")) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label),
                      after_stat(rr.label),
                      sep = " ~~~ ")),
    formula = y ~ x,
    parse = TRUE)+
  labs(x = "Expected", y = "Observed") +
  annotate("text",x = -Inf, y = Inf, label = "B",hjust = 2,
           vjust = 1.5, size= 6,fontface= "bold")+
  coord_cartesian(clip = "off")+
  theme(plot.margin = margin(10, 10, 10, 20))+  # extra space on left
  theme_minimal()

FigureS1 <-ggarrange(FigS1A, FigS1B, nrow=1, ncol = 2)
ggsave(here("./output/FigureS1.png"),FigureS1, width=12, height=6, dpi=600)
 ?ggarrange
plot(ErobALDO$EXPECTED, ErobALDO$OBSERVED.M, 
     xlab = "Mass added", ylab = "Mass recovered",
     pch = 2, frame.plot = FALSE, xlim=c(0,4000 ), ylim=c(0,4000))
abline(lm(OBSERVED.M ~ EXPECTED, data = ErobALDO, na.action=na.omit), col = "gray23", lty="dashed")
summary(lm(OBSERVED.F ~ EXPECTED, data = ErobALDO, na.action=na.omit))


# check the original data and run dilution linearity

GWM_ALDO<- drm(RESPONSE~DOSE, CURVE,
               fct=LL.4 (names=c("Slope", "Lower", "Upper", "ED50")),data=ErobALDO)
plot(GWM_ALDO)

# # Save coefficients so you can extract needed quantities:
# CF <- summary(GWM_ALDO)$coef   
# (d.slope <- CF[1,1] - CF[2,1])    # Difference in slopes (STD -Pool)
# (se.d <- sqrt(CF[1,2]^2 + CF[2,2]^2))    # Standard error of difference
# (t.stat<- abs(d.slope / se.d))
# # df= n of observation total (STD + DILUTION) - p (parameters , for each logistic curve so 4 *2)
# # df=5 (11-8)
# 2*(1 - pt(t.stat, df=3 )) # 1





# validation graphics


ALDO_ALL <- read.csv(here("./data/ALDO_V_AS.csv"), header=TRUE)
colnames(ALDO_ALL) <- c("STD", "B_BO", "STD2", "GWF", "GWM", "BWF", "BWM")
ggplot(ALDO_ALL, aes(x=log(STD), y=OD)) +
        xlab("Log(relative dose)")+ ylab("Optical density")+
        geom_point(size=2)+ geom_line() +
        geom_point(aes(x=log(STD.1), y=GWF_OD), shape= 15, size=2)+ 
        geom_line(aes(x=log(STD.1), y=GWF_OD ), linetype="twodash")+
        geom_point(aes(x=log(STD2)+4, y=GWM), shape= 17, size=2)+
        geom_line(aes(x=log(STD2)+4, y=GWM ), linetype="longdash")+
        theme_bw()+
        theme(legend.position = "none") + annotate("text", x = 25 , y=18, label = "C", size=7, fontface= "bold")+
        theme(panel.border= element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 0.7, linetype = "solid",
                                       colour = "black"))+
  theme(axis.text.x= element_text(size=14, face= "bold")) + 
        theme(axis.text.y=element_text(size=14, face= "bold")) +
        theme(axis.title.y=element_text(size=14)) +
        theme(axis.title.x=element_text(size=14))






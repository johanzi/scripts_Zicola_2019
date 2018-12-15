

#flowering time new RNAi line culture 34615


#libraries to load
########################################################
library(ggplot2)#to make plots with colored factors
library(multcomp)#Do multicomparison
library(agricolae)#Load TukeyHSD package
########################################################

df <- read.table("figure_3d_data.txt", header=TRUE, sep="\t")

# Reorder factors to display in proper order
df$line <- factor(df$line, levels = c("0T", "ft-10", "2-15-12-T5","1-10","2-14","2-30")) 

#plot individual number
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}
#Plot with ggplot2
ggplot(df, aes(line, rosette)) + 
  geom_boxplot(width=0.8, fill="gray90") + 
  stat_summary(fun.data = give.n, geom = "text") +
  theme_light() +
  expand_limits(y=0) +
  scale_y_continuous(breaks=seq(0,45,10)) +
  ylab("") + 
  xlab("") #+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 


#Normality
qqnorm(df$rosette)
qqline(df$rosette)
shapiro.test(df$rosette)

#ANOVA
attach(df)

fit= aov(rosette~line)
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="greater"))

detach(df)

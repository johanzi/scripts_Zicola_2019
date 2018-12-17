#Flowering time new RNAi line culture 40889

#libraries to load
########################################################
library(multcomp)#Do multicomparison
library(agricolae)#Load TukeyHSD package
########################################################

df <-  read.table("figure_3e_data.txt", stringsAsFactors=TRUE,sep="\t", na.strings = c("NA",""), header=TRUE)

# Reorder factors to display in proper order
df$line <- factor(df$line, levels = c("0T", "ft-10", "C-15","E-16","E-18","E-27","E-29"))

#plot individual number
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}
#Plot with ggplot2
ggplot(df, aes(line, rosette)) + 
  geom_boxplot(width=0.8, fill="gray90") + 
  stat_summary(fun.data = give.n, geom = "text") +
  expand_limits(y=0) +
  scale_y_continuous(breaks=seq(0,60,10)) +
  theme_light() +
  ylab("") + 
  xlab("")# +
  theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 


#Test equality of variances
bartlett.test(rosette~line, data=df)
#If p<5%, we reject hypothesis that variances are equal


#Test normality of the residuals
qqnorm(df$rosette)
qqline(df$rosette)
with(df, shapiro.test(rosette))
#If p<5%, we reject hypothesis that residues follow normal distribution


#ANOVA + Dunnett's test (assuming the 2 previous hypotheses of homoscedasticity and normality of residuals)
fit = aov(rosette~line, data=df)
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="greater")) #We assume a one-sided risk (WT flowers earlier)




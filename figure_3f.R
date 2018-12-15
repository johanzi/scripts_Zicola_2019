#libraries to load


# Samples from cultures 44052, 44082, 44127

library(ggplot2)#to make plots with colored factors
library(multcomp)#Do multicomparison
library(agricolae)

#FT expression for Block C IR construct

df = read.table("figure_3f_data.txt", 
                stringsAsFactors=TRUE,sep = "\t",
                dec = ".",
                na.strings = c("NA",""), 
                header=TRUE)

df$ratio <- as.numeric(df$ratio)

bartlett.test(df, ratio~line)
qqnorm(df$ratio)
qqline(df$ratio)
with(df, shapiro.test(ratio))


# Plot

# Reorder factors to display in proper order
df$line <- factor(df$line, levels = c("0T", "ft-10", "15-2-2","16-4","18-4","27-1","29-11"))

# Plot FT expression
ggplot(df, aes(line, ratio)) +
  stat_summary(fun.y = mean, geom = "bar", position="dodge", fill="gray70") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position="dodge", size=0.5) + 
  geom_jitter(aes(x = line), width=0.1, height=NULL) + 
  theme_light() +
  ylab("") + 
  xlab("")#+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 



###################################################
###################################################

# Parametric test

#ANOVA
attach(df)

fit = aov(ratio~line, df)
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="less"))





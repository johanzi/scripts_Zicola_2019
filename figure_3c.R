#libraries to load

library(ggplot2)#to make plots with colored factors
library(multcomp)#Do multicomparison


#FT expression for Block C IR construct

df = read.table("figure_3c_data.txt", 
                stringsAsFactors=TRUE,sep = "\t",
                dec = ".",
                na.strings = c("NA",""), 
                header=TRUE)

# Replace line30 values by NA
df$ratio[df$line=="line30"] <- NA



# Reorder factors to display in proper order
df$line <- factor(df$line, levels = c("0T", "ft-10", "line15-2","line10","line14","line30"))

# Plot FT expression
ggplot(df, aes(line, ratio)) +
  stat_summary(fun.y = mean, geom = "bar", position="dodge", fill="gray70") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position="dodge", size=0.5) + 
  geom_point(aes(x = line), shape = 19, position = position_dodge(width = 0.9)) + 
  theme_light()  + 
  ylab("") + 
  xlab("")#+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 

# Remove #30 data
df <- subset(df, df$line != "line30")

bartlett.test(df, ratio~line)
qqnorm(df$ratio)
qqline(df$ratio)
with(df, shapiro.test(ratio))

#ANOVA
fit = aov(ratio~line, df)
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="less"))



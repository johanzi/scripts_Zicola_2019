#libraries to load
########################################################
library(ggplot2)#to make plots with colored factors
library(multcomp)#Do multicomparison
library(agricolae)#Load TukeyHSD package
########################################################

# Cultures involved: 34081 and 34073

df = read.table("extended_figure_9b_data.txt", stringsAsFactors=TRUE,sep="\t", na.strings = c("",""), header=TRUE)

col_names = names(df)[1:2]

#Turn the columns in factors
df[,col_names] <- lapply(df[,col_names] , factor)


# Reorder factors to display in proper order
df$line <- factor(df$line, levels = c("0T", "ft-10", "27-4","1-2","21-4","45-1","47-10"))

# Plot data -------------------------

ggplot(df, aes(line, rosette)) + 
  geom_boxplot(width=0.8, fill="gray90") + 
  stat_summary(fun.data = give.n, geom = "text") +
  expand_limits(y=0) +
  scale_y_continuous(breaks=seq(0,60,10)) +
  theme_light() +
  ylab("") + 
  xlab("")+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 


# Statistics ---------------------

boxplot(rosette~line, data = df)
bartlett.test(df, rosette~line)
qqnorm(df$rosette)
qqline(df$rosette)
with(df, shapiro.test(rosette))


# ANOVA with all generations pooled with ft-10
fit = aov(rosette~line, df)
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="greater"))


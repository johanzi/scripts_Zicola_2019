#libraries to load
########################################################
library(ggplot2)#to make plots with colored factors
library(multcomp)#Do multicomparison
library(agricolae)#Load TukeyHSD package
########################################################

#Import data frame with all data

#The data sets contains different columns with different factors and 2 columns for response variable which are number of rosette laves and number of cauline leaves
df = read.table("extended_figure_1_data.txt", 
                stringsAsFactors=TRUE,sep = "\t", 
                na.strings = c("NA",""), 
                header=TRUE)


#Create list with the desired column names
col_names = names(df)[1:3]

#Turn the columns in factors
df[,col_names] <- lapply(df[,col_names] , factor)


give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}


ggplot(df, aes(line, rosette, fill=generation)) + 
  geom_boxplot() + stat_summary(fun.data = give.n, geom = "text") +
  theme_light() +
  expand_limits(y=0) +
  ylab("") + 
  xlab("") +
  theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 

# Statistics---------------

#T3
T3_line = subset(df, generation == "T3" | line == "0T")
bartlett.test(T3_line$rosette~T3_line$line, data = T3_line)
qqnorm(T3_line$rosette)
qqline(T3_line$rosette)

fit_df = aov(rosette~line, T3_line)
summary(glht(fit_df, linfct=mcp(line="Dunnett"), alternative="greater"))


#T4
T4_line = subset(df, generation == "T4" | line == "0T")
bartlett.test(T4_line$rosette~T4_line$line, data = T4_line)
qqnorm(T4_line$rosette)
qqline(T4_line$rosette)

fit_df = aov(rosette~line, T4_line)
summary(glht(fit_df, linfct=mcp(line="Dunnett"), alternative="greater"))


#T5
T5_line = subset(df, generation == "T5" | line == "0T")
bartlett.test(T5_line$rosette~T5_line$line, data = T5_line)
qqnorm(T5_line$rosette)
qqline(T5_line$rosette)

fit_df = aov(rosette~line, T5_line)
summary(glht(fit_df, linfct=mcp(line="Dunnett"), alternative="greater"))






#Flowering time Block C IR line culture 31660

#Grown in MD conditions (12h light / 12h dark)

#libraries to load
########################################################
library(ggplot2)#to make plots with colored factors
library(multcomp)#Do multicomparison
library(agricolae)#Load TukeyHSD package
########################################################


df = read.table("figure_1c_data.txt", stringsAsFactors=TRUE,sep="\t", na.strings = c("NA",""), header=TRUE)

col_names = names(df)[1:2]

#Turn the columns in factors
df[,col_names] <- lapply(df[,col_names] , factor)


# Plot data-------------------

# Reorder to have the lines matching to the picture in fig. 1d
df$line <- factor(df$line, levels =  c("0T", "ft.10", "2.15", "3.15", "3.27", "4.27"))


# Plot all lines and separate by generation
ggplot(data = df, aes(line, rosette, fill=generation)) +   geom_boxplot(lwd=0.2, outlier.size = 0.5) + theme(axis.text.y=element_text(size=rel(1), angle=0, colour="black")) + ylab("") + xlab("")+ theme_bw() + theme(legend.text=element_text(size=8), legend.title=element_text(size=8))


# Statistics-----------------

#df
boxplot(rosette~line, data = df)
bartlett.test(df, rosette~line)
qqnorm(df$rosette)
qqline(df$rosette)
with(df, shapiro.test(rosette))

#Variance test For each generation
bartlett.test(subset(df, generation == "T3" | line == "0T"), rosette~line)
bartlett.test(subset(df, generation == "T4" | line == "0T"), rosette~line)
bartlett.test(subset(df, generation == "T5" | line == "0T"), rosette~line)
bartlett.test(subset(df, generation == "T6" | line == "0T"), rosette~line)
#In no case the Ho is not rejected, I do not have homoscedasticity

#Normality test for each generation
with(subset(df, generation == "T3" | line == "0T"), shapiro.test(rosette))
with(subset(df, generation == "T4" | line == "0T"), shapiro.test(rosette))
with(subset(df, generation == "T5" | line == "0T"), shapiro.test(rosette))
with(subset(df, generation == "T6" | line == "0T"), shapiro.test(rosette))
#Only for T3 I do not respect normality of residuals (probably because of 4#27)


# ANOVA with all generations pooled w/o ft-10
fit = aov(rosette~line, subset(df, line != "ft.10"))
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="greater"))

# ANOVA with all generations pooled with ft-10
fit = aov(rosette~line, df)
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="greater"))


#ANOVA + Dunnett's test Per generation
fitT3 = aov(rosette~line, subset(df, generation == "T3" | line == "0T"))
summary(glht(fitT3, linfct=mcp(line="Dunnett"), alternative="greater"))

fitT4 = aov(rosette~line, subset(df, generation == "T4" | line == "0T"))
summary(glht(fitT4, linfct=mcp(line="Dunnett"), alternative="greater"))

fitT5 = aov(rosette~line, subset(df, generation == "T5" | line == "0T"))
summary(glht(fitT5, linfct=mcp(line="Dunnett"), alternative="greater"))

fitT6 = aov(rosette~line, subset(df, generation == "T6" | line == "0T"))
summary(glht(fitT6, linfct=mcp(line="Dunnett"), alternative="greater"))



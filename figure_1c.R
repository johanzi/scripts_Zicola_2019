#Flowering time new RNAi line culture 31660


#libraries to load
########################################################
library(ggplot2)#to make plots with colored factors
library(multcomp)#Do multicomparison
library(agricolae)#Load TukeyHSD package
########################################################


df31660 = read.table("figure_1c_data.txt", stringsAsFactors=TRUE,sep="\t", na.strings = c("NA",""), header=TRUE)

col_names = names(df31660)[1:7]

#Turn the columns in factors
df31660[,col_names] <- lapply(df31660[,col_names] , factor)


# Plot data-------------------

# Reorder to have the lines matching to the picture in fig. 1d
df31660$line <- factor(df31660$line, levels =  c("0T", "ft.10", "2.15", "3.15", "3.27", "4.27"))


# Plot all lines and separate by generation
ggplot(data = df31660, aes(line, rosette, fill=generation)) +   geom_boxplot(lwd=0.2, outlier.size = 0.5) + theme(axis.text.y=element_text(size=rel(1), angle=0, colour="black")) + ylab("") + xlab("")+ theme_bw() + theme(legend.text=element_text(size=8), legend.title=element_text(size=8))


# Statistics-----------------

#df31660
boxplot(rosette~line, data = df31660)
bartlett.test(df31660, rosette~line)
qqnorm(df31660$rosette)
qqline(df31660$rosette)
with(df31660, shapiro.test(rosette))

#Variance test For each generation
bartlett.test(subset(df31660, generation == "T3" | line == "0T"), rosette~line)
bartlett.test(subset(df31660, generation == "T4" | line == "0T"), rosette~line)
bartlett.test(subset(df31660, generation == "T5" | line == "0T"), rosette~line)
bartlett.test(subset(df31660, generation == "T6" | line == "0T"), rosette~line)
#In no case the Ho is not rejected, I do not have homoscedasticity

#Normality test for each generation
with(subset(df31660, generation == "T3" | line == "0T"), shapiro.test(rosette))
with(subset(df31660, generation == "T4" | line == "0T"), shapiro.test(rosette))
with(subset(df31660, generation == "T5" | line == "0T"), shapiro.test(rosette))
with(subset(df31660, generation == "T6" | line == "0T"), shapiro.test(rosette))
#Only for T3 I do not respect normality of residuals (probably because of 4#27)


# ANOVA with all generations pooled w/o ft-10
fit = aov(rosette~line, subset(df31660, line != "ft.10"))
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="greater"))

# ANOVA with all generations pooled with ft-10
fit = aov(rosette~line, df31660)
summary(glht(fit, linfct=mcp(line="Dunnett"), alternative="greater"))


#ANOVA + Dunnett's test Per generation
fitT3 = aov(rosette~line, subset(df31660, generation == "T3" | line == "0T"))
summary(glht(fitT3, linfct=mcp(line="Dunnett"), alternative="greater"))

fitT4 = aov(rosette~line, subset(df31660, generation == "T4" | line == "0T"))
summary(glht(fitT4, linfct=mcp(line="Dunnett"), alternative="greater"))

fitT5 = aov(rosette~line, subset(df31660, generation == "T5" | line == "0T"))
summary(glht(fitT5, linfct=mcp(line="Dunnett"), alternative="greater"))

fitT6 = aov(rosette~line, subset(df31660, generation == "T6" | line == "0T"))
summary(glht(fitT6, linfct=mcp(line="Dunnett"), alternative="greater"))



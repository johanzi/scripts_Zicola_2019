#libraries to load

library(ggplot2)#to make plots with colored factors
library(multcomp)#Do multicomparison
library(agricolae)

#FT expression for Block C IR construct

df = read.table("figure_1e_data.txt", 
                stringsAsFactors=TRUE,sep = "\t",
                dec = ".",
                na.strings = c("NA",""), 
                header=TRUE)

# Plot FT expression
ggplot(df, aes(line, ratio, fill = generation)) +
                stat_summary(fun.y = mean, geom = "bar", position="dodge") + 
                stat_summary(fun.data = mean_se, geom = "errorbar", position="dodge", size=0.5) + 
                geom_point(aes(x = line), shape = 20, position = position_dodge(width = 0.9)) + 
                theme_light()  + 
                scale_y_continuous(breaks = round(seq(min(df$ratio), max(df$ratio), by = 0.5),1))+
                ylab("") + 
                xlab("")# +
                theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 



bartlett.test(df, ratio~line)
qqnorm(df$ratio)
qqline(df$ratio)
with(df, shapiro.test(ratio))


#T3
fitT3 = aov(ratio~line, subset(df, generation=="T3"))
DunnettT3 <- glht(fitT3, linfct=mcp(line="Dunnett"), alternative="less")
summary(DunnettT3)

#T5
# Remove outsider value 3.27 T5
df_sub <- subset(df, df$line != "3.27")
fitT5 = aov(ratio~line, subset(df_sub, generation=="T5"))
DunnettT5 <- glht(fitT5, linfct=mcp(line="Dunnett"), alternative="less")
summary(DunnettT5)


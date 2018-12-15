
library(ggplot2)
library(agricolae)
library(multcomp)


# Analysis 48031 ----------------------------------------------------------

df <- read.table("figure_3g_data.txt", header=TRUE)

df$total <- df$rosette + df$cauline


# Plot data ---------------------------------------------------------------
df$cross <- factor(df$cross, levels = c("0T", "ft-10", "Block_C_27","Block_C_15","Block_E_18","Block_E_16","Col-0xBlock_C_27","Col-0xBlock_C_15","Block_C_27xCol-0","Block_C_15xCol-0","Block_E_16xCol-0","Col-0xBlock_E_18","Block_E_18xCol-0","Block_E_16xBlock_E_27","Block_E_16xBlock_E_15","Block_E_18xBlock_E_27","Block_E_18xBlock_E_15","Block_C_27xBlock_E_18","Block_C_27xBlock_E_16","Block_C_15xBlock_E_18","Block_C_15xBlock_E_16")) 

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

#Plot with ggplot2
ggplot(df, aes(cross, rosette)) + 
  geom_boxplot(width=0.7, fill="gray90") + 
  stat_summary(fun.data = give.n, geom = "text" ) +
  theme_light() +
  expand_limits(y=0) +
  scale_y_continuous(breaks=seq(0,45,10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ylab("") + 
  xlab("") #+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 


# Normality and homoscedasticity ------------------------------------------

# Rosette leaf count
bartlett.test(df, rosette~construct)
qqnorm(df$rosette)
qqline(df$rosette)
with(df, shapiro.test(rosette))

# total leaf count
bartlett.test(df, total~construct)
qqnorm(df$total)
qqline(df$total)
with(df, shapiro.test(total))


# Kruska Wallis test -------------------------------------------------------------------

# Kruskal-Wallis multiple comparison with agricolae
kru <- kruskal(df$rosette, df$construct, alpha=0.05, p.adj=c("bonferroni"), group=TRUE)
kru


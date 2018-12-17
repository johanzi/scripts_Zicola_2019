#analysis Block C ChIP data

#Data analysis H3K9me2

#libraries to load
########################################################
library(ggplot2)#to make plots with colored factors
########################################################

#Load data in a dataframe
data = read.table("figure_2e_data.txt", header=TRUE, dec=",", sep="\t")

# Plot data
ggplot(data, aes(line, value, fill=rep)) +
  stat_summary(fun.y = mean, geom = "bar", position="dodge") + 
  theme_light()  + 
  ylab("") + 
  xlab("") #+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=8), legend.title=element_text(size=8)) 



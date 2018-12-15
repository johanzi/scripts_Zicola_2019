library(ggplot2)

df <- read.table("figure_4b_data.txt", header=TRUE)

# Plot data
ggplot(data=df, aes(y=value, x=construct)) +  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_bw() +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, hjust=1)) 


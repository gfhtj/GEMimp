#plot
library(ggplot2)

# create a dataframe
data <- data.frame(
  Method = rep(c("Maaslin2", "corncob", "deseq"), each = 6),
  Indicator = rep(c("Precision", "Recall", "F1-score"), times = 6),
  Category = rep(c("Normal DA", "GE-impute DA"), each = 3, times = 3),
  Value = c(0, 0.06667, 0.04444, NaN, 0.5, 0.33333,
            0.02222, 0.08889, 0.04444, 0.5, 0.44444, 0.66667,
            NaN, 0.11765, 0.07842, 0.04255, 0.14815, 0.08333)
)

indicator_order <- c("Precision", "Recall", "F1-score")
data$Indicator <- factor(data$Indicator, levels = indicator_order)

color_scale <- c( "black","gray")


p <- ggplot(data, aes(x = Indicator, y = Value, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_grid(cols = vars(Method), scales = "free_x", switch = "x") +  
  scale_fill_manual(values = color_scale) +  
  theme_minimal() +
  theme(legend.position = "bottom",  
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(angle = 0),
        strip.background = element_blank(),  
        axis.title = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank(),  
        axis.line = element_line(color = "black")) +  
  labs(fill = "Category") +  
  scale_x_discrete(limits = indicator_order)  

print(p)


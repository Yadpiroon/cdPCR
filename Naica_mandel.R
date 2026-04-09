# Naica analysis

######ILS analysis####
library(readxl)
library(ILS)

dataframe <- read_excel("dataframe.xlsx", 
                        col_types = c("numeric", "text", "text", 
                                      "text"))
View(dataframe)

qcdata <- lab.qcdata(dataframe)
summary(qcdata)
head(qcdata)

# Fig 1
plot(qcdata, ylab = "Methods", xlab = "ARGs concentration in river water samples")

# Fig 2
qcdata <- lab.qcdata(
  data = dataframe,
  var.index = 1,         
  replicate.index = 2,   
  material.index = 3,    
  laboratory.index = 4,  
  data.name = "dataframe"  
)
head(qcdata)

library(dplyr)
variability <- qcdata %>%
  group_by(material) %>%
  summarize(
    S_r = sd(x),               # Within-method variability
    S_B = mean(x),             # Between-method variability
    S_R = sqrt(S_r^2 + S_B^2)  # Total variability
  )
head(variability)

library(tidyr)
variability_long <- variability %>%
  pivot_longer(cols = c(S_r, S_B, S_R), names_to = "Measure", values_to = "Deviation")
View(variability_long)
arg_names <- c("ampC", "mcr", "mecA", "NDM", "vanA")

# Plot with ARG names on the y-axis
plot(
  variability_long$Deviation,
  as.numeric(as.factor(variability_long$material)),
  pch = 16,
  xlab = "Deviations",
  ylab = "ARGs",
  main = "Measures of Variability",
  col = as.numeric(as.factor(variability_long$Measure)),
  yaxt = "n" # Suppress default y-axis labels
)

# Add custom y-axis labels
axis(2, at = 1:5, labels = arg_names)

# Add legend
legend(
  "topright",
  legend = unique(variability_long$Measure),
  col = 1:length(unique(variability_long$Measure)),
  pch = 16
)

# Cochran's Test
cochran_results <- cochran.test(qcdata)
print(cochran_results) #Critical value: 0.4595731

# k statistic 
k <- k.qcs(qcdata, alpha = 0.005)
plot(k)
summary(k) #Critical value:  1.208986
print(k)
head(k)

# h statistics 
h <- h.qcs(qcdata, alpha = 0.005)
summary(h) #Critical value:  1.154665
plot(h)
head(h)

h_stats <- h.qcs(qcdata, alpha = 0.005)
plot(h_stats)  # Visualize results
summary(h_stats)

# Grubbs' Test
grubbs_results <- grubbs.test(qcdata)
print(grubbs_results) #Critical value: 2.086079

aov <- ils.aov(qcdata)

plot(h_stats, ylab = "Methods", xlab = "ARGs concentration in river samples")

library(readxl)
library(ILS)

dataframe2 <- read_excel("dataframe2.xlsx", 
                        col_types = c("numeric", "text", "text", 
                                      "text"))
View(dataframe2)

qcdata2 <- lab.qcdata(dataframe2)
summary(qcdata2)
head(qcdata2)

# Fig 1
plot(qcdata2, ylab = "Methods", xlab = "ARGs concentration in river water samples")

# Fig 2
qcdata66 <- lab.qcdata(
  data = dataframe2,
  var.index = 1,         
  replicate.index = 2,   
  material.index = 3,    
  laboratory.index = 4,  
  data.name = "dataframe"  
)

head(qcdata66)
View(qcdata66)

library(dplyr)
variability2 <- qcdata2 %>%
  group_by(material) %>%
  summarize(
    S_r = sd(x),               # Within-method variability
    S_B = mean(x),             # Between-method variability
    S_R = sqrt(S_r^2 + S_B^2)  # Total variability
  )
head(variability2)

library(tidyr)
variability_long2 <- variability2 %>%
  pivot_longer(cols = c(S_r, S_B, S_R), names_to = "Measure", values_to = "Deviation")

arg_names <- c("ampC", "mcr-1", "mecA", "NDM", "TEM", "vanA")

# Plot with ARG names on the y-axis
plot(
  variability_long2$Deviation,
  as.numeric(as.factor(variability_long2$material)),
  pch = 16,
  xlab = "Deviations",
  ylab = "ARGs",
  main = "Measures of Variability",
  col = as.numeric(as.factor(variability_long2$Measure)),
  yaxt = "n" # Suppress default y-axis labels
)

axis(2, at = 1:1, labels = arg_names)

# Add legend
legend(
  "topright",
  legend = unique(variability_long2$Measure),
  col = 1:length(unique(variability_long2$Measure)),
  pch = 16
)


# Cochran's Test
cochran_results2 <- cochran.test(qcdata2)
print(cochran_results2) #Critical value: 0.4595731

# k statistic 
k2 <- k.qcs(qcdata2, alpha = 0.005)
plot(k2)
summary(k2) #Critical value:  1.208986
print(k2)
head(k2)

# h statistics 
h2 <- h.qcs(qcdata2, alpha = 0.005)
summary(h2) #Critical value:  1.154665
plot(h2)
head(h2)

h_stats2 <- h.qcs(qcdata2, alpha = 0.005)
plot(h_stats2)  # Visualize results
summary(h_stats2)

# Grubbs' Test
grubbs_results2 <- grubbs.test(qcdata2)
print(grubbs_results2) #Critical value: 2.086079

# gene compare

# library
getwd()
library(ggplot2)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)

library(readr)
threemethodx <- read_csv("C:/Users/Yadpiroon Siri/OneDrive - University of Yamanashi/PhD_1drive/PhD year 1/Naica6_YS21jun24/Analysis_Naica/threemethodx.csv", 
                         col_types = cols(type = col_character(), 
                                          site = col_character(), qPCR = col_number(), 
                                          cdPCR = col_number(), 
                                          gene = col_character()))
head(threemethodx)

library(readr)
allextem <- read_csv("allextem.csv", col_types = cols(type = col_character(), 
                                                      site = col_character(), x = col_number(), 
                                                      gene = col_character(), method = col_character(), 
                                                      qdc = col_character()))
head(allextem)

# Method 
one.way.anova_event <- aov(allextem$x ~ allextem$qdc, data = allextem)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)

###################qPCR##################
# Site 
one.way.anova_event <- aov(allextem$site ~ allextem$x, data = allextem)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)

# Type 
one.way.anova_event <- aov(threemethodx$qPCR ~ threemethodx$type, data = threemethodx)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)


# Gene 
one.way.anova_event <- aov(threemethodx$qPCR ~ threemethodx$gene, data = threemethodx)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)

###################dPCR##################
# Site 

up <- allextem[1:414, ]
dw <- allextem[415:828, ]

summary(up$x)
sd(up$x)
shapiro.test(up$x)

summary(dw$x)
sd(dw$x)
shapiro.test(dw$x)

t.test(up$x, dw$x)

# Type 
one.way.anova_event <- aov(threemethodx$dPCR ~ threemethodx$type, data = threemethodx)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)

# Gene 
one.way.anova_event <- aov(threemethodx$dPCR ~ threemethodx$gene, data = threemethodx)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)

###################cdPCR##################
# Site 
one.way.anova_event <- aov(threemethodx$cdPCR ~ threemethodx$site, data = threemethodx)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)

# Type 
one.way.anova_event <- aov(threemethodx$cdPCR ~ threemethodx$type, data = threemethodx)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)

# Gene 
one.way.anova_event <- aov(threemethodx$cdPCR ~ threemethodx$gene, data = threemethodx)
summary(one.way.anova_event)
TukeyHSD(one.way.anova_event)


#############r #################
# new corr plot
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(reshape2)

library(readr)
corrrrrr <- read_csv("corrrrrr.csv", col_types = cols(qPCR_mecA = col_number(), 
                                                      qPCR_vanA = col_number(), qPCR_TEM = col_number(), 
                                                      qPCR_NDM = col_number(), qPCR_mcr1 = col_number(), 
                                                      qPCR_ampC = col_number(), dPCR_NDM = col_number(), 
                                                      dPCR_ampC = col_number(), dPCR_mcr1 = col_number(), 
                                                      dPCR_mecA = col_number(), dPCR_vanA = col_number(), 
                                                      cdPCR_mecA = col_number(), cdPCR_vanA = col_number(), 
                                                      cdPCR_TEM = col_number(), cdPCR_NDM = col_number(), 
                                                      cdPCR_mcr1 = col_number(), cdPCR_ampC = col_number()))
View(corrrrrr)

corr <- round(cor(corrrrrr, method = "spearman"), 2)

p.df <- as.data.frame(ggcorrplot::cor_pmat(corrrrrr))

labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")
}

p.labs = p.df %>%
  mutate_all(labs.function)

p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

cor_plot = ggcorrplot(corr, hc.order = F, type = "lower",
                      lab = T, ggtheme = ggplot2::theme_gray, colors = c("#3EDBF0", "white", "#FF75A0")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 10), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 2),
        panel.background = element_rect(fill = "#F1F1F1", colour = NA),
        panel.grid.minor = element_line(colour = "white", size = 0.2),
        panel.grid.major = element_line(colour = "white", size = 0.2))

p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                  paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot.labs = cor_plot +
  geom_text(aes(x = p.labs$Var1,
                y = p.labs$Var2),
            label = p.labs$lab,
            nudge_y = 0.25,
            size = 5)

cor.plot.labs


# cut heg gene meca vana mcr1

library(readr)
corr3 <- read_csv("corr3.csv", col_types = cols(qPCR_tem = col_number(), 
                                                qPCR_ampC = col_number(), qPCR_ndm = col_number(), 
                                                cdPCR_tem = col_number(), cdPCR_ampC = col_number(), 
                                                cdPCR_ndm = col_number()))
head(corr3)

corr <- round(cor(corr3, method = "spearman"), 2)

p.df <- as.data.frame(ggcorrplot::cor_pmat(corr3))

labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 ~ "*")
}

p.labs = p.df %>%
  mutate_all(labs.function)

p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

cor_plot = ggcorrplot(corr, hc.order = F, type = "lower",
                      lab = T, ggtheme = ggplot2::theme_gray, colors = c("#3EDBF0", "white", "#FF75A0")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 10), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 2),
        panel.background = element_rect(fill = "#F1F1F1", colour = NA),
        panel.grid.minor = element_line(colour = "white", size = 0.2),
        panel.grid.major = element_line(colour = "white", size = 0.2))

p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                  paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot.labs = cor_plot +
  geom_text(aes(x = p.labs$Var1,
                y = p.labs$Var2),
            label = p.labs$lab,
            nudge_y = 0.25,
            size = 5)

cor.plot.labs

# qPCR and cd tem pos corr

# cut ampC

library(readr)
corr3 <- read_csv("corr3.csv", col_types = cols(qPCR_tem = col_number(), 
                                                qPCR_ampC = col_number(), qPCR_ndm = col_number(), 
                                                cdPCR_tem = col_number(), cdPCR_ampC = col_number(), 
                                                cdPCR_ndm = col_number()))
head(corr3)

new3 <- corr3 %>% select(qPCR_tem, qPCR_ndm, cdPCR_tem, cdPCR_ndm)
head(new3)

corr <- round(cor(new3, method = "spearman"), 2)

p.df <- as.data.frame(ggcorrplot::cor_pmat(new3))

# Rename columns and rows for proper subscripts in labels
colnames(corr) <- gsub("ndm", "bla[NDM]", colnames(corr))
colnames(corr) <- gsub("tem", "bla[TEM]", colnames(corr))
rownames(corr) <- colnames(corr)  # Ensure row names match

# Rename columns in p-value dataframe for consistency
colnames(p.df) <- gsub("ndm", "bla[NDM]", colnames(p.df))
colnames(p.df) <- gsub("tem", "bla[TEM]", colnames(p.df))
rownames(p.df) <- colnames(p.df)  # Ensure row names match


labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 ~ "*")
}

p.labs = p.df %>%
  mutate_all(labs.function)

p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

cor_plot = ggcorrplot(corr, hc.order = F, type = "lower",
                      lab = T, ggtheme = ggplot2::theme_gray, colors = c("#3EDBF0", "white", "#FF75A0")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 10), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 2),
        panel.background = element_rect(fill = "#F1F1F1", colour = NA),
        panel.grid.minor = element_line(colour = "white", size = 0.2),
        panel.grid.major = element_line(colour = "white", size = 0.2))

p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                  paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot.labs = cor_plot +
  geom_text(aes(x = p.labs$Var1,
                y = p.labs$Var2),
            label = p.labs$lab,
            nudge_y = 0.25,
            size = 5)

cor.plot.labs

# only tem in corr graph


ggsave(file="cor.plot.labsa.jpeg", cor.plot.labs,
       width= 90, height = 70, units = "mm", dpi=600)

####amc####
library(readr)
library(ggplot2)
library(ggcorrplot)
library(reshape2)
library(dplyr)

# Load Data
corr3 <- read_csv("corr3.csv", col_types = cols(qPCR_tem = col_number(), 
                                                qPCR_ampC = col_number(), qPCR_ndm = col_number(), 
                                                cdPCR_tem = col_number(), cdPCR_ampC = col_number(), 
                                                cdPCR_ndm = col_number()))

# Select Relevant Columns
new3 <- corr3 %>% select(qPCR_tem, qPCR_ndm, cdPCR_tem, cdPCR_ndm)

# Compute Correlation Matrix
corr <- round(cor(new3, method = "spearman"), 2)

# Compute p-values
p.df <- as.data.frame(ggcorrplot::cor_pmat(new3))

# Assign Unique and Properly Formatted Labels
new_colnames <- c(
  "bla[TEM]~(qPCR)",  # qPCR_tem
  "bla[NDM]~(qPCR)",  # qPCR_ndm
  "bla[TEM]~(cdPCR)", # cdPCR_tem
  "bla[NDM]~(cdPCR)"  # cdPCR_ndm
)

# Apply to Correlation Matrix and p-value Matrix
colnames(corr) <- rownames(corr) <- parse(text = new_colnames)
colnames(p.df) <- parse(text = new_colnames)
rownames(p.df) <- parse(text = new_colnames)

# Function to Label Significance
labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 ~ "*")
}

p.labs = p.df %>%
  mutate_all(labs.function)

p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

# Generate Correlation Plot with Custom Labels
cor_plot = ggcorrplot(corr, hc.order = F, type = "lower",
                      lab = T, ggtheme = ggplot2::theme_gray(), colors = c("#3EDBF0", "white", "#FF75A0")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 10),
        legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", family = "Arial", size = 10),
        axis.text.y = element_text(color = "black", family = "Arial", size = 10),
        panel.background = element_rect(fill = "#F1F1F1", colour = NA),
        panel.grid.minor = element_line(colour = "white", size = 0.2),
        panel.grid.major = element_line(colour = "white", size = 0.2))

# Add Significance Labels
p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                  paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot.labs = cor_plot +
  geom_text(aes(x = p.labs$Var1,
                y = p.labs$Var2),
            label = p.labs$lab,
            nudge_y = 0.25,
            size = 5)

# Display the Plot
cor.plot.labs


####analysis####


library(readr)
detection <- read_csv("detection.csv", col_types = cols(qPCR_positive = col_number(), 
                                                        qt = col_number(), cdPCR_positive = col_number(), 
                                                        ct = col_number()))
head(detection)

library(dplyr)


detection %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(
      matrix(
        c(qPCR_positive, qt - qPCR_positive,
          cdPCR_positive, ct - cdPCR_positive),
        nrow = 2
      )
    )$p.value
  ) -> fisher_results

print(fisher_results)


library(readr)
concentrationKW <- read_csv("concentrationKW.csv", 
                            col_types = cols(qPCR_Concentration = col_number(), 
                                             cdPCR_Concentration = col_number()))
View(concentrationKW)

# Wilcoxon Signed-Rank Test for Concentration Data
wilcoxon_results <- concentrationKW %>%
  group_by(ARG) %>%
  summarise(
    p_value = wilcox.test(qPCR_Concentration, cdPCR_Concentration, paired = TRUE)$p.value
  )

print(wilcoxon_results)

######Graph qPCR and cdPCR####

getwd()
library(tidyverse)
library(gapminder)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(corrplot)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)
library(ggsignif)
library(readxl)
library(readxl)
library(patchwork)
library(cowplot)
library(lattice)
library(PASWR)
library(factoextra)


####CP round 1 *4 ####

Round1 <- read_excel("A2.xlsx", col_types = c("text", 
                                              "text", "numeric", "numeric", "numeric", 
                                              "numeric", "numeric", "numeric"))

View(Round1)

A <- Round1[1:20, ]
head(A)
B <- Round1[21:26, ]
head(B)

xmin <- 0
xmax <- 6
ymin <- 0
ymax <- 6

p1 <- ggplot(Round1, aes(x = temq, y = temc)) +
  geom_point(alpha = 1, size = 1.5, colour = "black", fill = "#00AFBB", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  
  ggtitle(expression(italic(bla)[TEM])) +
  
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p1
save_plot("p1.jpeg", p1)

p2 <- ggplot(Round1, aes(x = ndmq, y = ndmc)) +
  geom_point(alpha = 1, size = 1.5, colour = "black", fill = "#E7B800", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  ggtitle(expression(italic(bla)[NDM])) +
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p2
save_plot("p2.jpeg", p2)

p3 <- ggplot(Round1, aes(x = ampcq, y = ampcc)) +
  geom_point(alpha = 1, size = 1.5, colour = "black", fill = "#FC4E07", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +

  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  ggtitle(expression(italic(AmpC))) +
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p3
save_plot("p3.jpeg", p3)



# Combine the plots into a single figure
fig1 <- plot_grid(p1, p2, p3, ncol = 1,
                           labels = c("a", "b", "c"), label_size = 11)

# Save the combined figure
ggsave(file="fig1.jpeg", fig1, width= 90, height = 210, units = "mm", dpi=600)


#### tip/mem only ####
Round1 <- read_excel("A2.xlsx", col_types = c("text", 
                                              "text", "numeric", "numeric", "numeric", 
                                              "numeric", "numeric", "numeric"))

head(Round1)

A <- Round1[1:20, ]
head(A)
B <- Round1[21:26, ]
head(B)

xmin <- 0
xmax <- 6
ymin <- 0
ymax <- 6

p1 <- ggplot(A, aes(x = temq, y = temc)) +
  geom_point(alpha = 1, size = 2, colour = "black", fill = "#00AFBB", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  
  ggtitle(expression("CP Select -" ~ italic(bla)[TEM])) +
  
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p1
save_plot("p1.jpeg", p1)

p2 <- ggplot(A, aes(x = ndmq, y = ndmc)) +
  geom_point(alpha = 1, size = 2, colour = "black", fill = "#E7B800", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  ggtitle(expression("CP Select -" ~ italic(bla)[NDM])) +
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p2
save_plot("p2.jpeg", p2)

p3 <- ggplot(A, aes(x = ampcq, y = ampcc)) +
  geom_point(alpha = 1, size = 2, colour = "black", fill = "#FC4E07", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  ggtitle(expression("CP Select -" ~ italic(AmpC))) +
  
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p3
save_plot("p3.jpeg", p3)

p4 <- ggplot(B, aes(x = temq, y = temc)) +
  geom_point(alpha = 1, size = 2, colour = "black", fill = "#00AFBB", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  
  ggtitle(expression("EMF -" ~ italic(bla)[TEM])) +
  
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p4
save_plot("p4.jpeg", p4)

p5 <- ggplot(B, aes(x = ndmq, y = ndmc)) +
  geom_point(alpha = 1, size = 2, colour = "black", fill = "#E7B800", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  ggtitle(expression("EMF -" ~ italic(bla)[NDM])) +
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p5
save_plot("p5.jpeg", p5)

p6 <- ggplot(B, aes(x = ampcq, y = ampcc)) +
  geom_point(alpha = 1, size = 2, colour = "black", fill = "#FC4E07", shape = 21, stroke = 0.5) +
  
  geom_segment(aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
               linetype = "dashed", color = "grey") +
  
  scale_x_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  scale_y_continuous(limits = c(0, 6),
                     labels = function(x) format(round(x, 1), nsmall = 1)) +
  
  xlab(expression("qPCR " * (Log[10] * " copies/μL"))) + 
  ylab(expression("cdPCR " * (Log[10] * " copies/μL"))) +
  ggtitle(expression("EMF -" ~ italic(AmpC))) +
  
  theme(axis.title.x = element_text(color = "black", size = 9, face = "bold"), # X-axis title
        axis.title.y = element_text(color = "black", size = 9, face = "bold"), # Y-axis title
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p6
save_plot("p6.jpeg", p6)

# Combine the plots into a single figure
fig1 <- plot_grid(p1, p4, p2, p5, p3, p6, ncol = 2,
                  labels = c("a", "b", "c", "d", "e", "f"), label_size = 11)

# Save the combined figure
ggsave(file="fig1.jpeg", fig1, width= 170, height = 230, units = "mm", dpi=600)

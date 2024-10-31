library(reshape2)
library(readr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(bbmle)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggpubr)
library(lubridate)
library(circular)
library(ggpubr)
library(lme4)
library(lmerTest)

period_2024 <- read_csv("./data/period_2024.csv")[,-1]
period_2022 <- read_csv("./data/period_2022.csv")[,-1]

period_all <- rbind(period_2022, period_2024)
unique(period_all$family)

period_all$family = as.factor(period_all$family) #group by father's family line
levels(period_all$family) = list(
  A = 'ah18',
  B = 'dh63',
  C = 'a17_h',
  D = 'ch27',
  E = c('CD', 'CH', 'CI', 'CJ', 'CK'),
  F = c('CE', 'CF', 'CG', 'CO')
)

period_all$period = as.factor(period_all$period)
period_all$period <- relevel(period_all$period, ref = "mutant")
period_all$ld <- as.factor(period_all$ld)
period_all$ld <- relevel(period_all$ld, ref = '12L12D')
period_all$family <- as.factor(period_all$family)
period_all$azt <- floor(period_all$azt)
period_all$azt <- circular(period_all$azt, units = 'hours', template='clock24') 


figure3b_l = ggplot(data = period_all %>% filter(ld == '12L12D'),
                    aes(x = pdd, y = family, colour = period, alpha = period)) + 
  geom_boxplot(outlier.shape = NA, size=0.25, alpha = 1) +
  geom_point(aes(x = pdd, y = family),
             position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0.5),
             size = 2, shape = 16) +
  scale_alpha_manual(values = c(0.5, 0.5)) +
  scale_colour_manual(values = c('mutant' = "#00CD00", 'wt' = "darkorchid"),
                      labels = c('wt' = expression(italic(per)^"+/+"),
                                 'mutant' = expression(italic(per)^"-/+")),
                      breaks = c('wt', 'mutant'))+
  scale_x_continuous(name = "PDD Time (days)",
                     limits = c(0,70), breaks = seq(0,70,by=10)) +
  scale_y_discrete(limits = rev, name = "Family") +
  theme_bw()  +
  theme(text = element_text(size = 14, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(alpha = 'none')

#Square root of weights to stabilize variance across period genotypes and family, so that no single family biases the overall mean by contributing more. 
#unbalanced data (uneven sampling)
weighted_contrasts <- function(factor_var) {
  n <- length(levels(factor_var))
  contrasts <- contr.sum(n)
  levels_weight <- table(factor_var) / length(factor_var)
  # Apply weights to each column of the contrasts matrix
  for (i in 1:ncol(contrasts)) {
    contrasts[, i] <- contrasts[, i] * sqrt(levels_weight) #Takes square root to stabilize variance when dealing with groups of different sizes (smaller families have larger variance)
  }
  return(contrasts)
}

contrasts(period_all$family) <- weighted_contrasts(period_all$family)
contrasts(period_all$period) <- weighted_contrasts(period_all$period)

m3 = lm(pdd ~ period*family, data = period_all) #(period*family)
m2 = lm(pdd ~ period+period:family, data = period_all)
m1 = lm(pdd ~ period+family, data = period_all) #best model (period+family)
m0 = lm(pdd ~ period, data = period_all)

ICtab(m3, m2, m1, m0) #m1 is the 'best' model
lmtest::lrtest(m1, m3)
anova(m1, m3)

mm1 = lmer(pdd ~ period + (1|family), data = period_all, REML=FALSE)
anova(mm1, m1) #m1 is the best fit model

summary(m1)
summary(contrast(emmeans::emmeans(m1, ~ period, type = 'response'), list('mutant - wt' = c(1, -1))), infer =c(TRUE, TRUE))
summary(contrast(emmeans::emmeans(m3, ~ period, type = 'response'), list('mutant - wt' = c(1, -1))), infer =c(TRUE, TRUE))

summary(m1)

summary(m3)
car::Anova(m3)
m3_emms = emmeans::emmeans(m3, ~ period | family , type = 'response')
m3_contrast = contrast(m3_emms, list('mutant - wt' = c(1, -1)), adjust = "BH", wts = NA)
contrast_summary = summary(m3_contrast, infer = c(TRUE, TRUE))
plot_data = as.data.frame(contrast_summary)

figure3b_r = ggplot(plot_data, aes(x = estimate, y = family, xmin = lower.CL, xmax = upper.CL)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
  labs(x = "Mutant PDD - Wildtype PDD", y = "Family") +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  theme_bw()  +
  theme(text = element_text(size = 14, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  scale_y_discrete(limits = rev) +
  annotate("text", x=-10.05, y=6.25, label= "**", size = 4)+
  annotate("text", x=-7.64, y=1.25, label= "***", size = 4)

figure3b_r


figure_3b = ggarrange(figure3b_l, figure3b_r, ncol = 2,
                      common.legend = TRUE, legend = 'top',
                      align = 'h')

ggsave("./plots/figure3b.png", plot = figure_3b, width = 6, height = 3, dpi = 600, units=c('in'),
       bg='white')
####################################################################################################################################################
#effect of rearing photoperiod and mutation on PDD time (interaction with diapause intensity?)
period_ld = period_all %>%
  group_by(family) %>%
  filter(n_distinct(ld) == n_distinct(period_all$ld)) %>%
  ungroup() %>% droplevels()

with(period_ld, table(family, period))
contrasts(period_ld$ld) <- weighted_contrasts(period_ld$ld)
contrasts(period_ld$family) <- weighted_contrasts(period_ld$family)
contrasts(period_ld$period) <- weighted_contrasts(period_ld$period)

period_ld %>% group_by(family, ld, period) %>%
  summarize(n = n(),
            pdd.mu = mean(pdd, na.rm=TRUE))

m3 = lm(pdd ~ period*ld*family, data = period_ld)
step(m3)
m2 = lm(pdd ~ period + ld + family + ld:family, data = period_ld) #best model
m2b = lm(pdd ~ period:ld + family, data = period_ld) #best model
m1 = lm(pdd ~ period +  ld + family, data = period_ld)
mm1 = lmer(pdd ~ period + (1|family) + (1|ld), data = period_ld, REML=FALSE)
m0 = lm(pdd ~ period +  ld, data = period_ld)
m0b = lm(pdd ~ period*ld, data = period_ld)
mm0 = lmer(pdd ~ period + ld + (1|family), data = period_ld, REML=FALSE)

anova(m2, m3)
ICtab(m3, m2, m2b, m1, mm1, m0, mm0)
anova(m1, m2) #m2 is the best model
anova(m0, m2) #m2 is the best model
ICtab(m0, m0b) #m0 is the more simple, best model
Anova(m0b)

###########
#effect of rearing photoperiod and mutation on PDD time (interaction with diapause intensity?)
car::Anova(m2)
m2_emms = emmeans::emmeans(m2, ~ period | ld*family , type = 'response') #results are averaged over the levels of family
multcomp_m2 = multcomp::cld(m2_emms, adjust = 'by', Letters = letters, decreasing=TRUE) %>%
  mutate(.group = gsub(" ", "", .group))

multcomp_m2$period = factor(multcomp_m2$period, c("wt", "mutant"))
period_ld$period = factor(period_ld$period, c("wt", "mutant"))

figureS4 = ggplot(data = multcomp_m2, aes(x=ld, y=emmean, color = period, group = period)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), position = position_dodge(width = 0.1),
                width = 0.2) +
  geom_point(position = position_dodge(width=0.1), size = 3) +
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 14, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  geom_text(aes(label = .group, x = ld, y = rep(50, 8), group=period),
            position=position_dodge(width=0.1),
            color = 'black', size = 6) +
  scale_colour_manual(values = c('mutant' = '#00CD00', 'wt' = "darkorchid"),
                      labels = c('wt' = expression(italic(per)^"+/+"),
                                 'mutant' = expression(italic(per)^"-/+")),
                      breaks = c('wt', 'mutant')) +
  scale_alpha_manual(values = c(0.5, 0.5)) +
  scale_y_continuous(name = "PDD Time (days)",
                     limits = c(0,45), breaks = seq(0,40,by=10)) +
  scale_x_discrete(breaks = c("12L12D", "15L09D"),
                   labels = c("12L12D" = "12L:12D",
                              "15L09D" = "15L:9D"),
                   name = "Photoperiod") +
  facet_wrap(~family) +
  annotate("text", x = 1, y = 42, label = "*", size = 4) +
  annotate("text", x = 2, y = 42, label = "*", size = 4)

figureS4
ggsave("./plots/figureS4.png", plot = figure4, width = 4, height = 4, dpi = 800, units=c('in'), bg='white')


##################################################################################################################
ecl_2024 = period_all %>% filter(year == 2024, !is.na(azt), ld == '12L12D')
ecl_2022 = period_all %>% filter(year == 2022, !is.na(azt), ld == '12L12D')

period_all %>% filter(!is.na(azt)) %>% 
  group_by(ld, year, period, family) %>% filter(n()>3) %>%
  summarise(
    eclosion_mu = mean.circular(azt, na.rm=TRUE)) %>% 
  pivot_wider(names_from = period, values_from = eclosion_mu) %>%
  mutate(eclosion_diff = wt - mutant)

period_all %>% filter(!is.na(azt)) %>% 
  group_by(year, period) %>% filter(n()>1) %>%
  summarise(
    eclosion_mu = mean.circular(azt, na.rm=TRUE)) %>% 
  pivot_wider(names_from = period, values_from = eclosion_mu) %>%
  mutate(eclosion_diff = wt - mutant)

period_all %>% filter(!is.na(azt)) %>% 
  group_by(year, period) %>%
  summarize(n = n())

#circle plot for eclosion grouping ALL samples together from each year
par(mfrow=c(1,2), mar = c(0, 0, 0, 0))
#plot 2022 data
plot.circular(period_all$azt[which(period_all$period == 'wt' & period_all$year == 2022)], stack = TRUE, pch = 20, sep = 0.05, cex=1.5, shrink = 1.3,
              ticks = FALSE, axes=FALSE, col = alpha('darkorchid', 0.50))
#title('DD\n(16L:8D)', line = -16)
arrows.circular(mean.circular(period_all$azt[which(period_all$period == 'wt' & period_all$year == 2022)], na.rm=TRUE), col = 'darkorchid')
points(period_all$azt[which(period_all$period == 'mutant' & period_all$year == 2022)], stack = TRUE, col = alpha('#00CD00', 0.5),
       pch = 20, sep = 0.05,cex=1.5)
arrows.circular(mean.circular(period_all$azt[which(period_all$period == 'mutant' & period_all$year == 2022)], na.rm=TRUE), col = alpha("#00CD00", 1))
axis.circular(at = circular(seq(0,20,by=4), units = 'hours', template = 'clock24'),
              labels = c("24/0", "4", "8", "12", "16", "20"), col = 'black')

#plot 2024 data
plot.circular(period_all$azt[which(period_all$period == 'wt' & period_all$year == 2024)], stack = TRUE, pch = 20, sep = 0.05, cex=1.5, shrink = 1.3,
              ticks = FALSE, axes=FALSE, col = alpha('darkorchid', 0.50))
#title('LD\n(12L:12D)', line = -16)
arrows.circular(mean.circular(period_all$azt[which(period_all$period == 'wt' & period_all$year == 2024)], na.rm=TRUE), col = 'darkorchid')
points(period_all$azt[which(period_all$period == 'mutant' & period_all$year == 2024)], stack = TRUE, col = alpha('#00CD00', 0.5),
       pch = 20, sep = 0.05,cex=1.5)
arrows.circular(mean.circular(period_all$azt[which(period_all$period == 'mutant' & period_all$year == 2024)], na.rm=TRUE), col = alpha("#00CD00", 1))
axis.circular(at = circular(seq(0,20,by=4), units = 'hours', template = 'clock24'),
              labels = c("24/0", "4", "8", "12", "16", "20"), col = 'black')

watson.two.test(ecl_2022$azt[ecl_2022$period == 'wt'],
                ecl_2022$azt[ecl_2022$period == 'mutant'])

watson.two.test(ecl_2024$azt[ecl_2024$period == 'wt'],
                ecl_2024$azt[ecl_2024$period == 'mutant'])
